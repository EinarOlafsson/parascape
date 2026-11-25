"""Infer E10: the spaCR OPS clustering ledger, repeated on ParaScape gene data.

The spaCR paper and spaCR master sessions benchmarked how to turn the OPS gene x feature table into clusters
(Funk plate A pilot, 2026-11-23; ledger at BMCspacr/reanalysis/a1_clustering_20261123). Their protocol:
  - change ONE step at a time against a baseline chain (a "ledger");
  - compare every arm at a MATCHED cluster count (Louvain resolution bisected to the target);
  - score known biology as EXCESS over a size-matched permutation null (cluster labels permuted among the
    clustered genes, so every cluster keeps its size);
  - split held-out sets complex-wise; keep a gene-shuffled profile as a negative control.
Their findings to test here: robust z collapses to the null; KNN vs median impute has no effect; whitening
(TVN) helps; SNF fusion of views hurts against concatenation; the method matters little at matched k;
consensus over replicate redraws is the biggest single gain.

ParaScape has no non-targeting controls and no guides, so two steps are translated:
  - TVN -> PCA whitening over all genes (sphering the typical variation);
  - guide-bootstrap consensus -> feature-bootstrap consensus (columns redrawn within each view, B redraws,
    co-association, average-linkage cut at the same k). A seed-only consensus separates redraw from averaging.

Scores (each as excess over the size-matched null, on the dev and the test half):
  - complexes: share of CORUM complexes (copied by orthology, >= 3 Toxoplasma genes clustered) with >= 50 % of
    their clustered members in one cluster; complexes are split into halves by complex;
  - STRING pairs: precision, recall and F1 of same-cluster calls (STRING copied by orthology, combined score
    >= 700, plus STRING measured in Toxoplasma); pairs are split in halves by connected component where
    possible, else by pair;
  - planted modules (positive control): the share of 20 planted 15-gene modules with >= 50 % in one cluster,
    run through the same chain.
Writes `ledger.parquet` (one row per arm x k x seed x split), `bootstrap.parquet` (paired differences against
the baseline with 95 % intervals), `summary.json`, and keeps a checkpoint per arm.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e3_preparation as E3  # noqa: E402
import e5_clustering as E5  # noqa: E402
import e9_strategies as E9  # noqa: E402

LABELS = "/media/carruthers/mnt2/parascape/labels/interaction_2026_11_25"
TARGET_K = (100, 300)
SEEDS = 3
N_NULL = 50
N_BOOT = 200
CONSENSUS_B = 10
STRING_MIN = 700
PLANT_STRENGTH = 2.0
SPLIT_SEED = 2026
MEASURED_VIEWS = ("transcript", "translation", "chromatin", "protein_abundance", "fitness", "place",
                  "interaction", "sequence", "structure_computed", "conservation")

ARMS = [
    # name, description
    ("baseline", "figure 6 chain: per view median fill, z, PCA <= 30; views side by side; cosine kNN 15 Louvain"),
    ("robust z", "per view median/IQR scaling instead of mean/SD (spaCR: collapsed to the null)"),
    ("rank-normal", "per view rank to normal scores instead of z"),
    ("KNN impute", "per view KNN fill instead of the median (spaCR: no effect)"),
    ("no per-view PCA", "z-scored columns side by side, no reduction"),
    ("equal view weight", "each view block scaled to the same total variance before joining"),
    ("whitened (TVN analogue)", "after joining: PCA to 99 % variance, each component to unit variance"),
    ("fused graph (SNF-like)", "one cosine kNN graph per view, averaged, Louvain on the fused graph (spaCR: SNF hurt)"),
    ("measured views only", "drop the sequence-derived ESM, 3Di and secondary-structure views"),
    ("k-means", "clustering: k-means at k"),
    ("Ward", "clustering: Ward linkage at k"),
    ("average linkage, cosine", "clustering: average linkage on cosine distance at k"),
    ("spectral, cosine kNN 15", "clustering: spectral on the cosine kNN graph at k (spaCR dev winner)"),
    ("consensus, seeds", f"{CONSENSUS_B} Louvain seeds, co-association, average-linkage cut at k"),
    ("consensus, feature bootstrap", f"{CONSENSUS_B} column redraws within views, Louvain each, co-association cut"),
    ("whitened + consensus, seeds", "combination: whitening, then 10-seed Louvain consensus"),
    ("no per-view PCA + whitened + consensus, seeds", "combination: z-scored columns, whitening, 10-seed consensus"),
    ("shuffled genes (negative control)", "every column permuted over genes independently, then the baseline"),
]


# ------------------------------------------------------------------ reference sets and the split
def reference(nodes: pd.DataFrame, log) -> dict:
    """Complexes (index arrays) and STRING pairs (index pairs), each with a dev/test half."""
    idx = {g: i for i, g in enumerate(nodes["gene_id"].astype(str))}
    cx = pd.read_parquet(os.path.join(LABELS, "corum_complexes_mapped.parquet"))
    comps = []
    for _, r in cx.iterrows():
        m = sorted({idx[g] for g in r["toxo_genes"] if g in idx})
        if len(m) >= 3:
            comps.append((str(r["complex_name"]), np.array(m)))
    # identical member sets (the same Toxoplasma genes behind several human complexes) count once
    seen, uniq = set(), []
    for name, m in comps:
        key = tuple(m)
        if key not in seen:
            seen.add(key); uniq.append((name, m))
    rng = np.random.default_rng(SPLIT_SEED)
    half = rng.permutation(len(uniq)) % 2
    s1 = pd.read_parquet(os.path.join(LABELS, "pairs_string_copied_by_orthology.parquet"))
    s2 = pd.read_parquet(os.path.join(LABELS, "pairs_string_measured_in_toxoplasma.parquet"))
    s = pd.concat([s1[s1["combined_score"] >= STRING_MIN], s2])
    a = s["gene_a"].map(idx); b = s["gene_b"].map(idx)
    ok = a.notna() & b.notna() & (a != b)
    P = np.unique(np.sort(np.c_[a[ok].astype(int), b[ok].astype(int)], axis=1), axis=0)
    n = len(nodes)
    g = csr_matrix((np.ones(len(P)), (P[:, 0], P[:, 1])), shape=(n, n))
    ncomp, lab = connected_components(g, directed=False)
    sizes = np.bincount(lab[P[:, 0]], minlength=ncomp)
    if sizes.max() <= 0.5 * len(P):
        comp_half = rng.permutation(ncomp) % 2
        p_half = comp_half[lab[P[:, 0]]]
        how = "by connected component"
    else:
        p_half = rng.permutation(len(P)) % 2
        how = f"by pair (largest component holds {sizes.max()} of {len(P)} pairs)"
    log(f"reference: {len(uniq)} complexes with >= 3 genes, {len(P)} STRING pairs, pairs split {how}")
    return {"complexes": uniq, "complex_half": half, "pairs": P, "pair_half": p_half, "pair_split": how}


# ------------------------------------------------------------------ scoring
def complex_share(c: np.ndarray, comps: list) -> float:
    """Share of complexes with >= 50 % of their clustered members in one cluster (unclustered = -1 excluded)."""
    hits = []
    for _, m in comps:
        cm = c[m]
        cm = cm[cm >= 0]
        if len(cm) < 3:
            continue
        hits.append(np.bincount(cm).max() / len(cm) >= 0.5)
    return float(np.mean(hits)) if hits else np.nan


def complex_hits(c: np.ndarray, comps: list) -> np.ndarray:
    out = []
    for _, m in comps:
        cm = c[m]
        cm = cm[cm >= 0]
        out.append(np.nan if len(cm) < 3 else float(np.bincount(cm).max() / len(cm) >= 0.5))
    return np.array(out)


def pair_scores(c: np.ndarray, P: np.ndarray) -> tuple:
    """STRING precision (same-cluster pairs that are STRING pairs, among clustered genes), recall, F1."""
    ok = (c[P[:, 0]] >= 0) & (c[P[:, 1]] >= 0)
    Pk = P[ok]
    same = c[Pk[:, 0]] == c[Pk[:, 1]]
    tp = int(same.sum())
    sizes = np.bincount(c[c >= 0])
    called = float((sizes * (sizes - 1) / 2).sum())
    prec = tp / called if called else np.nan
    rec = tp / len(Pk) if len(Pk) else np.nan
    f1 = 2 * prec * rec / (prec + rec) if prec and rec else 0.0
    return prec, rec, f1


def size_matched_null(c: np.ndarray, seed: int) -> np.ndarray:
    """Cluster labels permuted among the clustered genes: every cluster keeps its size."""
    out = c.copy()
    on = np.flatnonzero(c >= 0)
    out[on] = np.random.default_rng(seed).permutation(c[on])
    return out


def score(c: np.ndarray, ref: dict, planted_y=None, n_null: int = N_NULL) -> list:
    """Excess over the size-matched null for each metric and split."""
    rows = []
    comps, half, P, ph = ref["complexes"], ref["complex_half"], ref["pairs"], ref["pair_half"]
    nulls = [size_matched_null(c, 9000 + s) for s in range(n_null)]
    for split, h in (("dev", 0), ("test", 1)):
        cs = [x for x, hh in zip(comps, half) if hh == h]
        Ps = P[ph == h]
        obs_c = complex_share(c, cs)
        null_c = np.nanmean([complex_share(z, cs) for z in nulls])
        pr, rc, f1 = pair_scores(c, Ps)
        nl = np.array([pair_scores(z, Ps) for z in nulls])
        rows.append({"split": split, "complexes": obs_c, "complexes_null": null_c,
                     "complexes_excess": obs_c - null_c,
                     "string_precision": pr, "string_recall": rc, "string_f1": f1,
                     "string_f1_null": float(np.nanmean(nl[:, 2])), "string_f1_excess": f1 - float(np.nanmean(nl[:, 2])),
                     "string_recall_excess": rc - float(np.nanmean(nl[:, 1]))})
    return rows


def planted_share(c: np.ndarray, y) -> tuple:
    groups = [np.flatnonzero(y == g) for g in pd.unique(y[pd.notna(y)])]
    comps = [("m", m) for m in groups]
    obs = complex_share(c, comps)
    null = np.nanmean([complex_share(size_matched_null(c, 9500 + s), comps) for s in range(20)])
    return obs, obs - null


# ------------------------------------------------------------------ the chain
def prep_view(X, cover, arm):
    if arm == "robust z":
        return E3.prepare(X, cover, scaling="robust")
    if arm == "rank-normal":
        return E3.prepare(X, cover, scaling="rank")
    if arm == "KNN impute":
        return E3.prepare(X, cover, missing="knn")
    if "no per-view PCA" in arm:
        Xc = X[cover].copy()
        fill = np.nanmedian(Xc, axis=0)
        fill = np.where(np.isfinite(fill), fill, 0.0)
        holes = np.isnan(Xc)
        Xc[holes] = np.take(fill, np.where(holes)[1])
        sd = Xc.std(axis=0); sd[sd == 0] = 1.0
        out = np.full((X.shape[0], X.shape[1]), np.nan)
        out[cover] = (Xc - Xc.mean(axis=0)) / sd
        return out
    return E3.prepare(X, cover)


def views_for(nodes, views, extra, arm, rng=None) -> dict:
    """(matrix, cover) per view; for the feature bootstrap, columns redrawn within each view."""
    raw = E3.views_raw(nodes, views, extra, frozenset())
    if arm == "measured views only":
        raw = {k: v for k, v in raw.items() if k in MEASURED_VIEWS}
    if rng is not None:
        raw = {k: (X[:, rng.integers(0, X.shape[1], X.shape[1])], cov) for k, (X, cov) in raw.items()}
    return raw


def blocks_for(raw: dict, arm: str) -> list:
    return [prep_view(X, cov, arm) for X, cov in raw.values()]


def join(blocks: list, arm: str):
    """Blocks side by side; genes with >= 50 % of the columns kept, missing blocks at 0 (the mean)."""
    bl = []
    for B in blocks:
        if arm == "equal view weight":
            tot = np.nansum(np.nanvar(B, axis=0))
            B = B / np.sqrt(tot) if tot > 0 else B
        bl.append(B)
    E = np.hstack(bl)
    covered = (~np.isnan(E)).mean(axis=1) >= 0.5
    E = np.nan_to_num(E[covered])
    if "whitened" in arm:
        E = PCA(n_components=0.99, whiten=True, random_state=0).fit_transform(E)
    return E, np.flatnonzero(covered)


def knn_graph(E, k=E5.K_GRAPH):
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(E)
    _, ind = nn.kneighbors(E)
    rows = np.repeat(np.arange(len(E)), k)
    A = csr_matrix((np.ones(len(rows)), (rows, ind[:, 1:].ravel())), shape=(len(E), len(E)))
    return ((A + A.T) > 0).astype(float)


def louvain_graph(A, resolution, seed):
    import networkx as nx
    from networkx.algorithms.community import louvain_communities
    g = nx.from_scipy_sparse_array(A)
    comms = louvain_communities(g, weight="weight", resolution=resolution, seed=seed)
    out = np.empty(A.shape[0], int)
    for c, m in enumerate(comms):
        out[list(m)] = c
    return out


def louvain_at_k(A, k, seed, log=None):
    """Louvain with the resolution bisected (in log space) to give k clusters (within 5 %)."""
    lo, hi = np.log(0.05), np.log(200.0)
    best = None
    for _ in range(18):
        mid = (lo + hi) / 2
        c = louvain_graph(A, float(np.exp(mid)), seed)
        nk = len(np.unique(c))
        if best is None or abs(nk - k) < abs(best[1] - k):
            best = (c, nk, float(np.exp(mid)))
        if abs(nk - k) <= 0.05 * k:
            break
        if nk < k:
            lo = mid
        else:
            hi = mid
    return best


def coassoc_cut(labelings: list, k: int) -> np.ndarray:
    """Average linkage on 1 - co-association, cut at k clusters."""
    n = len(labelings[0])
    C = np.zeros((n, n), np.float32)
    for lab in labelings:
        C += (lab[:, None] == lab[None, :])
    D = 1.0 - C / len(labelings)
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    return fcluster(Z, t=k, criterion="maxclust") - 1


def cluster_arm(nodes, views, extra, arm, k, seed):
    """Cluster labels for every gene (-1 = not clustered), the cluster count and the resolution."""
    n = len(nodes)
    if arm == "consensus, feature bootstrap":
        labs, rows0 = [], None
        for b in range(CONSENSUS_B):
            E, rows = join(blocks_for(views_for(nodes, views, extra, arm, np.random.default_rng(500 + 31 * seed + b)), arm), arm)
            full = np.full(n, -1)
            full[rows] = louvain_at_k(knn_graph(E), k, seed + b)[0]
            labs.append(full)
        on = np.flatnonzero(np.all([lab >= 0 for lab in labs], axis=0))
        c = coassoc_cut([lab[on] for lab in labs], k)
        full = np.full(n, -1); full[on] = c
        return full, len(np.unique(c)), np.nan
    raw = views_for(nodes, views, extra, arm)
    if arm == "fused graph (SNF-like)":
        blocks = blocks_for(raw, arm)
        E, rows = join(blocks, arm)
        cov_all = np.zeros(n, bool); cov_all[rows] = True
        A = None
        for B in blocks:
            okb = ~np.isnan(B).any(axis=1)
            sub = np.flatnonzero(okb & cov_all)
            G = knn_graph(np.nan_to_num(B[sub]))
            pos = np.searchsorted(rows, sub)
            Gi = csr_matrix((G.data, (pos[G.nonzero()[0]], pos[G.nonzero()[1]])), shape=(len(rows), len(rows)))
            A = Gi if A is None else A + Gi
        A = A / len(blocks)
        c, nk, res = louvain_at_k(A, k, seed)
    else:
        E, rows = join(blocks_for(raw, arm), arm)
        if arm == "k-means":
            c = KMeans(n_clusters=k, n_init=4, random_state=seed).fit_predict(E); nk, res = k, np.nan
        elif arm == "Ward":
            c = AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(E); nk, res = k, np.nan
        elif arm == "average linkage, cosine":
            c = AgglomerativeClustering(n_clusters=k, linkage="average", metric="cosine").fit_predict(E); nk, res = k, np.nan
        elif arm == "spectral, cosine kNN 15":
            c = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=seed,
                                   assign_labels="cluster_qr").fit_predict(knn_graph(E)); nk, res = k, np.nan
        elif "consensus, seeds" in arm:
            A = knn_graph(E)
            labs = [louvain_at_k(A, k, seed * 100 + b)[0] for b in range(CONSENSUS_B)]
            c = coassoc_cut(labs, k); nk, res = len(np.unique(c)), np.nan
        else:
            c, nk, res = louvain_at_k(knn_graph(E), k, seed)
    full = np.full(n, -1)
    full[rows] = c
    return full, nk, res


def shuffled_nodes(nodes, views, extra, seed):
    """Every feature column permuted over genes independently (the gene-shuffled negative control)."""
    rng = np.random.default_rng(seed)
    out = nodes.copy()
    for cols in views.values():
        for col in cols:
            out[col] = rng.permutation(out[col].to_numpy())
    ex = {k: X[rng.permutation(len(X))] for k, X in extra.items()}
    return out, ex


# ------------------------------------------------------------------ paired bootstrap against the baseline
def paired_bootstrap(hits: dict, ref: dict) -> pd.DataFrame:
    """Difference in complex share (paired over complexes) and STRING recall (paired over pairs), test half."""
    rows = []
    rng = np.random.default_rng(11)
    half = ref["complex_half"]
    for (k, seed), base in hits.get("baseline", {}).items():
        bc, bp = base
        for arm, d in hits.items():
            if arm == "baseline" or (k, seed) not in d:
                continue
            ac, ap = d[(k, seed)]
            m = (half == 1) & np.isfinite(ac) & np.isfinite(bc)
            diff_c = ac[m] - bc[m]
            mp = ref["pair_half"] == 1
            diff_p = ap[mp].astype(float) - bp[mp].astype(float)
            bc_ = [diff_c[rng.integers(0, len(diff_c), len(diff_c))].mean() for _ in range(N_BOOT)] if len(diff_c) else [np.nan]
            bp_ = [diff_p[rng.integers(0, len(diff_p), len(diff_p))].mean() for _ in range(N_BOOT)]
            rows.append({"arm": arm, "k": k, "seed": seed,
                         "complexes_diff": float(diff_c.mean()) if len(diff_c) else np.nan,
                         "complexes_lo": float(np.nanpercentile(bc_, 2.5)), "complexes_hi": float(np.nanpercentile(bc_, 97.5)),
                         "string_recall_diff": float(diff_p.mean()),
                         "string_lo": float(np.percentile(bp_, 2.5)), "string_hi": float(np.percentile(bp_, 97.5))})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=SEEDS)
    ap.add_argument("--k", default=",".join(map(str, TARGET_K)))
    ap.add_argument("--arms", default="", help="';'-separated arm names (default: all)")
    ap.add_argument("--plant-strength", type=float, default=PLANT_STRENGTH)
    ap.add_argument("--calibrate", default="", help="comma-separated planted strengths: baseline only, k 100, then exit")
    ap.add_argument("--quick", action="store_true", help="1 seed, k 100, 5 nulls, 3 consensus draws: a code check")
    args = ap.parse_args(argv)
    global N_NULL, CONSENSUS_B
    if args.quick:
        N_NULL, CONSENSUS_B = 5, 3
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    names = [a for a, _ in ARMS]
    unknown = [a for a in args.arms.split(";") if a and a not in names]
    if unknown:
        raise SystemExit(f"unknown arms {unknown}; known: {names}")
    arms = [a for a in names if not args.arms or a in args.arms.split(";")]
    ks = [int(x) for x in args.k.split(",")]
    seeds = 1 if args.quick else args.seeds
    nodes, views, extra, _L = H.load()
    ref = reference(nodes, log)
    if args.calibrate:
        cal = []
        for st in [float(x) for x in args.calibrate.split(",")]:
            pn, py, _ = E9.planted(nodes, views, st)
            for k in ks:
                cp = cluster_arm(pn, views, extra, "baseline", k, 0)[0]
                obs, ex = planted_share(cp, py)
                cal.append({"strength": st, "k": k, "planted": obs, "planted_excess": ex})
                log(f"calibration: strength {st} k {k}: planted modules recovered {obs:.2f} (excess {ex:+.2f})")
        pd.DataFrame(cal).to_csv(os.path.join(args.out, "calibration.csv"), index=False)
        return 0
    pn, py, pcols = E9.planted(nodes, views, args.plant_strength)
    log(f"planted control: {len(pd.unique(py[pd.notna(py)]))} modules in {len(pcols)} columns at {args.plant_strength} SD")
    part = os.path.join(args.out, "ledger.partial.parquet")
    rows = pd.read_parquet(part).to_dict("records") if os.path.exists(part) else []
    done = {(r["arm"], r["k"], r["seed"]) for r in rows}
    hits_path = os.path.join(args.out, "hits.partial.npz")
    hits: dict = {}
    if os.path.exists(hits_path):
        z = np.load(hits_path, allow_pickle=True)
        hits = z["hits"].item()
    for arm in arms:
        for k in ks:
            for seed in range(seeds):
                if (arm, k, seed) in done:
                    continue
                t1 = time.time()
                if arm.startswith("shuffled genes"):
                    sn, sx = shuffled_nodes(nodes, views, extra, 700 + seed)
                    c, nk, res = cluster_arm(sn, views, sx, "baseline", k, seed)
                    cp = None
                else:
                    c, nk, res = cluster_arm(nodes, views, extra, arm, k, seed)
                    cp = cluster_arm(pn, views, extra, arm, k, seed)[0]
                sc = score(c, ref)
                pl, pl_ex = planted_share(cp, py) if cp is not None else (np.nan, np.nan)
                for r in sc:
                    rows.append({"arm": arm, "k": k, "seed": seed, "clusters": int(nk), "resolution": res,
                                 "genes_clustered": int((c >= 0).sum()), "planted": pl, "planted_excess": pl_ex,
                                 "seconds": round(time.time() - t1, 1), **r})
                same = (c[ref["pairs"][:, 0]] == c[ref["pairs"][:, 1]]) & (c[ref["pairs"][:, 0]] >= 0)
                hits.setdefault(arm, {})[(k, seed)] = (complex_hits(c, ref["complexes"]), same)
                pd.DataFrame(rows).to_parquet(part, index=False)
                np.savez(hits_path, hits=np.array(hits, dtype=object))
                t = [r for r in rows if r["arm"] == arm and r["k"] == k and r["seed"] == seed and r["split"] == "test"][0]
                log(f"{arm} k {k} seed {seed}: {nk} clusters; test complexes excess {t['complexes_excess']:+.3f}, "
                    f"STRING F1 excess {t['string_f1_excess']:+.4f}, planted {pl:.2f} ({time.time() - t1:.0f}s)")
    L = pd.DataFrame(rows)
    L.to_parquet(os.path.join(args.out, "ledger.parquet"), index=False)
    B = paired_bootstrap(hits, ref)
    B.to_parquet(os.path.join(args.out, "bootstrap.parquet"), index=False)
    agg = (L.groupby(["arm", "k", "split"])[["clusters", "complexes_excess", "string_f1_excess",
                                             "string_recall_excess", "planted", "planted_excess"]]
           .mean().round(4).reset_index())
    agg.to_csv(os.path.join(args.out, "ledger_mean.csv"), index=False)
    json.dump({"arms": dict(ARMS), "k": ks, "seeds": seeds, "n_null": N_NULL, "consensus_b": CONSENSUS_B,
               "string_min": STRING_MIN, "plant_strength": args.plant_strength, "pair_split": ref["pair_split"],
               "complexes": len(ref["complexes"]), "pairs": int(len(ref["pairs"])),
               "source_protocol": "spaCR a1_clustering_20261123 ledger + cluster_bench (spaCR paper, spaCR master)",
               "table": agg.to_dict("records")},
              open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    for f in (part, hits_path):
        if os.path.exists(f):
            os.remove(f)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
