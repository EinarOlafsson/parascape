"""Infer E11: optimise, per label, how features are extracted and how genes are clustered.

Einar 2026-11-25: "optimize clustering for each label, and optimize the extraction of important features; test
strategies, use known information as a positive control and optimize cluster quality; test the hyperparameter
space, the strategy space and the algorithm space."

Search space (random sampling, then successive halving):
  extraction  all columns | top k (30/100/300) by Kruskal-Wallis, mutual information, random-forest importance or
              L1 stability selection | Kruskal-Wallis significance set | L1 stability set (frequency >= 0.5);
              sequence/structure views (ESM, 3Di, secondary structure) in or out
  prepare     per view PCA off/10/30; equal view weight on/off; PCA whitening (99 %) on/off
  supervised  none | LDA scores appended | NCA metric learning (10 or 30 dims) - fitted on training genes only
  graph       kNN k 5/10/15/30/50, cosine or euclidean
  cluster     Louvain (resolution log-uniform 0.2-8; consensus over 1 or 5 seeds) | k-means | Ward (k log-uniform
              10-400) | HDBSCAN (min cluster 5/10/20/40, on 20 principal components)
Objective: held-out class recovery F1. Clusters enriched for a class among training genes (hypergeometric,
Benjamini-Hochberg over clusters x classes, q < 0.05) predict that class; F1 on validation genes, macro over
classes (binary labels: the positive class only).

Protocol (no information from the untouched genes reaches any choice):
  - genes split by gene family: development 80 % (E7 final split) and untouched 20 %;
  - inside development, 3 family folds: stage 1 scores N random configurations on fold 1, the best third goes on
    to folds 2-3, the best 8 by mean F1 get a stability score (5 redraws of the columns within views + new seed;
    size-weighted best Jaccard);
  - the chosen configuration maximises mean F1 x stability among those whose largest cluster holds <= 10 % of the
    genes; it is refitted on all development genes and graded once on the untouched genes;
  - references on the untouched genes: the E8 reference (all views, per view PCA 30, kNN 15 cosine, Louvain 0.5)
    and the E10 default (no per view PCA, whitening, Louvain consensus, k about 300);
  - negative control (--shuffled): the whole search run with development labels permuted; its untouched-gene F1
    shows how much an optimiser gains on nothing;
  - positive controls: ribosome / proteasome classes inside the labels (per-class F1), and (--planted) 20 planted
    modules of 15 genes at 8 column SD as their own label.

Writes trials.parquet (every evaluation), final.json (choice, test scores, references), per_class.csv; resumes
from trials.parquet.

    python scripts/infer_bench/e11_optimize.py --label "hyperLOPIT niche" --out DIR [--trials 120] [--shuffled]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom, kruskal
from sklearn.cluster import HDBSCAN, AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors, NeighborhoodComponentsAnalysis

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e7_inference as E7  # noqa: E402
import e9_strategies as E9  # noqa: E402
import e10_ops_ledger as E10  # noqa: E402

IM = H.IM
EXTRA_SEQ = ("esm", "secondary_structure", "3di")
TOP_K = (30, 100, 300)
RANKERS = ("kruskal", "mutual information", "forest", "l1 stability")


# ------------------------------------------------------------------ data
class Data:
    def __init__(self, nodes, views, extra, drop):
        self.nodes, self.views, self.extra = nodes, views, extra
        self.n = len(nodes)
        self.pool = []
        seen = set()
        for v, cols in views.items():
            for c in cols:
                if c not in drop and c not in seen and (~np.isnan(nodes[c].to_numpy(float))).sum() >= 100:
                    self.pool.append((v, c)); seen.add(c)
        self.view_of = dict((c, v) for v, c in self.pool)
        self.cols = [c for _v, c in self.pool]
        self.filled = E9.filled(nodes, self.cols)
        self._blocks: dict = {}

    def block(self, view, cols, pca, extra=False):
        """One view prepared: median fill, z-score, optional PCA; rows outside coverage NaN. Cached."""
        key = (view, extra, pca, hashlib.md5("|".join(cols).encode()).hexdigest())
        if key in self._blocks:
            return self._blocks[key]
        if extra:
            X = self.extra[view]
            cover = ~np.isnan(X).any(axis=1)
            Xc = X[cover].copy()
        else:
            X = self.nodes[cols].to_numpy(float)
            cover = (~np.isnan(X)).mean(axis=1) >= 0.5
            Xc = X[cover].copy()
            med = np.nanmedian(Xc, axis=0)
            med = np.where(np.isfinite(med), med, 0.0)
            holes = np.isnan(Xc)
            Xc[holes] = np.take(med, np.where(holes)[1])
        sd = Xc.std(axis=0); sd[sd == 0] = 1.0
        Xc = (Xc - Xc.mean(axis=0)) / sd
        if pca and Xc.shape[1] > pca and cover.sum() > pca:
            Xc = PCA(n_components=pca, random_state=0).fit_transform(Xc)
        out = np.full((self.n, Xc.shape[1]), np.nan, dtype=np.float32)
        out[cover] = Xc
        self._blocks[key] = out
        return out


# ------------------------------------------------------------------ feature rankings (training genes only)
def rankings(D: Data, y, train, seed=0) -> pd.DataFrame:
    lab = train & pd.notna(y)
    X = D.filled[lab]
    yl = np.asarray(y[lab]).astype(str)
    kr = []
    for j in range(X.shape[1]):
        try:
            kr.append(kruskal(*[X[yl == c, j] for c in np.unique(yl)]).pvalue)
        except ValueError:
            kr.append(1.0)
    mi = mutual_info_classif(X, yl, random_state=seed)
    fi = RandomForestClassifier(n_estimators=300, min_samples_leaf=5, n_jobs=2,
                                random_state=seed).fit(X, yl).feature_importances_
    Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1)
    rng = np.random.default_rng(seed)
    freq = np.zeros(X.shape[1])
    B = 20
    for b in range(B):                                     # stability selection (Meinshausen & Buehlmann 2010)
        idx = rng.choice(len(yl), len(yl) // 2, replace=False)
        if len(np.unique(yl[idx])) < 2:
            continue
        try:
            m = LogisticRegression(penalty="l1", C=0.05, solver="saga", max_iter=100, tol=1e-2)
        except TypeError:
            m = LogisticRegression(l1_ratio=1.0, C=0.05, solver="saga", max_iter=100, tol=1e-2)
        m.fit(Xs[idx], yl[idx])
        freq += (np.abs(m.coef_).max(axis=0) > 1e-8)
    t = pd.DataFrame({"view": [D.view_of[c] for c in D.cols], "feature": D.cols, "kruskal_p": kr,
                      "mutual information": mi, "forest": fi, "l1 stability": freq / B})
    t["kruskal"] = -np.log10(np.clip(t["kruskal_p"], 1e-300, 1))
    t["kruskal_q"] = E1.bh(t["kruskal_p"].to_numpy())
    return t


def select(R: pd.DataFrame, cfg) -> list:
    m = cfg["extract"]
    if m == "all":
        return list(R["feature"])
    if m == "kruskal significant":
        return R.loc[R["kruskal_q"] < 0.05, "feature"].tolist()
    if m == "l1 stability set":
        return R.loc[R["l1 stability"] >= 0.5, "feature"].tolist()
    ranker, k = m.rsplit(" top ", 1)
    return R.sort_values(ranker, ascending=False)["feature"].head(int(k)).tolist()


# ------------------------------------------------------------------ embedding
def embed(D: Data, R, cfg, y, train, col_rng=None):
    """Embedding of all genes; supervised steps fitted on `y` over `train` only. col_rng: keep a random 80 %
    of the selected columns (the stability redraw)."""
    feats = select(R, cfg)
    if col_rng is not None and len(feats) > 5:
        feats = list(col_rng.choice(feats, int(0.8 * len(feats)), replace=False))
    if len(D._blocks) > 400:                               # bound the block cache (memory)
        D._blocks.clear()
    by: dict = {}
    for c in feats:
        by.setdefault(D.view_of[c], []).append(c)
    blocks = []
    for v, cols in by.items():
        if len(cols) >= 1:
            blocks.append(D.block(v, cols, cfg["pca"] if len(cols) > 1 else 0))
    if cfg["seq_views"]:
        blocks += [D.block(v, [], cfg["pca"], extra=True) for v in D.extra]
    if not blocks:
        return None, None
    if cfg["equal_views"]:
        blocks = [b / np.sqrt(max(np.nansum(np.nanvar(b, axis=0)), 1e-9)) for b in blocks]
    E = np.hstack(blocks)
    rows = np.flatnonzero((~np.isnan(E)).mean(axis=1) >= 0.5)
    E = np.nan_to_num(E[rows]).astype(np.float64)
    if cfg["whiten"] and E.shape[1] > 2:
        E = PCA(n_components=0.99, whiten=True, random_state=0).fit_transform(E)
    sup = cfg["supervised"]
    if sup != "none":
        tr = train[rows] & pd.notna(y[rows])
        yl = np.asarray(y[rows][tr]).astype(str)
        if len(np.unique(yl)) >= 2 and tr.sum() > 20:
            if sup == "lda":
                lda = LinearDiscriminantAnalysis(n_components=min(len(np.unique(yl)) - 1, 20)).fit(E[tr], yl)
                P = lda.transform(E)
                P = P / max(P.std(), 1e-9) * np.sqrt(E.shape[1] / max(P.shape[1], 1)) * 0.5
                E = np.hstack([E, P])
            else:
                dims = int(sup.split()[1])
                base = PCA(n_components=min(50, E.shape[1]), random_state=0).fit(E)
                Eb = base.transform(E)
                nca = NeighborhoodComponentsAnalysis(n_components=min(dims, Eb.shape[1]), max_iter=50,
                                                     random_state=0).fit(Eb[tr], yl)
                E = nca.transform(Eb)
    return E, rows


# ------------------------------------------------------------------ clustering
def knn_graph(E, k, metric):
    nn = NearestNeighbors(n_neighbors=k + 1, metric=metric).fit(E)
    _, ind = nn.kneighbors(E)
    r = np.repeat(np.arange(len(E)), k)
    A = csr_matrix((np.ones(len(r)), (r, ind[:, 1:].ravel())), shape=(len(E), len(E)))
    return ((A + A.T) > 0).astype(float)


def cluster(E, cfg, seed=0):
    m = cfg["method"]
    if m == "louvain":
        A = knn_graph(E, cfg["knn"], cfg["metric"])
        labs = [E10.louvain_graph(A, cfg["resolution"], seed + 97 * b) for b in range(cfg["consensus"])]
        if len(labs) == 1:
            return labs[0]
        k = int(np.median([len(np.unique(lab)) for lab in labs]))
        return E10.coassoc_cut(labs, max(k, 2))
    if m == "kmeans":
        return KMeans(n_clusters=cfg["k"], n_init=3, random_state=seed).fit_predict(E)
    if m == "ward":
        return AgglomerativeClustering(n_clusters=cfg["k"], linkage="ward").fit_predict(E)
    if m == "hdbscan":
        Z = PCA(n_components=min(20, E.shape[1]), random_state=0).fit_transform(E)
        return HDBSCAN(min_cluster_size=cfg["min_size"], copy=True).fit_predict(Z)
    raise ValueError(m)


def sample_cfg(rng) -> dict:
    ext = ["all", "kruskal significant", "l1 stability set"] + [f"{r} top {k}" for r in RANKERS for k in TOP_K]
    cfg = {"extract": str(rng.choice(ext)), "seq_views": bool(rng.random() < 0.7),
           "pca": int(rng.choice([0, 5, 10, 30])), "equal_views": bool(rng.random() < 0.3),
           "whiten": bool(rng.random() < 0.5),
           "supervised": str(rng.choice(["none", "none", "lda", "nca 10", "nca 30"])),
           "method": str(rng.choice(["louvain", "louvain", "louvain", "kmeans", "ward", "hdbscan"])),
           "knn": int(rng.choice([5, 10, 15, 30, 50])), "metric": str(rng.choice(["cosine", "euclidean", "correlation"])),
           "resolution": float(np.exp(rng.uniform(np.log(0.2), np.log(8.0)))),
           "consensus": int(rng.choice([1, 1, 5])), "k": int(np.exp(rng.uniform(np.log(10), np.log(400)))),
           "min_size": int(rng.choice([5, 10, 20, 40]))}
    return cfg


# ------------------------------------------------------------------ scoring
def recovery_f1(c_full, y, kind, train, val) -> tuple:
    """Macro F1 over classes on validation genes of 'in a cluster enriched (training genes) for my class'."""
    tr = train & pd.notna(y) & (c_full >= 0)
    counts = pd.Series(y[tr]).value_counts()
    if len(counts) < 2:
        return np.nan, {}
    classes = [counts.index[-1]] if kind == "binary" else list(counts.index)
    N = int(tr.sum())
    keys, ps = [], []
    for k in np.unique(c_full[tr]):
        m = tr & (c_full == k)
        n = int(m.sum())
        for cls in classes:
            x = int((y[m] == cls).sum())
            if x:
                keys.append((k, cls)); ps.append(hypergeom.sf(x - 1, N, int(counts[cls]), n))
    q = E1.bh(np.array(ps)) if ps else np.array([])
    enr: dict = {}
    for (k, cls), qq in zip(keys, q):
        if qq < 0.05:
            enr.setdefault(cls, set()).add(k)
    va = val & pd.notna(y)
    per = {}
    for cls in classes:
        truth = va & (y == cls)
        pred = va & np.isin(c_full, list(enr.get(cls, set())))
        tp = int((truth & pred).sum())
        p = tp / pred.sum() if pred.sum() else 0.0
        r = tp / truth.sum() if truth.sum() else np.nan
        per[cls] = 2 * p * r / (p + r) if (p + r) and np.isfinite(r) else 0.0
    vals = [v for v in per.values() if np.isfinite(v)]
    return (float(np.mean(vals)) if vals else np.nan), per


def run_once(D, R, cfg, y, kind, train, val, seed=0, y_fit=None, col_rng=None):
    """Embed (supervised parts fitted on y_fit, default y), cluster, score F1 of y on val from train calls."""
    E, rows = embed(D, R, cfg, y if y_fit is None else y_fit, train, col_rng)
    if E is None or len(rows) < 50:
        return np.nan, {}, None
    c = cluster(E, cfg, seed)
    full = np.full(D.n, -1)
    full[rows] = c
    f1, per = recovery_f1(full, y, kind, train, val)
    on = full[full >= 0]
    sizes = np.bincount(on) if len(on) else np.array([0])
    return f1, per, {"clusters": int(len(np.unique(on))), "largest_share": float(sizes.max() / max(len(on), 1)),
                     "clustered": int(len(on)), "full": full}


def stability(D, R, cfg, y, train, ref_full, redraws=5):
    """Size-weighted best Jaccard of the reference clusters over column redraws (within views) + new seeds."""
    import cluster_quality as Q
    labs = []
    for b in range(redraws):
        E, rows = embed(D, R, cfg, y, train, col_rng=np.random.default_rng(1000 + b))
        if E is None:
            continue
        c = cluster(E, cfg, seed=11 + b)
        full = np.full(D.n, -1); full[rows] = c
        labs.append(full)
    if not labs:
        return np.nan
    ref = ref_full.copy()
    ok = ref >= 0
    per, _ = Q.stability(np.where(ok, ref, -1), labs)
    sizes = np.bincount(ref[ok], minlength=len(per))
    return float(np.average(per[:len(sizes)], weights=sizes)) if sizes.sum() else np.nan


# ------------------------------------------------------------------ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--trials", type=int, default=120)
    ap.add_argument("--shuffled", action="store_true")
    ap.add_argument("--planted", action="store_true", help="the label is 20 planted modules at 8 column SD")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quick", action="store_true", help="code check: 6 trials, 2 finalists, 2 redraws")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    if args.planted:
        nodes, y, _cols = E9.planted(nodes, views, 8.0)
        kind, drop = "categorical", set()
    else:
        y, kind, tcols = L[args.label]
        y = np.asarray(y, dtype=object)
        drop = set(IM.exclusions(tcols, views))
    if kind == "continuous":
        raise SystemExit("continuous labels have no classes to cluster")
    fam = H.families(nodes)
    fin = E7.final_split(fam)
    dev = ~fin
    y_search = y.copy()
    if args.shuffled:                                      # negative control: labels permuted inside development
        m = dev & pd.notna(y)
        y_search[m] = np.random.default_rng(4242).permutation(y_search[m])
    inner = H.family_folds(fam, 3, 7)
    folds = [(dev & (inner != f), dev & (inner == f)) for f in range(3)]
    D = Data(nodes, views, extra, drop)
    log(f"{args.label}{' [shuffled]' if args.shuffled else ''}: {len(D.cols)} features, "
        f"{int((pd.notna(y) & dev).sum())} development / {int((pd.notna(y) & fin).sum())} untouched labelled genes")
    Rf = [rankings(D, y_search, tr, seed=f) for f, (tr, _va) in enumerate(folds)]
    Rdev = rankings(D, y_search, dev, seed=9)
    log("rankings done")
    tpath = os.path.join(args.out, "trials.parquet")
    T = pd.read_parquet(tpath).to_dict("records") if os.path.exists(tpath) else []
    done = {(r["trial"], r["fold"]) for r in T}
    rng = np.random.default_rng(args.seed)
    n_trials = 6 if args.quick else args.trials
    cfgs = [sample_cfg(rng) for _ in range(n_trials)]

    def evaluate(i, f):
        if (i, f) in done:
            return
        tr, va = folds[f]
        t1 = time.time()
        try:
            f1, _per, info = run_once(D, Rf[f], cfgs[i], y_search, kind, tr, va, seed=f)
            info = {k: v for k, v in (info or {}).items() if k != "full"}
        except Exception as e:                              # a configuration that cannot run scores nothing
            f1, info = np.nan, {"error": str(e)[:200]}
        T.append({"trial": i, "fold": f, "f1": f1, "seconds": round(time.time() - t1, 1),
                  **{f"cfg_{k}": v for k, v in cfgs[i].items()}, **info})
        pd.DataFrame(T).to_parquet(tpath, index=False)

    for i in range(len(cfgs)):                              # stage 1: every configuration on fold 1
        evaluate(i, 0)
        if i % 10 == 0:
            log(f"stage 1: {i + 1}/{len(cfgs)}")
    t = pd.DataFrame(T)
    s1 = t[t["fold"] == 0].sort_values("f1", ascending=False)
    keep = s1["trial"].head(max(len(cfgs) // 3, 5)).tolist()
    for i in keep:                                          # stage 2: the best third on folds 2-3
        evaluate(i, 1); evaluate(i, 2)
    log("stage 2 done")
    t = pd.DataFrame(T)
    m = t[t["trial"].isin(keep)].groupby("trial").agg(f1=("f1", "mean"), f1_sd=("f1", "std"),
                                                       largest=("largest_share", "max"),
                                                       clusters=("clusters", "mean")).reset_index()
    top = m.sort_values("f1", ascending=False).head(2 if args.quick else 8)
    finals = []
    for _, r in top.iterrows():                             # stage 3: stability of the 8 best, fitted on development
        cfg = cfgs[int(r["trial"])]
        f1d, _p, info = run_once(D, Rdev, cfg, y_search, kind, dev, dev)
        st = stability(D, Rdev, cfg, y_search, dev, info["full"], redraws=2 if args.quick else 5) if info else np.nan
        finals.append({"trial": int(r["trial"]), "cv_f1": r["f1"], "cv_f1_sd": r["f1_sd"], "largest": r["largest"],
                       "clusters": r["clusters"], "stability": st, "score": r["f1"] * (st if np.isfinite(st) else 0)})
        log(f"  final candidate {int(r['trial'])}: cv F1 {r['f1']:.3f}, stability {st:.2f}, largest {r['largest']:.2f}")
    Fn = pd.DataFrame(finals)
    ok = Fn[(Fn["largest"] <= 0.10) & (Fn["stability"] >= 0.5)]       # balanced and not dissolved (Hennig 2007)
    if not len(ok):
        ok = Fn[Fn["largest"] <= 0.10]
    best = (ok if len(ok) else Fn).sort_values("score", ascending=False).iloc[0]
    cfg = cfgs[int(best["trial"])]
    # grade once on the untouched genes (true labels there, always)
    test_f1, per, info = run_once(D, Rdev, cfg, y, kind, dev, fin, y_fit=y_search)
    refs = {}
    ref_cfgs = {"E8 reference": {**sample_cfg(np.random.default_rng(0)), "extract": "all", "seq_views": True,
                                 "pca": 30, "equal_views": False, "whiten": False, "supervised": "none",
                                 "method": "louvain", "knn": 15, "metric": "cosine", "resolution": 0.5,
                                 "consensus": 1},
                "E10 default (resolution 4)": {**sample_cfg(np.random.default_rng(0)), "extract": "all", "seq_views": True,
                                "pca": 0, "equal_views": False, "whiten": True, "supervised": "none",
                                "method": "louvain", "knn": 15, "metric": "cosine", "resolution": 4.0,
                                "consensus": 5}}
    for name, rc in ref_cfgs.items():
        f_, p_, i_ = run_once(D, Rdev, rc, y, kind, dev, fin)
        refs[name] = {"test_f1": f_, "clusters": i_["clusters"] if i_ else None,
                      "largest_share": i_["largest_share"] if i_ else None, "per_class": p_}
    out = {"label": args.label, "shuffled": args.shuffled, "planted": args.planted, "trials": len(cfgs),
           "chosen_trial": int(best["trial"]), "chosen_config": cfg, "cv_f1": float(best["cv_f1"]),
           "stability": float(best["stability"]), "test_f1": test_f1,
           "test_clusters": info["clusters"] if info else None,
           "test_largest_share": info["largest_share"] if info else None, "references": refs,
           "finalists": Fn.to_dict("records"), "minutes": round((time.time() - t0) / 60, 1)}
    json.dump(out, open(os.path.join(args.out, "final.json"), "w"), indent=1, default=str)
    pd.DataFrame([{"class": k, "chosen_f1": v, **{f"{n}_f1": refs[n]["per_class"].get(k, np.nan) for n in refs}}
                  for k, v in per.items()]).to_csv(os.path.join(args.out, "per_class.csv"), index=False)
    log(f"chosen trial {int(best['trial'])}: cv F1 {best['cv_f1']:.3f}, stability {best['stability']:.2f}; "
        f"untouched F1 {test_f1:.3f} (E8 reference {refs['E8 reference']['test_f1']:.3f}, "
        f"E10 default {refs['E10 default (resolution 4)']['test_f1']:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
