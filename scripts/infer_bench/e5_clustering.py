"""E5 of the Infer benchmark: which clustering puts together the genes that belong together.

Design: parascape-paper `infer_benchmark/DESIGN.md`, section 7. On all views together (figure 6
preparation, cosine), for every categorical or binary catalogue label:

- **methods**: graph communities on the kNN graph (Louvain from networkx; Leiden is the design's
  choice but `leidenalg` is not in the env, and envs are not changed unasked), at resolution 0.5,
  1, 2 and 4 with k = 15; HDBSCAN at minimum cluster size 5, 10 and 20; k-means and a Gaussian
  mixture at k = classes x 1, 2 and 4; Ward agglomerative at the same k.
- **scores**:
  - held-out recovery: clusters are made on all genes (they use no label), labels of the training
    genes of each family fold find the clusters enriched for each class (hypergeometric, BH
    q < 0.05); the score is the share of held-out labelled genes that land in a cluster enriched
    for their own class, AVERAGED OVER CLASSES (chance is near 0, not the commonest class share);
  - each label is clustered on the views left after its exclusions (no place view for a place
    label), exactly as E1-E4 score it; the first run clustered every label on all views (fixed
    2026-11-25);
  - adjusted mutual information between clusters and the label, on labelled genes;
  - planted complexes: a separate run (planting changes the data), not in this script yet;
  - stability: adjusted Rand index between two seeds.

Writes `clustering.parquet` and `summary.json`.

    python scripts/infer_bench/e5_clustering.py --out DIR [--labels "a;b"]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from sklearn.cluster import HDBSCAN, AgglomerativeClustering, KMeans
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e3_preparation as E3  # noqa: E402

IM = H.IM
K_GRAPH = 15


def embedding(nodes, views, extra, drop=frozenset()) -> tuple:
    """All views prepared as figure 6 and side by side; (matrix of covered genes, their rows)."""
    raw = E3.views_raw(nodes, views, extra, drop)
    E = np.hstack([E3.prepare(X, c) for X, c in raw.values()])
    # Genes missing a view keep their other views: the missing blocks are set to 0 (the mean).
    covered = (~np.isnan(E)).mean(axis=1) >= 0.5
    E = np.nan_to_num(E)
    return E[covered], np.flatnonzero(covered)


def louvain(E: np.ndarray, resolution: float, seed: int) -> np.ndarray:
    """Communities of the cosine kNN graph."""
    import networkx as nx
    from networkx.algorithms.community import louvain_communities
    nn = NearestNeighbors(n_neighbors=K_GRAPH + 1, metric="cosine").fit(E)
    _, ind = nn.kneighbors(E)
    g = nx.Graph()
    g.add_nodes_from(range(len(E)))
    g.add_edges_from((i, int(j)) for i, row in enumerate(ind) for j in row[1:])
    comms = louvain_communities(g, resolution=resolution, seed=seed)
    out = np.empty(len(E), int)
    for c, members in enumerate(comms):
        out[list(members)] = c
    return out


def cluster(E: np.ndarray, method: str, setting, seed: int) -> np.ndarray:
    """Cluster labels for every row (HDBSCAN noise is -1)."""
    if method == "louvain":
        return louvain(E, setting, seed)
    if method == "hdbscan":
        return HDBSCAN(min_cluster_size=int(setting), copy=True).fit_predict(E)
    if method == "kmeans":
        return KMeans(n_clusters=int(setting), n_init=4, random_state=seed).fit_predict(E)
    if method == "gmm":
        return GaussianMixture(n_components=int(setting), covariance_type="diag",
                               random_state=seed).fit(E).predict(E)
    if method == "ward":
        return AgglomerativeClustering(n_clusters=int(setting), linkage="ward").fit_predict(E)
    raise ValueError(method)


def heldout_recovery(clusters: np.ndarray, y: np.ndarray, folds: np.ndarray) -> float:
    """Share of held-out labelled genes in a cluster enriched (BH q < 0.05) for their own class."""
    lab = pd.notna(y)
    per_class: dict = {}
    for k in np.unique(folds):
        tr = lab & (folds != k)
        te = lab & (folds == k)
        if not te.any():
            continue
        N = int(tr.sum())
        rows, ps = [], []
        for c in np.unique(clusters[clusters >= 0]):
            in_c = tr & (clusters == c)
            n = int(in_c.sum())
            if n == 0:
                continue
            for cls, K_ in pd.Series(y[tr]).value_counts().items():
                x = int((y[in_c] == cls).sum())
                if x:
                    rows.append((c, cls)); ps.append(hypergeom.sf(x - 1, N, K_, n))
        q = E1.bh(np.array(ps)) if ps else np.array([])
        enriched = {r for r, qq in zip(rows, q) if qq < 0.05}
        for i in np.flatnonzero(te):
            per_class.setdefault(y[i], []).append((clusters[i], y[i]) in enriched)
    # Averaged over classes: one giant cluster enriched for the commonest class otherwise "recovers"
    # most genes while finding nothing (HDBSCAN on effect on host scored 0.93 that way, 2026-11-25).
    return float(np.mean([np.mean(v) for v in per_class.values()])) if per_class else np.nan


def main(argv=None) -> int:
    """Cluster with every method and setting; score them against every label."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--labels", default="")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    nodes, views, extra, L = H.load()
    L.update(E1.literature_labels(nodes))
    L = {k: v for k, v in L.items() if v[1] != "continuous"}
    wanted = [x for x in args.labels.split(";") if x]
    if wanted:
        L = {k: v for k, v in L.items() if k in wanted}
    fam = H.families(nodes)
    ids = nodes["gene_id"].astype(str).to_numpy()
    by_drop = {}
    for lname, (_y, _k, tcols) in L.items():
        by_drop.setdefault(frozenset(IM.exclusions(tcols, views)), []).append(lname)
    n_classes = max(len(pd.unique(pd.Series(v[0]).dropna())) for v in L.values())
    settings = ([("louvain", r) for r in (0.5, 1, 2, 4)] + [("hdbscan", m) for m in (5, 10, 20)]
                + [(m, n_classes * f) for m in ("kmeans", "gmm", "ward") for f in (1, 2, 4)])
    out = []
    for drop, names in by_drop.items():
      E, rows = embedding(nodes, views, extra, drop)
      log(f"embedding without {len(drop)} columns for {names}: {E.shape}")
      for method, setting in settings:
        try:
            c0 = cluster(E, method, setting, 0)
            c1 = cluster(E, method, setting, 1)
        except Exception as error:                          # one failing method must not stop the run
            log(f"{method} {setting}: {error}")
            continue
        full = np.full(len(nodes), -1)
        full[rows] = c0
        stability = adjusted_rand_score(c0, c1)
        n_clusters = len(set(c0[c0 >= 0]))
        for lname in names:
            y, kind, _t = L[lname]
            yr = np.where(np.isin(np.arange(len(nodes)), rows), y, None).astype(object)
            rec = np.nanmean([heldout_recovery(full, yr, H.family_folds(fam, H.N_FOLDS, s)) for s in range(3)])
            null = np.nanmean([heldout_recovery(full, H.shuffled(yr, 100 + s), H.family_folds(fam, H.N_FOLDS, s))
                               for s in range(3)])
            lab = pd.notna(yr) & (full >= 0)
            ami = adjusted_mutual_info_score(yr[lab].astype(str), full[lab]) if lab.sum() > 10 else np.nan
            out.append({"method": method, "setting": setting, "label": lname, "clusters": n_clusters,
                        "heldout_recovery": rec, "shuffled_recovery": null, "recovery_gain": rec - null,
                        "ami": ami, "stability_ari": stability})
        log(f"{method} {setting}: {n_clusters} clusters, stability {stability:.2f}")
    t = pd.DataFrame(out)
    t.to_parquet(os.path.join(args.out, "clustering.parquet"), index=False)
    best = t.sort_values("recovery_gain").groupby("label").tail(1)
    summary = {"embedding_genes": int(len(rows)),
               "best_per_label": best[["label", "method", "setting", "heldout_recovery", "shuffled_recovery",
                                       "recovery_gain", "ami"]].to_dict("records"),
               "mean_over_labels": t.groupby(["method", "setting"])[["recovery_gain", "ami", "stability_ari"]]
               .mean().round(3).reset_index().to_dict("records")}
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
