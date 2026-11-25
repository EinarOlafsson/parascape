"""E4 of the Infer benchmark: the similarity and the reduction used to find neighbours.

Design: parascape-paper `infer_benchmark/DESIGN.md`, section 6. Every catalogue label x every view
(and all views together), scored as in E1 (family-held-out 5-fold, `SEEDS` seeds, kNN k = 15),
under:

- **similarity**: `cosine` (figure 6), `euclidean`, `correlation`;
- **reduction**: `none` (all prepared columns), PCA to `10`, `30` (figure 6) or `50` components, or
  to `90%` of the variance; for all views together also `per_view` (each view reduced to at most 10
  components, then side by side).

Preparation is figure 6 (median fill, z-score), the E3 baseline; when E3 chooses another, rerun
with it. spaCR max-projection over views (items 992 / 1012) is added when 1012 lands
(one-implementation rule); until then the cross-view comparison is concatenation only.

Writes `similarity.parquet` and `summary.json`.

    python scripts/infer_bench/e4_similarity.py --out DIR [--seeds 3] [--labels "a;b"]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e3_preparation as E3  # noqa: E402

IM = H.IM
SEEDS = 3
METRICS = ("cosine", "euclidean", "correlation")
REDUCTIONS = ("none", "10", "30", "50", "90%")


def reduce(X: np.ndarray, cover: np.ndarray, how: str) -> np.ndarray:
    """A view prepared as figure 6 does, then reduced as asked."""
    full = E3.prepare(X, cover) if how == "30" else _prepared_unreduced(X, cover)
    if how in ("none", "30"):
        return full
    Z = full[cover]
    if how == "90%":
        p = PCA(random_state=0).fit(Z)
        k = int(np.searchsorted(np.cumsum(p.explained_variance_ratio_), 0.90) + 1)
    else:
        k = int(how)
    k = min(k, Z.shape[1], Z.shape[0] - 1)
    if Z.shape[1] > k:
        Z = PCA(n_components=k, random_state=0).fit_transform(Z)
    out = np.full((X.shape[0], Z.shape[1]), np.nan)
    out[cover] = Z
    return out


def _prepared_unreduced(X, cover):
    """Median fill and z-score, no PCA."""
    Xc = X[cover].copy()
    med = np.nanmedian(Xc, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    holes = np.isnan(Xc)
    Xc[holes] = np.take(med, np.where(holes)[1])
    sd = Xc.std(axis=0); sd[sd == 0] = 1.0
    Xc = (Xc - Xc.mean(axis=0)) / sd
    out = np.full((X.shape[0], Xc.shape[1]), np.nan)
    out[cover] = Xc
    return out


def knn_predict(E, y, kind, train, test, classes, metric):
    """IM.knn_predict with the metric as a choice."""
    nn = NearestNeighbors(n_neighbors=min(IM.K, train.sum() - 1), metric=metric)
    nn.fit(E[train])
    _, ind = nn.kneighbors(E[test])
    ytr = y[train]
    if kind == "continuous":
        return ytr[ind].mean(axis=1)
    return np.stack([(ytr[ind] == c).mean(axis=1) for c in classes], axis=1)


def cv_score(E, y, kind, folds, metric) -> float:
    """H.cv_score with the metric as a choice."""
    if metric == "cosine":
        return H.cv_score(E, y, kind, folds)[0]
    ok = ~np.isnan(E).any(axis=1) & pd.notna(y)
    idx = np.flatnonzero(ok)
    if len(idx) < 40:
        return np.nan
    yv = y[idx]
    classes = np.unique(yv) if kind != "continuous" else None
    preds = np.full((len(idx),) if kind == "continuous" else (len(idx), len(classes)), np.nan)
    f = folds[idx]
    for k in np.unique(f):
        tr, te = f != k, f == k
        if tr.sum() < IM.K + 1 or te.sum() == 0:
            continue
        m_tr = np.zeros(len(E), bool); m_tr[idx[tr]] = True
        m_te = np.zeros(len(E), bool); m_te[idx[te]] = True
        preds[te] = knn_predict(E, y, kind, m_tr, m_te, classes, metric)
    done = ~np.isnan(preds if preds.ndim == 1 else preds[:, 0])
    return IM.score(yv[done], preds[done], kind)


def main(argv=None) -> int:
    """Score every metric x reduction for every label x view, and all views together."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=SEEDS)
    ap.add_argument("--labels", default="")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    wanted = [x for x in args.labels.split(";") if x]
    if wanted:
        L = {k: v for k, v in L.items() if k in wanted}
    fam = H.families(nodes)
    rows = []

    def run(lname, vname, E, y, kind, metric, how):
        s = [cv_score(E, y, kind, H.family_folds(fam, H.N_FOLDS, k), metric) for k in range(args.seeds)]
        rows.append({"label": lname, "view": vname, "metric": metric, "reduction": how,
                     "dims": E.shape[1], "score": np.nanmean(s), "score_sd": np.nanstd(s, ddof=1)})

    for lname, (y, kind, tcols) in L.items():
        raw = E3.views_raw(nodes, views, extra, IM.exclusions(tcols, views))
        for how in REDUCTIONS:
            red = {v: reduce(X, c, how) for v, (X, c) in raw.items()}
            for metric in METRICS:
                for vname, E in red.items():
                    run(lname, vname, E, y, kind, metric, how)
                run(lname, "all views", np.hstack(list(red.values())), y, kind, metric, how)
        per_view = np.hstack([reduce(X, c, "10") for X, c in raw.values()])
        for metric in METRICS:
            run(lname, "all views", per_view, y, kind, metric, "per_view")
        log(f"{lname}: done")
    t = pd.DataFrame(rows)
    base = t[(t.metric == "cosine") & (t.reduction == "30")].set_index(["label", "view"])["score"]
    t["change_from_figure6"] = t["score"] - [base.get((a, b), np.nan) for a, b in zip(t["label"], t["view"])]
    t.to_parquet(os.path.join(args.out, "similarity.parquet"), index=False)
    agg = t.groupby(["metric", "reduction"]).agg(
        median_change=("change_from_figure6", "median"),
        better_beyond_noise=("change_from_figure6", lambda c: int((c > 2 * t.loc[c.index, "score_sd"]).sum())),
        worse_beyond_noise=("change_from_figure6", lambda c: int((c < -2 * t.loc[c.index, "score_sd"]).sum())),
        cells=("change_from_figure6", "size")).reset_index()
    best_all = (t[t.view == "all views"].sort_values("score").groupby("label").tail(1)
                [["label", "metric", "reduction", "score"]])
    json.dump({"variants": agg.to_dict("records"), "best_all_views": best_all.to_dict("records")},
              open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
