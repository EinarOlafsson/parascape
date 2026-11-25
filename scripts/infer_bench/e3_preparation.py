"""E3 of the Infer benchmark: how a view is prepared before neighbours are found.

Design: parascape-paper `infer_benchmark/DESIGN.md`, section 5. Every catalogue label x every view,
scored as in E1 (family-held-out 5-fold, `SEEDS` seeds, kNN k = 15 on cosine), under:

- **missing values**: `median` (figure 6), `mean`, `knn` (5 nearest genes within the view), and
  `mask` (median plus one indicator column per partly missing column);
- **scaling**: `z` (figure 6), `robust` (median and IQR), `rank` (quantile to normal);
- **view weighting**, for all views at once: `per_view` (each view's block scaled to unit total
  variance, so a 320-column ESM view counts as much as a 2-column view) against `per_column`
  (blocks concatenated as they are).

The preparation is fitted on all genes but uses no label, so it cannot leak the label. Each
variant changes one choice from figure 6 at a time (missing x z, then median x scaling).

Writes `preparation.parquet` and `summary.json`.

    python scripts/infer_bench/e3_preparation.py --out DIR [--seeds 3] [--labels "a;b"]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import iqr, norm
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402

IM = H.IM
SEEDS = 3
MISSING = ("median", "mean", "knn", "mask")
SCALING = ("z", "robust", "rank")


def prepare(X: np.ndarray, cover: np.ndarray, missing: str = "median", scaling: str = "z") -> np.ndarray:
    """One view prepared as figure 6 does (`IM.embed`), with one choice changed; rows outside `cover` NaN."""
    Xc = X[cover].copy()
    holes = np.isnan(Xc)
    if missing == "knn" and holes.any():
        Xc = KNNImputer(n_neighbors=5).fit_transform(Xc)
    else:
        fill = np.nanmean(Xc, axis=0) if missing == "mean" else np.nanmedian(Xc, axis=0)
        fill = np.where(np.isfinite(fill), fill, 0.0)
        Xc[holes] = np.take(fill, np.where(holes)[1])
    if scaling == "rank":
        r = pd.DataFrame(Xc).rank(pct=True).to_numpy()
        Xc = norm.ppf(np.clip(r - 0.5 / len(r), 1e-6, 1 - 1e-6))
    elif scaling == "robust":
        med, spread = np.median(Xc, axis=0), iqr(Xc, axis=0)
        spread[spread == 0] = 1.0
        Xc = (Xc - med) / spread
    else:
        sd = Xc.std(axis=0)
        sd[sd == 0] = 1.0
        Xc = (Xc - Xc.mean(axis=0)) / sd
    if missing == "mask":
        partly = holes.any(axis=0) & ~holes.all(axis=0)
        if partly.any():
            Xc = np.c_[Xc, holes[:, partly].astype(float)]
    k = min(30, Xc.shape[1], Xc.shape[0] - 1)
    if Xc.shape[1] > k:
        Xc = PCA(n_components=k, random_state=0).fit_transform(Xc)
    out = np.full((X.shape[0], Xc.shape[1]), np.nan)
    out[cover] = Xc
    return out


def views_raw(nodes, views, extra, drop) -> dict:
    """(matrix, cover) per view, before preparation, as figure 6 selects them."""
    out = {}
    for name, cols in views.items():
        keep = [c for c in cols if c not in drop]
        if len(keep) < 2:
            continue
        X = nodes[keep].to_numpy(float)
        cover = (~np.isnan(X)).mean(axis=1) >= 0.5
        if cover.sum() >= 100:
            out[name] = (X, cover)
    for name, X in extra.items():
        cover = ~np.isnan(X).any(axis=1)
        if cover.sum() >= 100:
            out[name] = (X, cover)
    return out


def combined(prepped: dict, weighting: str) -> np.ndarray:
    """All views side by side, on the genes every view covers."""
    blocks = []
    for E in prepped.values():
        B = E.copy()
        if weighting == "per_view":
            total = np.nansum(np.nanvar(B, axis=0))
            B = B / np.sqrt(total) if total > 0 else B
        blocks.append(B)
    return np.hstack(blocks)


def main(argv=None) -> int:
    """Score every preparation variant for every label x view, and the two weightings."""
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
    variants = [(m, "z") for m in MISSING] + [("median", s) for s in SCALING if s != "z"]
    rows = []
    for lname, (y, kind, tcols) in L.items():
        raw = views_raw(nodes, views, extra, IM.exclusions(tcols, views))
        for missing, scaling in variants:
            prepped = {v: prepare(X, c, missing, scaling) for v, (X, c) in raw.items()}
            for vname, E in prepped.items():
                s = [H.cv_score(E, y, kind, H.family_folds(fam, H.N_FOLDS, k))[0] for k in range(args.seeds)]
                rows.append({"label": lname, "view": vname, "missing": missing, "scaling": scaling,
                             "weighting": "", "score": np.nanmean(s), "score_sd": np.nanstd(s, ddof=1)})
            if (missing, scaling) == ("median", "z"):
                for weighting in ("per_view", "per_column"):
                    E = combined(prepped, weighting)
                    s = [H.cv_score(E, y, kind, H.family_folds(fam, H.N_FOLDS, k))[0] for k in range(args.seeds)]
                    rows.append({"label": lname, "view": "all views", "missing": missing, "scaling": scaling,
                                 "weighting": weighting, "score": np.nanmean(s), "score_sd": np.nanstd(s, ddof=1)})
        log(f"{lname}: done")
    t = pd.DataFrame(rows)
    base = t[(t.missing == "median") & (t.scaling == "z") & (t.weighting == "")].set_index(["label", "view"])["score"]
    t["change_from_figure6"] = t["score"] - [base.get((a, b), np.nan) for a, b in zip(t["label"], t["view"])]
    t.to_parquet(os.path.join(args.out, "preparation.parquet"), index=False)
    per = t[t.weighting == ""].groupby(["missing", "scaling"]).agg(
        median_change=("change_from_figure6", "median"),
        better_beyond_noise=("change_from_figure6", lambda c: int((c > 2 * t.loc[c.index, "score_sd"]).sum())),
        worse_beyond_noise=("change_from_figure6", lambda c: int((c < -2 * t.loc[c.index, "score_sd"]).sum())),
        cells=("change_from_figure6", "size")).reset_index()
    w = t[t.weighting != ""].pivot_table(index="label", columns="weighting", values="score")
    summary = {"variants": per.to_dict("records"), "weighting": w.round(4).to_dict()}
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
