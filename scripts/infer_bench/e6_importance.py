"""E6 of the Infer benchmark: do the ways of ranking what matters agree, and are they stable?

Design: parascape-paper `infer_benchmark/DESIGN.md`, section 3 ("Importance methods"). Per label,
the views are ranked by four methods (the per-feature ranking is E1's table):

1. **alone**: each view scored by itself (as E1, family folds);
2. **leave one view out**: the all-views score minus the score without that view;
3. **permutation**: in the all-views neighbour model, the held-out score lost when that view's
   block is shuffled across held-out genes (training genes untouched);
4. **random forest** (cross-check only): impurity importance of a forest on all prepared columns,
   summed per view.

Each method is run with `SEEDS` seeds. Reported: the Spearman agreement between methods (mean
ranks), and per method the stability, the mean Spearman between the rankings of two seeds. The
method kept is the most stable one whose top views also give the best held-out recovery (E7).

Writes `importance.parquet` (label x view x method x seed) and `summary.json`.

    python scripts/infer_bench/e6_importance.py --out DIR [--seeds 3] [--labels "a;b"]
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e3_preparation as E3  # noqa: E402

IM = H.IM
SEEDS = 3


def all_views(prepped: dict, skip: str = "") -> np.ndarray:
    """The prepared views side by side; genes missing a view keep the others (missing = 0)."""
    blocks = [E for v, E in prepped.items() if v != skip]
    E = np.hstack(blocks)
    covered = (~np.isnan(E)).mean(axis=1) >= 0.5
    E = np.nan_to_num(E)
    E[~covered] = np.nan
    return E


def permuted_loss(prepped, y, kind, folds, seed) -> dict:
    """Held-out score lost when one view's block is shuffled among the held-out genes."""
    names = list(prepped)
    E = all_views(prepped)
    widths = [prepped[v].shape[1] for v in names]
    starts = np.cumsum([0] + widths)
    ok = ~np.isnan(E).any(axis=1) & pd.notna(y)
    idx = np.flatnonzero(ok)
    yv = y[idx]
    classes = np.unique(yv) if kind != "continuous" else None
    rng = np.random.default_rng(seed)
    base_pred, perm_pred = {}, {v: {} for v in names}
    f = folds[idx]
    for k in np.unique(f):
        tr, te = idx[f != k], idx[f == k]
        m_tr = np.zeros(len(E), bool); m_tr[tr] = True
        m_te = np.zeros(len(E), bool); m_te[te] = True
        base_pred[k] = IM.knn_predict(E, y, kind, m_tr, m_te, classes)
        for v, a, b in zip(names, starts[:-1], starts[1:]):
            Ep = E.copy()
            Ep[te, a:b] = Ep[rng.permutation(te), a:b]
            perm_pred[v][k] = IM.knn_predict(Ep, y, kind, m_tr, m_te, classes)
    order = np.concatenate([idx[f == k] for k in np.unique(f)])
    pos = {g: i for i, g in enumerate(idx)}
    truth = yv[[pos[g] for g in order]]
    base = IM.score(truth, np.concatenate([base_pred[k] for k in np.unique(f)]), kind)
    return {v: base - IM.score(truth, np.concatenate([perm_pred[v][k] for k in np.unique(f)]), kind)
            for v in names}


def forest(prepped, y, kind, seed) -> dict:
    """Random-forest impurity importance summed per view (all labelled genes; a cross-check)."""
    names = list(prepped)
    E = all_views(prepped)
    ok = ~np.isnan(E).any(axis=1) & pd.notna(y)
    model = (RandomForestRegressor if kind == "continuous" else RandomForestClassifier)(
        n_estimators=200, min_samples_leaf=5, n_jobs=2, random_state=seed)
    model.fit(E[ok], y[ok].astype(float) if kind == "continuous" else y[ok].astype(str))
    imp = model.feature_importances_
    starts = np.cumsum([0] + [prepped[v].shape[1] for v in names])
    return {v: float(imp[a:b].sum()) for v, a, b in zip(names, starts[:-1], starts[1:])}


def main(argv=None) -> int:
    """Rank views by every method and seed; report agreement and stability."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=SEEDS)
    ap.add_argument("--labels", default="")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    nodes, views, extra, L = H.load()
    L.update(E1.literature_labels(nodes))
    wanted = [x for x in args.labels.split(";") if x]
    if wanted:
        L = {k: v for k, v in L.items() if k in wanted}
    fam = H.families(nodes)
    rows = []
    for lname, (y, kind, tcols) in L.items():
        raw = E3.views_raw(nodes, views, extra, IM.exclusions(tcols, views))
        prepped = {v: E3.prepare(X, c) for v, (X, c) in raw.items()}
        for s in range(args.seeds):
            folds = H.family_folds(fam, H.N_FOLDS, s)
            full = H.cv_score(all_views(prepped), y, kind, folds)[0]
            for v, E in prepped.items():
                rows.append({"label": lname, "view": v, "method": "alone", "seed": s,
                             "importance": H.cv_score(E, y, kind, folds)[0]})
                rows.append({"label": lname, "view": v, "method": "leave one view out", "seed": s,
                             "importance": full - H.cv_score(all_views(prepped, skip=v), y, kind, folds)[0]})
            for v, loss in permuted_loss(prepped, y, kind, folds, 100 + s).items():
                rows.append({"label": lname, "view": v, "method": "permutation", "seed": s, "importance": loss})
            for v, imp in forest(prepped, y, kind, s).items():
                rows.append({"label": lname, "view": v, "method": "random forest", "seed": s, "importance": imp})
        log(f"{lname}: done")
    t = pd.DataFrame(rows)
    t.to_parquet(os.path.join(args.out, "importance.parquet"), index=False)
    stab, agree = [], []
    for (lname, method), g in t.groupby(["label", "method"]):
        w = g.pivot(index="view", columns="seed", values="importance")
        rs = [spearmanr(w[a], w[b]).statistic for a, b in itertools.combinations(w.columns, 2)]
        stab.append({"label": lname, "method": method, "stability": float(np.nanmean(rs))})
    for lname, g in t.groupby("label"):
        m = g.groupby(["method", "view"])["importance"].mean().unstack(0)
        for a, b in itertools.combinations(m.columns, 2):
            agree.append({"label": lname, "a": a, "b": b, "spearman": float(spearmanr(m[a], m[b]).statistic)})
    S, A = pd.DataFrame(stab), pd.DataFrame(agree)
    summary = {"stability_by_method": S.groupby("method")["stability"].median().round(3).to_dict(),
               "agreement_between_methods": A.groupby(["a", "b"])["spearman"].median().round(3)
               .reset_index().to_dict("records")}
    S.to_csv(os.path.join(args.out, "stability.tsv"), sep="\t", index=False)
    A.to_csv(os.path.join(args.out, "agreement.tsv"), sep="\t", index=False)
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
