"""E2 of the Infer benchmark: which confounder correction, per feature x label.

Design: parascape-paper `infer_benchmark/DESIGN.md`, section 4. For every (label, feature) pair that
E1 found to carry information, the held-out score is compared under:

- **raw**: the feature as it is (E1);
- **residualised, linear**: the feature minus a linear fit on the confounder;
- **residualised, smooth**: the feature minus a cubic fit on the confounder's rank;
- **stratified**: scored within confounder quintiles, then averaged.

The residual fit is made on the TRAINING genes of each fold only and applied to the held-out genes,
so the choice cannot flatter the held-out score (DESIGN.md rule). Family-held-out 5-fold, `SEEDS`
seeds, kNN k = 15 on cosine, as in E0 and E1.

**Confounders** (node columns): protein length (`length`), peak expression (`expr_max`, log),
measured share (`heldout_sets.tsv`), papers naming the gene (focal + substantive + incidental),
codon GC3 (`codon_gc3`) and conservation breadth (`has_pf_ortholog` + `has_cp_ortholog`).
Cell-cycle periodicity has no column in this release and is left out.

**A label that depends on the confounder** (the confounder alone scores above `DEPENDS` on it) is
flagged: removing the confounder can then remove real signal, so both versions are reported and the
pair is not silently corrected.

Writes `confounders.parquet` (one row per label x feature x confounder x method) and `summary.json`.

    python scripts/infer_bench/e2_confounders.py --e1 E1_DIR --out DIR [--seeds 3] [--max-pairs N]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402

IM = H.IM
SEEDS = 3
QUINTILES = 5
#: A confounder that alone scores this far above chance on a label means the label depends on it
#: (AUROC above 0.5 + DEPENDS, or |rho| above DEPENDS).
DEPENDS = 0.10
HELDOUT = "/media/carruthers/mnt2/parascape/infer_bench/heldout/heldout_sets.tsv"


def confounders(nodes: pd.DataFrame) -> dict:
    """Each confounder as an array aligned to the node table."""
    ids = nodes["gene_id"].astype(str)
    out = {"protein length": nodes["length"].to_numpy(float),
           "peak expression": np.log1p(nodes["expr_max"].clip(lower=0).to_numpy(float)),
           "papers": (nodes["n_papers_focal"].fillna(0) + nodes["n_papers_substantive"].fillna(0)
                      + nodes["n_papers_incidental"].fillna(0)).to_numpy(float),
           "codon GC3": nodes["codon_gc3"].to_numpy(float),
           "conservation breadth": (nodes["has_pf_ortholog"].fillna(False).astype(float)
                                    + nodes["has_cp_ortholog"].fillna(False).astype(float)).to_numpy()}
    if os.path.exists(HELDOUT):
        h = pd.read_csv(HELDOUT, sep="\t").set_index("gene_id")["measured_share"]
        out["measured share"] = ids.map(h).to_numpy(float)
    return out


def _residual(x_tr, c_tr, x_all, c_all, smooth: bool):
    """x minus a fit on c, the fit made on the training genes only."""
    ok = np.isfinite(x_tr) & np.isfinite(c_tr)
    if ok.sum() < 20:
        return x_all
    if smooth:
        order = np.sort(c_tr[ok])
        rank = lambda c: np.searchsorted(order, c) / max(len(order), 1)
        # Polynomial.fit scales the input to [-1, 1] first, so a confounder with many ties
        # (most genes have no papers) does not make the cubic ill-conditioned.
        degree = int(min(3, len(np.unique(c_tr[ok])) - 1))   # 3 values hold at most a quadratic
        if degree < 1:
            return x_all
        fit = np.polynomial.Polynomial.fit(rank(c_tr[ok]), x_tr[ok], degree)(
            rank(np.nan_to_num(c_all, nan=np.nanmedian(c_tr))))
    else:
        coef = np.polyfit(c_tr[ok], x_tr[ok], 1)
        fit = np.polyval(coef, np.nan_to_num(c_all, nan=np.nanmedian(c_tr)))
    return x_all - fit


def cv_corrected(x, conf, y, kind, folds, method) -> float:
    """Held-out score of one feature after a correction fitted inside each training fold."""
    ok = np.isfinite(x) & pd.notna(y)
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
        xi = x[idx].copy()
        if method in ("linear", "smooth"):
            xi = _residual(xi[tr], conf[idx][tr], xi, conf[idx], smooth=(method == "smooth"))
        E = E1.one_column(np.where(np.isfinite(xi), xi, np.nan)[:, None])
        Efull = np.full((len(x), E.shape[1]), np.nan)
        Efull[idx] = E
        m_tr = np.zeros(len(x), bool); m_tr[idx[tr]] = True
        m_te = np.zeros(len(x), bool); m_te[idx[te]] = True
        m_tr &= ~np.isnan(Efull).any(axis=1)
        preds[te] = IM.knn_predict(Efull, y, kind, m_tr, m_te, classes)
    done = ~np.isnan(preds if preds.ndim == 1 else preds[:, 0])
    return IM.score(yv[done], preds[done], kind)


def stratified(x, conf, y, kind, folds) -> float:
    """Mean held-out score within confounder quintiles (genes compared only with similar genes)."""
    q = pd.qcut(pd.Series(conf).rank(method="first"), QUINTILES, labels=False).to_numpy()
    scores, weights = [], []
    for b in range(QUINTILES):
        keep = q == b
        yb = np.where(keep, y, None).astype(object) if kind != "continuous" else np.where(keep, y, np.nan)
        s = cv_corrected(x, conf, yb, kind, folds, "raw")
        if np.isfinite(s):
            scores.append(s); weights.append(keep.sum())
    return float(np.average(scores, weights=weights)) if scores else np.nan


def dependence(conf, y, kind, fam) -> float:
    """How far the confounder alone predicts the label above chance (AUROC - 0.5, or |rho|)."""
    s, _ = H.cv_score(E1.one_column(conf[:, None]), y, kind, H.family_folds(fam, H.N_FOLDS, 0))
    if not np.isfinite(s):
        return np.nan
    return abs(s) if kind == "continuous" else s - 0.5


def main(argv=None) -> int:
    """Score the corrections for every informative (label, feature) pair of E1."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--e1", required=True, help="the E1 output folder (information.parquet)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, default=SEEDS)
    ap.add_argument("--max-pairs", type=int, default=0, help="a quick check: the first N pairs only")
    ap.add_argument("--part", default="", help="i/n: score only pairs at positions i, i+n, ... (parallel parts)")
    ap.add_argument("--merge", action="store_true", help="join the part checkpoints and write the outputs")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    info = pd.read_parquet(os.path.join(args.e1, "information.parquet"))
    pairs = info[(info["level"] == "feature") & info["has_information"].fillna(False)][["label", "feature", "kind"]]
    if args.max_pairs:
        pairs = pairs.head(args.max_pairs)
    nodes, views, extra, L = H.load()
    L.update(E1.added_labels(nodes))
    fam = H.families(nodes)
    C = confounders(nodes)
    dep = {(lab, c): dependence(conf, L[lab][0], L[lab][1], fam)
           for lab in pairs["label"].unique() for c, conf in C.items()}
    seed_part = os.path.join(args.out, "confounders.partial.parquet")
    if args.merge:
        import glob as _g
        files = sorted(_g.glob(os.path.join(args.out, "*.partial.parquet")))
        t = pd.concat([pd.read_parquet(f) for f in files]).drop_duplicates(["label", "feature", "confounder", "method"])
        have = set(map(tuple, t[["label", "feature"]].drop_duplicates().to_numpy()))
        want = set(map(tuple, pairs[["label", "feature"]].to_numpy()))
        if want - have:
            raise SystemExit(f"merge: {len(want - have)} of {len(want)} pairs missing; parts not finished")
        rows = t.to_dict("records")
        args.part = ""
    part = seed_part
    if args.part:
        i_, n_p = (int(x) for x in args.part.split("/"))
        pairs = pairs.iloc[i_::n_p]
        part = os.path.join(args.out, f"part_{i_}_of_{n_p}.partial.parquet")
    if not args.merge:
        rows = pd.read_parquet(part).to_dict("records") if os.path.exists(part) else []
    done = {(d["label"], d["feature"]) for d in rows}
    if os.path.exists(seed_part) and part != seed_part:
        done |= set(map(tuple, pd.read_parquet(seed_part)[["label", "feature"]].drop_duplicates().to_numpy()))
    if done:
        log(f"resuming: {len(done)} pairs already scored in {part}")
    for n_, r in enumerate(pairs.itertuples(index=False)):
        if (r.label, r.feature) in done:
            continue
        y, kind = L[r.label][0], L[r.label][1]
        x = nodes[r.feature].to_numpy(float)
        for cname, conf in C.items():
            for method in ("raw", "linear", "smooth", "stratified"):
                s = [stratified(x, conf, y, kind, H.family_folds(fam, H.N_FOLDS, k)) if method == "stratified"
                     else cv_corrected(x, conf, y, kind, H.family_folds(fam, H.N_FOLDS, k), method)
                     for k in range(args.seeds)]
                rows.append({"label": r.label, "feature": r.feature, "kind": kind, "confounder": cname,
                             "method": method, "score": np.nanmean(s), "score_sd": np.nanstd(s, ddof=1),
                             "label_depends": dep[(r.label, cname)]})
        if n_ % 25 == 0:
            pd.DataFrame(rows).to_parquet(part, index=False)
            log(f"{n_ + 1} of {len(pairs)} pairs (checkpoint saved)")
    pd.DataFrame(rows).to_parquet(part, index=False) if rows and not args.merge else None
    if args.part:
        log(f"part {args.part} done: {len(rows)} rows in {part}")
        return 0
    t = pd.DataFrame(rows)
    raw = t[t["method"] == "raw"].set_index(["label", "feature", "confounder"])["score"]
    t["change_from_raw"] = t["score"] - t.set_index(["label", "feature", "confounder"]).index.map(raw)
    t["label_depends_on_confounder"] = t["label_depends"] > DEPENDS
    t.to_parquet(os.path.join(args.out, "confounders.parquet"), index=False)
    corr = t[t["method"] != "raw"]
    summary = {"pairs": int(len(pairs)), "seeds": args.seeds,
               "median_change_by_confounder_and_method": corr.groupby(["confounder", "method"])["change_from_raw"]
               .median().round(4).unstack().to_dict(),
               "labels_depending_on_a_confounder": sorted({f"{a} ~ {b}" for (a, b), v in dep.items()
                                                          if np.isfinite(v) and v > DEPENDS})}
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out} (part checkpoints kept)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
