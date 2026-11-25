"""E7 of the Infer benchmark: the chosen pipeline per label, calibrated, checked on untouched genes.

Design: parascape-paper `infer_benchmark/DESIGN.md`, sections 7-8. For each categorical or binary
catalogue label:

1. **final split**: 20% of gene FAMILIES (seed 2026) are put aside before anything else; no
   earlier stage chose a setting on them;
2. **calibration** on the other 80%: family 5-fold neighbour votes (all views, the preparation,
   similarity and reduction passed in; figure 6 by default) give each held-out gene the share of
   its neighbours voting for the top class; an isotonic fit maps that share to how often the top
   class was right;
3. **the untouched check**: trained on the 80%, the final 20% are called and graded: accuracy,
   how often right per calibrated band, and the expected calibration error, each against CHANCE
   (the share of the commonest class);
4. **candidates**: every unlabelled gene, called from all labelled genes, with its calibrated how
   often right; kept when that beats chance.

Writes `final_check.tsv`, `calibration.tsv` (per label and band) and `candidates.parquet`.

    python scripts/infer_bench/e7_inference.py --out DIR [--labels "a;b"] [--metric cosine]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e1_information as E1  # noqa: E402
import e4_similarity as E4  # noqa: E402
import e6_importance as E6  # noqa: E402
import e3_preparation as E3  # noqa: E402

IM = H.IM
FINAL_SEED = 2026
FINAL_SHARE = 0.20
BANDS = (0.0, 0.5, 0.7, 0.9, 1.01)


def final_split(fam: np.ndarray) -> np.ndarray:
    """True for genes whose family is in the untouched 20%."""
    uniq = pd.unique(fam)
    rng = np.random.default_rng(FINAL_SEED)
    chosen = set(rng.choice(uniq, int(round(FINAL_SHARE * len(uniq))), replace=False))
    return np.array([f in chosen for f in fam])


def votes(E, y, train, test, metric):
    """(top class, share of neighbours voting for it) for each test gene."""
    classes = np.unique(y[train])
    p = E4.knn_predict(E, y, "categorical", train, test, classes, metric)
    return classes[p.argmax(axis=1)], p.max(axis=1)


def ece(conf: np.ndarray, right: np.ndarray, bins: int = 10) -> float:
    """Expected calibration error: |how often right - stated| averaged over equal-width bins."""
    edges = np.linspace(0, 1, bins + 1)
    out = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        m = (conf >= a) & (conf < b) if b < 1 else (conf >= a)
        if m.any():
            out += m.mean() * abs(right[m].mean() - conf[m].mean())
    return float(out)


def main(argv=None) -> int:
    """Calibrate, check on the untouched genes, and write candidates."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--labels", default="")
    ap.add_argument("--metric", default="cosine", choices=E4.METRICS)
    ap.add_argument("--e3", default="", help="E3 output folder: view weighting chosen per label")
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
    final = final_split(fam)
    weighting = {}
    if args.e3:
        w = pd.read_parquet(os.path.join(args.e3, "preparation.parquet"))
        w = w[w["weighting"] != ""].pivot_table(index="label", columns="weighting", values="score")
        weighting = {lab: ("per_view" if r["per_view"] > r["per_column"] else "per_column")
                     for lab, r in w.iterrows()}
    ids = nodes["gene_id"].astype(str).to_numpy()
    checks, bands, cands = [], [], []
    for lname, (y, kind, tcols) in L.items():
        raw = E3.views_raw(nodes, views, extra, IM.exclusions(tcols, views))
        prepped = {v: E3.prepare(X, c) for v, (X, c) in raw.items()}
        if weighting.get(lname) == "per_view":          # E3: each view's block to unit variance
            for v, B in prepped.items():
                total = np.nansum(np.nanvar(B, axis=0))
                if total > 0:
                    prepped[v] = B / np.sqrt(total)
        E = E6.all_views(prepped)
        usable = ~np.isnan(E).any(axis=1)
        lab = pd.notna(y) & usable
        dev, fin = lab & ~final, lab & final
        # 2. calibration on the 80%, family folds
        folds = H.family_folds(fam, H.N_FOLDS, 0)
        share, right = np.full(len(y), np.nan), np.zeros(len(y), bool)
        for k in np.unique(folds[dev]):
            tr, te = dev & (folds != k), dev & (folds == k)
            if tr.sum() <= IM.K or not te.any():
                continue
            top, s = votes(E, y, tr, te, args.metric)
            share[te], right[te] = s, top == y[te]
        ok = dev & np.isfinite(share)
        iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip").fit(share[ok], right[ok])
        # 3. the untouched check
        top, s = votes(E, y, dev, fin, args.metric)
        conf = iso.predict(s)
        hit = top == y[fin]
        counts = pd.Series(y[dev]).value_counts()
        chance = float(counts.iloc[0] / counts.sum())
        auroc = np.nan
        if kind == "binary":                      # unbalanced: rank of the rarer class, not accuracy
            from sklearn.metrics import roc_auc_score
            rare = counts.index[-1]
            classes = np.unique(y[dev])
            p = E4.knn_predict(E, y, "categorical", dev, fin, classes, args.metric)
            truth = (y[fin] == rare).astype(int)
            if 0 < truth.sum() < len(truth):
                auroc = float(roc_auc_score(truth, p[:, list(classes).index(rare)]))
        checks.append({"label": lname, "weighting": weighting.get(lname, "per_column"),
                       "auroc_rarer_class": auroc,
                       "development_genes": int(dev.sum()), "final_genes": int(fin.sum()),
                       "accuracy": float(hit.mean()), "chance": chance,
                       "lift_over_chance": float(hit.mean() - chance),
                       "stated_mean": float(conf.mean()), "ece": ece(conf, hit)})
        for a, b in zip(BANDS[:-1], BANDS[1:]):
            m = (conf >= a) & (conf < b)
            if m.any():
                bands.append({"label": lname, "band": f"{a:.1f}-{min(b, 1):.1f}", "genes": int(m.sum()),
                              "stated": float(conf[m].mean()), "right": float(hit[m].mean()), "chance": chance})
        # 4. candidates: unlabelled genes, called from every labelled gene
        unl = ~pd.notna(y) & usable
        if unl.any():
            top_u, s_u = votes(E, y, lab, unl, args.metric)
            c_u = iso.predict(s_u)
            keep = c_u > chance
            if kind == "binary":                  # "other" is no finding; only the rarer class is
                keep &= top_u != counts.index[0]
            cands.append(pd.DataFrame({"label": lname, "gene_id": ids[unl][keep], "call": top_u[keep],
                                       "how_often_right": c_u[keep], "chance": chance}))
        log(f"{lname}: accuracy {hit.mean():.3f} vs chance {chance:.3f}, ECE {checks[-1]['ece']:.3f}")
    pd.DataFrame(checks).to_csv(os.path.join(args.out, "final_check.tsv"), sep="\t", index=False)
    pd.DataFrame(bands).to_csv(os.path.join(args.out, "calibration.tsv"), sep="\t", index=False)
    C = pd.concat(cands, ignore_index=True) if cands else pd.DataFrame()
    C.to_parquet(os.path.join(args.out, "candidates.parquet"), index=False)
    json.dump({"metric": args.metric, "final_share": FINAL_SHARE, "final_seed": FINAL_SEED,
               "checks": checks, "candidates_per_label": C.groupby("label").size().to_dict() if len(C) else {}},
              open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
