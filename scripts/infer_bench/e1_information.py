"""E1 of the Infer benchmark: the feature information table (raw, before any confounder correction).

Every feature (node column of a view) and every view, against every catalogue label, scored on
held-out genes with the harness of E0 (`e0_harness.py`):

- family-held-out 5-fold cross-validation, `SEEDS` seeds; the shuffled-label baseline with the same
  folds and seeds (labels permuted among the labelled genes, so class sizes are kept);
- score: AUROC of the neighbour vote (binary and categorical, one against the rest, averaged) or
  Spearman rho (continuous), kNN k = 15 on cosine, exactly as figure 6;
- **has information** (decided in E0, 2026-11-25): effect = mean real - mean shuffled; effect greater
  than 2 x max(real SD, shuffled SD) AND Benjamini-Hochberg q < 0.05 across all features of the
  label, from a one-sided z = effect / sqrt(real SD^2/n + shuffled SD^2/n).

Features WITHOUT information are kept, so a later release can show a feature gaining or losing it.

The catalogue is the figure 6 labels plus, from the literature reading (`literature_affinage`):
- **compartment (papers)**: a gene's single experimentally shown compartment (genes with a tie are
  left out); classes under 10 genes are not scored;
- **host-facing (papers)**: shown host-facing genes against genes with paper findings that are not
  host-facing (observed negatives, not "everything else").
No view column is derived from paper counts or co-mention, so neither literature label has a
circular feature (checked 2026-11-25); the exclusion is still declared through the provenance.

Writes `information.parquet` (one row per release x label x feature) and `summary.json`.

    python scripts/infer_bench/e1_information.py --out DIR [--labels "a;b"] [--seeds 3]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402

IM = H.IM
METHOD_VERSION = "e1-v1 (kNN k=15 cosine, family 5-fold, raw)"
MIN_CLASS = 10
SEEDS = 3


def literature_labels(nodes: pd.DataFrame) -> dict:
    """The two literature labels, aligned to the node table, in the figure 6 label format."""
    from parascape import literature_findings as LF
    t = LF.labels(log=lambda *a: None).reindex(nodes["gene_id"].astype(str))
    comp = t["lit_compartment"].fillna("").astype(str)
    comp = comp.where(comp.ne("") & ~comp.str.contains(r"\|"))
    counts = comp.value_counts()
    comp = comp.where(comp.map(counts).fillna(0) >= MIN_CLASS)
    studied = t["lit_n_papers"].fillna(0) > 0
    hf = np.where(t["lit_host_facing"].fillna(False).astype(bool), "host-facing",
                  np.where(studied, "not host-facing", None)).astype(object)
    return {"compartment (papers)": (comp.to_numpy(object), "categorical", ["__place__"]),
            "host-facing (papers)": (hf, "binary", [])}


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    v = p[ok]
    order = np.argsort(v)
    ranked = v[order] * len(v) / (np.arange(len(v)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(v)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def one_column(X: np.ndarray) -> np.ndarray:
    """A single column embedded for cosine kNN: standardised value plus a constant axis."""
    cover = ~np.isnan(X[:, 0])
    E = IM.embed(X, cover)
    E = np.c_[E, np.ones(len(E))]
    E[~cover] = np.nan
    return E


def score_row(E, y, kind, fam, seeds) -> dict:
    """Real and shuffled scores over seeds, family folds."""
    runs = [H.cv_score(E, y, kind, H.family_folds(fam, H.N_FOLDS, s)) for s in range(seeds)]
    real, genes = [r[0] for r in runs], runs[0][1]
    null = [H.cv_score(E, H.shuffled(y, 100 + s), kind, H.family_folds(fam, H.N_FOLDS, s))[0]
            for s in range(seeds)]
    rm, rs = np.nanmean(real), np.nanstd(real, ddof=1)
    nm, ns = np.nanmean(null), np.nanstd(null, ddof=1)
    eff = rm - nm
    se = np.sqrt(rs ** 2 / seeds + ns ** 2 / seeds)
    z = eff / se if se > 0 else np.nan
    return {"genes": genes, "score": rm, "score_sd": rs, "score_ci_low": rm - 1.96 * rs / np.sqrt(seeds),
            "score_ci_high": rm + 1.96 * rs / np.sqrt(seeds), "null_mean": nm, "null_sd": ns,
            "effect": eff, "z": z, "p": float(norm.sf(z)) if np.isfinite(z) else np.nan,
            "noise_floor": 2 * np.nanmax([rs, ns])}


def main(argv=None) -> int:
    """Score every feature and view against every label and write the table."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--labels", default="")
    ap.add_argument("--seeds", type=int, default=SEEDS)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    log = lambda s: print(f"[{time.time() - t0:7.0f}s] {s}", flush=True)
    nodes, views, extra, L = H.load()
    L.update(literature_labels(nodes))
    wanted = [x for x in args.labels.split(";") if x]
    if wanted:
        missing = set(wanted) - set(L)
        if missing:
            raise SystemExit(f"unknown labels: {sorted(missing)}; known: {list(L)}")
        L = {k: v for k, v in L.items() if k in wanted}
    fam = H.families(nodes)
    release = os.path.basename(IM.CACHE.rstrip("/"))
    rows = []
    for lname, (y, kind, tcols) in L.items():
        drop = IM.exclusions(tcols, views)
        labelled = int(pd.notna(y).sum())
        for vname, cols in views.items():
            for c in cols:
                if c in drop:
                    rows.append({"label": lname, "feature": c, "view": vname, "level": "feature",
                                 "excluded": "shares a source with the label"})
                    continue
                X = nodes[[c]].to_numpy(float)
                if (~np.isnan(X[:, 0])).sum() < 100:
                    continue
                r = score_row(one_column(X), y, kind, fam, args.seeds)
                rows.append({"label": lname, "feature": c, "view": vname, "level": "feature",
                             "excluded": "", **r})
        for vname, E in H.matrices(nodes, views, extra, drop).items():
            r = score_row(E, y, kind, fam, args.seeds)
            rows.append({"label": lname, "feature": f"view:{vname}", "view": vname, "level": "view",
                         "excluded": "", **r})
        log(f"{lname}: {labelled} labelled genes")
    t = pd.DataFrame(rows)
    t["q"] = np.nan
    for (lname, level), part in t[t["excluded"] == ""].groupby(["label", "level"]):
        t.loc[part.index, "q"] = bh(part["p"].to_numpy())
    t["has_information"] = (t["q"] < 0.05) & (t["effect"] > t["noise_floor"])
    t["kind"] = t["label"].map({k: v[1] for k, v in L.items()})
    t["release"], t["correction"], t["method_version"] = release, "raw", METHOD_VERSION
    t["seeds"], t["folds"] = args.seeds, f"family {H.N_FOLDS}-fold"
    t.to_parquet(os.path.join(args.out, "information.parquet"), index=False)
    f = t[(t["level"] == "feature") & (t["excluded"] == "")]
    v = t[t["level"] == "view"]
    summary = {"release": release, "method_version": METHOD_VERSION, "labels": list(L),
               "features_scored": int(f["feature"].nunique()), "rows": int(len(t)),
               "features_with_information_per_label": f.groupby("label")["has_information"].sum().astype(int).to_dict(),
               "views_with_information_per_label": v.groupby("label")["has_information"].sum().astype(int).to_dict(),
               "best_view_per_label": v.sort_values("effect").groupby("label").tail(1).set_index("label")["view"].to_dict()}
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1, default=str)
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
