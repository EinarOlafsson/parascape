"""Supplementary figures 4a, 10a and 12a: the design choices behind figures 4, 10 and 1/12.

4a  similarity of contrasts: (A) minimum shared genes, (B) signed vs absolute rho,
    (C) Spearman vs Pearson, (D) without the two studies with the most contrasts.
10a host response: (A) the cut for a "shared" host gene (5, 10, 15 effectors), (B) overlap within
    the translocon, translocon to cargo, and within cargo, (C) single-cell knockouts: cells per
    knockout against the evidence it gives (power).
12a coverage: (A) columns set aside because they hold a value for every gene, per kind,
    (B) coverage per kind with and without that rule, (C) datasets with and without a citation
    year, by level of biology.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import paper_panels as PP                               # noqa: E402

CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")


def pair_rows(labels, ref, val, xlabel, legend, log=False):
    """Bars (the value), with the reference as a tick when there is one."""
    FS.use()
    fig, ax = FS.ref_bars(list(val), list(labels), refs=None if ref is None else list(ref),
                          highlight=list(labels), ref_label=legend[0].replace(" (ref.)", ""),
                          ylabel=xlabel, log=log)
    return fig


def supp4():
    from parascape import perturbation_atlas as PA
    F = "supplementary_figure_4a"
    store = PA.load(); lfc, reg = store["lfc"], store["registry"]
    study = reg["study"].astype(str).to_dict()
    cols = list(lfc.columns)
    s = np.array([study[c] for c in cols])
    same = (s[:, None] == s[None, :])[np.triu_indices(len(cols), 1)]
    srcs = [os.path.join(ROOT, "results", "what_structures_contrasts_2026_11_18", "pairs.tsv")]

    def medians(r, absval=True):
        v = r.to_numpy()[np.triu_indices(len(cols), 1)]
        v = np.abs(v) if absval else v
        ok = ~np.isnan(v)
        return float(np.median(v[ok & same])), float(np.median(v[ok & ~same]))

    rows = []
    for mn in (500, 1000, 2000, 4000):
        r = lfc.corr(method="spearman", min_periods=mn)
        a, b = medians(r); rows.append({"variant": f">= {mn:,} genes", "same": a, "other": b})
        if mn == 1000:
            base = r
    a, b = medians(base, absval=False); rows.append({"variant": "signed rho", "same": a, "other": b})
    pr = lfc.corr(method="pearson", min_periods=1000)
    a, b = medians(pr); rows.append({"variant": "Pearson", "same": a, "other": b})
    big = pd.Series(s).value_counts().index[:2]
    keep = [c for c in cols if study[c] not in set(big)]
    s2 = np.array([study[c] for c in keep])
    same2 = (s2[:, None] == s2[None, :])[np.triu_indices(len(keep), 1)]
    v = np.abs(base.loc[keep, keep].to_numpy()[np.triu_indices(len(keep), 1)]); ok = ~np.isnan(v)
    rows.append({"variant": "no 2 largest", "same": float(np.median(v[ok & same2])),
                 "other": float(np.median(v[ok & ~same2]))})
    d = pd.DataFrame(rows); d["ratio"] = d.same / d.other
    PP.panel(F, "a", fig=pair_rows(list(d.variant), d.other, d.same, "Median similarity",
                                   ["other study (ref.)", "same study"]),
             data=d, stats=d, sources=srcs, what="same vs other study under each variant")
    print(d.round(3).to_string(index=False))
    PP.compose(F, "supplementary", "The laboratory effect of figure 4 holds under every variant "
               "of the similarity measure.")


def supp10():
    F = "supplementary_figure_10a"
    HR = os.path.join(ROOT, "results", "host_response_shared_2026_11_18")
    epg = pd.read_csv(os.path.join(HR, "effectors_per_host_gene.tsv"), sep="\t")
    summ = json.load(open(os.path.join(HR, "summary.json")))
    ht = pd.read_parquet(os.path.join(ROOT, "results", "hosttx_per_host_gene_2026_11_18",
                                      "hosttx_host_gene_effects.parquet"))
    srcs = [os.path.join(HR, f) for f in ("effectors_per_host_gene.tsv", "summary.json")]
    rows = []
    for cut in (5, 10, 15):
        sh = set(epg.loc[epg.effectors_moving_it >= cut, "host_id"])
        sub = ht[ht.host_id.isin(sh)]
        sign = sub.groupby("host_id")["hosttx_log2fc"].apply(lambda v: max((v > 0).mean(),
                                                                          (v < 0).mean()))
        rows.append({"cut": f">= {cut} effectors", "shared_genes": len(sh),
                     "same_direction_pct": float(sign.mean() * 100) if len(sign) else np.nan})
    a = pd.DataFrame(rows)
    PP.panel(F, "a", fig=pair_rows(list(a.cut), None, a.same_direction_pct,
                                   "Effectors moving it the same way (%)", ["shared host genes"]),
             data=a, stats=a, sources=srcs + [os.path.join(ROOT, "results",
                                                           "hosttx_per_host_gene_2026_11_18")],
             what="the cut for a shared host gene")
    b = pd.DataFrame({"pairs": ["within translocon", "translocon to cargo", "within cargo"],
                      "overlap": [summ["overlap_translocon_within"],
                                  summ["overlap_translocon_to_cargo"],
                                  summ["overlap_cargo_within"]]})
    FS.use()
    fb, ax = FS.figure(3, vertical=True)
    FS.bars(ax, b.overlap.to_numpy() * 100, list(b.pairs), highlight=["within translocon"])
    ax.set_ylabel("Host genes in common (%)"); ax.set_xlim(-0.6, 2.4)
    PP.panel(F, "b", fig=fb, data=b, stats=b, sources=srcs, what="overlap by effector class")
    DP = os.path.join(ROOT, "results", "dual_perturb_single_cell_2026_11_18")
    dp = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(DP, "*_targets.tsv"))])
    FS.use()
    fc, ax = FS.square_figure()
    ax.plot(dp.cells, -np.log10(dp.parasite_q.clip(lower=1e-12)), "o", color=FS.ACCENT, ms=3,
            alpha=FS.ALPHA, lw=0, label="knockout")
    ax.axhline(-np.log10(0.05), color=FS.INK, lw=0.6, ls="--", label="q = 0.05")
    ax.set_xlabel("Cells with the knockout"); ax.set_ylabel("Parasite change ($-$log$_{10}$ q)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "c", fig=fc, data=dp, stats={"median cells per knockout": float(dp.cells.median()),
             "knockouts": len(dp)}, sources=sorted(glob.glob(os.path.join(DP, "*"))),
             what="cells per knockout against the evidence it gives")
    print(a.to_string(index=False))
    PP.compose(F, "supplementary", "The design choices behind figure 10, varied or checked.")


def supp12():
    F = "supplementary_figure_12a"
    import build_figure01_v2 as F1
    from parascape import information_categories as IC
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    derived = set(getattr(IC, "DERIVED_FROM", {}))
    kept, all_cols, set_aside = {}, {}, {}
    for c in nodes.columns:
        for cat in IC.categories_of(c):
            if cat in F1.NOT_KNOWLEDGE:
                continue
            all_cols.setdefault(cat, []).append(c)
            if c in derived or not nodes[c].isna().any():
                set_aside[cat] = set_aside.get(cat, 0) + 1
            else:
                kept.setdefault(cat, []).append(c)
    rows = []
    for cat in all_cols:
        with_rule = nodes[kept[cat]].notna().any(axis=1).mean() * 100 if cat in kept else 0.0
        without = nodes[all_cols[cat]].notna().any(axis=1).mean() * 100
        rows.append({"kind": F1.CATEGORY_NAMES.get(cat, cat), "columns": len(all_cols[cat]),
                     "set_aside": set_aside.get(cat, 0), "coverage_with_rule_pct": with_rule,
                     "coverage_without_rule_pct": without})
    d = pd.DataFrame(rows).sort_values("coverage_with_rule_pct", ascending=False)
    srcs = [os.path.join(CACHE, "nodes.parquet"), os.path.join(ROOT, "parascape",
                                                                "information_categories.py")]
    d["short"] = [F1.SHORT_KIND.get(k, k) for k in d.kind]
    PP.panel(F, "a", fig=pair_rows(list(d.short), None, d.set_aside, "Columns set aside",
                                   ["value for every gene"]),
             data=d, stats={"columns set aside": int(d.set_aside.sum())}, sources=srcs,
             what="columns set aside because they hold a value for every gene")
    PP.panel(F, "b", fig=pair_rows(list(d.short), d.coverage_without_rule_pct,
                                   d.coverage_with_rule_pct, "Genes measured (%)",
                                   ["without the rule (ref.)", "with the rule (figure 1)"]),
             data=d, stats=d, sources=srcs, what="coverage with and without the rule")
    yrs = F1.dataset_years()
    from parascape import datasets as DS
    yrs["level"] = [d_.level for d_ in DS.REGISTRY]
    c = yrs.assign(dated=yrs.year.notna()).groupby("level")["dated"].agg(["sum", "count"])
    c["undated"] = c["count"] - c["sum"]
    PP.panel(F, "c", fig=pair_rows(list(c.index), c["count"], c["sum"], "Datasets",
                                   ["all (ref.)", "with a citation year"]),
             data=c.rename_axis("level"), stats={"dated": int(c["sum"].sum()),
                                                  "all": int(c["count"].sum())},
             sources=[os.path.join(ROOT, "parascape", "datasets.py")],
             what="datasets with a citation year, by level")
    print(d.round(1).to_string(index=False))
    PP.compose(F, "supplementary", "The coverage rule behind figures 1 and 12: a column with a "
               "value for every gene is not a measurement.")


if __name__ == "__main__":
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    which = sys.argv[1:] or ["4", "10", "12"]
    for w in which:
        {"4": supp4, "10": supp10, "12": supp12}[w]()
