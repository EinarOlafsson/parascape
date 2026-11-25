"""Figure 12, one claim: the most valuable new dataset for the field is a genome-wide host-response
screen. Every kind of knowledge that covers the genome today got there through one genome-wide
dataset; host response has none, and it cannot be predicted from what is measured.

Einar, 2026-11-19: placeholders unhelpful. All panels real, from the dataset registry (each
dataset's own columns, its citation year, its categories of knowledge) and the gene table.

  A  per kind of knowledge: largest single dataset vs all datasets combined     (one screen does it)
  B  genes measured over time, by kind (cumulative, by citation year)           (fitness: one jump)
  C  genes per dataset, by kind (median, IQR)                                   (host datasets small)
  D  host-response coverage in gene groups                                      (who is missing)
  E  the least-measured tenth of genes against all genes                        (what they are)
  F  unmeasured genes vs what prediction fills at a validated rate              (prediction no fix)
"""
from __future__ import annotations

import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from parascape import organisms as _ORG              # noqa: E402
TG, PF = _ORG.TOXOPLASMA, _ORG.FALCIPARUM
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import datasets as DS                    # noqa: E402
from parascape import figstyle as FS                    # noqa: E402
from parascape import information_categories as IC      # noqa: E402
import build_figure01_v2 as F1                          # noqa: E402

OUT = os.path.join(ROOT, "results", "figure12_v2")
PNL = os.path.join(OUT, "panels")
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
FM = os.path.join(ROOT, "results", "filling_the_map_2026_11_18", "best_per_target.tsv")
NAMES = F1.CATEGORY_NAMES
SHOWN = ["fitness", "transcript", "place", "protein_abundance", "interaction", "phenotype",
         "host_response", "function"]
SECRETED = F1.SECRETED


def save(fig, name):
    os.makedirs(PNL, exist_ok=True)
    p = os.path.join(PNL, f"{name}.png")
    fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return p


def rows_fig(nr, left=1.35, value=1.6, right=1.6, rh=0.2):
    w, h = left + value + right, 0.15 + nr * rh + 0.6
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([left / w, 0.6 / h, value / w, nr * rh / h])
    ax.set_ylim(nr - 0.5, -0.5); ax.tick_params(axis="y", length=0)
    return fig, ax


def per_dataset(nodes, knowledge):
    """Each Toxoplasma dataset: its year, its categories, and the genes it measured."""
    rows = []
    kcols = {c for v in knowledge.values() for c in v}
    for d in DS.REGISTRY:
        cols = [c for c in d.columns if c in nodes.columns and c in kcols]
        if not cols:
            continue
        cats = set()
        for c in cols:
            cats |= {x for x in IC.categories_of(c) if x in knowledge}
        y = re.findall(r"\b(19[5-9]\d|20[0-2]\d)\b", str(d.citation or ""))
        genes = nodes[cols].notna().any(axis=1).to_numpy()
        rows.append({"key": d.key, "year": int(y[0]) if y else np.nan, "cats": cats,
                     "genes": genes, "n": int(genes.sum())})
    return rows


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    knowledge = F1.knowledge_columns(nodes)
    ds = per_dataset(nodes, knowledge)
    G = len(nodes)
    cover = {c: nodes[v].notna().any(axis=1).to_numpy() for c, v in knowledge.items()}
    shown = [c for c in SHOWN if c in cover]

    # A
    rows = []
    for c in shown:
        mine = [d for d in ds if c in d["cats"]]
        rows.append((c, max((d["n"] for d in mine), default=0) / G * 100,
                     cover[c].mean() * 100, len(mine)))
    a = pd.DataFrame(rows, columns=["cat", "largest", "all", "n"]).sort_values("all",
                                                                             ascending=False)
    FS.use()
    na = [F1.kind_name(c) for c in a["cat"]]
    figa, ax = FS.ref_bars(a["largest"].to_numpy(), na, refs=a["all"].to_numpy(),
                           highlight=[F1.kind_name(c) for c in ("host_response", "phenotype")],
                           ref_label="all datasets combined",
                           ylabel="Genes in largest dataset (%)")
    ax.set_ylim(0, 105)
    print(a.round(1).to_string(index=False))

    # B: genes measured by the datasets published up to 2015, 2020 and now, per kind
    styles = {"fitness": None, "place": None, "phenotype": None, "host_response": None}
    years = np.arange(2005, 2027)
    cuts = [2015, 2020, 2026]
    series = []
    for yr in cuts:
        vals = []
        for c in styles:
            got = np.zeros(G, bool)
            for d in ds:
                if c in d["cats"] and not np.isnan(d["year"]) and d["year"] <= yr:
                    got |= d["genes"]
            vals.append(got.mean() * 100)
        series.append(vals)
    FS.use()
    figb, ax = FS.grouped_bars(series, [f"by {y}" for y in cuts],
                               [F1.kind_name(c) for c in styles],
                               colors=[FS.BACKGROUND, FS.CONTROL, FS.ACCENT],
                               ylabel="Genes measured (%)")
    ax.set_ylim(0, 100)

    # C: genes per dataset by kind (median bar, one point per dataset)
    FS.use()
    figc, ax = FS.dist_bars([[d_["n"] for d_ in ds if c in d_["cats"]] for c in a["cat"]], na,
                            highlight=[F1.kind_name(c) for c in ("host_response", "phenotype")],
                            ylabel="Genes per dataset", log=True)
    ax.set_yticks([10, 100, 1000, 10000]); ax.set_yticklabels(["10", "100", "1,000", "10,000"])
    ax.minorticks_off()

    # D: host-response coverage by gene group
    hr = cover["host_response"]
    fit = nodes["fit_invitro_hff"]
    prod = nodes["product"].astype(str).str.lower()
    groups = {"all genes": np.ones(G, bool),
              "secreted (hyperLOPIT)": nodes["lopit_map"].astype(str).isin(SECRETED).to_numpy(),
              "essential (lowest tenth)": (fit <= fit.quantile(0.1)).to_numpy(),
              "hypothetical": prod.str.contains("hypothetical").to_numpy(),
              "no localisation": nodes["lopit_map"].isna().to_numpy()}
    dvals = {k: hr[m].mean() * 100 for k, m in groups.items()}
    GS = {"all genes": "all genes", "secreted (hyperLOPIT)": "secreted",
          "essential (lowest tenth)": "essential", "hypothetical": "hypothetical",
          "no localisation": "no location"}
    FS.use()
    figd, ax = FS.ref_bars(list(dvals.values()), [GS[k] for k in dvals],
                           highlight=["secreted"], ylabel="Host response measured (%)")
    ax.set_ylim(0, 100)
    print({k: round(v, 1) for k, v in dvals.items()})

    # E: least-measured tenth vs all genes
    n_cat = np.column_stack([cover[c] for c in cover]).sum(axis=1)
    dark = n_cat <= np.quantile(n_cat, 0.1)
    props = {"hypothetical": groups["hypothetical"], "secreted (hyperLOPIT)":
             groups["secreted (hyperLOPIT)"], "essential (lowest tenth)":
             groups["essential (lowest tenth)"]}
    # "no localisation" is left out here: lacking a localisation is one of the counts that
    # DEFINES the least-measured tenth, so its 99.5 % would be circular.
    allv = [m.mean() * 100 for m in props.values()]
    darkv = [m[dark].mean() * 100 for m in props.values()]
    FS.use()
    fige, ax = FS.grouped_bars([allv, darkv], ["all genes", "least-measured tenth"],
                               [GS[k] for k in props], ylabel="Genes (%)")
    ax.set_ylim(0, 100)
    print(f"least-measured tenth: {int(dark.sum())} genes, <= {int(np.quantile(n_cat, .1))} "
          f"kinds; {dict(zip(props, np.round(darkv, 1)))}")

    # F: unmeasured vs filled by prediction at a validated rate
    fm = pd.read_csv(FM, sep="\t")
    fm = fm[fm.organism == TG]
    pick = {"place": "lopit_unified", "phenotype": "screen_any_phenotype"}
    rows = []
    for c in ["place", "phenotype", "host_response", "function"]:
        missing = int((~cover[c]).sum())
        filled = 0.0
        if c in pick:
            r = fm[fm.target == pick[c]].iloc[0]
            if r.lift_over_baseline > 0.05:
                col = nodes[pick[c]]
                miss_t = int(col.isna().sum())
                filled = miss_t * r.answer_rate * r.right_when_answered
        rows.append((c, missing, filled))
    f = pd.DataFrame(rows, columns=["cat", "missing", "filled"])
    FS.use()
    figf, ax = FS.grouped_bars([f["missing"].to_numpy(), f["filled"].where(f["filled"] > 0)
                                .to_numpy()],
                               ["never measured", "filled correctly by prediction"],
                               [F1.kind_name(c) for c in f["cat"]], ylabel="Genes", log=True)
    ax.set_yticks([1, 10, 100, 1000, 10000]); ax.minorticks_off()
    ax.set_yticklabels(["1", "10", "100", "1,000", "10,000"]); ax.set_ylim(0.7, 2e4)
    print(f.round(0).to_string(index=False))

    # ---- Figure 8 (storyboard 2026-11-20): what is missing ----
    import paper_panels as PP
    F = "figure_8__final"
    NODES = os.path.join(CACHE, "nodes.parquet")
    REG = os.path.join(ROOT, "parascape", "datasets.py")
    kinds = ["fitness", "place", "phenotype", "host_response"]
    kcol = {"fitness": FS.DARK_GREY, "place": FS.TEAL, "phenotype": FS.CONTROL,
            "host_response": FS.ACCENT}

    # A: saturation, genes covered as datasets are added largest first
    FS.use()
    figA, ax = FS.square_figure()
    sat = []
    for c in kinds:
        mine = sorted([d_ for d_ in ds if c in d_["cats"]], key=lambda d_: -d_["n"])
        got = np.zeros(G, bool); ys = [0.0]
        for d_ in mine:
            got |= d_["genes"]; ys.append(got.mean() * 100)
        ax.plot(range(len(ys)), ys, "-o", ms=2.5, lw=1.2, color=kcol[c], label=F1.kind_name(c))
        sat += [{"data_type": F1.kind_name(c), "datasets": k, "genes_pct": v}
                for k, v in enumerate(ys)]
    ax.set_xlabel("Datasets, largest first"); ax.set_ylabel("Genes measured (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    PP.panel(F, "a", fig=figA, data=pd.DataFrame(sat), stats={"genes": G}, sources=[REG, NODES],
             what="genes covered as datasets are added, largest first")

    # B: coverage by year, lines
    FS.use()
    figB, ax = FS.square_figure()
    cy = []
    for c in kinds:
        mine = [d_ for d_ in ds if c in d_["cats"] and not np.isnan(d_["year"])]
        ys = []
        for yr in years:
            got = np.zeros(G, bool)
            for d_ in mine:
                if d_["year"] <= yr:
                    got |= d_["genes"]
            ys.append(got.mean() * 100)
        ax.plot(years, ys, lw=1.2, color=kcol[c], label=F1.kind_name(c))
        cy += [{"data_type": F1.kind_name(c), "year": int(y_), "genes_pct": v}
               for y_, v in zip(years, ys)]
    ax.set_xlabel("Year of publication"); ax.set_ylabel("Genes measured (%)"); ax.set_ylim(0, 100)
    ax.set_xticks([2005, 2015, 2025])
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    PP.panel(F, "b", fig=figB, data=pd.DataFrame(cy), stats={"datasets with a year":
             int(sum(not np.isnan(d_["year"]) for d_ in ds))}, sources=[REG, NODES],
             what="genes measured over time, per data type")

    # C: what the least-measured tenth is, 100% stacked against all genes
    prodl = nodes["product"].astype(str).str.lower()
    ann = np.where(prodl.str.contains("hypothetical"), "hypothetical protein",
                   np.where(nodes["lopit_map"].astype(str).isin(F1.SECRETED), "secreted",
                            "annotated"))
    comp_ = pd.DataFrame({"group": np.where(dark, "least-measured tenth", "all genes"),
                          "annotation": ann})
    tabs = {"all genes": pd.Series(ann).value_counts(normalize=True) * 100,
            "least-measured tenth": pd.Series(ann[dark]).value_counts(normalize=True) * 100}
    cats_ = ["hypothetical protein", "annotated", "secreted"]
    FS.use()
    figC, ax = FS.plot_figure(2 * FS.SLOT_IN * 25.4 * 1.6, FS.PLOT_H_MM)
    bot = np.zeros(2)
    for cc, colr in zip(cats_, (FS.ACCENT, FS.CONTROL, FS.TEAL)):
        v = np.array([tabs[g].get(cc, 0.0) for g in tabs])
        ax.bar([0, 1], v, bottom=bot, width=0.6, color=colr, alpha=0.85, lw=0, label=cc)
        bot += v
    ax.set_xticks([0, 1]); ax.set_xticklabels(["all genes", "least known"], rotation=45,
                                              ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0); ax.set_ylim(0, 100); ax.set_ylabel("Genes (%)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    PP.panel(F, "c", fig=figC, data=pd.DataFrame(tabs).fillna(0).rename_axis("annotation"),
             stats={"least-measured tenth, genes": int(dark.sum()),
                    "data types at most": int(np.quantile(n_cat, .1))},
             sources=[(NODES, ["product", "lopit_map"])],
             what="annotation of all genes and of the least-measured tenth")

    # D: per data type, measured / filled by prediction / unknown
    fill = dict(zip(f["cat"], f["filled"]))
    rows_d = []
    for c in ["place", "function", "phenotype", "host_response"]:
        meas = cover[c].mean() * 100
        pr = fill.get(c, 0.0) / G * 100
        rows_d.append({"data_type": F1.kind_name(c), "measured": meas, "predicted": pr,
                       "unknown": 100 - meas - pr})
    dd_ = pd.DataFrame(rows_d)
    FS.use()
    figD, ax = FS.plot_figure(len(dd_) * FS.SLOT_IN * 25.4 * 1.2, FS.PLOT_H_MM)
    bot = np.zeros(len(dd_)); x_ = np.arange(len(dd_))
    for col_, c_ in (("measured", FS.DARK_GREY), ("predicted", FS.ACCENT),
                     ("unknown", "#E6E6E6")):
        ax.bar(x_, dd_[col_], bottom=bot, width=0.6, color=c_, lw=0, label=col_)
        bot += dd_[col_].to_numpy()
    ax.set_xticks(x_); ax.set_xticklabels(dd_.data_type, rotation=45, ha="right",
                                          rotation_mode="anchor")
    ax.tick_params(axis="x", length=0); ax.set_ylim(0, 100); ax.set_ylabel("Genes (%)")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3)
    PP.panel(F, "d", fig=figD, data=dd_, stats={"prediction": "best strategy where it beats "
             "the majority-class baseline by > 5 points (Fig. 4)"}, sources=[FM, NODES],
             what="measured, predicted and unknown per data type")
    # E and F: what papers show by experiment (the literature reading, instruction 102 / Affinage;
    # forms chosen by Einar 2026-11-25 from the candidate sheet: 3 stacked counts, 5 stacked area).
    sys.path.insert(0, os.path.join(ROOT, "scripts", "infer_bench"))
    import affinage_panel_data as AP
    from parascape import literature_findings as LF
    ev = AP.evidence(nodes, LF.labels(log=lambda *a: None).reset_index())
    rows_e = []
    for named, label in ((True, "named"), (False, "hypothetical protein")):
        g_ = ev[ev["named"] == named]
        rows_e.append({"product": label, "genes": len(g_),
                       "function observed": int(g_["observed_function"].sum()),
                       "none, needed in fibroblasts": int((~g_["observed_function"] & g_["needed_in_fibroblasts"]).sum()),
                       "none, not needed": int((~g_["observed_function"] & ~g_["needed_in_fibroblasts"]).sum())})
    de_ = pd.DataFrame(rows_e)
    FS.use()
    figE, ax = FS.plot_figure(2 * FS.SLOT_IN * 25.4 * 1.6, FS.PLOT_H_MM)
    bot = np.zeros(2)
    for col_, c_ in (("function observed", FS.ACCENT), ("none, needed in fibroblasts", FS.AMBER),
                     ("none, not needed", FS.CONTROL)):
        ax.bar([0, 1], de_[col_], bottom=bot, width=0.6, color=c_, alpha=FS.ALPHA, lw=0, label=col_)
        bot += de_[col_].to_numpy()
    ax.set_xticks([0, 1]); ax.set_xticklabels(["named", "hypoth."], rotation=45, ha="right",
                                              rotation_mode="anchor")
    ax.tick_params(axis="x", length=0); ax.set_ylabel("Genes")
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02))
    PP.panel(F, "e", fig=figE, data=de_,
             stats={"needed in fibroblasts": f"fibroblast phenotype score < {AP.FITNESS_CUT}",
                    "observed function": "a paper finding (activity, phenotype, host effect, complex), "
                                         "an EC number or an imaging phenotype"},
             sources=[(NODES, ["product", "ec_number", "fit_invitro_hff", "screen_any_phenotype"]),
                      os.path.join(ROOT, "parascape", "data", "literature_affinage",
                                   "literature_gene_labels.parquet")],
             what="genes with an observed function, named against hypothetical protein")

    yr = AP.by_year(LF.findings())
    kinds_f = ["localisation", "phenotype", "activity", "host_facing", "partner", "complex", "structure"]
    names_f = {"host_facing": "host effect"}
    FS.use()
    figF, ax = FS.square_figure()
    ax.stackplot(yr["year"], [yr.get(k, 0) for k in kinds_f],
                 colors=[FS.ACCENT, FS.TEAL, FS.NAVY, FS.MAGENTA, FS.ACCENT_LIGHT, FS.AMBER, FS.CONTROL],
                 alpha=FS.ALPHA, labels=[names_f.get(k, k) for k in kinds_f], lw=0)
    ax.set_xlabel("Year of publication"); ax.set_ylabel("Findings"); ax.set_ylim(0, None)
    ax.set_xticks([2000, 2010, 2020])
    FS.legend(ax, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    PP.panel(F, "f", fig=figF, data=yr,
             stats={"findings": int(yr["all"].sum()), "last year": "partial (to 2026-11-24)"},
             sources=[os.path.join(ROOT, "parascape", "data", "literature_affinage",
                                   "literature_findings.parquet")],
             what="findings shown by experiment, per year of publication and kind")
    PP.compose(F, "main", "Function, mutant phenotype and host response are missing for most "
               "Toxoplasma genes, and prediction does not fill them.")


if __name__ == "__main__":
    main()
