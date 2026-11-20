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
    fm = fm[fm.organism == "Tg"]
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

    import paper_panels as PP
    F = "figure_12"
    NODES = os.path.join(CACHE, "nodes.parquet")
    REG = os.path.join(ROOT, "parascape", "datasets.py")
    cum = []
    for c in styles:
        mine = [d_ for d_ in ds if c in d_["cats"] and not np.isnan(d_["year"])]
        for yr in years:
            got = np.zeros(G, bool)
            for d_ in mine:
                if d_["year"] <= yr:
                    got |= d_["genes"]
            cum.append({"kind": NAMES[c], "year": int(yr), "genes_measured_pct": got.mean() * 100})
    PP.panel(F, "a", fig=figa, data=a.rename(columns={"largest": "largest_single_pct",
             "all": "all_combined_pct", "n": "datasets"}), stats={"genes": G},
             sources=[REG, NODES], what="largest single dataset vs all combined")
    PP.panel(F, "b", fig=figb, data=pd.DataFrame(cum), stats={"datasets with a citation year":
             int(sum(not np.isnan(d_["year"]) for d_ in ds))}, sources=[REG, NODES],
             what="genes measured over time")
    PP.panel(F, "c", fig=figc, data=pd.DataFrame([{"dataset": d_["key"], "genes": d_["n"],
             "kinds": ";".join(sorted(d_["cats"]))} for d_ in ds]), stats={"datasets": len(ds)},
             sources=[REG, NODES], what="genes per dataset by kind")
    PP.panel(F, "d", fig=figd, data=pd.Series(dvals, name="host_response_measured_pct")
             .rename_axis("group").reset_index(), stats={k: int(m_.sum()) for k, m_ in
             groups.items()}, sources=[(NODES, ["lopit_map", "fit_invitro_hff", "product"])],
             what="host-response coverage by gene group")
    PP.panel(F, "e", fig=fige, data=pd.DataFrame({"group": list(props), "all_genes_pct": allv,
             "least_measured_tenth_pct": darkv}), stats={"least-measured tenth, genes":
             int(dark.sum()), "kinds of knowledge at most": int(np.quantile(n_cat, .1))},
             sources=[NODES], what="the least-measured tenth of genes")
    PP.panel(F, "f", fig=figf, data=f, stats={"prediction": "best strategy from figure 11, "
             "only where it beats guessing by > 5 points"}, sources=[FM, NODES],
             what="unmeasured genes vs what prediction fills")
    out = PP.compose(F, "main", "The most valuable new dataset for the field is a genome-wide host-response screen.", rows=[3, 3])
    print("figure 12 ->", out)


if __name__ == "__main__":
    main()
