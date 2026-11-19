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
    figa, ax = rows_fig(len(a))
    y = np.arange(len(a))
    ax.hlines(y, a["largest"], a["all"], color=FS.CONTROL, lw=0.6)
    ax.plot(a["all"], y, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA,
            label="all datasets combined (ref.)")
    ax.plot(a["largest"], y, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA,
            label="largest single dataset")
    ax.set_yticks(y); ax.set_yticklabels([NAMES[c] for c in a["cat"]])
    ax.set_xlim(-3, 103); ax.set_xlabel("Genes measured (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(a.round(1).to_string(index=False))

    # B: cumulative coverage by year
    years = np.arange(2005, 2027)
    FS.use()
    figb, ax = FS.square_figure()
    styles = {"fitness": ("-", FS.DARK_GREY), "place": ("--", FS.DARK_GREY),
              "phenotype": (":", FS.DARK_GREY), "host_response": ("-", FS.ACCENT)}
    for c, (ls, col) in styles.items():
        mine = [d for d in ds if c in d["cats"] and not np.isnan(d["year"])]
        cum = []
        for yr in years:
            got = np.zeros(G, bool)
            for d in mine:
                if d["year"] <= yr:
                    got |= d["genes"]
            cum.append(got.mean() * 100)
        ax.plot(years, cum, ls=ls, color=col, lw=1.2, label=NAMES[c])
    ax.set_xlabel("Year of publication"); ax.set_ylabel("Genes measured (%)")
    ax.set_ylim(0, 100)
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))

    # C: genes per dataset by kind
    FS.use()
    figc, ax = rows_fig(len(a), right=0.3)
    for yy, c in zip(y, a["cat"]):
        v = np.array([d["n"] for d in ds if c in d["cats"]], float)
        if not len(v):
            continue
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        col = FS.ACCENT if c == "host_response" else FS.CONTROL
        ax.hlines(yy, max(q1, 1), q3, color=col, lw=1.2, alpha=FS.ALPHA)
        ax.plot(max(med, 1), yy, "o", color=col, ms=3.5)
    ax.set_yticks(y); ax.set_yticklabels([NAMES[c] for c in a["cat"]])
    ax.set_xscale("log"); ax.set_xticks([10, 100, 1000, 10000]); ax.minorticks_off()
    ax.set_xticklabels(["10", "100", "1,000", "10,000"])
    ax.set_xlabel("Genes per dataset")

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
    FS.use()
    figd, ax = rows_fig(len(groups), right=0.3)
    yd = np.arange(len(groups))
    ax.hlines(yd, 0, list(dvals.values()), color=FS.ACCENT, lw=0.6, alpha=FS.ALPHA)
    ax.plot(list(dvals.values()), yd, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA)
    ax.set_yticks(yd); ax.set_yticklabels(list(dvals))
    ax.set_xlim(0, 100); ax.set_xlabel("Host response measured (%)")
    print({k: round(v, 1) for k, v in dvals.items()})

    # E: least-measured tenth vs all genes
    n_cat = np.column_stack([cover[c] for c in cover]).sum(axis=1)
    dark = n_cat <= np.quantile(n_cat, 0.1)
    props = {"hypothetical": groups["hypothetical"], "secreted (hyperLOPIT)":
             groups["secreted (hyperLOPIT)"], "essential (lowest tenth)":
             groups["essential (lowest tenth)"]}
    # "no localisation" is left out here: lacking a localisation is one of the counts that
    # DEFINES the least-measured tenth, so its 99.5 % would be circular.
    FS.use()
    fige, ax = rows_fig(len(props))
    ye = np.arange(len(props))
    allv = [m.mean() * 100 for m in props.values()]
    darkv = [m[dark].mean() * 100 for m in props.values()]
    ax.hlines(ye, allv, darkv, color=FS.CONTROL, lw=0.6)
    ax.plot(allv, ye, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA, label="all genes (ref.)")
    ax.plot(darkv, ye, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA,
            label="least-measured tenth")
    ax.set_yticks(ye); ax.set_yticklabels(list(props))
    ax.set_xlim(0, 100); ax.set_xlabel("Genes (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
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
    figf, ax = rows_fig(len(f))
    yf = np.arange(len(f))
    ax.hlines(yf, f["filled"].clip(lower=1), f["missing"], color=FS.CONTROL, lw=0.6)
    ax.plot(f["missing"], yf, "o", color=FS.CONTROL, ms=3.5, alpha=FS.ALPHA,
            label="genes never measured (ref.)")
    ax.plot(f["filled"].clip(lower=1), yf, "o", color=FS.ACCENT, ms=3.5, alpha=FS.ALPHA,
            label="filled correctly by prediction")
    ax.set_yticks(yf); ax.set_yticklabels([NAMES[c] for c in f["cat"]])
    ax.set_xscale("log"); ax.set_xticks([1, 10, 100, 1000, 10000]); ax.minorticks_off()
    ax.set_xticklabels(["0", "10", "100", "1,000", "10,000"])
    ax.set_xlabel("Genes")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    print(f.round(0).to_string(index=False))

    panels = [("a", save(figa, "a_largest_dataset")), ("b", save(figb, "b_over_time")),
              ("c", save(figc, "c_genes_per_dataset")), ("d", save(figd, "d_host_by_group")),
              ("e", save(fige, "e_least_measured")), ("f", save(figf, "f_prediction_gap"))]
    import assemble_paper_figures as A_
    out = A_.compose("12", [{"png": p, "letter": l} for l, p in panels], OUT, rows=[3, 3])
    print("figure 12 ->", out)


if __name__ == "__main__":
    main()
