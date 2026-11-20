"""Supplementary figure 3a: the design choices behind figure 3, each tested with data.

Figure 3 cuts genes into 24 co-expression modules (Ward linkage) and names a module after a
compartment only when it holds more of that compartment than random genes MATCHED on how well
studied they are. Each choice is varied here:

  A  how many modules (chosen 24): share of modules named, share of genes in named modules
  B  the baseline: matched random genes vs the genome share of the compartment, per module
  C  what the null is matched on (chosen: kinds of knowledge per gene) vs not matched at all
  D  module sizes, and the 20-gene minimum below which a module is not named
  E  linkage method (chosen Ward) vs average linkage
  F  stability of the naming across five random seeds of the null
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from parascape import figstyle as FS                    # noqa: E402
import build_coexpression_modules as CM                 # noqa: E402
import paper_panels as PP                               # noqa: E402

F = "supplementary_figure_3a"
CACHE = os.environ.get("PARASCAPE_CACHE",
                       "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots")
#: As figure 3. NOT fewer: with 300 draws the smallest p is 1/301 and the Bonferroni q over 23
#: modules can never fall below 0.05, so no module could be named in any variant (found
#: 2026-11-20 -- every row of this figure read 0 named).
CM.N_NULL = 1000


def summary(table, assignment):
    named = table[table["named"]]
    genes_named = assignment.isin(named["module"]).mean() * 100
    return {"modules": len(table), "named": int(len(named)),
            "named_pct": len(named) / max(1, len(table)) * 100, "genes_in_named_pct": genes_named}


def rows_fig(labels, values, chosen, xlabel, ref=None, ref_label=None):
    FS.use()
    nr = len(labels); rh = 0.2
    w, h = 1.3 + 1.6 + 1.6, 0.15 + nr * rh + 0.6
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([1.3 / w, 0.6 / h, 1.6 / w, nr * rh / h])
    y = np.arange(nr)
    for yy, lab, v in zip(y, labels, values):
        col = FS.ACCENT if lab == chosen else FS.CONTROL
        ax.hlines(yy, 0, v, color=col, lw=0.6, alpha=FS.ALPHA)
        ax.plot(v, yy, "o", color=col, ms=4, alpha=FS.ALPHA)
    if ref is not None:
        for yy, r in zip(y, ref):
            ax.plot(r, yy, "o", color=FS.DARK_GREY, ms=3, mfc="white")
        ax.plot([], [], "o", color=FS.DARK_GREY, mfc="white", label=ref_label)
    ax.set_yticks(y); ax.set_yticklabels(labels); ax.set_ylim(nr - 0.5, -0.5)
    ax.tick_params(axis="y", length=0); ax.set_xlabel(xlabel)
    ax.plot([], [], "o", color=FS.ACCENT, label="choice in figure 3")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    return fig


def main():
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    comps = pd.read_parquet(CM.COMPONENTS)
    nodes = pd.read_parquet(os.path.join(CACHE, "nodes.parquet"))
    nodes = nodes.set_index(nodes["gene_id"].astype(str))
    shared = comps.index.intersection(nodes.index)
    comps, nodes = comps.loc[shared], nodes.loc[shared]
    placed = comps.notna().all(axis=1)
    comps, nodes = comps[placed], nodes[placed]
    labels = nodes["compartment_best"].astype(str)
    measured = nodes.notna().sum(axis=1).astype(float)
    srcs = [CM.COMPONENTS, (os.path.join(CACHE, "nodes.parquet"), ["compartment_best"])]

    # A: number of modules
    ks, rA = (12, 18, 24, 30, 36), []
    for k in ks:
        asg = CM.modules(comps, n=k)
        rA.append({"k": k, **summary(CM.name_modules(asg, labels, measured), asg)})
        print("A", rA[-1])
    a = pd.DataFrame(rA)
    PP.panel(F, "a", fig=rows_fig([f"{k} modules" for k in ks], a["named_pct"], "24 modules",
                                  "Modules named (%)", ref=a["genes_in_named_pct"],
                                  ref_label="genes in named modules (%)"),
             data=a, stats=a, sources=srcs, what="how many modules")

    asg = CM.modules(comps, n=24)
    tab = CM.name_modules(asg, labels, measured)

    # B: matched baseline vs genome prior
    prior = labels.value_counts(normalize=True)
    b = tab.dropna(subset=["compartment"]).copy()
    b = b[b["compartment"] != ""]
    b["genome_share"] = b["compartment"].map(prior)
    FS.use()
    figb, ax = FS.square_figure()
    mx = float(max(b.genome_share.max(), b.matched_baseline.max())) * 100 * 1.1
    ax.plot([0, mx], [0, mx], color=FS.BACKGROUND, lw=0.6)
    ax.plot(b.genome_share * 100, b.matched_baseline * 100, "o", color=FS.ACCENT, ms=4,
            alpha=FS.ALPHA, lw=0, label="module")
    ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_xlabel("Genome share of the compartment (%)")
    ax.set_ylabel("Matched random genes (%)")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    above = int((b.matched_baseline > b.genome_share).sum())
    PP.panel(F, "b", fig=figb, data=b[["module", "compartment", "share", "matched_baseline",
                                       "genome_share"]],
             stats={"modules": len(b), "matched baseline above genome share": above,
                    "median ratio matched / genome": float((b.matched_baseline /
                                                            b.genome_share).median())},
             sources=srcs, what="matched baseline vs genome share of the compartment")

    # C: matched vs unmatched null
    flat = pd.Series(1.0, index=measured.index)          # every gene equally measured
    tu = CM.name_modules(asg, labels, flat)
    rC = [{"null": "matched on knowledge", **summary(tab, asg)},
          {"null": "not matched", **summary(tu, asg)}]
    c = pd.DataFrame(rC)
    PP.panel(F, "c", fig=rows_fig(list(c["null"]), c["named"], "matched on knowledge",
                                  "Modules named"),
             data=c, stats=c, sources=srcs, what="what the null is matched on")

    # D: module sizes
    sizes = asg.value_counts().sort_values(ascending=False)
    FS.use()
    figd, ax = FS.square_figure()
    ax.hist(sizes.values, bins=20, color=FS.CONTROL, alpha=FS.ALPHA, lw=0, label="modules")
    ax.axvline(CM.MIN_MODULE, color=FS.ACCENT, lw=1.2, label=f"minimum ({CM.MIN_MODULE} genes)")
    ax.set_xlabel("Genes per module"); ax.set_ylabel("Modules")
    FS.legend(ax, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    PP.panel(F, "d", fig=figd, data=sizes.rename("genes").rename_axis("module").reset_index(),
             stats={"modules": len(sizes), f"below {CM.MIN_MODULE} genes":
                    int((sizes < CM.MIN_MODULE).sum()), "median size": float(sizes.median())},
             sources=srcs, what="module sizes and the minimum")

    # E: linkage method
    v = comps.to_numpy(float); v = (v - v.mean(0)) / (v.std(0) + 1e-12)
    avg = pd.Series(fcluster(linkage(pdist(v), method="average"), t=24, criterion="maxclust"),
                    index=comps.index, name="module")
    rE = [{"linkage": "Ward", **summary(tab, asg)},
          {"linkage": "average", **summary(CM.name_modules(avg, labels, measured), avg)}]
    e = pd.DataFrame(rE)
    PP.panel(F, "e", fig=rows_fig(list(e["linkage"]), e["genes_in_named_pct"], "Ward",
                                  "Genes in named modules (%)"),
             data=e, stats=e, sources=srcs, what="linkage method")

    # F: seed stability
    rF = [{"seed": s, **summary(CM.name_modules(asg, labels, measured, seed=s), asg)}
          for s in range(5)]
    f = pd.DataFrame(rF)
    PP.panel(F, "f", fig=rows_fig([f"seed {s}" for s in f.seed], f["named"], "seed 0",
                                  "Modules named"),
             data=f, stats={"named, min": int(f.named.min()), "named, max": int(f.named.max())},
             sources=srcs, what="stability across null seeds")
    print(c.to_string(index=False)); print(e.to_string(index=False)); print(f.to_string(index=False))
    PP.compose(F, "supplementary", "Every design choice behind figure 3, varied with the rest "
               "held fixed.")


if __name__ == "__main__":
    main()
