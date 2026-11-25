"""Figure of E8: does clustering on the features E1 selects beat all views and random sets of the same size?

DEVELOPMENT FIGURE made by a script from an E8 output folder (`e8_selected_clustering.py`). Panels:

A. label x embedding (all views = reference, E1-selected features, random sets of the same size = mean of 3
   draws) at Louvain resolution 0.5 and 1.0: held-out recovery gain over shuffled development labels;
B. the same embeddings: kNN accuracy on untouched genes minus the majority-class rate (kNN lift).

    python scripts/infer_bench/fig_e8.py E8_DIR OUT_DIR
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parascape import figstyle as F  # noqa: E402
import fig_e9 as FE9  # noqa: E402

GROUPS = ("all views", "selected", "random")
GROUP_NAME = {"all views": "all (ref.)", "selected": "selected", "random": "random"}
SHORT = {"interaction module (STRING, measured in Toxoplasma)": "STRING module, measured",
         "interaction module (STRING, copied by orthology)": "STRING module, orthology",
         "complex (CORUM human, copied by orthology)": "CORUM complex, orthology",
         "enzymatic function (EC class)": "enzyme class (EC)"}
MM = 1 / 25.4


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src, out = argv
    os.makedirs(out, exist_ok=True)
    F.use()
    t = pd.read_parquet(os.path.join(src, "results.parquet"))
    t["group"] = t["embedding"].str.replace(r"random \d+", "random", regex=True)
    labels = list(pd.unique(t["label"]))
    A = t.groupby(["label", "resolution", "group"])["recovery_gain"].mean().unstack(["resolution", "group"])
    A = A.reindex(index=labels, columns=[(r, g) for r in sorted(t["resolution"].unique()) for g in GROUPS])
    res_list = sorted(t["resolution"].unique())
    A.columns = [f"{GROUP_NAME[g]} @{r:g}" for r, g in A.columns]
    Bk = t.groupby(["label", "group"])["knn_lift"].mean().unstack().reindex(index=labels, columns=list(GROUPS))
    Bk.columns = [GROUP_NAME[g] for g in Bk.columns]
    A.index = Bk.index = [SHORT.get(x, x) for x in labels]
    va = max(float(np.nanmax(np.abs(A.to_numpy(float)))), 0.05)
    vb = max(float(np.nanmax(np.abs(Bk.to_numpy(float)))), 0.05)
    cm = FE9.diverging()
    W, labw, cell = 128.0, 46.0, FE9.CELL_MM
    aw, ah, bw = cell * A.shape[1], cell * A.shape[0], cell * Bk.shape[1]
    ticks = 1.3 * max(len(c) for c in GROUP_NAME.values()) + 3 + 9
    H_ = 6 + ah + ticks + 4
    fig = plt.figure(figsize=(W * MM, H_ * MM))

    def axes(x, y_top, w, h):
        return fig.add_axes([x / W, 1 - (y_top + h) / H_, w / W, h / H_])

    def letter(x, t_):
        fig.text(x / W, 1 - 3 / H_, t_, fontsize=F.LETTER_PT, weight="bold", va="top")

    axA = axes(labw, 6, aw, ah)
    FE9.heat(axA, A, TwoSlopeNorm(vmin=-va, vcenter=0, vmax=va), cm)
    axA.set_xticklabels([GROUP_NAME[g] for _ in res_list for g in GROUPS], rotation=90)
    for i, r in enumerate(res_list):
        x0, x1 = (labw + cell * 3 * i + 0.6) / W, (labw + cell * 3 * (i + 1) - 0.6) / W
        yl = 1 - (6 + ah + ticks - 9) / H_
        fig.add_artist(plt.Line2D([x0, x1], [yl, yl], color=F.INK, lw=0.6))
        fig.text((x0 + x1) / 2, yl - 0.8 / H_, f"{r:g}", ha="center", va="top", fontsize=F.FONT_PT)
    fig.text((labw + aw / 2) / W, 1 - (6 + ah + ticks - 1) / H_, "Louvain resolution", ha="center", va="bottom",
             fontsize=F.FONT_PT)
    letter(2, "A")
    caxA = axes(labw + aw + 2, 6, 2.5, ah)
    fig.colorbar(plt.cm.ScalarMappable(TwoSlopeNorm(vmin=-va, vcenter=0, vmax=va), cm), cax=caxA, alpha=F.ALPHA)
    caxA.set_ylabel("Recovery gain", labelpad=2)
    caxA.yaxis.set_label_position("right")
    xb = labw + aw + 2 + 2.5 + 22
    axB = axes(xb, 6, bw, ah)
    FE9.heat(axB, Bk, TwoSlopeNorm(vmin=-vb, vcenter=0, vmax=vb), cm)
    axB.set_xticklabels(Bk.columns, rotation=90)
    axB.set_yticklabels([])
    letter(xb - 6, "B")
    caxB = axes(xb + bw + 2, 6, 2.5, ah)
    fig.colorbar(plt.cm.ScalarMappable(TwoSlopeNorm(vmin=-vb, vcenter=0, vmax=vb), cm), cax=caxB, alpha=F.ALPHA)
    caxB.set_ylabel("kNN accuracy lift", labelpad=2)
    caxB.yaxis.set_label_position("right")
    fig.text(0.5, 0.004, "DEVELOPMENT FIGURE - made by a script (fig_e8.py)", ha="center", va="bottom",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_e8.{ext}"), dpi=300)
    clusters = t.groupby("resolution")["clusters"].agg(["min", "max"])
    open(os.path.join(out, "fig_e8_legend.md"), "w").write(
        "**Figure E8. Clustering on E1-selected features does not beat random feature sets of the same size.** "
        "PRELIMINARY (development figure). "
        "Features were selected by the E1 rule (Benjamini-Hochberg q < 0.05 and effect above the noise floor, "
        "5 seeds) on development genes only (80 % of gene families); grading used the untouched 20 % "
        "(final split, seed 2026). "
        "(A) Recovery gain: class-averaged share of untouched genes placed in a cluster enriched for their class, "
        "minus the same with development labels shuffled (5 shuffles); Louvain on the cosine kNN graph (k = 15) "
        f"at resolution 0.5 ({clusters.loc[0.5, 'min']}-{clusters.loc[0.5, 'max']} clusters) and 1.0 "
        f"({clusters.loc[1.0, 'min']}-{clusters.loc[1.0, 'max']} clusters). "
        "(B) kNN lift: accuracy of a cosine kNN vote from development to untouched genes minus the majority-class "
        "rate. Columns: all views (reference); selected features; random feature sets of the selected size "
        "(mean of 3 draws); in A they are grouped by Louvain resolution. Recovery gain averages over classes, so for "
        "the two binary host labels the large negative class inflates it (compare their kNN lift of about 0 in B). "
        "White = 0 (no gain over shuffled labels or chance). "
        "Abbreviations: kNN, k nearest neighbours; EC, Enzyme Commission; CORUM, comprehensive resource of "
        "mammalian protein complexes; STRING, search tool for the retrieval of interacting genes/proteins.\n")
    print(os.path.join(out, "fig_e8.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
