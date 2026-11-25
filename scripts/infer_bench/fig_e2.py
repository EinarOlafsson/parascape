"""Figure of E2: how much of each label's predictability is explained by a confounder.

DEVELOPMENT FIGURE from an E2 output folder. Heatmap: label x confounder, median change of the per-feature score
after the smooth (spline) correction for that confounder, over the label's informative features (E1); the
raw median score is the first column's colour bar companion in the legend. Cells marked with a dot: the label
itself depends on the confounder (E2 dependence test).

    python scripts/infer_bench/fig_e2.py E2_DIR OUT_DIR
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

MM = 1 / 25.4


def main(argv=None) -> int:
    src, out = (argv if argv is not None else sys.argv[1:])
    os.makedirs(out, exist_ok=True)
    F.use()
    t = pd.read_parquet(os.path.join(src, "confounders.parquet"))
    sm = t[t["method"] == "smooth"]
    M = sm.pivot_table(index="label", columns="confounder", values="change_from_raw", aggfunc="median")
    dep = sm.groupby(["label", "confounder"])["label_depends_on_confounder"].first().unstack().reindex_like(M)
    M.index = [FE9.FE8_SHORT.get(i, i) for i in M.index]
    dep.index = M.index
    v = max(float(np.nanmax(np.abs(M.to_numpy()))), 0.05)
    norm = TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)
    cell, labw = FE9.CELL_MM * 1.4, 62.0
    aw, ah = cell * M.shape[1], cell * M.shape[0]
    W, H_ = labw + aw + 30, 8 + ah + 34
    fig = plt.figure(figsize=(W * MM, H_ * MM))
    ax = fig.add_axes([labw / W, 1 - (8 + ah) / H_, aw / W, ah / H_])
    FE9.heat(ax, M, norm, FE9.diverging())
    ax.set_xticklabels(M.columns, rotation=90)
    for (i, j), d in np.ndenumerate(dep.to_numpy()):
        if d is True or d == 1:
            ax.plot(j, i, "o", ms=2.2, color=F.INK)
    cax = fig.add_axes([(labw + aw + 3) / W, 1 - (8 + ah) / H_, 2.5 / W, ah / H_])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=FE9.diverging()), cax=cax)
    cb.solids.set_alpha(F.ALPHA)
    cb.set_label("Score change", labelpad=2)
    fig.text(0.995, 0.995, "DEVELOPMENT FIGURE - made by a script (fig_e2.py)", ha="right", va="top",
             color=F.DARK_GREY, fontsize=F.NOTE_PT)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_e2.{ext}"), dpi=300)
    raw = t[t["method"] == "raw"].groupby("label")["score"].median()
    open(os.path.join(out, "fig_e2_legend.md"), "w").write(
        "**Figure E2. How much of each label's predictability a confounder explains.** PRELIMINARY (development "
        "figure). For every informative label x feature pair of E1 (1,626 pairs), the per-feature kNN score "
        "(E1 score, family-split 5-fold, 3 seeds) was recomputed after removing the confounder from the feature "
        "with a cross-validated smoothing spline; cells: median change from the uncorrected score over the "
        "label's pairs (magenta = the confounder explains part of the information; white = no change). Dots: the "
        "label itself depends on the confounder (dependence above the E2 threshold). Confounders: number of "
        "papers on the gene (study attention), share of datasets that measured the gene, peak expression, protein "
        "length, codon GC3 and conservation breadth. Median uncorrected scores: "
        + "; ".join(f"{k} {v:.2f}" for k, v in raw.items()) +
        ". Chance is 0.5 for class labels and 0 for continuous scores.\n")
    print(os.path.join(out, "fig_e2.png"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
