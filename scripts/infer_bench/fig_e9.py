"""Figure of E9: which features to cluster on and how to cluster, with positive and negative controls.

DEVELOPMENT FIGURE made by a script from an E9 output folder (`e9_strategies.py`). Panels:

A. stage A heatmap: label x feature strategy, held-out recovery gain over shuffled labels
   (Louvain 0.5); the best real strategy per label is boxed; random sets are the mean of 3 draws;
B. positive control: planted modules, recovery gain against planted strength, per ranking method
   (top 30), with all views (reference, grey) and random 30 (dashed grey);
C. positive control: share of the 20 planted features each ranking method finds in its top 30
   and in its significance set;
D. stage B heatmap: label x clustering method, on all views and on the stage A winner;
E. biological positive controls: recovery of the ribosome and proteasome classes.

    python scripts/infer_bench/fig_e9.py E9_DIR OUT_DIR
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parascape import figstyle as F  # noqa: E402
import e9_strategies as E9  # noqa: E402

METHODS = ("knn", "mutual information", "kruskal", "forest")
METHOD_NAME = {"knn": "neighbour score", "mutual information": "mutual information",
               "kruskal": "Kruskal-Wallis", "forest": "random forest"}
METHOD_COLOUR = dict(zip(METHODS, (F.ACCENT, F.TEAL, F.NAVY, F.MAGENTA)))
CELL_MM = 4.2
MM = 1 / 25.4


def diverging():
    """Pale magenta below 0, white at 0 (chance), blue above."""
    return LinearSegmentedColormap.from_list("gain", ["#C878A8", "#FFFFFF", F.ACCENT])


def strategy_columns(a: pd.DataFrame) -> list:
    """Stage A columns in reading order: reference, each method (top 10/30/100, significant), random."""
    cols = ["all views"]
    for m in METHODS:
        cols += [f"{m} top {k}" for k in E9.TOPK]
        if f"{m} significant" in set(a["family"]):
            cols.append(f"{m} significant")
    cols += [f"random {k}" for k in E9.TOPK]
    return [c for c in cols if c in set(a["family"])]


def column_label(c: str) -> str:
    """Short x tick text."""
    for m in METHODS:
        if c.startswith(m):
            return f"{METHOD_NAME[m]}, {c[len(m) + 1:]}"
    return c


def heat(ax, M: pd.DataFrame, norm, cmap, box=None):
    """Heatmap of M (rows x columns), cells at 60 %, NaN cells marked with a thin black x."""
    ax.imshow(M.to_numpy(float), cmap=cmap, norm=norm, aspect="equal", alpha=F.ALPHA, interpolation="none")
    for (i, j), v in np.ndenumerate(M.to_numpy(float)):
        if not np.isfinite(v):
            ax.plot(j, i, marker="x", color=F.INK, ms=3, mew=0.5)
    if box is not None:
        for i, j in box:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=F.INK, lw=0.8))
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([column_label(c) for c in M.columns], rotation=90)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(M.index)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)


def main(argv=None) -> int:
    """Draw the figure into OUT_DIR as PNG (300 dpi) and SVG, plus the legend as markdown."""
    argv = argv if argv is not None else sys.argv[1:]
    if len(argv) != 2:
        raise SystemExit(__doc__)
    src, out = argv
    os.makedirs(out, exist_ok=True)
    F.use()
    a = pd.read_parquet(os.path.join(src, "stage_a.parquet"))
    b = pd.read_parquet(os.path.join(src, "stage_b.parquet"))
    pc = pd.read_parquet(os.path.join(src, "per_class.parquet"))
    ps = pd.read_parquet(os.path.join(src, "planted_selection.parquet"))
    a["family"] = a["strategy"].str.replace(r" #\d+$", "", regex=True)
    real = [x for x in pd.unique(a["label"]) if not x.startswith("planted")]
    planted = sorted([x for x in pd.unique(a["label"]) if x.startswith("planted")],
                     key=lambda s: float(s.rsplit(" ", 1)[1]))
    cols = strategy_columns(a)
    A = a.groupby(["label", "family"])["recovery_gain"].mean().unstack().reindex(index=real + planted, columns=cols)
    boxes = []
    for i, lab in enumerate(A.index):
        row = A.loc[lab, [c for c in cols if c != "all views" and not c.startswith("random")]]
        if row.notna().any() and row.max() > 0:          # a tie at 0 boxes nothing
            boxes.append((i, cols.index(row.idxmax())))
    vals = A.to_numpy(float)
    vmax = max(float(np.nanmax(np.abs(vals))) if np.isfinite(vals).any() else 0.0, 0.05)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    b = b.copy()
    b["col"] = np.where(b["method"].isin(["kmeans", "gmm", "ward"]), b["method"] + ", k = 2 x classes",
                        b["method"] + " " + b["setting"].map(lambda s: f"{float(s):g}"))
    Bpiv = b.pivot_table(index="label", columns=["features", "col"], values="recovery_gain")
    Bpiv = Bpiv.reindex(index=[x for x in real + planted if x in Bpiv.index])
    Dm = Bpiv.copy()
    Dm.columns = [f"{'winner' if f == 'stage A winner' else f}: {c}" for f, c in Dm.columns]
    # ---------------------------------------------------------------- layout (mm), three rows
    W, labw, sq, gap = 180.0, 62.0, 34.0, 6.0
    tick_mm = lambda names: 1.3 * max(len(str(n)) for n in names) + 3     # rotated 8 pt tick labels
    aw, ah = CELL_MM * A.shape[1], CELL_MM * A.shape[0]
    a_ticks = tick_mm([column_label(c) for c in A.columns])
    dw, dh = CELL_MM * Dm.shape[1], CELL_MM * Dm.shape[0]
    d_ticks = tick_mm(Dm.columns)
    ctrl = pc[pc["class"].isin(E9.CONTROL_CLASSES) & (pc["stage"] == "A")]
    y1 = 6.0
    y2 = y1 + ah + a_ticks + 10
    y3 = y2 + (sq + 16 if planted else 0)
    H_ = y3 + max(dh, CELL_MM * max(1, ctrl["class"].nunique())) + d_ticks + 12
    fig = plt.figure(figsize=(W * MM, H_ * MM))

    def axes(x, y_top, w, h):
        return fig.add_axes([x / W, 1 - (y_top + h) / H_, w / W, h / H_])

    def letter(x, y_top, t):
        fig.text(x / W, 1 - y_top / H_, t, fontsize=F.LETTER_PT, weight="bold", va="top")

    axA = axes(labw, y1, aw, ah)
    heat(axA, A, norm, diverging(), boxes)
    letter(2, y1 - 3, "A")
    cax = axes(labw + aw + 3, y1, 2.5, max(ah, 14))
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=diverging()), cax=cax)
    cb.solids.set_alpha(F.ALPHA)
    cb.set_label("Recovery gain")
    if planted:
        strengths = [float(s.rsplit(" ", 1)[1]) for s in planted]
        axB = axes(labw - 40, y2, sq, sq)
        for m in METHODS:
            v = [A.loc[p, f"{m} top 30"] if f"{m} top 30" in A else np.nan for p in planted]
            axB.plot(strengths, v, "-o", color=METHOD_COLOUR[m], alpha=F.ALPHA, lw=1.0, ms=3, label=METHOD_NAME[m])
        axB.plot(strengths, [A.loc[p, "all views"] for p in planted], "-o", color=F.CONTROL, alpha=F.ALPHA,
                 lw=1.0, ms=3, label="all views (ref.)")
        if "random 30" in A:
            axB.plot(strengths, [A.loc[p, "random 30"] for p in planted], "--o", color=F.CONTROL,
                     alpha=F.ALPHA, lw=1.0, ms=3, label="random 30 (ref.)")
        axB.set_xscale("log", base=2)
        axB.set_xticks(strengths)
        axB.set_xticklabels([f"{s:g}" for s in strengths])
        axB.set_xlabel("Planted strength (column SD)")
        axB.set_ylabel("Recovery gain")
        axB.set_ylim(0, 1)
        letter(2, y2 - 3, "B")
        xC = labw - 40 + sq + 22
        axC = axes(xC, y2, sq, sq)
        for m in METHODS:
            v = [ps.loc[(ps.label == p) & (ps.strategy == f"{m} top 30"), "recall"].mean() for p in planted]
            axC.plot(strengths, np.array(v) * 100, "-o", color=METHOD_COLOUR[m], alpha=F.ALPHA, lw=1.0, ms=3)
        for m in ("knn", "kruskal"):
            v = [ps.loc[(ps.label == p) & (ps.strategy == f"{m} significant"), "recall"].mean() for p in planted]
            axC.plot(strengths, np.array(v) * 100, ":s", color=METHOD_COLOUR[m], alpha=F.ALPHA, lw=1.0, ms=3)
        v = [ps.loc[(ps.label == p) & ps.strategy.str.startswith("random 30"), "recall"].mean() for p in planted]
        axC.plot(strengths, np.array(v) * 100, "--o", color=F.CONTROL, alpha=F.ALPHA, lw=1.0, ms=3)
        axC.set_xscale("log", base=2)
        axC.set_xticks(strengths)
        axC.set_xticklabels([f"{s:g}" for s in strengths])
        axC.set_xlabel("Planted strength (column SD)")
        axC.set_ylabel("Planted features found (%)")
        axC.set_ylim(0, 105)
        letter(xC - 14, y2 - 3, "C")
        handles, labels = axB.get_legend_handles_labels()
        handles.append(plt.Line2D([], [], color=F.INK, ls=":", marker="s", ms=3, lw=1.0))
        labels.append("significance set")
        fig.legend(handles, labels, loc="upper left", frameon=False,
                   bbox_to_anchor=((xC + sq + 4) / W, 1 - y2 / H_))
    axD = axes(labw, y3, dw, dh)
    heat(axD, Dm, norm, diverging())
    axD.set_xticklabels(list(Dm.columns), rotation=90)
    letter(2, y3 - 3, "D")
    if len(ctrl):
        win = {A.index[i]: A.columns[j] for i, j in boxes}
        keep = ["all views"] + sorted({win[l] for l in ctrl["label"].unique() if l in win})
        Em = ctrl.groupby(["class", "strategy"])["recovery"].mean().unstack().reindex(columns=keep)
        Em = Em.loc[:, Em.notna().any()]
        xE = W - 4 - CELL_MM * Em.shape[1]
        axE = axes(xE, y3, CELL_MM * Em.shape[1], CELL_MM * Em.shape[0])
        heat(axE, Em, plt.Normalize(0, 1), plt.cm.Blues)
        letter(labw + dw + 4, y3 - 3, "E")
    fig.text(0.5, 1.5 / H_, "DEVELOPMENT FIGURE - made by a script (fig_e9.py)", ha="center", color=F.DARK_GREY)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(out, f"fig_e9.{ext}"), dpi=300)
    with open(os.path.join(out, "fig_e9_legend.md"), "w") as f:
        f.write(LEGEND)
    print(os.path.join(out, "fig_e9.png"))
    return 0


LEGEND = """**Figure E9. Choosing, from the data, which features to cluster on and how to cluster.**
Every choice is made on development genes (80% of gene families) and graded on the untouched 20% (gene
families held back, seed 2026). Recovery: the share of untouched labelled genes that fall in a cluster
enriched (hypergeometric test, Benjamini-Hochberg q < 0.05) for their own class among development genes,
averaged over classes; recovery gain: recovery minus the recovery with development labels shuffled (5
shuffles), so 0 is chance.
**(A)** Recovery gain per label (rows) and feature strategy (columns), Louvain communities at resolution
0.5 on the cosine 15-nearest-neighbour graph. Features ranked on development labels by the neighbour
score of E1 (single-feature 15-nearest-neighbour score against shuffled labels), mutual information,
the Kruskal-Wallis test, or random-forest importance; top 10, 30 or 100 features, or the method own
significance set (neighbour score: q < 0.05 and a gain above twice the larger seed-to-seed SD;
Kruskal-Wallis: q < 0.05). Extra views (sequence models and structure) enter whole when their view
score passes, for every strategy alike. All views: the reference. Random: as many features drawn at
random, mean of 3 draws (the effect of fewer columns alone). Boxed: the best non-reference strategy of
each row. Rows "planted modules at s": positive control, 20 groups of 15 genes planted along one random
direction per group in 20 features (10 transcript, 10 fitness), strength s column SDs per feature.
x: a cell that could not be computed.
**(B)** Positive control: recovery gain against planted strength for the top 30 features of each
ranking method, all views (grey) and 30 random features (grey, dashed).
**(C)** Positive control: share of the 20 planted features found in the top 30 of each method (solid)
and in the significance sets (dotted squares); random 30 in grey, dashed.
**(D)** Recovery gain per label and clustering method, on all views and on the stage A winner of the
label: Louvain at resolution 0.5, 1 and 2; HDBSCAN with a minimum cluster of 10; k-means, Gaussian
mixture (diagonal covariance) and Ward at k = 2 x classes.
**(E)** Biological positive controls: recovery of the ribosome and proteasome classes (hyperLOPIT
niche; CORUM complexes copied by orthology), on all views and on the stage A winner.
"""

if __name__ == "__main__":
    sys.exit(main())
