"""Figure 1A: ParaScape as a flowchart -- input data -> what it does -> what it is used for.

Einar, 2026-11-20 (second rejection of 1A): "it should be a flowchart of the components of
parascape and the data that the package takes as input, then what it does with the data and then
what it uses the data for". Boxes and arrows, few words per box, no lists of text, no charts.
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from parascape import figstyle as FS                    # noqa: E402

INPUTS = ["RNA-seq", "Proteomics,\nhyperLOPIT", "CRISPR\nscreens", "Interactions\nXL-MS, IP-MS",
          "Structures\nAlphaFold", "Single-cell\nRNA-seq", "Plasmodium\n& host data"]
STEPS = [("Re-quantify raw reads", "salmon, kallisto"),
         ("Map to one gene ID", "VEuPathDB, orthology"),
         ("Remove confounders", "cell cycle, laboratory"),
         ("Gene table and graph", "traits, 12 evidence types"),
         ("Predict,\nheld-out test", "localisation,\nphase, function")]
USES = ["Look up a gene\nor gene list", "Find co-expression\nmodules", "Localise proteins",
        "Transfer to\nPlasmodium", "Rank host-facing\nproteins and targets",
        "Choose the next\nexperiment"]


def box(ax, x, y, w, h, text, fc, ec, bold=False, sub=None):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=1.2", lw=0.6, ec=ec, fc=fc))
    if sub:
        ax.text(x, y + h / 2 - 1.2, text, ha="center", va="top", fontsize=FS.FONT_PT,
                fontweight="bold" if bold else "normal", color=FS.INK, linespacing=1.0)
        ax.text(x, y - h / 2 + 1.0, sub, ha="center", va="bottom", fontsize=FS.FONT_PT,
                color=FS.DARK_GREY, linespacing=1.0)
    else:
        ax.text(x, y, text, ha="center", va="center", fontsize=FS.FONT_PT,
                fontweight="bold" if bold else "normal", color=FS.INK, linespacing=1.0)


def arrow(ax, p, q, color=FS.INK):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=6, lw=0.7, color=color,
                                 shrinkA=0, shrinkB=0))


def draw():
    """Three stages, each under a coloured header bar; equal boxes; grey connectors."""
    FS.use()
    W, H = 180.0, 70.0
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    pale, grey, line = "#E6F0F8", "#F0F0F0", "#7A7A7A"
    cols = {"in": (2, 36), "mid": (46, 134), "out": (144, 178)}

    def header(key, text):
        x0, x1 = cols[key]
        ax.add_patch(FancyBboxPatch((x0, H - 7), x1 - x0, 5.2,
                                    boxstyle="round,pad=0,rounding_size=1.0", lw=0,
                                    fc=FS.ACCENT if key == "mid" else FS.DARK_GREY))
        ax.text((x0 + x1) / 2, H - 4.4, text, ha="center", va="center", color="white",
                fontsize=FS.FONT_PT, fontweight="bold")
    header("in", "Published data"); header("mid", "ParaScape"); header("out", "Uses")

    def column(key, items, fc, ec):
        x0, x1 = cols[key]
        n = len(items); top, bot = H - 9.5, 2
        hh = (top - bot) / n
        ys = []
        for k, t in enumerate(items):
            y = top - hh * (k + 0.5)
            ax.add_patch(FancyBboxPatch((x0, y - hh / 2 + 0.6), x1 - x0, hh - 1.2,
                                        boxstyle="round,pad=0,rounding_size=1.0", lw=0.5,
                                        ec=ec, fc=fc))
            ax.text((x0 + x1) / 2, y, t, ha="center", va="center", fontsize=FS.FONT_PT,
                    linespacing=1.0)
            ys.append(y)
        return ys
    yin = column("in", INPUTS, grey, line)
    yout = column("out", USES, "white", line)

    # pipeline: 2 x 2 grid of steps plus prediction, read left-to-right, top-to-bottom
    bw, bh = 38, 15
    grid = [(65, H - 20), (115, H - 20), (65, H - 41), (115, H - 41)]
    for (x, y), (t, sub) in zip(grid, STEPS[:4]):
        ax.add_patch(FancyBboxPatch((x - bw / 2, y - bh / 2), bw, bh,
                                    boxstyle="round,pad=0,rounding_size=1.2", lw=0.6,
                                    ec=FS.ACCENT, fc=pale))
        ax.text(x, y + 2.6, t.replace("\n", " "), ha="center", va="center",
                fontsize=FS.FONT_PT, fontweight="bold")
        ax.text(x, y - 3.0, sub.replace("\n", " "), ha="center", va="center",
                fontsize=FS.FONT_PT, color=FS.DARK_GREY)
    px, py = 90, H - 61
    ax.add_patch(FancyBboxPatch((px - 44, py - bh / 2), 88, bh,
                                boxstyle="round,pad=0,rounding_size=1.2", lw=0.6, ec=FS.ACCENT,
                                fc=FS.ACCENT_LIGHT))
    ax.text(px, py + 2.6, "Predict and test on held-out genes", ha="center", va="center",
            fontsize=FS.FONT_PT, fontweight="bold")
    ax.text(px, py - 3.0, "localisation, cell-cycle phase, function; your gene list",
            ha="center", va="center", fontsize=FS.FONT_PT, color=FS.INK)

    def arr(p, q):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=6, lw=0.7,
                                     color=line, shrinkA=0, shrinkB=0,
                                     connectionstyle="arc3,rad=0"))
    # inputs converge into step 1
    xj = 40.5
    ax.plot([xj, xj], [yin[-1], yin[0]], color=line, lw=0.7)
    for y in yin:
        ax.plot([cols["in"][1], xj], [y, y], color=line, lw=0.7)
    arr((xj, grid[0][1]), (grid[0][0] - bw / 2, grid[0][1]))
    arr((grid[0][0] + bw / 2, grid[0][1]), (grid[1][0] - bw / 2, grid[1][1]))
    arr((grid[1][0], grid[1][1] - bh / 2), (grid[3][0], grid[3][1] + bh / 2))
    arr((grid[3][0] - bw / 2, grid[3][1]), (grid[2][0] + bw / 2, grid[2][1]))
    arr((grid[2][0], grid[2][1] - bh / 2), (grid[2][0], py + bh / 2))
    # outputs fan out from the pipeline
    xo = 139.5
    ax.plot([xo, xo], [yout[-1], yout[0]], color=line, lw=0.7)
    arr((px + 44, py), (xo, py)) if False else ax.plot([px + 44, xo], [py, py], color=line, lw=0.7)
    ax.plot([grid[1][0] + bw / 2, xo], [grid[1][1], grid[1][1]], color=line, lw=0.7)
    for y in yout:
        arr((xo, y), (cols["out"][0], y))
    return fig


if __name__ == "__main__":
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    out = os.path.join(ROOT, "results", "figure01_flowchart")
    os.makedirs(out, exist_ok=True)
    f = draw()
    f.savefig(os.path.join(out, "panel_a.png"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    f.savefig(os.path.join(out, "panel_a.pdf"), bbox_inches="tight", pad_inches=0.02)
    print(out)
