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
STEPS = [("Re-quantify\nraw reads", "salmon, kallisto"),
         ("One gene\nidentity", "VEuPathDB IDs,\northology"),
         ("Gene table\n+ gene graph", "traits, 12 kinds\nof link"),
         ("Remove\nconfounders", "cell cycle,\nlaboratory"),
         ("Predict,\nheld-out test", "localisation,\nphase, function")]
USES = ["Look up a gene\nor a gene list", "Find modules\nand organelles", "Place proteins",
        "Transfer to\nPlasmodium", "Rank host-facing\nproteins, drug targets",
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
    FS.use()
    W, H = 180.0, 78.0
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    light = "#EEF4FA"; grey = "#F2F2F2"

    # column headers
    for x, t in ((17, "Input: published data"), (90, "ParaScape"), (163, "Used to")):
        ax.text(x, H - 2, t, ha="center", va="top", fontsize=FS.FONT_PT, fontweight="bold")

    # inputs
    n = len(INPUTS); top, bot = H - 9, 4
    step = (top - bot) / n
    ys = [top - step * (i + 0.5) for i in range(n)]
    for y, t in zip(ys, INPUTS):
        box(ax, 17, y, 30, step - 1.6, t, grey, FS.DARK_GREY)
    # bracket of inputs into the pipeline
    xm = 36
    ax.plot([xm, xm], [ys[-1], ys[0]], color=FS.INK, lw=0.7)
    for y in ys:
        ax.plot([32, xm], [y, y], color=FS.INK, lw=0.7)

    # pipeline: two rows, left to right then down, as a flow
    px = [52, 76, 100, 124]
    box_w, box_h = 21, 21
    y1 = H / 2 + 11
    arrow(ax, (xm, H / 2), (px[0] - box_w / 2, y1))
    pos = [(px[0], y1), (px[1], y1), (px[2], y1), (px[3], y1)]
    for (x, y), (t, sub) in zip(pos, STEPS[:4]):
        box(ax, x, y, box_w, box_h, t, light, FS.ACCENT, bold=True, sub=sub)
    for (a, _), (b, _) in zip(pos[:-1], pos[1:]):
        arrow(ax, (a + box_w / 2, y1), (b - box_w / 2, y1))
    # prediction under the table, fed by table and corrected data
    y2 = H / 2 - 14
    x5 = (px[1] + px[2]) / 2
    box(ax, x5 + 12, y2, 30, box_h, STEPS[4][0], light, FS.ACCENT, bold=True, sub=STEPS[4][1])
    arrow(ax, (px[3], y1 - box_h / 2), (x5 + 22, y2 + box_h / 2))
    # user's own data enters too
    box(ax, px[0], y2, 21, box_h, "Your gene\nlist or data", "white", FS.ACCENT)
    arrow(ax, (px[0] + 10.5, y2), (x5 - 3, y2))
    # frame around ParaScape
    ax.add_patch(FancyBboxPatch((40, 5), 96, H - 13, boxstyle="round,pad=0,rounding_size=2",
                                lw=0.6, ec=FS.ACCENT, fc="none", ls=(0, (3, 2))))

    # uses
    m = len(USES); step2 = (top - bot) / m
    ys2 = [top - step2 * (i + 0.5) for i in range(m)]
    xo = 144
    ax.plot([xo, xo], [ys2[-1], ys2[0]], color=FS.INK, lw=0.7)
    arrow(ax, (136, H / 2), (xo, H / 2))
    for y, t in zip(ys2, USES):
        arrow(ax, (xo, y), (148, y))
        box(ax, 163, y, 30, step2 - 1.6, t, "white", FS.DARK_GREY)
    return fig


if __name__ == "__main__":
    os.environ.setdefault("PARASCAPE_NO_PANEL_LETTERS", "1")
    out = os.path.join(ROOT, "results", "figure01_flowchart")
    os.makedirs(out, exist_ok=True)
    f = draw()
    f.savefig(os.path.join(out, "panel_a.png"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    f.savefig(os.path.join(out, "panel_a.pdf"), bbox_inches="tight", pad_inches=0.02)
    print(out)
