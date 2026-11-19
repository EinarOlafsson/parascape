"""Figure 1a: what ParaScape can DO, as a schematic.

Einar, 2026-11-19: "your flow chart schematic is not only ugly it is not informative and does not
cover the capabilities of parascape at all", and -- about my reply, which talked about box heights
and arrow geometry -- "shows me that you have completely missed the point".

He is right, and the point is CONTENT. The old panel was a four-box left-to-right pipeline:
sources, one gene space, made comparable, modules. That describes how the corpus was assembled. It
says nothing about what a person can DO with it, which is what a capability schematic is for, and
tidying its geometry did not change that by one word.

WHAT THIS ONE SHOWS. The three things you can bring in (a gene, a list, your own measurements), the
one map they join, and the four things you can do with it -- look it up, explore it, ASK it, and
take it away. The Ask cluster is drawn largest and in the accent because it is the largest band
(8 of the 25 modules) and because inference is the point of the program: its glyph is a gene with
three facts measured and a fourth, missing, being filled in.

HOW IT OBEYS THE SKILL. Illustrated, not described: every capability is a glyph and the only text
is one short label beside each, nine words in all, under the limit of ten. Flat shapes, no
gradients or shadows, light-grey bodies with thin darker-grey outlines, ONE accent for the element
under study. Tight bounding box, no larger than one data panel.
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import datasets as DS       # noqa: E402
from parascape import figstyle as FS       # noqa: E402
from parascape import modules as MOD       # noqa: E402

#: A schematic is no larger than one data panel (skill, 2026-11-18).
WIDTH_MM, HEIGHT_MM = 92.0, 42.0

#: Schematic bodies: light grey with a thin darker outline (skill).
BODY = "#E4E4E4"
BODY_DARK = "#D0D0D0"
EDGE = "#9A9A9A"


def body(ax, x, y, w, h, face=BODY, lw=0.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.002,rounding_size=0.012",
                                facecolor=face, edgecolor=EDGE, linewidth=lw, zorder=2))


def word(ax, x, y, text, colour=None, weight="normal"):
    ax.text(x, y, text, ha="center", va="top", fontsize=FS.FONT_PT,
            color=colour or FS.INK, fontweight=weight, zorder=5)


def arrow(ax, x0, y0, x1, y1, colour=None):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=5,
                                 linewidth=0.6, color=colour or "#000000", zorder=4,
                                 shrinkA=0, shrinkB=0))


def glyph_gene(ax, cx, cy, s):
    """A gene: a bar with a darker coding block."""
    body(ax, cx - s, cy - s * 0.24, s * 2, s * 0.48)
    ax.add_patch(Rectangle((cx - s * 0.30, cy - s * 0.24), s * 0.60, s * 0.48,
                           facecolor=BODY_DARK, edgecolor="none", zorder=3))


def glyph_list(ax, cx, cy, s):
    """A gene list: stacked rules."""
    for i, dy in enumerate((0.34, 0.0, -0.34)):
        ax.add_patch(Rectangle((cx - s * 0.85, cy + s * dy - s * 0.09), s * 1.7, s * 0.18,
                               facecolor=BODY if i % 2 == 0 else BODY_DARK,
                               edgecolor=EDGE, linewidth=0.5, zorder=3))


def glyph_upload(ax, cx, cy, s):
    """Your own measurements: a sheet with an arrow going in."""
    body(ax, cx - s * 0.70, cy - s * 0.50, s * 1.4, s * 1.0)
    arrow(ax, cx, cy - s * 0.26, cx, cy + s * 0.24)


def glyph_map(ax, cx, cy, s):
    """The map: a node-link graph, which is what the corpus becomes."""
    nodes = [(-0.72, 0.32), (0.02, 0.66), (0.72, 0.28), (-0.46, -0.44), (0.42, -0.52), (0.0, 0.0)]
    for a in range(len(nodes)):
        for b in range(a + 1, len(nodes)):
            ax.plot([cx + nodes[a][0] * s, cx + nodes[b][0] * s],
                    [cy + nodes[a][1] * s, cy + nodes[b][1] * s],
                    color=EDGE, linewidth=0.4, zorder=2)
    for i, (dx, dy) in enumerate(nodes):
        ax.add_patch(Circle((cx + dx * s, cy + dy * s), s * 0.14,
                            facecolor=FS.CONTROL if i == len(nodes) - 1 else BODY_DARK,
                            edgecolor=EDGE, linewidth=0.5, zorder=3))


def glyph_card(ax, cx, cy, s):
    """Look up: a record card with filled fields."""
    body(ax, cx - s * 0.72, cy - s * 0.52, s * 1.44, s * 1.04)
    for i, dy in enumerate((0.26, 0.0, -0.26)):
        ax.add_patch(Rectangle((cx - s * 0.50, cy + s * dy - s * 0.07),
                               s * (0.86 - i * 0.16), s * 0.14,
                               facecolor=BODY_DARK, edgecolor="none", zorder=3))


def glyph_explore(ax, cx, cy, s):
    """Explore: a cloud of genes with one group picked out."""
    pts = [(-0.62, 0.34), (-0.28, 0.60), (-0.66, -0.12), (-0.20, 0.04), (0.06, 0.46),
           (0.46, 0.30), (0.66, -0.06), (0.30, -0.36), (-0.10, -0.56), (0.60, 0.52)]
    for i, (dx, dy) in enumerate(pts):
        ax.add_patch(Circle((cx + dx * s, cy + dy * s), s * 0.12,
                            facecolor=BODY_DARK if i < 7 else FS.CONTROL,
                            edgecolor="none", zorder=3))


def glyph_ask(ax, cx, cy, s):
    """Ask: a gene with three facts measured and a fourth, missing, being filled in.

    The one element under study, so it carries the only accent. This glyph is the argument of the
    paper: the map exists so that a fact nobody measured can be predicted from the ones that were.
    """
    for dy in (0.54, 0.20, -0.14):
        ax.add_patch(Rectangle((cx - s * 0.66, cy + s * dy - s * 0.10), s * 1.32, s * 0.20,
                               facecolor=BODY_DARK, edgecolor=EDGE, linewidth=0.5, zorder=3))
    ax.add_patch(Rectangle((cx - s * 0.66, cy - s * 0.58), s * 1.32, s * 0.20,
                           facecolor=FS.ACCENT, edgecolor=FS.ACCENT, linewidth=0.5,
                           alpha=FS.ALPHA, zorder=3))
    arrow(ax, cx + s * 1.16, cy - s * 0.48, cx + s * 0.72, cy - s * 0.48, colour=FS.ACCENT)


def glyph_export(ax, cx, cy, s):
    """Export: a sheet with an arrow leaving it."""
    body(ax, cx - s * 0.70, cy - s * 0.50, s * 1.4, s * 1.0)
    arrow(ax, cx, cy + s * 0.24, cx, cy - s * 0.26)


def main() -> None:
    FS.use()
    bands = Counter(m.band for m in MOD.MODULES)

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    s = 0.070

    # ---- what you bring in
    for cy, draw, text in ((0.82, glyph_gene, "Gene"),
                           (0.52, glyph_list, "List"),
                           (0.22, glyph_upload, "Your data")):
        draw(ax, 0.080, cy, s)
        word(ax, 0.080, cy - s * 0.74, text)
        arrow(ax, 0.152, cy, 0.212, 0.52)

    # ---- the one map everything joins
    glyph_map(ax, 0.300, 0.52, s * 1.40)
    word(ax, 0.300, 0.52 - s * 1.40 * 0.78, "Map")

    # ---- what you can do with it: three capabilities OUT of the map, all feeding Export.
    # The first draft had Look up and Explore pointing INTO Ask and Export pointing backwards
    # into it, which says something false about the program.
    glyph_card(ax, 0.600, 0.85, s * 0.95)
    word(ax, 0.600, 0.85 - s * 0.72, "Look up")
    glyph_explore(ax, 0.600, 0.17, s * 0.95)
    word(ax, 0.600, 0.17 - s * 0.72, "Explore")
    glyph_ask(ax, 0.610, 0.52, s * 1.30)
    word(ax, 0.610, 0.52 - s * 1.30 * 0.76, "Ask", colour=FS.ACCENT, weight="bold")
    glyph_export(ax, 0.905, 0.52, s * 0.95)
    word(ax, 0.905, 0.52 - s * 0.72, "Export")

    # map -> each capability
    arrow(ax, 0.392, 0.60, 0.528, 0.83)
    arrow(ax, 0.392, 0.52, 0.516, 0.52)
    arrow(ax, 0.392, 0.44, 0.528, 0.21)
    # each capability -> export
    arrow(ax, 0.672, 0.83, 0.836, 0.60)
    arrow(ax, 0.706, 0.52, 0.836, 0.52)
    arrow(ax, 0.672, 0.21, 0.836, 0.44)

    out = os.path.join(ROOT, "results", f"figure1_schematic_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    for extension in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"panel_schematic.{extension}"),
                    dpi=600 if extension == "png" else None,
                    bbox_inches="tight", pad_inches=0.01)
    labels = ["Gene", "List", "Your data", "Map", "Look up", "Explore", "Ask", "Export"]
    n_words = sum(len(w.split()) for w in labels)
    json.dump({"datasets": len(DS.REGISTRY), "modules": len(MOD.MODULES), "bands": dict(bands),
               "words_in_schematic": n_words,
               "shows": "what can be brought in, the map they join, and what can be done with it"},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print(f"wrote {out}  ({n_words} words in the schematic, limit 10)")
    print("bands:", dict(bands))


if __name__ == "__main__":
    main()
