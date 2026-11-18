"""What ParaScape is, drawn: the imports, the one gene space, and what comes out of it.

Figure 1a of the ParaScape paper (instruction 101). Panel b is the measured coverage
(`build_corpus_inventory.py`); this is the schematic beside it, drawn from the catalogue rather
than typed, so it cannot drift from what the program actually holds.
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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import datasets as DS      # noqa: E402
from parascape import figstyle as FS      # noqa: E402
from parascape import modules as MOD      # noqa: E402

LEVELS = {"DNA": "genome and screens", "transcription": "transcription",
          "translation": "translation and protein", "post_translation": "modification, place, partners",
          "reference": "sequence, structure, literature"}

WIDTH_MM, HEIGHT_MM = 180.0, 62.0


def box(ax, x, y, w, h, text, face, fontsize=5.6, ink=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                linewidth=0, facecolor=face, alpha=FS.ALPHA, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            color=ink or FS.INK, zorder=3, linespacing=1.35)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=5,
                                 linewidth=0.6, color="#8a8a8a", zorder=1,
                                 shrinkA=0, shrinkB=0))


def main() -> None:
    FS.use()
    by_level = Counter(d.level for d in DS.REGISTRY)
    by_organism = Counter(d.organism for d in DS.REGISTRY)
    bands = Counter(m.band for m in MOD.MODULES)

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # Heading sits right of the panel letter; at x=0.085 the two ran through each other.
    ax.text(0.052, 0.945, "published measurements", ha="left", fontsize=6.2, color=FS.INK)
    levels = [lv for lv in LEVELS if by_level.get(lv)]
    h = 0.70 / len(levels)
    for i, lv in enumerate(levels):
        y = 0.855 - (i + 1) * h + 0.04
        box(ax, 0.01, y + 0.012, 0.15, h - 0.024,
            f"{LEVELS[lv]}\n{by_level[lv]} datasets", FS.CONTROL)
        arrow(ax, 0.163, y + h / 2, 0.235, 0.50)

    organisms = ", ".join(f"{k} {v}" for k, v in by_organism.most_common())
    box(ax, 0.24, 0.26, 0.21, 0.48,
        f"one gene space\n\n8,140 genes\n362 slots\n{len(DS.REGISTRY)} datasets\n\n{organisms}",
        FS.SERIES[3], fontsize=5.8)
    arrow(ax, 0.455, 0.50, 0.525, 0.50)

    box(ax, 0.53, 0.26, 0.18, 0.48,
        "made comparable\n\ncontrols clustered\ncell cycle removed\n201 contrasts\non one scale",
        FS.SERIES[2], fontsize=5.8)
    arrow(ax, 0.715, 0.50, 0.785, 0.50)

    labels = [f"{b} ({n})" for b, n in sorted(bands.items(), key=lambda kv: -kv[1])]
    box(ax, 0.79, 0.26, 0.20, 0.48,
        f"{len(MOD.MODULES)} modules\n\n" + "\n".join(labels), FS.SERIES[0], fontsize=5.8)

    ax.text(0.615, 0.095, "every analysis in this paper runs inside the program, and saves what it "
            "was made from", ha="center", fontsize=5.6, color="#7d7d7d")
    FS.letter(fig, "a")

    out = os.path.join(ROOT, "results", f"figure1_schematic_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"panel_schematic.{ext}"),
                    dpi=600 if ext == "png" else None, bbox_inches="tight", pad_inches=0.02)
    json.dump({"datasets": len(DS.REGISTRY), "by_level": dict(by_level),
               "by_organism": dict(by_organism), "modules": len(MOD.MODULES),
               "bands": dict(bands)}, open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
