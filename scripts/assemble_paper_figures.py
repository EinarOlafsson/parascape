"""Assemble the panels into FIGURES, one image per figure.

Einar, 2026-11-19: "dont show me single panels show me figures ... white space should be minimized
in graphs and between panels."

So this lays each figure's panels out side by side at their true size, with a fixed small gutter
between them and nothing else. A panel is placed at the millimetre size it was drawn at -- it is
never rescaled, because rescaling a panel changes its font size, and the whole point of sizing
everything in ems is that 8 pt is 8 pt on every panel of every figure.

Panels whose measurement does not exist yet are drawn with PLACEHOLDER DATA in the real graph form
(see `scripts/build_paper_figure_sheet.py`), so the choice of graph can be judged before the data
arrives, rather than with a hatched box that shows nothing.

Writes `results/paper_figures_<date>/figureNN.png` and `.pdf`.
"""
from __future__ import annotations

import glob
import json
import os
import sys
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import figstyle as FS      # noqa: E402

#: The gutter between panels, in ems, so it scales with the font like everything else.
GUTTER_EM = 1.2
GUTTER_MM = GUTTER_EM * FS.EM_MM

#: Panels are rendered at this many dots per millimetre when composed.
DPMM = 600 / 25.4


def newest_sheet() -> str:
    found = sorted(glob.glob(os.path.join(ROOT, "results", "paper_figure_sheet_*")))
    if not found:
        raise SystemExit("run scripts/build_paper_figure_sheet.py first")
    return found[-1]


def compose(figure_id: str, panels: list, out: str) -> str:
    """One figure from its panels, laid out in a row at their true millimetre size."""
    images = [mpimg.imread(p["png"]) for p in panels]
    # True size in millimetres: the panel PNGs are written at 600 dpi.
    widths = [img.shape[1] / DPMM for img in images]
    heights = [img.shape[0] / DPMM for img in images]
    total_w = sum(widths) + GUTTER_MM * (len(images) - 1)
    total_h = max(heights)

    fig = plt.figure(figsize=(total_w / 25.4, total_h / 25.4))
    x = 0.0
    for img, w, h in zip(images, widths, heights):
        # Panels in a row share their TOP and BOTTOM edges (skill checklist, item 3). A panel is
        # never rescaled -- scaling changes its font size -- so a shorter panel is centred in the
        # row height and the difference is white space, not a different type size.
        y = (total_h - h) / 2.0
        ax = fig.add_axes([x / total_w, y / total_h, w / total_w, h / total_h])
        ax.imshow(img)
        ax.set_axis_off()
        x += w + GUTTER_MM
    path = os.path.join(out, f"figure{figure_id}")
    for extension in ("png", "pdf"):
        fig.savefig(f"{path}.{extension}", dpi=600 if extension == "png" else None,
                    bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)
    return f"{path}.png"


def main() -> None:
    sheet = newest_sheet()
    summary = json.load(open(os.path.join(sheet, "summary.json")))
    out = os.path.join(ROOT, "results", f"paper_figures_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)

    by_figure: dict[str, list] = {}
    for row in summary["rows"]:
        png = os.path.join(sheet, f"fig{row['figure']}{row['letter']}.png")
        if not os.path.exists(png):
            print(f"  missing panel image: {png}")
            continue
        by_figure.setdefault(row["figure"], []).append({**row, "png": png})

    for figure_id in sorted(by_figure):
        panels = sorted(by_figure[figure_id], key=lambda r: r["letter"])
        path = compose(figure_id, panels, out)
        letters = "".join(p["letter"] for p in panels)
        real = sum(1 for p in panels if p["real"])
        note = "" if real == len(panels) else f"  ({len(panels) - real} with placeholder data)"
        print(f"  figure {figure_id}  panels {letters}{note}  -> {os.path.basename(path)}")

    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
