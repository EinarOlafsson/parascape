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

#: A strip above each panel holds its letter, so the letter never sits on an axis title.
LETTER_STRIP_MM = 4.5

#: Panels are rendered at this many dots per millimetre when composed.
DPMM = 600 / 25.4


def newest_sheet() -> str:
    found = sorted(glob.glob(os.path.join(ROOT, "results", "paper_figure_sheet_*")))
    if not found:
        raise SystemExit("run scripts/build_paper_figure_sheet.py first")
    return found[-1]


#: Panels per row before wrapping. A single row of twelve is 8,000 px wide and unusable; a
#: journal column takes three or four across.
PER_ROW = 4


def compose(figure_id: str, panels: list, out: str, rows: list | None = None) -> str:
    """One figure from its panels, wrapped into a GRID at their true millimetre size.

    Panels are never rescaled -- scaling changes the font size, which is the one thing sizing in
    ems exists to prevent -- so a row is as tall as its tallest panel and shorter panels are
    centred in it. Rows are as wide as the widest row, and a short row is centred.
    """
    images = [mpimg.imread(p["png"]) for p in panels]
    # `scale` is for a panel drawn to ANOTHER paper's type size: the bradyzoite paper sets 6 pt,
    # so its panels are enlarged by 8/6 to bring their text to this paper's size. Scaling a panel
    # drawn here would change its font, which is what the rule against rescaling exists to stop.
    scales = [float(p.get("scale", 1.0)) for p in panels]
    widths = [img.shape[1] / DPMM * k for img, k in zip(images, scales)]
    heights = [img.shape[0] / DPMM * k for img, k in zip(images, scales)]

    if rows:
        # Explicit row sizes, e.g. [1, 2, 3]: a full-width schematic on a row of its own.
        sizes, rows, i = rows, [], 0
        for n in sizes:
            rows.append(list(range(i, min(i + n, len(images)))))
            i += n
    else:
        rows = [list(range(i, min(i + PER_ROW, len(images))))
                for i in range(0, len(images), PER_ROW)]
    row_w = [sum(widths[i] for i in r) + GUTTER_MM * (len(r) - 1) for r in rows]
    row_h = [max(heights[i] for i in r) + LETTER_STRIP_MM for r in rows]
    total_w = max(row_w)
    total_h = sum(row_h) + GUTTER_MM * (len(rows) - 1)

    fig = plt.figure(figsize=(total_w / 25.4, total_h / 25.4))
    y_top = total_h
    for r, rw, rh in zip(rows, row_w, row_h):
        x = (total_w - rw) / 2.0            # a short row is centred, not left-ragged
        for i in r:
            w, h = widths[i], heights[i]
            y = y_top - LETTER_STRIP_MM - h
            ax = fig.add_axes([x / total_w, y / total_h, w / total_w, h / total_h])
            ax.imshow(images[i])
            ax.set_axis_off()
            # The letter comes from the FIGURE spec, stamped here, top left of the panel cell.
            letter = panels[i].get("letter")
            if letter:
                fig.text(x / total_w, (y_top) / total_h, letter.upper(),
                         fontsize=FS.LETTER_PT, fontweight="bold", color=FS.INK,
                         ha="left", va="top")
            x += w + GUTTER_MM
        y_top -= rh + GUTTER_MM
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
