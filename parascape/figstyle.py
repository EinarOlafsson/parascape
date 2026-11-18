"""The house figure style: the palette, the geometry and the rules that keep panels symmetrical.

One module so that every panel of every paper is drawn by the same code rather than by the same
intentions. Agreed with Einar on 2026-11-18 and copied from the reference figures named in the
`publication-figures` skill (Nature Communications s41467-026-73079-z figures 2 and 4; Nature
Microbiology s41564-024-01754-2 figures 1 and 3).

The rules this exists to enforce, because every one of them has been got wrong by hand:

* bars, violins and boxes take the SAME width and the SAME gap, in every graph of a figure;
* a panel with three conditions is the same width as any other panel with three conditions;
* the control or reference is grey, never a colour;
* elements are drawn at 60 per cent opacity;
* more than 20 points per group is a violin or a box, never a bar with whiskers.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib import font_manager

#: Categorical series, assigned in this order and never cycled.
SERIES = ("#B02418", "#E8912A", "#35A7DE", "#2B4E9B",
          "#1E7D3C", "#6B3F9E", "#C2185B", "#7B2D26")

#: Structural greys. The control is CONTROL and nothing else.
CONTROL = "#9E9E9E"
SECOND_GREY = "#CFCFCF"
DARK_GREY = "#4D4D4D"
INK = "#1A1A1A"

#: Opacity of every drawn element.
ALPHA = 0.6

#: One category occupies this much width, of which the bar takes BAR_FRACTION.
#: Fixing the SLOT rather than the bar is what makes two panels with three conditions identical.
SLOT_MM = 7.0
BAR_FRACTION = 0.62

#: A category laid out ALONG THE X AXIS needs more room than one laid out down the y axis,
#: because its name is written horizontally under it rather than beside it. Seven millimetres is
#: ample for a horizontal bar and far too narrow for a box with a two-word label under it, which
#: is how the first draft of figure 4 came out with its three labels written through each other.
#: The symmetry rule still holds, within an orientation: every vertical panel uses this width per
#: category and every horizontal one uses SLOT_MM.
SLOT_MM_VERTICAL = 17.0

#: Above this many points per group, a bar with whiskers is refused.
MAX_POINTS_FOR_A_BAR = 20

#: Millimetres of panel outside the plotting area, left and bottom, for labels.
MARGIN_MM = (19.0, 13.0, 3.0, 4.0)      # left, bottom, right, top


def use() -> None:
    """Set the rcParams of the house style. Call once before drawing."""
    if "Open Sans" in {f.name for f in font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = "Open Sans"
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.size": 6.2, "axes.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": False, "axes.edgecolor": INK,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
    })


def panel_width_mm(n_categories: int, vertical: bool = False) -> float:
    """How wide a panel showing `n_categories` must be, so two such panels match exactly."""
    left, _bottom, right, _top = MARGIN_MM
    return left + right + (SLOT_MM_VERTICAL if vertical else SLOT_MM) * n_categories


def figure(n_categories: int, height_mm: float = 52.0, vertical: bool = False):
    """A figure and axes sized by the number of categories, with the house margins.

    The plotting area is exactly one slot per category, so a three-condition panel drawn here is
    the same width as any other three-condition panel drawn here, which is the symmetry Einar
    asked for and the thing eyeballed margins never deliver.
    """
    left, bottom, right, top = MARGIN_MM
    width_mm = panel_width_mm(n_categories, vertical=vertical)
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax = fig.add_axes([left / width_mm, bottom / height_mm,
                       (width_mm - left - right) / width_mm,
                       (height_mm - bottom - top) / height_mm])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return fig, ax


def bar_width() -> float:
    """The width of one bar, in category units, so bars and violins match."""
    return BAR_FRACTION


def colours(names, control=()) -> list:
    """One colour per category: grey for a control, then the series palette in order."""
    control = {str(c) for c in control}
    out, i = [], 0
    for name in names:
        if str(name) in control:
            out.append(CONTROL)
        else:
            out.append(SERIES[i % len(SERIES)])
            i += 1
    return out


def refuse_a_bar(n_points: int) -> bool:
    """Whether this many points per group must be a violin or a box rather than a bar."""
    return n_points > MAX_POINTS_FOR_A_BAR


def letter(fig, text: str) -> None:
    """The panel letter: bold, top-left, outside the plotting area."""
    fig.text(0.012, 0.985, text, fontsize=7.5, fontweight="bold", color=INK,
             ha="left", va="top")


#: Where a label may sit relative to its point, in points, tried in this order.
LABEL_CANDIDATES = ((4.5, 0.0, "left", "center"), (-4.5, 0.0, "right", "center"),
                    (4.5, 5.0, "left", "bottom"), (-4.5, 5.0, "right", "bottom"),
                    (4.5, -5.5, "left", "top"), (-4.5, -5.5, "right", "top"),
                    (0.0, 6.5, "center", "bottom"), (0.0, -7.0, "center", "top"))


def _text_box(text, renderer):
    """The bounding box of a label, measured on the TEXT and nothing else.

    `Annotation.get_window_extent` includes the leader line, which measures a box several times
    too large and makes every candidate position look occupied.
    """
    from matplotlib import text as mtext
    return mtext.Text.get_window_extent(text, renderer=renderer)


def place_labels(ax, fig, points, labels, colours, fontsize=4.9):
    """Label every point in the first candidate position that hits nothing already placed.

    Guessed offsets collide as soon as the data moves: figure 7 put `struct` under
    `all layers, weighted` and ran three interaction layers through each other. Measuring is cheap
    at a few dozen labels and stays correct when the numbers change.
    """
    renderer = fig.canvas.get_renderer()
    taken = [ax.transData.transform((x, y)) for x, y in points]
    placed = []
    width, height = fig.get_size_inches() * fig.dpi
    for (x, y), name, colour in zip(points, labels, colours):
        for dx, dy, ha, va in LABEL_CANDIDATES:
            text = ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                               fontsize=fontsize, color=colour, ha=ha, va=va,
                               annotation_clip=False)
            fig.canvas.draw()
            box = _text_box(text, renderer).expanded(1.06, 1.16)
            clash = (box.x0 < 2 or box.y0 < 2 or box.x1 > width - 2 or box.y1 > height - 2
                     or any(box.overlaps(other) for other in placed)
                     or any(box.x0 - 1 <= px <= box.x1 + 1 and box.y0 - 1 <= py <= box.y1 + 1
                            for px, py in taken))
            if not clash:
                placed.append(box)
                break
            text.remove()
