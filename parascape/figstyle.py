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

#: THE BASE FONT, and the unit every other size is derived from. Einar, 2026-11-19: "always use
#: open sans, and in the figure graphs should be sized relative to the font size of 8". So nothing
#: below is a millimetre chosen by eye -- every length is a multiple of the em, and changing
#: FONT_PT rescales the whole figure without changing a single proportion.
FONT_PT = 8.0
EM_MM = FONT_PT * 25.4 / 72.0            # one em, in millimetres

#: Text sizes, all relative to the base. Minimising the text on a graph is a house rule, so there
#: are only three: the base for axis labels and category names, a smaller one for tick numbers and
#: in-plot notes, and the panel letter.
TICK_PT = FONT_PT * 0.80
NOTE_PT = FONT_PT * 0.80
LETTER_PT = FONT_PT * 1.15

#: One category occupies this much of the CATEGORY axis, of which the bar takes BAR_FRACTION.
#: Fixing the SLOT rather than the bar is what makes two panels with three conditions identical.
SLOT_EM = 3.6
SLOT_MM = SLOT_EM * EM_MM
BAR_FRACTION = 0.62

#: A category laid out ALONG THE X AXIS carries its name written underneath, so it needs more room
#: than one laid out down the y axis. Names longer than this wrap or tilt rather than widen the
#: slot, because widening it would break the rule that two panels of the same category count are
#: the same size.
# 3.2 em, not 5.4: the names are TILTED under a vertical bar, so they no longer need a slot wide
# enough to hold them flat. At 5.4 a fourteen-category panel came out 213 mm wide beside a
# one-category panel of 40 mm, which is the opposite of "white space should be minimized".
SLOT_EM_VERTICAL = 3.2
SLOT_MM_VERTICAL = SLOT_EM_VERTICAL * EM_MM

#: A panel with many conditions narrows its slot AND its bars (Einar, 2026-11-19: "if a graph has
#: that many conditions the bars and bar spasing needs to be smaller"). Above CROWDED_ABOVE
#: categories the slot shrinks towards SLOT_EM_MIN, so a fourteen-bar panel is not simply a
#: three-bar panel repeated until it runs off the page. The bar keeps the same FRACTION of its
#: slot, so bars and gaps shrink together and the panel stays symmetrical -- what changes is the
#: pitch, not the proportion.
CROWDED_ABOVE = 6
SLOT_EM_MIN = 1.9


def slot_em(n_categories: int, vertical: bool = True) -> float:
    """The slot one category gets, in ems, narrowing once a panel is crowded."""
    base = SLOT_EM_VERTICAL if vertical else SLOT_EM
    if n_categories <= CROWDED_ABOVE:
        return base
    # Linear from the base at CROWDED_ABOVE down to the minimum at twice that, then flat: a panel
    # with 30 categories is already as tight as it can be read at 8 pt.
    span = max(CROWDED_ABOVE, 1)
    t = min(1.0, (n_categories - CROWDED_ABOVE) / span)
    return base + (SLOT_EM_MIN - base) * t

#: Length of the VALUE axis, in ems. The same in every panel, so two panels differ only in how
#: many categories they have.
VALUE_EM = 20.0

#: Above this many points per group, a bar with whiskers is refused.
MAX_POINTS_FOR_A_BAR = 20

#: Panel outside the plotting area, in ems: left, bottom, right, top. Whitespace is minimised, so
#: these are the room the tick numbers and the axis title actually need and nothing more. The left
#: margin is the one that has to hold a rotated axis title plus a number.
MARGIN_EM = (4.2, 3.0, 0.4, 0.8)
MARGIN_MM = tuple(e * EM_MM for e in MARGIN_EM)


def use() -> None:
    """Set the rcParams of the house style. Call once before drawing.

    Open Sans is required, not preferred: a figure that silently falls back to DejaVu is a figure
    that does not match the others. If it is missing this says so rather than drawing it wrong.
    """
    if "Open Sans" in {f.name for f in font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = "Open Sans"
    else:
        raise RuntimeError("Open Sans is not installed: the house style requires it")
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.size": FONT_PT, "axes.linewidth": 0.5,
        "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
        "axes.labelsize": FONT_PT, "legend.fontsize": NOTE_PT,
        "xtick.major.size": 2.0, "ytick.major.size": 2.0,
        "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
        "axes.labelpad": 1.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": False, "axes.edgecolor": INK,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
    })


def panel_width_mm(n_categories: int, vertical: bool = True) -> float:
    """How wide a VERTICAL panel showing `n_categories` must be, so two such panels match.

    Only meaningful for a vertical panel, where the categories run along x. A horizontal panel
    takes its width from the value axis and its HEIGHT from the categories; `panel_height_mm`
    answers that one.
    """
    left, _bottom, right, _top = MARGIN_MM
    return left + right + slot_em(n_categories, vertical) * EM_MM * n_categories


def panel_height_mm(n_categories: int) -> float:
    """How tall a HORIZONTAL (`barh`) panel showing `n_categories` must be."""
    _left, bottom, _right, top = MARGIN_MM
    return bottom + top + slot_em(n_categories, False) * EM_MM * n_categories


#: Length of the value axis, in millimetres, derived from the em like everything else.
VALUE_MM = VALUE_EM * EM_MM


def figure(n_categories: int, height_mm: float = 52.0, vertical: bool = False,
           value_mm: float = VALUE_MM):
    """A figure and axes sized by the number of categories, with the house margins.

    The plotting area is exactly one slot per category ALONG THE CATEGORY AXIS, so a
    three-condition panel drawn here is the same size as any other three-condition panel drawn
    here, which is the symmetry Einar asked for and the thing eyeballed margins never deliver.

    Which axis that is depends on the orientation, and getting this wrong is what made the first
    set of panels inconsistent. `vertical=True` lays the categories along x (violins, boxes,
    vertical bars), so the slots set the WIDTH and `height_mm` sets the value axis.
    `vertical=False` lays them down y (`barh`), so the slots set the HEIGHT and `value_mm` sets
    the value axis. The old code put the slots on the width in BOTH cases, so every horizontal
    panel took its bar thickness from a hand-picked height -- 48 mm here, 60 mm there -- and no
    two of them had bars of the same thickness, which is rule 1.
    """
    left, bottom, right, top = MARGIN_MM
    if vertical:
        width_mm = left + right + slot_em(n_categories, True) * EM_MM * n_categories
    else:
        width_mm = left + right + value_mm
        height_mm = bottom + top + slot_em(n_categories, False) * EM_MM * n_categories
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


#: At or below this many points in a group, the individual points are drawn over the summary, as
#: the reference figures draw replicates over their bars. Above it, a violin IS the points and
#: overplotting 1,200 markers only hides the shape it exists to show.
MAX_POINTS_DRAWN_OVER = 60


def bars(ax, values, names, control=(), highlight=(), vertical=True, rotate=None):
    """Bars in the house geometry, VERTICAL by default.

    Einar, 2026-11-19: "horizontal bar graphs are ok to use in rare occasion when they are the best
    way to show the data. but you use them way to much ... exchange them." Nine panels of the paper
    were horizontal because a long category name fits beside a bar without thought; that is a
    reason to shorten the name, not to turn the graph on its side. A horizontal bar is for when the
    category axis genuinely cannot be read any other way.

    `control` names the bars drawn grey, `highlight` the ones drawn in the first series colour; a
    name in neither is grey too, since grey is the house reference. `rotate` tilts the category
    names, and defaults to tilting only when the names are too long to sit flat.
    """
    positions = list(range(len(values)))
    if highlight:
        fills = [SERIES[0] if n in set(highlight) else CONTROL for n in names]
    else:
        fills = colours(names, control=control)
    if vertical:
        ax.bar(positions, values, width=bar_width(), linewidth=0, color=fills, alpha=ALPHA)
        ax.set_xticks(positions)
        longest = max((len(str(n)) for n in names), default=0)
        tilt = (30 if longest > 10 else 0) if rotate is None else rotate
        ax.set_xticklabels(names, rotation=tilt,
                           ha="right" if tilt else "center",
                           rotation_mode="anchor" if tilt else None)
        ax.tick_params(axis="x", length=0)
    else:
        ax.barh(positions, values, height=bar_width(), linewidth=0, color=fills, alpha=ALPHA)
        ax.set_yticks(positions)
        ax.set_yticklabels(names)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    return fills


def violin(ax, groups, names, control=(), vertical=True, draw_points=None):
    """A violin per group, at exactly the width and spacing a bar would have had.

    This exists because the rule "more than twenty points is never a bar" was written down and then
    not followed: three panels of the paper summarised distributions of hundreds of genes as one
    bar each. A violin here is not a different style of bar -- it occupies the same slot, takes
    `bar_width()` like everything else, and carries the same opacity, so a violin panel and a bar
    panel of the same category count are still dimensionally identical.

    `groups` is one array of values per group; `names` labels them; `control` names the groups
    drawn grey. Individual points are drawn over a group with at most `MAX_POINTS_DRAWN_OVER`
    values unless `draw_points` says otherwise.
    """
    fills = colours(names, control=control)
    positions = list(range(len(groups)))
    # `orientation=` since matplotlib 3.11; `vert=` is deprecated and goes in 3.13.
    parts = ax.violinplot([list(g) for g in groups], positions=positions,
                          widths=bar_width(), showextrema=False, showmedians=False,
                          orientation="vertical" if vertical else "horizontal")
    for body, fill in zip(parts["bodies"], fills):
        body.set_facecolor(fill)
        body.set_alpha(ALPHA)
        body.set_linewidth(0)
    # The median as a short rule, in ink: a violin with no median cannot be read against another.
    for position, values, in zip(positions, groups):
        values = [v for v in values if v == v]
        if not values:
            continue
        median = sorted(values)[len(values) // 2]
        span = bar_width() * 0.42
        if vertical:
            ax.plot([position - span, position + span], [median, median],
                    color=INK, linewidth=0.9, zorder=4, solid_capstyle="butt")
        else:
            ax.plot([median, median], [position - span, position + span],
                    color=INK, linewidth=0.9, zorder=4, solid_capstyle="butt")
    for position, values, fill in zip(positions, groups, fills):
        values = [v for v in values if v == v]
        show = (len(values) <= MAX_POINTS_DRAWN_OVER) if draw_points is None else draw_points
        if not show or not values:
            continue
        jitter = _spread(len(values), bar_width() * 0.30)
        at = [position + j for j in jitter]
        xs, ys = (at, values) if vertical else (values, at)
        ax.plot(xs, ys, "o", markersize=1.9, markerfacecolor=fill, markeredgecolor=INK,
                markeredgewidth=0.25, alpha=ALPHA, linestyle="none", zorder=5)
    return parts


def _spread(n: int, half_width: float) -> list:
    """Deterministic offsets across a slot, so a redraw puts the points back where they were."""
    if n == 1:
        return [0.0]
    step = (2.0 * half_width) / (n - 1)
    return [-half_width + i * step for i in range(n)]


def placeholder(what: str, why: str = "", letter_text: str = "a", n_categories: int = 6,
                height_mm: float = 54.0):
    """A panel that is NOT a result: a marked, hatched box standing in for one not yet measured.

    Einar, 2026-11-19: "we will make the paper figure with real data and with placeholders. then we
    will fill in the data for the paper." So a placeholder has to be assemblable into the page at
    the right size, and it has to be impossible to mistake for a finding -- a blank panel in a
    draft becomes a figure nobody remembers is empty.

    It takes the same geometry as the panel that will replace it, so the page does not reflow when
    the real one arrives, and it carries hatching plus the word PLACEHOLDER because the reference
    figures mark an absent cell with hatching rather than leaving it blank, which is the same rule
    as reporting `n.m.` instead of nothing.
    """
    use()
    fig, ax = figure(n_categories, height_mm=height_mm)
    ax.add_patch(plt.Rectangle((0.02, 0.02), 0.96, 0.96, transform=ax.transAxes,
                               facecolor=SECOND_GREY, edgecolor=CONTROL, linewidth=0.6,
                               hatch="////", alpha=ALPHA, zorder=0))
    ax.text(0.5, 0.60, "PLACEHOLDER", transform=ax.transAxes, fontsize=FONT_PT, fontweight="bold",
            color=DARK_GREY, ha="center", va="center")
    ax.text(0.5, 0.44, what, transform=ax.transAxes, fontsize=NOTE_PT, color=INK,
            ha="center", va="center", wrap=True)
    if why:
        ax.text(0.5, 0.26, why, transform=ax.transAxes, fontsize=NOTE_PT, color=DARK_GREY,
                ha="center", va="center", style="italic")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    letter(fig, letter_text)
    return fig


def letter(fig, text: str) -> None:
    """The panel letter: bold, top-left, outside the plotting area."""
    # ABOVE the axes, not at the figure corner: a long rotated y-axis label reaches the top of the
    # figure, and the letter landed on the word "way".
    ax = fig.axes[0] if fig.axes else None
    if ax is not None:
        box = ax.get_position()
        fig.text(max(box.x0 - 0.10, 0.002), min(box.y1 + 0.085, 0.995), text,
                 fontsize=LETTER_PT, fontweight="bold", color=INK, ha="left", va="top")
    else:
        fig.text(0.002, 0.995, text, fontsize=LETTER_PT, fontweight="bold", color=INK,
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


def place_labels(ax, fig, points, labels, colours, fontsize=None):
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
                               fontsize=fontsize or NOTE_PT, color=colour, ha=ha, va=va,
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
