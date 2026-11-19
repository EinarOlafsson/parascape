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

import os

import matplotlib.pyplot as plt
from matplotlib import font_manager

#: The reference-paper palette, measured from the figure pixels by the maintainer and recorded in
#: the `spacr-graphs` skill. These REPLACE the set this module used to carry, which was a different
#: palette I had assembled from two of the papers.
ACCENT = "#0A77B8"          # primary accent: the thing the reader should look at
TEAL = "#0A9878"            # second series
NAVY = "#283868"            # third series
ACCENT_LIGHT = "#78AED2"    # paired lighter tint of the accent
MAGENTA = "#C878A8"         # fourth series
AMBER = "#F0A830"           # fifth series only; never darkened towards brown
SERIES = (ACCENT, TEAL, NAVY, ACCENT_LIGHT, MAGENTA, AMBER)

#: Structural greys. The control is CONTROL and nothing else; background points are lighter.
CONTROL = "#8C8C8C"
BACKGROUND = "#C8C8C8"
SECOND_GREY = "#C8C8C8"
DARK_GREY = "#4D4D4D"
INK = "#000000"

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
TICK_PT = FONT_PT          # 8 pt minimum for EVERY text element
NOTE_PT = FONT_PT          # 8 pt minimum for EVERY text element
LETTER_PT = 10.0           # bold CAPITAL, 10 pt, top left, outside the axes

#: One condition = one slot of SLOT_IN = 0.25 in; the bar takes 0.6 of a slot and the gap 0.4.
#: These are the maintainer's numbers, not a proportion chosen here.
SLOT_IN = 0.25

#: No categorical panel is narrower than this many slots: at 0.25 in a one-condition violin came
#: out 6 mm wide, which is a sliver and not a graph.
MIN_SLOTS = 4
SLOT_MM = SLOT_IN * 25.4
SLOT_EM_VERTICAL = SLOT_MM / EM_MM
SLOT_MM_VERTICAL = SLOT_MM
BAR_FRACTION = 0.6

#: A crowded panel still uses the fixed slot: the slot IS the rule, so a 14-condition panel is
#: 14 x 0.25 in wide and nothing narrows. (2026-11-19: my earlier shrinking rule is withdrawn --
#: the maintainer's spec fixes the slot, and "smaller bars when crowded" is satisfied by the slot
#: being 0.25 in rather than the 5.4 em I had invented.)
CROWDED_ABOVE = 10 ** 9
SLOT_EM_MIN = SLOT_EM_VERTICAL


def slot_em(n_categories: int, vertical: bool = True) -> float:
    """The slot one condition gets, in ems. Fixed, by the maintainer's spec."""
    return SLOT_MM / EM_MM


#: Length of the value axis. A categorical panel is as small as 8 pt text allows.
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
        "font.size": FONT_PT,
        "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
        "axes.labelsize": FONT_PT, "legend.fontsize": NOTE_PT,
        "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
        "axes.labelpad": 2.0,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": False, "axes.edgecolor": INK,
        "axes.spines.top": True, "axes.spines.right": True,    # FULL BOX FRAME (4/5 papers)
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
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
        # A panel never narrower than MIN_SLOTS slots: at 0.25 in per condition a one-condition
        # violin came out 6 mm wide, which is a sliver and not a graph.
        slots = max(n_categories, MIN_SLOTS)
        width_mm = left + right + slot_em(n_categories, True) * EM_MM * slots
    else:
        width_mm = left + right + value_mm
        height_mm = bottom + top + slot_em(n_categories, False) * EM_MM * n_categories
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax = fig.add_axes([left / width_mm, bottom / height_mm,
                       (width_mm - left - right) / width_mm,
                       (height_mm - bottom - top) / height_mm])
    return fig, ax


def square_figure(side_mm: float = 46.0):
    """A figure whose PLOTTING AREA is square, for every continuous-axis plot.

    The skill: "Continuous-axis plots (histograms, density plots, scatter plots, line graphs,
    sweeps, curves) are SQUARE by default: plot area width = height." A scatter drawn in a
    categorical panel's shape is wide for no reason and wastes the page.
    """
    left, bottom, right, top = MARGIN_MM
    width_mm = left + right + side_mm
    height_mm = bottom + top + side_mm
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    ax = fig.add_axes([left / width_mm, bottom / height_mm,
                       side_mm / width_mm, side_mm / height_mm])
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
    # One accent per comparison, everything else grey (skill, 5/5 reference papers). With no
    # highlight named, every bar is grey: the palette in order made five colours mean nothing.
    fills = [ACCENT if n in set(highlight) else CONTROL for n in names]
    if vertical:
        ax.bar(positions, values, width=bar_width(), linewidth=0, color=fills, alpha=ALPHA)
        ax.set_xticks(positions)
        longest = max((len(str(n)) for n in names), default=0)
        # Tilt whenever a name is wider than its slot. A 0.25 in slot holds about 4 characters at
        # 8 pt, so "Toxoplasma" written flat ran into "Plasmodium" beside it.
        tilt = (30 if longest > 4 else 0) if rotate is None else rotate
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
    # ONE accent colour, everything else grey -- the same rule as `bars`. Cycling the palette made
    # four placeholder violins four colours meaning nothing. Every name NOT listed as a control
    # takes the single accent; with no control named, all are grey.
    control = {str(c) for c in control}
    fills = [ACCENT if (control and str(n) not in control) else CONTROL for n in names]
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


def legend(ax, loc: str = "best", **kw):
    """The house legend: no frame, ink text, tight. The ONLY place text other than an axis title
    is allowed to sit, per Einar 2026-11-19 ("there should only be text on x and y axes or in
    legends"). A reference line is identified here, never by a label floating beside it.
    """
    leg = ax.legend(frameon=False, loc=loc, handlelength=1.0, borderaxespad=0.2,
                    labelspacing=0.25, fontsize=NOTE_PT, **kw)
    for text in leg.get_texts():
        text.set_color(INK)
    return leg


def pad_limits(ax, values, axis: str = "y", frac: float = 0.08) -> None:
    """Give the data room so nothing touches the frame.

    Violins drawn to the exact data range sat on the top and bottom spines, which reads as data
    running off the panel. The pad is a fraction of the span, so it scales with the quantity.
    """
    lo, hi = float(min(values)), float(max(values))
    span = (hi - lo) or abs(hi) or 1.0
    if axis == "y":
        ax.set_ylim(lo - span * frac, hi + span * frac)
    else:
        ax.set_xlim(lo - span * frac, hi + span * frac)


def letter(fig, text: str) -> None:
    """The panel letter: bold, top-left, outside the plotting area."""
    # ABOVE the axes, not at the figure corner: a long rotated y-axis label reaches the top of the
    # figure, and the letter landed on the word "way".
    if not text:
        return                       # a figure with ONE panel has no panel letter
    # Panels drawn for a multi-panel FIGURE carry no letter of their own: the assembler stamps
    # the letter the figure spec assigns. A letter baked into a reused panel is the letter it had
    # in some earlier figure, which is how one figure came out lettered A B C D B B.
    if os.environ.get("PARASCAPE_NO_PANEL_LETTERS"):
        return
    # Top left of the FIGURE, outside the axes: placing it relative to the axes box put it on the
    # top tick label of a narrow panel.
    fig.text(0.0, 1.0, text.upper(), fontsize=LETTER_PT, fontweight="bold", color=INK,
             ha="left", va="top")


#: Where a label may sit relative to its point, in points, tried in this order.
LABEL_CANDIDATES = ((4.5, 0.0, "left", "center"), (-4.5, 0.0, "right", "center"),
                    (4.5, 5.0, "left", "bottom"), (-4.5, 5.0, "right", "bottom"),
                    (4.5, -5.5, "left", "top"), (-4.5, -5.5, "right", "top"),
                    (0.0, 6.5, "center", "bottom"), (0.0, -7.0, "center", "top"),
                    (9.0, 11.0, "left", "bottom"), (-9.0, 11.0, "right", "bottom"),
                    (9.0, -11.0, "left", "top"), (-9.0, -11.0, "right", "top"))


def _text_box(text, renderer):
    """The bounding box of a label, measured on the TEXT and nothing else.

    `Annotation.get_window_extent` includes the leader line, which measures a box several times
    too large and makes every candidate position look occupied.
    """
    from matplotlib import text as mtext
    return mtext.Text.get_window_extent(text, renderer=renderer)


def _leader(dy, colour):
    """A thin line from a point to a label set above or below it. Beside the point (dy == 0) a
    label is unambiguous; above or below, with neighbours close, it read as naming the wrong
    point (figure 6D: 'rhoptries' looked like the micronemes point)."""
    if dy == 0:
        return None
    return dict(arrowstyle="-", lw=0.4, color=colour, shrinkA=0, shrinkB=2.5)


def place_labels(ax, fig, points, labels, colours, fontsize=None):
    """Label every point in the first candidate position that hits nothing already placed.

    Guessed offsets collide as soon as the data moves: figure 7 put `struct` under
    `all layers, weighted` and ran three interaction layers through each other. Measuring is cheap
    at a few dozen labels and stays correct when the numbers change.
    """
    renderer = fig.canvas.get_renderer()
    taken = [ax.transData.transform((x, y)) for x, y in points]
    placed = []
    # Labels stay INSIDE the axes: bounding them by the figure let 'dense granules' run across the
    # left spine (figure 6D). A label that fits nowhere goes in its least-crowded candidate rather
    # than being dropped -- 'rhoptries' silently vanished from the same panel.
    frame = ax.get_window_extent(renderer)
    for (x, y), name, colour in zip(points, labels, colours):
        fallback = None
        for dx, dy, ha, va in LABEL_CANDIDATES:
            text = ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                               fontsize=fontsize or NOTE_PT, color=colour, ha=ha, va=va,
                               annotation_clip=False, arrowprops=_leader(dy, colour))
            fig.canvas.draw()
            box = _text_box(text, renderer).expanded(1.06, 1.16)
            outside = (box.x0 < frame.x0 + 2 or box.y0 < frame.y0 + 2
                       or box.x1 > frame.x1 - 2 or box.y1 > frame.y1 - 2)
            hits = (sum(box.overlaps(other) for other in placed)
                    + sum(box.x0 - 1 <= px <= box.x1 + 1 and box.y0 - 1 <= py <= box.y1 + 1
                          for px, py in taken))
            if not outside and not hits:
                placed.append(box)
                break
            if not outside and (fallback is None or hits < fallback[0]):
                fallback = (hits, dx, dy, ha, va)
            text.remove()
        else:
            if fallback is not None:
                _, dx, dy, ha, va = fallback
                text = ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                                   fontsize=fontsize or NOTE_PT, color=colour, ha=ha, va=va,
                                   annotation_clip=False, arrowprops=_leader(dy, colour))
                fig.canvas.draw()
                placed.append(_text_box(text, renderer).expanded(1.06, 1.16))
