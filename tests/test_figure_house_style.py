"""The house figure style, enforced rather than intended.

Einar, 2026-11-19: "there should also be a test for figures."

The rules were written into the `publication-figures` skill on 2026-11-18 and then broken within a
day, in ways no reviewer caught by looking: a panel whose bars were the means of 1,253 genes, a
panel whose highlight was on the two best categories rather than the two worst, and -- the one that
mattered -- every horizontal panel taking its bar thickness from a hand-picked height, so no two
`barh` panels in the paper had bars of the same thickness.

So these are the rules as assertions. They are deliberately about GEOMETRY and PALETTE, which are
checkable, rather than about whether a form suits its data, which is a judgement. The one judgement
that IS checkable is the twenty-point rule, and `refuse_a_bar` exists so a drawing script can ask.

What this cannot catch, and what still needs a person: whether the form fits the measurement,
whether the highlight is on the right end, and whether a label collides. Render it and look at it.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from parascape import figstyle as FS

MM = 25.4


@pytest.fixture(autouse=True)
def _close_figures():
    """Every test here makes figures; left open they warn and hold memory."""
    yield
    plt.close("all")


def size_mm(fig) -> tuple:
    """A figure's size in millimetres, the unit every house rule is written in."""
    # Six decimals, not three: every length is now a multiple of the em (FONT_PT * 25.4 / 72),
    # which is not a round number of millimetres, and rounding to 3 made an exact match look wrong.
    width, height = fig.get_size_inches()
    return round(width * MM, 6), round(height * MM, 6)


# ----------------------------------------------------------------- rule 4: equal counts, equal size
@pytest.mark.parametrize("n", [2, 3, 5, 9, 14])
def test_two_horizontal_panels_with_the_same_count_are_identical(n):
    """Rule 4: "graphs with 3 conditions should always be the same width and so on."

    For a horizontal panel the categories run down y, so it is the HEIGHT that must match, and the
    value axis is the same in both. This is the rule the old code broke: it sized the WIDTH by the
    category count in both orientations and let a caller pick the height.
    """
    first, _ = FS.figure(n)
    second, _ = FS.figure(n)
    assert size_mm(first) == size_mm(second)
    assert size_mm(first)[1] == pytest.approx(FS.panel_height_mm(n))


@pytest.mark.parametrize("n", [2, 3, 5, 9])
def test_two_vertical_panels_with_the_same_count_are_identical(n):
    """The same rule for a vertical panel, where the categories run along x."""
    first, _ = FS.figure(n, vertical=True)
    second, _ = FS.figure(n, vertical=True)
    assert size_mm(first) == size_mm(second)
    assert size_mm(first)[0] == pytest.approx(FS.panel_width_mm(n, vertical=True))


def test_a_horizontal_panel_grows_by_exactly_one_slot_per_category():
    """Rule 1: the same gap and the same bar width, ACROSS panels and not only within one.

    If adding a category grew a panel by anything other than one slot, two panels of different
    counts would have bars of different thickness, which is what the paper had.
    """
    heights = [size_mm(FS.figure(n)[0])[1] for n in (3, 4, 5, 6)]
    steps = [b - a for a, b in zip(heights, heights[1:])]
    for step in steps:
        assert step == pytest.approx(FS.SLOT_MM)


def test_the_value_axis_of_a_horizontal_panel_does_not_depend_on_the_category_count():
    """Two panels with different counts still share a value axis, so their bars are comparable."""
    widths = {size_mm(FS.figure(n)[0])[0] for n in (2, 7, 13)}
    assert len(widths) == 1


# ----------------------------------------------------------------- rule 3: violins match the bars
def test_a_violin_takes_the_same_width_as_a_bar():
    """Rule 3: "all boxes, and violins should ... have the same width as the bars"."""
    fig, ax = FS.figure(3, vertical=True)
    groups = [np.random.default_rng(i).normal(size=200) for i in range(3)]
    parts = FS.violin(ax, groups, ["a", "b", "c"], control=["a"])
    for body in parts["bodies"]:
        xs = body.get_paths()[0].vertices[:, 0]
        assert (xs.max() - xs.min()) <= FS.bar_width() + 1e-9


def test_a_violin_is_drawn_at_the_house_opacity_and_the_control_is_grey():
    """Rules 5 and 6: control grey, every element at 60 per cent."""
    fig, ax = FS.figure(3, vertical=True)
    groups = [np.random.default_rng(i).normal(size=50) for i in range(3)]
    parts = FS.violin(ax, groups, ["control", "b", "c"], control=["control"])
    for body in parts["bodies"]:
        assert body.get_alpha() == pytest.approx(FS.ALPHA)
    first = parts["bodies"][0].get_facecolor()
    assert matplotlib.colors.to_hex(first[0][:3]).lower() == FS.CONTROL.lower()


def test_the_series_palette_is_assigned_in_order_and_never_cycled_before_it_is_spent():
    """The palette rule: assign in order, never generate a hue."""
    names = [f"s{i}" for i in range(len(FS.SERIES))]
    assert FS.colours(names) == list(FS.SERIES)


def test_a_control_takes_grey_wherever_it_sits_in_the_order():
    chosen = FS.colours(["a", "wild type", "b"], control=["wild type"])
    assert chosen[1] == FS.CONTROL
    assert chosen[0] == FS.SERIES[0] and chosen[2] == FS.SERIES[1]


# ----------------------------------------------------------------- rule 2: twenty points
@pytest.mark.parametrize("n,refused", [(1, False), (20, False), (21, True), (1253, True)])
def test_more_than_twenty_points_may_not_be_a_bar(n, refused):
    """Rule 2: "if there are more that 20 data points show a violin plot or box plot instead"."""
    assert FS.refuse_a_bar(n) is refused


def test_the_points_drawn_over_a_group_are_placed_deterministically():
    """A redraw must put the points back where they were, or a figure changes without its data."""
    assert FS._spread(7, 0.3) == FS._spread(7, 0.3)
    spread = FS._spread(7, 0.3)
    assert spread[0] == pytest.approx(-0.3) and spread[-1] == pytest.approx(0.3)
    assert spread == sorted(spread)


# ----------------------------------------------------------------- the rest of the house
def test_no_panel_carries_gridlines_and_the_top_and_right_spines_are_off():
    """"no gridlines" is in the list of what the reference figures do NOT do."""
    FS.use()
    _fig, ax = FS.figure(4)
    assert matplotlib.rcParams["axes.grid"] is False
    assert not ax.spines["top"].get_visible()
    assert not ax.spines["right"].get_visible()


def test_text_stays_text_in_the_pdf():
    """`pdf.fonttype = 42`, so Einar can edit the labels when he assembles the figures."""
    FS.use()
    assert matplotlib.rcParams["pdf.fonttype"] == 42
    assert matplotlib.rcParams["ps.fonttype"] == 42


def test_the_opacity_and_the_control_are_the_values_einar_set():
    """Pinned so a later edit to the module cannot quietly move them."""
    assert FS.ALPHA == 0.6
    assert FS.CONTROL == "#9E9E9E"
    assert FS.MAX_POINTS_FOR_A_BAR == 20


def test_every_series_colour_is_a_six_digit_hex_and_none_repeats():
    assert len(set(FS.SERIES)) == len(FS.SERIES)
    for colour in FS.SERIES + (FS.CONTROL, FS.SECOND_GREY, FS.DARK_GREY, FS.INK):
        assert matplotlib.colors.is_color_like(colour)
        assert colour.startswith("#") and len(colour) == 7


def test_a_violin_group_that_is_all_unknown_does_not_raise():
    """A category nobody measured is drawn as nothing, not as an exception or as a zero."""
    fig, ax = FS.figure(2, vertical=True)
    FS.violin(ax, [np.array([np.nan, np.nan]), np.array([1.0, 2.0, 3.0])], ["empty", "full"])
