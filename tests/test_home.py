"""The home screen: what it offers, and what it admits it cannot do.

This replaces the Gene/List/Label node model as the front door. The tests that
matter here are not about layout -- they are about honesty. A home screen that
hides what is missing teaches the user the feature never existed; one that greys
it out with the reason teaches them what to go and get.
"""
import pytest

from parascape import modules as MOD


# --------------------------------------------------------------------- catalogue


def test_every_module_belongs_to_a_declared_band():
    declared = {title for title, _ in MOD.BANDS}
    for m in MOD.MODULES:
        assert m.band in declared, f"{m.key} sits in an undeclared band {m.band!r}"


def test_the_bands_answer_what_the_person_came_to_do():
    """Bands are tasks, not data shapes. The node model named them after the
    matrix -- a row, a block of rows, a column -- which is correct about the data
    and useless as a front door: nobody arrives wanting a column."""
    assert [t for t, _ in MOD.BANDS] == [
        "Look up", "Explore", "Ask", "Bring your own", "Data"]
    for _, note in MOD.BANDS:
        assert note.endswith("."), "a band says what you have, in a sentence"


def test_every_module_says_what_you_get_not_what_it_is():
    for m in MOD.MODULES:
        assert m.note, f"{m.key} has no note; a one-word tile is unidentifiable"
        assert len(m.name.split()) <= 2, f"{m.name!r} is not a one or two word title"
        assert not m.note[0].isupper() or m.note.split()[0] in ("CSV", "ParaScape"), \
            f"{m.key}: the note is a continuation of the name, not a new sentence"


def test_module_keys_are_unique():
    keys = [m.key for m in MOD.MODULES]
    assert len(keys) == len(set(keys))


def test_get_raises_rather_than_returning_none():
    """A missing module must fail where it is asked for, not later where it is used."""
    with pytest.raises(KeyError):
        MOD.get("not_a_module")
    assert MOD.get("gene").name == "Gene"


# ------------------------------------------------------------------ availability


def test_fill_is_unavailable_by_computation_not_by_hand():
    """The one that must never be hardcoded. Fill opens when a question has cleared
    its bar, and whether one has is a measurement that changes as evidence arrives,
    so a stored 'coming soon' would be wrong in both directions."""
    fill = MOD.get("fill")
    assert MOD.unavailable_reason(fill, fillable=0)
    assert "90%" in MOD.unavailable_reason(fill, fillable=0)
    assert MOD.unavailable_reason(fill, fillable=1) == "", \
        "a question clearing the bar must open Fill without a code change"


def test_a_module_with_an_unbuilt_view_says_so():
    for key in ("runs", "list"):
        reason = MOD.unavailable_reason(MOD.get(key))
        assert reason and reason.startswith("needs "), f"{key}: {reason!r}"


def test_a_module_whose_dock_is_absent_names_the_dock():
    reason = MOD.unavailable_reason(MOD.get("gene"), docks=("maps", "jobs"))
    assert "evidence" in reason


def test_a_module_whose_dock_is_present_is_available():
    assert MOD.unavailable_reason(MOD.get("gene"), docks=("evidence",)) == ""


# ------------------------------------------------------------------------ search


def test_search_ranks_a_named_module_above_a_mentioned_one():
    hits = MOD.search("map")
    assert hits and hits[0].name == "Map", [h.name for h in hits]


def test_search_is_empty_for_nothing():
    assert MOD.search("") == ()
    assert MOD.search("   ") == ()


# -------------------------------------------------------------------- the widget


@pytest.fixture
def home(qtbot):
    from parascape.home import HomeScreen
    widget = HomeScreen(counts={"T. gondii genes": ("8,140", "ME49")})
    qtbot.addWidget(widget)
    return widget


def test_every_module_gets_a_tile(home):
    assert len(home._tiles) == len(MOD.MODULES)


def test_an_unavailable_tile_is_greyed_with_its_reason_and_never_hidden(home):
    missing = home.unavailable()
    assert missing, "something is always unfinished; a home screen claiming otherwise is lying"
    for key, reason in missing.items():
        tile = home.tile(key)
        assert tile.isVisible() or not tile.isHidden(), f"{key} is hidden rather than greyed"
        assert not tile.isEnabled()
        assert reason and reason in tile.toolTip()


def test_pressing_an_available_tile_reports_its_module(home):
    seen = []
    home.module_chosen.connect(seen.append)
    assert home.open_module("gene") is True
    assert seen == ["gene"]


def test_an_unavailable_tile_cannot_be_opened(home):
    seen = []
    home.module_chosen.connect(seen.append)
    assert home.open_module("fill") is False
    assert seen == [], "a greyed tile must not fire, or the app opens an empty view"


def test_the_window_can_override_availability(qtbot):
    """The window knows which docks it actually built; the catalogue only knows what
    exists in principle."""
    from parascape.home import HomeScreen
    widget = HomeScreen(reasons={"map": "no OpenGL context on this machine"})
    qtbot.addWidget(widget)
    assert widget.unavailable()["map"] == "no OpenGL context on this machine"
    assert not widget.tile("map").isEnabled()


def test_the_umap_and_the_cell_drawing_are_modules_not_the_landing_page():
    """Einar's correction: they are the best thing in the application and they are
    not a front door. They moved into Explore rather than being discarded."""
    assert MOD.get("map").band == "Explore"
    assert MOD.get("cell").band == "Explore"
    assert MOD.get("map").opens and MOD.get("cell").opens


def test_the_page_is_not_painted_white(home):
    """Caught by rendering the page to an image rather than reading the code: a
    QScrollArea paints its viewport from the Base role, which is near-white, so the
    band column came out white and every heading on it vanished."""
    from parascape import theme as TH
    sheet = home.styleSheet()
    assert "HomeScrollViewport" in sheet and "HomeBands" in sheet
    assert TH.palette_for("dark")["bg"] in sheet
