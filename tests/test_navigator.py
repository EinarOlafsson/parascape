"""The Navigator: one field over the animated field, answering with anything.

Each test below pins something that went wrong while building it, found by
rendering the page rather than reading the code.
"""
import pytest
from PySide6 import QtCore

from parascape import navigator as NV

GENES = [("TGME49_208830", "dense granule protein GRA16"),
         ("TGME49_223920", "rhoptry neck protein RON3"),
         ("TGME49_214320", "facilitative glucose transporter GT1")]


@pytest.fixture
def index():
    from parascape import help_index as HI
    return HI.static_entries()


# ------------------------------------------------------------------ the ranking


def test_a_module_named_outright_comes_first(index):
    assert NV.suggestions("map", index=index)[0] == NV.suggestions("map", index=index)[0]
    first = NV.suggestions("map", index=index)[0]
    assert (first.kind, first.key) == ("module", "map")


def test_a_gene_is_found_by_its_product(index):
    hits = NV.suggestions("GRA16", genes=GENES, index=index)
    assert ("gene", "TGME49_208830") in {(s.kind, s.key) for s in hits}


def test_a_gene_is_found_by_its_accession(index):
    hits = NV.suggestions("TGME49_223920", genes=GENES, index=index)
    assert hits[0].kind == "gene" and hits[0].key == "TGME49_223920"


def test_questions_are_still_searchable(index):
    """Dropping the old Gene/List/Label cards first dropped the question cards with
    them, and "glucose" stopped finding the glucose question."""
    kinds = [s.kind for s in NV.suggestions("glucose", index=index)]
    assert "question" in kinds


def test_methods_and_datasets_come_from_the_same_ranking(index):
    kinds = {s.kind for s in NV.suggestions("conformal", index=index)}
    assert "strategy" in kinds


def test_the_old_node_cards_do_not_come_back(index):
    """The modules replaced Gene/List/Label as the front door."""
    for s in NV.suggestions("map", index=index):
        assert not (s.kind == "gene" and s.title == "Map"), "the old Gene ▸ Map card resurfaced"


def test_an_empty_query_suggests_nothing():
    assert NV.suggestions("") == [] and NV.suggestions("   ") == []


def test_a_broken_index_does_not_break_the_field():
    def boom():
        raise RuntimeError("index unavailable")
    assert [s.kind for s in NV.suggestions("map", index=boom)][:1] == ["module"]


# ---------------------------------------------------------------- card text


def test_a_dataset_subtitle_drops_its_breadcrumb():
    """On a two-line card the slot-tree path used both lines and said nothing."""
    raw = "Slot tree ▸ Tg ▸ relational measurements ▸ host effect — host effect human K562; invasion"
    assert NV._subtitle(raw).startswith("host effect human K562")


def test_a_plain_subtitle_is_left_alone():
    assert NV._subtitle("facilitative glucose transporter GT1") == "facilitative glucose transporter GT1"


def test_clip_shortens_at_a_word_with_an_ellipsis():
    out = NV._clip("one two three four five six seven", 14)
    assert out.endswith("…") and len(out) <= 15 and " " not in out[-2:]


# --------------------------------------------------------------------- widget


@pytest.fixture
def page(qtbot, index):
    widget = NV.NavigatorPage(genes=GENES, index=index)
    qtbot.addWidget(widget)
    widget.resize(1200, 800)
    widget.show()
    return widget


def type_and_wait(page, qtbot, text):
    page.field.setText(text)
    page.refresh()


def test_typing_builds_a_row_of_cards_and_a_details_pane(page, qtbot):
    type_and_wait(page, qtbot, "rhoptry")
    assert page._cards and page.row.isVisible() and page.details.isVisible()
    assert page.current() is page._cards[0].suggestion


def test_return_opens_the_first_suggestion(page, qtbot):
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "map")
    page.open_current()
    assert seen == [("module", "map")]


def test_return_before_the_debounce_still_opens_what_was_typed(page, qtbot):
    """Typing then pressing Return immediately must not open the PREVIOUS row."""
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    page.field.setText("map")                  # debounce running, row not built yet
    page.open_current()
    assert seen == [("module", "map")]


def test_an_unavailable_module_explains_rather_than_opens(page, qtbot):
    seen = []
    page.chosen.connect(lambda kind, key: seen.append((kind, key)))
    type_and_wait(page, qtbot, "fill")
    assert page.current().key == "fill"
    page.open_current()
    assert seen == [], "Fill has nothing to write; opening it would show an empty view"
    assert "not available" in page.feedback.text()


def test_nothing_found_says_so(page, qtbot):
    type_and_wait(page, qtbot, "zzqxnotathing")
    assert not page._cards and "Nothing here matches" in page.feedback.text()


def test_ask_sits_inside_the_field_on_the_right(page, qtbot):
    """A child of a QLineEdit is drawn at (0, 0) unless something positions it; the
    first render put Ask above the text at the top left."""
    type_and_wait(page, qtbot, "rhoptry")
    page._place_ask()
    geo = page.ask_button.geometry()
    assert geo.x() > page.field.width() // 2
    assert 0 <= geo.y() and geo.bottom() <= page.field.height()


def test_ask_sends_the_text_to_the_assistant(page, qtbot):
    asked = []
    page.ask.connect(asked.append)
    page.field.setText("which GRAs reach the host nucleus?")
    page.ask_button.click()
    assert asked == ["which GRAs reach the host nucleus?"]


def test_escape_clears_the_field(page, qtbot):
    page.field.setText("rhoptry")
    qtbot.keyClick(page.field, QtCore.Qt.Key.Key_Escape)
    assert page.field.text() == ""


# --------------------------------------------------------------------- landing


@pytest.fixture
def landing(qtbot, index, tmp_path, monkeypatch):
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "ambient.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))
    from parascape.landing import Landing
    widget = Landing(genes=GENES, index=index)
    qtbot.addWidget(widget)
    return widget


def test_the_landing_opens_on_the_navigator(landing):
    assert landing.page == "Navigator"


def test_the_switch_reaches_the_modules(landing):
    landing.show_page("Modules")
    assert landing.page == "Modules" and landing.stack.currentWidget() is landing.home


def test_one_backdrop_serves_both_pages(landing):
    """Installing one per page would restart the animation on every switch and run
    two timers for one picture."""
    assert landing.animated
    assert landing.home.animated is False, "the tiles page must not install its own"
    before = landing._ambient
    landing.show_page("Modules")
    landing.show_page("Navigator")
    assert landing._ambient is before


def test_a_module_from_either_page_arrives_the_same_way(landing):
    seen = []
    landing.module_chosen.connect(seen.append)
    landing.navigator.field.setText("map")
    landing.navigator.open_current()
    landing.home.open_module("gene")
    assert seen == ["map", "gene"]


def test_a_gene_from_the_navigator_is_an_item_not_a_module(landing):
    seen = []
    landing.item_chosen.connect(lambda kind, key: seen.append((kind, key)))
    landing.navigator.field.setText("TGME49_223920")
    landing.navigator.open_current()
    assert seen == [("gene", "TGME49_223920")]


def test_the_theme_can_be_changed_live(landing):
    from parascape.backdrops import engines as E
    landing.set_theme("blobs")
    assert landing.animated
    landing.set_theme(E.NO_ANIMATION)
    assert not landing.animated
    # And the page paints its own colour again. Left transparent with no backdrop
    # behind it, it would show through to whatever is underneath.
    assert landing.property("ambient") == "off"
    assert landing.home.property("ambient") == "off"
    assert not landing.stack.testAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground)
