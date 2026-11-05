"""The home screen in a real window: the entry point, and every card's destination.

Two kinds of check. The ones that keep the promise -- the home dock is what is up on startup and
every other dock is still there and still reachable, because phase 1 removes nothing -- and the ones
that keep it honest: every card resolves to something this window actually has, an unavailable card
says why and changes nothing when clicked, every control carries a tooltip (a release failed its
checks for one missing tooltip), and a dropped file reports the identifiers it could not resolve
rather than quietly shortening the list.
"""
import os

import pandas as pd
import pytest
from PyQt6 import QtCore, QtGui, QtWidgets

from parascape import home_panel as HP
from parascape import nodes as ND

#: The docks that existed before the home screen, and must still exist after it.
EXISTING_DOCKS = ("find", "evidence", "analysis", "strategies", "start here", "star map",
                  "discoveries", "compare", "information map", "maps", "gallery", "jobs",
                  "console")


@pytest.fixture(scope="module")
def window():
    pytest.importorskip("PyQt6")
    from parascape import app as A
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                      # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    win.close()
    qapp.processEvents()


@pytest.fixture
def panel():
    QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    return HP.HomePanel(known_ids=["TGME49_208830", "TGME49_233460"])


def _docks(window):
    return {d.windowTitle(): d for d in window.findChildren(QtWidgets.QDockWidget)}


def _raised_dock(window) -> str:
    """The title of the dock tab that is up.

    Not `isVisible()`: the suite never shows the window, so every dock reports invisible. The tab
    bar Qt builds for a tabified dock area is what says which one of them is raised.
    """
    titles = set(_docks(window))
    for bar in window.findChildren(QtWidgets.QTabBar):
        labels = {bar.tabText(i) for i in range(bar.count())}
        if bar.count() and labels <= titles:
            return bar.tabText(bar.currentIndex())
    return ""


# --------------------------------------------------------------------------- the entry point
def test_the_home_dock_is_installed_and_raised_on_startup(window):
    assert isinstance(getattr(window, "home", None), HP.HomePanel)
    dock = window.home_dock
    assert dock.windowTitle() == "home"
    assert _raised_dock(window) == "home", "the home dock is not the tab that is up on startup"
    tabbed = window.tabifiedDockWidgets(dock)
    assert tabbed, "the home dock should share the tab group with the other right-hand docks"


def test_nothing_was_removed_to_make_room_for_the_home_screen(window):
    """Phase 5 retires a dock only once each has a node; this is phase 1."""
    present = _docks(window)
    missing = [name for name in EXISTING_DOCKS if name not in present]
    assert not missing, f"docks lost: {missing}"
    for name in EXISTING_DOCKS:
        present[name].show()
        present[name].raise_()
        assert present[name].widget() is not None


def test_every_node_resolves_to_a_destination_this_window_has(window):
    """The check that would have caught a card pointing at a panel that does not exist."""
    docks = _docks(window)
    unresolved = []
    for node in window.home.nodes:
        action = ND.resolve(node)
        if action.kind == "unavailable":
            assert node.unavailable_reason
        elif action.kind == "dock":
            if action.target not in docks:
                unresolved.append(f"{node.id}: no {action.target} dock")
            elif action.tab:
                tabs = docks[action.target].findChild(QtWidgets.QTabWidget)
                labels = [tabs.tabText(i) for i in range(tabs.count())] if tabs else []
                if action.tab not in labels:
                    unresolved.append(f"{node.id}: {action.target} has no {action.tab} tab")
        elif action.kind == "window":
            if not callable(getattr(window, action.target, None)):
                unresolved.append(f"{node.id}: window has no {action.target}()")
        elif action.kind == "workflow":
            assert callable(getattr(window, "_open_workflows", None))
        elif action.kind == "dialog":
            from parascape.data_download_dialog import DataDownloadDialog
            assert DataDownloadDialog is not None
        else:
            assert action.kind in ("catalogue", "question", "drop")
    assert not unresolved, unresolved


def test_opening_every_card_reports_something_and_raises_nothing(window):
    for node in ND.by_category(window.home.nodes):
        pass
    for card in window.home.cards:
        message = window.home.choose(card.node)
        assert message, f"{card.node.id} said nothing"
        if not card.node.available:
            assert card.node.unavailable_reason in message
    window.home_dock.raise_()


def test_an_unavailable_card_states_why_and_does_not_open_anything(window):
    node = ND.get("list_pulldown", window.home.nodes)
    before = window.home.pages.currentIndex()
    message = window.home.choose(node)
    assert "background frequency" in message
    assert window.home.pages.currentIndex() == before
    assert "background frequency" in window.home.note.text()


def test_a_dock_card_raises_its_dock_on_the_named_tab(window):
    message = HP.open_node(window, ND.get("gene_expression", window.home.nodes))
    assert "Gene expression" in message
    tabs = _docks(window)["compare"].findChild(QtWidgets.QTabWidget)
    assert tabs.tabText(tabs.currentIndex()) == "Gene expression"
    assert _raised_dock(window) == "compare"
    message = HP.open_node(window, ND.get("gene_network", window.home.nodes))
    assert "star map" in message
    assert _raised_dock(window) == "star map"
    window.home_dock.raise_()


def test_a_missing_destination_is_reported_rather_than_raised(window):
    node = ND.Node("gone", "Gone", "panel", "gene", "opens a dock that is not here",
                   ND.Destination("dock", "no such dock"))
    assert "not available" in HP.open_node(window, node)
    node = ND.Node("gone2", "Gone", "panel", "gene", "calls a method that is not here",
                   ND.Destination("window", "no_such_method"))
    assert "cannot be opened" in HP.open_node(window, node)


def test_the_search_field_is_the_one_from_help_search_and_finds_nodes(window):
    from parascape import help_index as HI
    from parascape import help_search as HS
    field = HS.field_of(window)
    results = field.type_and_search("pulldown_background")
    node_rows = [e for e in results if e.kind == "node"]
    assert node_rows, "the home nodes are not in the search index"
    assert HI.KIND_LABEL["node"]
    assert "background frequency" in HS.open_entry(window, node_rows[0])
    window.home_dock.raise_()


# --------------------------------------------------------------------------- the cards
def test_the_three_categories_appear_in_order_with_their_cards(panel):
    titles = [label.text() for label in panel.board.findChildren(QtWidgets.QLabel)
              if label.objectName() == "HomeSection"]
    assert [t.split(" — ")[0] for t in titles] == ["Gene", "List", "Label"]
    assert all(" — " in t for t in titles), "a section must state its purpose"
    assert len(panel.cards) == len([n for n in panel.nodes if n.level == 2])


def test_the_label_cards_are_drawn_in_the_settled_order(panel):
    titles = [c.node.title for c in panel.cards if c.node.category == "label"]
    assert titles == ["Catalogue", "Infer", "Strategies", "Scorecard", "Map", "Runs", "Fill",
                      "Compare"]


def test_every_interactive_control_has_a_tooltip(panel):
    """A release failed CI for one missing tooltip; a card with no tooltip is a card nobody can
    ask about."""
    panel.show_page(HP.PAGE_DROP)
    panel.show_page(HP.PAGE_CATALOGUE)
    kinds = (QtWidgets.QAbstractButton, QtWidgets.QLineEdit, QtWidgets.QComboBox,
             QtWidgets.QTableWidget, QtWidgets.QListWidget, QtWidgets.QTextBrowser)
    missing = []
    for kind in kinds:
        for widget in panel.findChildren(kind):
            if isinstance(widget.parent(), QtWidgets.QLineEdit) \
                    or widget.objectName().startswith("qt_"):
                continue          # parts Qt puts inside a line edit or a table, not our controls
            if not widget.toolTip().strip():
                missing.append(f"{type(widget).__name__} {widget.objectName() or widget.text()!r}"
                               if hasattr(widget, "text") else type(widget).__name__)
    assert not missing, f"no tooltip: {missing}"
    panel.show_page(HP.PAGE_HOME)


def test_a_card_tooltip_states_the_claim_rule_and_the_scoreboard(panel):
    facts = next(c for c in panel.cards if c.node.id == "gene_facts")
    tip = facts.tip_text()
    assert "probability" in tip and "uncertainty" in tip and "source" in tip
    questions = ND.children("label_catalogue", panel.nodes)
    scored = next(n for n in questions if n.scoreboard and n.scoreboard.scored)
    card = HP.NodeCard(scored)
    assert "when the hidden genes rotate" in card.tip_text() and "on hidden genes" in card.tip_text()
    assert "not a negative" in card.tip_text()


def test_an_unavailable_card_is_visibly_distinct_and_carries_its_reason(panel):
    card = next(c for c in panel.cards if c.node.id == "list_pulldown")
    assert not card.node.available
    assert any("unavailable:" in line for line in card.lines())
    assert "Not available" in card.tip_text()
    card.resize(HP.CARD_WIDTH, HP.CARD_HEIGHT)
    card.grab()                                 # painting the dashed, dimmed state must not raise


def test_search_filters_and_ranks_the_cards(panel):
    assert [n.title for n in panel.filter_cards("network")] == ["Network"]
    assert [c.node.title for c in panel.shown_cards()] == ["Network"]
    assert [n.title for n in panel.filter_cards("pulldown_background")] == ["Pulldown"]
    assert "result(s) for" in panel.note.text()
    assert len(panel.filter_cards("")) == len(panel.cards)
    panel.filter_cards("glucose")
    assert panel.shown_cards() == [], "a question is not a card"
    assert "Catalogue" in panel.note.text(), "it must say where the question is"
    panel.search.setText("")


def test_the_cards_are_keyboard_navigable_and_return_opens_one(panel):
    opened = []
    panel.node_chosen.connect(lambda node: opened.append(node.id))
    panel.resize(900, 700)
    panel.show()                        # focus needs a shown widget, even offscreen
    QtWidgets.QApplication.processEvents()
    first, second = panel.cards[0], panel.cards[1]
    first.setFocus(QtCore.Qt.FocusReason.OtherFocusReason)
    assert first.hasFocus()
    panel.keyPressEvent(QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Right,
                                        QtCore.Qt.KeyboardModifier.NoModifier))
    assert second.hasFocus(), "an arrow key must move between cards"
    second.keyPressEvent(QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Return,
                                         QtCore.Qt.KeyboardModifier.NoModifier))
    assert opened == [second.node.id]


# --------------------------------------------------------------------------- the catalogue
def test_the_label_catalogue_shows_both_numbers_and_a_state_per_question(panel):
    panel.choose(ND.get("label_catalogue", panel.nodes))
    assert panel.pages.currentIndex() == HP.PAGE_CATALOGUE
    view = panel.catalogue_view
    assert "Effective independent axes" in view.header.text()
    assert list(view.COLUMNS).count("right when rotated") == 1
    assert "right on hidden genes" in view.COLUMNS
    row = next(i for i, n in enumerate(view.nodes)
               if n.question_id == "required_without_glucose")
    values = [view.table.item(row, c).text() for c in range(len(view.COLUMNS))]
    assert any(v.startswith("≥") for v in values), values
    states = [view.table.item(i, len(view.COLUMNS) - 1).text() for i in range(len(view.nodes))]
    assert "cannot be answered yet" in states, "a declared gap must render as a state, not a page"
    view.show_detail(row)
    assert "Scoreboard" in view.detail.toHtml() and "Coverage" in view.detail.toHtml()
    assert "Forbidden evidence" in view.detail.toHtml()
    panel.show_page(HP.PAGE_HOME)


def test_opening_one_question_selects_it_in_the_catalogue(panel):
    node = next(n for n in ND.children("label_catalogue", panel.nodes)
                if n.question_id == "required_without_glucose")
    message = panel.choose(node)
    assert panel.pages.currentIndex() == HP.PAGE_CATALOGUE
    assert "when the hidden genes rotate" in message and "on hidden genes" in message
    assert panel.catalogue_view.nodes[panel.catalogue_view.table.currentRow()].question_id \
        == "required_without_glucose"
    panel.show_page(HP.PAGE_HOME)


def test_fill_shows_only_what_cleared_the_bar_or_says_there_is_nothing(panel):
    node = ND.get("label_fill", panel.nodes)
    if node.available:
        panel.choose(node)
        assert all(n.scoreboard.may_fill for n in panel.catalogue_view.nodes)
    else:
        assert "frozen-test" in panel.choose(node)
    panel.show_page(HP.PAGE_HOME)


# --------------------------------------------------------------------------- a dropped file
def _written(tmp_path, name):
    table = pd.DataFrame({"gene": ["TGME49_208830", "TgME49.233460", "nonsense"],
                          "enrichment": [12.0, 3.0, 1.0]})
    path = tmp_path / name
    if name.endswith(".xlsx"):
        table.to_excel(path, index=False, sheet_name="hits")
    else:
        table.to_csv(path, index=False)
    return str(path)


@pytest.mark.parametrize("name", ["hits.csv", "hits.xlsx"])
def test_a_dropped_file_reaches_the_drop_node(panel, tmp_path, name):
    assert panel.load_file(_written(tmp_path, name))
    assert panel.wait_for_file()
    QtWidgets.QApplication.processEvents()
    assert panel.pages.currentIndex() == HP.PAGE_DROP
    assert panel.drop is not None and panel.drop.name == name
    assert "TGME49_208830" in panel.genes()
    assert os.path.basename(name) in panel.drop_view.headline.text()
    assert panel.drop_view.column.currentText() == "gene"
    assert panel.drop_view.extras.count() >= 1


def test_an_unresolved_identifier_is_shown_and_not_discarded(panel, tmp_path):
    panel.load_file(_written(tmp_path, "hits.csv"))
    panel.wait_for_file()
    QtWidgets.QApplication.processEvents()
    table = panel.drop_view.report
    inputs = [table.item(r, 0).text() for r in range(table.rowCount())]
    assert inputs == ["TGME49_208830", "TgME49.233460", "nonsense"]
    states = {table.item(r, 0).text(): table.item(r, 2).text() for r in range(table.rowCount())}
    assert states["nonsense"] == "unresolved"
    assert table.item(inputs.index("nonsense"), 3).text()
    assert "not discarded" in panel.drop_view.summary.text()


def test_the_list_actions_wake_up_once_a_list_is_loaded(panel, tmp_path):
    assert not panel.drop_view.buttons["list_answer"].isEnabled()
    panel.load_file(_written(tmp_path, "hits.csv"))
    panel.wait_for_file()
    QtWidgets.QApplication.processEvents()
    assert panel.drop_view.buttons["list_answer"].isEnabled()
    assert panel.drop_view.buttons["list_map"].isEnabled()
    assert not panel.drop_view.buttons["list_pulldown"].isEnabled(), "there is no background frequency yet"
    assert "background frequency" in panel.drop_view.buttons["list_pulldown"].toolTip()
    assert "Not available" not in panel.drop_view.buttons["list_answer"].toolTip(), \
        "an enabled button must stop saying it needs a list"


def test_an_unsupported_file_says_what_is_accepted(panel, tmp_path):
    path = tmp_path / "hits.docx"
    path.write_text("not a table")
    assert panel.load_file(str(path)) is False
    for extension in HP.GD.READABLE:
        assert extension in panel.note.text()
    assert ".docx" in panel.note.text()


def test_the_window_accepts_a_file_dropped_anywhere(window, tmp_path):
    path = _written(tmp_path, "dropped.csv")
    mime = QtCore.QMimeData()
    mime.setUrls([QtCore.QUrl.fromLocalFile(path)])
    event = QtGui.QDropEvent(QtCore.QPointF(10, 10), QtCore.Qt.DropAction.CopyAction, mime,
                             QtCore.Qt.MouseButton.LeftButton,
                             QtCore.Qt.KeyboardModifier.NoModifier)
    assert window._home_drops.eventFilter(window, event)
    assert window.home.wait_for_file()
    QtWidgets.QApplication.processEvents()
    assert _raised_dock(window) == "home", "a dropped file must bring the home screen up"
    assert window.home.pages.currentIndex() == HP.PAGE_DROP
    assert window.home.drop is not None and window.home.drop.name == "dropped.csv"
    assert window.home.genes()


def test_a_loaded_list_makes_the_list_cards_available(panel, tmp_path):
    panel.load_file(_written(tmp_path, "hits.csv"))
    panel.wait_for_file()
    QtWidgets.QApplication.processEvents()
    loader = panel._loader
    if loader is not None:
        loader.wait(30000)
        QtWidgets.QApplication.processEvents()
    assert ND.get("list_answer", panel.nodes).available, \
        "once a list is loaded its cards must stop being greyed"
    assert ND.get("list_pulldown", panel.nodes).available is False


# --------------------------------------------------------------------------- the navigator
def test_one_search_bar_reaches_genes_and_the_index(window):
    """Einar's navigator: one field for everything, and a gene result lands on THAT gene's card."""
    panel = window.home
    gene = str(window.nodes["gene_id"].iloc[0])
    found = panel.filter_cards(gene)
    assert found and found[0].destination.kind == "gene"
    assert found[0].destination.target == gene
    card = next(c for c in panel.cards if c.node.id == f"gene:{gene}")
    assert not card.isHidden(), "a found gene must be drawn as a card"
    message = panel.choose(found[0])
    assert gene in message
    assert str(window.search.text()) == gene, "the gene must be selected, not just announced"
    assert gene in window.detail.toPlainText(), "its card must be the one showing"
    rows = panel.filter_cards("hdbscan")
    assert any(n.destination.kind == "entry" for n in rows), "strategies must come from one bar too"
    assert "could not open" not in panel.choose(next(n for n in rows
                                                     if n.destination.kind == "entry"))
    panel.search.setText("")
    panel.filter_cards("")
    window.home_dock.raise_()


def test_down_from_the_search_box_moves_into_the_results(panel):
    panel.resize(900, 700)
    panel.show()
    QtWidgets.QApplication.processEvents()
    panel.filter_cards("network")
    panel.search.setFocus(QtCore.Qt.FocusReason.OtherFocusReason)
    QtWidgets.QApplication.processEvents()
    down = QtGui.QKeyEvent(QtCore.QEvent.Type.KeyPress, QtCore.Qt.Key.Key_Down,
                           QtCore.Qt.KeyboardModifier.NoModifier)
    assert panel.eventFilter(panel.search, down) is True
    assert panel.shown_cards_widgets()[0].hasFocus()


def test_a_found_section_appears_and_disappears_with_the_query(panel):
    panel.gene_lookup = [("TGME49_208830", "GRA16 dense granule protein")]
    before = len(panel.cards)
    panel.filter_cards("GRA16")
    sections = [label.text() for label in panel.board.findChildren(QtWidgets.QLabel)
                if label.objectName() == "HomeSection"]
    assert sections[0].startswith("Found"), sections
    assert len(panel.cards) > before
    panel.filter_cards("")
    sections = [label.text() for label in panel.board.findChildren(QtWidgets.QLabel)
                if label.objectName() == "HomeSection"]
    assert not any(s.startswith("Found") for s in sections)
    assert len(panel.cards) == before


def test_the_cards_follow_spacr_s_tile_treatment(panel):
    """spaCR's rule: a launcher on the home pane earns an edge only when it is interactive."""
    from parascape import theme as TH
    assert HP.CARD_RADIUS == TH.RADIUS["tile"]
    assert HP.GAP is TH.SPACING
    panel.resize(900, 700)
    panel.show()                       # focus is only real on a shown widget, even offscreen
    QtWidgets.QApplication.processEvents()
    card = panel.cards[0]
    card.resize(HP.CARD_WIDTH, HP.CARD_HEIGHT)
    resting = card.grab().toImage()
    card.setFocus(QtCore.Qt.FocusReason.OtherFocusReason)
    QtWidgets.QApplication.processEvents()
    assert card.grab().toImage() != resting, "focus must be visible on the card"
