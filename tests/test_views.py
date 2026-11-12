"""Every view the catalogue declares is real, named and reachable -- and the reverse.

Instruction 98, WP1 counted what WP10 had to fix: 29 tabs inside five docks and four
further windows that **the navigator could not name at all**. That is the whole of Einar
complaint that "everything behind the navigator seems to be unchanged": the modules were
in the right shell and the views inside them had no names anybody could type.

WP10 declares each of them on the module whose panel holds it
(:class:`parascape.modules.View`). These tests hold the two directions that matter:

1. a declared view exists in a real window -- the tab is there, or the window method is;
2. a tab that exists is declared, or is a module in its own right, or is a signpost left
   where a module took its controls.

The second is the direction the first version of the module catalogue got wrong: it was
written from memory and pointed six modules at views that did not exist. A table checked
only against itself is the bug spaCR records as #16i.
"""
import pytest
from PySide6 import QtWidgets

from parascape import modules as MOD


# ------------------------------------------------------- what can be checked with no Qt

def test_every_view_says_what_it_shows_and_where_it_is():
    for module, view in MOD.views():
        assert view.name.strip(), f"{module.key} has a view with no name"
        assert len(view.shows.split()) >= 6, \
            f"{module.key} ▸ {view.name}: {view.shows!r} says nothing about the data"
        assert view.opens.strip(), f"{module.key} ▸ {view.name} opens nothing"


def test_a_view_is_declared_on_exactly_one_module():
    """Two modules claiming one tab is two places to look for it, which is the thing the
    one catalogue exists to stop."""
    seen = {}
    for module, view in MOD.views():
        key = (view.opens, view.name.lower())
        assert key not in seen, f"{view.name} is on {seen.get(key)} and on {module.key}"
        seen[key] = module.key


def test_no_view_repeats_the_destination_of_the_module_it_is_on():
    """The module already opens its own destination; declaring it again as a view would
    put the same thing on screen under two names."""
    for module, view in MOD.views():
        assert view.opens != module.opens, \
            f"{module.key} ▸ {view.name} is the module itself, not a view of it"


def test_a_view_is_found_by_its_name_and_by_what_it_shows():
    module, view = MOD.find_view("Condition comparison")
    assert module.key == "expression" and view.opens == "compare:Condition comparison"
    by_name = MOD.search_views("coexpression")
    assert by_name and by_name[0][1].name == "Coexpression"
    by_what = MOD.search_views("slot")
    assert any(v.name == "Slot tree" for _m, v in by_what)


def test_an_unknown_view_raises_rather_than_answering_none():
    with pytest.raises(KeyError):
        MOD.find_view("a view nobody declared")


def test_the_views_cover_the_tabs_the_inventory_found_orphaned():
    """WP1 named them. Each had to end up routed, a module, or deleted, and these are
    the ones that were routed, so the count is held here rather than in prose."""
    opens = {v.opens for _m, v in MOD.views()}
    for tab in ("1 · Data", "4 · Inference", "5 · Search", "6 · Validation", "7 · Discover"):
        assert f"analysis:{tab}" in opens, tab
    for tab in ("Studies & samples", "Coexpression", "Published expression",
                "Deposited matrices", "Other measurements", "Sparse measurements",
                "Original workbooks", "Reference identifiers", "Dataset coverage",
                "Condition comparison"):
        assert f"compare:{tab}" in opens, tab
    for tab in ("Labels and functions", "Inferred claims", "Functional tests",
                "Strategy coverage", "Host source status"):
        assert f"discoveries:{tab}" in opens, tab
    assert {"information map:One gene", "information map:One trait"} <= opens
    assert {"strategies:Guide", "strategies:Settings", "strategies:Results"} <= opens
    # The two docks that had no module at all, and the window reached only by a link.
    assert {"maps", "start here", "star map"} <= opens
    assert {"action:open_slot_tree", "action:open_gene_evidence",
            "action:open_label_browser", "action:open_question_scoreboard"} <= opens


# ------------------------------------------------------------------- in a real window

@pytest.fixture(scope="module")
def window():
    from parascape import app as A
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                   # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    from parascape.map_module import wait_for_projections
    wait_for_projections()
    win.close()
    qapp.processEvents()


def _docks(win):
    return {d.windowTitle(): d for d in win.findChildren(QtWidgets.QDockWidget)}


def _current_tab(win, dock_title: str) -> str:
    """The tab on top of the panel of one dock, wherever that panel currently lives.

    A borrowed panel sits under the module shell while the module is on screen, so asking
    the dock alone would answer about an empty frame.
    """
    dock = _docks(win)[dock_title]
    panel = dock.widget()
    for shell in win.tool_modules.values():
        borrowed = getattr(getattr(shell, "module", None), "main", None)
        if getattr(borrowed, "dock", None) is dock and getattr(borrowed, "panel", None):
            panel = borrowed.panel
    tabs = panel if isinstance(panel, QtWidgets.QTabWidget) else \
        panel.findChild(QtWidgets.QTabWidget)
    return tabs.tabText(tabs.currentIndex()) if tabs else ""


def _tabs_of(dock):
    tabs = dock.findChild(QtWidgets.QTabWidget)
    return [tabs.tabText(i) for i in range(tabs.count())] if tabs else []


@pytest.mark.parametrize("pair", list(MOD.views()),
                         ids=lambda p: f"{p[0].key}:{p[1].name}")
def test_every_declared_view_exists_in_the_window(window, pair):
    module, view = pair
    if view.opens.startswith("action:"):
        method = view.opens.split(":", 1)[1]
        assert callable(getattr(window, method, None)), \
            f"{module.key} ▸ {view.name}: no window.{method}"
        return
    dock, _, tab = view.opens.partition(":")
    docks = _docks(window)
    assert dock in docks, f"{module.key} ▸ {view.name}: no {dock!r} dock"
    if tab:
        assert tab in _tabs_of(docks[dock]), \
            f"{module.key} ▸ {view.name}: no {tab!r} tab in {dock}"


def test_every_tab_in_the_window_is_a_view_a_module_or_a_signpost(window):
    """The direction a catalogue cannot answer on its own: a tab nobody can name.

    A tab is allowed to be unnamed only when a module TOOK its controls and a signpost to
    that module was left in its place -- which is what "move, do not copy" leaves behind
    (`map_module._replace_tab`).
    """
    moved = {"analysis:2 · Map", "analysis:3 · Clusters", "analysis:8 · Questions"}
    named = {v.opens for _m, v in MOD.views()} | {m.opens for m in MOD.MODULES} | moved
    docks = _docks(window)
    orphans = []
    for title in ("analysis", "compare", "discoveries", "information map", "strategies"):
        for tab in _tabs_of(docks[title]):
            if f"{title}:{tab}" not in named:
                orphans.append(f"{title}:{tab}")
    assert not orphans, f"tabs nobody can name: {orphans}"


def test_every_dock_is_a_module_or_holds_a_declared_view(window):
    """The same for the docks. `evidence`, `find`, `gene list`, `console`, `jobs`,
    `gallery` and `assistant` are modules; `maps`, `start here` and `star map` had no
    module at all, which is why they are declared views now."""
    named = {m.opens.partition(":")[0] for m in MOD.MODULES}
    named |= {v.opens.partition(":")[0] for _m, v in MOD.views()}
    named |= {"assistant"}          # folded into every module, and in the Tools menu
    missing = [t for t in _docks(window) if t not in named]
    assert not missing, f"docks nobody can name: {missing}"


def test_opening_a_view_reports_what_it_opened_and_raises_nothing(window):
    from parascape.landing_install import open_view
    for module, view in MOD.views():
        if view.opens.startswith("action:"):
            continue                 # these open windows; `test_app_menus` exercises those
        message = open_view(window, module.key, view.name)
        assert message, f"{module.key} ▸ {view.name} reported nothing"
        assert "not available" not in message, f"{module.key} ▸ {view.name}: {message}"


def test_a_view_inside_the_panel_of_its_module_opens_in_the_module(window):
    """The point of routing rather than raising a dock: Einar asked for every module to
    look alike, so a tab of a panel a module already shows arrives inside that shell."""
    from parascape.landing_install import open_view
    message = open_view(window, "expression", "Coexpression")
    assert "Coexpression" in message
    assert window.landing.page == "Module"
    assert window.landing.stack.currentWidget() is window.tool_modules["expression"]
    assert _current_tab(window, "compare") == "Coexpression"


def test_a_view_the_search_field_finds_opens_through_the_one_route(window):
    """The search field and the navigator reach a view the same way, because both go
    through `help_search.open_entry` onto `landing_install.open_view`."""
    from parascape import help_index as HI
    from parascape import help_search as HS
    rows = [e for e in HI.panel_entries(window) if e.data.get("view") == "One trait"]
    assert rows, "the search index has no row for the One trait view"
    HS.open_panel(window, rows[0])
    assert _current_tab(window, "information map") == "One trait"


def test_asking_for_a_view_a_module_does_not_have_says_so(window):
    """It reports, never raises: this is reached from a click on a search result, and a
    stale result naming a view that has since been renamed must not take the window."""
    from parascape.landing_install import open_view
    assert "no view" in open_view(window, "gene", "a view nobody declared")
    assert "no" in open_view(window, "not_a_module", "Slot tree")
