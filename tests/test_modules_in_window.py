"""The landing view inside a real window: every tile and result goes somewhere real.

The first version of the module catalogue was written from memory and pointed six
modules at tabs that did not exist -- Cell at a `maps` tab when the drawing lives in
`find`, Questions at `discoveries` when it is tab 8 of `analysis`, Datasets and Export
at tabs when they are window actions, and Structure at a view that has not been built.
Every one of those would have been a tile that opened nothing. These tests check the
catalogue against the window it runs in, so it cannot drift back.
"""
import pytest
from PySide6 import QtWidgets

from parascape import modules as MOD


@pytest.fixture(scope="module")
def window():
    from parascape import app as A
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                      # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    win.close()
    qapp.processEvents()


def _docks(win):
    return {d.windowTitle(): d for d in win.findChildren(QtWidgets.QDockWidget)}


@pytest.mark.parametrize("module", [m for m in MOD.MODULES if m.opens], ids=lambda m: m.key)
def test_every_destination_exists_in_the_window(window, module):
    if module.opens.startswith("action:"):
        method = module.opens.split(":", 1)[1]
        assert callable(getattr(window, method, None)), f"{module.key}: no window.{method}"
        return
    dock, _, tab = module.opens.partition(":")
    docks = _docks(window)
    assert dock in docks, f"{module.key}: no {dock!r} dock (have {sorted(docks)})"
    if tab:
        tabs = docks[dock].findChild(QtWidgets.QTabWidget)
        names = [tabs.tabText(i) for i in range(tabs.count())] if tabs else []
        assert tab in names, f"{module.key}: no {tab!r} tab in {dock} (have {names})"


def test_the_window_opens_on_the_landing_with_nothing_under_it_changed(window):
    """An overlay, not a replacement. The stack design moved the map out of the central
    widget; the dock-hiding design took docks out of their tab groups and made their
    children stop being visible to the window. This checks neither happens."""
    from parascape.landing_install import landing_up
    assert window.landing is not None, "the landing failed to install"
    assert landing_up(window)
    assert window.centralWidget() is window.view, "the map stays the central widget"
    assert window.landing.parentWidget() is window
    assert window.tabifiedDockWidgets(window.analysis_dock), "docks keep their tab groups"


def test_the_landing_covers_everything_between_the_menu_and_the_status_bar(window):
    from parascape.landing_install import fit_landing
    window.resize(1400, 900)
    fit_landing(window)
    geo = window.landing.geometry()
    assert geo.width() == window.width()
    assert geo.top() >= window.menuBar().geometry().bottom()


def test_the_map_module_takes_the_landing_away_and_home_brings_it_back(window):
    from parascape.landing_install import go_home, landing_up, open_module
    open_module(window, "map")
    assert not landing_up(window)
    go_home(window)
    assert landing_up(window)


def test_a_menu_command_takes_the_landing_away(window):
    """A panel opened from a menu must not appear underneath the overlay."""
    from parascape.landing_install import go_home, landing_up
    go_home(window)
    view = next(a.menu() for a in window.menuBar().actions()
                if a.text().replace("&", "") == "View")
    command = next(a for a in view.actions()
                   if not a.isSeparator() and a.menu() is None and a is not window.landing_home_action)
    view.triggered.emit(command)
    assert not landing_up(window)
    go_home(window)


def test_home_and_the_animation_choices_keep_the_landing(window):
    from parascape.landing_install import go_home, landing_up
    go_home(window)
    view = next(a.menu() for a in window.menuBar().actions()
                if a.text().replace("&", "") == "View")
    view.triggered.emit(window.landing_home_action)
    assert landing_up(window)
    first = window.landing_animation_menu.actions()[0]
    window.landing_animation_menu.triggered.emit(first)
    assert landing_up(window)


def test_the_dock_watcher_ignores_the_window_appearing(window):
    """Every dock reports becoming visible when the window first shows. Read as "a panel
    was opened", that took the landing away at launch -- invisible to any test that never
    shows the window, so this one shows it."""
    from parascape.landing_install import go_home, landing_up
    go_home(window)
    window.hide()
    window.show()
    QtWidgets.QApplication.processEvents()
    assert landing_up(window), "showing the window removed the landing"
    window.hide()


def test_a_dock_module_raises_its_dock_and_tab(window):
    from parascape.landing_install import go_home, open_module
    go_home(window)
    message = open_module(window, "questions")
    assert message == "analysis ▸ 8 · Questions"
    tabs = _docks(window)["analysis"].findChild(QtWidgets.QTabWidget)
    assert tabs.tabText(tabs.currentIndex()) == "8 · Questions"


def test_an_unavailable_module_reports_rather_than_opens(window):
    from parascape.landing_install import open_module
    message = open_module(window, "structure")
    assert "needs" in message


def test_home_is_reachable_from_the_view_menu_and_a_shortcut(window):
    action = window.landing_home_action
    # Not Ctrl+Shift+H: that is Help ▸ Search ParaScape, and taking it broke the search shortcut.
    assert action.shortcut().toString() == "Alt+Home"
    view = next(a.menu() for a in window.menuBar().actions()
                if a.text().replace("&", "") == "View")
    assert action in view.actions()


def test_the_animation_menu_offers_spacrs_themes(window):
    from parascape.backdrops import engines as E
    labels = [a.text() for a in window.landing_animation_menu.actions()]
    assert len(labels) == len(E.ANIMATION_CHOICES)
    checked = [a for a in window.landing_animation_menu.actions() if a.isChecked()]
    assert len(checked) == 1, "exactly one theme is the current one"


def test_the_window_still_builds_every_dock_it_built_before(window):
    """Nothing was removed: the landing is in front, not instead."""
    for title in ("find", "analysis", "maps", "strategies", "star map", "discoveries",
                  "compare", "information map", "evidence", "home", "console", "jobs",
                  "assistant", "gallery", "start here"):
        assert title in _docks(window), title


def test_tiles_for_docks_built_late_are_not_greyed(window):
    """Installed from `_home()`, the landing ran before `_tool_docks()` and greyed Jobs and
    Console as "not available in this window" -- docks that existed a moment later."""
    unavailable = window.landing.home.unavailable()
    for key in ("jobs", "console", "infer", "network", "datasets", "export"):
        assert key not in unavailable, f"{key} greyed: {unavailable.get(key)}"
    assert set(unavailable) == {"structure", "runs", "fill", "list"}, unavailable
