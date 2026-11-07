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
    from parascape.map_module import wait_for_projections
    wait_for_projections()
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


def test_the_map_module_opens_its_shell_and_the_3d_map_is_one_action_away(window):
    """Instruction 86: Map is a module in spaCR shell; its "Full 3D view" action is
    what takes the landing away now, and Home comes back to the home page."""
    from parascape.landing_install import go_home, landing_up, open_module
    open_module(window, "map")
    assert landing_up(window) and window.landing.page == "Module"
    window.map_module.actions["Full 3D view"].click()
    assert not landing_up(window)
    go_home(window)
    assert landing_up(window) and window.landing.page != "Module"


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


def test_a_tool_module_lends_its_panel_with_the_right_tab_forward(window):
    from parascape.landing_install import go_home, open_module
    window.show()
    try:
        go_home(window)
        assert open_module(window, "expression") == "Expression"
        QtWidgets.QApplication.processEvents()
        lent = window.tool_modules["expression"].module.main.panel
        tabs = [lent] + lent.findChildren(QtWidgets.QTabWidget)
        tabs = next(t for t in tabs if isinstance(t, QtWidgets.QTabWidget))
        assert tabs.tabText(tabs.currentIndex()) == "Gene expression"
        go_home(window)
    finally:
        window.hide()


def test_questions_and_infer_open_the_infer_module(window):
    from parascape.landing_install import go_home, open_module
    for key in ("questions", "infer"):
        go_home(window)
        assert open_module(window, key) == "Infer"
        assert window.landing.page == "Module"


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
                  "compare", "information map", "evidence", "gene list", "console", "jobs",
                  "assistant", "gallery", "start here"):
        assert title in _docks(window), title


def test_tiles_for_docks_built_late_are_not_greyed(window):
    """Installed from `_home()`, the landing ran before `_tool_docks()` and greyed Jobs and
    Console as "not available in this window" -- docks that existed a moment later."""
    unavailable = window.landing.home.unavailable()
    for key in ("jobs", "console", "infer", "network", "datasets", "export"):
        assert key not in unavailable, f"{key} greyed: {unavailable.get(key)}"
    assert set(unavailable) == {"structure", "runs", "fill", "list"}, unavailable


def test_the_navigator_can_be_turned_off_from_view_and_help(window, tmp_path, monkeypatch):
    """Einar: "I couldn't turn off the navigator in View or Tools like I can in spaCR"."""
    from PySide6 import QtCore
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "landing.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))
    action = window.landing_navigator_action
    for menu_name in ("View", "Help"):
        menu = next(a.menu() for a in window.menuBar().actions()
                    if a.text().replace("&", "") == menu_name)
        assert action in menu.actions(), f"no Navigator switch in {menu_name}"
    action.setChecked(False)
    assert window.landing.page == "Modules"
    assert not any(b.isVisibleTo(window.landing) for b in window.landing._switch.buttons())
    action.setChecked(True)
    assert window.landing.page == "Navigator"


def test_leaving_the_landing_never_shows_a_second_home(window):
    """Einar: "I used the navigator, typed map, clicked map and was sent to the old home screen
    with tiles." The old node-card dock was the raised right-hand tab; it is now the gene list,
    not raised, and draws no cards."""
    from parascape.landing_install import go_home, open_module
    from parascape.home_panel import DOCK_TITLE, PAGE_DROP
    assert "home" not in _docks(window), "the old home dock is still installed"
    go_home(window)
    open_module(window, "map")
    dock = window.home_dock
    assert dock.windowTitle() == DOCK_TITLE
    assert window.home.with_cards is False
    assert window.home.search.isHidden(), "the old navigator field is still drawn"
    assert window.home.pages.currentIndex() == PAGE_DROP


def test_the_drop_module_opens_the_gene_list(window):
    from parascape.landing_install import go_home, open_module
    go_home(window)
    assert open_module(window, "drop") == "Drop"
    assert window.tool_modules["drop"].module.main.dock is window.home_dock


def test_the_dock_lists_home_and_every_module_and_opens_them(window):
    """Einar: "there should also be a dock like there is in spaCR" (instruction 95)."""
    from parascape import modules as MOD
    from parascape.landing_install import go_home, landing_up
    dock = window.landing.dock
    assert dock is not None, "the dock failed to build"
    keys = [row.key for row in dock._rows]
    assert keys[0] == "home" and set(keys[1:]) == {m.key for m in MOD.MODULES}
    go_home(window)
    seen = []
    window.landing.module_chosen.connect(lambda key: seen.append(key))
    dock.nav_selected.emit("gene")
    assert seen == ["gene"]
    dock.nav_selected.emit("map")
    assert window.landing.page == "Module", "the dock opens the Map module"
    dock.nav_selected.emit("home")
    assert landing_up(window) and window.landing.page != "Module"


def test_every_tile_with_a_panel_opens_a_module_in_the_same_shell(window):
    """Einar: "all modules need to be thematically the same"."""
    from parascape import modules as MOD
    from parascape.landing_install import go_home, open_module
    from parascape.shell import ModuleShell
    for m in MOD.MODULES:
        if MOD.unavailable_reason(m) or m.opens.startswith("action:"):
            continue
        go_home(window)
        open_module(window, m.key)
        current = window.landing.stack.currentWidget()
        assert isinstance(current, ModuleShell), f"{m.key} did not open a module shell"


def test_a_lent_panel_goes_back_to_its_dock(window):
    from parascape.landing_install import go_home, open_module
    window.show()
    try:
        open_module(window, "jobs")
        QtWidgets.QApplication.processEvents()
        lent = window.tool_modules["jobs"].module.main.panel
        assert lent is not None and window.tool_modules["jobs"].isAncestorOf(lent)
        go_home(window)
        QtWidgets.QApplication.processEvents()
        assert window.jobs_dock.widget() is lent
    finally:
        window.hide()


def test_partners_reads_the_network_by_binding_evidence(window):
    from parascape.landing_install import open_module
    open_module(window, "partners")
    module = window.network_module.module
    assert set(module.layers()) <= set(module.BINDING_LAYERS) and module.layers()
    module._relayout()
