"""spaCR's module shell (instruction 86) and the Map module, the first one moved onto it.

Einar: "spaCR has its settings in a container on the left, a console, AI chatbox, and an
actions section. There is a live section to live test parameters and a hyperparameter
search ... spaCR containers are also collapsible so after setting the settings you can
give more space to the content."
"""
import pytest
from PySide6 import QtCore, QtWidgets


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """Folds and pane sizes go to a scratch file, not the user's settings."""
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "shell.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))
    return path


def test_the_vendored_panes_import_without_spacr(qapp):
    from parascape.panes import collapsible_section, collapsible_splitter, foldable  # noqa: F401
    import sys
    assert not any(name == "spacr" or name.startswith("spacr.") for name in sys.modules)


def test_settings_fold_away_to_give_the_content_the_room(qapp, settings):
    from parascape.shell import ModuleShell
    shell = ModuleShell("t", "Test", "what it is for")
    shell.add_setting("Group", QtWidgets.QLabel("a setting"))
    shell.add_section("Live preview", QtWidgets.QLabel("content"))
    shell.resize(1200, 700)
    shell.show()
    qapp.processEvents()
    before = shell.body.sizes()
    assert before[0] > 0
    shell.body.toggle_pane("Settings", by_user=True)
    qapp.processEvents()
    after = shell.body.sizes()
    assert after[0] == 0 and after[1] > before[1], (before, after)
    shell.close()


def test_actions_appear_with_the_first_and_report_under_the_row(qapp, settings):
    from parascape.shell import ModuleShell
    shell = ModuleShell("t2", "Test")
    assert "Actions" not in shell.sections
    clicked = []
    shell.add_action("go", lambda: clicked.append(1), tip="do it", primary=True)
    assert "Actions" in shell.sections
    shell.actions["go"].click()
    shell.say("done")
    assert clicked == [1] and shell.status_line.text() == "done"


def test_a_borrowed_panel_goes_back_to_its_dock(qapp, settings):
    from parascape.shell import ModuleShell
    window = QtWidgets.QMainWindow()
    dock = QtWidgets.QDockWidget("console", window)
    panel = QtWidgets.QLabel("the console")
    dock.setWidget(panel)
    shell = ModuleShell("t3", "Test")
    shell.borrow("Console", [("Console", dock)])
    shell.show()
    qapp.processEvents()
    assert dock.widget() is not panel and shell.isAncestorOf(panel)
    shell.hide()
    qapp.processEvents()
    assert dock.widget() is panel
    shell.show()
    qapp.processEvents()
    assert shell.isAncestorOf(panel), "a second visit borrows it again"
    shell.hide()


# ------------------------------------------------------------------ in the window
@pytest.fixture(scope="module")
def window(qapp):
    from parascape import app as A
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                      # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    from parascape.map_module import wait_for_projections
    wait_for_projections()
    win.close()
    qapp.processEvents()


def test_the_map_controls_moved_into_the_module_not_copied(window):
    shell = window.map_module
    assert shell is not None, "the Map module failed to build"
    panel = window.panel
    for name in ("nn", "md", "seed", "sample", "nn_grid", "md_grid", "mcs_grid", "emb_name"):
        assert shell.isAncestorOf(getattr(panel, name)), f"{name} is not in the Map module"
    assert panel.spec().n_neighbors == panel.nn.value(), "the panel still reads the moved control"


def test_the_old_map_tab_points_to_the_module(window):
    tabs = window.panel.findChild(QtWidgets.QTabWidget)
    index = next(i for i in range(tabs.count()) if tabs.tabText(i) == "2 · Map")
    texts = [b.text() for b in tabs.widget(index).findChildren(QtWidgets.QPushButton)]
    assert "Open the Map module" in texts


def test_the_map_module_is_built_around_its_goal(window):
    """Instruction 96, agreed with Einar: your genes and the colouring first, what makes
    genes similar in plain words, the UMAP settings folded under Advanced; on the right the
    map itself, the neighbourhood reading, and finding the best map for a label."""
    shell = window.map_module
    sections = list(shell.sections)
    for part in ("Your genes", "Colour by", "What makes genes similar", "Advanced",
                 "Map", "Your genes and their neighbours", "Find the best map for a label",
                 "Actions", "Console and AI chat"):
        assert part in sections, f"no {part}"
    assert sections.index("Your genes") < sections.index("Advanced")
    assert not shell.sections["Advanced"].is_expanded(), "the UMAP settings start folded"
    assert shell.sections["Console and AI chat"].folder.shut, "the console starts folded"
    actions = set(shell.actions)
    assert {"Build this map", "Predict from neighbours", "Find the best map",
            "Full 3D view"} <= actions


def test_the_map_shown_is_the_current_map_not_a_sample(window):
    module = window.map_module.module
    module.show_current()
    assert module.mode == "map"
    assert len(module.view.ids) == int(window.visible_mask().sum())


def test_a_typed_gene_is_put_on_the_map_and_read(window):
    module = window.map_module.module
    gid = str(window.nodes.gene_id.iloc[0])
    module.gene_field.setText(gid)
    module._genes_typed()
    assert module.selected == [0] and 0 in module.view.highlight
    html = module.reading.toHtml()
    assert gid in html and "nearest genes" in html


def test_a_list_file_is_put_on_the_map(window, tmp_path):
    module = window.map_module.module
    ids = window.nodes.gene_id.astype(str).tolist()[:5]
    path = tmp_path / "hits.csv"
    path.write_text("gene\n" + "\n".join(ids + ["NOT_A_GENE"]) + "\n")
    module.load_list(str(path))
    assert module.selected == [0, 1, 2, 3, 4]
    assert "5 genes" in module.reading.toHtml()


def test_ticking_a_kind_of_measurement_ticks_its_blocks(window):
    module = window.map_module.module
    tree = window.panel.blocks
    box, item = next((b, it) for b, it in module.similarity.items() if not b.isChecked())
    module._similarity_clicked(item, True)
    keys = {leaf.data(0, QtCore.Qt.ItemDataRole.UserRole) for leaf in tree.leaves(item)
            if leaf.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled}
    assert keys <= set(tree.checked()) and box.isChecked()
    module._similarity_clicked(item, False)
    assert not keys & set(tree.checked()) and not box.isChecked()


def test_colour_by_is_the_window_choice(window):
    module = window.map_module.module
    other = next(t for t in window.color_sources() if t != module.colour_box.currentText())
    module.colour_box.setCurrentText(other)
    assert window.category_box.currentText() == other


def test_the_console_is_in_the_module_while_it_is_up_and_back_after(window):
    """Borrowing follows the shell on and off screen, so the window has to be shown."""
    from parascape.landing_install import go_home, open_module
    window.show()
    open_module(window, "map")
    QtWidgets.QApplication.processEvents()
    try:
        assert window.map_module.isAncestorOf(window.console)
        go_home(window)
        QtWidgets.QApplication.processEvents()
        assert window.console_dock.widget() is window.console
        assert window.chat_dock.widget() is window.chat
    finally:
        window.map_module.module.live.setChecked(False)
        window.hide()
