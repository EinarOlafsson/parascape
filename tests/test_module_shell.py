"""spaCR's module shell (instruction 86) and the Map module, the first one moved onto it.

Einar: "spaCR has its settings in a container on the left, a console, AI chatbox, and an
actions section. There is a live section to live test parameters and a hyperparameter
search ... spaCR containers are also collapsible so after setting the settings you can
give more space to the content."
"""
import numpy as np
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


def test_groups_are_found_kept_and_read_against_a_label(window, qtbot):
    """Einar: "I need to be able to cluster, and see if I can learn any labels based on my
    clustering, and any cluster's precision and recall with a label"."""
    module = window.map_module.module
    window.panel.mcs.setValue(40)
    module.find_groups()
    qtbot.waitUntil(lambda: module.groups is not None, timeout=120000)
    assert len(module.groups) == window.n
    assert module.colour_box.currentText().startswith("clustering: "), "the map shows the groups"
    assert module.group_table.rowCount() == len(module._group_rows) > 0
    assert module.label_table.rowCount() > 0
    assert "by chance" in module.groups_summary.text()
    assert module.label_name() == module.compare_box.currentText(), \
        "coloured by groups, neighbours are read by the compared label"
    module._group_clicked(0, 0)
    group = module._group_rows[0].group
    assert set(module.selected) == set(np.flatnonzero(module.groups == group).tolist())


def test_dropdowns_drop_below_and_scroll(window):
    module = window.map_module.module
    for box in (module.colour_box, module.compare_box):
        assert "combobox-popup: 0" in box.styleSheet()
        assert box.maxVisibleItems() >= 15


def test_the_old_clusters_tab_points_to_the_module(window):
    tabs = window.panel.findChild(QtWidgets.QTabWidget)
    index = next(i for i in range(tabs.count()) if tabs.tabText(i) == "3 · Clusters")
    texts = [b.text() for b in tabs.widget(index).findChildren(QtWidgets.QPushButton)]
    assert "Open the Map module" in texts


def test_a_deleted_stop_button_does_not_break_starting_a_job(window):
    """The Map tab's stop button was deleted with its page; every job start then raised
    'Internal C++ object already deleted'."""
    window.panel._set_stoppable(True)
    window.panel._set_stoppable(False)
    from shiboken6 import isValid
    assert all(isValid(b) for b in window.panel._stop_btns)


# ------------------------------------------------------------------ the Gene module
def test_the_gene_module_opens_and_reads_a_typed_gene(window, qtbot):
    """Einar chose "synthesis first" for the gene card (instruction 93)."""
    from parascape.landing_install import open_module
    shell = window.gene_module
    assert shell is not None, "the Gene module failed to build"
    open_module(window, "gene")
    assert window.landing.page == "Module"
    module = shell.module
    gid = str(window.nodes.gene_id.iloc[0])
    module.field.setText(gid)
    module._typed()
    assert module.gene == gid
    assert gid in module.known.toHtml()
    sections = list(shell.sections)
    for part in ("Your gene", "What to read", "Advanced", "What is known", "The evidence",
                 "What is not known, and what can be inferred", "Actions",
                 "Console and AI chat"):
        assert part in sections, part
    assert sections.index("What is known") < sections.index("The evidence")
    assert not shell.sections["Advanced"].is_expanded()


def test_unticking_a_kind_of_knowledge_removes_it_from_the_reading(window):
    module = window.gene_module.module
    if module.reading is None:
        module.show(str(window.nodes.gene_id.iloc[0]))
    module.include["fitness"].setChecked(False)
    assert "needs it" not in module.known.toPlainText()
    module.include["fitness"].setChecked(True)


def test_a_list_steps_gene_by_gene(window, tmp_path):
    module = window.gene_module.module
    ids = window.nodes.gene_id.astype(str).tolist()[:3]
    path = tmp_path / "three.csv"
    path.write_text("gene\n" + "\n".join(ids) + "\n")
    module.load_list(str(path))
    assert module.gene == ids[0]
    module._step(1)
    assert module.gene == ids[1]


# ------------------------------------------------------------------ the Network module
def test_the_network_shows_every_gene_and_focuses_a_typed_one(window):
    """Einar: "I want to see the entire network with the gene's region zoomed in, with the
    ability to zoom the network and change the edge definitions of the network" (instr. 87)."""
    from parascape.landing_install import open_module
    shell = window.network_module
    assert shell is not None, "the Network module failed to build"
    open_module(window, "network")
    module = shell.module
    assert len(module.view.xy) == window.n, "every gene is drawn, not only the neighbourhood"
    gid = str(window.nodes.gene_id.iloc[0])
    module.field.setText(gid)
    module._typed()
    assert module.focus == [0]
    assert module.view.dist[0] == 0
    assert gid in module.reading.toHtml()
    for part in ("Your gene", "Colour by", "Which links count", "Advanced", "The network",
                 "Your gene and its links", "Actions", "Console and AI chat"):
        assert part in shell.sections, part


def test_the_kinds_of_link_and_their_strength_define_the_network(window):
    module = window.network_module.module
    before = len(module._links)
    for box in module.layer_boxes.values():
        box.setChecked(False)
    assert len(module._links) == 0
    first = next(iter(module.layer_boxes))
    module.layer_boxes[first].setChecked(True)
    some = len(module._links)
    module.strength.setValue(90)
    module.refresh()
    assert len(module._links) < some
    module.strength.setValue(50)
    for layer, box in module.layer_boxes.items():
        from parascape.network_module import DEFAULT_LAYERS
        box.setChecked(layer in DEFAULT_LAYERS)
    assert len(module._links) == before


def test_the_view_zooms_about_the_pointer_and_travels_to_a_gene(window, qtbot):
    module = window.network_module.module
    view = module.view
    view.resize(800, 600)
    view.fit()
    view.travel_to(0, zoom=5.0)
    qtbot.waitUntil(lambda: abs(view.zoom - 5.0) < 1e-6, timeout=3000)
    assert np.allclose(view.center, view.xy[0])


# ------------------------------------------------------------------ the Infer module
def test_the_infer_module_holds_the_question_controls_moved_not_copied(window):
    from parascape.landing_install import open_module
    shell = window.infer_module
    assert shell is not None, "the Infer module failed to build"
    open_module(window, "infer")
    panel = window.panel
    for name in ("question_choice", "question_detail", "question_table", "question_report",
                 "question_pdf"):
        assert shell.isAncestorOf(getattr(panel, name)), name
    for part in ("Your question", "Your genes", "The question", "The answer",
                 "How far to trust it", "Actions", "Console and AI chat"):
        assert part in shell.sections, part
    assert {"Answer this question", "Check what it may use"} <= set(shell.actions)


def test_the_infer_module_copies_a_methods_text_cited_from_the_registry(window, qtbot):
    """The methods text cites the sources of the columns the question may use, via the clipboard."""
    from PySide6 import QtWidgets
    shell = window.infer_module
    buttons = [b for b in shell.findChildren(QtWidgets.QPushButton) if b.text() == "Copy the methods text"]
    assert len(buttons) == 1 and buttons[0].toolTip()
    QtWidgets.QApplication.clipboard().setText("")
    window.panel.copy_methods_text()
    qtbot.waitUntil(lambda: "References" in QtWidgets.QApplication.clipboard().text()
                    or "refused" in window.panel.question_report.toPlainText(), timeout=60000)
    text = QtWidgets.QApplication.clipboard().text()
    if text:
        assert text.startswith("Analyses were run in ParaScape")
        assert "BibTeX" in text


def test_the_infer_module_saves_an_answer_as_a_bundle(window, tmp_path):
    """The named genes leave marked inferred, with the methods text and a manifest a reload verifies."""
    import pandas as pd
    from types import SimpleNamespace
    from parascape import export_bundle
    ids = window.nodes.gene_id.astype(str).tolist()
    recipe = SimpleNamespace(question="q?", holdout="compartment_best", validation_holdout="ec_number",
                             min_precision=0.3)
    window.panel._question_result = SimpleNamespace(
        recipe=recipe, inference=pd.DataFrame({"gene_id": ids[:2], "label": ["a", "b"]}),
        per_cluster=pd.DataFrame(), recovery=pd.DataFrame({"label": ["a"], "recall": [0.5]}),
        closure=SimpleNamespace(columns=tuple(window.nodes.columns[:5])))
    folder = window.infer_module.module.save_bundle(str(tmp_path))
    manifest = export_bundle.load(folder)
    assert manifest["files"]["named_genes.csv"]["kind"] == "inferred"
    assert manifest["query"]["holdout"] == "compartment_best"
    window.panel._question_result = None
    reopened = window.infer_module.module.open_bundle(folder)
    assert reopened["query"]["question"] == "q?"
    assert len(window.infer_module.module.loaded_bundle) == 2
    (tmp_path / "named_genes.csv").write_text("gene_id,label\nx,y\n")
    assert window.infer_module.module.open_bundle(folder) == {}      # an edited bundle is refused


def test_map_gene_and_network_save_a_bundle_a_reload_verifies(window, tmp_path, qtbot):
    """Every gene-reading module leaves the program the same way: labelled tables, a manifest."""
    from parascape import export_bundle
    from parascape.landing_install import open_module
    gid = str(window.nodes.gene_id.iloc[0])
    gene = window.gene_module.module
    gene.show(gid)
    folder = gene.save_bundle(str(tmp_path / "gene"))
    manifest = export_bundle.load(folder)
    assert manifest["query"]["gene"] == gid
    other = str(window.nodes.gene_id.iloc[1])
    gene.show(other)
    assert gene.open_bundle(folder)["query"]["gene"] == gid and gene.gene == gid
    assert manifest["files"]["what_is_measured.csv"]["kind"] == "measured"
    open_module(window, "map")
    qtbot.wait(200)
    m = window.map_module.module
    if getattr(window, "xyz", None) is not None and len(window.xyz):
        manifest = export_bundle.load(m.save_bundle(str(tmp_path / "map")))
        assert manifest["files"]["map_coordinates.csv"]["kind"] == "computed"
        assert "map.svg" in manifest["files"]
    net = window.network_module.module
    if net.focus and net._links is not None:
        manifest = export_bundle.load(net.save_bundle(str(tmp_path / "net")))
        assert manifest["files"]["links.csv"]["kind"] == "measured"


def test_the_infer_data_step_holds_out_the_labels_kind_and_lets_you_leave_out_others(window):
    """Instruction 74: the label's own kind of evidence is shown held out; unticking another kind
    reaches the recipe, and closure removes those columns with the reason."""
    from dataclasses import replace
    from parascape.recipes import close
    module = window.infer_module.module
    module.refresh_uses()
    boxes = module.uses_boxes
    assert boxes, "the question's inputs hold no categorised evidence"
    free = [c for c, b in boxes.items() if b.isEnabled()]
    assert free
    boxes[free[0]].setChecked(False)
    assert window.panel.question_leave_out == (free[0],)
    recipe = window.panel._question_recipe(window.panel.current_question())
    assert recipe.leave_out_categories == (free[0],)
    closure = close(window.nodes, recipe)
    assert any("left out by you" in why for why in closure.excluded.values())
    module.refresh_uses()
    assert window.panel.question_leave_out == ()


def test_the_gene_module_lists_what_each_method_said_held_out(window):
    module = window.gene_module.module
    module.show(str(window.nodes.gene_id.iloc[0]))
    assert "What each method said about it, held out" in window.gene_module.sections
    assert "methods were tested on this gene" in module.methods.toPlainText() \
        or "could not be read" in module.methods.toPlainText()


def test_your_genes_are_read_against_the_answer(window):
    import pandas as pd
    module = window.infer_module.module
    ids = window.nodes.gene_id.astype(str).tolist()
    module.field.setText(f"{ids[0]}, {ids[1]}")
    module._typed()
    assert module.mine == ids[:2]

    class Result:
        inference = pd.DataFrame({"gene_id": [ids[0], ids[5]], "label": ["dense granules", "x"],
                                  "control_agrees": [True, False]})
    window.panel._question_result = Result()
    module._after_answer()
    text = module.answer_note.text()
    assert "2</b> genes named" in text and "1</b> are named" in text and "dense granules" in text
    window.panel._question_result = None


def test_the_gene_module_keeps_a_map_suggestion_and_forgets_it(window, tmp_path, monkeypatch):
    from parascape import kept_claims as KC
    monkeypatch.setenv("PARASCAPE_KEPT_CLAIMS", str(tmp_path / "kept.jsonl"))
    module = window.gene_module.module
    gid = str(window.nodes.gene_id.iloc[0])
    module.show(gid)
    claim = module.keep_map_suggestion()
    if claim is None:
        pytest.skip("no labelled neighbours on the test map")
    assert "Kept by you" in module.inferred.toHtml()
    from PySide6.QtCore import QUrl
    module._inferred_link(QUrl(f"forget:{claim.claim_id}"))
    assert KC.claims(gid).empty
