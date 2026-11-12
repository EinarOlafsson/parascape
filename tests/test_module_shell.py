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


# -------------------------------------------- the gene card to the mould (98, WP8)
def test_the_gene_card_is_laid_out_in_the_order_the_mould_asks_for(window):
    """Masthead, then the subject, the lens, the choice that changes the meaning,
    Advanced folded and Save folded on the left; the thing itself first on the right,
    then the synthesis, then what can be inferred, then the next step."""
    shell = window.gene_module
    sections = list(shell.sections)
    for part in ("Your gene", "Where it is, read by", "What to read", "Advanced", "Save",
                 "Where it is in the parasite", "What is known",
                 "What you can do with this gene next", "Actions"):
        assert part in sections, part
    assert sections.index("Your gene") < sections.index("Where it is, read by")
    assert sections.index("Where it is, read by") < sections.index("What to read")
    assert sections.index("What to read") < sections.index("Advanced") < sections.index("Save")
    assert sections.index("Where it is in the parasite") < sections.index("What is known")
    assert sections.index("What is known") \
        < sections.index("What is not known, and what can be inferred")
    assert sections.index("What is not known, and what can be inferred") \
        < sections.index("What you can do with this gene next")
    assert sections.index("What you can do with this gene next") < sections.index("Actions")
    assert not shell.sections["Advanced"].is_expanded()
    assert not shell.sections["Save"].is_expanded()
    assert shell.sections["Console and AI chat"].shut, "the console must start folded"
    assert shell.title_label.text() == "Gene card"


def test_the_parasite_schematic_survived_the_restructure(window):
    """Einar named the parasite schematic as a thing the overhaul must keep. It is the
    thing the card is about, so it is the first section of the runtime column, and it
    shows the real reading rather than a sample."""
    from parascape.celldiagram import CellDiagram
    module = window.gene_module.module
    assert isinstance(module.parasite, CellDiagram), "no parasite drawing on the gene card"
    placed = window.nodes[window.nodes["compartment"].astype(str) == "dense granules"]
    gene = str(placed.gene_id.iloc[0])
    module.show(gene)
    place = module.placement()
    assert place.compartment == "dense granules"
    assert "dense granules" in module.parasite_note.text()
    assert module.parasite.selected == "dense granules"
    # The window's own palette, never a second one.
    assert module.parasite.color_of == window.color_of
    assert module.parasite.fills(), "the organelle is not filled in"


def test_a_gene_with_no_measured_compartment_says_unknown_rather_than_nothing(window):
    unplaced = window.nodes[window.nodes["compartment"].astype(str) == "unassigned"]
    module = window.gene_module.module
    module.show(str(unplaced.gene_id.iloc[0]))
    assert "unknown" in module.parasite_note.text()
    assert not module.parasite.fills(), "an unknown place must colour nothing"


def test_the_parasite_is_read_by_the_same_label_the_window_is_coloured_by(window):
    """One source of truth: the card and the map must not show two different
    compartments for one gene."""
    module = window.gene_module.module
    columns = module._placement_columns()
    assert "compartment" in columns
    other = next(c for c in columns if c != "compartment")
    window.category_box.setCurrentText(other)
    assert module.place_box.currentText() == other
    module.place_box.setCurrentText("compartment")
    assert window.category_box.currentText() == "compartment"


def test_the_trust_line_reports_against_chance_with_n(window):
    module = window.gene_module.module
    module.show(str(window.nodes.gene_id.iloc[0]))
    trust = module.parasite_trust.text()
    assert "measured for" in trust and "commonest answer" in trust
    assert "by chance" in trust or "no evidence at all" in trust
    assert f"{module.placement().measured:,}" in trust


def test_the_card_offers_every_module_that_reads_a_gene_pre_filled(window):
    """Einar: "so all the modules that take 1 gene and do something should come
    downstream of the gene card module"."""
    from parascape import ports as P
    module = window.gene_module.module
    gene = str(window.nodes.gene_id.iloc[0])
    module.show(gene)
    offered = module.next_steps.offers
    assert {o.module for o in offered} == set(P.consumers_of(P.GENE)) - {"gene"}
    assert all(o.value == gene for o in offered)
    blocked = [o for o in offered if not o.ok]
    assert blocked and all(o.blocked and o.fix for o in blocked)
    for offer in offered:
        assert module.next_steps.buttons[offer.module].isEnabled() is offer.ok


def test_clearing_the_card_empties_it_and_the_offers(window):
    """Through the button, not the method: `clicked` carries a checked flag, and a slot
    taking no arguments fails only when somebody presses it."""
    module = window.gene_module.module
    module.show(str(window.nodes.gene_id.iloc[0]))
    clear = next(b for b in module.shell.settings.findChildren(QtWidgets.QPushButton)
                 if b.text() == "Clear")
    clear.click()
    assert module.gene == "" and module.reading is None
    assert module.field.text() == "" and module.matches.count() == 0
    assert "Type a gene" in module.parasite_note.text()
    assert all(o.value == "" for o in module.next_steps.offers)


def test_a_next_step_opens_the_module_with_the_gene_in_it(window):
    """The pre-filling half of an offer: a next step that opened an empty module would
    leave the gene to be typed a second time."""
    module = window.gene_module.module
    gene = str(window.nodes.gene_id.iloc[1])
    module.show(gene)
    module.next_steps.buttons["network"].click()
    assert window.landing.page == "Module"
    assert window.search.text() == gene, "the window was not moved onto the gene"


def test_typing_a_gene_in_the_navigator_lands_on_the_gene_card(window):
    """Einar: "now when i put in a gene and pressed the modual that came up for the gene,
    i was brought to the old version of parascape". It opened the evidence dock."""
    from parascape.landing_install import open_item
    gene = str(window.nodes.gene_id.iloc[2])
    open_item(window, "gene", gene)
    assert window.landing.page == "Module"
    assert window.landing.stack.currentWidget() is window.gene_module
    assert window.gene_module.module.gene == gene


def test_the_gene_route_never_reaches_the_pre_navigator_workflow_window(window):
    """WP10 retires `workflows.py`; WP8 only has to stop the gene route reaching it."""
    from parascape.landing_install import open_item
    open_item(window, "gene", str(window.nodes.gene_id.iloc[3]))
    assert not hasattr(window, "workflows_dialog"), \
        "opening a gene built the pre-navigator Workflows window"


def test_the_gene_module_lists_links_to_other_organisms_in_their_own_spaces(window):
    """Instruction 64.24: a link's other end is resolved in its own space, never guessed."""
    from parascape import bridges
    module = window.gene_module.module
    sections = list(window.gene_module.sections)
    assert sections.index("The evidence") < sections.index("Links to other organisms")
    ortho = bridges.orthology_links()
    if len(ortho):
        gene = str(ortho.iloc[0]["source_id"])
        module.show(gene)
        html = module.links.toHtml()
        assert "orthology (shared OrthoMCL group)" in html and "open:" in html
        assert any(c["space"] == "Pf" for c in module._counterparts.values())
    lit = bridges.host_links()
    lit = lit[lit["target_organism"] == bridges.HOST_UNRECORDED]
    if len(lit) and str(lit.iloc[0]["source_id"]) in set(window.nodes.gene_id.astype(str)):
        module.show(str(lit.iloc[0]["source_id"]))
        assert "does not say which host species" in module.links.toPlainText()


def test_the_species_menu_lists_host_spaces_with_their_state(window, tmp_path, monkeypatch):
    """Instruction 64.23: host gene packs are listed with an explicit state, never as calibrated."""
    from parascape import paths, space_catalogue
    actions = window.host_space_actions
    assert set(actions) == set(space_catalogue.HOST_PACK_SPACES)
    for act in actions.values():
        assert "not calibrated" in act.statusTip()
    if space_catalogue.local_archives("Mm"):
        monkeypatch.setenv(paths.ENV_STATE, str(tmp_path / "state"))
        assert window.install_host_space("Mm")
        assert "installed" in actions["Mm"].text() and not actions["Mm"].isEnabled()


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


def test_the_infer_data_step_says_how_many_columns_circularity_removed(window):
    """Instruction 73: the kinds are not enough -- the COUNT and the columns themselves are shown.

    "place is held out" does not tell a reader whether that cost two columns or two hundred, and the answer
    changes how much to trust what is left. The names and reasons go in the tooltip, because the note must stay
    one line.
    """
    module = window.infer_module.module
    module.refresh_uses()
    note = module.uses_note.text()
    assert "left out for that reason" in note
    circular = module._circular
    assert isinstance(circular, dict)
    assert str(len(circular)) in note, "the note must carry the number of columns held out"
    tooltip = module.uses_note.toolTip()
    if circular:
        first = sorted(circular)[0]
        assert first in tooltip and circular[first][:12] in tooltip
    else:
        assert "No column shares" in tooltip

def test_the_gene_module_lists_what_each_method_said_held_out(window):
    module = window.gene_module.module
    module.show(str(window.nodes.gene_id.iloc[0]))
    assert "What each method said about it, held out" in window.gene_module.sections
    assert "strategy:feature_knn" in module.methods.toHtml()
    assert module.open_strategy("feature_knn")
    assert window.strategy_panel.current.key == "feature_knn"
    assert not module.open_strategy("no_such_strategy")
    assert "methods were tested on this gene" in module.methods.toPlainText() \
        or "could not be read" in module.methods.toPlainText()


@pytest.mark.parametrize("width", [1024, 1366, 1920])
def test_the_gene_module_reads_a_measured_and_a_barely_measured_gene_at_each_width(window, qtbot, width):
    """Instruction 64.36: real widgets at supported widths, including the missing-data case."""
    from PySide6 import QtWidgets
    nodes = window.nodes
    numeric = nodes.select_dtypes("number")
    sparse = str(nodes.gene_id.iloc[int(numeric.notna().sum(axis=1).to_numpy().argmin())])
    rich = str(nodes.gene_id.iloc[int(numeric.notna().sum(axis=1).to_numpy().argmax())])
    window.resize(width, 900)
    window.show()
    try:
        module = window.gene_module.module
        for gene in (rich, sparse):
            module.show(gene)
            QtWidgets.QApplication.processEvents()
            assert gene in module.known.toPlainText()
            assert module.unusual.toPlainText().strip()          # either the unusual list or why not
            assert module.methods.toPlainText().strip()
        assert window.gene_module.width() <= width
    finally:
        window.hide()


def test_going_home_and_opening_another_gene_never_shows_the_last_one(window):
    """Instruction 64.36: back and search keep a valid context; no stale result from another entity."""
    from PySide6 import QtWidgets
    from parascape.gene_module import open_in
    from parascape.landing_install import go_home
    a, b = (str(g) for g in window.nodes.gene_id.iloc[[10, 20]])
    open_in(window, a)
    QtWidgets.QApplication.processEvents()
    go_home(window)
    open_in(window, b)
    QtWidgets.QApplication.processEvents()
    module = window.gene_module.module
    assert module.gene == b
    text = module.known.toPlainText()
    assert b in text and a not in text
    assert module.reading.gene == b


def test_a_gene_with_a_host_target_shows_its_host_counterpart(window):
    """Instruction 64.36: the host side is reachable from the gene -- a parasite protein with a recorded
    host target names it in what the gene binds."""
    import pandas as pd
    from parascape.paths import data_dir
    import os
    h = pd.read_parquet(os.path.join(data_dir(), "host_interactions.parquet"))
    h = h[h["gene_id"].astype(str).isin(set(window.nodes.gene_id.astype(str)))]
    if not len(h):
        pytest.skip("no host interaction for a gene of this table")
    gene, target = str(h.iloc[0]["gene_id"]), str(h.iloc[0]["host_target"])
    module = window.gene_module.module
    module.show(gene)
    assert target in module.known.toPlainText() or target in module.reading.partners.to_string()


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


def test_a_value_scorecard_states_its_units(window):
    """Instruction 64.11: MAE and RMSE are read in the unit of the predicted column, so the card says it."""
    from parascape import strategies as S
    from parascape import scorecard as SC
    panel = window.strategy_panel
    value = next(s for s in S.catalog() if s.task == SC.T_VALUES)
    html = panel.scorecard_html(value)
    assert "MAE and RMSE are in the units of" in html


def test_the_infer_module_answers_with_a_chosen_method_and_shows_its_grade(window, monkeypatch):
    """Instruction 74, Infer step 3: the same question answered by a harness method, graded on held-out genes."""
    import pandas as pd
    from parascape import harness_answer as HA
    shell = window.infer_module
    module = shell.module
    sections = list(shell.sections)
    assert sections.index("What the answer may use") < sections.index("How to answer") < sections.index("Your genes")
    keys = [module.method_choice.itemData(i) for i in range(module.method_choice.count())]
    assert keys[0] == "map" and set(HA.METHODS) <= set(keys)
    ids = window.nodes.gene_id.astype(str).tolist()
    fake = {"refusal": "", "grades": pd.DataFrame([dict(label="x = yes", labelled=500, frozen_auc=0.81,
                                                        average_precision=0.4, share_of_yes=0.1, may_fill=False)]),
            "calls": pd.DataFrame({"gene_id": ids[:3], "answer": ["x = yes"] * 3, "probability": [0.9, 0.6, 0.2],
                                   "may_fill": [False] * 3})}
    monkeypatch.setattr(HA, "answer", lambda *a, **k: fake)
    monkeypatch.setattr(window.panel, "_run", lambda job, done, name="", on_error=None: done(job(lambda *_: None)))
    module.method_choice.setCurrentIndex(keys.index("hist_gbt"))
    module.answer_question()
    report = window.panel.question_report.toPlainText()
    assert "graded on genes it never saw" in report and "0.81" in report and "rank genes only" in report
    assert window.panel.question_table.rowCount() == 3
    module.method_choice.setCurrentIndex(0)
