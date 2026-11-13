"""The GENE LIST CARD in a real window (instruction 98, WP9).

Einar's words, which are the acceptance criteria: *"if i type in a list of genes, the same
logic but this time not gene card instead a gene list module and after clicking it all the
moduals that come downstream of gene lists."*

The same shape as the gene card tests in `test_module_shell.py`, and deliberately a file of
its own: the window fixture there is module-scoped and already builds every module, so a
second file costs one more window and keeps a long file from growing by a third.
"""
import pytest
from PySide6 import QtWidgets

from parascape import modules as MOD


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


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


def a_list(window, n=24):
    """A list worth reading: genes that share a compartment, as a real list would."""
    placed = window.nodes[window.nodes["compartment"].astype(str) == "dense granules"]
    return [str(g) for g in placed.gene_id.astype(str).head(n)]


def test_the_gene_list_card_is_laid_out_in_the_order_the_mould_asks_for(window):
    """Masthead, then the subject, the lens, the choice that changes the meaning, Advanced
    folded and Save folded on the left; the thing itself first on the right, then where the
    genes are, then how much is measured, then what is unknown, then the next step."""
    shell = window.list_module
    assert shell is not None, "the gene list module failed to build"
    sections = list(shell.sections)
    for part in ("Your genes", "Where they are, read by", "What to read", "Advanced", "Save",
                 "What the genes of the list share",
                 "Where the genes of the list are in the parasite",
                 "How much of the list is measured",
                 "What is unknown for most of the list, and what can be inferred",
                 "The genes of the list", "What you can do with this list next", "Actions"):
        assert part in sections, part
    assert sections.index("Your genes") < sections.index("Where they are, read by")
    assert sections.index("Where they are, read by") < sections.index("What to read")
    assert sections.index("What to read") < sections.index("Advanced") < sections.index("Save")
    assert sections.index("What the genes of the list share") \
        < sections.index("Where the genes of the list are in the parasite")
    assert sections.index("How much of the list is measured") \
        < sections.index("What is unknown for most of the list, and what can be inferred")
    assert sections.index("What is unknown for most of the list, and what can be inferred") \
        < sections.index("What you can do with this list next") < sections.index("Actions")
    assert not shell.sections["Advanced"].is_expanded()
    assert not shell.sections["Save"].is_expanded()
    assert shell.sections["The genes of the list"].shut, "the identifiers are not the card"
    assert shell.sections["Console and AI chat"].shut, "the console must start folded"
    # The masthead is the catalogue row since instruction 98, WP11: it said "Gene list card"
    # while the tile said "Gene list", which is one fact in two tables and the bug this whole
    # instruction exists to remove.
    assert shell.title_label.text() == "Gene list"
    assert shell.title_label.text() == MOD.get("list").name


def test_a_list_card_is_not_a_table_of_identifiers(window):
    """The mould: make the information about a list accessible and SYNTHESISED. What opens
    first is what the list SHARES, read against every other gene; the identifiers are the
    folded section at the bottom."""
    module = window.list_module.module
    module.set_genes(a_list(window), "a test")
    assert module.shared, "nothing was found that the list shares"
    assert module.shared_view.rowCount() == len(module.shared)
    for row in range(module.shared_view.rowCount()):
        for column in range(module.shared_view.columnCount()):
            assert module.shared_view.item(row, column).toolTip()
    lead = module.shared_lead.text()
    assert "read against every other gene" in lead
    assert str(module.resolved.count) in lead


def test_every_share_the_card_shows_carries_what_it_is_a_share_of(window):
    """Against chance, with n, as whole numbers -- in the table, not only in the tooltip."""
    module = window.list_module.module
    module.set_genes(a_list(window), "a test")
    top = module.shared[0]
    assert f"{top.in_list:,} of {top.list_measured:,}" in module.shared_view.item(0, 2).text()
    assert f"of {top.outside_measured:,}" in module.shared_view.item(0, 3).text()
    assert "random lists" in module.shared_view.item(0, 4).text()


def test_the_parasite_schematic_is_on_the_list_card_too(window):
    """Einar named the schematic as a thing the overhaul must keep. A list of genes sits
    somewhere in the cell exactly as one gene does."""
    from parascape.celldiagram import CellDiagram
    module = window.list_module.module
    module.set_genes(a_list(window), "a test")
    assert isinstance(module.parasite, CellDiagram), "no parasite drawing on the list card"
    assert module.spread().commonest == "dense granules"
    assert "dense granules" in module.parasite_note.text()
    assert module.parasite.selected == "dense granules"
    assert module.parasite.color_of == window.color_of, "one palette, never a second"
    assert module.parasite.fills(), "the organelle is not filled in"
    assert "of the" in module.parasite_trust.text()


def test_a_list_with_no_measured_compartment_says_unknown_rather_than_nothing(window):
    unplaced = window.nodes[window.nodes["compartment"].astype(str) == "unassigned"]
    module = window.list_module.module
    module.set_genes([str(g) for g in unplaced.gene_id.astype(str).head(12)], "a test")
    assert "unknown" in module.parasite_note.text()
    assert not module.parasite.fills(), "an unknown place must colour nothing"


def test_how_much_of_the_list_is_measured_is_said_against_the_other_genes(window):
    module = window.list_module.module
    module.set_genes(a_list(window), "a test")
    wanted = [c for c in module.coverage if module.wanted(c.category)]
    assert module.coverage_view.rowCount() == len(wanted) > 1
    lead = module.coverage_lead.text()
    assert "kinds of knowledge" in lead and "unknown for most of them" in lead
    # Least measured first: the emptiest kinds decide what to do next.
    shares = [int(module.coverage_view.item(r, 1).text().split("(")[1].split(" %")[0])
              for r in range(module.coverage_view.rowCount())]
    assert shares == sorted(shares)
    assert "%" in module.coverage_view.item(0, 2).text(), "the other genes, for comparison"


def test_what_is_unknown_for_most_of_the_list_points_at_inference(window):
    module = window.list_module.module
    module.set_genes(a_list(window), "a test")
    said = module.unknown_view.toPlainText()
    assert "unknown" in said or "no gap here" in said
    assert "Predict what is missing" in module.shell.actions


def test_the_card_offers_every_module_that_reads_a_list_pre_filled(window):
    """Einar: "after clicking it all the moduals that come downstream of gene lists"."""
    from parascape import ports as P
    module = window.list_module.module
    genes = a_list(window)
    module.set_genes(genes, "a test")
    offered = module.next_steps.offers
    assert {o.module for o in offered} == set(P.consumers_of(P.GENE_LIST)) - {"list"}
    assert all(o.value == ",".join(genes) for o in offered)
    blocked = [o for o in offered if not o.ok]
    assert blocked, "Infer needs a question, so one offer must be greyed with its reason"
    assert all(o.blocked and o.fix for o in blocked)
    for offer in offered:
        assert module.next_steps.buttons[offer.module].isEnabled() is offer.ok
    # Said as a count, never as twenty-four accessions.
    assert f"{len(genes)} genes" in module.next_steps.lead.text()
    assert not any(g in module.next_steps.lead.text() for g in genes)


def test_clearing_the_list_card_empties_it_and_the_offers(window):
    """Through the button, not the method: `clicked` carries a checked flag, and a slot
    taking no arguments fails only when somebody presses it."""
    module = window.list_module.module
    module.set_genes(a_list(window), "a test")
    clear = next(b for b in module.shell.settings.findChildren(QtWidgets.QPushButton)
                 if b.text() == "Clear")
    clear.click()
    assert module.genes == [] and module.shared == ()
    assert module.box.toPlainText() == "" and module.members.rowCount() == 0
    assert "Type or choose a list" in module.parasite_note.text()
    assert all(o.value == "" for o in module.next_steps.offers)


def test_a_typed_list_and_a_list_file_both_reach_the_card(window, tmp_path):
    module = window.list_module.module
    genes = a_list(window, 6)
    module.box.setPlainText(", ".join(genes))
    module._typed()
    assert module.genes == genes
    path = tmp_path / "hits.csv"
    path.write_text("gene\n" + "\n".join(genes[:4] + ["NOT_A_GENE"]) + "\n")
    module.load_list(str(path))
    assert module.genes == genes[:4]
    assert "hits.csv" in module.list_note.text()


def test_the_lens_of_the_list_card_is_synced_with_the_window(window):
    """One source of truth: the card and the map must not mean two different things."""
    module = window.list_module.module
    columns = module._placement_columns()
    assert "compartment" in columns
    other = next(c for c in columns if c != "compartment")
    window.category_box.setCurrentText(other)
    assert module.place_box.currentText() == other
    module.place_box.setCurrentText("compartment")
    assert window.category_box.currentText() == "compartment"


def test_a_next_step_from_the_list_card_opens_the_module_with_the_list_in_it(window):
    """The pre-filling half of an offer. The gene card steps THROUGH a list, which its own
    port has always said it does."""
    module = window.list_module.module
    genes = a_list(window, 5)
    module.set_genes(genes, "a test")
    module.next_steps.buttons["gene"].click()
    assert window.landing.page == "Module"
    assert window.gene_module.module.genes == genes
    assert window.gene_module.module.gene == genes[0]


def test_one_gene_of_the_list_opens_on_its_own(window):
    module = window.list_module.module
    genes = a_list(window, 5)
    module.set_genes(genes, "a test")
    assert module.members.rowCount() == 5
    module._open_member(2)
    assert window.gene_module.module.gene == genes[2]


def test_carrying_a_list_into_the_window_gates_it(window):
    """The gate is what the application already means by "this set of genes is in hand":
    the map marks it, the strategy panel fills a gene box from it, and Export writes it."""
    from parascape import landing_install as LI
    from parascape import ports as P
    genes = a_list(window, 7)
    assert LI.carry_subject(window, P.GENE_LIST, ",".join(genes))
    gated = [str(g) for g in window.nodes.gene_id.astype(str).to_numpy()[window.gated]]
    assert set(gated) == set(genes)
    window.clear_gate()


def test_a_list_saved_as_a_bundle_can_be_read_back(window, tmp_path):
    module = window.list_module.module
    genes = a_list(window, 6)
    module.set_genes(genes, "a test")
    folder = module.save_bundle(str(tmp_path / "bundle"))
    assert folder
    module.clear()
    assert module.genes == []
    module.open_bundle(folder)
    assert module.genes == genes
