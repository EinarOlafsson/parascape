"""The tree of what ParaScape can do, and the rail that draws it (instruction 98, WP7).

Einar: "there should be a tree of things that parascape can do, and this progression
should guide the user to the modual they want. just like in spacr."

Every test here pins something that has to stay true for the tree to be computed rather
than written down, because the failure spaCR records as bug #16i -- two tables of the same
fact, drifting apart silently -- is exactly what a hand-written tree would be.
"""
import pytest
from PySide6 import QtWidgets

from parascape import modules as MOD
from parascape import ports as PORTS
from parascape import tree as T


# ------------------------------------------------------------------ what it computes


def test_the_six_headings_are_einars_own_in_his_order():
    """His decision of 2026-11-12: ParaScape keeps its own headings rather than
    mirroring the sections of spaCR, which describe an imaging pipeline this has not."""
    assert T.SECTIONS == ("Start", "Look up", "Explore", "Ask", "Bring your own", "Data")


def test_the_five_bands_come_from_the_catalogue_rather_than_a_second_list():
    assert T.SECTIONS[1:] == tuple(MOD.SECTIONS)


def test_start_is_not_a_band_of_the_catalogue():
    """It holds no module -- Home is not one and a subject is not one -- and a declared
    band must have modules (`test_no_section_is_empty`). spaCR keeps SECTION_HELP out of
    its SECTION_ORDER for the same reason and records that putting it in was reverted."""
    assert T.START not in MOD.SECTIONS


def test_every_module_is_filed_exactly_once():
    filed = T.modules_in()
    assert sorted(filed) == sorted(m.key for m in MOD.MODULES)
    assert len(filed) == len(set(filed))


def test_a_module_is_filed_under_the_heading_its_catalogue_row_declares():
    for module in MOD.MODULES:
        assert T.section_of(module.key) == module.band, module.key


def test_every_declared_view_is_in_the_tree():
    """"Include the 32 declared views so nothing behind the navigator is unreachable."
    A view missing here is a panel the rail cannot name, which is the complaint."""
    assert sorted(T.views_in()) == sorted((m.key, v.name) for m, v in MOD.views())
    assert len(T.views_in()) == 32


def test_a_module_leads_to_what_reads_what_it_hands_on():
    """The next steps are `modules.downstream`, which is computed from the declared
    subjects -- the same function the gene card and the module network read."""
    nodes = T.tree()
    for _depth, node in T.walk(nodes):
        if node.kind != "module":
            continue
        steps = [c.key for c in node.children if c.kind == "step"]
        assert steps == list(MOD.downstream(node.key)), node.key


def test_a_module_row_says_what_it_takes_and_what_it_hands_on():
    gene = T.find("gene")
    assert "one gene" in gene.note
    assert MOD.takes_words("gene") in gene.note
    assert MOD.hands_on_words("gene") in gene.note


def test_the_progression_is_over_the_subject_you_have():
    """ParaScape has no single pipeline to progress along, so the progression is the
    subject in hand and what becomes possible once you hold one."""
    brought = T.subjects_you_can_bring()
    assert PORTS.GENE in brought and PORTS.GENE_LIST in brought
    for kind in brought:
        assert kind in PORTS.SUBJECT_SOURCES, kind
        assert MOD.downstream_of_subject(kind), kind


def test_nothing_offers_a_compound_because_nothing_reads_one():
    """Declared in ports and consumed by no module, which is the inventory finding of
    WP1 rather than an oversight. A heading for it would be a dead row."""
    assert PORTS.COMPOUND not in T.subjects_you_can_bring()


def test_a_subject_offers_the_modules_that_read_it_ready_ones_first():
    node = next(n for _d, n in T.walk() if n.kind == "subject" and n.key == PORTS.GENE)
    assert [c.key for c in node.children] == list(MOD.downstream_of_subject(PORTS.GENE))
    blocked = [i for i, c in enumerate(node.children) if c.blocked]
    ready = [i for i, c in enumerate(node.children) if not c.blocked]
    assert not ready or not blocked or min(blocked) > max(ready)


def test_a_module_under_a_subject_is_a_reference_not_a_second_filing():
    """`test_module_catalogue.py::test_no_module_is_in_two_sections` is untouched and
    still passes, because the catalogue has not changed. The rows under Start point at
    the modules filed elsewhere and say so, which is what keeps both true at once."""
    for _depth, node in T.walk():
        if node.kind == "step":
            assert node.reference
            assert T.section_of(node.key) != T.START


def test_a_module_that_cannot_open_keeps_its_row_with_the_reason():
    """Einar asked the tiles to be an honest inventory rather than a menu of promises."""
    structure = T.find("structure")
    assert structure is not None and structure.blocked
    assert "structure view" in structure.blocked


def test_a_window_reason_beats_the_catalogue_reason():
    nodes = T.tree({"gene": "the evidence view is not available in this window"})
    assert "not available in this window" in T.find("gene", nodes).blocked


def test_the_tree_is_four_levels_and_no_deeper():
    assert max(depth for depth, _n in T.walk()) == 2


def test_find_never_returns_a_reference():
    assert T.find("partners").kind == "module" and not T.find("partners").reference
    assert T.find("nothing-is-called-this") is None


# ------------------------------------------------------------------------ the rail


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture
def rail(qapp):
    from parascape import iconset
    from parascape.module_tree import ModuleTree
    widget = ModuleTree(icon_for=lambda key: iconset.module_icon(key))
    yield widget
    widget.deleteLater()


def test_the_rail_draws_one_row_for_every_row_of_the_tree(rail):
    assert len(rail._entries) == len(T.walk())


def test_the_rail_lists_home_and_every_module_once(rail):
    keys = [row.key for row in rail.rows()]
    assert keys[0] == "home"
    assert sorted(keys[1:]) == sorted(m.key for m in MOD.MODULES)


def test_every_declared_view_has_a_row_in_the_rail(rail):
    assert len(rail.view_rows()) == 32
    named = {(row.node.owner, row.node.key) for row in rail.view_rows()}
    assert named == {(m.key, v.name) for m, v in MOD.views()}


def shown(rows) -> list:
    """The rows of the rail that are not folded away.

    `isVisible` is False for every child until the rail itself is on screen, and these
    tests build the rail alone, so the question asked is whether anything HID the row --
    which is `isHidden`, the same distinction `landing_install.fit_landing` records.
    """
    return [row for row in rows if not row.isHidden()]


def test_only_the_progression_is_open_when_the_rail_appears(rail):
    """Six headings collapsed told nobody anything. Start open is the one that says
    what to do: here is what you can have, and here is what it lets you do."""
    names = {w.text() for w in rail.visible_rows() if hasattr(w, "text")}
    assert "Home" in names and "One gene" in names
    assert not shown([row for row in rail.rows() if row.key == "map"])


def test_a_heading_opens_and_shuts(rail):
    assert rail.toggle_section("Look up") is True
    assert shown([row for row in rail.rows() if row.key == "gene"])
    assert rail.toggle_section("Look up") is False
    assert not shown([row for row in rail.rows() if row.key == "gene"])


def test_clicking_a_module_opens_it_and_shows_what_is_behind_it(rail):
    opened = []
    rail.nav_selected.connect(opened.append)
    rail.toggle_section("Look up")
    next(r for r in rail.rows() if r.key == "expression").click()
    assert opened == ["expression"]
    assert shown([v for v in rail.view_rows() if v.node.key == "Coexpression"])
    assert shown([s for s in rail.step_rows() if s.node.owner == "module:expression"])


def test_only_one_module_of_a_heading_stays_unfolded(rail):
    """Nine modules each unfolding ten views put the rail past any screen. Opening one
    folds the last, so a heading shows its modules and the inside of at most one."""
    rail.toggle_section("Look up")
    next(r for r in rail.rows() if r.key == "expression").click()
    next(r for r in rail.rows() if r.key == "gene").click()
    assert not shown([v for v in rail.view_rows() if v.node.owner == "expression"])
    assert shown([v for v in rail.view_rows() if v.node.owner == "gene"])


def test_shutting_a_heading_forgets_which_module_was_unfolded(rail):
    rail.toggle_section("Look up")
    next(r for r in rail.rows() if r.key == "gene").click()
    rail.toggle_section("Look up")
    rail.toggle_section("Look up")
    assert not shown([v for v in rail.view_rows() if v.node.owner == "gene"])


def test_clicking_a_view_names_the_module_it_belongs_to(rail):
    chosen = []
    rail.view_selected.connect(lambda key, name: chosen.append((key, name)))
    rail.open_to("expression")
    next(v for v in rail.view_rows() if v.node.key == "Coexpression").click()
    assert chosen == [("expression", "Coexpression")]


def test_clicking_a_next_step_opens_the_module_it_leads_to(rail):
    opened = []
    rail.nav_selected.connect(opened.append)
    rail.open_to("gene")
    step = next(s for s in rail.step_rows() if s.node.owner == "module:gene")
    step.click()
    assert opened[-1] == step.node.key
    assert step.node.key in MOD.downstream("gene")


def test_a_step_under_a_subject_is_not_a_step_under_the_module_of_that_name(rail):
    """The module `gene` and the subject `gene` are both called "gene", and the second
    reads as the first with everything else that reads a gene added to it."""
    under_module = {s.node.key for s in rail.step_rows()
                    if s.node.owner == "module:gene"}
    under_subject = {s.node.key for s in rail.step_rows()
                     if s.node.owner == "subject:gene"}
    assert under_module == set(MOD.downstream("gene"))
    assert under_subject == set(MOD.downstream_of_subject("gene"))
    assert under_module != under_subject


def test_opening_a_module_from_elsewhere_opens_the_tree_down_to_it(rail):
    rail.open_to("scorecard")
    assert shown([r for r in rail.rows() if r.key == "scorecard"])
    assert shown([v for v in rail.view_rows() if v.node.owner == "scorecard"])


def test_a_module_that_folds_nothing_still_opens_its_heading(rail):
    """Export hands nothing on and has no view, so there is nothing to unfold under it
    -- and the heading above it was being opened and never drawn."""
    rail.open_to("export")
    assert shown([r for r in rail.rows() if r.key == "export"])


def test_every_row_says_what_it_does_on_hover(rail):
    """Tooltip on every control, and in the rail the tooltip is the only place what a
    module takes and what follows it can be read without opening anything."""
    for row in rail.all_rows():
        assert row.toolTip().strip(), row.key
        assert row.toolTip().startswith(row.node.name), row.key
    for heading in rail._headers.values():
        assert heading.toolTip().strip(), heading.section


def test_a_name_too_long_for_the_rail_is_still_readable_on_hover(rail):
    """The rail is 220 to 320 px wide and "Hide labels and see if they come back" does
    not fit in it at any width the dock allows, so the tooltip opens with the name."""
    rail.open_to("scorecard")
    row = next(v for v in rail.view_rows()
               if v.node.key == "Hide labels and see if they come back")
    assert row.toolTip().startswith("Hide labels and see if they come back")
    assert row.accessibleName() == row.node.name


def test_no_row_carries_a_cross(rail):
    """Einar: "I don't want to see that anywhere." The only marks in the rail are the
    disclosure glyphs the folding headings already use, and an arrow on a next step."""
    from parascape import no_cross
    assert no_cross.find_crosses(rail) == []
    for row in rail.all_rows():
        assert row.text().strip() not in no_cross.CROSS_TEXT, row.text()


def test_every_row_has_a_symbol(rail):
    """Including the views and the next steps, which the flat rail had no rows for. A
    view is drawn with the symbol of the module it is behind."""
    rail.open_to("expression")
    for row in rail.all_rows():
        assert not row.icon().isNull(), row.key


def test_the_rail_is_wide_enough_for_whatever_is_unfolded(rail):
    """A rail measured on the module names alone elided "Condition comparison" the
    moment a module unfolded."""
    rail.open_to("expression")
    rail.setFixedWidth(rail.fitting_width())
    rail.resize(rail.fitting_width(), 900)
    assert rail.fitting_width() >= rail.WIDTH_MIN
