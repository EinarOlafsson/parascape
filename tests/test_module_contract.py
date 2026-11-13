"""THE CONTRACT EVERY MODULE KEEPS, stated once, as a table.

Instruction 98, WP11. Einar asked for one thing above all the rest: *"every module should be
thematically consistent. every analasys should be thematically consistent and intuitive to use
and get to."* Ten work packages built that. This file is what stops it decaying the next time a
module is added, because a contract nobody checks is a contract that lasts until the next
commit.

WHAT THIS FILE IS, AND WHAT IT DELIBERATELY IS NOT
==================================================
:data:`CONTRACT` is the whole contract, one clause per line. Each clause either carries a
``check`` that runs HERE against every catalogue row, or names the test that already holds it.
Nothing is asserted twice: two tests of one rule are two things to keep in step, they drift,
and the one that drifts is the one nobody reads. A clause with neither a check nor a holder
fails :func:`test_every_clause_is_either_checked_here_or_held_elsewhere`, and a holder that has
been renamed away fails :func:`test_every_test_this_contract_leans_on_still_exists` -- so the
contract cannot quietly lose a clause either.

The table-driven test is parametrised by module AND by clause, so a failure reads

    test_the_contract[jobs-it offers a control only for a subject it reads]

which names the module and the missing part without anybody reading the assertion.

WHAT IS CHECKED HERE, AND WHY IT WAS NOT CHECKED BEFORE
=======================================================
Four of the clauses had no test anywhere, and each of the four is a failure that has actually
happened on this instruction rather than one imagined for a test:

* **a second table of a fact the catalogue holds.** spaCR records this as bug #16i. ParaScape
  has had it three times: `iconset.MODULE_ICON` against the icon field, `modules.DOWNSTREAM`
  against `ports.next_modules`, and -- found by WP11 -- the masthead of every shell module
  against its own catalogue row, where five of seven had already drifted ("Gene" on the tile,
  "Gene card" on the masthead). Two are closed by name; this file stops a fourth.
* **a control for a subject the module does not read.** Every tool module was given a "Your
  gene" field, so Jobs asked for a gene identifier above the list of what is running, and
  Console, Gallery, Drop, Strategies and Scorecard did the same.
* **greyed for want of something it does not use.** Network was greyed whenever the star map
  dock failed to build -- a dock the Network module does not use.
* **a recipe surface wherever there are settings.** Map and Network had one and the other six
  shells did not, although each has settings of its own. The distinction is the one Einar drew: a gene or
  a question is a SUBJECT and not a setting, so a module whose Settings column holds only the
  subject needs no recipe. That is why the shell records which group is which
  (`ModuleShell.subject_groups`) instead of this file guessing from the headings.
"""
import ast
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import pytest
from PySide6 import QtWidgets

from parascape import modules as MOD
from parascape import ports as P

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent / "parascape"

#: The controls a recipe can keep. A group holding one of these is settings; a group holding
#: none of them is prose, a drawing or a pair of buttons, and there is nothing to save.
KEEPABLE = (QtWidgets.QComboBox, QtWidgets.QCheckBox, QtWidgets.QSpinBox,
            QtWidgets.QDoubleSpinBox, QtWidgets.QSlider, QtWidgets.QLineEdit)

#: The subjects the application can carry INTO a module (`landing_install.carry_subject`).
#: A module that requires one of these and has nowhere to put it would take an offer from the
#: gene card or the gene list card and open empty, which is the one thing an offer exists to
#: stop.
CARRIED = (P.GENE, P.GENE_LIST)


# --------------------------------------------------------------------------- the clauses

@dataclass(frozen=True)
class Clause:
    """One line of the contract.

    :param name: what the clause asks for, as it appears in a failing test id.
    :param check: ``check(module, window) -> str``; the empty string passes and anything else
        is what is missing, reported beside the module key. ``None`` when the clause is held
        by a test somewhere else.
    :param held_by: ``("file::test", ...)`` naming the tests that already hold this clause.
        Read by :func:`test_every_test_this_contract_leans_on_still_exists`, so a rename over
        there is a failure here rather than a clause that silently stops being checked.
    :param why: one line on why the clause is in the contract at all.
    """

    name: str
    why: str
    check: Optional[Callable] = None
    held_by: tuple = field(default_factory=tuple)


def _shell_of(window, module):
    """The shell this module opens in, or None when it has no shell in this window.

    Several keys may share one shell -- Partners is the Network module read by binding
    evidence, Questions is the Infer module opened on its catalogue -- so the shell carries
    the key of the module it BELONGS to. A clause about the shell is checked against the
    owner, which is why most of the checks below start by comparing the two keys.
    """
    built_by = MOD.shell_module(module.key)
    if built_by:
        return getattr(window, built_by, None)
    return getattr(window, "tool_modules", {}).get(module.key)


def _owns_its_shell(window, module):
    """The shell, when this module is the one the shell was built for; otherwise None."""
    shell = _shell_of(window, module)
    if shell is None or getattr(shell, "key", None) != module.key:
        return None
    return shell


# ---- reachable by typing its own name

def _findable_by_name(module, _window) -> str:
    """Typing the name of a module finds that module, ahead of everything else.

    `test_navigator` checks this for one module. It is a property of every row: a module
    somebody cannot reach by typing what it is called is a module they will not find, whatever
    else is true of it.
    """
    hits = MOD.search_all(module.name)
    if not hits:
        return f"typing {module.name!r} in the navigator finds nothing"
    kind, found, view = hits[0]
    if kind != "module" or found.key != module.key:
        what = f"the {view.name!r} view of {found.key}" if view is not None else found.key
        return f"typing {module.name!r} offers {what} first, not {module.key!r}"
    return ""


# ---- a recipe surface wherever there are settings

def _settings_that_a_recipe_could_keep(shell) -> list:
    """The groups of the Settings column that hold settings rather than the subject."""
    out = []
    for name, section in shell.settings_groups.items():
        if name in shell.subject_groups or name == "Recipes":
            continue
        if any(section.findChildren(cls) for cls in KEEPABLE):
            out.append(name)
    return out


def _has_recipes_where_it_has_settings(module, window) -> str:
    shell = _owns_its_shell(window, module)
    if shell is None:
        return ""
    groups = _settings_that_a_recipe_could_keep(shell)
    if groups and shell.recipes is None:
        return (f"it has settings ({', '.join(groups)}) and no Recipes group, so there is no "
                f"way to keep them under a name as every other module allows")
    if not groups and shell.recipes is not None:
        return ("it has a Recipes group and nothing but the subject to put in it, which is a "
                "promise the module cannot keep")
    if shell.recipes is not None and not shell.recipes.settings:
        return "its Recipes group keeps no setting at all"
    return ""


# ---- a control only for a subject it reads

def _offers_only_subjects_it_reads(module, window) -> str:
    shell = _owns_its_shell(window, module)
    if shell is None:
        return ""
    reads = set(module.consumes)
    wrong = {name: kind for name, kind in shell.subject_groups.items() if kind not in reads}
    if wrong:
        said = "; ".join(f"{name!r} takes a {kind}" for name, kind in sorted(wrong.items()))
        return (f"{said}, which its ports do not declare (it reads {sorted(reads)}) -- a tile "
                f"that asks for something it does not use")
    return ""


def _has_somewhere_to_put_what_it_requires(module, window) -> str:
    shell = _owns_its_shell(window, module)
    if shell is None:
        return ""
    required = {p.kind for p in module.ports.consumes if p.required and p.kind in CARRIED}
    if not required:
        return ""
    # `one_of` says two ports meet one need, so filling either is enough.
    offered = set(shell.subject_groups.values())
    if not (required & offered):
        return (f"it requires {sorted(required)} and its Settings column has no group for one, "
                f"so an offer carrying it would open the module empty")
    return ""


# ---- greyed only for a reason of its own

def _greyed_only_for_its_own_reason(module, window) -> str:
    """A module may be greyed for what its own row declares, and for nothing else.

    The two honest reasons are the `needs` written on the row, and the destination the row
    names in `opens` failing to build in this window. Anything else is the Network failure:
    it was greyed whenever the star map dock did not build, a dock it does not open.
    """
    from parascape.landing_install import _reasons
    reason = _reasons(window).get(module.key, "")
    if not reason or module.needs:
        return ""
    destination = module.opens.partition(":")[0]
    if destination and destination in reason:
        return ""                    # its own view did not build, which is its own reason
    return (f"the window greys it -- {reason!r} -- and that is neither what its row says it "
            f"needs nor the {destination or 'nothing'!r} it says it opens, which is how "
            f"Network came to be greyed for want of a dock it does not use")


# ---- the masthead is the catalogue row

def _masthead_is_the_catalogue_row(module, window) -> str:
    shell = _owns_its_shell(window, module)
    if shell is None:
        return ""
    if shell.title_label.text() != module.name:
        return (f"the masthead says {shell.title_label.text()!r} and the tile says "
                f"{module.name!r}")
    if shell.description_label.text() != module.note:
        return (f"the masthead describes it as {shell.description_label.text()!r} and the tile "
                f"as {module.note!r}")
    if shell.instruction_label.text() != module.instruction:
        return (f"the masthead says to {shell.instruction_label.text()!r} and the tile says "
                f"to {module.instruction!r}")
    return ""


#: THE CONTRACT. Every module keeps every line of it.
CONTRACT = (
    Clause("declares its ports",
           "a module with no declared subject is invisible to the navigator, the tree and "
           "every readiness answer, while looking perfectly fine on Home",
           held_by=("test_module_catalogue::test_every_module_in_the_catalogue_has_ports_declared",
                    "test_module_catalogue::test_the_two_tables_are_keyed_the_same")),
    Clause("has an icon that exists and renders",
           "a symbol that cannot be loaded is a blank square on a tile, and a symbol built "
           "some other way drifts from spaCR unnoticed",
           held_by=("test_module_catalogue::test_every_module_has_an_icon_name",
                    "test_module_catalogue::test_the_symbol_of_a_module_is_declared_once",
                    "test_symbols_and_no_cross::test_every_module_has_a_shipped_spacr_symbol",
                    "test_symbols_and_no_cross"
                    "::test_every_symbol_renders_visible_pixels_in_the_ink_of_the_theme")),
    Clause("appears in exactly one section",
           "a heading drawn over nothing is a promise the application does not keep, and a "
           "module in two sections is two places to look for one thing",
           held_by=("test_module_catalogue::test_no_section_is_empty",
                    "test_module_catalogue::test_no_module_is_in_two_sections",
                    "test_module_catalogue::test_the_grouping_covers_every_module_exactly_once")),
    Clause("is reachable in the interface",
           "a catalogue row nobody can open, and a module somebody can open that has no row, "
           "are the same bug seen from the two ends",
           held_by=("test_module_catalogue::test_every_catalogue_module_is_reachable_in_the_interface",
                    "test_module_catalogue::test_every_module_reachable_in_the_interface_is_in_the_catalogue",
                    "test_modules_in_window::test_every_destination_exists_in_the_window")),
    Clause("is reachable from the tree",
           "the rail is the answer to 'what can ParaScape do'; a module missing from it is a "
           "module that exists only for somebody who already knew its name",
           held_by=("test_tree::test_every_module_is_filed_exactly_once",
                    "test_tree::test_the_rail_lists_home_and_every_module_once")),
    Clause("is reachable by typing its own name",
           "the navigator is the front door Einar asked for, and a module its search field "
           "cannot name is behind that door rather than through it",
           check=_findable_by_name),
    Clause("its views exist and open",
           "WP10 gave every panel behind the navigator a name; a declared view that is not "
           "there, and a tab nobody declared, are both how that decays",
           held_by=("test_views::test_every_declared_view_exists_in_the_window",
                    "test_views::test_every_tab_in_the_window_is_a_view_a_module_or_a_signpost",
                    "test_views::test_opening_a_view_reports_what_it_opened_and_raises_nothing",
                    "test_tree::test_every_declared_view_is_in_the_tree")),
    Clause("has a recipe surface wherever it has settings",
           "spaCR lets every screen keep its settings under a name; a module here that cannot "
           "is a module whose working set-up has to be rebuilt by hand each time",
           check=_has_recipes_where_it_has_settings),
    Clause("offers a control only for a subject it reads",
           "Jobs asked for 'a gene id or a word of its product' above the list of what is "
           "running, because every tool module was given a gene field whether it read one or "
           "not",
           check=_offers_only_subjects_it_reads),
    Clause("has somewhere to put what it requires",
           "a module that needs a gene and has no field for one takes an offer carrying a gene "
           "and opens empty, which is the thing the offer exists to stop",
           check=_has_somewhere_to_put_what_it_requires),
    Clause("says the same thing on its tile and on its masthead",
           "the tile said 'Gene' and the masthead said 'Gene card'; five of the seven shells "
           "had drifted from their own row this way, which is bug #16i a third time",
           check=_masthead_is_the_catalogue_row),
    Clause("is greyed only for a reason its own row declares",
           "Network was greyed whenever the star map dock failed to build -- a dock the "
           "Network module does not use",
           check=_greyed_only_for_its_own_reason),
    Clause("carries a docstring",
           "the generated reference IS the documentation of this project, so a public name "
           "with no docstring is a hole in it",
           held_by=("test_docstrings::test_every_module_explains_itself",
                    "test_docstrings::test_every_public_function_and_class_has_a_docstring",
                    "test_docstrings::test_public_methods_of_public_classes_have_docstrings")),
    Clause("says what to do first, in plain words",
           "twelve mastheads said what the module IS and never what to DO with it",
           held_by=("test_module_catalogue::test_every_module_says_what_to_do_first",
                    "test_home::test_every_module_says_what_you_get_not_what_it_is")),
    Clause("has a tooltip on every control",
           "85 controls had none when that audit was first run, almost all of them where "
           "choosing wrongly costs the most",
           held_by=("test_app_menus::test_no_control_anywhere_is_left_without_an_explanation",)),
    Clause("shows no small cross",
           "Einar asked for this one by name -- 'I don't want to see that anywhere' -- and it "
           "creeps back with every new search field and every new dock",
           held_by=("test_symbols_and_no_cross::test_the_real_window_has_no_small_cross_anywhere",
                    "test_symbols_and_no_cross::test_no_cross_survives_opening_every_module")),
)

BY_NAME = {clause.name: clause for clause in CONTRACT}


# --------------------------------------------------------- the contract is itself complete

def test_every_clause_is_either_checked_here_or_held_elsewhere():
    """A clause with neither is a line of the contract nobody enforces, which is worse than
    not writing it down: it reads as covered."""
    idle = [c.name for c in CONTRACT if c.check is None and not c.held_by]
    assert idle == [], f"clauses nothing checks: {idle}"
    both = [c.name for c in CONTRACT if c.check is not None and c.held_by]
    assert both == [], f"clauses asserted twice, which is what this file exists to avoid: {both}"


def test_every_clause_says_why_it_is_in_the_contract():
    """Each line is here because something went wrong. A clause with no reason attracts the
    next person to delete it when it gets in the way."""
    for clause in CONTRACT:
        assert len(clause.why.split()) >= 8, clause.name


def test_every_test_this_contract_leans_on_still_exists():
    """The other half of not duplicating: if a clause says `test_views` holds it and that test
    is renamed or deleted, the clause stops being checked and nothing says so. Reading the
    names out of the files is what makes leaning on them safe."""
    missing = []
    for clause in CONTRACT:
        for reference in clause.held_by:
            module_name, _, test_name = reference.partition("::")
            path = HERE / f"{module_name}.py"
            if not path.is_file():
                missing.append(f"{reference}: no {path.name}")
                continue
            names = {node.name for node in ast.walk(ast.parse(path.read_text(encoding="utf8")))
                     if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
            if test_name not in names:
                missing.append(f"{reference}: no such test")
    assert missing == [], f"the contract leans on tests that are gone: {missing}"


# ------------------------------------------------------ no second table keyed by module

def _module_keyed_tables() -> list:
    """Every module-level collection outside the catalogue whose keys are module keys.

    THE GUARD AGAINST A FOURTH #16i. `iconset.MODULE_ICON` and `modules.DOWNSTREAM` were each
    a table of module keys sitting beside the catalogue and answering a question the catalogue
    could answer, and each drifted. Both are closed by name in `test_module_catalogue`; this
    finds the next one by shape rather than by name, because the next one will have a
    different name.

    Four keys is the bar: two or three module keys in a tuple is a list of special cases with
    a comment beside it, and those are honest. A collection naming four or more of the
    twenty-three is a table ABOUT the modules.
    """
    keys = {m.key for m in MOD.MODULES}
    found = []
    for path in sorted(ROOT.glob("*.py")):
        if path.name in ("modules.py", "ports.py"):
            continue                        # the two that are allowed to be keyed by module
        tree = ast.parse(path.read_text(encoding="utf8"))
        for node in tree.body:
            if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue
            value = node.value
            if isinstance(value, ast.Dict):
                literal = [k.value for k in value.keys
                           if isinstance(k, ast.Constant) and isinstance(k.value, str)]
            elif isinstance(value, (ast.Tuple, ast.List, ast.Set)):
                literal = [e.value for e in value.elts
                           if isinstance(e, ast.Constant) and isinstance(e.value, str)]
            else:
                continue
            if len(literal) >= 4 and set(literal) <= keys:
                targets = getattr(node, "targets", [getattr(node, "target", None)])
                name = next((t.id for t in targets if isinstance(t, ast.Name)), "?")
                found.append(f"{path.name}:{node.lineno} {name}")
    return found


def test_nothing_keeps_a_second_table_about_the_modules():
    assert _module_keyed_tables() == [], (
        "a collection keyed by module key, outside the catalogue and the ports: it answers a "
        "question modules.py can answer, and it will drift from it, which is the bug spaCR "
        "records as #16i. Compute it from the catalogue instead")


# ---------------------------------------------------------------------- the real window

@pytest.fixture(scope="module")
def window():
    """The real window, built offscreen, with every module opened once.

    Opening them all first matters: a shell builds its panels when it is first shown, so a
    clause read off an unopened module would be read off half a module.
    """
    from parascape import app as A
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                   # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    win.show()
    qapp.processEvents()
    from parascape.landing_install import open_module
    for module in MOD.MODULES:
        if MOD.unavailable_reason(module) or module.opens.startswith("action:"):
            continue
        try:
            open_module(win, module.key)
        except Exception:                        # a module that will not open fails its own
            continue                             # clause below, not the fixture
        qapp.processEvents()
    yield win
    from parascape.map_module import wait_for_projections
    wait_for_projections()
    win.close()
    qapp.processEvents()


@pytest.mark.parametrize("name", [c.name for c in CONTRACT if c.check is not None])
@pytest.mark.parametrize("key", [m.key for m in MOD.MODULES])
def test_the_contract(window, key, name):
    """Every clause of the contract, against every row of the catalogue.

    The failure names the module and the part it is missing, which is the whole point: a
    module added without an icon, without a section, or with a field for a subject it does not
    read fails here, by name, rather than being noticed later by Einar.
    """
    missing = BY_NAME[name].check(MOD.get(key), window)
    assert not missing, f"{key}: {name} -- {missing}"


def test_a_control_built_when_a_subject_arrives_also_carries_its_explanation(window):
    """The audit in `test_app_menus` walks the window as it is built. The cards build more
    controls when a subject arrives -- the matches of a gene, the row of next steps -- and
    those are the ones somebody is looking at when they most need the explanation."""
    from parascape.gene_module import open_in as gene_in
    from parascape.list_module import open_in as list_in
    gene = str(window.nodes.gene_id.iloc[0])
    genes = ", ".join(str(g) for g in window.nodes.gene_id.iloc[:12])
    gene_in(window, gene)
    QtWidgets.QApplication.processEvents()
    list_in(window, genes)
    QtWidgets.QApplication.processEvents()
    missing = []
    for shell in (window.gene_module, window.list_module):
        standard = {b for box in shell.findChildren(QtWidgets.QDialogButtonBox)
                    for b in box.buttons()}
        for cls in (QtWidgets.QComboBox, QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox,
                    QtWidgets.QCheckBox, QtWidgets.QPushButton):
            for widget in shell.findChildren(cls):
                if widget in standard or widget.toolTip().strip():
                    continue
                missing.append(f"{shell.key}: {cls.__name__} {widget.text()!r}")
    assert missing == [], f"controls with no explanation once a subject is in hand: {missing}"


def test_a_module_with_nothing_to_set_gives_the_room_to_what_it_shows(window):
    """Found by rendering the Jobs module (instruction 98, WP11). Six tool modules stopped
    offering a gene field they had no business offering, and the Settings column was then a
    third of the screen holding an empty box, with the panel squeezed into the rest."""
    from parascape.landing_install import open_module
    bare = []
    for module in MOD.MODULES:
        shell = _owns_its_shell(window, module)
        if shell is None or shell.settings_groups:
            continue
        open_module(window, module.key)
        QtWidgets.QApplication.processEvents()
        if shell.settings_scroll.isVisible():
            bare.append(module.key)
    assert bare == [], f"an empty Settings column is still taking the room: {bare}"


def test_a_squeezed_action_is_shortened_with_an_ellipsis_rather_than_cut(qapp):
    """Found by rendering the Network module (instruction 98, WP11).

    Eight actions in one row left "Open a saved bundle…" drawn as "pen a saved bundle" and
    "Save everything as a bundle…" as "everything as a bur" -- cut at both ends, so neither
    said which button it was and both read as a different, wrong action. A plain
    `QPushButton` does that; `ActionButton` shortens with an ellipsis instead.

    Compared against a plain button rather than against the elided string, because the string
    is what the code under test computes: the question is what is DRAWN.
    """
    from parascape.shell import ActionButton
    label = "Save everything as a bundle…"
    ours, plain = ActionButton(label), QtWidgets.QPushButton(label)
    for button in (ours, plain):
        button.resize(90, 28)
        button.show()
    qapp.processEvents()
    assert ours.grab().toImage() != plain.grab().toImage(), \
        "the action button draws the same as a plain one, so it is not shortening anything"
    wide = ActionButton(label)
    wide.resize(ours.fontMetrics().horizontalAdvance(label) + 80, 28)
    wide.show()
    qapp.processEvents()
    plain.resize(wide.size())
    qapp.processEvents()
    assert wide.grab().toImage() == plain.grab().toImage(), \
        "a label that fits must be drawn exactly as any other button draws it"
    for button in (ours, plain, wide):
        button.close()


def test_every_action_keeps_its_whole_label_where_it_can_still_be_read(window):
    """The other half: `text()` is still the whole label, because that is what the module
    keyed the button by, and the tooltip opens with it, because a shortened label on screen
    makes the tooltip the only place the action is named in full."""
    from parascape.shell import ActionButton
    for module in MOD.MODULES:
        shell = _owns_its_shell(window, module)
        if shell is None:
            continue
        for label, button in shell.actions.items():
            assert isinstance(button, ActionButton), f"{module.key}: {label!r} is a plain button"
            assert button.text() == label, \
                f"{module.key}: the button reports {button.text()!r} rather than {label!r}"
            assert label in button.toolTip(), \
                f"{module.key}: the tooltip of {label!r} does not carry the whole label"


def test_every_table_says_what_its_columns_are_before_anything_is_loaded(window):
    """Found by rendering the modules and looking (instruction 98, WP11).

    A table sized its columns only once something had been loaded into it, so the empty table
    a person meets FIRST drew each heading centred in a box too narrow for it and cut it at
    both ends: the Drop module said "NTIFIER IN THE F" for "identifier in the file", and the
    Scorecard said "NFERRED CLAIMS". Neither says which column it is, and neither says that
    anything has been cut, which is worse than an ellipsis by exactly the amount a reader
    trusts it.

    Checked on the EMPTY tables, because a table with rows in it has been resized to its
    content and the fault is hidden.
    """
    from parascape.landing_install import open_module
    from parascape.panes.tables import HEADING_MARGIN
    narrow = []
    for module in MOD.MODULES:
        shell = _owns_its_shell(window, module)
        if shell is None:
            continue
        open_module(window, module.key)
        for _ in range(4):
            QtWidgets.QApplication.processEvents()
        for table in shell.findChildren(QtWidgets.QTableWidget):
            header = table.horizontalHeader()
            if not header.isVisible() or table.rowCount():
                continue
            for column in range(table.columnCount()):
                item = table.horizontalHeaderItem(column)
                if item is None:
                    continue
                wanted = header.fontMetrics().horizontalAdvance(item.text()) + HEADING_MARGIN
                if header.sectionSize(column) < wanted:
                    narrow.append(f"{module.key}: {item.text()!r} has "
                                  f"{header.sectionSize(column)} px and needs {wanted}")
    assert narrow == [], f"headings cut rather than shown: {narrow}"


def test_the_subject_groups_of_every_module_name_a_declared_subject(window):
    """Not that the module reads it -- `test_the_contract` holds that -- but that the kind is
    a subject at all. A typo in a tag would otherwise read as "this group is for nothing"."""
    for module in MOD.MODULES:
        shell = _owns_its_shell(window, module)
        if shell is None:
            continue
        for name, kind in shell.subject_groups.items():
            assert kind in P.ALL_SUBJECTS, f"{module.key} ▸ {name}: {kind!r} is not a subject"


def test_no_module_promises_a_subject_in_words_that_its_ports_leave_out(window):
    """What a module SAYS it takes, against what it declares.

    The Screen tile said it joined your screen and opened a dock of RNA-seq published by somebody else;
    its ports said nothing about a screen until WP10. The first instruction is where a module
    promises a subject -- "Name a gene", "Choose a question" -- so that line is read against
    the declared ports. Only the plain forms are looked for, because a note is prose and
    "genes it was not shown" is not a promise to take a list.
    """
    promises = (
        (r"\ba gene\b|\bone gene\b|\bthe gene\b", P.GENE),
        (r"\ba list of genes\b|\byour genes\b|\bgene list\b", P.GENE_LIST),
        (r"\ba question\b|\bthe question\b", P.QUESTION),
        (r"\byour (?:measured )?table\b|\byour screen\b", P.SCREEN),
        (r"\ba trait\b", P.TRAIT),
    )
    wrong = []
    for module in MOD.MODULES:
        said = module.instruction.lower()
        declared = set(module.consumes)
        for pattern, kind in promises:
            if re.search(pattern, said) and kind not in declared:
                wrong.append(f"{module.key}: says {module.instruction!r} and does not read a "
                             f"{kind}")
    assert wrong == [], wrong
