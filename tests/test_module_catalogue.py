"""One catalogue: the invariants that stop it becoming several again.

spaCR records a bug, #16i, in which Home grouped the applications by one table
and the tabs by another, so a section could acquire a tab it had no tiles for.
ParaScape had the same shape of problem in four pieces: `modules.MODULES`, the
icon map in `navigator.py`, the shell routes in `landing_install.open_module`
and the tool specs in `dock_module.SPECS`, which named the same ten modules in
different words. These tests are what keep them one.

The Qt tests at the foot need a real window, because the question "is every
module somebody can reach in the application in the catalogue?" cannot be
answered from the catalogue. That is the direction the first version of the
catalogue got wrong, and `tests/test_modules_in_window.py` records how: it was
written from memory and pointed six modules at views that did not exist.
"""
import pytest
from PySide6 import QtWidgets

from parascape import modules as MOD
from parascape import ports as P


# ------------------------------------------------------------------- sections

def test_no_section_is_empty():
    """spaCR has this one by name. A heading drawn over nothing is a promise
    the application does not keep, and it appears the day the last module in a
    section is retired, which is exactly when nobody is looking at the section."""
    for section in MOD.SECTIONS:
        assert MOD.section_members(section), f"{section} has no modules"


def test_no_module_is_in_two_sections():
    """`band` is one string, so a module cannot literally be in two. What CAN
    happen is a module appearing twice in the table under different sections,
    which is the same bug wearing a different hat."""
    seen = {}
    for section in MOD.SECTIONS:
        for module in MOD.section_members(section):
            assert module.key not in seen, (
                f"{module.key} is in {seen.get(module.key)} and {section}")
            seen[module.key] = section
    assert len(seen) == len(MOD.MODULES)


def test_every_module_is_in_a_declared_section():
    for module in MOD.MODULES:
        assert module.band in MOD.SECTIONS, f"{module.key}: {module.band!r}"
        assert module.section == module.band


def test_the_sections_are_computed_from_the_band_table():
    """Two lists of section names, drifting apart, is bug #16i itself."""
    assert MOD.SECTIONS == tuple(title for title, _note in MOD.BANDS)


def test_the_bands_and_the_sections_are_one_grouping_under_two_names():
    """The home screen says bands and the instructions say sections. One
    implementation, because two would be the thing this is here to prevent."""
    assert MOD.bands() == MOD.sections()
    for section, _note, members in MOD.sections():
        assert members == MOD.section_members(section)


def test_the_grouping_covers_every_module_exactly_once():
    listed = [m.key for _s, _n, members in MOD.sections() for m in members]
    assert sorted(listed) == sorted(m.key for m in MOD.MODULES)
    assert len(listed) == len(set(listed))


# ---------------------------------------------------------------------- ports

def test_every_module_in_the_catalogue_has_ports_declared():
    """A module with no declared subject is invisible to the navigator, to the
    tree and to every readiness answer, while looking perfectly fine on Home."""
    for module in MOD.MODULES:
        spec = module.ports
        assert spec.key == module.key
        assert spec.consumes or spec.produces, module.key


def test_every_declared_port_belongs_to_a_module_in_the_catalogue():
    """The other direction: ports for a key nothing can open."""
    keys = {m.key for m in MOD.MODULES}
    for key in P.known_modules():
        assert key in keys, f"ports declared for {key}, which is not a module"


def test_the_two_tables_are_keyed_the_same():
    assert set(P.known_modules()) == {m.key for m in MOD.MODULES}


def test_a_module_reports_its_subjects_without_a_second_table():
    """`consumes` and `produces` read `ports`; they are not copied here."""
    module = MOD.get("fill")
    assert P.CLAIM in module.consumes and P.SCORECARD in module.consumes
    assert module.produces == tuple(p.kind for p in module.ports.produces)


def test_the_modules_downstream_of_one_gene_are_all_in_the_catalogue():
    """Instruction 98 point 6: clicking a gene card offers these, so every one
    of them has to be a module somebody can actually open."""
    keys = {m.key for m in MOD.MODULES}
    assert set(P.consumers_of(P.GENE)) <= keys
    assert set(P.consumers_of(P.GENE_LIST)) <= keys


def test_nothing_hand_writes_what_follows_a_module_any_more():
    """Instruction 98, what WP8 had to close first: `DOWNSTREAM` was twenty-two
    hand-written lists beside a computed `ports.next_modules`, and the
    hand-written one was what the navigator offered. It is gone, and the
    catalogue is not allowed to grow another."""
    assert not hasattr(MOD, "DOWNSTREAM"), \
        "what follows a module is computed from the declared subjects, not listed"


def test_what_follows_a_module_is_what_reads_what_it_hands_on():
    """Every offer is a module that declares it reads one of the subjects the
    first one declares it produces -- which is the whole reason the graph can be
    checked rather than trusted."""
    for module in MOD.MODULES:
        produced = set(module.produces)
        for target in MOD.downstream(module.key):
            reads = set(MOD.get(target).consumes)
            assert reads & produced, \
                f"{module.key} offers {target}, which reads none of {sorted(produced)}"


def test_the_modules_downstream_of_one_gene_are_the_ones_that_take_a_gene():
    """Einar: "all the modules that take 1 gene and do something should come
    downstream of the gene card module"."""
    offered = MOD.downstream_of_subject(P.GENE, without="gene")
    assert set(offered) == set(P.consumers_of(P.GENE)) - {"gene"}
    assert set(offered) == {"cell", "expression", "information", "network",
                            "partners", "structure"}


def test_an_offer_nobody_can_take_is_never_the_first_one():
    """spaCR sorts its next steps ready first, blocked last. Structure has no
    view built, and before this it was the first thing offered under a gene."""
    offered = MOD.downstream_of_subject(P.GENE, without="gene")
    blocked = [k for k in offered if MOD.unavailable_reason(MOD.get(k))]
    assert offered[0] not in blocked
    assert offered[-len(blocked):] == tuple(blocked) if blocked else True


def test_the_symbol_of_a_module_is_declared_once():
    """Instruction 98: `iconset.MODULE_ICON` and the icon field of a catalogue
    row held the same fact twice. The catalogue feeds the icon set."""
    from parascape import iconset
    icons = iconset._load_module_icons()
    for module in MOD.MODULES:
        assert icons[module.key] == module.icon, module.key


# ----------------------------------------------------------- the declared row

def test_every_module_says_what_to_do_first():
    """The masthead instruction. It lived in `dock_module.SPECS` for ten
    modules and nowhere for the other twelve, so twelve mastheads said what the
    module IS and never what to DO with it."""
    for module in MOD.MODULES:
        assert module.instruction, f"{module.key} has no instruction"
        assert module.instruction.endswith("."), module.key
        assert module.instruction[0].isupper(), module.key


def test_every_module_has_an_icon_name():
    for module in MOD.MODULES:
        assert module.icon, f"{module.key} has no icon"


def test_a_builder_is_one_of_the_three_declared_kinds():
    for module in MOD.MODULES:
        assert module.builder in ("", "tool") or module.builder.startswith("shell:"), \
            f"{module.key}: {module.builder!r}"


def test_a_tool_module_names_a_dock_to_borrow():
    """`builder="tool"` means "wrap the dock named in opens", so a tool row
    with no destination would build a shell around nothing."""
    for module in MOD.tool_modules():
        assert module.opens, f"{module.key} is a tool with nothing to open"
        assert not module.opens.startswith("action:"), module.key


def test_the_tool_specs_are_computed_from_the_catalogue():
    """`dock_module.SPECS` was a second list of the same ten modules with its
    own wording for each. One table, so a tile and the masthead of the module
    behind it cannot describe it differently."""
    from parascape import dock_module as DM
    assert [s.key for s in DM.SPECS] == [m.key for m in MOD.tool_modules()]
    for spec in DM.SPECS:
        module = MOD.get(spec.key)
        assert spec.title == module.name
        assert spec.description == module.note
        assert spec.instruction == module.instruction
        dock, _, tab = module.opens.partition(":")
        assert (spec.dock, spec.tab) == (dock, tab)


def test_only_the_console_declines_to_lend_itself_a_console():
    from parascape import dock_module as DM
    assert [s.key for s in DM.SPECS if not s.console] == ["console"]


def test_the_shell_builders_are_read_from_the_catalogue():
    """Six keys, four builders: Partners is the Network module read by binding
    evidence, and Questions is the Infer module opened on its catalogue."""
    assert MOD.shell_module("map") == "map_module"
    assert MOD.shell_module("partners") == "network_module"
    assert MOD.shell_module("questions") == "infer_module"
    assert MOD.shell_module("drop") == "", "a tool is not built by a shell module"
    assert set(MOD.shell_modules()) == {"map_module", "gene_module",
                                        "network_module", "infer_module"}
    assert len(MOD.shell_modules()) == len(set(MOD.shell_modules()))


def test_a_module_is_built_by_a_shell_or_a_tool_or_opens_something():
    """A row with no builder and nothing to open is a tile that does nothing.
    Four of them are honest about it through `needs`; the rest must not be."""
    for module in MOD.MODULES:
        if module.builder or module.opens:
            continue
        assert module.needs, f"{module.key} builds nothing and opens nothing"


def test_shell_module_raises_on_a_key_that_is_not_a_module():
    with pytest.raises(KeyError):
        MOD.shell_module("not_a_module")


# ------------------------------------------------------------ in a real window

@pytest.fixture(scope="module")
def window():
    """The real window, built offscreen. The same fixture
    `tests/test_modules_in_window.py` uses, for the same reason: a catalogue
    can only be checked against the application it describes."""
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


def test_every_module_reachable_in_the_interface_is_in_the_catalogue(window):
    """The dock lists what somebody can open. Anything in it that is not a
    catalogue row has no section, no icon and no declared subject, and would
    therefore be missing from Home, from the navigator and from the tree while
    still being one click away in the dock."""
    keys = {m.key for m in MOD.MODULES} | {"home"}
    dock = window.landing.dock
    assert dock is not None, "the dock failed to build"
    for row in dock._rows:
        assert row.key in keys, f"{row.key} is reachable and not in the catalogue"


def test_every_catalogue_module_is_reachable_in_the_interface(window):
    """And the other way: a row nobody can get to."""
    reachable = {row.key for row in window.landing.dock._rows}
    for module in MOD.MODULES:
        assert module.key in reachable, f"{module.key} is in the catalogue and unreachable"


def test_every_shell_builder_actually_built_something(window):
    """`shell_modules` names the Python modules the window builds its shells
    from, and the window keeps each under the same name. A builder named in the
    catalogue and absent from the window is a tile that opens the old route."""
    for name in MOD.shell_modules():
        assert hasattr(window, name), f"{name} was never built"


def test_every_tool_module_found_its_dock(window):
    docks = {d.windowTitle() for d in window.findChildren(QtWidgets.QDockWidget)}
    for module in MOD.tool_modules():
        dock = module.opens.partition(":")[0]
        assert dock in docks, f"{module.key} borrows {dock!r}, which does not exist"
        assert module.key in window.tool_modules, module.key
