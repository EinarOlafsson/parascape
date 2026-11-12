"""The declared subjects: the graph, the readiness answers, and the purity rule.

These exist because everything else in the interface restructure is computed
from `parascape/ports.py`. The navigator offers the modules downstream of a
gene, the tree walks the same graph, and a greyed tile carries the reason this
module produced. A wrong declaration here is not a wrong declaration in one
place; it is a wrong answer everywhere, and none of the places it goes wrong in
would obviously look wrong.

The purity test at the foot is the one that is easy to lose: `ports` is
imported by the module catalogue, and the module catalogue is imported while
the application is starting. The day somebody adds `import pandas` here for one
convenience, every launch pays for it and no test fails.
"""
import subprocess
import sys

import pytest

from parascape import ports as P


# ------------------------------------------------------------- the vocabulary

def test_every_subject_has_a_name_and_a_meaning():
    """A subject with no words for it reaches a user as a bare key."""
    for kind in P.ALL_SUBJECTS:
        assert kind in P.SUBJECT_NOTES, kind
        assert kind in P.SUBJECT_MEANINGS, kind
        assert P.SUBJECT_NOTES[kind], kind
        assert P.SUBJECT_MEANINGS[kind], kind


def test_the_short_name_of_a_subject_fits_in_a_sentence():
    """No commas and no trailing stop: these are dropped into "no X in hand"."""
    for kind, name in P.SUBJECT_NOTES.items():
        assert "," not in name, f"{kind}: {name!r} would break the sentence"
        assert not name.endswith("."), f"{kind}: {name!r}"
        assert name == name.strip()


def test_an_unknown_subject_describes_itself_rather_than_raising():
    assert P.describe_subject("no_such_subject") == "no_such_subject"
    assert P.explain_subject("no_such_subject") == "no_such_subject"


def test_every_declared_port_uses_a_subject_in_the_vocabulary():
    for key in P.known_modules():
        spec = P.module_ports(key)
        for port in spec.consumes + spec.produces:
            assert port.kind in P.ALL_SUBJECTS, f"{key}.{port.role}: {port.kind}"


def test_the_seven_subjects_einar_named_are_all_declared():
    """Instruction 98 names them: one gene, a gene list, a gene pair, a space,
    a screen, a trait, a compound. The six beyond them are additions and may be
    argued with; these seven are the requirement."""
    for kind in (P.GENE, P.GENE_LIST, P.GENE_PAIR, P.SPACE, P.SCREEN,
                 P.TRAIT, P.COMPOUND):
        assert kind in P.ALL_SUBJECTS


def test_nothing_takes_a_compound_yet_and_that_is_recorded_not_hidden():
    """A compound is in the vocabulary because Einar named it, and no module
    reads one because drug sensitivity is in ParaScape only as per-gene columns.
    If a compound module ever arrives this test is what says to update the
    docstring beside the constant rather than leave it lying."""
    assert P.consumers_of(P.COMPOUND) == ()
    assert P.producers_of(P.COMPOUND) == ()


# ------------------------------------------------------------- declarations

def test_every_module_declares_something():
    """A module with neither an input nor an output is not in the graph at all."""
    for key in P.known_modules():
        spec = P.module_ports(key)
        assert spec.consumes or spec.produces, key
        assert spec.summary, f"{key} has no summary"


def test_every_port_says_what_it_is_for():
    for key in P.known_modules():
        for port in P.module_ports(key).consumes + P.module_ports(key).produces:
            assert port.description, f"{key}.{port.role}"
            assert port.role, key


def test_roles_are_unique_inside_one_module():
    """Two ports of one kind in one module are told apart by role and nothing
    else, so a repeated role silently merges them."""
    for key in P.known_modules():
        spec = P.module_ports(key)
        roles = [p.role for p in spec.consumes + spec.produces]
        assert len(roles) == len(set(roles)), f"{key}: {roles}"


def test_module_ports_raises_on_an_unknown_key():
    with pytest.raises(P.UnknownModule):
        P.module_ports("no_such_module")


def test_registering_a_second_declaration_for_one_key_is_refused():
    """Two declarations claiming one key, last one quietly winning, is how a
    plugin would blank a module without anybody noticing."""
    with pytest.raises(ValueError):
        P.register_module_ports(P.ModulePorts(key="gene"))


def test_registering_a_declaration_with_a_repeated_role_is_refused():
    with pytest.raises(ValueError):
        P.register_module_ports(P.ModulePorts(
            key="made_up_for_the_test",
            consumes=(P.Port(P.GENE, "same"), P.Port(P.TRAIT, "same"))))


def test_a_new_module_can_join_the_graph_without_editing_the_table():
    """The seam. A module written later registers rather than being added here."""
    key = "test_only_module"
    try:
        P.register_module_ports(P.ModulePorts(
            key=key, summary="only in this test",
            consumes=(P.Port(P.SPACE, "map", description="a map to read"),),
            produces=(P.Port(P.CLAIM, "out", description="what it claims"),)))
        assert key in P.known_modules()
        assert key in P.consumers_of(P.SPACE)
        assert key in P.next_modules("map")
    finally:
        P.PORTS.pop(key, None)


def test_overwrite_replaces_a_declaration_when_it_is_asked_for():
    original = P.module_ports("structure")
    try:
        P.register_module_ports(
            P.ModulePorts(key="structure", summary="replaced",
                          consumes=(P.Port(P.GENE, "gene", description="x"),)),
            overwrite=True)
        assert P.module_ports("structure").summary == "replaced"
    finally:
        P.register_module_ports(original, overwrite=True)


def test_port_lookup_by_role():
    spec = P.module_ports("infer")
    assert spec.port("question").kind == P.QUESTION
    with pytest.raises(KeyError):
        spec.port("no_such_role")


def test_required_and_optional_kinds_are_reported_separately():
    spec = P.module_ports("infer")
    assert spec.required_kinds() == (P.QUESTION,)
    assert set(spec.optional_kinds()) == {P.GENE_LIST, P.TRAIT}


# ------------------------------------------------------------------ the graph

def test_the_modules_downstream_of_one_gene():
    """Instruction 98, point 6: "click it once -> the modules downstream of one
    gene". This is that list, and it is the whole reason the file exists."""
    downstream = P.consumers_of(P.GENE)
    assert set(downstream) == {"cell", "expression", "gene", "information",
                               "network", "partners", "structure"}


def test_the_modules_downstream_of_a_list():
    """Instruction 98, point 7."""
    downstream = P.consumers_of(P.GENE_LIST)
    for key in ("infer", "list", "map", "network", "export"):
        assert key in downstream, key


def test_required_only_narrows_the_consumers():
    """Map takes a list and works without one; List cannot."""
    assert "list" in P.consumers_of(P.GENE_LIST, required_only=True)
    assert "map" in P.consumers_of(P.GENE_LIST)
    assert "map" not in P.consumers_of(P.GENE_LIST, required_only=True)


def test_the_modules_that_run_work_are_what_produce_a_job():
    """Jobs and Console were a dead end: they consume a job, nothing declared
    producing one, and the only advice `check_ready` could give was that this
    "cannot be reached from here" -- while every long build in the application
    was handing work to the same runner."""
    assert set(P.producers_of(P.JOB)) == {"map", "infer", "predict", "strategies"}
    assert "jobs" in P.next_modules("map")
    fix = P.check_ready("jobs").errors[0].fix
    assert "build a map" in fix


def test_clicking_a_gene_is_a_way_of_producing_one():
    """The first version declared no producer of a gene at all, so a reader with
    no gene chosen was told one could not be reached -- with the map in front of
    them handing out genes on every click."""
    assert "map" in P.producers_of(P.GENE)
    assert "network" in P.producers_of(P.GENE)


def test_what_follows_the_map():
    following = P.next_modules("map")
    assert "gallery" in following, "a map goes to the gallery"
    assert "gene" in following, "a gene clicked on the map goes to Gene"
    assert "map" not in following, "a module never follows itself"


def test_export_is_a_sink():
    """What leaves ParaScape is a file, and a file is not a subject. The tree
    has to stop somewhere and this is where."""
    assert P.module_ports("export").produces == ()
    assert P.next_modules("export") == ()


def test_drop_needs_nothing_and_starts_a_chain():
    """A file on disk is where subjects come from, not a subject itself."""
    assert P.module_ports("drop").consumes == ()
    assert "list" in P.next_modules("drop")


def test_upstream_of_fill_names_who_makes_what_it_needs():
    upstream = P.upstream_modules("fill")
    assert "infer" in upstream, "claims come from Infer"
    assert "scorecard" in upstream, "the track record comes from the Scorecard"


def test_the_graph_functions_refuse_an_unknown_module():
    for call in (P.next_modules, P.upstream_modules, P.describe_ports):
        with pytest.raises(P.UnknownModule):
            call("no_such_module")


# ------------------------------------------------------------------ the route

def test_a_route_ends_at_the_module_asked_for():
    steps = P.route(P.GENE, "gallery")
    assert steps and steps[-1] == "gallery"
    assert "map" in steps, "a map is what the gallery holds"


def test_a_route_covers_everything_the_module_needs_not_just_the_first_thing():
    """The failure this is written against: a search over one subject at a time
    reached a claim, declared Fill reachable, and left out the track record Fill
    also requires. A half route reads as an instruction that will work."""
    steps = P.route(P.GENE, "fill")
    assert steps and steps[-1] == "fill"
    held = {P.GENE}
    for key in steps[:-1]:
        spec = P.module_ports(key)
        assert set(spec.required_kinds()) <= held, f"{key} cannot run at that point"
        held |= {p.kind for p in spec.produces}
    assert set(P.module_ports("fill").required_kinds()) <= held


def test_a_module_that_needs_nothing_is_already_reachable():
    assert P.route(P.GENE, "export") == ("export",)
    assert P.route("", "datasets") == ("datasets",)


def test_a_route_that_does_not_exist_is_empty_rather_than_invented():
    """Nothing produces a compound, so nothing can be routed to a module that
    would need one. Returning a plausible-looking chain would be worse."""
    P.register_module_ports(P.ModulePorts(
        key="needs_a_compound", summary="only in this test",
        consumes=(P.Port(P.COMPOUND, "compound", description="a compound"),)))
    try:
        assert P.route(P.GENE, "needs_a_compound") == ()
    finally:
        P.PORTS.pop("needs_a_compound", None)


def test_a_route_is_the_same_on_every_call():
    """An offer that reorders itself between two openings of one page reads as
    a bug to whoever is looking at it."""
    assert P.route(P.GENE, "fill") == P.route(P.GENE, "fill")


# --------------------------------------------------------------- is it ready?

def test_a_module_with_what_it_needs_is_ready():
    ready = P.check_ready("gene", P.Subject(P.GENE, "TGME49_215220"))
    assert ready
    assert ready.ok and "gene" in ready.satisfied
    assert "can work on this" in ready.reason


def test_a_bare_kind_is_accepted_as_well_as_a_value():
    """The navigator asks "would Gene take a gene?" before anybody has named
    one; a module about to run asks about the gene in hand. Both are this call."""
    assert P.check_ready("gene", P.GENE)
    assert P.check_ready("gene", P.Subject(P.GENE, "TGME49_215220"))


def test_a_missing_subject_blocks_and_says_where_to_get_one():
    answer = P.check_ready("structure")
    assert not answer
    assert answer.missing == (P.GENE,)
    problem = answer.errors[0]
    assert problem.level == P.ERROR and problem.is_error
    assert "one gene" in problem.message
    assert "navigator" in problem.fix, "the way in without a module is named first"


def test_an_empty_subject_is_a_different_answer_from_no_subject():
    """"A question is what I would hand you, but I have not chosen one" has a
    different fix from "I have no question": choose one, versus go and get one."""
    empty = P.check_ready("infer", P.Subject(P.QUESTION, ""))
    none = P.check_ready("infer")
    assert not empty and not none
    assert "none is chosen" in empty.errors[0].message
    assert "Choose" in empty.errors[0].fix
    assert "none is in hand" in none.errors[0].message


def test_a_zero_is_a_chosen_value_and_not_emptiness():
    """A trait whose value is 0 was chosen. Treating it as empty would make a
    legitimate choice look like no choice at all."""
    assert bool(P.Subject(P.TRAIT, 0))
    assert not bool(P.Subject(P.TRAIT, ""))
    assert not bool(P.Subject(P.TRAIT, []))
    assert not bool(P.Subject(P.TRAIT, None))


def test_an_optional_subject_warns_and_does_not_block():
    answer = P.check_ready("map")
    assert answer.ok, "Map works with no genes chosen; it shows every gene"
    assert answer.warnings
    assert all(not p.is_error for p in answer.warnings)


def test_a_subject_the_module_does_not_read_warns_rather_than_blocks():
    """Handing Map a question does not stop Map. Saying nothing would leave
    somebody believing their question had been taken into account."""
    answer = P.check_ready("map", P.Subject(P.QUESTION, "host_facing"))
    assert answer.ok
    unused = [p for p in answer.warnings if p.subject == P.QUESTION]
    assert unused and "does not read" in unused[0].message
    assert "infer" in unused[0].fix


def test_several_subjects_can_be_offered_at_once():
    answer = P.check_ready("fill", P.Subject(P.CLAIM, ["TGME49_215220"]),
                           have=[P.Subject(P.SCORECARD, {"right": 0.9})])
    assert answer.ok, P.format_readiness(answer)
    assert set(answer.satisfied) == {"claims", "record"}


def test_have_accepts_a_single_subject_as_well_as_a_list():
    answer = P.check_ready("fill", P.Subject(P.CLAIM, ["x"]),
                           have=P.Subject(P.SCORECARD, {"right": 0.9}))
    assert answer.ok


def test_fill_will_not_accept_a_claim_without_its_track_record():
    """The one check in this file that is about the science rather than the
    plumbing. A claim is a guess; the track record is what says how often a
    guess of that kind is right. Instruction 98: no calibration may be dropped
    in the move, and writing an unmeasured claim into the map beside measured
    ones drops the distinction between them."""
    answer = P.check_ready("fill", P.Subject(P.CLAIM, ["TGME49_215220"]))
    assert not answer
    assert P.SCORECARD in answer.missing
    assert "track record" in answer.errors[0].message


def test_every_problem_carries_a_fix():
    """A problem somebody cannot act on was reported to the wrong person."""
    for key in P.known_modules():
        for problem in P.check_ready(key).problems:
            assert problem.fix, f"{key}: {problem.message}"
            assert problem.message


def test_check_ready_refuses_an_unknown_module():
    with pytest.raises(P.UnknownModule):
        P.check_ready("no_such_module", P.GENE)


def test_the_reason_is_one_line_and_names_the_first_problem():
    answer = P.check_ready("fill")
    assert "\n" not in answer.reason
    assert answer.reason.startswith("fill cannot:")
    assert "more" in answer.reason, "the count of the rest is in the line"


def test_the_full_report_pairs_each_problem_with_its_fix():
    text = P.format_readiness(P.check_ready("fill"))
    assert text.startswith("NOT READY: fill")
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith(("error", "warning")):
            assert lines[i + 1].strip().startswith("fix:"), line
    assert str(P.check_ready("fill")) == text


def test_every_module_can_be_checked_without_a_window_or_any_data():
    """A readiness answer has to cost less than the thing it protects."""
    for key in P.known_modules():
        answer = P.check_ready(key)
        assert answer.module == key
        assert isinstance(answer.ok, bool)


# ------------------------------------------------------------------ describing

def test_describe_ports_names_both_halves():
    text = P.describe_ports("infer")
    assert "takes:" in text and "hands on:" in text
    assert "[optional]" in text


def test_a_module_that_takes_nothing_says_so():
    text = P.describe_ports("drop")
    assert "(nothing)" in text


def test_the_subject_listing_reads_as_a_sentence():
    assert P._listing([]) == ""
    assert P._listing(["a"]) == "a"
    assert P._listing(["a", "b"]) == "a or b"
    assert P._listing(["a", "b", "c"]) == "a, b or c"


# ----------------------------------------------------------------- the purity

def test_importing_ports_pulls_in_no_qt_and_no_scientific_stack():
    """Run in a fresh interpreter, because by the time this file runs under the
    full suite pandas and Qt are long since imported by something else and the
    check would pass whatever `ports` did.

    It matters because the module catalogue imports this, and the catalogue is
    read while the application is starting. One convenience import here is paid
    for on every launch, and nothing else would fail.
    """
    code = ("import sys, parascape.ports as P; "
            "assert P.GENE; "
            "bad = [m for m in ('PySide6', 'numpy', 'pandas', 'sklearn', 'umap', "
            "'matplotlib', 'pyqtgraph') if m in sys.modules]; "
            "print(','.join(bad))")
    done = subprocess.run([sys.executable, "-c", code],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "", f"ports dragged in {done.stdout.strip()}"


def test_ports_imports_nothing_else_out_of_parascape():
    """`parascape.modules` imports THIS module, not the other way round. An
    import back the other way is a cycle, and the cycle would only show up at
    application start where it is hardest to read."""
    code = ("import sys, parascape.ports; "
            "print(','.join(sorted(m for m in sys.modules "
            "if m.startswith('parascape.') and m != 'parascape.ports')))")
    done = subprocess.run([sys.executable, "-c", code],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "", f"ports imported {done.stdout.strip()}"
