"""THE ANALYSES BEHIND THE PAPER, declared once and runnable from the program.

Instruction 101, items B2 and B3: every analysis in the paper must be runnable inside ParaScape by
somebody who did not write it, and must leave its result beside the record of what it was made
from. These tests hold the catalogue honest -- a finding naming a script that does not exist, or a
claim that describes a method instead of the parasite, is caught here rather than by a reader.
"""
import os
import re

import pytest

from parascape import modules as MOD
from parascape import paper_findings as PF


def test_every_finding_names_a_script_that_exists():
    missing = [f.key for f in PF.FINDINGS if not os.path.exists(PF.script_path(f))]
    assert not missing, f"findings naming a script that is not there: {missing}"


def test_every_figure_of_the_paper_has_at_least_one_finding():
    """Twelve figures; 1a is a schematic and carries no analysis of its own."""
    figures = {f.figure.rstrip("abc") for f in PF.FINDINGS}
    assert figures >= {str(n) for n in range(2, 13)} - {"6"} | {"6"}, figures


def test_keys_and_scripts_are_unique():
    keys = [f.key for f in PF.FINDINGS]
    assert len(keys) == len(set(keys))
    scripts = [f.script for f in PF.FINDINGS]
    assert len(scripts) == len(set(scripts))


@pytest.mark.parametrize("finding", PF.FINDINGS, ids=lambda f: f.key)
def test_the_claim_is_about_the_parasite_and_carries_a_number(finding):
    """Einar's test for every figure: it ends in a sentence about the parasite.

    A claim with no number is an impression, and a claim that names its own method is a methods
    note. Both have been written here and both are caught by this.
    """
    assert finding.claim.endswith("."), finding.key
    assert re.search(r"\d", finding.claim), f"{finding.key}: the claim carries no number"
    method_words = ("cross-validated", "we ran", "the script", "this analysis computes")
    said = finding.claim.lower()
    assert not any(w in said for w in method_words), f"{finding.key}: that is a methods note"


def test_a_finding_that_has_never_run_reports_nothing_rather_than_guessing():
    made_up = PF.Finding("nope", "0", "x", "Nothing, 0 times.", "build_corpus_inventory.py",
                         "there_is_no_such_result_*")
    assert PF.latest(made_up) == ""
    assert PF.record(made_up) == {}


def test_a_record_says_where_it_came_from():
    """The B3 half: a result must carry the directory it was written to."""
    ran = [f for f in PF.FINDINGS if PF.latest(f)]
    assert ran, "no finding has produced a result in this checkout"
    for finding in ran:
        record = PF.record(finding)
        assert "result_directory" in record, finding.key
        assert os.path.isdir(record["result_directory"]), finding.key


def test_the_module_is_in_the_catalogue_and_sits_with_the_corpus():
    module = MOD.get("findings")
    assert module.band == "Data", "the analyses are about the corpus, not about a gene"
    assert MOD.shell_module("findings") == "findings_module"


def test_typing_what_it_is_called_finds_it():
    kind, found, view = MOD.search_all("Findings")[0]
    assert (kind, found.key, view) == ("module", "findings", None)
