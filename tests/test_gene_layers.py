"""The three layers of the gene card: observed, predicted (with how often right, against chance), still open."""
from __future__ import annotations

import os
import sys
from types import SimpleNamespace

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parascape import gene_layers as GL  # noqa: E402


def _findings():
    return pd.DataFrame([
        {"pmid": "1", "gene_id": "G1", "kind": "localisation", "finding": "G1 is in the rhoptries.",
         "method": "IFA", "method_tier": 2, "status": "established", "year": 2010,
         "named_in_abstract": True, "open_question": "Its partner is unknown."},
        {"pmid": "2", "gene_id": "G1", "kind": "phenotype", "finding": "G1 is not needed for invasion.",
         "method": "knockout", "method_tier": 1, "status": "refuted", "year": 2015,
         "named_in_abstract": True, "open_question": ""},
        {"pmid": "3", "gene_id": "G1", "kind": "activity", "finding": "A family member is a kinase.",
         "method": "assay", "method_tier": 1, "status": "established", "year": 2016,
         "named_in_abstract": False, "open_question": ""},
        {"pmid": "4", "gene_id": "G2", "kind": "localisation", "finding": "G2 is in the IMC.",
         "method": "IFA", "method_tier": 2, "status": "established", "year": 2020,
         "named_in_abstract": True, "open_question": ""},
    ])


def test_observed_lines_are_the_genes_own_findings_with_their_paper():
    lines = GL.observed_lines("G1", _findings())
    text = " ".join(lines)
    assert "rhoptries" in text and "pubmed.ncbi.nlm.nih.gov/1/" in text
    assert "G2" not in text
    assert lines[0].startswith("<b>Where it is</b>"), "place comes first"


def test_a_finding_whose_protein_the_abstract_does_not_name_is_left_out():
    assert "kinase" not in " ".join(GL.observed_lines("G1", _findings()))


def test_a_refutation_reads_as_one():
    text = " ".join(GL.observed_lines("G1", _findings()))
    assert "<b>Ruled out:</b> G1 is not needed for invasion." in text


def test_the_summary_names_the_strongest_kind_of_evidence():
    s = GL.observed_summary("G1", _findings())
    assert s.startswith("2 findings shown by experiment in 2 papers")
    assert "genetic loss of function" in s


def test_a_gene_with_no_findings_has_no_observed_layer():
    assert GL.observed_lines("G9", _findings()) == []
    assert GL.observed_summary("G9", _findings()) == ""


def test_chance_is_the_share_of_the_commonest_true_answer():
    ledger = pd.DataFrame({"target": ["place"] * 4 + ["other"], "gene_id": ["a", "b", "c", "d", "a"],
                           "truth": ["ER", "ER", "ER", "IMC", "x"]})
    share, answer = GL.chance_rate(ledger, "place")
    assert answer == "ER" and abs(share - 0.75) < 1e-9
    assert GL.chance_rate(ledger, "absent") == (None, "")


def test_each_agreement_line_is_set_against_chance():
    reading = SimpleNamespace(predicted=None, unknown=[])
    lines = GL.predicted_lines(reading, [("<b>place</b>: 5 of 6 methods said ER.", "place")],
                               {"place": (0.75, "ER")})
    assert "would be right 75% of the time" in lines[0]


def test_still_open_holds_what_is_unmeasured_and_what_papers_leave_open():
    reading = SimpleNamespace(unknown=[("phenotype", "what losing it does")])
    lines = GL.open_lines(reading, "G1", _findings())
    assert lines[0] == "Not measured: what losing it does."
    assert any("Its partner is unknown." in s and "pubmed.ncbi.nlm.nih.gov/1/" in s for s in lines)


def test_the_shipped_findings_load_and_read_for_a_known_effector():
    """ROP18 (TGME49_205250) was missing from its papers before the second reading pass."""
    lines = GL.observed_lines("TGME49_205250")
    assert lines, "ROP18 has no observed findings in the shipped table"
