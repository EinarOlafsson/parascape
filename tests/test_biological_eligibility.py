"""The biological admission step of instruction 65.04.

Every test here is a trap the first version of the script fell into. They are kept because each one
produced a confident WRONG verdict rather than an error, which is the failure mode that matters for
an audit: a wrong number looks exactly like a right one.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import verify_biological_eligibility as E  # noqa: E402

from parascape import datasets, slots as S  # noqa: E402


# ------------------------------------------------------------------ the assay vocabulary


def test_one_assay_spelled_three_ways_folds_to_one():
    """`Ribo-seq`, `RiboSeq` and `riboseq` are one assay; so are `IP-MS` and `IPMS`.

    13 registry entries carry these spellings. Anything that groups by assay without folding them
    reports one assay as two or three, and the audit then compares groups that were never apart.
    """
    assert E.normalise_kind("Ribo-seq") == E.normalise_kind("RiboSeq") == E.normalise_kind("riboseq")
    assert E.normalise_kind("IP-MS") == E.normalise_kind("IPMS")


def test_the_registry_still_spells_those_assays_more_than_one_way():
    """A live reading of the registry, so that tidying it up makes this test tell you, not fail you."""
    spellings = {}
    for e in datasets.REGISTRY:
        spellings.setdefault(E.normalise_kind(e.kind), set()).add(e.kind)
    drifted = {k: v for k, v in spellings.items() if len(v) > 1}
    assert set(drifted) <= {"riboseq", "ipms"}, f"a new assay spelling has drifted: {drifted}"


def test_every_assay_in_the_registry_is_in_the_capability_table():
    """An assay the table does not know becomes `unknown`, which queues it for reading rather than
    failing it -- but silently growing that queue is how an audit stops being an audit."""
    missing = sorted({e.kind for e in datasets.REGISTRY
                      if E.normalise_kind(e.kind) not in E.MEASURES})
    assert not missing, f"assays with no declared capability: {missing}"


# ------------------------------------------------------------------ the capability table


@pytest.mark.parametrize("kind,axis", [
    ("BioID", "localization"),            # a bait in an outer membrane reports what is near it
    ("proteomics", "chemistry"),          # thermal proteome profiling is target engagement
    ("in_vitro_evolution", "sequence"),   # evolved clones are sequenced whole-genome
    ("RNAseq", "regulation"),             # decay after actinomycin D is regulation
    ("RNAseq", "sequence"),               # a long-read reannotation gives UTR architecture
    ("nanopore", "regulation"),           # which isoform is used is regulation
    ("literature", "chemistry"),          # the only route to a Toxoplasma resistance allele
])
def test_an_assay_that_really_does_measure_this_axis_is_not_called_a_mismatch(kind, axis):
    """The first table was too NARROW and reported 28 correct pairings as biologically ineligible.

    A too-strict admission table rejects good evidence with a confident verdict, which is the same
    failure instruction 65 exists to prevent, pointing the other way.
    """
    assert axis in E.MEASURES[E.normalise_kind(kind)], \
        f"{kind} does measure {axis}; calling it a mismatch rejects a correct pairing"


def test_a_ribosome_profiling_rna_column_answers_a_transcription_slot():
    """The channel is chosen by the COLUMN, not by the assay label.

    A ribosome-profiling deposit ships the footprint channel and the matched RNA channel. Reading
    the assay alone called five correct candidates a mismatch -- among them `riboseq_mrna_`, which
    says what it is in its own name.
    """
    for pattern in ("rna99395_intracellular_", "rna129869_", "riboseq_mrna_",
                    "rna245775_parent_tachy_"):
        assert E._channel_axis(pattern) == "transcription", pattern
    for pattern in ("rpf99395_intracellular_", "rpf245775_parent_tachy_", "te302107_tachy_"):
        assert E._channel_axis(pattern) == "", pattern


# ------------------------------------------------------------------ the life stage


def test_an_oocyst_column_is_not_read_as_a_tissue_cyst():
    """"oocyst" contains "cyst". A bradyzoite stem that matched it would claim every oocyst column
    as a tissue cyst, which is a different stage in a different host."""
    assert E._matches("bradyzoite", "expr_cyst")
    assert not E._matches("bradyzoite", "expr_oocyst")
    assert E._matches("oocyst", "expr_oocyst")


def test_a_stage_named_only_in_prose_is_not_a_declaration(monkeypatch):
    """`gse245775` mentions "bradyzoite" in its note in order to DENY it: the arms are labelled
    `pre-bradyzoites`, "named for that and not for bradyzoite". A word search over prose finds the
    word inside the sentence written to exclude it, so prose counts as something to read, never as
    a fit."""
    entry = next(e for e in datasets.REGISTRY if e.key == "gse245775")
    assert "brady" in str(entry.note).lower(), "the fixture no longer exercises the trap"
    assert E.declares_stage(entry, "bradyzoite", (), "rpf245775_eif12ko_") == "prose"
    assert E.declares_stage(entry, "tachyzoite", (), "rna245775_parent_tachy_") == "column"


def test_a_candidate_does_not_inherit_a_stage_from_a_sibling_column():
    """`rna108740_Day` was credited with "bradyzoite" because `rna108740_Tissue_cysts` sits in the
    same deposit. A neighbouring column is not evidence about this one -- which is the whole reason
    slots are filled by pattern rather than by deposit."""
    entry = next(e for e in datasets.REGISTRY if e.key == "gse108740")
    siblings = tuple(entry.columns or ())
    assert E.declares_stage(entry, "bradyzoite", siblings, "rna108740_Day") != "column"
    assert E.declares_stage(entry, "bradyzoite", siblings, "rna108740_Tissue_cysts") == "column"


# ------------------------------------------------------------------ the slots that are exempt


def test_a_slot_that_is_never_a_feature_is_not_asked_whether_its_quantity_fits():
    """Asking whether a source fits an axis literally called "NEVER a feature" is a category error,
    not a finding."""
    assert "NEVER a feature" in E.EXEMPT_AXES
    assert {"never", "metadata"} <= E.EXEMPT_ROLES


def test_a_host_slot_admits_a_host_organism_and_refuses_the_parasite():
    host = next(s for s in S.all_slots() if s.unit == "host_gene")
    human = next(e for e in datasets.REGISTRY if e.organism == "Hs")
    parasite = next(e for e in datasets.REGISTRY if e.organism == "Tg")
    assert E.check(host, human, ())[0]["organism"] == "fits"
    assert E.check(host, parasite, ())[0]["organism"] == "mismatch"


# ------------------------------------------------------------------ the verdict itself


def test_an_unknown_queues_the_row_and_never_passes_it():
    """Instruction 65: unknown counts stay unknown rather than becoming zero or a guess. The same
    rule holds for biology, where a guess is worse."""
    assert E.admission({"organism": "fits", "stage": "unknown"}) == "needs reading"
    assert E.admission({"organism": "mismatch", "stage": "fits"}) == "INELIGIBLE"
    assert E.admission({"organism": "fits", "stage": "fits",
                        "evidence_grade": "inferred"}) == "admitted"


def test_being_inferred_rather_than_measured_does_not_by_itself_block_a_source():
    """Curated evidence is the only route to a Toxoplasma resistance allele. It must be marked as
    inferred and it must still be admissible; the grade is a label, not a verdict."""
    curated = next(e for e in datasets.REGISTRY if e.key == "curated_resistance_alleles")
    slot = next(s for s in S.all_slots() if s.key == "Tg_resistance_conferring_mutation")
    verdicts, _ = E.check(slot, curated, tuple(curated.columns or ()))
    assert verdicts["evidence_grade"] == "inferred"
    assert E.admission(verdicts) != "INELIGIBLE"
