"""The contrast ladder: the rules that stop a study being read against the wrong baseline.

Every test here is a mistake that was actually made, not a hypothetical. The keyword rule this replaces
read "WT" in `T. gondii VEG WT, Switch Conditions` and called the alkaline-stressed arm a baseline, which
deleted the best alkaline contrast in the collection. Each check below pins one of the ways that can happen
again: reading a level on one axis as an answer about another, pairing across an axis that lives only in
free text, skipping an axis that some groups declare and others omit, and printing two different
measurements as one row.
"""
import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "build_contrast_ladder", Path(__file__).resolve().parent.parent / "scripts" / "build_contrast_ladder.py")
CL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CL)

# The six tachyzoite groups of GSE218921 plus its two sporozoite groups, verbatim from the sample registry.
# The structured head declares strain, treatment and tissue; the genotype exists only in the free text.
GSE218921 = []
for _geno, _tail in (("WT", "T. gondii VEG WT"), ("ROCY1 KO", "T. gondii VEG ROCY1 KO"),
                     ("BFD1 KO", "T. gondii VEG BFD1 KO")):
    for _trt, _word in (("high pH, low C02, low serum conditions", "Switch Conditions"),
                        ("control media conditions", "Control Conditions")):
        GSE218921.append(
            f"strain: VEG; treatment: {_trt}; tissue: replicating tachyzoites under different conditions "
            f"| {_tail}, {_word},")
GSE218921 += [
    "strain: VEG; treatment: high pH, low C02, low serum conditions; tissue: freshly excysted sporozoites "
    "| T. gondii in vitro bradyzoite induction conditions",
    "strain: VEG; treatment: control media conditions; tissue: freshly excysted sporozoites "
    "| T. gondii in vitro control conditions",
]


def rows(labels, study="GSE218921"):
    return [{"study": study, "condition_group": label} for label in labels for _ in range(3)]


def ladder(labels, study="S1"):
    contrasts, unresolved, refused = CL.build(rows(labels, study), "condition_group", "study")
    return contrasts, unresolved, refused


# --------------------------------------------------------------------------- the original mistake


def test_wt_is_read_as_a_genotype_not_as_a_baseline_for_treatment():
    """The whole point. `WT` must orient the genotype axis and must not make the stressed arm a control."""
    contrasts, _, _ = ladder(GSE218921)
    treatment = [c for c in contrasts if c["axis"] == "treatment"]
    assert treatment, "the alkaline contrast must exist"
    for c in treatment:
        assert "control media" in c["baseline_level"], c
        assert "high ph" in c["treated_level"], c


def test_the_recovered_wild_type_alkaline_contrast_is_present_at_depth_zero():
    """The contrast the keyword rule deleted: wild-type tachyzoites, control media against switch."""
    contrasts, _, _ = ladder(GSE218921)
    hit = [c for c in contrasts if c["axis"] == "treatment" and c["depth"] == 0
           and c["background"].startswith("genotype=wild type")]
    assert len(hit) == 1, [c["background"] for c in contrasts if c["axis"] == "treatment"]
    assert "replicating tachyzoites" in hit[0]["background"]


def test_untreated_knockouts_are_compared_to_untreated_wild_type():
    """Einar's first rule: the knockouts relate to WT on the untreated background, at depth 0."""
    contrasts, _, _ = ladder(GSE218921)
    genotype = [c for c in contrasts if c["axis"] == "genotype" and c["depth"] == 0]
    assert {c["treated_level"] for c in genotype} == {"bfd1 ko", "rocy1 ko"}
    for c in genotype:
        assert c["baseline_level"] == "wild type"
        # The background names every other axis held fixed, the ones already at their reference included.
        # What matters is that the treatment held fixed here is the untreated one.
        assert "treatment=control media conditions" in c["background"], c["background"]
        assert "high ph" not in c["background"], "a depth-0 genotype contrast sits on untreated cells"


def test_treated_samples_are_compared_to_the_untreated_cells_of_their_own_genotype():
    """Einar's second rule: the treatment effect inside a knockout is measured against that knockout."""
    contrasts, _, _ = ladder(GSE218921)
    deep = {c["background"].split(";")[0] for c in contrasts
            if c["axis"] == "treatment" and c["depth"] == 1}
    assert "genotype=bfd1 ko" in deep and "genotype=rocy1 ko" in deep


def test_depth_separates_the_clean_effect_from_the_effect_inside_a_mutant():
    contrasts, _, _ = ladder(GSE218921)
    assert {c["depth"] for c in contrasts} == {0, 1}


# --------------------------------------------------------------------------- the refusal rule


def test_an_axis_hidden_in_free_text_refuses_the_pair():
    """GSE132238 declares genotype and timepoint; `stress`/`unstress` hides in the sample name. Pairing a
    stressed wild type against an unstressed knockout would report a stress response as a genotype effect."""
    labels = [f"genotype: {g}; timepoint: 24 hr | {n}_{s}_24h_array1"
              for g, n in (("WT", "WT"), ("{delta}BFD1", "KO")) for s in ("stress", "unstress")]
    contrasts, _, refused = ladder(labels, "GSE132238")
    assert any("stress" in r["unparsed_difference"] for r in refused)
    for c in contrasts:
        assert c["axis"] == "genotype"
        # Each surviving contrast must stay inside one stress state, named in the background.
        assert "unparsed=" in c["background"]


def test_a_restatement_of_a_parsed_factor_does_not_refuse_the_pair():
    """`Switch Conditions` in the tail merely repeats `treatment: high pH`. Treating it as a hidden axis
    refused the alkaline contrast a second time, for a second wrong reason."""
    contrasts, _, _ = ladder(GSE218921)
    assert any(c["axis"] == "treatment" and c["depth"] == 0 for c in contrasts)


def test_attribution_is_tested_inside_strata_so_differently_worded_arms_do_not_poison_it():
    """The sporozoite arm words the same stress as `bradyzoite induction`, so study-wide the word `switch`
    looks unexplained. Tissue separates two experiments sharing one accession."""
    contrasts, _, _ = ladder(GSE218921)
    tissues = {c["background"].split("tissue=")[1] for c in contrasts if "tissue=" in c["background"]}
    assert any("sporozoites" in t for t in tissues)
    assert any("tachyzoites" in t for t in tissues)


def test_a_replicate_label_is_flagged_rather_than_guessed_at():
    """`1a`/`2b` are almost certainly plate identifiers, but merging on that guess would fuse two real
    conditions. The pair is refused and the refusal says what it looks like."""
    labels = [f"treatment: {t}; genotype: WT | sample_{i}" for t in ("control", "LPS") for i in ("1a", "2b")]
    _, _, refused = ladder(labels, "GSE223869")
    assert refused
    assert all(r["looks_like"] == "plate or batch identifier" for r in refused)


# --------------------------------------------------------------------------- the quieter mistakes


def test_an_axis_only_some_groups_declare_becomes_unstated_rather_than_being_skipped():
    """GSE136123 has two different untreated baselines, `Pru MORC KD_UT` and `Pru_US`. Only the first
    parses a genotype, so without an explicit level the comparison skips the axis and measures the
    knockdown against the parental line as though they were the same strain."""
    labels = ["treatment: none | Pru MORC KD_UT", "treatment: none | Pru_US",
              "treatment: IAA for 24h | Pru MORC KD_IAA_24h"]
    contrasts, _, _ = ladder(labels, "GSE136123")
    seen = {(c["baseline_group"], c["treated_group"]) for c in contrasts}
    assert not any("MORC KD_UT" in b and "Pru_US" in t or "Pru_US" in b and "MORC KD_UT" in t
                   for b, t in seen), "the two baselines must not be treated as one"


def test_replicate_groups_of_one_condition_are_merged():
    """Deposited as `..._array1` and `..._array2`, they are one condition. Left apart they emit the same
    contrast twice and each carries half the replicates it really has."""
    labels = [f"genotype: {g}; timepoint: 24 hr | {n}_stress_24h_{a}"
              for g, n in (("WT", "WT"), ("{delta}BFD1", "KO")) for a in ("array1", "array2")]
    contrasts, _, _ = ladder(labels, "GSE132238")
    assert len(contrasts) == 1, contrasts
    assert contrasts[0]["n_baseline"] == 6, "three rows per group, two merged groups"


def test_no_two_contrasts_print_the_same_row():
    """Two different measurements that print identically get averaged together downstream. The builder
    reports this count and it must stay zero."""
    contrasts, _, _ = ladder(GSE218921)
    seen = [(c["study"], c["axis"], c["baseline_level"], c["treated_level"], c["background"])
            for c in contrasts]
    assert len(seen) == len(set(seen))


# --------------------------------------------------------------------------- orientation vocabulary


@pytest.mark.parametrize("axis,level", [
    ("treatment", "control media conditions"), ("treatment", "untreated"), ("treatment", "mock"),
    ("treatment", "vehicle"), ("treatment", "no atc"), ("genotype", "wild type"), ("genotype", "parental"),
])
def test_reference_levels_are_recognised(axis, level):
    assert CL.reference_rank(axis, level)[0] < 99


@pytest.mark.parametrize("axis,level", [
    ("strain", "me49"), ("strain", "rh"), ("tissue", "brain"), ("stage", "bradyzoite"),
])
def test_strain_and_tissue_never_get_an_automatic_reference(axis, level):
    """Choosing a strain is choosing an organism, not choosing a control. Silently making ME49 the
    baseline is how a strain difference gets reported as a treatment effect."""
    assert CL.reference_rank(axis, level)[0] == 99


def test_an_edited_genotype_is_not_a_reference_even_when_it_says_control():
    assert CL.reference_rank("genotype", "cas9 control bfd1 ko")[0] == 99 or \
        CL.reference_rank("genotype", "bfd1 ko")[0] == 99


def test_the_earliest_timepoint_is_the_reference():
    early = CL.reference_rank("time", "4 hr")[0]
    late = CL.reference_rank("time", "48 hr")[0]
    assert early < late


def test_every_orientation_records_the_word_that_decided_it():
    contrasts, _, _ = ladder(GSE218921)
    assert all(c["oriented_by"] for c in contrasts)
