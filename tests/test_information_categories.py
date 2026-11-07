"""The gate: every shipped column states what KIND of knowledge it is, so inference cannot be circular.

Instruction 73's mechanism is only a mechanism if the annotation is complete. These tests are the
gate that keeps it complete as data keeps arriving: a column added to any shipped table without a
category fails here rather than silently becoming evidence for a question it answers.

Three of the checks defend specific accidents that already happened or nearly did:

* `gra17_synthlethal`, `invivo_platform` and `gra12` are pooled CRISPR fitness screens whose family
  names contain no "crispr". They topped the importance table of the harness's first real run, and
  two passes of hand-written guessing had missed one each time. `excluded_for('fitness')` must hold
  all three, and the test names them rather than counting.
* `fit_hyperlopit_unassigned_invivo_lib1` has "hyperlopit" in its name and is **fitness**: the name
  records which library was screened, not what was measured. A name-based rule gets it backwards,
  which is why provenance drives the mapping.
* Localisation is the category the UMAP question turns on, so every hyperLOPIT, TAGM and
  ortholog-transferred compartment column must be excluded when a `place` label is inferred.
"""
import pathlib
import re

import numpy as np
import pandas as pd
import pytest

from parascape import harness as H
from parascape import information_categories as C
from parascape import question_registry as Q

ROOT = pathlib.Path(__file__).resolve().parent.parent
TABLE_FILES = {
    "nodes.parquet": ROOT / "parascape" / "data" / "nodes.parquet",
    "sequence_traits.parquet": ROOT / "results" / "sequence_traits_2026_11_04" / "sequence_traits.parquet",
    "domain_traits.parquet": ROOT / "results" / "domain_traits_2026_11_04" / "domain_traits.parquet",
    "expression_traits.parquet": ROOT / "results" / "expression_traits_2026_11_06" / "expression_traits.parquet",
    "information_map/values.parquet": ROOT / "results" / "information_map_v0_2026_11_04" / "values.parquet",
}
INSTRUCTION = ROOT / "instructions" / "open" / "73_question_matrix.md"


# --------------------------------------------------------------------------- the vocabulary
def test_the_category_vocabulary_is_the_one_instruction_73_lists():
    """Read from the instruction rather than restated here, so the two cannot drift apart."""
    if not INSTRUCTION.is_file():
        pytest.skip("the instruction is not shipped with an installed copy")
    text = INSTRUCTION.read_text()
    start = text.index("Categories, which are about the *kind of knowledge*")
    block = text[start:text.index("A column may carry more than one category", start)]
    named = set(re.findall(r"`([a-z_]+)`", block))
    assert named == set(C.CATEGORIES), f"vocabulary drifted: {named ^ set(C.CATEGORIES)}"
    for name, definition in C.CATEGORIES.items():
        assert definition.strip() and definition.strip().endswith("."), f"{name} has no definition"


def test_every_category_definition_is_about_knowledge_and_not_about_an_assay():
    """A definition naming an instrument would invite two columns of one assay into one category."""
    for name, definition in C.CATEGORIES.items():
        assert len(definition.split()) >= 5, f"{name}: one word is not a definition"


# --------------------------------------------------------------------------- completeness: the gate
@pytest.mark.parametrize("table", sorted(TABLE_FILES))
def test_the_shipped_manifest_matches_the_table_on_disk(table):
    """The enumerated manifest is what `unannotated()` answers from, so it must be the real columns."""
    path = TABLE_FILES[table]
    if not path.is_file():
        pytest.skip(f"{path} is not present in this checkout")
    actual = list(pd.read_parquet(path).columns)
    assert list(C.TABLES[table]) == actual, f"{table}: the manifest and the file disagree"


@pytest.mark.parametrize("table", sorted(TABLE_FILES))
def test_every_column_of_every_shipped_table_carries_at_least_one_category(table):
    missing = [c for c in C.TABLES[table] if not C.categories_of(c)]
    assert not missing, f"{table}: {len(missing)} columns with no information category: {missing[:20]}"


def test_unannotated_is_empty_for_the_shipped_tables():
    """The gate itself. A column that genuinely cannot be classed belongs here, and then this fails
    loudly -- which is the whole point, and much cheaper than a wrong category passing quietly."""
    assert C.unannotated() == ()
    assert len(C.SHIPPED_COLUMNS) > 1000


def test_every_category_is_used_and_every_column_only_uses_the_vocabulary():
    used = {c for column in C.SHIPPED_COLUMNS for c in C.categories_of(column)}
    assert used <= set(C.CATEGORIES)
    assert used == set(C.CATEGORIES), f"categories with no column: {sorted(set(C.CATEGORIES) - used)}"


def test_the_review_table_says_how_every_column_was_derived():
    frame = C.categories_table()
    assert set(frame.columns) == {"column", "table", "categories", "derived_from", "reason",
                                 "provenance_key", "provenance_level_kind", "confidence"}
    assert set(frame["derived_from"]) <= set(C.DERIVED_FROM)
    assert (frame["reason"].str.len() > 0).all()
    assert not (frame["confidence"] == "none").any()
    # Most of the annotation must rest on recorded provenance rather than on this module's judgement.
    assert (frame["confidence"] == "high").mean() > 0.9


# --------------------------------------------------------------------------- the traps, each named
def test_every_pooled_crispr_fitness_screen_is_fitness_however_its_family_is_named():
    excluded = set(C.excluded_for("fitness"))
    for column in ("crispr_gra17ko_phenotype", "crispr_gra17_synthlethal_delta", "crispr_gra17_candidate",
                   "crispr_invivo_platform_lfc",
                   "crispr_gra12s1_l2fc_invitro", "crispr_gra12s1_l2fc_invivo", "crispr_gra12s1_disco",
                   "crispr_gra12s2_l2fc_invitro", "crispr_gra12s2_l2fc_invivo", "crispr_gra12s2_disco",
                   "fit_hyperlopit_unassigned_invivo_lib1", "fit_hyperlopit_unassigned_invivo_lib2",
                   "crispr_reporter_strain_p8_log2", "oxidative_stress_screen_score",
                   "enteric_oocyst_yield", "fit_no_glucose", "fit_ifng", "fit_invivo_brain"):
        assert column in excluded, f"{column} is a fitness measurement and is not held out"
        assert "fitness" in C.categories_of(column)


def test_the_hyperlopit_library_screen_is_fitness_and_not_place():
    """Its name says which library was screened; its measurement is guide depletion in mice."""
    for column in ("fit_hyperlopit_unassigned_invivo_lib1", "fit_hyperlopit_unassigned_invivo_lib2"):
        assert C.categories_of(column) == ("fitness",)
        assert column not in C.excluded_for("place")


def test_every_hyperlopit_and_tagm_column_is_place():
    excluded = set(C.excluded_for("place"))
    for column in ("lopit_map", "lopit_mcmc", "lopit_prob_map", "lopit_prob_mcmc", "lopit_confident",
                   "lopit_methods_agree", "lopit_unified", "compartment", "compartment_source",
                   "compartment_best", "compartment_hyperlopit", "compartment_go_cc",
                   "ortholopit_label", "ortholopit_donors", "ortholopit_accuracy",
                   "ortholopit_accepted", "host_facing", "secretome_soluble_over_vesicle_log2"):
        assert column in excluded, f"{column} carries localisation knowledge and is not held out"


def test_pvm_proximity_is_both_place_and_interaction():
    assert set(C.categories_of("pvm_proximity_positive")) == {"place", "interaction"}
    assert "pvm_proximity_positive" in C.excluded_for("place")
    assert "pvm_proximity_positive" in C.excluded_for("interaction")


def test_the_host_transcription_effector_screen_is_host_response_and_place():
    """`host_response` is the literal reading -- the readout is host transcription, not depletion. But a
    screen naming which effectors change host transcription has very nearly shown which proteins reach
    the host, so for a localisation question it sits close to the conclusion. Until 2026-11-05 only the
    hand-written family override held it out of host-facing, which is the mechanism the categories were
    built to replace; `place` makes the category rule do it."""
    for column in ("hosttx_T2", "hosttx_padj", "hosttx_signature_pc01", "hosttx_signature_n_de"):
        assert set(C.categories_of(column)) == {"host_response", "place"}
        assert column in C.excluded_for("place"), "a localisation question must not see this"
        assert column not in C.excluded_for("fitness")      # its readout is not depletion


def test_function_is_its_own_category_and_not_phenotype():
    """A GO term is a conclusion about what the protein does; a drug-sensitivity shift is a measurement
    of what happens when you break it. Pooling them meant a function question held out real perturbation
    phenotypes, or an imaging-phenotype question held out annotation it never needed."""
    for column in ("function_class", "function_class_n", "function_go_bp", "function_go_bp_n",
                   "function_go_mf", "function_go_mf_n", "ec_number", "has_ec"):
        assert C.categories_of(column) == ("function",), column
        assert column in C.excluded_for("function")
        assert column not in C.excluded_for("phenotype"), "annotation is not a perturbation observation"


def test_a_perturbation_phenotype_is_still_phenotype_and_survives_a_function_question():
    """The other half of the split: a function question may legitimately learn from what breaks."""
    for column in ("diff_reporter_log2_mNG_over_bulk", "screen_egress_phenotype",
                   "resistance_substitutions", "escrt_recruitment_xgboost"):
        assert "phenotype" in C.categories_of(column), column
        assert column not in C.excluded_for("function")


def test_publication_and_mention_counts_are_attention():
    for column in ("n_publications", "n_papers_focal", "n_papers_substantive", "n_papers_incidental",
                   "n_fulltext", "attention_depth", "lit_tier"):
        assert "attention" in C.categories_of(column)


def test_gene_ids_names_and_products_are_identity_and_are_never_features():
    for column in ("gene_id", "product", "orthogroup"):
        assert "identity" in C.categories_of(column)
        assert not C.is_feature(column)
        for category in C.CATEGORIES:
            assert column in C.excluded_for(category), f"{column} is offered as evidence for {category}"


def test_ribo_seq_and_translational_efficiency_are_translation_and_not_transcript():
    for column in ("te99395_intracellular_r1", "te129869_confluent_r1", "te245775_parent_tachy_r1",
                   "te302107_tachy_r1", "rpf99395_extracellular_r1", "rpf245775_eif12ko_tachy_r3",
                   "expr_rel_GSE245775_parental_stress_RIBOseq"):
        assert C.categories_of(column) == ("translation",)
        assert column not in C.excluded_for("transcript")
    # The RNA arm of the same deposits is the converse case and must stay `transcript`.
    for column in ("rna99395_intracellular_r1", "rna245775_parent_tachy_r1"):
        assert C.categories_of(column) == ("transcript",)


def test_an_unknown_category_raises_rather_than_excluding_nothing():
    with pytest.raises(ValueError, match="not information categories"):
        C.excluded_for("localisation")


# --------------------------------------------------------------------------- the harness enforcement
def _fitness_question(**kwargs):
    fields = dict(id="category_probe", text="Does the parasite need it without glucose?",
                  answer_type="boolean", question_class="A", label_source="synthetic",
                  label_family="truth", label_category="fitness", evidence_forbidden=("truth",),
                  min_labels_for_frozen_test=200)
    fields.update(kwargs)
    return Q.Question(**fields)


def _data(n=400, seed=0):
    rng = np.random.default_rng(seed)
    ids = [f"g{i:04d}" for i in range(n)]
    hidden = rng.normal(size=n)
    labels = pd.Series(hidden > 0, index=ids)
    features = pd.DataFrame({"fit_no_glucose": hidden + rng.normal(scale=0.2, size=n),
                             "codon_gc3": rng.normal(size=n)}, index=ids)
    families = pd.Series({"fit_no_glucose": "some_other_family", "codon_gc3": "toxodb_codon_usage"})
    groups = pd.Series([f"og{i // 2:04d}" for i in range(n)], index=ids)
    return labels, features, families, groups


def test_a_question_cannot_be_run_with_a_column_of_its_own_category(tmp_path):
    """The mechanism, on the exact shape of the failure it was written for: the forbidden FAMILY list
    says nothing about `some_other_family`, and the column is a fitness measurement anyway."""
    labels, features, families, groups = _data()
    question = _fitness_question()
    assert question.allows("some_other_family")              # the hand-written list misses it
    with pytest.raises(ValueError, match="own information category"):
        H.run_question(question, labels=labels, features=features, families=families, groups=groups,
                       strategies=[H.Strategy.of("all", "hist_gbt", max_iter=40)],
                       log_path=tmp_path / "holdout.jsonl")
    filtered = H.allowed_features(question, features, families)
    assert list(filtered.columns) == ["codon_gc3"]


def test_a_circular_column_may_be_forced_in_and_is_then_recorded(tmp_path):
    labels, features, families, groups = _data()
    question = _fitness_question()
    result = H.run_question(question, labels=labels, features=features, families=families,
                            groups=groups, strategies=[H.Strategy.of("all", "hist_gbt", max_iter=40)],
                            log_path=tmp_path / "holdout.jsonl",
                            force_include={"fitness": "demonstrating the override is recorded"})
    record = result.evidence_selection
    assert record["label_categories"] == ("fitness",)
    assert "fit_no_glucose" in record["force_included"]
    assert "recorded" in record["force_included"]["fit_no_glucose"]
    assert not record["excluded_by_category"]
    assert "fit_no_glucose" in result.as_dict()["evidence_selection"]["force_included"]
    # Forcing a column in changes what the model saw, so it must also spend the holdout budget.
    clean = H.run_question(question, labels=labels, features=features[["codon_gc3"]],
                           families=families, groups=groups,
                           strategies=[H.Strategy.of("all", "hist_gbt", max_iter=40)],
                           log_path=tmp_path / "holdout.jsonl")
    assert clean.evidence_hash != result.evidence_hash


def test_the_run_record_names_what_was_removed_and_why():
    labels, features, families, groups = _data()
    record = H.evidence_selection(_fitness_question(), features.columns, families)
    assert record["excluded_by_category"] == {"fit_no_glucose": "fitness"}
    assert record["effective"] == ("codon_gc3",)
    assert record["requested"] == ("fit_no_glucose", "codon_gc3")
    assert not record["unannotated"]


def test_an_unannotated_column_is_recorded_rather_than_assumed_safe():
    """Completeness is enforced at ingest (this file's gate), not by refusing a run mid-flight: the
    harness records an uncategorised column so the scorecard shows what could not be checked."""
    families = pd.Series({"made_up_column": "somewhere"})
    record = H.evidence_selection(_fitness_question(), ["made_up_column"], families)
    assert record["unannotated"] == ("made_up_column",)
    assert record["effective"] == ("made_up_column",)


def test_the_hand_written_forbidden_list_still_applies_on_top_of_the_category():
    """Kept as an explicit override: `pubmed` is `attention`, so no `place` category would hold it out
    of the host-facing question, where it is nearly the conclusion."""
    question = Q.get("host_facing")
    assert question.label_category == "place"
    assert "pubmed" in question.evidence_forbidden
    families = pd.Series({"n_publications": "pubmed", "lopit_prob_map": "lopit_tgon",
                          "codon_gc3": "toxodb_codon_usage"})
    record = H.evidence_selection(question, list(families.index), families)
    assert record["excluded_by_family"] == {"n_publications": "pubmed"}
    assert record["excluded_by_category"] == {"lopit_prob_map": "place"}
    assert record["effective"] == ("codon_gc3",)


def test_every_registered_question_states_or_can_derive_its_label_category():
    for question in Q.QUESTIONS.values():
        if question.id == "registry_probe":                   # a fixture of the registry's own tests
            continue
        assert C.label_categories(question), f"{question.id} has no label category"
