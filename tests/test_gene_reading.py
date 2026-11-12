"""Reading one gene (the Gene module, instructions 93/96): the synthesis says what is known with
its sources, keeps predictions out of it, and says what is not known."""
import pandas as pd

from parascape import gene_reading as GR


class FakeMap:
    """The three PublishedMap calls the reading uses, over a hand-made gene."""

    def gene_facts(self, gene):
        return pd.DataFrame([
            dict(column="host_facing", tier=1, table="information map", observed=True,
                 description="faces the host", value=True, probability=0.95,
                 source_tier="observed-direct"),
            dict(column="compartment_hyperlopit", tier=1, table="information map", observed=True,
                 description="hyperLOPIT niche", value="dense granules", probability=0.9,
                 source_tier="observed-curated"),
            dict(column="dom_IPR000008", tier=2, table="domains", observed=True,
                 description="C2 domain", value=False, probability=float("nan"),
                 source_tier="measured"),
            dict(column="essential_call", tier=1, table="information map", observed=True,
                 description="essential", value="yes", probability=0.6,
                 source_tier="filled: strategy 07"),
        ])

    def gene_evidence(self, gene):
        return pd.DataFrame([
            dict(gene_id=gene, trait="host_facing", value="yes", probability=0.95,
                 source_tier="observed-direct", reference="PMID:23601110", dataset="review"),
        ])

    def gene_partners(self, gene):
        return pd.DataFrame([
            dict(partner="PP2A-B55", side="host", probability=0.9,
                 source_tier="observed-curated", reference="PMID:23739762"),
            dict(partner="TGME49_270250", side="parasite", probability=0.97,
                 source_tier="observed-direct", reference="PMC12505969"),
        ])


def _read():
    return GR.read(FakeMap(), "TGME49_208830", "dense granule protein GRA16")


def test_where_it_is_names_the_place_and_its_sources():
    text = dict(GR.synthesis(_read()))["place"]
    assert "faces the host" in text and "dense granules" in text
    assert "pubmed.ncbi.nlm.nih.gov/23601110" in text, "each claim links its source"


def test_partners_are_counted_by_side_strongest_first():
    text = dict(GR.synthesis(_read()))["interaction"]
    assert "1 in the host, 1 in the parasite" in text
    assert text.index("TGME49_270250") < text.index("PP2A-B55")


def test_a_prediction_is_never_phrased_as_knowledge():
    reading = _read()
    assert "essential" not in " ".join(s for _c, s in GR.synthesis(reading))
    inferred = " ".join(GR.inferences(reading))
    assert "Predicted, not measured" in inferred and "60% certain" in inferred


def test_an_absent_domain_is_not_something_known():
    assert "domain" not in dict(GR.synthesis(_read()))


def test_key_traits_without_a_value_are_listed_as_not_known():
    class Empty(FakeMap):
        def gene_facts(self, gene):
            return pd.DataFrame(columns=["column", "tier", "table", "description", "observed",
                                         "value", "probability", "source_tier"])
    reading = GR.read(Empty(), "TGME49_000000")
    assert len(GR.inferences(reading)) == len(GR.KEY_TRAITS)


def test_gene_table_measurements_join_the_published_facts():
    row = pd.Series({"gene_id": "TGME49_208830", "fit_invitro_hff": 2.28, "empty": float("nan")})
    reading = GR.read(FakeMap(), "TGME49_208830", node_row=row)
    fitness = dict(GR.synthesis(reading)).get("fitness", "")
    assert "2.28" in fitness and "Sidik" not in fitness or "PMID" in fitness


def test_sequences_are_shown_by_length():
    assert GR.value_text("seq_full", "MAAAK") == "5 residues"


# --------------------------------------------------------------------------- what is unusual (instruction 60)
def _unusual_nodes():
    import numpy as np
    rng = np.random.default_rng(0)
    n = 800
    frame = pd.DataFrame({"gene_id": [f"g{i}" for i in range(n)],
                          "expr_a": rng.normal(size=n), "expr_b": rng.normal(size=n),
                          "flag": rng.random(n) > 0.5, "thin": [1.0] * 10 + [None] * (n - 10)})
    frame.loc[0, "expr_a"] = 99.0       # g0 is the highest on expr_a
    frame.loc[0, "expr_b"] = 0.0        # and ordinary on expr_b
    return frame


def test_the_most_extreme_measurement_comes_first_as_an_exact_rank():
    frame = GR.unusual(_unusual_nodes(), "g0")
    assert frame.iloc[0]["column"] == "expr_a" and frame.iloc[0]["rank_high"] == 1
    line = GR.unusual_sentences(frame.head(1))[0]
    assert "the highest of the 800 genes measured" in line and "100" not in line


def test_flags_and_thinly_measured_columns_are_not_ranked():
    frame = GR.unusual(_unusual_nodes(), "g0")
    assert "flag" not in set(frame["column"]) and "thin" not in set(frame["column"])


def test_an_unknown_gene_reads_as_nothing_unusual():
    assert GR.unusual(_unusual_nodes(), "nope").empty


def test_every_catalogue_method_is_accounted_for_in_a_gene_profile():
    """Instruction 64.30: all strategies listed; those without a per-gene record say why."""
    from parascape import strategies as S
    profile = GR.strategy_profile("TGME49_259020", ledger=pd.DataFrame(columns=[
        "gene_id", "mode", "strategy", "target", "prediction", "truth", "correct", "abstained"]))
    assert len(profile) == len(S.catalog())
    assert (profile["status"] != "").all() and not (profile["status"] == "held-out calls").any()


def test_combined_answers_say_whether_they_were_admitted():
    """64.35: an admitted label shows the outer-test call; a failed admission says so, never a call."""
    lines = GR.meta_lines("TGME49_259020")
    if not lines:
        return
    text = " ".join(lines)
    assert "not admitted" in text                      # cell-cycle phase failed its admission
    assert "made with this gene's answer hidden" in text


def test_a_record_built_on_an_older_answer_is_flagged_stale():
    """64.30: when the table's value no longer matches what the hold-out record was scored against."""
    import os
    from parascape import track_record as T
    ledger = T.shipped("Tg")
    row = ledger[(ledger["mode"].astype(str) == "together") & (ledger["target"].astype(str) == "dtm_class")].iloc[0]
    nodes = pd.read_parquet(os.path.join(os.path.dirname(os.path.dirname(__file__)), "parascape", "data", "nodes.parquet"))
    assert GR.stale_records(str(row["gene_id"]), nodes) == []
    changed = nodes.copy()
    changed.loc[changed["gene_id"].astype(str) == str(row["gene_id"]), "dtm_class"] = "SOMETHING_ELSE"
    assert any("dtm class" in s for s in GR.stale_records(str(row["gene_id"]), changed))


def test_a_gene_in_a_frozen_claim_set_shows_its_frozen_prediction():
    """64.38: the frozen prediction is visible on the gene's card, with its freeze date and hash."""
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results", "prospective_claims_2026_11_06", "claims.tsv")
    if not os.path.exists(path):
        return
    gene = pd.read_csv(path, sep="\t")["gene_id"].iloc[0]
    line = GR.prospective_lines(gene)[0]
    assert "Frozen prediction" in line and "rank 1 of" in line and "SHA-256" in line


# --------------------------------------- where in the parasite the gene is (98, WP8)
#
# The gene card draws the parasite with the organelle the gene is found in filled in --
# the schematic Einar named by name. These are the pure functions behind that drawing,
# tested on a made-up table so the chance baseline can be checked by counting.

DRAWABLE = ("dense granules", "rhoptries 1", "nucleus - chromatin", "apicoplast")


def placement_table():
    """A gene table with one placement column, one that is not, and one unplaced gene."""
    return pd.DataFrame({
        "gene_id": ["A", "B", "C", "D", "E"],
        "compartment": ["dense granules", "nucleus - chromatin", "nucleus - chromatin",
                        "unassigned", "rhoptries 1"],
        "compartment_source": ["hyperLOPIT", "hyperLOPIT", "orthoLOPIT", "unknown",
                               "hyperLOPIT"],
        "fitness": [1.0, 2.0, 3.0, 4.0, 5.0],
        "product": ["a", "b", "c", "d", "e"],
    })


def test_a_placement_column_is_found_by_what_is_in_it_not_by_its_name():
    """`compartment_best` looks like a placement and holds coarser groupings the artwork
    has no shape for; picking columns by name would have drawn nothing."""
    nodes = placement_table()
    nodes["compartment_best"] = ["nucleus", "nucleus", "nucleus", "unassigned", "apical"]
    found = GR.placement_columns(nodes, DRAWABLE)
    assert found == ["compartment"]


def test_a_column_the_drawing_mostly_cannot_colour_is_not_a_placement():
    """"Does anything in it match" kept `compartment_best`, which has a cytosol and an
    apicoplast among four fifths of groupings the artwork has no shape for -- so choosing
    it drew an empty parasite for most genes."""
    nodes = placement_table()
    nodes["mostly_coarse"] = ["apicoplast", "nucleus", "nucleus", "apical", "pellicle"]
    assert "mostly_coarse" not in GR.placement_columns(nodes, DRAWABLE)


def test_a_gene_with_no_measured_compartment_is_unknown_and_never_implied_no():
    place = GR.placement(placement_table(), "D", "compartment", "compartment_source")
    assert not place.known and place.compartment == ""
    assert "unknown" in GR.placement_sentence(place)
    assert "not the same as it being nowhere" in GR.placement_sentence(place)


def test_a_placement_names_the_organelle_and_what_measured_it():
    place = GR.placement(placement_table(), "A", "compartment", "compartment_source")
    assert place.known and place.compartment == "dense granules"
    assert place.source == "hyperLOPIT"
    assert "dense granules" in GR.placement_sentence(place)
    assert "hyperLOPIT" in GR.placement_sentence(place)


def test_the_chance_baseline_is_the_commonest_answer_with_its_n():
    """Instruction 98: every suggestion reported against chance, with n, as whole numbers."""
    place = GR.placement(placement_table(), "A", "compartment", "compartment_source")
    assert (place.measured, place.total) == (4, 5)          # D is unassigned
    assert place.commonest == "nucleus - chromatin" and place.commonest_n == 2
    assert place.coverage_percent == 80 and place.chance_percent == 50
    trust = GR.placement_trust(place)
    assert "4</b> of the 5" in trust and "50 %" in trust and "80 %" in trust
    assert "." not in trust.split("%")[0][-4:], "percentages are whole numbers"


def test_a_table_with_no_placement_at_all_reads_as_unknown_rather_than_raising():
    nodes = placement_table().drop(columns=["compartment"])
    assert GR.placement_columns(nodes, DRAWABLE) == []
    place = GR.placement(nodes, "A", "compartment")
    assert not place.known and GR.placement_trust(place) == ""


def test_the_window_colouring_wins_when_it_is_a_placement():
    """The card and the map must never show two different compartments for one gene."""
    columns = ["compartment", "lopit_map"]
    assert GR.preferred_placement_column(columns, "lopit_map") == "lopit_map"
    assert GR.preferred_placement_column(columns, "fitness") == "compartment"
    assert GR.preferred_placement_column(["lopit_map"], "") == "lopit_map"
    assert GR.preferred_placement_column([], "anything") == ""
