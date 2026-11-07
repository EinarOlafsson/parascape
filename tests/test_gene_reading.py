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
