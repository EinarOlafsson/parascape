"""Which labels may stand beside which (instruction 73, rung 5)."""
from __future__ import annotations

from parascape import label_pairing as LP, question_registry as Q


def test_two_labels_of_the_same_category_may_not_be_learned_together():
    """A `place` target helped by another `place` label is the same question twice."""
    niches = [q for q in Q.QUESTIONS.values() if q.label_category == "place"]
    assert len(niches) >= 2
    allowed, reason = LP.may_stand_beside(niches[0], niches[1])
    assert not allowed and "place" in reason


def test_a_forbidden_family_is_refused_in_either_direction():
    """A question forbids what would restate its answer, so either side's refusal is enough."""
    frame = LP.pairs()
    refused = frame[~frame.allowed & frame.reason.str.contains("forbids")]
    assert len(refused) >= 2, "the domain question and the function question forbid each other"
    for row in refused.itertuples():
        assert "interpro" in row.reason or "forbids" in row.reason


def test_different_organisms_never_pair():
    """A joint model shares one row space, so two species cannot be learned together."""
    a = Q.Question(id="a", text="?", answer_type="boolean", question_class="B",
                   label_source="nodes.parquet:x", label_family="fa", label_category="fitness",
                   evidence_forbidden=("fa",), status="scored")
    b = Q.Question(id="b", text="?", answer_type="boolean", question_class="B",
                   label_source="pf_nodes.parquet:y", label_family="fb", label_category="place",
                   evidence_forbidden=("fb",), status="scored")
    # The organism is READ from the table named in label_source: the registry has no organism field.
    assert LP.organism_of(a) == "Tg" and LP.organism_of(b) == "Pf"
    allowed, reason = LP.may_stand_beside(a, b)
    assert not allowed and "organism" in reason


def test_every_pair_carries_a_reason_and_the_record_is_not_empty():
    """A refusal that says nothing cannot be argued with, so each row explains itself."""
    frame = LP.pairs()
    assert len(frame) > 100
    assert frame.reason.str.len().min() > 20
    assert frame.allowed.any() and not frame.allowed.all()
    assert not (frame.target == frame.auxiliary).any()


def test_auxiliaries_exclude_the_question_itself_and_its_own_category():
    """The practical call a joint model makes: who may help this question."""
    helpers = LP.auxiliaries_for("host_facing")
    assert helpers and "host_facing" not in helpers
    place = {q.id for q in Q.QUESTIONS.values() if q.label_category == "place"}
    assert not (set(helpers) & place), "a place question must not be helped by another place label"
