"""The instruction 73 question registry refuses the definitions that cost us a run.

`parascape.question_registry`, not `parascape.questions` -- that module is instruction 45's shipped
catalogue of what the table can be asked, and `tests/test_questions.py` covers it.

Three of these checks defend a lesson paid for in work rather than a coding convention: a question
must hold out the family its own labels came from (the first host-facing run's top feature was
nearly the conclusion); a question that cannot afford a frozen partition must say so in the registry
rather than have it decided per run; and a feature masquerading as a question must be visible as
rejected rather than quietly absent.
"""
import pytest

from parascape import question_registry as Q


def test_unknown_answer_type_and_class_are_refused():
    base = dict(id="x", text="Is it so?", label_source="s", label_family="f",
                evidence_forbidden=("f",))
    with pytest.raises(ValueError, match="answer_type"):
        Q.Question(answer_type="probability", question_class="B", **base)
    with pytest.raises(ValueError, match="question_class"):
        Q.Question(answer_type="boolean", question_class="E", **base)
    with pytest.raises(ValueError, match="label_tier"):
        Q.Question(answer_type="boolean", question_class="B", label_tier="guess", **base)
    with pytest.raises(ValueError, match="status"):
        Q.Question(answer_type="boolean", question_class="B", status="maybe", **base)
    assert Q.Question(answer_type="boolean", question_class="B", **base).answer_type == "boolean"


def test_a_question_may_not_be_answered_from_its_own_label_source():
    with pytest.raises(ValueError, match="must be in evidence_forbidden"):
        Q.Question(id="x", text="Is it so?", answer_type="boolean", question_class="B",
                   label_source="screen", label_family="splitcas9_imaging_screen",
                   evidence_forbidden=("pubmed",))


def test_a_family_cannot_be_allowed_and_forbidden_at_once():
    with pytest.raises(ValueError, match="both allowed and forbidden"):
        Q.Question(id="x", text="Is it so?", answer_type="boolean", question_class="B",
                   label_source="screen", label_family="f",
                   evidence_allowed=("f", "g"), evidence_forbidden=("f",))


def test_allows_reads_both_lists_and_defaults_to_everything_not_forbidden():
    open_question = Q.Question(id="x", text="?", answer_type="boolean", question_class="B",
                               label_source="s", label_family="f", evidence_forbidden=("f", "g"))
    assert open_question.allows("anything") and not open_question.allows("g")
    narrow = Q.Question(id="y", text="?", answer_type="boolean", question_class="B",
                        label_source="s", label_family="f", evidence_allowed=("h",),
                        evidence_forbidden=("f",))
    assert narrow.allows("h") and not narrow.allows("anything")


def test_the_frozen_test_threshold_is_a_property_of_the_question():
    question = Q.Question(id="x", text="?", answer_type="boolean", question_class="B",
                          label_source="s", label_family="f", evidence_forbidden=("f",),
                          min_labels_for_frozen_test=400)
    assert not question.frozen_test_possible(318)
    assert question.frozen_test_possible(400)
    with pytest.raises(ValueError, match="min_labels_for_frozen_test"):
        Q.Question(id="x", text="?", answer_type="boolean", question_class="B", label_source="s",
                   label_family="f", evidence_forbidden=("f",), min_labels_for_frozen_test=-1)


def test_the_three_admission_tests_are_recorded_and_combine():
    useless = Q.get("has_domain_ipr000719")
    assert not useless.useful and not useless.admitted and useless.status == "retired"
    assert Q.get("required_without_glucose").admitted


def test_register_refuses_a_different_definition_under_the_same_id():
    question = Q.Question(id="registry_probe", text="?", answer_type="boolean", question_class="B",
                          label_source="s", label_family="f", evidence_forbidden=("f",))
    Q.register(question)
    Q.register(question)                       # idempotent for the identical definition
    with pytest.raises(ValueError, match="already registered differently"):
        Q.register(Q.Question(id="registry_probe", text="Something else?", answer_type="boolean",
                              question_class="B", label_source="s", label_family="f",
                              evidence_forbidden=("f",)))
    del Q.QUESTIONS["registry_probe"]


def test_a_gap_question_declares_itself_unanswerable_rather_than_being_predicted_into():
    gap = Q.get("favourable_immune_response_human")
    assert gap.question_class == "C" and gap.status == "gap"
    assert not gap.answerable and not gap.admitted and gap.label_tier == "none"
    with pytest.raises(ValueError, match="not answerable"):
        Q.Question(id="z", text="?", answer_type="boolean", question_class="C", label_source="none",
                   label_family="f", label_tier="none", evidence_forbidden=("f",), answerable=True)


def test_the_catalogue_is_a_table_with_one_row_per_question_and_every_class_present():
    table = Q.catalogue()
    assert len(table) == len(Q.QUESTIONS)
    assert set(table["id"]) == set(Q.QUESTIONS)
    assert set(Q.QUESTION_CLASSES) <= set(table["question_class"])
    for column in ("id", "text", "answer_type", "label_source", "label_tier", "evidence_allowed",
                   "evidence_forbidden", "question_class", "min_labels_for_frozen_test", "status"):
        assert column in table.columns
    assert table["question_class"].is_monotonic_increasing


def test_every_registered_question_holds_out_its_own_label_family():
    for question in Q.QUESTIONS.values():
        assert question.label_family in question.evidence_forbidden
        assert not question.allows(question.label_family)
