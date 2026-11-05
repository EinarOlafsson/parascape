"""The question registry: one row per question the information matrix asks of the data.

Not to be confused with `parascape.questions`, the shipped catalogue of instruction 45, which asks
what this *table* can be asked and runs each candidate through the leakage closure. This module is
the instruction 73 registry: the questions the information matrix answers per gene, each with its
labels, its forbidden evidence and whether it can afford a frozen test. The two will likely meet --
a question admitted here is a candidate there -- but they are different objects and the names are
kept apart deliberately.

Instruction 73 replaces the framing in which host-facing, function and partners were *the* traits:
they are three questions among hundreds. A question is therefore a first-class, declared object
rather than a notebook written from scratch, and this module is where it is declared.

Two things are recorded here that cost us real work when they were left implicit:

* **What a question may not be answered from.** `evidence_forbidden` always contains the family the
  labels themselves came from, because a question answered from its own label source measures
  nothing. The first host-facing run learnt this the expensive way -- its top feature was a screen
  that is nearly the conclusion -- so the registry refuses a question whose `label_family` is not
  held out.
* **Whether a frozen test partition is affordable at all.** A held-back partition costs labels. A
  98-gene question cannot spare them, so below `min_labels_for_frozen_test` the harness keeps nested
  cross-validation only, flags the question `no_frozen_test`, and the question may rank and
  prioritise but may **not** fill cells. That is a property of the question, recorded once, not a
  judgement made per run.

The three admission tests of instruction 73 (answerable, non-redundant, useful) are fields too, so
the catalogue can be audited rather than argued about.

    register(Question(id='x', text='...', answer_type='boolean', question_class='B', ...))
    catalogue()        # the whole registry as a table
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import pandas as pd

from . import information_categories as IC

ANSWER_TYPES = ("boolean", "category", "ranked_list", "number")
QUESTION_CLASSES = ("A", "B", "C", "D")
LABEL_TIERS = ("observed-direct", "observed-curated", "observed-screen", "observed-proximity", "none")
STATUSES = ("proposed", "admitted", "scored", "gap", "retired")

# A one-sided Clopper-Pearson lower bound of 0.90 needs about 22 held-out calls even if every call
# is correct, and about 65 at a true precision of 0.95. A test partition is 15 % of the labels and
# only some of those genes are called, so a question needs several hundred labels before a frozen
# partition can demonstrate anything at all. 400 is that arithmetic, rounded, and it is a default:
# each question may state its own.
DEFAULT_MIN_LABELS_FOR_FROZEN_TEST = 400


@dataclass(frozen=True)
class Question:
    """One question, with its labels, the evidence it may use and whether it can afford a frozen test.

    `evidence_allowed` empty means "every family except the forbidden ones"; listing families
    narrows it further. `label_family` must appear in `evidence_forbidden`, which is the registry's
    one structural guarantee: no question is answered from the source its own labels came from.

    `label_category` is the instruction 73 information category of the labels themselves, and it is
    the primary exclusion: the harness holds out **every column of that category**, whatever its
    family or study, which is what the hand-written `evidence_forbidden` list failed to do. The list
    stays, as a recorded additional override for what a category cannot express -- `pubmed` is
    `attention` and would not be excluded from a `place` question by category alone, though for
    host-facing it is nearly the conclusion. Left empty, the category is read off the column named in
    `label_source` (`information_categories.label_categories`), so forgetting to declare one does not
    silently remove the mechanism.
    """

    id: str
    text: str
    answer_type: str
    question_class: str
    label_source: str                       # where the labels are read from, e.g. 'nodes.parquet:fit_no_glucose'
    label_family: str                       # the evidence family that produced them; must be held out
    label_tier: str = "observed-screen"
    label_category: str = ""                # the information category of the LABELS (instruction 73)
    evidence_allowed: tuple = ()            # empty = everything not forbidden
    evidence_forbidden: tuple = ()           # the hand-written override, kept IN ADDITION to the category
    min_labels_for_frozen_test: int = DEFAULT_MIN_LABELS_FOR_FROZEN_TEST
    status: str = "proposed"
    answerable: bool = True                 # admission test 1: a truth set exists or can be curated
    non_redundant: bool = True              # admission test 2: not a restatement of a question already here
    useful: bool = True                     # admission test 3: the answer would change a decision
    redundant_with: tuple = ()
    note: str = ""

    def __post_init__(self):
        if not isinstance(self.id, str) or not self.id or not isinstance(self.text, str) or not self.text:
            raise ValueError("A question needs an explicit identifier and the question in words")
        if self.answer_type not in ANSWER_TYPES:
            raise ValueError(f"{self.id}: answer_type must be one of {list(ANSWER_TYPES)}")
        if self.question_class not in QUESTION_CLASSES:
            raise ValueError(f"{self.id}: question_class must be one of {list(QUESTION_CLASSES)}")
        if self.label_tier not in LABEL_TIERS:
            raise ValueError(f"{self.id}: label_tier must be one of {list(LABEL_TIERS)}")
        if self.status not in STATUSES:
            raise ValueError(f"{self.id}: status must be one of {list(STATUSES)}")
        if type(self.min_labels_for_frozen_test) is not int or self.min_labels_for_frozen_test < 0:
            raise ValueError(f"{self.id}: min_labels_for_frozen_test must be a non-negative integer")
        if not self.label_source or not self.label_family:
            raise ValueError(f"{self.id}: name the label source and the evidence family it belongs to")
        if self.label_category and self.label_category not in IC.CATEGORIES:
            raise ValueError(f"{self.id}: label_category must be one of {sorted(IC.CATEGORIES)}")
        if self.label_family not in self.evidence_forbidden:
            raise ValueError(f"{self.id}: the label family {self.label_family!r} must be in evidence_forbidden")
        overlap = set(self.evidence_allowed) & set(self.evidence_forbidden)
        if overlap:
            raise ValueError(f"{self.id}: families both allowed and forbidden: {sorted(overlap)}")
        if self.answerable and self.label_tier == "none":
            raise ValueError(f"{self.id}: a question with no label tier is not answerable; say so")

    @property
    def admitted(self) -> bool:
        """Whether the question passes all three admission tests of instruction 73."""
        return bool(self.answerable and self.non_redundant and self.useful)

    def allows(self, family: str) -> bool:
        """Whether a named evidence family may be used to answer this question."""
        if family in self.evidence_forbidden:
            return False
        return not self.evidence_allowed or family in self.evidence_allowed

    def frozen_test_possible(self, n_labels: int) -> bool:
        """Whether this many labels can afford a partition held back from everything."""
        return int(n_labels) >= self.min_labels_for_frozen_test

    def as_row(self) -> dict:
        """The registry row for the catalogue, with tuples flattened for printing."""
        row = asdict(self)
        for key in ("evidence_allowed", "evidence_forbidden", "redundant_with"):
            row[key] = ", ".join(row[key])
        row["admitted"] = self.admitted
        return row


QUESTIONS: dict[str, Question] = {}


def register(question: Question) -> Question:
    """Add a question; registering the same id again with a different definition is an error."""
    if question.id in QUESTIONS and QUESTIONS[question.id] != question:
        raise ValueError(f"{question.id} is already registered differently")
    QUESTIONS[question.id] = question
    return question


def get(question_id: str) -> Question:
    """One registered question, by id."""
    if question_id not in QUESTIONS:
        raise KeyError(f"no such question: {question_id!r}")
    return QUESTIONS[question_id]


def catalogue(questions=None) -> pd.DataFrame:
    """The registry as a table, one row per question, ordered by class then id."""
    rows = [q.as_row() for q in (QUESTIONS.values() if questions is None else questions)]
    frame = pd.DataFrame(rows)
    if not len(frame):
        return frame
    return frame.sort_values(["question_class", "id"]).reset_index(drop=True)


# The CRISPR fitness screens all measure depletion of a knockout from a pool. They are correlated
# measurements of the same kind of thing, so a conditional-importance question is held out from all
# of them rather than only from the column its own labels were thresholded from -- otherwise
# `fit_complete_medium_2025` answers `fit_no_glucose` and nothing has been learnt.
# The list is enumerated by a rule rather than by eye, because two passes of guessing missed a
# family each time: it is EVERY family in `nodes.parquet` contributing a column whose name begins
# with `crispr` or `fit_` -- a pooled knockout screen read out as depletion -- plus the curated
# enteric fitness measurements. Four of these do not have "crispr" in their family name
# (`gra12`, `gra17_synthlethal`, `invivo_platform`, `gse253884_5`) and the harness's own
# permutation importance and best-single-feature baseline are what surfaced them, one per run.
_CRISPR_FITNESS = ("crispr_glucose_limitation", "crispr_invitro", "crispr_invivo_composite",
                   "crispr_macrophage", "crispr_parasite_density", "crispr_serum_restriction",
                   "crispr_young2019", "second_background_fitness", "curated_enteric_fitness",
                   "gra12", "gra17_synthlethal", "invivo_platform", "gse253884_5")

for _question in (
        Question(
            id="required_without_glucose",
            text="Which genes does the parasite need in order to grow when glucose is withdrawn?",
            answer_type="boolean", question_class="A",
            label_source="nodes.parquet:fit_no_glucose <= -2.0",
            label_family="crispr_glucose_limitation", label_tier="observed-screen",
            label_category="fitness",
            evidence_forbidden=_CRISPR_FITNESS,
            status="admitted",
            note="Rung 0 of the strategy ladder carries the measurement for the 7,393 measured genes; "
                 "this question is the prediction problem for the ~750 the screen missed, and it is "
                 "what the harness is demonstrated on."),
        Question(
            id="screen_egress_phenotype",
            text="Which genes, when disrupted, give an egress phenotype in the imaging screen?",
            answer_type="boolean", question_class="B",
            label_source="nodes.parquet:screen_egress_phenotype",
            label_family="splitcas9_imaging_screen", label_tier="observed-screen",
            label_category="phenotype",
            evidence_forbidden=("splitcas9_imaging_screen",),
            status="admitted",
            note="318 labelled genes, below the frozen-test threshold: this question may rank and "
                 "prioritise, and may not fill cells."),
        Question(
            id="host_facing",
            text="Does this protein reach beyond the parasite plasma membrane, into the host cell?",
            answer_type="boolean", question_class="B",
            label_source="information_map:host_facing", label_family="host_facing_truth",
            label_tier="observed-curated", label_category="place",
            evidence_forbidden=("host_facing_truth", "pvm_proximity", "hosttx_effectors",
                                "organelle_surface_turboid", "pubmed", "lopit_tgon"),
            status="scored",
            note="98 reviewed genes, 65 positives. Scored twice (2026-11-04) as a cross-validated "
                 "estimate with no frozen test; the measured ceiling, 0.822 precision, is below the "
                 "0.90 bar, so the binding constraint is the evidence and not the label count."),
        Question(
            id="function_class",
            text="What coarse functional class does this protein belong to?",
            answer_type="category", question_class="B",
            label_source="results/function_truth_2026_11_04", label_family="function_truth",
            label_tier="observed-curated", label_category="phenotype",
            evidence_forbidden=("function_truth", "quickgo", "interpro", "toxodb_ec_numbers"),
            status="scored",
            note="242 genes over 15 classes; the learning curve is still rising, so here the label "
                 "count binds and curation helps."),
        Question(
            id="binding_partners",
            text="Which proteins is this protein most likely to bind?",
            answer_type="ranked_list", question_class="D",
            label_source="starpath crosslink pairs (1 - local FDR)", label_family="starpath_xlms",
            label_tier="observed-direct", label_category="interaction",
            evidence_forbidden=("starpath_xlms", "ipms_baits", "ipms_corpus_membership"),
            min_labels_for_frozen_test=800,
            status="scored",
            note="Cold-node pair partitions, both endpoints unseen. Cold-node scoring keeps only "
                 "~22 % of labelled pairs, so the frozen-test threshold is higher than the default."),
        Question(
            id="favourable_immune_response_human",
            text="Which genes provoke a favourable immune response in humans?",
            answer_type="boolean", question_class="C",
            label_source="none: no human screen exists", label_family="human_immune_screen",
            label_tier="none", label_category="fitness",
            evidence_forbidden=("human_immune_screen",),
            status="gap", answerable=False,
            note="Declared gap. The missing measurement is the human twin of the Wang design: a "
                 "genome-wide screen in IFN-gamma-activated human macrophages against matched "
                 "unactivated human cells. Listed so the absence is visible rather than predicted into."),
        Question(
            id="autophagy_involvement",
            text="Is this protein involved in autophagy?",
            answer_type="boolean", question_class="C",
            label_source="none: fewer than 10 experimentally annotated genes",
            label_family="quickgo", label_tier="none", label_category="phenotype",
            evidence_forbidden=("quickgo",),
            status="gap", answerable=False,
            note="Declared gap: curation (the host-facing review is the template) would have to come first."),
        Question(
            id="has_domain_ipr000719",
            text="Does this protein carry InterPro domain IPR000719?",
            answer_type="boolean", question_class="B",
            label_source="nodes.parquet:domain flags", label_family="interpro",
            label_tier="observed-direct", label_category="domain",
            evidence_forbidden=("interpro",),
            status="retired", useful=False,
            note="Fails admission test 3. It is measured, not predicted, and the answer changes no "
                 "decision: it is a feature and belongs on the evidence side. Kept registered as the "
                 "worked example of a rejected question."),
):
    register(_question)
