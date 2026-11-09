"""Which labels may stand beside which, so a joint model can be tested (instruction 73, rung 5).

Rung 5 of the strategy ladder is a model that learns several questions at once. It was recorded as untested for
one reason, in `borrowed_models.UNTESTED`: the harness excludes every column of a target's own information
category, and an **auxiliary label** of that category would bring the excluded evidence back in by another door.
Nobody had written down which labels may stand beside which, so the strategy could not be graded honestly.

This module is that record, and it is derived rather than asserted. Two labels may be learned together when:

* their information categories do not intersect -- a `place` target may not be helped by another `place` label,
  because the second label is the first question in different words;
* neither label's dataset family is on the other's forbidden list -- the hand-written `evidence_forbidden` of a
  question names the families that would restate its answer, and a label drawn from one of them is that answer;
* they are about the same organism, because a joint model shares one row space -- and the organism is read from
  the table the label comes from, since the registry records no organism field.

Everything else is allowed, and the pair carries the reason it was allowed. A pair that fails says which rule it
failed, so the record is reviewable rather than a verdict.
"""
from __future__ import annotations

import pandas as pd

from . import question_registry as Q


def organism_of(question) -> str:
    """Which space a question's label lives in, read from the table its label comes from.

    The registry records no organism field, so the table named in `label_source` is the only statement of it:
    `pf_nodes.parquet` is Plasmodium, anything else Toxoplasma. Derived rather than assumed, because defaulting
    every question to one organism would let two species pair silently.
    """
    from . import organisms

    source = str(getattr(question, "label_source", "") or "").lower()
    declared = str(getattr(question, "organism", "") or "")
    if declared:
        return declared
    # Longest table name first: "nodes.parquet" is a substring of "pf_nodes.parquet", so matching in
    # registration order called every Plasmodium question Toxoplasma.
    for code, space in sorted(organisms.SPACES.items(), key=lambda kv: -len(kv[1].nodes or "")):
        if space.nodes and space.nodes.lower() in source:
            return code
    return organisms.TOXOPLASMA


def _categories(question) -> frozenset:
    """The information categories a question's label carries, from the question or its column."""
    from . import information_categories as IC

    declared = str(getattr(question, "label_category", "") or "")
    if declared:
        return frozenset(part.strip() for part in declared.split(",") if part.strip())
    source = str(getattr(question, "label_source", "") or "")
    column = source.split(":")[-1] if ":" in source else ""
    return frozenset(IC.categories_of(column, organism_of(question))) if column else frozenset()


def _families(question) -> frozenset:
    """The dataset families a question refuses as evidence, including its own label family."""
    forbidden = set(getattr(question, "evidence_forbidden", ()) or ())
    own = str(getattr(question, "label_family", "") or "")
    if own:
        forbidden.add(own)
    return frozenset(forbidden)


def may_stand_beside(target, auxiliary) -> tuple:
    """(may it, why) for one ordered pair of questions.

    `target` is the question being answered and `auxiliary` the label offered alongside it. The order matters:
    a question forbids families that would restate ITS answer, so both directions are checked and either
    refusal is enough.
    """
    if target.id == auxiliary.id:
        return False, "a label cannot stand beside itself"
    if organism_of(target) != organism_of(auxiliary):
        return False, (f"different organisms ({organism_of(target)} and {organism_of(auxiliary)}): "
                       f"no shared row space")
    shared = _categories(target) & _categories(auxiliary)
    if shared:
        return False, (f"both labels are {', '.join(sorted(shared))}: the auxiliary label is the target "
                       f"question in different words")
    if str(auxiliary.label_family or "") in _families(target):
        return False, (f"the auxiliary label comes from {auxiliary.label_family!r}, which the target forbids "
                       f"as evidence")
    if str(target.label_family or "") in _families(auxiliary):
        return False, (f"the target's label comes from {target.label_family!r}, which the auxiliary question "
                       f"forbids as evidence")
    target_categories = ", ".join(sorted(_categories(target))) or "uncategorised"
    auxiliary_categories = ", ".join(sorted(_categories(auxiliary))) or "uncategorised"
    return True, (f"{target_categories} beside {auxiliary_categories}: different categories, neither family "
                  f"forbidden by the other")


def pairs(questions=None) -> pd.DataFrame:
    """Every ordered pair of registry questions, with whether they may be learned together and why.

    One row per (target, auxiliary). `allowed` is the answer and `reason` is always filled, so a refusal can be
    read and argued with rather than guessed at.
    """
    questions = list((questions or Q.QUESTIONS).values()) if isinstance(questions or Q.QUESTIONS, dict) \
        else list(questions)
    rows = []
    for target in questions:
        for auxiliary in questions:
            if target.id == auxiliary.id:
                continue
            allowed, reason = may_stand_beside(target, auxiliary)
            rows.append(dict(target=target.id, auxiliary=auxiliary.id,
                             target_category=", ".join(sorted(_categories(target))),
                             auxiliary_category=", ".join(sorted(_categories(auxiliary))),
                             allowed=allowed, reason=reason))
    return pd.DataFrame(rows)


def auxiliaries_for(question_id: str, questions=None) -> list:
    """The question ids that may stand beside this one, best used as the auxiliary tasks of a joint model."""
    frame = pairs(questions)
    return sorted(frame.loc[(frame.target == question_id) & frame.allowed, "auxiliary"])
