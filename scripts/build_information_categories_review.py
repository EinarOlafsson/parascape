#!/usr/bin/env python3
"""Write the reviewable derivation of every column's information category (instruction 73).

The annotation itself lives in `parascape/information_categories.py`, because it has to be importable
by the harness and by the data picker. This script is how a human reads it: one row per column per
shipped table, with the categories, how they were derived, the dataset key and level/kind behind
that, and whether the derivation rests on recorded provenance or on a judgement made in that module.

It writes `results/information_categories_<date>/columns.tsv`, `counts.tsv`, `judgements.tsv` and
`README.md`, and it prints the same counts, so a reviewer can disagree with a specific line rather
than with the idea. Run it again after adding a column to any shipped table: the manifest in
`information_categories` is checked against the real tables by
`tests/test_information_categories.py`, and this table is what says why each entry is what it is.

    python scripts/build_information_categories_review.py [--date 2026-11-05]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.getcwd())

import pandas as pd                                                           # noqa: E402

from parascape import information_categories as C                             # noqa: E402
from parascape import question_registry as Q                                  # noqa: E402

#: Judgement calls that a reviewer should look at first, each with what would change the answer.
#: These are statements about the vocabulary and about composite columns, not about single names, so
#: they belong in prose rather than in the table.
OPEN_JUDGEMENTS = """
* **`phenotype` is carrying functional annotation.** Instruction 73's vocabulary has no `function`
  category, so `function_class`, `function_go_bp`, `function_go_mf`, `ec_number` and `has_ec` are
  annotated `phenotype` -- what the protein is held to do. The exclusion this produces is right
  (inferring function excludes GO and EC) and slightly wide (inferring the egress screen phenotype
  also excludes GO). If Einar wants them separated, the vocabulary needs a seventeenth category and
  the gate test's vocabulary check is what will say so.
* **`n_holes` carries four categories.** It is a degree count in the graph of pairs that
  co-expression and co-fitness support while the literature does not co-mention them, so it is
  declared `interaction`, `attention`, `transcript` and `fitness`. A composite column is where one
  category would be a lie, and this is the only one in the map.
* **The host-transcription effector screen is `host_response` alone**, as instruction 73 says. It is
  near-conclusive for host-facing as well, and the host-facing question still holds it out through
  its hand-written `evidence_forbidden` list -- which is the clearest demonstration of why that list
  is kept as an override rather than deleted.
* **`domclass_rhoptry_kinase_ropk` is `domain` and `place`.** Its second arm matches
  `rhoptry|ROP\\d+` in the PRODUCT NAME, so the column reads a localisation word out of the
  annotation. Everything else in the domain table is `domain` alone.
* **Pooled CRISPR screens with a non-depletion readout** -- the differentiation reporter and the host
  ESCRT recruitment screen -- are annotated `phenotype` and also `fitness`. The ratio they report
  largely cancels depletion, so `fitness` is conservative rather than literal; it is kept because
  these are the kind of column that has bitten this project before.
* **Strain SNP counts are `conservation`, not `sequence`.** They measure constraint within the
  species. A learner wanting raw sequence features is not excluded by a `conservation` label and
  vice versa, which is the distinction intended.
* **`orthogroup` is `identity` and `conservation`.** It is the protected group key of every split,
  so it must never be a feature; it is also homology knowledge, and both are declared.
"""


def markdown(frame: pd.DataFrame) -> str:
    """A GitHub-flavoured table from a frame, written here because `tabulate` is not a dependency."""
    head = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    rows = ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join([head, rule, *rows])


def main() -> int:
    """Write the review directory and print the counts; returns a process exit status."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default=time.strftime("%Y_%m_%d"))
    parser.add_argument("--root", default="results")
    args = parser.parse_args()

    out = os.path.join(args.root, f"information_categories_{args.date}")
    os.makedirs(out, exist_ok=True)

    table = C.categories_table()
    table.to_csv(os.path.join(out, "columns.tsv"), sep="\t", index=False)
    counts = C.counts()
    counts.to_csv(os.path.join(out, "counts.tsv"), sep="\t", index=False)
    judged = table[table["confidence"] == "judged"]
    judged.to_csv(os.path.join(out, "judgements.tsv"), sep="\t", index=False)

    unannotated = C.unannotated()
    per_table = (table.groupby("table").size().rename("columns").reset_index())
    routes = table["derived_from"].value_counts().rename("columns").rename_axis("derived_from").reset_index()
    questions = pd.DataFrame([{"question": q.id, "label_category": q.label_category or "(derived)",
                               "derived": ",".join(C.label_categories(q)),
                               "evidence_forbidden (override)": ", ".join(q.evidence_forbidden)}
                              for q in Q.QUESTIONS.values()])
    questions.to_csv(os.path.join(out, "questions.tsv"), sep="\t", index=False)

    readme = f"""# Information categories: every data column annotated with the kind of knowledge it carries

Built by `scripts/build_information_categories_review.py` from `parascape/information_categories.py`
on {args.date.replace('_', '-')}. Instruction 73, "Information categories: the mechanism that stops
circular inference".

## The rule

**Every column in the map carries an information category: what kind of knowledge it embodies.**
Inferring a label of category X **excludes every column of category X by default**, whatever its
name, family or study. Including one is possible but requires an explicit reason that is recorded in
the run and shown on the scorecard, so a circular run is visible rather than impossible. A column may
carry more than one category where it genuinely does, and the exclusion then applies if **any** of
its categories matches the label's. `identity` columns are never offered as features at all.

Why it exists: a UMAP built from "all the data" and scored against held-out localisation labels is
circular if any input column carries `place`, and the circularity is in the *input selection*, not in
the model, so nothing downstream can see it. The per-question `evidence_forbidden` list that this
generalises already failed once -- `gra17_synthlethal`, `invivo_platform` and `gra12` are pooled
CRISPR fitness screens whose family names contain no "crispr", and a human reading an importance
table is what caught them.

## Where it is enforced

* **The harness.** `harness.evidence_selection` computes the exclusion from categories, not from the
  hand-written family list, and returns the requested and the effective selection with the reason for
  every removal. `run_question` refuses a circular column rather than dropping it quietly;
  `force_include={{'<column or category>': '<reason>'}}` includes one on purpose, records it in
  `HarnessResult.evidence_selection`, and hashes it into the evidence identity so it also spends the
  frozen-partition budget.
* **The hand-written list is kept as an additional, recorded override.** It expresses what a category
  cannot: `pubmed` is `attention`, so no `place` category would hold it out of the host-facing
  question, where it is nearly the conclusion. It is no longer the mechanism.
* **The gate.** `tests/test_information_categories.py` fails if any column of any shipped table has
  no category, if the vocabulary drifts from instruction 73, if the enumerated manifest and the real
  tables disagree, or if any of the named traps regresses.

## Counts

{len(C.SHIPPED_COLUMNS):,} distinct columns over {len(C.TABLES)} shipped tables, none unannotated
(`information_categories.unannotated()` is `{unannotated}`).

### Per table

{markdown(per_table)}

### Per category

A column with two categories counts once under each, so the total exceeds the column count; that is
the honest way to report it, since the exclusion fires on any match.

{markdown(counts)}

### How each column's category was derived

`provenance:level_kind` and `provenance:dataset` come from `datasets.provenance` -- the dataset's
on-disk level and kind (`transcription/RNAseq`, `DNA/CRISPR_screen`, `post_translation/LOPIT`), or
its own registry entry where the level is true but not decisive. `release_registry` comes from the
per-column registry shipped beside a trait table (`domain_lookup.tsv`, `column_registry.tsv`, the
`sequence_traits` build rules), which is evidence of the same kind kept beside the table.
`column_override` and `name_override` are judgements made in `information_categories.py`, each with
its reason in the table; every one of them is listed in `judgements.tsv`.

{markdown(routes)}

{len(table) - len(judged):,} of {len(table):,} rows rest on recorded provenance and {len(judged):,}
on a judgement.

## The traps, verified rather than assumed

| Trap | Verified |
|---|---|
| every pooled CRISPR fitness screen is `fitness` however its family is named | all 16 `DNA/CRISPR_screen` datasets were enumerated from the registry, not guessed: `gra17_synthlethal`, `invivo_platform` and `gra12` are `fitness`, as are `gse253884_5`, `second_background_fitness`, `oxidative_stress_screen`, `curated_enteric_fitness` and the seven `crispr_*` keys |
| hyperLOPIT and every TAGM-derived column is `place` | the ten `lopit_tgon` columns, plus the four `ortholopit_*` columns that no dataset entry claims and that a name rule would miss |
| `pvm_proximity_positive` is both `place` and `interaction` | yes, and the test asserts both |
| the host-transcription effector screen is `host_response` | the 24 `hosttx_*` columns, although their dataset's kind is `CRISPR_screen` |
| publication, mention and attention counts are `attention` | `n_publications`, `n_papers_focal`, `n_papers_substantive`, `n_papers_incidental`, `n_fulltext`, `attention_depth`, `lit_tier` |
| gene ids, names and products are `identity` | `gene_id`, `product`, `orthogroup`; none is ever a feature |
| translational-efficiency and Ribo-seq columns are `translation`, not `transcript` | every `te*` and `rpf*` column, and the three `expr_*_RIBOseq` columns the expression `column_registry.tsv` identifies; the **RNA arm** of the same deposits stays `transcript`, which is the converse error |
| `fit_hyperlopit_unassigned_invivo_lib*` | **`fitness`, not `place`**: the name records which library was screened, the measurement is guide depletion in mice. A name-based rule gets this one exactly backwards |

## Judgement calls, which are the part to argue with
{OPEN_JUDGEMENTS}
## The registered questions and their label categories

{markdown(questions)}

## Files

* `columns.tsv` -- one row per column per table: column, table, categories, derived_from, reason,
  provenance key, provenance level/kind, confidence.
* `judgements.tsv` -- the subset whose categories rest on a judgement rather than on provenance.
* `counts.tsv` -- columns per category, split by whether they came from provenance.
* `questions.tsv` -- each registered question's label category and its remaining hand-written override.
"""
    with open(os.path.join(out, "README.md"), "w") as stream:
        stream.write(readme)
    print(counts.to_string(index=False))
    print(routes.to_string(index=False))
    print(f"unannotated: {unannotated}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
