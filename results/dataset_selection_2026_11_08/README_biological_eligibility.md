# 65.04 — the biological eligibility read

Run 13 November 2026 by `scripts/verify_biological_eligibility.py`. Offline and deterministic: it
reads the shipped registry, the slot declarations and the audit tables already in this folder, so
re-running it reproduces the table exactly. Output: `biological_eligibility.tsv`, one row for each
of the 156 candidates that currently fill a slot, covering 66 sources across 60 slots.

## Why this step comes before the ranking

Instruction 65 ranks datasets by citation rate, newness and comprehensiveness, and it says plainly
that **biological admission comes first**: organism, host or cell type, stage, assay, quantity,
units, measured against inferred evidence, mapping, quality checks and data access must fit the
question. Ranking a source that does not answer the question is worse than not ranking it, because
the number that comes out looks like a comparison.

The check compares **two declarations** rather than my reading of a paper:

* what the slot asks — `axis` (the quantity), `unit` (what one row is), `context` and the three
  facet paths, and `role`;
* what the source declares — `organism`, `kind` (the assay), `level` (the quantity class),
  `coverage`, `accession` and its column names.

Where the two disagree, that is a finding whatever any paper says. Where they are silent, the
verdict is **unknown** and the row is queued for reading. It is never guessed and never counted as a
pass — the rule instruction 65 sets for dates and counts, which matters more for biology.

## The answer

| | candidates | share |
|---|---:|---:|
| admitted on declared grounds | 87 | 56 per cent |
| needs reading | 69 | 44 per cent |
| **biologically ineligible** | **0** | **0 per cent** |

**No candidate now filling a slot is biologically ineligible on declared grounds.** Organism, unit
and mapping fit for all 156. That is a real result and it is the one worth stating first.

**31 of the 156 are inferred rather than measured** — curated readings of papers, computed features,
annotation, orthology. They are admitted as evidence and labelled, never counted as measurements.
Being inferred is a grade, not a verdict: curated reading is the only route to a Toxoplasma
resistance allele, because this organism has no resistome to parse.

### What the 69 need, and why

| what is unknown | rows | what has to be read |
|---|---:|---|
| quantity | 28 | the source records `level: reference`, which states a resource rather than a quantity, so the axis cannot be checked from the declaration |
| stage | 27 | the slot asks about a life stage that the column does not name |
| access | 7 | a measurement with no recorded accession |
| access and stage | 7 | both of the above |

## What I got wrong first, and what it cost

The first version of this check reported **36 candidates as biologically ineligible. Every one was a
correct pairing and the table was wrong.** The capability model was too narrow:

* BioID was allowed to report only relations — but a bait anchored in an outer membrane reports what
  is near it, which is localization evidence about the outside of an organelle.
* Proteomics was not allowed to answer a chemistry question — but thermal proteome profiling is
  target engagement, and that is chemistry.
* Whole-genome sequencing of evolved clones was not allowed to be sequence evidence.
* Curated reading of papers was not allowed to carry chemistry or fitness.
* RNA-seq was confined to steady-state abundance, so decay after transcriptional arrest and the
  architecture of an untranslated region both came out as mismatches.

A too-strict admission table rejects good evidence with a confident verdict. That is the same
failure instruction 65 exists to prevent, pointing the other way, and it is harder to notice because
a rejection reads as rigour. Every widening in the script now names the source that justified it.

Three further traps, each of which produced a confident wrong answer before it was found:

1. **The channel is chosen by the column, not by the assay.** A ribosome-profiling deposit ships the
   footprint channel and the matched RNA channel. Reading the assay alone called five correct
   candidates a mismatch, among them `riboseq_mrna_`, which says what it is in its own name.
2. **Prose that mentions a stage may be denying it.** `gse245775` matched "bradyzoite" from its
   note — and that note exists to say the opposite, that the arms are labelled `pre-bradyzoites`,
   "named for that and not for bradyzoite". A word search found the word inside the sentence written
   to exclude it. A column name now counts as a declaration; prose counts as something to read.
3. **A candidate must not inherit a stage from a sibling column.** `rna108740_Day` was credited with
   "bradyzoite" because `rna108740_Tissue_cysts` sits in the same deposit. A neighbouring column is
   not evidence about this one, which is the whole reason slots are filled by pattern rather than by
   deposit.

All four are regression tests in `tests/test_biological_eligibility.py` (18 tests), and each was
confirmed to fail when the guard is removed.

## One finding about the registry itself

Two assays are spelled more than one way across 13 entries: `Ribo-seq` / `RiboSeq` / `riboseq`, and
`IP-MS` / `IPMS`. Anything that groups or filters by assay without folding them treats one assay as
two or three, and then compares groups that were never apart. The script folds them and records the
drift in `assay_spelling_drift.json`; a test fails if a third assay starts to drift.

## What this does not cover

* It admits on **declared** grounds. A source whose declaration is accurate and whose data is not
  passes here; that is what the 69 queued readings and the separate origin-paper and supplementary
  table checks are for.
* It examines the 156 candidates that currently fill a slot, not every possible pairing.
* Cp declares 63 slots and the registry holds no Cp source, so no Cp candidate appears in this
  table. That is a gap in coverage, not a judgement about eligibility.
* It says nothing about which of several eligible sources is best. That is the ranking step, and it
  only runs on what passes here.
