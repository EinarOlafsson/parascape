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
| admitted on declared grounds | 115 | 74 per cent |
| needs reading | 41 | 26 per cent |
| **biologically ineligible** | **0** | **0 per cent** |

**No candidate now filling a slot is biologically ineligible on declared grounds.** Organism, unit
and mapping fit for all 156. That is a real result and it is the one worth stating first.

**31 of the 156 are inferred rather than measured** — curated readings of papers, computed features,
annotation, orthology. They are admitted as evidence and labelled, never counted as measurements.
Being inferred is a grade, not a verdict: curated reading is the only route to a Toxoplasma
resistance allele, because this organism has no resistome to parse.

### What the 41 need, and why

| what is unknown | rows | what has to be read |
|---|---:|---|
| stage | 27 | the slot asks about a life stage that the column does not name |
| access | 7 | a measurement with no recorded accession |
| access and stage | 7 | both of the above |

A further **28 rows record `level: reference`**, which states a resource rather than a point on the
central dogma. The level field does not speak to the quantity for those, but the assay does, and for
all 28 the assay establishes it — so they read **carried by the assay** and are admitted rather than
queued. An earlier version of this table queued them, which put 28 rows in front of a reader with
nothing specific to read. A queue padded with rows that need no work hides the rows that do.

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

## Asking the deposit, and the one real defect it found

A declaration check cannot resolve a stage the column does not name, so `resolve_stage_from_deposit.py`
asks the authority that knows: the submitter, through the per-sample characteristics in GEO. Same
principle the origin-paper check settled on — ask the deposit before the paper — applied to biology
instead of provenance.

It asks about **both** kinds of row: the ones with no stage, and the 12 that claim a stage from
their column name. Checking only the unknowns would have left 12 claims believed on the strength of
a column name, and a column named `tachy` whose samples say otherwise is worse than a column that
says nothing — the first is trusted, the second is queued.

| what the deposit said | rows |
|---|---:|
| confirmed the stage, from the samples the column selects | 10 |
| **declared a different stage than the slot asks for** | **1** |
| declares no stage in any sample, so the paper is the only authority left | 9 |
| several stages and the column cannot be mapped to samples | 4 |
| no deposit to ask — a paper supplementary table, or nothing recorded | 26 |

Ten stage claims are now backed by the submitter rather than by a column name, among them the two
`gse206344` columns, which the deposit settles as oocyst. Nine are unverifiable this way: GSE22258
names its columns `tachyzoite` and `bradyzoite`, but no sample characteristic in that series
declares a stage at all, so those two remain column-name claims.

**And it found one disagreement, which is a real defect.**


> `Tg_transcription_bradyzoite_tissue_cyst` includes the column pattern `rna108740_Day`. That
> pattern selects 14 samples of GSE108740, and the submitter declares every one of them
> **`developmental stage: Enteropithelial stages`**, host **`Felis catus`**. Those are the stages in
> the gut of the cat, not bradyzoites in a tissue cyst in a mouse.

The slot averages 12 columns. **Eight of the 12 are the cat enteroepithelial samples**; only
`rna108740_Tissue_cysts_A/B`, `expr_cyst` and `rna22258_bradyzoite` are the stage the slot asks for.

What it costs, measured against the same slot with the cat columns removed, over 7,865 genes:

| | |
|---|---|
| agreement of the two versions | Spearman 0.79 |
| genes moving more than 10 percentile points | 3,105 (39 per cent) |
| genes moving more than 25 percentile points | 973 (12 per cent) |

Two things say the inclusion was not deliberate. There is **no recorded rationale** anywhere for
putting the Day columns in a bradyzoite slot. And `expr_cyst`, the curated summary of the same
deposit, is documented as "the mean **tissue-cyst** FPKM" — the hand-built column excludes exactly
what the pattern pulls in. The pattern `rna108740_Day` simply over-matches.

**I have not changed the slot.** Removing the pattern changes the feature values every downstream
model sees for that slot, which is a scientific change with consequences for the frozen benchmarks,
so it is a decision to take deliberately rather than overnight. The recommendation is to drop
`rna108740_Day` from the patterns of `Tg_transcription_bradyzoite_tissue_cyst`, and — if the cat
stages are wanted as evidence, which they reasonably might be — to give them a slot of their own,
since no slot currently asks about the enteroepithelial stages.

## The same mistake, found a second way: patterns that over-match

A slot pattern is a prefix. It is written with particular columns in mind and keeps matching
whatever the deposit gains afterwards, so a slot quietly acquires columns nobody chose for it. None
of this is visible in the slot declaration, in the registry, or in any test.
`scripts/audit_pattern_overmatch.py` looks at what each pattern **actually selects** and asks
whether those columns describe one condition. Five cases:

**`Tg_protein_abundance_tachyzoite`, pattern `proteome_`, selects 15 columns.** Its facets say
`tachyzoite` and `unperturbed`, and the pattern takes:

* 3 untreated columns — what the slot asks for;
* **9 columns of AP2XII-1-depleted parasites**, which PRIDE titles "In vitro production of
  cat-restricted Toxoplasma pre-sexual stages by epigenetic reprogramming" — perturbed, and on the
  way to a different stage;
* **3 log2 fold-change columns**, averaged in with log2 abundances. The registry note for this very
  source says the abundance and fold-change columns "are distinct quantification types".

Each pattern is collapsed by row mean before the policy runs, so all 15 become one number.

**The measured cost is small, and saying so matters.** Against the same pattern restricted to the
untreated columns, over 3,005 genes: **Spearman 0.987**, with 126 genes (4 per cent) moving more
than 10 percentile points. The columns are already normalized and centred near zero, so mixing them
does not move the scale the way it would if they were raw. This is a real defect in what the slot
claims to be, and a minor one in what it returns — both halves of that are worth stating. It is not
in the same class as the cat-stage finding above, which moves 39 per cent of genes.

Separately, **102 genes take their tachyzoite protein abundance from `protein_ibaq_log2`**, the Pru
immunoprecipitation column, whose own note says "Enrichment, NOT a deep proteome; do not report as
proteome-wide". The `fill` policy reaches for it wherever the total proteome is absent.

**Three more**, all in the eIF1.2 depletion slots: `rpf245775_eif12ko_`, `te245775_eif12ko_` and
`rna245775_eif12ko_` each select three tachyzoite and three pre-bradyzoite columns, so the slot
averages across two stages. Whether that is intended is a reading, not a defect — it is listed
because a slot is normally meant to hold the stage constant.

**What this tool cannot do.** It catches only the half the column names give away. The cat-stage
defect is invisible to it: those columns are named `rna108740_Day3_CS4_FPKM`, and "Day" is a
timepoint, not a stage. Nothing in that name says cat, or enteroepithelial, or anything a word list
could match — only the submitter knew. The two scripts are complementary and neither replaces the
other.

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

---

## Both fixes, tried rather than described (2026-11-13)

`decisions/apply_65_04_slot_fixes.py` applies either or both, idempotently, and refuses rather than
guesses if a slot has moved on. I applied them to a scratch copy, ran the suite, measured, and
reverted. Two things came out that neither the report nor I had noticed, and both change the
decision:

**1. Fix A is not finished by itself.** Dropping `rna108740_Day` leaves its eight columns claimed by
no slot at all, and `test_the_column_partition_holds_for_every_species` fails: instruction 41
requires that no column is claimed twice **and none by nobody**. So the cat enteroepithelial columns
need a home in the same change -- either a slot of their own, which is what they deserve since no
slot currently asks about the enteroepithelial stages, or an explicit exclusion. The suite is
enforcing the thing I recommended in prose.

**2. The downloaded data release must not be patched in place.** `slots.CATALOG` is
`data_dir()/slots.json`, not the repository copy, so editing only the repo changes what ships at the
next release and nothing that runs now -- including the tests, which read the downloaded copy. My
first attempt edited the repo file, re-ran 141 tests, watched them pass and had proved nothing: the
slots under test were the unmodified ones. The obvious correction, writing both, is also wrong:
hand-editing the downloaded copy fails
`test_packaged_layout_matches_the_shared_matrix_and_gene_order`, because the release pins the
checksum of what it packaged. The release is meant to be rebuilt, not patched.

So the route is: change the repository source, give the orphaned columns a home, rebuild the data
release, then regenerate the atlas because the counts this slot feeds will move.

**Both fixes now carry their companion slots, and the prepared change passes.** Fix B orphaned
twelve columns for the same reason fix A orphaned eight, so it has the same obligation: the nine
AP2XII-1-depleted abundance columns and the three fold-changes go to TWO new slots rather than one,
because putting an abundance and a log2 fold change in a single slot is precisely the defect fix B
removes -- each pattern is collapsed by row mean, and a mean over two quantification types is not a
quantity.

With fix A, fix B and all three companion slots applied: **196 passed, 1 failed**, and the single
failure is `test_packaged_layout_matches_the_shared_matrix_and_gene_order` -- which fails because I
patched the downloaded copy in place to test it, not because of the change. The proper route
rebuilds the data release, which recomputes that checksum.

So the decision in front of you is a complete, tested change rather than a one-line edit with two
hidden consequences. That is what it took to find out.
