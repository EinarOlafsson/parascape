# 73 · The question matrix: 100–1000 trait questions, each with an answer, a probability and an uncertainty

**Status: open (requested 2026-11-05).** Supersedes the framing of instruction 72, which treated host-facing,
function and partners as *the* three traits. They are three questions among hundreds.

## The request (Einar, 2026-11-05)
* "I do not only care about host-facing parasite proteins. I care about systematically generating as much
  knowledge with the information available today on Toxoplasma from the Toxoplasma literature. Whether a protein
  is host-facing or not is just one trait. There are hundreds or thousands of information columns I would like to
  systematically predict, with assigned uncertainty to each prediction."
* His examples: importance for nutrient acquisition, and under low nutrient conditions; cell division; host-cell
  invasion; egress; switching to bradyzoite and from sporozoite; dissemination; survival in the liver;
  bradyzoite survival; a favourable immune response in mice, and in humans; the cell cycle; translation; protein
  scavenging; transcription; autophagy; methylation; acetylation; phosphorylation; lipid biogenesis and
  scavenging; the most likely binding partners; what a protein is expressed with; what it is expressed in
  response to. "A very small sample."
* **Answer types differ**: some boolean, some a category, some a ranked list (partners). Each needs a
  probability *and* an uncertainty.
* "The models will keep getting refined the more data comes in, but the actual architecture we should be able to
  set up already."

## The finding that reshapes the problem
**286 of the 411 numeric columns in `nodes.parquet` are measured for 5,000 or more genes.** Most of his questions
are therefore not prediction problems for most genes — the answer is already measured. The questions divide
into four classes, and only one of them is hard:

| Class | Examples from his list | Labels available | What the architecture must do |
|---|---|---|---|
| **A · Measured genome-wide** | nutrient acquisition (`fit_no_glucose`, `fit_no_glutamine`, 7,393 genes); low-nutrient (serum and lipid restriction, 7,395); survival in liver, lung, spleen, brain, heart, peritoneum (`fit_invivo_*`, 7,460 each); dissemination (the cross-tissue pattern of the same set); favourable immune response **in mice** (`fit_ifng` vs `fit_naive_bmdm`, 7,402); cell division (in vitro growth, 7,325) | ~7,400 genes each | **Carry, do not predict.** The work is turning a continuous score into a calibrated answer with an uncertainty — which needs each screen's published FDR, the thing found missing for the Wang screen. Prediction is needed only for the ~700 genes a screen missed |
| **B · Measured for a subset — a well-posed prediction problem** | egress, F-actin, apicoplast, replication (`screen_*`, 318 genes); localisation (716 markers); function class (242); phospho/acetyl/methyl/ubiquitin sites (717) | 242–717 | Standard scoring, and **far better posed than host-facing was with 65** |
| **C · Few or no labels** | a favourable immune response **in humans** (no human screen exists — identified 2026-11-04); autophagy; protein scavenging; sporozoite-to-tachyzoite switching | < 100, sometimes 0 | Curation (the host-facing review is the template) or a new experiment. Say so rather than predicting into a void |
| **D · Not boolean** | most likely binding partners; what a protein is expressed with; what it is expressed in response to | crosslinks 1,630 genes; co-expression 8,052; 266 expression columns | **Ranked lists with a probability per item**, not a single cell |

## Architecture
Five parts. Only the second is per-question work.

1. **A question registry.** One row per question: `id`, the question in Einar's own words, `answer_type`
   (boolean / category / ranked_list / number), `label_source`, the tier of those labels, `evidence_allowed`
   (families), `evidence_forbidden` (`depends_on`, so a question is never answered from its own label source),
   `class` (A–D above) and `status`. This generalises the `Trait` registry in `parascape/information_map.py`,
   which already carries tier, family, normalization, role and depends_on.
2. **A label loader per question** — the only bespoke code. For class A it reads a screen column and a
   threshold; for B a curated table; for D a pair or neighbour table.
3. **One scoring harness, written once.** Orthogroup-grouped folds, the stated baselines, isotonic calibration,
   the Clopper–Pearson lower bound, the fill rule, the no-knowledge split, permutation importance by family.
   This has now been written three times (host-facing, function, partners); factor it out rather than a fourth.
4. **Two numbers per cell, not one.** `probability` = calibrated P(answer), and `uncertainty` = the width of its
   interval (a credible interval for a boolean, the spread of the fold predictions, the standard error for a
   measured score). A confident 0.5 and an ignorant 0.5 are different cells and must look different.
5. **Question generation, not question writing.** 100–1000 questions are not typed by hand; they are derived
   from the label sources: 26 hyperLOPIT niches, 15 function classes, ~45 fitness and condition columns, every
   GO term with at least 10 experimental genes, 266 expression columns as "responds to X", 313 domain families
   as "has domain X" (measured, not predicted). That reaches several hundred immediately, and each arrives with
   its label source already attached.

## The honest caveat that must ship with the matrix
**A thousand questions is not a thousand independent pieces of knowledge.** Many will be variations of a few
underlying axes: the seven in-vivo tissue columns are largely one dissemination signal, and the serum, lipid and
density screens share a growth axis. The matrix must therefore publish a **redundancy report** — questions
clustered by the correlation of their *labels*, stating how many independent axes the set really contains.
Without it the result is an impressive-looking matrix carrying five facts. The same rule as everywhere else
applies: report each question against the hardest honest comparison, and "cannot be answered yet" is a result.

## Order of work
1. The question registry and the refactored scoring harness (the architecture itself).
2. Class A done properly: every genome-wide screen carried as a calibrated answer with uncertainty, which alone
   answers a dozen of his questions for ~90 % of genes. Needs the published per-gene statistics for each screen.
3. Class B scored through the harness: egress and the other imaging phenotypes, localisation niches, function.
4. Class D as ranked lists: partners, co-expression neighbours, expression response.
5. Class C listed as gaps with what each would need, so the missing experiments are visible rather than implied.

## The plan: choosing the questions, and choosing a strategy for each (Einar, 2026-11-05)
"What questions do we ask the data to answer, and what strategy do we use to answer each systematically?"

### Part 1 — How a question earns its place
A question is admitted only if it passes three tests, each recorded in the registry so the catalogue can be
audited rather than argued about:
1. **Answerable.** A truth set exists, or can be made by curation. Without one the question is a wish, and it is
   listed as a gap with what it would need (the host-facing run proved the cost of discovering this late).
2. **Non-redundant.** It is not a restatement of a question already present. Tested by the correlation of
   *labels*, not of predictions: the seven in-vivo tissue columns are largely one dissemination axis, and the
   serum, lipid and density screens share a growth axis. The catalogue reports **effective independent axes**,
   not question count.
3. **Useful.** The answer would change a decision -- which experiment to run, how to read a screen hit. "Has
   domain IPR000719" fails this test as a *question*: it is a feature, and it belongs on the evidence side.

**A distinction that governs the whole catalogue: "important for X" is not "involved in X".** A gene required
for growth without glucose may be a glycolytic enzyme (involved) or simply essential everywhere (not
specifically involved). So every condition question is posed as **conditional** importance -- the interaction,
not the absolute score -- and computed from each screen's own published per-gene statistics rather than by
subtracting two scores, which the mouse-immune attempt showed amplifies noise.

**The catalogue, built by enumeration rather than by hand** (counts are what exists today):

| Source | Questions it yields | Labelled genes each |
|---|---|---|
| Fitness and condition screens (nutrient, lipid, serum, density, glucose, glutamine) | ~16 conditional-importance questions | ~7,395 |
| In vivo tissue screens | 7 tissue questions + 1 dissemination pattern | ~7,460 |
| Mouse macrophage, naive and IFN-gamma | 2 immune questions | ~7,402 |
| Imaging screen phenotypes | egress, F-actin, apicoplast, replication | 318 |
| hyperLOPIT niches | 26 localisation questions | 716 markers |
| Function classes | 15 | 242 |
| Experimental GO terms with >= 10 genes | dozens | >= 10 each |
| PTM site tables | is phosphorylated / acetylated / methylated / ubiquitinated | 717 |
| Cell-cycle periodicity | is cell-cycle regulated, and peak phase | ~5,377 periodic |
| Expression columns | 266 "responds to X" questions | 8,140 |
| Pairs | most likely partners; expressed with | 1,630 / 8,052 |

That is **over 400 questions with labels already in hand**, which is how 100-1000 is reached without inventing
any. Einar's list maps onto it directly, except for four genuine gaps that must be declared, not fudged: a
favourable immune response **in humans** (no human screen exists), **autophagy** involvement, **protein
scavenging**, and **sporozoite-to-tachyzoite** switching.

### Part 2 — The strategy is chosen by the question, not fixed
A ladder. Each rung is used only when the one below is insufficient, and every rung reports the same numbers:
grouped folds, the stated baselines, isotonic calibration, the Clopper-Pearson lower bound, coverage at the
90 % bar, the no-knowledge split, and permutation importance by family.

| Rung | When it applies | What it does |
|---|---|---|
| **0 · Carry the measurement** | class A, ~7,400 genes measured | No model at all. Convert the measured score to an answer plus an uncertainty from the source's own statistics. Predict only the few hundred genes the screen missed. This alone answers a dozen of Einar's questions for ~90 % of genes |
| **1 · Baselines** | every question, always | Majority, prevalence, best single feature, ortholog copy. Nothing is reported except as improvement over the best of these |
| **2 · Per-family model, then stacking** | class B | One model per evidence family, combined by a meta-learner on out-of-fold predictions. Makes each family's contribution measurable instead of hiding 86 of them in one model |
| **3 · Guilt by association and network propagation** | questions about complexes, pathways, co-membership | Label propagation over co-expression, crosslinks and co-fitness, with the cofitness rule: high correlation **and** top rank **and** support from several studies |
| **4 · Matrix completion** | the matrix as a whole | Low-rank factorisation predicts a cell from the pattern of co-occurrence across all questions. The right frame when no single feature carries signal, which is what every run here has shown |
| **5 · Multi-task joint model** | last | Shared representation across correlated questions; replaces per-question models only where it beats them on the same folds |

Three rules bind every rung: positive-unlabelled weighting, because "not recorded" is not "not true";
inverse-propensity weighting, because a well-studied gene is more likely to be measured at all; and a question
is never answered using the family its own labels came from.

### Part 3 — Phases, each with a gate
1. **The registry and the harness.** Question registry (id, text, answer type, label source, tier, allowed and
   forbidden families, class, status) and the scoring harness factored out of the three runs already written.
   *Gate: one question reproduces its published scoreboard through the new harness.*
2. **Rung 0 across class A.** Every genome-wide screen carried as a calibrated answer with uncertainty.
   *Gate: a dozen of Einar's questions answered for ~90 % of genes, with each screen's own statistics cited.*
   Blocker to clear first: several screens are shipped as effect sizes only; the published per-gene FDR tables
   must be acquired.
3. **The catalogue and the redundancy report.** All ~400 questions enumerated, with effective independent axes.
   *Gate: Einar sees the list and prunes it.*
4. **Rungs 1-3 on class B**, in order of label count: localisation niches (716), function (242), imaging
   phenotypes (318), PTM (717), host-facing (98). *Gate: the per-question scoreboard, and the count of cells
   fillable at 90 %.*
5. **Rungs 4-5**, only if 1-3 leave cells unfillable. *Gate: beats the per-question models on the same folds.*
6. **Freeze and date the filled matrix** for the prospective test against papers published afterwards.

### What makes the whole thing honest
Four numbers published per release, never one: **questions with enough labels to score**, **effective
independent axes**, **cells fillable at 90 %**, and **prospective accuracy** once the first frozen matrix is old
enough to test. A matrix of a thousand questions carrying five facts is the failure mode to avoid, and only the
second number reveals it.

## Validation, automated: nested cross-validation AND a frozen test partition (Einar, 2026-11-05)
"Whenever predictions are made there should always be a holdout set with known answers that the model predicts,
and the quality of the strategy measured that way. Did you bake in a mechanism like this, or cross validation?
Preferably both."

### Honest answer: cross-validation yes, a frozen holdout no — and the primitive for it already existed
All three scoring runs (host-facing v1 and v2, function, partners) used **5-fold cross-validation grouped by
orthogroup**, with calibration cross-fitted. Out-of-fold predictions on known labels are a held-out measurement,
so the estimates are not circular. **But there was no partition held back from everything.** The same folds that
produced the reported number also chose the feature set: 14 feature-set variants were tried for host-facing and
10 for function. Picking the best by cross-validation and then quoting its cross-validated score is optimistic,
and the honest numbers (0.727, 0.617, 0.891) are therefore **upper bounds on what a frozen test would give**.
Worse, `parascape/splits.py` already provides exactly the right primitive and the runs did not use it:
`ROLES = ("train", "tune", "calibration", "test")`, `make_split` splitting whole protected groups, and
`guard_fit`, which **raises** if a stage touches a partition it may not see (`FIT_ACCESS` lets model fitting and
feature selection see only `train`, setting selection `train`+`tune`, calibration only `calibration`). That
guard is the mechanism; it was simply not wired in.

### What the harness must enforce, for every question, automatically
1. **Four partitions by `splits.make_split`**, grouped by orthogroup (whole-source where a question's labels come
   from one study; cold-node for pair questions), seed recorded. `train` fits, `tune` selects settings,
   `calibration` fits the isotonic curve, **`test` is touched once**.
2. **Nested cross-validation inside `train` + `tune`** for every model and feature choice. This is where the
   searching happens and where the comparison between strategies is made.
3. **`guard_fit` called at each stage**, so a leak raises instead of being caught in review. A question cannot
   reach the fill step without its guard calls recorded.
4. **A holdout budget, recorded per question.** Each strategy version may consume the `test` partition **once**;
   the use is logged with the strategy's version hash and the result. A second look at the same partition with a
   changed strategy makes it another tuning set, so it is refused and the question is marked as needing a fresh
   partition. Without this the frozen test silently becomes a third cross-validation.
5. **Both numbers published, never one**: the nested cross-validated estimate (how the strategy compares, on
   many folds) and the frozen-test result (what it actually achieves, once). **Filling a cell requires the
   frozen-test lower bound to clear the bar**, not the cross-validated one.
6. **Small-label questions are exempt and must say so.** A frozen partition costs labels that a 98-gene question
   cannot spare. Below a recorded threshold the question keeps nested cross-validation only, is flagged
   `no_frozen_test`, and **may not fill cells** -- it can rank and prioritise, which is what those questions are
   for. The threshold is a property of the question in the registry, not a judgement made per run.
7. **The prospective test is the final holdout**: the filled matrix is frozen and dated, then scored against
   papers published afterwards. No partition can be reused that way, and it is the only test immune to every
   form of selection on today's data.

### Consequence for what has already been reported
The three published scoreboards stay as they are, relabelled as **cross-validated estimates without a frozen
test**, which is what they are. Nothing was filled from them, so no claim has to be withdrawn -- but the
architecture work must not repeat the omission, and re-running the three questions through the harness with a
frozen partition is part of phase 1's gate.

## Phase 1 built (2026-11-05): the registry, the harness, and what the harness immediately caught
`parascape/question_registry.py`, `parascape/harness.py`, `tests/test_question_registry.py`,
`tests/test_harness.py`, `notebooks/harness_demo_2026_11_05.ipynb`,
`results/harness_demo_2026_11_05/`, `results/harness_holdout_log.jsonl`.

* **Named apart from `parascape/questions.py`**, which is instruction 45's catalogue of what the
  *table* can be asked. The instruction 73 registry is `question_registry`; the two are different
  objects and the first attempt at this work overwrote the older module, which the suite caught.
* **Registry fields**: `id`, `text`, `answer_type`, `question_class` (A-D), `label_source`,
  `label_family`, `label_tier`, `evidence_allowed`, `evidence_forbidden`,
  `min_labels_for_frozen_test` (default 400, the arithmetic of a 0.90 Clopper-Pearson bound through
  a 15 % partition), `status`, the three admission tests (`answerable`, `non_redundant`, `useful`)
  and `redundant_with`. One structural guarantee: **registering a question whose `label_family` is
  not in `evidence_forbidden` raises**. Eight questions registered, at least one per class,
  including the two declared gaps and `has_domain_ipr000719` kept as the worked example of a
  question rejected for failing the usefulness test.
* **Harness**: four `splits.make_split` partitions grouped by orthogroup (`cold_node_pairs` for
  ranked lists), `StageGuard` calling `split.guard_fit` at every stage and recording each call,
  nested cross-validation inside train + tune where the strategy and the fill threshold are chosen,
  isotonic calibration on the `calibration` partition only through `claims.CertaintyModel`, one
  scored look at `test`, and both numbers returned with Clopper-Pearson lower bounds, coverage at
  the bar, calibration error, baselines and permutation importance by family. `harness.may_fill`
  reads the **frozen-test** bound alone.
* **Two API frictions worth recording.** `FIT_ACCESS["model_fit"]` permits `train` only, so nested
  cross-validation inside train + tune fits under the `refit` stage, which is the permission
  `splits` already defines for exactly that; and a cold-node pair question is guarded with its
  endpoint **nodes**, because that is the unit `splits` partitions.
* **The budget hashes the evidence as well as the strategy.** Widening `evidence_forbidden` changes
  what a model sees while leaving the model, its settings and its version identical, so the strategy
  hash alone called a changed run a repeat. The log entry carries `strategy_hash` **and**
  `evidence_hash`, and a change in either refuses the partition.
* **Demonstrated on `required_without_glucose`** (class A, 7,393 labels, 1,765 positives, 371
  features in 77 families): cross-validated lower bound **0.7722**, frozen-test lower bound
  **0.5938** (8 calls of 1,082 test genes), ROC-AUC 0.83 on both. Nothing may be filled, and nothing
  about glucose biology is claimed: rung 0 carries the measurement for these genes and the
  conditional form still needs the screen's published per-gene statistics.
  `screen_egress_phenotype` (318 labels) is flagged `no_frozen_test` and cannot fill, as designed.
* **What the harness caught on its first use, which is the point of it.** Permutation importance and
  the best-single-feature baseline found three further pooled CRISPR fitness screens whose family
  names do not contain "crispr" -- `gra17_synthlethal`, `invivo_platform`, `gra12` (and
  `gse253884_5`) -- each answering the question with the same kind of measurement. The forbidden
  list is now enumerated by a rule (every family contributing a column named `crispr*` or `fit_*`,
  plus the curated enteric set), and the correction cost two frozen partitions: seeds 0 and 1 are
  spent and logged, and the reported run is seed 2. That is the budget working, and the honest price
  of a corrected definition.
* **Still open from phase 1's gate:** re-running host-facing, function and partners through the
  harness with a frozen partition. Host-facing (98 reviewed genes) and function (242) are below the
  frozen-test threshold as registered, so phase 1's gate for them is a `no_frozen_test` flag rather
  than a number -- which is itself the finding, and the reason the three published scoreboards stay
  labelled as cross-validated estimates.

## Information categories: the mechanism that stops circular inference (Einar, 2026-11-05)
"If the user builds a UMAP from a subset of data and does clustering, then measures cluster-to-held-out-label
precision and recall, it is important that **each data column in ParaScape is annotated with an information
category**, otherwise the UMAP inference will often be circular."

This is the generalisation of the per-question `evidence_forbidden` list, and it is needed because the
hand-written version already failed once: on the harness's first real run, the family importance exposed three
further pooled CRISPR fitness screens whose family names contain no "crispr" (`gra17_synthlethal` was top of the
table), each answering the question with the same kind of measurement as its labels. A per-question list caught
that only because someone read the importance table. **A per-column category catches it automatically.**

### The rule
**Every column in the map carries an information category: what kind of knowledge it embodies.** Inferring a
label of category X **excludes every column of category X by default**, whatever its name, family or study.
Including one is possible but requires an explicit reason that is recorded in the run and shown on the
scorecard, so a circular run is visible rather than impossible.

Categories, which are about the *kind of knowledge*, not the assay:
`place` · `fitness` · `transcript` · `translation` · `protein_abundance` · `modification` · `interaction` ·
`sequence` · `structure` · `domain` · `conservation` · `chromatin` · `phenotype` ·
`function` (what the protein is held to DO: GO terms, EC numbers, functional class) · `host_response` ·
`attention` (bookkeeping: publication counts, study effort — a predictor of nothing biological) ·
`identity` (names and accessions — never a feature).

A column may carry more than one category where it genuinely does; the exclusion then applies if **any** of its
categories matches the label's.

### Why this is load-bearing for the UMAP strategies specifically
An embedding built from "all the data" and then evaluated against held-out localisation labels is circular if
any input column carries `place` -- and 8 of the map's columns do, plus every hyperLOPIT-derived value. Without
the annotation the user cannot tell, because the circularity is in the *input selection*, not in the model. With
it, the data picker shows which categories are in and which were removed and why, and the scorecard records it.

### Where it is enforced
1. **At ingest.** A new column cannot enter the map without an information category. This is the gate that keeps
   the annotation complete as data keeps arriving, rather than a catalogue that rots.
2. **In the data picker** (`Label → Infer`): the user chooses *some or all* data by category, family and tier,
   and the excluded-by-circularity set is shown, named, and counted rather than silently dropped.
3. **In the harness**: the exclusion is computed from categories, not from a hand-written family list, and the
   run record stores both the requested selection and the effective one.
4. **On the scorecard**: every result states which categories it used, which were excluded automatically, and
   any that were force-included with the stated reason.

## Saved inference runs (Einar, 2026-11-05)
"It should be easy to save inference runs and they should hold all the metadata necessary to know which data
they used, from which data, and so on." A run record is therefore reproducible or it is not saved, and holds:
the question and its version; the strategy, its version and settings; the **requested and effective column
selection with categories**; the data release version and the sha256 of each table read; the split identity,
seed and partition sizes; every guard stage recorded; the cross-validated and frozen-test numbers with their
bounds; the calibration curve; the holdout-log entry; which claims were kept and by whom; and the date.
Re-running a saved record on the same release must reproduce it exactly; on a newer release it must say what
changed.

## Information categories built (2026-11-05): every column annotated, and the exclusion derived from it
`parascape/information_categories.py`, `tests/test_information_categories.py`,
`scripts/build_information_categories_review.py`, `results/information_categories_2026_11_05/`.

* **1,056 columns over five shipped tables, none unannotated.** `nodes.parquet` (443), the
  sequence (22), domain (313) and expression (266) trait tables and the map's own `values.parquet`
  (20). `information_categories.unannotated()` is empty and the gate test keeps it empty; a column
  that genuinely cannot be classed is to be left out of the mapping so the gate fails loudly rather
  than a wrong category passing quietly.
* **Derived from provenance, not from names.** `datasets.provenance(column, organism)` gives the
  level and kind, and that pair is the primary rule; a dataset whose level is true but not decisive
  gets a record of its own justified from its registry entry; a single column that differs from its
  dataset gets an exact override. 1,003 of the 1,064 rows rest on recorded provenance (the dataset
  registry, or the per-column registry shipped beside a trait release) and 60 on a judgement, each
  listed with its reason in `judgements.tsv`.
* **The traps, verified by enumerating the registry rather than by guessing.** All 16
  `DNA/CRISPR_screen` datasets are `fitness` -- `gra17_synthlethal`, `invivo_platform` and `gra12`
  among them -- except the host-transcription effector screen, which is `host_response` as the
  instruction says. `fit_hyperlopit_unassigned_invivo_lib1/2` are **`fitness`, not `place`**: the
  name records which library was screened and the measurement is guide depletion, so a name-based
  rule gets that one exactly backwards. hyperLOPIT, TAGM and the four `ortholopit_*` columns that no
  dataset entry claims are `place`; `pvm_proximity_positive` is `place` and `interaction`; every
  `te*`, `rpf*` and `expr_*_RIBOseq` column is `translation` while the RNA arm of the same deposits
  stays `transcript`; `gene_id`, `product` and `orthogroup` are `identity` and are never features.
* **The harness exclusion is now computed from categories.** `harness.evidence_selection` returns
  the requested and the effective selection with the reason for every removal;
  `run_question(force_include={'<column or category>': '<reason>'})` is the recorded escape hatch,
  and it is hashed into the evidence identity, so including a circular column on purpose spends the
  frozen partition exactly as widening the forbidden list does. The hand-written `evidence_forbidden`
  list is **kept as an additional, recorded override**: `pubmed` is `attention`, so no `place`
  category would hold it out of the host-facing question, where it is nearly the conclusion. It is no
  longer the mechanism. Each registered question now also carries `label_category`, and a question
  that leaves it empty has it read off the column named in `label_source`.
* **The gap in the vocabulary, flagged and then closed.** The vocabulary listed no `function`
  category, so functional annotation (`function_class`, the GO BP and MF columns, `ec_number`) was
  annotated `phenotype` -- right for a function question, slightly wide for an imaging-phenotype one.
  Einar took the call on 2026-11-05: `function` is now the seventeenth category and holds those
  columns. See "A seventeenth category" below for why `phenotype` was the wrong home.
* **Still open from this piece:** the data picker (`Label -> Infer`) does not yet show the
  excluded-by-circularity set, and the ingest gate that refuses a new column with no category is the
  test rather than the loader.


## A seventeenth category: `function` (Einar, 2026-11-05)

The vocabulary shipped with sixteen categories and no `function`, so GO terms, EC numbers and coarse
functional class were filed under `phenotype` -- recorded at the time as a judgement, with the reason that
`phenotype` ("what the protein is held to do") was the nearest available kind of knowledge.

It was the wrong home, and the cost is specific. `phenotype` also holds what is *observed when the gene is
perturbed*: imaging phenotypes, drug sensitivity shifts, differentiation readouts. Those are measurements. A
GO molecular-function term is a conclusion somebody already drew about what the protein does -- which is the
thing a function question is trying to predict. Pooling the two means a function question either holds out
real perturbation phenotypes it could legitimately learn from, or an imaging-phenotype question holds out
annotation it never needed to. One of those is always wrong, and which one depends on the question.

`function` is therefore its own category: **what the protein is held to do**, as distinct from `phenotype`
(what happens when you break it) and `place` (where it is). It takes `function_class`, `function_class_n`,
`function_go_bp`, `function_go_bp_n`, `function_go_mf`, `function_go_mf_n`, `ec_number`, `has_ec` and
`toxodb_ec_numbers`.

`phenotype` keeps its definition narrowed to match: "what is observed when the gene is perturbed --
imaging, drug sensitivity, differentiation".

Consequence to expect: the function question now holds these columns out by category rather than by a
hand-written family list, so its evidence hash changes and the holdout budget will demand a fresh seed.
That is the budget working as designed.

## The host-transcription effector screen also carries `place` (Einar, 2026-11-05)

`hosttx_effectors` was `host_response` alone, which is literally right -- the readout is the host
transcriptional response, not parasite depletion. But a screen that finds which parasite effectors change
host transcription is very nearly a list of proteins that reach the host, so for a `place` question it sits
close to the conclusion. It was being held out of host-facing only by the hand-written override, which is
exactly the mechanism information categories were built to replace.

It now carries `place` as well, so the category rule excludes it from any localisation question whatever a
family list says. This is the conservative direction: the cost of the extra category is one feature a
`place` question may not use, and the cost of omitting it is a host-facing model reading its own answer.

## Does the category rule reach the map-and-cluster questions? Measured 2026-11-06

The harness holds out a label's own information category; the Infer module's map-and-cluster questions
(`questions.py` recipes) use the older association-based closure (`recipes.close`, `excluded_detail`). Measured
over the 21 shipped questions: **for 18 the closure already removes every input of the label's category.**
Three would change under the category rule, all with a `transcript` label:

| question (holdout) | inputs of the label category / all inputs | examples |
|---|---|---|
| persist in the mouse brain at 14 d (`brain168465_14d_lfc`) | **156 / 174** | `expr_tachy`, `rna108740_Tachyzoites_*` |
| what the BFD2 switch holds (`bfd2_rip_log2_ip_over_input`) | 1 / 27 | `novel_transcript_models` |
| held silent in tachyzoites (`ees_vs_tachyzoite_log2`) | 4 / 51 | `n_m6a_peaks`, mRNA decay after actinomycin |

**Decision for Einar (not applied):** is predicting one transcriptome contrast from other transcriptome columns
circular? Under the category rule the brain question loses nearly all its evidence and would be refused; under
the closure it stands. Until decided, the shipped answers are unchanged, and this table is the record of the gap.

## Decided and applied (Einar, 2026-11-07): strict

Einar chose the strict rule: predicting a label from columns of the same information category is circular. It
is now in `recipes.close` (`recipes.same_category`): every input sharing a category with the **label** (identity
aside) is held out, with the reason "same information category as <label> (<category>)". Two deliberate limits:

* **Inputs only, never the control check.** A control is meant to measure the same kind of thing
  independently; applied to the control-independence check the rule refused 9 questions (one fitness screen
  "not independent" of another). It is applied after that check.
* **The label's category only.** Applied for the control's category too, it emptied the expression inputs of
  two fitness questions whose control is an expression column -- beyond what was decided.

Effect on the 21 map questions: 19 ship (was 21). Refused: the mouse-brain question (`brain168465_14d_lfc`,
156 of 174 inputs were expression, as measured) and "held silent in tachyzoites" (`ees_vs_tachyzoite_log2`,
whose 4 transcript-category columns were a whole input block). The BFD2 RIP question loses one column and
stands. Category lookups are cached per session (a closure of all questions: 92 s -> 44 s). Test in
`tests/test_recipes.py`; 84 recipe and question tests pass.
