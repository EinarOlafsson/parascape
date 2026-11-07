# 70 · One inference scoreboard, one deployment rule, and a question-first navigator

**Status: open (requested 2026-10-31).**

**Paused 2026-11-01** while Einar's dataset goal ("finish all downloads and implement and normalize all
datasets") takes priority. Direct-tier evidence acquired so far: QuickGO experimental cellular-component
annotations for taxon 5811 (78 annotations, 37 proteins, 28 matrix genes; snapshot
`toxoplasma_projects/datasets/sources/quickgo_tgondii_cc_experimental_2026_10_31/`). Einar chose: score
first on the 718 curation-grade markers, report the ~37 direct genes as a check, and build a larger
literature-mined direct set for his review.

## The request (Einar, 2026-10-31)
* The core goal has always been the same: a comprehensive knowledge matrix of everything known about
  Toxoplasma genes, used to fill in the blanks -- function, localization, expression profiles, binding
  partners.
* Method: take known information, hold some out, apply a strategy and measure how well it recovers the
  held-out labels. If a strategy works, cross-validate it; if it holds, deploy it for inference.
* "We have a ton of strategies in ParaScape but I have not really figured out how to gauge the tests,
  and I don't know how good they are." "As it is now, I get lost right away."
* Borrow the spaCR navigator and protocols: a search field for modules and other things, and predefined
  questions that guide the user to what ParaScape can do.
* Decisions: build the scoreboard and the question tree (approved); **localization is the first target**.

## Why the current tests cannot answer "how good is it"
* The 39 strategy self-tests (strategies.py `label_transfer_test`, `judge`) hide 25 % of a label once per
  seed and compare against the same procedure on shuffled labels. PASS means "better than chance", so a
  barely useful strategy and a strongly useful one both pass.
* The localization truth used so far (`compartment`, `lopit_unified`) is a classifier output
  (`truth_grade: prediction`, `benchmark_admitted: false` in results/label_knn_pilot_2026_10_08_v5).
  Agreement with it measures agreement with another model, not with biology.
* Strategies are not compared on one frozen split against the same simple baselines, there is no
  source/study hold-out, and there is no explicit rule for deployment.

## Plan
1. **Graded truth per target.** Reuse `results/spatial_truth_source_review_2026_10_08` (GT-SPATIAL-01),
   which already audited the hyperLOPIT tables; do not repeat it. Localization tiers:
   * **Direct experiment** (new acquisition): GO cellular-component annotations for *T. gondii* with
     experimental evidence codes (EXP/IDA and relatives) from UniProt-GOA/QuickGO, each with its PMID,
     plus the 9 Figure-3 IFA outcomes of Barylyuk 2020 that are absent from both marker fields (the six
     unresolved attempts stay unknown, never negative).
   * **Curation**: the 718 final hyperLOPIT markers (mmc5.xls `markers`, 26 niches; 716 exact IDs, 2
     unmapped retained). The paper states they combine 656 prior markers -- selected from previous
     knowledge, inferred function and profile consistency -- with all 62 new microscopy outcomes, so marker
     status is *not* independent experimental truth. Snapshot re-downloaded 2026-10-31 to
     `toxoplasma_projects/datasets/sources/barylyuk2020_mmc5/mmc5.xls` (sha256
     `a58a8778714c42e344099de9429daa3be4e9c70021b0aeab3de65578c165b8f2`, MD5 `92ef73be6dd910b33a4d7eb8b20b1377`
     equals the published PMC MD5); 686 of 716 equal the installed `compartment` value.
   * **Prediction** (reported separately, never as accuracy): confident TAGM calls for non-marker proteins.
   Report both the 26 niches and a coarse level (merging numbered sub-niches such as rhoptries 1/2); figure
   and GO labels are coarse and are not translated to refined niches automatically.
2. **Whole-source hold-out.** Every column derived from the hyperLOPIT experiment (the `lopit_tgon`
   provenance family, including MAP/MCMC posteriors and `lopit_methods_agree`) and every transfer built
   from it (orthoLOPIT) is removed from the features. The question is whether the other evidence places
   a protein. TAGM itself cannot be scored on its own training markers; it is listed, not ranked.
3. **One frozen split for every strategy.** Five folds grouped by orthogroup, fixed seed and pinned
   cohort digest (ground_truth.cohort_digest); tuning only on inner folds (splits.guard_fit).
4. **Simple baselines in the same table.** Training majority, prevalence-weighted random, kNN on a single
   evidence block, and ortholog-label copy. Every strategy is reported as improvement over the *best*
   baseline, not over shuffled labels.
5. **Deployment rule (initial, to be confirmed on the first results).** Deploy for a target when the
   strategy beats the best baseline in at least 4 of 5 folds *and* its confident calls reach at least
   80 % precision on direct-evidence genes, with the confidence threshold chosen on inner folds.
   Otherwise the strategy is "exploratory" for that target. Deployed strategies fill gaps with calibrated
   confidence; nothing else is presented as an inference.
6. **Scoreboard view.** Strategy x target table coloured deploy / exploratory / fails, each cell opening
   the folds, baselines, truth grade and excluded sources.
7. **Question-first navigator.** Promote the existing ranked search (help_search.py, Ctrl+Shift+H) to the
   start page and add a protocol tree in the spaCR style (module_workflows.json `protocol_tree`): "Where
   is my gene?", "What does it do?", "Is it essential?", "Which genes behave like mine?", "Does a
   strategy work for label X?", "Fill in the blanks for my gene list". Leaves open a module with the
   context filled in; strategies sit behind questions rather than in the first menu.
8. **Expression evidence as scored strategies.** The Codex expression proposals enter as candidates:
   coexpression neighbours/modules -> localization, function and partner targets; the sample/context
   registry and replication browser precede any RNA–protein comparison; transformers only after PCA,
   regularized regression and a small autoencoder are beaten with whole-study hold-outs.

## Rules that apply
* Unlabeled is not negative; missing measurements are not zeros; repeated runs are not replicates.
* A score against prediction-grade truth is never presented as biological accuracy.
* Leakage: register every new evidence block with its provenance family before scoring.

## Progress
| Action | Item | Status |
|---|---|---|
| 70.01 | Localization truth tiers: curation markers (done in GT-SPATIAL-01), direct GO-experimental + 9 IFA acquired and pinned | 30% |
| 70.02 | Whole-LOPIT feature exclusion verified by leakage audit | 0% |
| 70.03 | Frozen grouped 5-fold split + baselines | 0% |
| 70.04 | All label strategies run on the frozen localization split | 0% |
| 70.05 | Deployment rule applied, scoreboard data written | 0% |
| 70.06 | Scoreboard view in the app | 0% |
| 70.07 | Question-tree start page with ranked search | 0% |
| 70.08 | Coexpression strategies scored on the localization board | 0% |

## Where each action stands (2026-11-06)

Most of this plan was since built generically -- for every question rather than for localization alone --
by the harness (instruction 73) and the navigator (84). Stated here so the table above is not read as 0 %.

| Action | Delivered by | What is still missing for localization |
|---|---|---|
| 70.02 whole-LOPIT exclusion | `harness.evidence_selection`: inferring a label of category `place` holds out every `place` column, whatever its family (`information_categories`); `lopit_tgon` is also named in `evidence_forbidden` where it applies | nothing for the mechanism |
| 70.03 frozen grouped split + baselines | `harness.partition_question` (orthogroup-grouped train / tune / calibration / test, seed recorded, budget-logged) and `harness._baselines` (training majority, prevalence-weighted random, best single feature, ortholog copy; the single-block kNN baseline proposed here is not among them) | a localization question in the registry |
| 70.04 all strategies on the split | the borrowed strategies (77) and boosted trees run through it | **the harness scores boolean labels only** (`labels.astype(bool)` in `run_question`); a 26-niche label needs a categorical path, or one-vs-rest per coarse niche |
| 70.05 deployment rule | `harness.may_fill`: the frozen-test Clopper-Pearson lower bound must clear 0.90 -- stricter than the 80 % proposed here, and decided on a partition touched once | -- |
| 70.06 scoreboard view | the Infer module "How far to trust it" and the strategy cards; the track record (62) per gene, class and label | a strategy x target table for the harness results |
| 70.07 question-first start page | the Navigator (84): search on Return, module cards, downstream network | -- |
| 70.08 coexpression strategies | not started | -- |

**The next concrete step** is a categorical path in the harness (or a one-vs-rest decomposition over the
coarse niches) so that `compartment` from the 718 curation-grade markers can be registered and scored.

## 70.04 first results (2026-11-06): localization ranks well, and cannot fill yet

`scripts/build_localization_questions_notebook.py` -> `notebooks/localization_questions_2026_11_06.ipynb`,
`results/localization_questions_2026_11_06/`. Localization asked one coarse niche at a time (13 niches with
at least 25 of the 718 curation-grade hyperLOPIT markers; every other marker a curated "no"), every `place`
column held out, boosted trees, nested cross-validation, frozen partition once (seed 0; the sensitivity run
seed 1). Built in the notebook, **not registered** in the application.

* **Frozen-test AUC 0.85-1.0 for all 13 niches** (median 0.964): ribosome 1.0, proteasome 1.0, apicoplast
  0.997, dense granules 0.993, rhoptries 0.99, micronemes 0.985, nucleus 0.964 ... IMC 0.846. The synchronised
  cell-cycle series (`gse19092`) is the top evidence family for most niches.
* **Nothing may fill**: a frozen partition holds 104 genes with 4-13 positives, so even perfect calls give a
  Clopper-Pearson lower bound well under 0.90. More markers is what lifts it, not a better model.
* **Not placed by who chose to measure.** `gra12` (a CRISPR screen of 236 secretory candidates) topped dense
  granules and rhoptries, and being measured at all encodes that choice. Rerun with only families measured on
  at least half the genes: median cross-validated AUC change 0.0 (dense granules 0.957 -> 0.937).
* **Caveat that remains**: the markers were partly chosen from prior knowledge, and for dense granules the
  genome-wide top family becomes InterPro, whose domain names can carry the location.

Next for localization: more independent positives (the GO experimental set of 70.01 and the 9 IFA outcomes) to
widen the frozen partitions, and Einar's decision on registering these questions in the application.

**70.08 (2026-11-06):** co-expression scored on the localization board -- see instruction 71, product 1 through
the harness: median AUC 0.897 from co-expression alone, matching all evidence for the division organelles.

## Registered, rank only (Einar, 2026-11-07)

The 13 localization niche questions are in the registry (`question_registry.LOCALIZATION_NICHES`,
`localization_<niche>`), status `scored`, each note carrying its frozen-test AUC, cross-validated AUC and the
co-expression-alone AUC, and saying it may rank and may not fill. Labels ship as
`parascape/data/lopit_markers_coarse.tsv` (718 markers from mmc5, sha256 `a58a8778...`, with the coarse niche).
They appear in the Navigator (rendered: a search for "rhoptries" finds "Localization Rhoptries", "Is this protein
in the rhoptries?"). 229 node, home, module, menu and search tests pass; one search test was widened because a
product search for "dense granule" now also, correctly, finds the dense-granule question.
