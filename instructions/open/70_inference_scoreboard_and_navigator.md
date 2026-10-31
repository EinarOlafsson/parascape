# 70 · One inference scoreboard, one deployment rule, and a question-first navigator

**Status: open (requested 2026-10-31).**

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
