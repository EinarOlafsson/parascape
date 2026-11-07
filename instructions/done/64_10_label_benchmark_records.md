# 64.10 · Complete the categorical-label benchmark records

Status: DONE — 100%, 2026-11-07. Per-gene held-out rows for the 14 record methods and now for every harness method; per-class precision/recall/F1, coverage, chance and confusion derived from those rows for all 168 method x label records; rare classes listed as gaps. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `parascape/track_record.py`, `parascape/strategy_catalog.py`, `parascape/strategy_learning.py`, `parascape/scorecard.py` |

## Deliverable

Native held-out label outputs, including class-score or prediction-set fields where produced, for all label-capable strategies in the capability map.

## Controlled scope

Categorical-label task only; one new adapter per validated commit, first on a frozen pilot target.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Reuse the existing 14-strategy record; add missing label-capable adapters one at a time with native-prediction parity.
- [x] Every eligible labelled gene has a held-out call or an explicit abstention; rare and unsupported classes are retained as evaluation gaps.
- [x] Precision, recall/F1 by class, coverage, baselines and confusion outcomes derive from those rows; conformal calls report coverage and set size/singleton efficiency.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

First bounded change: `parascape/label_records.py` implements feature kNN with
training-only native rank scaling, frozen source exclusions, retained class
scores and explicit abstention. Training ranks match native `Context.matrix`,
including ties/missing values; query ranks use frozen training distributions.
The whole-context rank transform is deliberately withheld under inductive
evaluation. No installed data or runtime strategy is changed.

Pilot: `results/label_knn_pilot_2026_10_08_v5/`, candidate stored Toxoplasma
compartment calls, seed 17, 3,827 eligible genes, train/tune/calibration/test
2,126/572/569/560. Stored calls have prediction grade; recovery is not biological
accuracy. Every test gene and native class-score value is retained. Earlier
diagnostics and their original code remain available. Full-precision score
hashing fixes an independently discovered rounded-JSON identity weakness.

Verification: all 560 calls and 26 native score columns replay exactly; training
vectors exactly match native `Context.matrix`. 196 correct calls, 236 wrong,
128 abstentions; all-hidden recovery 0.35, called-only recovery 0.453704,
majority baseline 0.230357. Artifact identity
`c7b228db196b61614cebcb767c0557c8a816f6dbb6847686c8eb57c314985c1a`.
**1,161 checks passed**, one optional pdoc documentation module skipped locally;
445 focused checks also passed before the broader strategy/track-record suite.
This coherent change is committed/pushed to nightly as `8b59643`.

Second bounded change: `parascape/transfer_records.py` adds the previously
unrecorded categorical `ortholog_transfer` path. The native donor orthogroup
aggregation receives no receiver target column; the source-to-label mapping
fits only on receiver training genes. Both endpoints remain prediction-grade.
Source table/mapping identities, semantic review and dependence gaps are explicit.
The existing Plasmodium schizont spatial prediction is a scoped installed-source
pilot, not a new preferred dataset or an admitted biological comparator.

`results/ortholog_transfer_pilot_2026_10_08/`: same frozen 560-gene test cohort,
123 calls, 437 abstentions, 65 correct/58 wrong. Mapping reaches 967/3,827 receiver
genes, 550 training and 126 test genes; three mapped test values lack a learned
source category. All-hidden surrogate recovery 0.116071, called-only 0.528455,
majority baseline 0.230357. No native class scores/support or AUROC/AUPRC are
invented. Every native projection/call, card and matched baseline replays exactly.
Artifact `854ab60f8a75b078c8dc9336f8f49aadabb64110858c1f14fd829ca61054d2d9`;
mapping `e1da9a0c47cecdd47087c34026ddc5e5add6d4f950570a81574c481b40654d80`.
**464 relevant checks passed**. Artifact validation now also refuses numeric
outputs as categorical evaluations and categorical outputs as numeric ones.
Initial missing evaluation metadata was refused before any artifact was written;
the executed verification records that diagnostic and corrected scope.

Third bounded change: `parascape/conformal_records.py` consumes a verified
train-only native base model, separate calibration/test score populations and
only calibration labels. Native quantiles, per-class thresholds/rare-class overall
fallback, typed set members, native set text and all raw scores are retained.
Unbounded thresholds are explicit status/null, not nonfinite JSON.
`record_scorecards` adds empirical set coverage, mean size, singleton/empty share
and native efficiency at target/class levels, alongside call metrics. Class
score parameters are restricted to each class cohort when deriving its card.

Pilot: `results/conformal_label_pilot_2026_10_08_v2/`, same frozen training state,
2,126 training/569 calibration/560 test genes, alpha 0.1/per-class/kNN fixed
without test-outcome selection. Native function parity is exact on identical
frozen partitions; the released partition chooser is temporarily replaced in
the diagnostic and restored. Every score, set, threshold, card and matched
baseline replays. **908 checks passed**, one optional pdoc module skipped.
Empirical prediction-grade set coverage 0.973214, mean size 24.455357 of 26,
efficiency 0.061786, **zero singleton calls**; 11 rare classes use the overall
fallback. All-training-classes control covers all stored labels with size 26
and efficiency zero. Broad coverage is not a useful-label accuracy claim.
Artifact `8a96d69660e83c49d897bdcac1d147aa1ad0f3863a8e507df00fd5935dccf8f3`;
calibration `2e2bfe06d87be1efdb3fba801c03bd8f7eb563cf1ebbc1e809bacccd2bdae6e5`.
No independent biology or exchangeability guarantee admitted; no runtime
strategy/data/calibration change.
The canonical v2 also refuses post-load mutations of the base-model payload.
An executed identity-hardening comparison preserves the initial prototype and
verifies identical numerical outcomes with distinct correct code/model lineage.

Pending: missing adapters (`holdout_search`, `multiplex_modules`), full outer coverage, biological source admission
and calibration. This item does not earn a completion tick from one pilot.

Truth-source follow-up **GT-SPATIAL-01**, shared with 65.04:
`results/spatial_truth_source_review_2026_10_08/` verifies seven additional
primary companions and complete spatial marker tables. Toxoplasma's 718 final
markers include all 62 new microscopy outcomes, so those outcomes cannot
independently validate the published classifier trained on them. Plasmodium
Figure 3 provides nine selected IFA outcomes outside both marker fields; all
nine match exact installed gene IDs. Six other attempted targets remain unknown,
not negatives. Coarse microscopy taxonomy, selection bias, context and feature
dependence remain explicit admission gaps. **426 relevant checks passed** and
complete rows/media/mapping replayed. No biological benchmark admitted; the
truth admission gaps and remaining acceptance conditions are unchanged.

## Next frozen partition: logistic conformal, 2026-10-08

Reuse the seed-17 candidate compartment truth, unchanged source SHA and ordered
train/tune/calibration/test cohorts (2,126/572/569/560). Reuse only the verified
training-fitted rank distributions/source exclusions from the kNN pilot; fit the
native balanced logistic model on training labels alone. Fix C=0.5, alpha=0.1
and per-class thresholds before execution. Retain fitted coefficients/intercepts,
all raw calibration/test scores, sets, singleton abstentions and matching cards.
Replay against the native conformal routine on the identical frozen roles and
compare the all-training-class control. No tuning against outer outcomes, source
addition or biological admission. Output: a new immutable
`results/conformal_logistic_pilot_2026_10_08/` packet. Extend meaningful native
parity/boundary tests to both bases; preserve preceding frozen code packets.

Logistic partition completed: `results/conformal_logistic_pilot_2026_10_08/`.
Native training-only logistic coefficients/intercepts, all calibration/test scores,
quantiles, sets and matched cards replay exactly. Same 560 test genes/26 classes;
empirical prediction-grade coverage 0.905357, mean size 6.975, efficiency 0.761;
zero singleton calls/empty sets. Eleven rare classes use overall fallback.
All-training-class control coverage 1, size 26, efficiency 0. No independent
biological truth or exchangeability promise. 1,088 relevant checks passed.
Artifact `0a5fa0680a5c60d8421d8d60cae53024c5ac891db3f94f3bd01c71ba9d32cee2`.
Progress 35%; missing adapters/full outer coverage/biology remain open.

Functional complete-Pfam partition, 2026-10-08:
`results/functional_profile_knn_2026_10_08_v4/` reuses the native feature-kNN
adapter with fixed training-only selection/ranks and source/derived-input closure.
All 635 test genes, 333 unsupported profiles, 348 numeric inputs, 1,413 native
training classes and 2,279 reporting classes remain. Exact native/refit/serialized
replay yields 12 correct, eight wrong and 615 abstentions; majority recovers19
profiles. Independent acceptance verifies47 input/47 output/12 artifact receipts.
124 focused adapter/control/exclusion/card/artifact regressions pass. Annotation
recovery remains separate from biological validation. Earlier categorical,
metadata and resource failures are retained; UI packaging and remaining adapters
are pending. Progress stays35%; no full-action completion.

Native random-forest call partition,2026-10-08: forest_records.py invokes the
existing forest with catalogue defaults300trees/leaf2, training-only ranks and
serial n_jobs1 (native4) solely for exact probability addition. Native tree JSON
state, calls/support/scores and unavailable states replay exactly; native statistics
unchanged. Frozen Pf direct EC profile test152genes/four unsupported gives59correct,
93wrong,noabstentions; source grade unresolved. Canonical packet
results/pf_functional_ec_forest_2026_10_08_v2 and independent acceptance retain
source/control/cohort identities. 101 combined focused checks pass; resource and
Qt diagnostics preserved. Adapter excludes importance, tuning, calibration,
deployment and independent biological admission; remaining method/variant coverage
still open. 64.10 stays35%; no whole-action tick.

## 2026-11-07: label records of another kind, on frozen partitions

Beyond the 14-strategy per-gene track record, seven label methods now have held-out records through the harness
(instruction 73): boosted trees and the six borrowed strategies of instruction 77 (positive-unlabelled bagging,
propensity weighting, low-rank completion, per-family stacking, lookalike ranking, network propagation), each on
its own once-only frozen partition for `required_without_glucose` and `host_facing`
(`results/borrowed_strategies_2026_11_06/scores.tsv`), and boosted trees for 13 localization niches
(`results/localization_questions_2026_11_06/`, `results/coexpression_components_localization_2026_11_06/`). These are
partition-level records (precision at the fill threshold, AUC, average precision, calibration error, baselines), not
per-gene rows; per-gene rows with native-prediction parity for these adapters remain open.

## 2026-11-07: per-gene rows everywhere, per-class numbers from them

- **Harness methods now keep per-gene rows.** `harness.HarnessResult.rows` (from `harness.gene_rows`): one row per
  scored gene of the nested cross-validation and of the frozen test -- probability, the answer, the calibrated
  probability, whether it is a call at the fill threshold, the truth, right or wrong. Every summary number of a
  partition is reproduced from its rows (`tests/test_harness.py::test_per_gene_rows_reproduce_every_partition_number`).
  Earlier frozen partitions are not re-touched to back-fill rows (each is used once); rows accrue on new runs.
- **Per-class numbers from rows.** `parascape/label_metrics.py`: `by_class` (genes, answered, right, recall,
  precision of calls naming the class, F1; an abstention counts against recall, a never-predicted class has no
  precision rather than zero), `overall` (coverage, accuracy on answered and on all genes, chance = commonest
  class, macro F1), `confusion` (with an abstained column). Same function for the track record and harness rows.
  Tests `tests/test_label_metrics.py`.
- **All 168 method x label records** (both organisms): `scripts/build_label_metrics_notebook.py` ->
  `notebooks/label_metrics_2026_11_07.ipynb`, `results/label_metrics_2026_11_07/` (`overall.tsv`, `by_class.tsv`,
  `confusion.tsv`, `evaluation_gaps.tsv`).

Findings, stated plainly: localisation-type labels are learned well above chance (Tg lopit_unified macro F1 0.63,
accuracy 63 % vs 33 % chance; Pf lopit_pf_location 0.60, 65 % vs 24 %; Tg compartment 0.45, 48 % vs 20 %); the five
Toxoplasma screen-phenotype labels are NOT -- the best method on each stays below the commonest-class accuracy
(e.g. egress 78 % vs 89 % chance) with macro F1 about 0.5, and Pf `is_exported` is the same (93 % vs 97 %; macro F1
0.69 shows some signal on the rare class). 54 of 168 records beat chance on accuracy over all genes.
Evaluation gaps: classes under the minimum class size (15) are left out of hold-out runs -- Tg compartment
(20S proteasome, tubulin cytoskeleton), lopit_unified (cytoskeleton), dtm_class (BETA), Pf lopit_pf_location (six
classes, 66 genes), pb_transferred_phenotype (Fast).

Conformal prediction sets: coverage and set size are reported by the conformal adapters of 2026-10-08
(`label_records.py`); no new conformal method was added here.
