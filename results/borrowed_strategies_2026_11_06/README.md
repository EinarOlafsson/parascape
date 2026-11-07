# The borrowed inference strategies, scored by the harness (instruction 77)

Notebook: `notebooks/borrowed_strategies_2026_11_06.ipynb` (executed). Builder:
`scripts/build_borrowed_strategies_notebook.py`. Models: `parascape/borrowed_models.py`, registered in
`harness.MODELS`. Each strategy ran through `harness.run_question` unchanged (nested 3 x 3 cross-validation,
guard_fit, isotonic calibration, baselines, family importance) and touched its own frozen partition once
(seeds 10-16, logged in `results/harness_holdout_log.jsonl`).

## Results (average precision; `scores.tsv` has every number)

| strategy | glucose: cross-validated | glucose: frozen test | host-facing: cross-validated | host-facing: frozen test | verdict (glucose / host-facing) |
|---|---|---|---|---|---|
| boosted trees (reference) | 0.581 | 0.538 | 0.721 | 0.962 | reference |
| positive-unlabelled bagging | 0.572 | 0.531 | **0.735** | 0.592 | **ships** (equal, not better) / worse |
| inverse-propensity weighting | 0.527 | **0.587** | 0.687 | 0.624 | worse on cross-validation, better on the frozen test |
| low-rank completion | 0.527 | 0.542 | 0.634 | 0.975 | worse |
| per-family stacking | 0.542 | 0.578 | 0.541 | 0.745 | no confident calls at the bar |
| lookalike ranking | 0.479 | 0.455 | 0.681 | 0.969 | worse |
| network propagation | 0.487 | 0.451 | 0.369 | 0.432 | worse |

Chance (prevalence): glucose about 0.24, host-facing about 0.11. **No strategy may fill a cell** on either
question: the frozen-test lower bound is below 0.90 for all 14.

## How to read the host-facing frozen column: the shuffled-label control

The same pipeline on randomly permuted labels (`shuffled_label_control.tsv`, `DRYRUN=1`) gave host-facing
frozen-test AUC from 0.40 to 0.85 and average precision from 0.08 to 0.46 -- by chance alone, because a
frozen partition holds about 10 host-facing positives. Differences between strategies on that column are
inside the noise; the cross-validated column is the steadier comparison. Glucose (about 260 positives per
frozen partition) varied only 0.47-0.51 AUC under shuffling.

## Reading

* **Nothing borrowed beats boosted trees.** Positive-unlabelled bagging matches them on glucose and is
  marginally best on host-facing cross-validation; it is the only one that clears the shipping rule.
* **Propensity weighting splits**: worse cross-validated, better on its frozen partition. Undecided until a
  second partition, and the cross-validated estimate is the one to trust.
* **Lookalike and propagation are worse** on glucose. The likely reason, not measured here: with ~370 evidence columns, distance in standardised
  space is dominated by columns that do not bear on the question, which trees ignore and a similarity cannot.
* **Stacking makes no confident calls** at the threshold the bar needs. Likely, not measured: a non-negative combination of 77 family models smooths its
  probabilities below it.
* **Untested, with the reason** (`untested.tsv`): selective prediction is the harness itself
  (coverage at fixed precision); multi-task learning waits on the registry recording which labels may stand
  beside which; empirical-Bayes shrinkage needs per-cell standard errors only the expression fold changes have.
