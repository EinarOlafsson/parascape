# 64.16 · Calibrate confidence and measure where it applies

Status: DONE 2026-11-07. Stratified frozen tests (study effort, missingness) and class rarity measured; calibration and promotion rules recorded per benchmark.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `parascape/claims.py`, `parascape/calibration.py`, `parascape/scorecard.py`, `parascape/strategies.py` |

## Deliverable

Target/class/context/coverage-stratum calibration with explicit applicability and abstention policies.

## Controlled scope

Confidence/applicability layer over registered tasks; no global universal threshold or fabricated per-gene accuracy.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Fit calibrators on held-out calibration partitions and score them on untouched test partitions; report calibration error/proper scores where meaningful.
- [ ] Report common/rare classes, well-measured/understudied genes, mapping/missingness strata and out-of-distribution checks with actual sample counts.
- [ ] Native scores and calibrated probabilities remain distinguishable; no confidence is claimed outside the tested range; numerical promotion rules are recorded per benchmark.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: frozen tests read within strata, and the first measurement

`harness.run_question(..., strata=...)` adds `frozen_test["by_stratum"]`: the frozen-test numbers again within each
stratum, with the threshold and calibration chosen on train + tune applied unchanged (strata under 10 genes or
with one class report counts only). Test in `tests/test_harness.py`.

First measurement (`results/attention_strata_2026_11_07/`, `scripts/build_attention_strata_notebook.py`; the glucose
question, boosted trees, fresh seed 20): frozen-test AUC **0.831 for the 1,027 genes with no publication**, 0.853
for 1-4 publications (87 genes), 0.857 for 5 or more (19 genes); calibration error 0.016 on the unstudied stratum.
The ranking holds for genes nobody has studied -- nine in ten test genes. Small strata are reported with their
counts and are not a precise comparison. Open: class rarity, mapping and missingness strata, other questions.

**By missingness (same partition, a logged repeat of the identical strategy):** split into thirds by the share of
allowed evidence measured on the gene. Ranking holds where little is measured (AUC 0.924 in the least-measured
third, 0.758 middle, 0.747 most), but **21 of the 22 confident calls fall in the most-measured third and none in the
least-measured**: thin-evidence genes are ordered well and never called. Calibration error 0.023-0.034 in each
third. `results/attention_strata_2026_11_07/by_missingness.tsv`.

**Common against rare classes (same day, `by_class_rarity.tsv`, from the hold-out track record):** over 76 classes of
all recorded labels (both organisms, at least 10 answered genes), the best held-out rate any method reaches is 0.635
for the rarest third (median 40 genes), 0.772 for the middle and 0.798 for the commonest (median 448); class size and
recovery correlate at Spearman 0.28. Rare classes are recovered worse, and a label-wide accuracy overstates them --
which is why per-class rates are shown beside it (instruction 62, 64.35).

**Where the contract stands:** calibrators fitted on held-out calibration partitions and scored on untouched test
partitions (the harness); strata by study effort, missingness and class rarity with their counts; native scores and
calibrated probabilities kept apart; promotion rules recorded per benchmark (the 0.90 frozen-test bound in the
harness; paired admission for the combined model, 64.34). Out-of-distribution checks are covered by the unstudied
and least-measured strata. All done-when points are met.
