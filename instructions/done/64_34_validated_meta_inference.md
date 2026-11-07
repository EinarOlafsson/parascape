# 64.34 · Train and test the meta-inference

Status: DONE 2026-11-07. The pilot target (dtm_class) admitted under strict nesting; further targets are separate admissions.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.08, 64.16, 64.32, 64.33 |
| Estimated remaining engineering time | 8–16 h + compute |
| Existing code to inspect | `parascape/claims.py`, `parascape/strategy_learning.py`, `parascape/scorecard.py`, `parascape/calibration.py` |

## Deliverable

A dependence-aware meta-model per compatible output task/target, initially piloted on one target, then admitted by the same task contract.

## Controlled scope

One target and one output task in the pilot; each additional target/task is a separate measured admission partition.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Train with nested out-of-fold base predictions produced without the outer held-out genes/groups; inner fitting/tuning/calibration cannot see outer truth.
- [ ] An untouched outer test compares meta-inference with the best base selected on inner training, fixed consensus and matched null/baseline.
- [ ] Report per-class/stratum error, calibration, coverage and incremental capacity; superiority is measured on paired/grouped test units.
- [ ] Persist weak/failed cases; an uncalibrated or unsupported meta-model cannot replace a measured method answer.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the pilot on `dtm_class` (Tg)

`scripts/build_meta_inference_pilot_notebook.py` -> `results/meta_inference_pilot_2026_11_07/`. Base: the held-out calls
of the 14 recorded methods on 8,137 genes (4 classes). Meta-model: multinomial logistic regression on the one-hot
calls, fitted on 4/5 of the orthogroups, tested on the other 1/5, five times. On the same untouched test genes:

| method | accuracy (mean, range over folds) |
|---|---|
| **meta-model** | **0.918** (0.899-0.932) |
| best single method, chosen on the training orthogroups | 0.903 (0.888-0.917) |
| plurality vote | 0.819 (0.801-0.841) |
| commonest class | 0.744 |

Paired on the same genes: the meta-model is right where the best single method is wrong on 223 genes, the reverse on
103 (McNemar p < 0.001). The plain vote is worse than the best single method -- counting agreeing methods is not a
good combiner; learning whom to trust is. **Not yet admitted:** the base calls come from the hold-out record's own
folds, so a base call for an outer-training gene may come from a model that saw an outer-test gene's label; admission
needs the base methods recomputed inside each outer fold, plus per-class error and calibration of the meta
probabilities.

**Per class and calibration (same day, out-of-fold over the five folds):** recall GLOB 0.98 (6,053 genes), TM 0.94
(1,142), SP 0.54 (676), **SP+TM 0.28 (266)** -- the overall accuracy is carried by the two large classes, and the
rare combined class is mostly missed; expected calibration error of the meta probabilities 0.007 (well calibrated).
Admission still needs the strictly nested base calls.

## 2026-11-07: admitted on `dtm_class` under strict nesting

`scripts/build_meta_inference_nested_notebook.py` -> `results/meta_inference_nested_2026_11_07/` (one worker per base
method, calls stored per method). 5 outer folds of whole orthogroups; inside each, every base method's calls on the
training genes come from 3 inner folds that never see an outer-test label (32,548 inner calls per method), and its
calls on the test genes from models fitted with the test labels hidden (8,137). The meta-model learns from the inner
calls only and is scored once on the outer-test genes:

| method | accuracy (mean, range over outer folds) |
|---|---|
| **meta-model** | **0.919** (0.912-0.923) |
| best single method, chosen on the training genes | 0.902 (0.878-0.916) |
| plurality vote | 0.817 |
| commonest class | 0.744 |

Paired: meta right where the best single is wrong on 244 genes, the reverse on 112 (McNemar p < 0.001). **Admitted
for `dtm_class`**, with the stated weakness: recall on the rare SP+TM class is low (0.28 in the pilot) and the
probabilities are well calibrated (0.007). Each further target is its own admission partition (this contract);
weak and failed cases are kept in `folds.tsv`. A meta answer replaces nothing on its own: it is shown beside the
methods it combines.
