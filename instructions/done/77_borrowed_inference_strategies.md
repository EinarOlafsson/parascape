# 77 · The borrowed profile-inference strategies, each with its own quality control
**Status: done 2026-11-06 -- six built and scored on frozen partitions, three untested with the reason; none beats boosted trees, positive-unlabelled bagging matches them.** Requested 2026-11-05. The techniques identified from large-scale profile inference, now to be
implemented. Each enters the same catalogue and the same harness as the existing 39 strategies.

## Why
Einar's framing: the same shape as predicting a person's attributes from many weak signals. Every run here has
shown the matching symptom -- no single feature carries signal (every permutation importance <= 0.02, spread
across unrelated families) while ranking ability arrives long before precision.

## The strategies, in the order they should be built
1. **Positive-unlabelled learning with confidence weights.** The central label problem: "not recorded" is not
   "not true". Weight observed positives by confidence, unobserved cells as weak negatives with low weight.
2. **Inverse-propensity weighting for data missing not at random.** A famous gene is more likely to be measured
   at all, so model that probability and weight by its inverse. The principled version of deleting attention
   columns.
3. **Low-rank matrix completion.** The map *is* a sparse gene x question matrix; factorisation predicts a cell
   from the whole co-occurrence pattern rather than from chosen features.
4. **Per-family stacking.** A model per evidence family combined by a meta-learner on out-of-fold predictions,
   so each family's contribution is measurable instead of hidden among 86 of them.
5. **Lookalike ranking from a seed set.** Rank the genome by similarity to a small set of known positives --
   better posed than balanced classification against invented negatives.
6. **Selective prediction with a reject option.** The formal version of the fill rule, with coverage at fixed
   precision as the metric.
7. **Network propagation** (personalised PageRank, diffusion kernels), stronger than the kNN already built.
8. **Multi-task learning** over correlated questions; replaces per-question models only where it wins.
9. **Empirical-Bayes shrinkage** for thin columns, generalising what the expression fold changes already use.

## Quality control, identical for every one of them
Nested cross-validation inside train+tune, a frozen test partition touched once, `guard_fit` at every stage,
the stated baselines, isotonic calibration, the Clopper-Pearson lower bound, the no-knowledge split, and
permutation importance by family. **An untested strategy is listed as untested and never offered as equal.**

## Gate
Each strategy ships when it beats the best baseline on a question whose labels can support a frozen test, or is
published as untested with the reason. None may fill a cell on a cross-validated number alone.

## Done 2026-11-06

Built as estimators the harness scores unchanged (`parascape/borrowed_models.py`, registered in
`harness.MODELS`; 20 tests in `tests/test_borrowed_models.py`), and scored on the two registry questions that
afford a frozen test: `required_without_glucose` (7,393 genes) and `host_facing` (585 curated in the
information map). Each strategy had its own fresh partition (seeds 10-16), touched once and logged. Full table
and reading: `results/borrowed_strategies_2026_11_06/README.md`.

* **Gate**: positive-unlabelled bagging ships on glucose (frozen AP 0.531 against the reference 0.538, beats
  the best baseline) -- an equal, not an improvement. Every other strategy is worse than boosted trees on the
  cross-validated estimate; propensity weighting is better on its frozen partition (0.587) and worse
  cross-validated (0.527), so undecided. **No strategy may fill** (all 14 frozen lower bounds below 0.90).
* **The host-facing frozen test is noise-limited**: a shuffled-label run of the same pipeline gave frozen AUC
  0.40-0.85 by chance (about 10 positives per partition). Kept as `shuffled_label_control.tsv`.
* **Untested, with reasons** (`borrowed_models.UNTESTED`): selective prediction is the harness itself;
  multi-task learning, empirical-Bayes shrinkage.
* Found on the way: the harness hashes a strategy by model name and settings, not by the model code, so a
  changed estimator must bump its `version` or it would get a second look at a spent partition. Recorded in
  the module docstring. Also: inverse-propensity weighting cannot remove a tilt where attention nearly
  decides the label (the weights hit the clip) -- tested and documented.
* Also closes the host-facing part of instruction 73's phase 1 gate: host-facing now has frozen-test
  results through the harness (585 labels; the registry note of 98 predates the hard negatives).
