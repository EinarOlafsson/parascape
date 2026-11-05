# 77 · The borrowed profile-inference strategies, each with its own quality control
**Status: open (requested 2026-11-05).** The techniques identified from large-scale profile inference, now to be
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
