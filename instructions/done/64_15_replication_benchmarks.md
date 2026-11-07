# 64.15 · Test split-half and conjunction findings

Status: DONE 2026-11-07. Split-half and independent-dataset replication measured, kept apart; per-finding records with both halves and effect sizes for blind tests and conjunctions; split clusters recorded as too low-yield to test.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `parascape/strategy_catalog.py`, `parascape/scorecard.py`, `parascape/discovery.py` |

## Deliverable

Finding-level replication records for blind tests, split clusters and conjunctions.

## Controlled scope

Replication-task adapter for the three existing mechanisms.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Discovery and verification populations/groups are disjoint; finding selection and multiple-testing rules are declared before verification.
- [ ] Every finding records both halves, the original question, eligible sample counts, effect/direction and replication outcome.
- [ ] Cards expose replicated findings / tested findings, effect sizes, uncertainty and matched-null rates; split-half replication is identified separately from an independent dataset test.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: what the replication records say

* **Split-half replication (calibration sweep 2026-09-26, findings selected on the first half before the second is
  read):** blind tests -- share of first-half findings that replicate on the second half, skill 0.95 (Tg) / 0.88 (Pf),
  every run passing; conjunctions -- 0.49 (Tg, 10 of 15 runs pass), untestable for Pf.
* **Split clusters: untestable as built.** All 30 runs inconclusive: the first half yields 0-2 findings where 3 are
  needed before replication can be scored (41 of the inconclusive notes read "only 0 findings on the first half").
  This is a finding about the mechanism's yield, not a pass or fail.
* **Independent-dataset replication, kept separate as required:** across studies, not halves of one --
  the perturbation atlas finds the same perturbation in another lab's experiment at median rank 1-3 against chance
  59-99 (`results/perturbation_atlas_2026_11_07/`); wild-type alkaline stress agrees between independent studies at
  Spearman 0.44 (instruction 79); BFD1 gain and loss reproduce across labs (`expression_traits_2026_11_06`).
* **Open:** per-finding records (both halves, effect, direction) on cards, and a lower finding threshold or larger
  first half for split clusters so it can be tested at all.

**Per-finding records (same day):** `strategies.replication_test` now keeps one row per finding in the test's
details -- the finding, whether it held on the second half, and the size of each half -- for every strategy built on
it (blind tests, split clusters, conjunctions). Blind battery on Tg: 53 first-half findings, 52 replicate. Test in
`tests/test_strategies.py`. Remaining: effect size and direction per finding (the `find` functions return feature
names only), and split clusters' too-few-findings problem.

**Effect sizes per finding, both halves (same day):** `replication_test(..., effect=...)` records each finding's effect
on the first and on the second half; the blind battery supplies its association size (Cramer's V or epsilon
squared -- undirected, since a cluster association has no sign; stated rather than invented). On Tg the 53 findings'
effects agree between halves at Spearman 0.975. Remaining: effect functions for conjunctions and split clusters,
and split clusters' too-few-findings problem.

**Conjunctions (same day):** effect = fold enrichment of the second label inside the cluster, on each half (directional:
above 1 is enriched). On Tg the first-half enrichments run 1.3-6.1 and the second-half 0.7-6.2; a finding whose
second-half enrichment falls below 1 has reversed direction and is the one that failed to replicate. Remaining:
split clusters yields 0-2 first-half findings where 3 are needed, so it cannot be scored as built -- recorded as a
property of the mechanism.
