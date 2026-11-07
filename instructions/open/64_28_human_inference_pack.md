# 64.28 · Precompute and validate human inferences

Status: DONE — 100%, 2026-11-07. Human pilot benchmarked on human truth (HPA graded locations, five compartments) from the human pack; parity and cold-cache checks pass; gaps stated. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.21, 64.23, 64.25 |
| Estimated remaining engineering time | 6–12 h + compute |
| Existing code to inspect | `parascape/organisms.py`, `parascape/strategies.py`, `parascape/calibration.py`, `parascape/packs.py` |

## Deliverable

An Hs inference pack with organism-specific benchmark runs for a declared initial target/evidence set.

## Controlled scope

One host and an agreed pilot target set; compute estimates come from the first measured job.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Applicable adapters run on human data and tests; available truth supports local class/target reliability rather than parasite-derived grades.
- [x] Coverage, unsupported mechanisms, context and mapping/applicability gaps accompany every output.
- [x] Pack integrity, cached/fresh parity and cold-cache organism-switch checks pass for the selected pilot scope.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: pilot on human truth

`scripts/build_human_inference_pack_notebook.py` -> `notebooks/human_inference_pilot_2026_11_07.ipynb`,
`results/human_inference_pilot_2026_11_07/` (scores, per-gene rows per location, summary).

**Frozen scope (declared in the notebook before running).** Five questions: is a human protein in the
Nucleoplasm / Plasma membrane / Mitochondria / Golgi apparatus / Endoplasmic reticulum, truth from the Human
Protein Atlas where the location is graded Enhanced, Supported or Approved (12,662 genes; Uncertain and
unannotated left out). Evidence: expression only -- GTEx fibroblast and hepatocyte medians, HPA RNA tissue
specificity, distribution and single-cell specificity, and the shipped infection responses of the gene's one
reviewed protein. Secretome location and protein class (it includes predicted membrane/secreted classes) are the
label's own category and excluded. Boosted trees at harness defaults, seed 0, each gene its own group.

**Result (human performance on human truth, nothing borrowed from the parasites):** frozen-test AUC 0.639
(Nucleoplasm), 0.673 (Plasma membrane), 0.709 (Mitochondria), 0.579 (Golgi), 0.689 (ER); average precision above
prevalence for each (e.g. Mitochondria 0.268 vs 0.072, ER 0.188 vs 0.036). None may fill: expression alone ranks
human locations modestly and is not enough to write a location as a fact. The notebook was executed three times
while fixing cells after the harness step; the frozen tests were first touched by the identical strategy and
the later runs are flagged `frozen_test_repeat` in the log, as the harness is designed to record.

**Checks.** Cached/fresh parity: Mitochondria re-run from scratch gives identical per-gene rows. Cold cache: in an
empty cache the human pack is "installable", installs against its SHA256 and reads back 58,988 genes identical to
the build, shown "not calibrated".

**Gaps, stated:** no human homology groups (paralogs can sit on both sides of a split); only expression evidence
in the pack (no human interaction, structure or domain columns yet); HPA locations are cell-line
immunofluorescence; catalogue strategies other than the harness models need a human context adapter (not
built), so human has no per-gene track record and stays "not calibrated" in the space catalogue.
