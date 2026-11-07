# 64.29 · Precompute and validate mouse inferences

Status: DONE — 100%, 2026-11-07. Mouse pilot benchmarked on mouse truth (BMDM surface proteomics) from the mouse pack; parity and cold-cache checks pass; gaps stated. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.22, 64.23, 64.25 |
| Estimated remaining engineering time | 6–12 h + compute |
| Existing code to inspect | `parascape/organisms.py`, `parascape/strategies.py`, `parascape/calibration.py`, `parascape/packs.py` |

## Deliverable

An Mm inference pack with its own measured target/context reliability for a declared initial evidence set.

## Controlled scope

One host and a fixed pilot target set; additional contexts are distinct partitions.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Applicable adapters run on mouse data and tests; human/parasite performance is never reused as measured mouse performance.
- [x] Tissue, infection/context, missingness and mapping differences remain visible in results and scorecards.
- [x] Pack integrity, cached/fresh parity and cold-cache organism-switch checks pass for the selected pilot scope.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: pilot on mouse truth

`scripts/build_mouse_inference_pack_notebook.py` -> `notebooks/mouse_inference_pilot_2026_11_07.ipynb`,
`results/mouse_inference_pilot_2026_11_07/`.

**Frozen scope (declared in the notebook before running).** Question `mm_bmdm_surface`: is a mouse protein
detected on the surface of bone-marrow-derived macrophages (surface proteomics shipped in
`mm_host_proteins.parquet`: 1,296 proteins tested, 150 detected; untested proteins are unknown, not negative).
Evidence: mouse expression only -- BMDM, brain and skeletal-muscle TPM and the Bgee summaries of the protein's gene
(one reviewed link, from the mouse pack). Excluded: the surface intensity (the label's own source) and every
location column, including the locations transferred from human orthologs (the label's own category, and
inferred). Boosted trees at harness defaults, seed 0, each protein its own group.

**Result (mouse performance measured on mouse truth, nothing borrowed):** cross-validated AUC 0.760, average
precision 0.270; frozen test (195 proteins, 22 detected, prevalence 0.11) AUC 0.798, average precision 0.279.
At the fill threshold the frozen test made 5 calls, 0 right: it ranks surface proteins well above chance and may
NOT fill. Per-gene rows saved (`rows.parquet`).

**Checks.** Cached/fresh parity: the question re-run from scratch gives identical per-protein rows (flagged a
repeat and logged). Cold cache: in an empty cache the mouse pack is "installable", installs against its SHA256,
then reads back 18,428 genes identical to the pack, shown "not calibrated".

**Gaps, stated:** one observed mouse label in the pilot; no mouse homology groups (paralogs can sit on both sides
of a split); BMDM surface proteomics measures one cultured cell type; catalogue strategies other than the harness
models need a mouse context adapter (not built), so mouse has no per-gene track record and stays "not
calibrated" in the space catalogue.
