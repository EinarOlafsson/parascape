# 64.27 · Precompute the Plasmodium inference atlas

Status: DONE — 100%, 2026-11-07. Plasmodium atlas manifest (all 39 methods with a role and a status), query cache by identity, and sampled cached records reproduced exactly by a fresh run. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.07, 64.17, 64.25, 64.26 |
| Estimated remaining engineering time | 4–8 h + compute |
| Existing code to inspect | `parascape/organisms.py`, `parascape/strategies.py`, `parascape/track_record.py`, `parascape/claims.py`, `parascape/packs.py` |

## Deliverable

A Pf atlas using its own targets, graph, available evidence, split/calibration artifacts and fixed settings manifest.

## Controlled scope

One parasite species and a frozen pilot manifest; expand targets in a later bounded partition.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] The same atlas admission/parity checks as Tg pass using Pf's own gene and target spaces.
- [x] Missing evidence/ground truth stays explicit; Tg metrics, defaults or models cannot silently supply Pf results.
- [x] Build runtime, RAM, size, coverage and incomplete expensive queries are recorded in the manifest.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the atlas manifest, canonical queries and a reproducibility check

`parascape/atlas.py`: `manifest(organism)` gives every one of the 39 catalogue methods its role (gene label,
numeric value, gene ranking, module, replication, gene set) and, for known genes, unknown genes and its card,
either what is stored or a specific status. Known genes are read only from held-out per-gene rows; unknown genes
only from deployment claims; a method with no stored deployment output says it runs on demand from its module.
`cached_or_run` serves a query (organism, method, label, genes, settings) from the atlas cache by its identity
and runs only what is missing; a different list is a different query. Tests `tests/test_atlas.py`.

Executed `notebooks/inference_atlas_2026_11_07.ipynb` (`results/inference_atlas_2026_11_07/`): manifests for
both parasites, canonical class-set queries already cached in the track record (`canonical_set_queries.tsv`), and
a reproducibility check -- 4 cached (method, label) records per organism drawn at random (seed 7) and rebuilt from
scratch with the production computation (`precompute_track_record`):
cluster_guilt, stratum_focus, random_forest and graph_convolution on lopit_pf_location -- 1,580 genes each,
every call identical, accuracy identical (0.144, 0.492, 0.646, 0.594).

Coverage on Plasmodium: 14 of 39 methods have held-out per-gene rows (on the three Pf labels), 2 have deployment
claims, all 39 have a measured card; the rest run on demand. Remaining coverage, stated: as for Toxoplasma
(64.26), plus six Pf location classes under the minimum class size (64.10 evaluation gaps).
