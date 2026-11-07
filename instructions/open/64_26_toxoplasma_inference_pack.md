# 64.26 · Precompute the Toxoplasma inference atlas

Status: DONE — 100%, 2026-11-07. Toxoplasma atlas manifest (all 39 methods with a role and a status), query cache by identity, and sampled cached records reproduced exactly by a fresh run. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.07, 64.17, 64.25 |
| Estimated remaining engineering time | 4–8 h + compute |
| Existing code to inspect | `parascape/organisms.py`, `parascape/strategies.py`, `parascape/track_record.py`, `parascape/claims.py`, `parascape/packs.py` |

## Deliverable

A Tg atlas over an explicitly frozen initial target/settings manifest and the reusable outputs each applicable strategy can precompute.

## Controlled scope

One organism and a declared pilot target/settings set; record the remaining target coverage rather than calling the atlas complete.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Each of the 39 strategies has an output or a specific per-query status; gene-label, numeric, relationship and module outputs preserve their roles.
- [x] Known genes use held-out records; unknown genes use separately generated deployment predictions; provenance and cards are attached.
- [x] Canonical class/list queries are cached by query identity; arbitrary new lists/parameters reuse base artifacts then run the remaining bounded work on demand.
- [x] Random sampled cached outputs match a fresh run of the same frozen configuration within declared tolerances and fit-time conditions.
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
layer_vote / screen_replication_phenotype, supervised_classifier / screen_actin_phenotype, map_neighbours /
screen_replication_phenotype and map_neighbours / screen_apicoplast_phenotype -- 318 genes each, every call
identical, accuracy identical (0.333, 0.679, 0.862, 0.805).

Coverage on Toxoplasma: 14 of 39 methods have held-out per-gene rows (the label methods on all nine labels), 2
have deployment claims (the calibrated recipes of `claim_recipes.parquet`), all 39 have a measured card; the
other 25 are ranking, value, module, replication and gene-set methods graded by their own held-out scorecards and
run on demand. Remaining coverage, stated: per-gene deployment output for every label method (only the best
calibrated recipe per label is precomputed), and per-gene rows for the non-label tasks.
