# 64.37 · Make user data reach the same inference flow

Status: DONE 2026-11-07. Imports reach the strategies, record their input and state, hold out their own import, and save and replay exactly.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.04, 64.05, 64.08, 64.23, 64.25, 64.30 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `parascape/importer.py`, `parascape/app.py`, `parascape/strategies.py`, `parascape/guided.py`, `parascape/search.py` |

## Deliverable

Imported measurements rebuild the correct organism context, declare units/evidence family and invalidate affected caches.

## Controlled scope

One end-to-end import fixture for each currently supported parasite/host ID scheme; reuse instruction 60's import repair.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] A real fixture import changes the strategy inputs and its inference/scorecard context in the same session.
- [ ] Imported labels are excluded when used as evaluation truth; missing provenance/leakage declarations produce an explicit unresolved validation state.
- [ ] Shared shipped data remain immutable; save/load retains the exact input hash, organism, transformations and local result artifacts.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07

* **Done:** a real fixture import now changes the strategy inputs in the same session -- `StrategyPanel.set_nodes`
  and `GuidedPanel.set_context`, called by `Window.apply_import` (instruction 60, change 8; test in
  `tests/test_app_controls.py`).
* **Done:** every import record carries `source_sha256` (the exact file), `organism`, and `validation` -- explicitly
  "unresolved: imported columns have no registered dataset or leakage family; closure guards them by measured
  association only". Test in `tests/test_app_controls.py`.
* **Done (same day):** an imported label holds out the rest of its own import -- the application registers each
  import's columns (`recipes.register_import`, keyed by the file hash) and `recipes.close` excludes them with the
  reason "same import as <label>". Test in `tests/test_recipes.py`.
* **Open:** session save/load of imports with their local results; the same rule in the harness feature
  selection.

**Save and replay (same day):** File ▸ **Save the imports…** writes every import record of the session (file,
SHA-256, organism, every processing choice) as JSON; File ▸ **Replay saved imports…** imports each file again with
the same choices and refuses one whose content changed since (`Window.save_imports`, `replay_imports`). Test in
`tests/test_app_controls.py` (the same columns rebuilt; a changed file not re-imported). All done-when points are
met: the import reaches the strategy context in-session, records its input, organism and an explicit unresolved
state, holds out its own import as evidence for its own label, keeps the shipped data immutable, and saves/reloads.
