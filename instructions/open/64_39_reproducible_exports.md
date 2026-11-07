# 64.39 · Export the evidence and inference trail

Status: OPEN — about 70 %, 2026-11-06. One bundle per Infer answer built and tested; adapters for the other modules open.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.08, 64.17, 64.30, 64.31, 64.35 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `parascape/results.py`, `parascape/report.py`, `parascape/datasets.py`, `parascape/strategies.py` |

## Deliverable

Exportable result tables, scorecards, test records, figures and a concise methods/citation/provenance bundle for the selected query.

## Controlled scope

One common export bundle and adapters for the new profiles; reuse instruction 60's methods/figure work.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] The bundle identifies organism/entity/context, source and artifact hashes, method/settings, benchmark/split and whether outputs are measured/inferred.
- [ ] CSV/table and vector figure exports reconcile with the visible selected data; resolved source citations trace to verified registry records.
- [ ] A reloaded export recovers the same query and scorecard evidence; missing sources stay explicit instead of generating invented citations.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-06: the bundle

`parascape/export_bundle.py`: `write(folder, ...)` puts each table as CSV, **marked measured or inferred** (a
table without the mark is refused), each figure as a vector SVG (`figure_export`), `methods.md` (methods
paragraph, references and BibTeX from the registry; unconfirmed sources flagged, never cited) and
`manifest.json` (query, organism, version, table fingerprint, evidence columns, SHA-256 of every file). `load`
refuses a bundle edited after it was written. In the application: Infer ▸ Save the answer ▸ **Save everything as
a bundle…** (named genes inferred, recovery measured, the Map view as a figure). Tests:
`tests/test_export_bundle.py` (3), `tests/test_module_shell.py` (the Infer bundle). Open: adapters for the Map,
Gene and Network modules, and reloading a bundle back into the application.
