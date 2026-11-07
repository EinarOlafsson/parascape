# 64.36 · Connect every evidence and inference level

Status: DONE 2026-11-07. The navigation chain checked on real widgets: scorecard in one click, widths, missing data, back/search, organism change, host counterpart.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.18, 64.19, 64.20, 64.23, 64.24, 64.30, 64.31, 64.33, 64.35 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `parascape/app.py`, `parascape/guided.py`, `parascape/guided_panel.py`, `parascape/strategy_panel.py`, `parascape/help_search.py` |

## Deliverable

A consistent navigation chain: organism/dataset ↔ label/class ↔ gene ↔ strategy/method ↔ meta-inference ↔ test record/source.

## Controlled scope

Navigation/integration and accessibility checks; preserve the existing layout where it can serve the route.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Each required scorecard is reachable within one action from its corresponding result, with deeper truth/split/evidence records from that card.
- [ ] Back/search and organism changes retain valid query context and never show stale results from another entity.
- [ ] Real-widget checks cover a gene, class, label, host counterpart and missing-data case at supported window widths; compact summaries expand on demand.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: gene -> method -> scorecard in one click

In the Gene module, every method of **What each method said about it, held out** is a link: one click opens the
Strategies module with that strategy's card (guide, calibration, scorecard, track-record line) shown
(`GeneModule.open_strategy`). Test in `tests/test_module_shell.py`. Also in place from earlier work: gene -> class
and label pages (track record), question -> its best-recovering strategies (Infer), result -> methods text and
bundle. Width and missing-data checks (same day): the Gene module reads the best- and the least-measured gene at 1024,
1366 and 1920 px on the real window, each section saying something (the unusual list or why there is none, the
method list) and the module never wider than the window (`tests/test_module_shell.py`). Back and search (same day): going Home and opening another gene shows that gene and nothing of the last one
(`tests/test_module_shell.py`); the Navigator resets on return (`tests/test_navigator.py`). Open: the host
counterpart (no host gene space ships yet, instruction 64.21/64.22), organism change.

**Organism change and host counterpart (same day):** switching organism opens the other species in a window of its
own (by design: the tables are never merged), and its Gene module reads that organism's genes with nothing of the
previous one (`tests/test_app_menus.py`); a parasite protein with a recorded host target (e.g. ROP16 -> STAT3) names
its host counterpart in what the gene binds (`tests/test_module_shell.py`). With the gene -> method -> scorecard link,
widths, missing data and back/search checks, all done-when points are met; a full host gene space (64.21/64.22)
will add host-side pages to the same chain.
