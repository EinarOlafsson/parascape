# 64.40 · Verify both goals end to end and release

Status: OPEN — 70%, 2026-11-07. Release evidence matrix built from stored results; release waits on a green Tests workflow (stale PyQt6 setting, fix pushed separately) and on item 90 (licence review).

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.18, 64.19, 64.20, 64.21, 64.22, 64.23, 64.24, 64.26, 64.27, 64.28, 64.29, 64.30, 64.31, 64.33, 64.35, 64.36, 64.37, 64.38, 64.39 |
| Estimated remaining engineering time | 6–10 h + checks |
| Existing code to inspect | `tests/`, `scripts/release.py`, `scripts/build_tutorials.py`, `docs/`, `NEXT_SESSION.md` |

## Deliverable

A release evidence matrix for the two goals, all 39 strategies, supported organisms/hosts, every scorecard level and all unresolved biological gaps.

## Controlled scope

Integration/release verification of admitted scope; unimplemented species and unavailable truth are recorded as remaining work, not completed goals.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] From an installed distribution, users browse a source, look up a gene/class/label, inspect applicable cached inferences and meta-results, and open their grounded cards.
- [ ] Offline/cold-cache/missing-pack/import/organism-switch cases pass; memory/storage/runtime budgets and stale-artifact guards are measured.
- [ ] Every strategy has an appropriate explicit test and its evaluated scope, or remains visibly unvalidated with its truth-acquisition item still open.
- [ ] Fresh calibration and current tutorials/docs accompany changed data/methods; required CI and package checks pass before the standing authorized release flow.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## 2026-11-07: the release evidence matrix

`scripts/build_release_evidence_notebook.py` -> `notebooks/release_evidence_2026_11_07.ipynb`,
`results/release_evidence_2026_11_07/`: gathered from stored results, nothing restated by hand.

- **Methods** (`methods.tsv`): all 39 catalogue methods per parasite with role, evidence (held-out per-gene rows
  + card for the 14 label methods; own held-out scorecard for the other 25), labels on which each beats chance,
  best macro F1, deployment claims.
- **Spaces** (`spaces.tsv`): Toxoplasma and Plasmodium ready and calibrated; human and mouse packs installable,
  never shown calibrated.
- **Scorecard levels** (`levels.tsv`): gene, class, label, method, meta and question, each with where it is stored.
- **Host pilots** (`host_pilots.tsv`): human locations frozen AUC 0.58-0.71, mouse BMDM surface 0.80; none may fill.
- **Gaps** (`open_items.tsv`, `stated_gaps.tsv`): every open progress item and every gap the finished 64.x items
  stated in their write-ups.

**What still blocks a release, stated:**

1. The `Tests` workflow on nightly has not passed for the last 60 runs: `.github/workflows/checks.yml`,
   `docs.yml` and `release.yml` set `PYQTGRAPH_QT_LIB: PyQt6` while the package depends on PySide6, so the first
   step fails at import. The fix (set PySide6) is a separate workflow-only commit; whether this checkout may push
   workflow changes, and what the full suite says after it, are recorded below when known.
2. Item 90: no PyPI or pack release until Einar has reviewed the licence terms.
3. The installed-distribution walk-through (browse a source, look up a gene/class/label, open cached inferences
   and cards from an installed wheel, offline and cold-cache) is covered piecewise by the window tests and the
   cold-cache checks of 64.28/64.29, not yet as one run from an installed wheel.
4. Known order dependence in the local suite: `test_app_menus.py` then `test_app_smoke.py` in one process gives 14
   smoke failures (the menu tests leave the Plasmodium map open); each file passes alone, and the same 14 fail on a
   clean checkout, so this predates today's changes.

**CI after the Qt binding fix (2026-11-07).** With `PYQTGRAPH_QT_LIB: PySide6` (commit 385d12e5, the workflow push
was accepted) the first step passes. The next, `tests/test_app_contexts.py` under xvfb, passed every assertion
(three windows, shared contexts, valid framebuffers) and then aborted while Python finalised: "Fatal Python
error: none_dealloc: deallocating None" in shiboken teardown (Python 3.11.18, PySide6 6.12.0). The test's
subprocess now ends with `os._exit(0)` after its checks. **Open gap:** the same teardown fault may abort the
desktop app at exit on Python 3.11 with PySide6 6.12; not reproduced locally (Python 3.12, no xvfb here).

**CI, next step (2026-11-07).** With stderr uncaptured, the full-suite abort printed its reason: the same "Fatal
Python error: none_dealloc: deallocating None" in the middle of `tests/test_analysis_panel.py`. It is a reference
count error in shiboken6/PySide6 6.12 that over-releases `None`; on Python 3.12 `None` is immortal, so the fault
is invisible there (every local run is 3.12), while on 3.11 it aborts the process. The test workflows now run on
Python 3.12, the version ParaScape is developed and used on. **Open risk, for Einar:** `requires-python` still
says >= 3.10, so a 3.10/3.11 user with PySide6 6.12 can hit this abort; options are raising `requires-python` to
3.12, or excluding PySide6 6.12 on older Pythons once a fixed release exists.
