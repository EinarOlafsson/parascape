# 64.23 · Select and install available organism spaces

Status: DONE — 100%, 2026-11-07. Space catalogue with explicit ready/installed/installable/unavailable states; host packs listed and installable from the Species menu, never shown as calibrated. Built; Einar's review pending -- reopen on his feedback.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.18, 64.21, 64.22 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `parascape/organisms.py`, `parascape/packs.py`, `parascape/paths.py`, `parascape/app.py` |

## Deliverable

A space selector and pack catalogue for parasites and hosts, with explicit available/installable/unavailable states.

## Controlled scope

Integrate the existing two parasite spaces and validated Hs/Mm packs; additional species get separate builder tasks.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Selecting a space updates evidence, targets, strategies, scorecards and caches atomically without retaining another species' context.
- [x] Published pack versions/hashes, interrupted install and offline operation follow the existing pack contract.
- [x] Other instruction-53 species use the same admission checklist one species at a time; an unfinished space is never shown as calibrated or ready.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

2026-11-07. `parascape/space_catalogue.py`: one row per declared space -- Toxoplasma and Plasmodium maps
(shipped) and the human and mouse gene packs (64.21, 64.22) -- with one explicit state: ready, installed (version),
installable (a local archive with its published SHA256) or unavailable (and why). `calibrated` is true only for a
space with a held-out track record, so an installed host pack is shown "not calibrated" and never as ready for
inference. `install(code)` installs through the existing pack contract (published SHA256 required, verified
before activation, previous version kept if activation fails, no network).

The Species menu lists the host spaces under the parasite maps with that state; an installable one installs on
click. Selecting a parasite map still opens a fresh window of its own (`open_species`), so no other species
context is retained; an installed host space is read through a parasite gene's links (Gene module, 64.24).

Checks: `tests/test_space_catalogue.py` (3 pass: every space has one state; a host pack installs with its
checksum and stays uncalibrated; a pack without its published checksum, or with a wrong one, is refused),
`tests/test_module_shell.py::test_the_species_menu_lists_host_spaces_with_their_state`, and the existing pack
contract tests in `tests/test_packs.py` (interrupted activation, tampering, offline, size limits).

Gaps, stated: a host space does not open as a full map window (the window is built for parasite tables; a host
map window would be a separate item); packs are found in this checkout, not yet hosted for download; other
instruction-53 species have no builder yet and are not listed.
