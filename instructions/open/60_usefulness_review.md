# 60 · What would make ParaScape most useful (review at 0.49.0)

**Status: open (reviewed 2026-09-30).** A read-only review of the whole program against the way a
working biologist uses it, asked for by the user: "think about how to make this tool as useful as
possible". Weeks 1-2 below are the agreed order of work.

## What is strong, and should be built on

* **The calibration layer**, and nothing else in this field has it: 7,640 self-test runs per
  strategy × organism × setting × held-out label × seed, against a *measured* null, with intervals
  and a grade. `docs/calibration.md` says in its own preamble that strategy 18's AUROC near 0.97 is
  an artefact of its null and that the honest figure is about 0.67. That paragraph is worth more
  than most of the UI.
* **The leakage closure** (`Context.banned`, `banned_layers`, block-level dropping in
  `Context.blocks`), which carries the record of the bug that forced it.
* **Provenance in one place**: 162 `Dataset` entries with PMID, accession, citation and columns;
  `datasets.provenance(column)`, `derived_sources`, `unresolved()`. **The single biggest unexploited
  capability in the tree.**
* **Scorecards**: 6 tasks, 51 metrics, each with its chance level and how to read it.
* **The 0.49.0 trio**: `guided.py` is Qt-free and ranks by the grade measured *for that organism*,
  refuses strategies the table cannot fill, and says so when nothing better than weak exists.
* **Absence is handled honestly** in the evidence panel ("unknown, not absent"; "absence of
  attention, not absence of function").

## Where a biologist still hits friction

* **The 0.48.0 dock/card bug** breaks the first thing anyone does. Everything else is behind it.
* **Eleven docks**, 39 strategies, 290 slots, 40 techniques. "Start here" fixed the entry problem
  but *added* a surface; nothing yet removes one.
* **A gene's page is a hand-picked 25 fields of 443 columns.** There is no answer to "of everything
  measured on this gene, what is unusual?", and no sight of what the shipped strategy runs already
  said about it — which the star-edge cache holds.
* **Your own screen never reaches the strategies.** `File > Import data` updates only the evidence
  panel, while `StrategyPanel` built its `Context` once at construction and the guided panel reuses
  that same object, so an imported column cannot be used in Strategies or Start here in that
  session. The importer is also Toxoplasma-only, and an imported column has no leakage family, so it
  would be unguarded if it did arrive.
* **Nothing can leave the program.** Export is CSV, GraphML and a screen-sized PNG. The one
  publication-grade artefact, `report.recipe_pdf`, is wired to exactly one sub-tab. There is no
  methods sentence and no citation list anywhere, though every citation a run needs is already
  resolvable from the columns it touched.
* **A user's own run is not reproducible**: `StrategyResult.save` records tables and settings but
  not the organism, the version, the table fingerprint, the seed or the resolved columns.
  `searches.py` already implements exactly the right manifest; the newer panels do not use it.
* **Finding versus artefact**: nothing reports when a result rides on median-filled missing values
  or on one dominant dataset, and the `all_*` gallery maps leave 51-62% of genes unclustered without
  saying so. `per_target` skill is in the shipped calibration for every strategy and is shown
  nowhere, so "this label is unlearnable by this strategy" is averaged away.

## The ten changes, ranked by usefulness per unit of effort

| # | Change | Why | Effort | Main files |
|---|---|---|---|---|
| 1 | Fix the dock/card bug on a real display, with a regression test | Everything is behind it | ½ day | `strategy_panel.py`, `strategy_card.py`, `glass.py`, `app.py` |
| 2 | **Methods paragraph + resolved citation list + BibTeX** for any run | The step that costs a day per figure and is where wrong citations enter | 2-3 days | new `methods_text.py`, `datasets.py`, `strategies.py` |
| 3 | Record provenance on results (organism, version, table fingerprint, seed, columns); session save/load | Prerequisite for 2, 4, 6, 10 | 2 days | `strategies.py`, `results.py`, `searches.py`, `app.py` |
| 4 | Figure export (PDF/SVG) from every panel via `report.py` | A screen-sized PNG is not a figure | 2 days | `report.py`, the three panels |
| 5 | **A gene digest**: coverage, the 10 most extreme measurements as percentiles, which strategies named it, which maps cluster it | A biologist's first question is currently unanswerable without four docks | 3 days | `app.py`, `star_edges.py`, `umap_gallery.py` |
| 6 | Surface `per_target` calibration ("what does the calibration say about MY label") and weight the guided ranking by it | Pure UI over data already shipped | 1-2 days | `guided.py`, `strategy_panel.py` |
| 7 | Missingness / dominance / noise warning on every result | The difference between a result and an artefact | 2 days | new `diagnostics.py`, `strategies.py` |
| 8 | Make "drop in my own screen" work end to end (rebuild contexts, Pf importer, require a leakage family, place it immediately) | The entry point the user asked for | 4 days | `app.py`, `importer.py`, `search.py`, `guided.py` |
| 9 | Many more gallery maps, kept by structure score rather than all | The 51-62% noise maps are the weakest visible artefact | 3 days + compute | `umap_gallery.py`, its build script |
| 10 | Generated screen-capture tutorial videos | They stay current because the script regenerates them | 2-3 days after 2 and 4 | `scripts/build_tutorials.py` |

## What NOT to build

* **A twelfth dock, or strategies 40+.** 39 strategies, 9 weak and 2 untestable on Pf, is already
  more than the calibration can vouch for. Hiding the "no skill" ones would add more than adding
  another family.
* **The assistant as a route to answers.** An LLM narrating over data whose whole value is the
  measured chance level is the one component that can manufacture the unvalidated claim this project
  exists to prevent.
* **More on-screen rendering realism.** Figure-quality *export* pays; PBR does not.
* **Cross-space bridges (WP10-12) before items 2, 3 and 7**, which would multiply the ways a leak
  can enter while a single-space result still has no provenance record.
* **Active learning**, already measured as no better than random (record 54).
* **A plugin/scripting layer**: `docs/API.md` plus notebooks already serve that.

## Order of work

| Week | Work | Lands in |
|---|---|---|
| 1 | Item 1, then item 3 as the substrate | 0.49.1, then 0.50.0 |
| 2 | Items 2 and 4 — "results that can leave the program" | 0.50.0 |
| 3 | Items 6, 7, 5 — "tell a finding from an artefact without the docs" | 0.51.0 |
| 4 | Items 8 and 9 — "start from the data you have" | 0.52.0 |
| tail | Item 10, then a fresh calibration sweep for anything whose inputs changed, then the organism spaces (record 53) | 0.53.0 |

Two standing conditions: any change to data, slots or strategies needs a fresh calibration sweep and
a tutorial rebuild before release, and no identifier in a generated citation list may be typed from
memory — the generator reads `datasets.py` and flags `unresolved()` rather than filling it.

## Change 2 built (2026-11-06): methods paragraph, references and BibTeX

`parascape/methods_text.py`: `sources(columns, organism)` resolves every column to its registry dataset and
follows declared derivations to the measurements behind them; `methods_paragraph`, `reference_list`, `bibtex`
and `full_text` write the text. **No identifier is typed from memory**: every citation is the registry entry
of a column used; a dataset with no recorded publication is listed as "not yet confirmed -- confirm before
citing" and kept out of the BibTeX; a layer ParaScape computed cites its inputs instead; a column with no
source is named. The BibTeX holds only the fields the registry records (the citation as recorded, PMID,
accession, URL) -- never a title, journal or DOI that would have to be guessed.

In the application: Infer ▸ Save the answer ▸ **Copy the methods text** -- for the selected question, from the
columns its closure allows (rendered and checked: 95 columns cited for the localisation question). Tests:
`tests/test_methods_text.py` (6), `tests/test_module_shell.py` (the button and the clipboard).

The registry has 103 entries without a confirmed publication; every methods text shows which of them it
touches. Confirming them is curation, and the text makes that list visible per run.

## Change 3 built (2026-11-06): provenance on every strategy result

`StrategyResult.provenance`, filled by `Strategy.run` through `strategies.run_provenance`: organism, ParaScape
version, the table it ran on (`searches.fingerprint`: gene identifiers hashed in order, gene and column
counts), the seed and the time; `save` writes it into the summary JSON. Test in `tests/test_strategies.py`.
Not done here: session save/load, and the resolved column list per strategy (strategies do not yet report
which columns they read; the methods text of change 2 uses the question closure instead).

## Change 6, where it stands (2026-11-06)

Half delivered by instruction 62: the guided recommendations (`guided.py`) append the measured hold-out rate for
the strategy on the label in hand (`track_record.record_phrase`), and the Infer module lists the best
recovering strategies for a question label. Still open: weighting the guided ranking itself by the per-label
calibration (`per_target`, shipped in `strategy_calibration.json` for every strategy, read nowhere yet).

**Done the same day:** `guided.label_weight` reads the per-label cell when the person chose a label: skill whose
lower bound is above zero adds up to 1.5 to the ranking score, measured no-skill on that label (upper bound at
or below zero) costs 2.0, and a label never held out moves nothing. Test in `tests/test_guided.py`.

## Change 4 built (2026-11-06): figure export from the module views

`parascape/figure_export.py`: `save(widget, path)` writes SVG or PDF as **vectors** -- the view is drawn again
onto the export device, so points and lines stay points and lines (checked: no embedded image in the Map or
Network SVG) -- or a PNG at three times the screen. A view may offer `paint_figure(painter)`; the Network view
does, because on screen its background is a cached image, and the figure draws it as vectors instead. Map and
Network each gain **Save the figure…** (SVG / PDF / PNG). Tests: `tests/test_figure_export.py` (5). Rendered and
checked on the real window: Map 1.9 MB SVG, Network 3.9 MB SVG with every link and gene.

## Change 7 built (2026-11-06): missingness and dominance warnings

`parascape/diagnostics.py`: `missingness(nodes, columns, genes)` compares the share of evidence missing for
the genes a result names with genes in general on the same columns, and warns at a gap of 15 points or more
("their places rest on filled-in values more than on measurements"); `dominance(importance)` warns when one
evidence family carries half or more of the positive importance. The Infer answer note shows the missingness
caution under the named genes; every harness result has `.dominance`. Tests: `tests/test_diagnostics.py` (6).
Not done: a noise check, and the warnings on the strategy cards (strategies do not report their columns).

## Change 5, the measurable half (2026-11-06): what is unusual about a gene

Gene module, new section **What is unusual about it** (`gene_reading.unusual` / `unusual_sentences`): of every
continuous column measured on at least 500 genes, the ten furthest from the middle, each as an exact rank
("the highest of the 8,131 genes measured") with its source, at most two per dataset, identity and attention
columns left out, filtered by "What to read". Checked on BAG1: the highest of every bradyzoite expression
measurement -- the textbook answer. Tests in `tests/test_gene_reading.py`; rendered on the real window.
Not done: "which strategies named it, which maps cluster it" (the star-edge cache), the other half of change 5.

## Change 8, the core (2026-11-07): an imported screen reaches the strategies

`StrategyPanel.set_nodes` rebuilds the strategy context over the table with the imported columns (same rows; the
graph is kept, it is checked against gene ids), and `GuidedPanel.set_context` hands the guided panel the same
context; `Window.apply_import` calls both. Before this an imported column reached the evidence panel and the
map, never a strategy. Test: `tests/test_app_controls.py`. Still open from change 8: a Plasmodium importer, and a
required leakage family for an imported column (today the measured-association closure is what guards it).
