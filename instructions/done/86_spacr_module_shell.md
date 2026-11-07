# 86 · Every module in spaCR's module shell

**Status: built 2026-11-06, waiting on Einar.** Every available tile opens a module in the shell:
Map, Gene, Network (and Partners = Network by binding evidence), Infer (and Questions) built with the
parascape-module mould; Expression, Screen, Cell, Gallery, Strategies, Scorecard, Information, Drop, Jobs and
Console as tool modules (`parascape/dock_module.py`: the panel lent into the shell while on screen, with the
shared masthead, Your gene, actions and folded console / AI chat). Datasets and Export open their dialogs;
Structure, Runs, Fill and List stay greyed with their reasons. Gate test:
`test_modules_in_window.py::test_every_tile_with_a_panel_opens_a_module_in_the_same_shell`. Later
refinements (a mould redesign per tool module) follow Einar's review.
- Step 1 done: spaCR's `collapsible_splitter`, `collapsible_section`, `foldable` and `collapse_arrow` are
  vendored in `parascape/panes/` (imports pointed at `panes/_compat.py`); `parascape/shell.py` is the shell:
  masthead with Home, Settings column that folds to its edge, runtime column of sections that fold to their
  headings, Actions row, and the window's Console and AI chat borrowed side by side while a module is up.
- Step 2 done for Map: `parascape/map_module.py`. The analysis panel's "2 · Map" controls are MOVED into it
  (same widgets, so `spec()`, `run_embed`, `run_umap_walk` are unchanged), plus a live 2D preview of a gene
  sample that redraws as n_neighbors / min_dist / seed change. The Map tile opens it; "show the 3D map" leaves
  to the 3D view; the old tab has a button to the module.
- Lesson: no pyqtgraph PlotWidget inside a shell -- deleting one mid-session segfaulted in
  `QGraphicsScene::~QGraphicsScene`; the preview paints its own points.
- 2026-11-06: Map rebuilt to the agreed mould (and Clusters moved into it as Groups); Gene, Network and
  **Infer** (the analysis Questions tab moved in; Questions and Infer tiles both open it; your genes read
  against the answer) built with the `parascape-module` skill. 5 of 22 modules on the shell.
- Next, in order: Strategies, Scorecard, Expression, Partners, Cell, Datasets, Information, Drop, Screen,
  Export, Gallery, Jobs, Console.

Requested 2026-11-05 (instruction 74, goals 6-7). Einar: "spaCR has its settings in a
container on the left, a console, AI chatbox, and an actions section. There is a live section to live test
parameters and a hyperparameter search. I'd like you to copy this structure whenever it makes sense and when
it doesn't build something that does ... spaCR containers are also collapsible so after setting the
settings you can give more space to the content."

## The plan
1. Port spaCR's `AppScreen` layout (`spacr/qt/screens/app_screen.py`): settings column left; console, AI
   chat and actions right; collapsible sections and splitter (`collapsible_section.py`,
   `collapsible_splitter.py`).
2. Live preview and hyperparameter search where a module has parameters worth sweeping -- in ParaScape that is
   inference runs and map building, not segmentation.
3. Migrate modules onto the shell one per change, each independently testable; nothing deleted until its
   replacement works.

## Gate
Each migrated module is reachable from its tile, keeps every control it had, and collapses its settings.
