---
name: parascape-module
description: Build or rebuild a ParaScape module (repo/starplast) in the agreed mould - spaCR module shell, goal-first layout (your genes -> what surrounds them -> how far to trust it -> next step to inference), plain wording, tests and render check. Use whenever a ParaScape module, tile, dock panel or Navigator entry is being created, moved onto the shell, or redesigned.
---

# The ParaScape module mould

Agreed with Einar on 2026-11-06 on the Map module (`parascape/map_module.py`, instruction 96).
Map is the reference implementation: when in doubt, do what it does.

## The two goals every module serves

1. Make information about **a gene or a list of genes** accessible and **synthesised**.
2. Make **inference** possible from labels or lists, to create knowledge.

Let these decide the layout and wording. Do NOT write them on screen. A module shows the
method only where it changes the answer, and never before the thing the user came for.

Before laying out anything, write one sentence: "This module lets you take <gene/list> and
see <synthesis>, so that you can <inference / next step>." If you cannot, ask Einar
(AskUserQuestion) before building.

## Layout (`parascape/shell.py` ModuleShell)

Masthead: the module name, a short description of what it IS ("every gene, placed by how alike
its measurements are"), and an instruction of what to DO ("Put your genes on the map and read
what surrounds them."). Home button on the right.

**Left column (settings), in this order:**
1. **Your genes** - the subject. A field (gene id, product word, several separated by commas;
   Return applies), "Choose a list…" (file dialog, `gene_drop.read_drop`), "Use the gene list"
   (`window.home.genes()`), "Clear", and a line saying what is selected. A file dropped on the
   module is read here (event filter on the shell, not the window-wide drop). Clicking the main
   view also selects.
2. **The lens** - what label the result is read by ("Colour by"), synced both ways with the
   window's own control (one source of truth).
3. **The choice that changes the meaning**, in plain words (for Map: "What makes genes similar"
   = top-level measurement categories as tick boxes driving the existing block tree).
4. **Advanced** (`expanded=False`) - method parameters, plain name first with the technical name
   in parentheses ("neighbours considered (n_neighbors)"), plus a link to the full detail panel.
5. **Save** (`expanded=False`) when the result can be kept.

**Right column (runtime sections), in this order:**
1. **The thing itself**, large (`stretch=3`, focus): the real current data, never a sample and
   never a button away. A one-line note says what is shown ("The current map · 8,140 genes ·
   coloured by compartment"); a legend; the selection marked. A preview of changed settings
   replaces it with a clear "Preview … Build to apply" note and a "Back to the current …" button.
2. **Your genes and <their reading>** - the synthesis, as sentences: what surrounds them (counts
   per label), the gene's own label, and when unknown, what the data suggest.
3. **Trust, on the same data** - every suggestion carries how often that kind of call is right
   on this data, AGAINST CHANCE (commonest label), with n. Warn when the label being read was
   one of the inputs (circularity). Pure functions in a testable module (`map_reading.py`).
4. **Find the best … for a label** (`folded=True`) - the optimiser, scored by the GOAL metric
   (agreement on the chosen label), not a method metric; ranked table; click a row to apply.
5. **Actions**: primary "Build/Run this" first; then the next step toward inference ("Predict
   from neighbours" -> Infer); then "Full … view", "Export". Stop button from the panel.
6. **Console and AI chat** - `shell.borrow(..., folded=True)`.

## Navigator and catalogue

- `parascape/modules.py`: the module entry (note = what it is, in plain words), `opens`, and a
  `DOWNSTREAM` entry ordered by the goals: synthesis next steps first, then the step to inference.
- `parascape/navigator.py` `_MODULE_ICON`: a spaCR icon from `data/icons/spacr/`.
- `landing_install.install`: build the module after every dock, `landing.stack.addWidget(shell)`
  at once (never leave it top-level), and route `open_module(window, key)` to it.

## Mechanics that have already bitten

- **Move, do not copy** existing controls: reparent the panel's own widgets so its methods keep
  reading them; replace the old tab's content with a note and an "Open the … module" button
  (with a tooltip).
- **No pyqtgraph PlotWidget and no second GL view** inside a module: draw in software
  (`map_module.ScatterView`: drag to turn, wheel to zoom, click to pick). A PlotWidget deleted
  mid-session segfaulted in `QGraphicsScene::~QGraphicsScene`.
- Worker threads: unparented `QThread` held in a module-level set until `finished`; stop timers
  in `hideEvent`; `aboutToQuit` waits for them.
- PySide6 ignores virtual overrides assigned on an instance: use an event filter or a subclass.
- Shell QSS rules are prefixed `QWidget#ModuleShell`; containers use the pane opacity
  (`backdrops._compat.get_pane_opacity`, default 60 %); fields fade (scoped by
  `parascapeFieldFade`). Sections that support the work start folded.
- Docstring on every public method; tooltip on every control (`tests/test_app_menus.py` checks).
- Foldable groups created inside a module go through `shell.section(...)` (heading in the theme
  colour; a bare CollapsibleSection heading is near-invisible). Every QTextBrowser with links gets
  `shell.link_style(browser)` (Qt's default link blue is unreadable on the dark containers).
- Facts about a gene: `gene_reading.read(published_map, gene, node_row=...)` -- the published
  information map (`window.information_map.map`, load it with `.load()` and wait) plus the
  window's gene-table row; categories from `information_categories` so a sentence and a
  circularity exclusion use the same words. A prediction is never phrased as knowledge.
- Second reference implementation: `parascape/gene_module.py` (one gene / stepping a list).

## Wording

Plain words over jargon (say what a term does; keep the identifier in parentheses at most).
No possessive 's (rephrase with "of"). Percentages as whole numbers for rates of being right.
Say "how often right, against chance" rather than a metric name. Absence is "unknown", grey.

## Tests and verification (all required before pushing)

- A pure-function test file for the reading/trust logic (synthetic data, chance baseline,
  unplaced/unknown handling, too-little-data returns None).
- Window tests (`tests/test_module_shell.py` pattern): sections present and in order; Advanced
  starts folded; console folded; the main view shows the real current data; a typed gene and a
  list file are selected and read; the plain-language choice drives the underlying control; the
  lens is synced with the window; controls moved not copied; the old tab points to the module.
- Run: the module tests + `test_modules_in_window`, `test_navigator`, `test_help_search`,
  `test_app_menus`, `test_app_smoke`, `test_docstrings`, `test_publishing`, inside
  `systemd-run --user --scope -p MemoryMax=14G nice -n 19 taskset -c 0-23`, `-p no:randomly`.
- Render the REAL `app.Window` offscreen with `XDG_CONFIG_HOME` set to a scratch directory (never
  touch Einar's `~/.config/ParaScape`), open the module, select a gene, `w.grab()`, and LOOK.
  Offscreen GL "Error while drawing item" tracebacks are an artifact; ignore them.

## Finishing

Update the module's instruction file in `instructions/open/` and its INDEX row (status, what was
built, tests). Commit by explicit paths as Einar Olafsson <einar.olafsson@gmail.com> with NO
co-author trailer; push `HEAD:nightly`. Then show Einar a render and say in plain words what the
module is for, and ask what would make it more intuitive.
