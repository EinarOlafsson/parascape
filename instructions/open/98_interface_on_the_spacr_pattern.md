# 98 · Restructure the ParaScape interface on the spaCR pattern

**Status: OPEN. Requested by Einar 2026-11-12, after installing 0.0.11.0.**

> "The new interface is really close to what I wanted in terms of the navigator screen. But everything behind
> it seems to be unchanged. So please take stock of each component and part of ParaScape then rearrange
> ParaScape so that the software and each module's structure resembles the structure in spaCR. ParaScape is
> already a powerful tool. Now we just need to restructure the user interface, so it is more intuitive and
> easier to use."

## What is actually being asked

Not a reskin. The navigator landing page already looks right; **everything behind it is the old application**.
Typing a gene and clicking the module it offers drops the user into a pre-navigator window. The request is to
make the structure underneath match the structure spaCR has, because that structure works.

Einar named these, and each is a requirement rather than a suggestion:

1. **A search bar in the navigator, like spaCR's.**
2. **No small "×" close affordance. Anywhere.** "I don't want to see that anywhere."
3. **The same symbols as spaCR**, and the **AI symbol behaving the same way**.
4. **Recipes** exist for ParaScape and work as they do in spaCR.
5. **A tree of the things ParaScape can do**, whose progression guides the user to the module they want.
6. **Type one gene → a gene card module** (overhauled to the module mould, carrying the information it carries
   today, *including the parasite schematic*). **Click it once → the modules downstream of one gene.**
7. **Type a list → a gene list module**, and clicking it reveals the modules downstream of a list.
8. **Every module and every analysis thematically consistent**, and intuitive to reach.

## The four things spaCR does that ParaScape does not, with where to read them

Study these before writing anything. They are the whole pattern.

| spaCR | what it gives | ParaScape has |
|---|---|---|
| `spacr/ports.py` | Each module declares what it **consumes** and **produces**, as typed ports with a *kind*. `next_modules`, `upstream_modules`, `producers_of`, `consumers_of`, and `check_ready` which answers **before a run starts** whether a module can run and, when not, why and what to fix. | nothing equivalent |
| `spacr/qt/app_catalog.py` | `DECLARED_APPS`: ONE canonical table of key, name, description, **section**, factory, stage, intro. Screens read their metadata from it instead of repeating strings. | `parascape/modules.py`, partially |
| `spacr/qt/app.py` `SECTION_*`, `section_members`, `home_categories`, `home_bands` | Home and every tab are **computed from that one table**, so "a section cannot acquire a tab it has no apps for or lose one it does". A documented bug (#16i) was Home grouping by one table and the tabs by another. | Home is partly hand-written |
| `spacr/qt/chaining.py` + `spacr/ports.py` | The strip above Run: what it will read, what is stale, **what comes next** — successors offered pre-filled, each checked, an offer that cannot work shown greyed **with the reason**. | nothing |

Also read `spacr/qt/recipes.py` (format, storage, module-scoping, version behaviour), `spacr/qt/iconset.py`
(the symbols), and `spacr/qt/command_palette.py` with `spacr/qt/help_search.py` (search).

**The translation to ParaScape.** spaCR's port *kind* is a folder of data on disk. ParaScape's equivalent is
the **subject** a module acts on: one gene, a gene list, a pair of genes, a space, a screen, a trait, a
compound. That substitution is the whole design. "The modules downstream of one gene" is exactly
`consumers_of("gene")`, and "downstream of a list" is `consumers_of("gene_list")`.

## Work packages

**WP1 — Take stock.** An inventory of every ParaScape component: each module, panel, tab and window; what it
takes as input, what it shows, what it produces, whether it is reachable from the navigator, and whether it
has been through the module mould. Record it as a table in this file. **Nothing else starts until this
exists**, because "rearrange ParaScape" cannot be scoped without it.

**WP2 — `parascape/ports.py`: typed subjects.** Every module declares the subject it consumes and what it
produces. `consumers_of`, `producers_of`, `next_modules`, `upstream_modules`, `check_ready(module, subject)`
returning a readiness with reasons. Pure, no Qt, fully tested.

**WP3 — One catalogue.** A single declared table carrying key, name, plain-words description, section, icon,
factory and the ports from WP2. Navigator, tree and every menu computed from it. A test that no section is
empty and that nothing appears in two sections — spaCR has `test_no_section_is_empty` by name.

**WP4 — Navigator search.** A search bar as spaCR's works. Typing a gene identifier, a product word, or a
comma-separated list resolves the **subject**, then offers the modules that consume that subject, ranked.

**WP5 — Symbols.** The spaCR icon set throughout (`parascape/data/icons/spacr/`), one icon per module declared
in the catalogue. The AI symbol behaves as spaCR's does. **A test that no "×" close affordance exists in any
ParaScape widget**, because Einar asked for it by name and it will otherwise creep back.

**WP6 — Recipes.** Named, reusable, shareable settings bundles, per module, carrying the ParaScape version
they were captured with, refusing to apply across modules, one JSON file each. Same shape as
`~/.spacr/recipes/<module>/<slug>.json`.

**WP7 — The tree.** A visible tree of what ParaScape can do, whose progression leads to a module. Computed
from WP2 and WP3, never written down twice.

**WP8 — The gene card.** Rebuilt on the module mould (`.claude/skills/parascape-module`), carrying everything
it carries today **including the parasite schematic**. Clicking it once reveals the single-gene modules
downstream, as spaCR's chaining offers successors.

**WP9 — The gene list module.** The same, for a list.

**WP10 — Every remaining module to the mould**, and the pre-navigator windows retired rather than left
reachable. This is the bulk of the work and the reason the request was made: the navigator is right and what
is behind it is not.

**WP11 — Consistency as tests.** Every module: declares its ports, has an icon, appears in exactly one
section, is reachable from the navigator, has a recipe surface, carries a docstring and tooltips, and shows no
"×". A module that fails any of these fails the suite rather than being noticed later by Einar.

## Rules that apply to all of it

* The module mould in `.claude/skills/parascape-module` governs every module: real data in the main view,
  trust reported against chance, the two goals, Advanced folded, console folded.
* **Move controls, do not copy them** — the mould records that reparenting is what keeps the old panel's
  methods working, and that a deleted `PlotWidget` segfaulted.
* Plain words over jargon; no possessive 's; one decimal on percentages.
* Verify by rendering the REAL `app.Window` offscreen with `XDG_CONFIG_HOME` pointed at a scratch directory,
  and LOOK at the result. Never Einar's `~/.config/ParaScape`.
* Commit by explicit paths as Einar Olafsson <einar.olafsson@gmail.com>, no co-author trailer, push
  `HEAD:nightly`.

## What must not be lost

ParaScape is already powerful and this is a restructuring, not a rewrite. No analysis, column, calibration or
frozen result may be dropped in the move. In particular the Plasmodium node table is pinned to the frozen
functional benchmarks (`deposits.HELD_OUT_OF_NODES`), and the gene card's parasite schematic was named by
Einar as something to keep.

## WP1 inventory

Taken 2026-11-12 against ParaScape 0.0.11.0 on `nightly`, by opening every file rather than reading
the names. Nothing here is guessed: where a module says what it takes, that is what its code reads.

**The count.** 22 modules, 15 docks, 29 tabs inside five of those docks, and the main window
plus 9 further windows and dialogs.
Of the 22 modules, **4 have been through the module mould** (`.claude/skills/parascape-module`):
Map, Gene, Network and Infer, each with its own file and its own tested reading functions. Two more
keys open one of those four shells (Partners is Network read by binding evidence; Questions is Infer
opened on its catalogue), so **6 of 22 keys reach a mould-built module**. **10 more sit in the
spaCR shell but have not been laid out to the mould** — they borrow an existing dock panel whole
(`dock_module.py`), so they look alike and are not yet goal-first inside. **6 have neither**: four
are honestly greyed out because the view does not exist (Structure, Runs, Fill, List) and two are
window actions rather than modules (Datasets, Export).

**What "takes" and "hands on" mean.** The subject, in the sense of WP2 — `parascape/ports.py`.
A `?` marks a subject the module can work without.

### The modules

#### Look up

| module | takes | hands on | opens | mould | in navigator |
|---|---|---|---|---|---|
| **Gene** (`gene`) | gene, gene_list? | gene?, gene_pair?, claim? | module shell | yes — `gene_module.py` | yes |
| **Structure** (`structure`) | gene | — | nothing built | no | greyed: no structure view; AF3 summaries are on the evidence of a gene for now |
| **Partners** (`partners`) | gene, gene_list? | gene_pair, gene_list? | module shell | yes — `network_module.py` | yes |
| **Expression** (`expression`) | gene, gene_list? | trait? | shell over the `compare` dock, tab "Gene expression" | shell only | yes |

#### Explore

| module | takes | hands on | opens | mould | in navigator |
|---|---|---|---|---|---|
| **Map** (`map`) | gene_list?, trait? | space, gene_list?, gene?, job? | module shell | yes — `map_module.py` | yes |
| **Network** (`network`) | gene?, gene_list?, space? | gene_pair, gene_list?, gene? | module shell | yes — `network_module.py` | yes |
| **Cell** (`cell`) | gene?, trait? | gene_list?, gene? | shell over the `find` dock | shell only | yes |
| **Gallery** (`gallery`) | space | space | shell over the `gallery` dock | shell only | yes |

#### Ask

| module | takes | hands on | opens | mould | in navigator |
|---|---|---|---|---|---|
| **Questions** (`questions`) | — | question | module shell | yes — `infer_module.py` | yes |
| **Infer** (`infer`) | question, gene_list?, trait? | claim, gene_list?, run, job? | module shell | yes — `infer_module.py` | yes |
| **Strategies** (`strategies`) | trait?, gene_list? | claim?, scorecard, job? | shell over the `strategies` dock | shell only | yes |
| **Scorecard** (`scorecard`) | trait?, claim? | scorecard | shell over the `discoveries` dock | shell only | yes |
| **Runs** (`runs`) | run | scorecard? | nothing built | no | greyed: no panel reads `results/harness_holdout_log.jsonl` |
| **Fill** (`fill`) | claim, scorecard | trait | nothing built | no | greyed: no question has cleared its 90 % bar |

#### Bring your own

| module | takes | hands on | opens | mould | in navigator |
|---|---|---|---|---|---|
| **Drop** (`drop`) | — | gene_list, screen? | shell over the `gene list` dock | shell only | yes |
| **Screen** (`screen`) | screen | gene_list, trait? | shell over the `compare` dock | shell only | yes |
| **List** (`list`) | gene_list | trait? | nothing built | no | greyed: no set-level enrichment view |
| **Export** (`export`) | gene_list?, space?, trait?, claim? | — | window action `export_visible` | no | yes |

#### Data

| module | takes | hands on | opens | mould | in navigator |
|---|---|---|---|---|---|
| **Datasets** (`datasets`) | — | dataset | window action `open_dataset_browser` | no | yes |
| **Information** (`information`) | gene?, trait?, dataset? | trait, claim?, gene? | shell over the `information map` dock | shell only | yes |
| **Jobs** (`jobs`) | job | — | shell over the `jobs` dock | shell only | yes |
| **Console** (`console`) | job | — | shell over the `console` dock | shell only | yes |

### The docks

Fifteen, all built by `app.Window`. The landing covers them rather than hiding them, so every one of
these is still reachable behind it. "Lent" means a module borrows the panel while it is on screen
and gives it back (`shell._Borrowed`).

| dock | built by | what it shows | lent to | reachable from the navigator |
|---|---|---|---|---|
| `find` | `app._left` | the search field, the colour-by control and the parasite drawing | Cell | yes, as Cell |
| `evidence` | `app._right` | one gene read as text, with links into its sources | — | yes, as Gene and Partners (which open the shells instead) |
| `analysis` | `app._analysis` + `analysis_panel.py` | eight tabs: data, map, clusters, inference, search, validation, discover, questions | Infer (tab 8 only) | tab 8 only, as Questions |
| `strategies` | `strategy_panel.py` | 39 named ways of inferring something, each graded and self-testing | Strategies | yes |
| `maps` | `umap_gallery_panel.py` | the pregenerated maps and their label scores | — | **no** |
| `star map` | `star_map.py` | one gene and its links, by source | — | the Network module replaced it; the dock is still installed |
| `discoveries` | `discoveries_panel.py` | five tabs: labels, inferred claims, functional tests, strategy coverage, host sources | Scorecard | yes, as Scorecard |
| `compare` | `compare_panel.py` | eleven tabs of measured values and study design | Expression (tab 1), Screen | yes, as Expression and Screen |
| `information map` | `information_map_panel.py` | two tabs: one gene, one trait — measured facts with source and certainty | Information | yes |
| `gene list` | `home_panel.py` | the gene list a dropped file becomes, and the question catalogue | Drop | yes, as Drop |
| `start here` | `guided_panel.py` | one question at a time, ending at the recommended strategies | — | **no** |
| `gallery` | `gallery.py` | a wall of thumbnails, one per map in a walk | Gallery | yes |
| `console` | `console.py` | what the program is doing, filterable | Console (and folded into every other module) | yes |
| `jobs` | `app._tool_docks` | what is running, what it costs, how to stop it | Jobs | yes |
| `assistant` | `chat.py` | the AI chat | folded into every module | the Ask button on the navigator |

### The tabs inside those docks

29 in all, across five of the fifteen docks; the other ten hold one panel each. These are where
most of the application actually lives, and all but three of them are
reachable only by finding the right dock and then the right tab — which is the "everything behind
the navigator is unchanged" Einar reported.

| dock | tabs | reachable by name |
|---|---|---|
| `analysis` | 1 · Data, 2 · Map, 3 · Clusters, 4 · Inference, 5 · Search, 6 · Validation, 7 · Discover, 8 · Questions | only "8 · Questions", as the Questions module |
| `compare` | Gene expression, Studies & samples, Coexpression, Published expression, Deposited matrices, Other measurements, Sparse measurements, Original workbooks, Reference identifiers, Dataset coverage, Condition comparison | only "Gene expression", as the Expression module |
| `discoveries` | Labels and functions, Inferred claims, Functional tests, Strategy coverage, Host source status | no |
| `information map` | One gene, One trait | no |
| `strategies` | Guide, Settings, Results | no |

Two tabs of the eight in `analysis` have been MOVED into modules rather than copied: "2 · Map"
into the Map module and "8 · Questions" into Infer, with the panel widgets reparented so
`AnalysisPanel.spec()`, `run_embed`, `run_umap_walk`, `check_question` and `run_question` still read
exactly what they read before.

### The windows and dialogs

| window | opened by | what it shows | reachable from the navigator |
|---|---|---|---|
| `app.Window` | the application | the 3D map as the central widget, every dock, the landing over all of it | it is the application |
| Dataset browser (`dataset_browser.py`) | `Window.open_dataset_browser` | all sources, their coverage and provenance | yes, as Datasets |
| Slot tree (`slot_tree.py`, a `QMainWindow`) | `Window.open_slot_tree` | two tabs: Slots, Datasets | **no** |
| Gene evidence browser (`gene_evidence_browser.py`) | `Window.open_gene_evidence` | one gene, every record behind it | **no** — reached by clicking a link in the evidence dock |
| Label space browser (`label_space_browser.py`) | `Window.open_label_browser` | two tabs: membership and sources, strategies and recorded tests | **no** |
| Workflows (`workflows.py`) | `Window._open_workflows` | three tabs: explore a gene, predict a trait, compare a screen | **no** — this is the pre-navigator front door and WP10 retires it |
| Question scoreboard (`question_scoreboard.py`) | `Window.open_question_scoreboard` | how each registered question has done | **no** |
| Preferences (`app.build_preferences`) | `Window.open_preferences` | three tabs: Appearance, Display, Window | the menus |
| Download dialogs (`data_download_dialog.py`, `component_download_dialog.py`) | the Data menu | what can be fetched and what is already here | the menus |

### What the inventory says about the rest of the work

1. **The modules are not the problem; the routes to them are.** Ten of the twenty-two already sit in
   the right shell, and four are properly to the mould. What is unreachable is the 29 tabs and
   five of the nine further windows, almost none of which the navigator can name.
2. **`Workflows` is the pre-navigator window Einar described.** Three tabs — explore a gene, predict
   a trait, compare a screen — that duplicate Gene, Infer and Screen. WP10 retires it rather than
   leaving it reachable from the menus.
3. **Two docks have no module at all**: `maps` (the pregenerated maps, which the Gallery module
   arguably covers) and `start here` (the guided panel). Both are candidates for the tree in WP7
   rather than for new tiles.
4. **`star map` is dead weight.** The Network module replaced what it did; the dock is still built
   and still installed. Deciding whether it goes is part of WP10, not of WP1.
5. **Nothing consumes a compound.** Drug sensitivity exists only as per-gene columns
   (`datasets.curated_drug_sensitivity`), so the subject is declared in `ports.py` and has no module
   to hand one to. That is recorded rather than quietly dropped, and
   `tests/test_ports.py::test_nothing_takes_a_compound_yet_and_that_is_recorded_not_hidden` says so.

### WP2 and WP3 as built

* **`parascape/ports.py`** declares 13 subjects — the seven Einar named, plus question, claim,
  scorecard, run, dataset and job, which the inventory above forced. It has `consumers_of`,
  `producers_of`, `next_modules`, `upstream_modules`, `route` and `check_ready`, is standard library
  only, and imports nothing else out of ParaScape. `tests/test_ports.py` has 50 tests including two
  that run a fresh interpreter to hold the purity line.
* **`parascape/modules.py`** is the one catalogue, extended rather than replaced: it gained
  `instruction`, `builder`, `section`, and `consumes`/`produces` read from `ports.py`.
  `dock_module.SPECS` and the shell routes in `landing_install` are now computed from it, so the
  four tables that described these modules are one. `tests/test_module_catalogue.py` has
  `test_no_section_is_empty` and `test_no_module_is_in_two_sections` by name, checks that every row
  has ports and that every reachable module in a real window is a row.

## WP5 status (symbols and the small cross) -- done 2026-11-12

* **One icon provider**: `parascape/iconset.py`. Every module tile, dock row and result card gets its
  symbol from it, by module key or kind of result, re-inked for the theme and cached. The artwork is
  byte-identical to spaCR (a test compares the files). A symbol that cannot be found raises
  `MissingIcon` under pytest and falls back to the settings symbol, with a logged warning, outside it.
  `home.png` was added for the Home row; the Home module tiles now carry their symbols (they had none).
* **The AI symbol**: `parascape/ai_button.py`, drawn as spaCR draws it (ring, provider-coloured dot only
  while on, pulse while typing, circling marker or glow while waiting, still when Reduce motion is on).
  Click switches it on or off; right-click, or the first click with no provider, chooses the provider;
  with it on, a search that finds nothing is answered in the panel below the field. The old Ask word
  that opened the assistant dock is gone, so the dock is reached from the Tools menu only.
* **No small cross**: removed the clear button from four search fields (Navigator, Search ParaScape,
  Home search, maps filter; Escape empties them), the closable flag from the console, jobs, assistant
  and gallery docks (the Tools menu shows and hides them), and the close mark from every glass window
  (Escape closes it). `parascape/no_cross.py` finds any cross by walking a real window; the test
  `tests/test_symbols_and_no_cross.py` runs it on the real `app.Window` and after opening every dock module.
* Left open: the native title bar of an undressed window is drawn by the window manager and has its own
  cross; it cannot be changed from inside the application.

## Status 2026-11-12: WP1, WP2, WP3, WP5 and WP6 built

| WP | state |
|---|---|
| WP1 take stock | **done** — 22 modules, 15 docks, 29 tabs, 9 further windows, inventoried by opening each file |
| WP2 typed subjects | **done** — `parascape/ports.py`, 13 subjects, 50 tests, no Qt import |
| WP3 one catalogue | **done** — `parascape/modules.py` extended; the tool specs and shell routes computed from it; 25 tests |
| WP5 symbols, no cross | **done** — `parascape/iconset.py`, `parascape/ai_button.py`, `parascape/no_cross.py`; 16 tests |
| WP6 recipes | **done** — `parascape/recipes.py` (named `module_recipes.py` until 2026-11-12), `parascape/recipes_panel.py`; 38 tests |
| WP8 the gene card | **done 2026-11-12** — see below |
| WP4, WP7, WP9, WP10, WP11 | not started; they are the ones Einar will actually see |

**Two tables that can now disagree, which is the exact failure spaCR documents as bug #16i.** Both must be
closed before WP7, because the tree is computed from them:

1. `iconset.MODULE_ICON` maps a module key to a symbol, while the catalogue row has its own icon field. One
   of the two has to feed the other.
2. `modules.DOWNSTREAM` is still hand-written and sits beside a computed `ports.next_modules()`. The
   hand-written one is what the navigator offers today. WP7 should delete it and read the computed answer,
   which is the whole point of WP2.

**What WP6 did not wire, and the wording risk (closed 2026-11-12, see "The word recipe" below).** `parascape/recipes.py` meant something else -- the
Infer biological-question recipes -- so the settings recipes lived in `module_recipes.py`. Map and Network have
a Recipes group; Gene, Infer and the dock modules do not, because a gene or a question is a subject rather
than a setting. **If both kinds ever appear in the Infer module the word "recipe" will mean two things on one
screen**, which is a naming decision for Einar rather than a thing to guess at.

**Three agents shared one working tree and that was a mistake of mine.** One agent's uncommitted edits made
another's test run fail (`test_the_dock_watcher_ignores_the_window_appearing`), which was then reported as a
possible defect. It passes on committed code. Later work packages should use separate worktrees.

## Einar's decisions, 2026-11-12

**1. The settings recipes keep the word "recipe"**, matching spaCR, where a recipe is always a saved settings
bundle. The Infer biological-question recipes are renamed -- "question templates" or a better word -- across
`parascape/recipes.py` and the five modules that import it, and `module_recipes.py` then takes the plain name.
The point of the restructure is that the word means one thing in both applications; two meanings on one screen
would undo it.

**2. The gene card path is built first** (WP8, then WP9), ahead of the tree and the search. It is the exact
journey that dropped Einar into the old window: type a gene, click the module offered, land in
`workflows.py`. It is therefore the thing he will notice fixed.

**3. The tree keeps ParaScape's own headings** -- Start / Look up / Explore / Ask / Bring your own / Data --
rather than mirroring spaCR's sections. They describe what the user is trying to do; spaCR's describe an
imaging pipeline that ParaScape does not have. There are no segmentation models or plates here, and spaCR's
own rule is that no declared section may be empty.

### What WP8 must close first

`modules.DOWNSTREAM` is hand-written and is what the navigator offers today. `ports.next_modules()` computes
the same answer from declared subjects. WP8 deletes the hand-written one and reads the computed one -- that is
what makes clicking the gene card reveal the right modules rather than a list someone maintained by hand. The
same applies to `iconset.MODULE_ICON` against the icon field of a catalogue row: one must feed the other.

## The word "recipe", settled 2026-11-12

**A recipe is a named bundle of the settings of one module, and nothing else**, here and in spaCR. It lives in
`parascape/recipes.py` (`Recipe`, `save_recipe`, `apply_recipe`), is shown by `parascape/recipes_panel.py`, and is
kept as one JSON file under `~/.parascape/recipes/<module>/`. The module was `module_recipes.py` until this day.

**The way of asking one biological question in the Infer module is a question template**: one question, the
measurements allowed to answer it, the label to hold out and the independent control. It lives in
`parascape/question_templates.py` (`QuestionTemplate`, `TemplateResult`, `TemplateStore`; `questions.as_template`
turns a catalogue entry into one; `report.question_pdf` writes the report of a run). It was chosen over "plan",
"test" and "design" because the sentences already written around it fit: a template is reusable, is stored
without data, and is applied to whichever table is loaded, and nothing else in the Infer code uses the word.
The screen still says "question" wherever it can, since that is what a person picks.

**What was kept.** Files written by the old `RecipeStore` have the same fields and still load, and the format of
the settings recipes (`"parascape_recipe": 1`) did not change. A name that used to be importable from
`parascape.recipes` (`run`, `close`, `RecipeStore` and the rest) now raises an `ImportError` that names the new
place. No `RecipeStore` folder is created by the application itself, so there is no stored question template
on any machine except the ones a person made by hand.

**Not renamed, on purpose.** An embedding keeps a "full recipe" (`EmbeddingSpec`, stored in a run or an embedding)
because that is a bundle of settings and the meaning is the settled one. Two other uses are a third meaning and
are left for a later decision: `claims.Recipe` and the `claim_recipes` table (which generator and verifiers a
label uses), and `baselines.BaselineRecipe`. The text of the catalogue `questions.json` still says "recipe" in
the reasons of 13 questions; it is hashed into the data release, so it was not edited, and none of those
reasons is shown on screen.

## WP8 as built, 2026-11-12

Einar's words, which were the acceptance criteria: *"if i type a gene into the navigator i get a few moduals
that pop up, one should be the gene card (overhauled to mach the spacr module structure but with the same
information as before including the parasite schematic), when i click the gene card once i should see the
moduals downstream of one gene (like the network in spacr). so all the modules that take 1 gene and do
something should come downstream of the gene card module."* And: *"now when i put in a gene and pressed the
modual that came up for the gene, i was brought to the old version of parascape. this also needs to be
overhauled."*

**Both tables that could disagree are closed.**

* `modules.DOWNSTREAM` is **deleted**. `modules.downstream(key)` is computed from the declared subjects: for
  each subject a module hands on, in the order it declares them, the modules that read that subject, ready
  ones first and blocked ones last -- which is the order spaCR sorts its own next steps in
  (`chaining.next_steps`: `(not ok, module)`). Optional readers count as well as required ones, which is the
  difference from `ports.next_modules`: Network reads a gene it can do without, and leaving Network out of
  what follows a gene would have been wrong in the one place Einar was looking.
  `modules.downstream_of_subject(kind)` answers the narrower question, and
  `downstream_of_subject("gene", without="gene")` is the six modules Einar asked for.
* `iconset.MODULE_ICON` is **built from the catalogue** (`iconset._load_module_icons`). The icon field of a
  row now holds the NAME OF THE spaCR ARTWORK (`profiler`, `umap`) rather than a second word for the module,
  which is the fact it was pretending to hold while the real one lived in the icon set.

**No port declaration needed correcting.** `consumers_of("gene")` returns cell, expression, gene, information,
network, partners and structure, and opening each of the seven confirmed it: every one reads a gene, and no
module that reads a gene is missing from it. Map was checked and deliberately left out -- it declares a gene
LIST and a trait, and the single gene the Gene card sends it is sent as a selection, not as a port.

**The gene card** (`parascape/gene_module.py`, masthead "Gene card"). Nothing was dropped; three things were
added and the order was put right:

* Left column, in the order the mould asks for: **Your gene** (field, Choose a list, Use the gene list,
  **Clear** -- new, the mould asks every module for a way back to nothing), **Where it is, read by** (the lens
  -- new), **What to read** (the choice that changes the meaning), **Advanced** folded, **Save** folded (the
  two bundle buttons moved here out of the action row).
* Right column: **Where it is in the parasite** first, large, with focus -- **the parasite schematic, which
  survived**; then What is known, What is unusual, The evidence, Links to other organisms (folded), What is
  not known and what can be inferred, What each method said held out (folded), **What you can do with this
  gene next** (new), Actions, Console and AI chat (folded).
* The schematic is a second `celldiagram.CellDiagram`, not the one from the `find` dock moved: the rule of the
  mould, "move, do not copy", is about CONTROLS whose methods must keep reading them, and the `find` dock is
  lent whole to the Cell module, which could not then borrow a widget the gene card had taken. It is fed from
  `window.color_of` -- the one palette the window itself draws with, never a second one -- and clicking an organelle ticks that
  compartment in the find panel, so it works in both directions as it does beside the map.
* Where the gene sits is a pure reading (`gene_reading.placement_columns`, `placed_at`, `placement`,
  `placement_sentence`, `placement_trust`, `preferred_placement_column`), tested on a made-up table. A column
  counts as a placement when the drawing can colour MORE THAN HALF of what it holds: "does anything in it
  match" admitted `compartment_best`, four fifths of which are coarser groupings the artwork has no shape for.
* Trust against chance, with n, under the drawing: "Where a gene is, is measured for 3,827 of the 8,140 genes
  in this table (47 %). The commonest answer is nucleus - chromatin, which covers 20 % of the measured genes
  -- so a guess with no evidence at all would be right about 20 % of the time." A gene with no measured
  compartment says **unknown**, in words, and the parasite is drawn grey.
* One palette: when the window is coloured by something the drawing has no shape for, the card moves the
  window onto the measurement it is reading and says so, rather than mixing in a colour of its own.

**What clicking it reveals** (`parascape/next_steps.py`, ported in shape from `spacr/qt/chaining.py`): the six
modules that read one gene -- Partners, Expression, Network, Cell, Information, and Structure greyed because
no structure view is built -- each pre-filled with the gene, each checked with `ports.check_ready`, and the
blocked one kept in the row with its reason and its fix in the tooltip rather than dropped. The same offers
appear twice, on purpose and from one function: as the row at the foot of the card, and as the first column of
the navigator module network when the gene card is clicked once.

**The route is fixed.** `landing_install.open_item(window, "gene", ...)` opened the gene as a NODE, which took
the landing away and raised the evidence dock -- the pre-navigator window Einar landed in. It now opens the
gene card with the gene in it. Typing a gene in the navigator offers the gene card first and the other
single-gene modules after it, each carrying the gene; `open_module` gained `subject` and `value`, and
`carry_subject` selects the gene in the window so every panel behind the card is reading the same gene.
`workflows.py` is untouched and still reachable from its menu -- retiring it is WP10 -- but nothing on the
gene route reaches it, and `test_module_shell.py` holds that.

**Tests**: `tests/test_next_steps.py` (9, new), the placement readings in `tests/test_gene_reading.py` (7,
new), the gene card in `tests/test_module_shell.py` (10, new), the reveal in `tests/test_navigator.py` (8,
new), and the two closed tables in `tests/test_module_catalogue.py` (5, new). Verified by rendering the real
`app.Window` offscreen at 1920x1080 with `XDG_CONFIG_HOME` in a scratch directory and looking at it.
