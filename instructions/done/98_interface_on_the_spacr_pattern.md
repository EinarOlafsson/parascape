# 98 · Restructure the ParaScape interface on the spaCR pattern

**Status: built 2026-11-12, waiting on Einar.** All eleven work packages are built: typed subjects
(`ports.py`), one catalogue (`modules.py`), the navigator search, the spaCR symbols with no small cross
anywhere, recipes, the tree down the left, the gene card and the gene list card with the modules downstream
of each, the pre-navigator window retired with every panel behind it named and routed, and the whole of that
held as a contract over every module (`tests/test_module_contract.py`), so a module added without a part of
it fails the suite by name rather than being noticed later by Einar.

**Two things are left and neither is work this instruction can finish.** (1) The star map: delete it once the
Network module draws the inferred and your-run links, or keep it -- a decision for Einar, with the code
untouched either way. (2) `docs/deck/` slide 4 still shows the retired workflows window; the slide images are
rendered from a running application, so bringing them up to date is a pass of its own.

*Requested by Einar 2026-11-12, after installing 0.0.11.0.*

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
| WP9 the gene list card | **done 2026-11-12** — see "WP9 as built" |
| WP10 retire the old windows, route every orphan | **done 2026-11-12** — see "WP10 as built" |
| WP7 the tree | **done 2026-11-12** — see "WP7 and WP4 as built" |
| WP4 navigator search | **done 2026-11-12** — see "WP7 and WP4 as built" |
| WP11 consistency as tests | **done 2026-11-12** — see "WP11 as built" |

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

## WP9 as built, 2026-11-12

Einar's words, which were the acceptance criteria: *"if i type in a list of genes, the same logic but this
time not gene card instead a gene list module and after clicking it all the moduals that come downstream of
gene lists."*

**The one sentence the mould asks for.** This module lets you take a list of genes and see what the genes of
it share that other genes do not -- and how much of the list is measured at all -- so that you can tell what
the list points to, and which of its gaps are worth inferring.

**The card is the `list` tile, built rather than greyed.** WP1 recorded List as one of the four modules
honestly greyed out, "needs a set-level enrichment view". That view is what WP9 is, so the row lost its
`needs`, gained `builder="shell:list_module"` and is named **Gene list** rather than List, to pair with the
gene card.

**What the card shows, and why.** A list card that is only a table of identifiers fails the first goal of the
mould: it makes the information no more accessible than the file the list arrived in, and synthesises nothing.
So the identifiers are the FOLDED section at the bottom, and what opens first is the only thing about a list
that cannot be read off the genes one at a time:

1. **What the genes of the list share** (the thing itself, large, with the focus). Every label any gene can
   carry, counted inside the list and outside it, kept where the list holds it at least 1.5 times as often as
   the rest AND rarely enough that a random list of the same size would do it less than once in 20 tries. Both
   bars are said out loud in the line above the table, because a bar that hides things without saying so is a
   bar nobody can argue with. Each row: the label, the measurement, how many of the list carry it out of how
   many of the list are measured for it, the same share among the other genes, and how surprising it is in
   plain words ("once in 2,122 random lists"). A share with nothing to compare it against cannot be read --
   half a list in the nucleus means nothing until you know what half the genome does.
2. **Where the genes of the list are in the parasite** -- the schematic Einar named, on the list as on one
   gene: a second `celldiagram.CellDiagram`, fed from `window.color_of` so the card and the map cannot
   disagree about what a colour means, coloured for the compartment most of the list is in, with how much of
   the list is in it against how much of the rest of the genome is. Clicking another organelle reads that one
   instead and ticks it in the find panel. A list with nothing placed says **unknown**, in words, and the
   parasite is drawn grey.
3. **How much of the list is measured**, per kind of knowledge, least measured first, each against the share
   of the genes OUTSIDE the list -- which is what tells a badly measured list apart from a kind of knowledge
   nobody has measured for anything. They are different problems with different answers.
4. **What is unknown for most of the list**, which is where inference is the only way forward and is therefore
   what the card hands to Infer -- the second goal of the mould.

Every sentence is computed by **`parascape/list_reading.py`**, which has no Qt in it: `split_identifiers`,
`resolve` (an identifier that names no gene is reported, never dropped), `shared_labels` with `merge_by_value`,
`coverage`, `unknown_for_most`, and `spread` for the drawing. 27 tests on a made-up table.

**What clicking it reveals**: the eight modules that read a list -- Gene, Partners, Expression, Map, Network,
Strategies, Export ready, and **Infer greyed with its reason** ("a question is needed and none is in hand"),
kept in the row rather than dropped. The same offers appear twice from one function (`next_steps.offers`, the
same one the gene card uses): as the row at the foot of the card, and as the first column of the navigator
module network when the gene list card is clicked once.

**The port declaration that was wrong, and the correction.** `consumers_of("gene_list")` returns nine modules
and all nine were right -- each was opened and each does read several genes. What was wrong was what
`check_ready` said about three of them. **Gene, Partners and Expression each declared one gene as REQUIRED and
a list as optional**, which reads as "it needs a gene and may also be given a list". What they actually do is
read either: the gene card steps THROUGH a list (its own port says so, "a list to step through, one gene at a
time") and has a button for taking the gene list whole; the Partners port says "several genes, whose partners
are shown together"; and the field of the expression panel says "Enter one or more exact reference gene IDs".
So somebody holding a list was told all three could not work on it, in the one place a list is the thing in
hand.

`ports.Port` gained **`one_of`**: the name of a set of ports of which ANY ONE will do. Only `check_ready`
reads it, `required_kinds` is unchanged -- what the module NEEDS has not changed; what changed is the
admission that two different things can meet the need. With nothing at all in hand nothing in the set is
filled, so those modules are still reported as not ready, now with a reason that names both ways of fixing it
("one gene or a list of genes is needed and none is in hand").

Two further corrections, both of wiring rather than of a declaration:

* **The `list` module produces the list itself**, declared FIRST, before the trait it was already declaring.
  `modules.downstream` walks the produced ports in declared order, so with the trait first the first column of
  the module network under the gene list card was the readers of a trait. It also produces one gene, because
  a gene of the list opens on its own from the folded table.
* **Export handed a list now writes the list.** Its catalogue row opens `export_visible`, which writes every
  gene passing the class filter -- so the offer "Export, with your 28 genes already in it" would have written
  eight thousand of them. `carry_subject` gates a list carried into the window (the gate is what the
  application already means by "this set of genes is in hand": the map marks it, the strategy panel fills a
  gene box from it, `describe_gate` reports it), and `open_module` routes Export to `export_gated` when a list
  is what is in hand.

**The navigator routes a typed list.** `genes_in_text` recognises several identifiers that name several genes
-- commas, semicolons, spaces, tabs or new lines, in any mixture, matched ignoring case and spelled as the
table spells them. One identifier is still a search, and three words that name no gene are still a search.
The resolution is one vectorised match through a `resolve_list` callable handed to the Navigator by
`landing_install`: the gene lookup beside it is a SEARCH taking the whole query, which asked for a pasted
column of 300 accessions looks for one gene whose product contains that whole string. The genes themselves
are kept as ONE group card -- one card each would be three hundred cards -- which opens to its members.
`open_item(window, "gene_list", ...)` opens the card, as `open_item(window, "gene", ...)` opens the gene card.

**Three things found by looking at the render that no test caught**, which is why the mould asks for it:

1. The card opened with FIVE rows all saying "dense granules" -- one finding said once per measurement that
   agrees with it, filling the top of the card and pushing everything else off the screen. `merge_by_value`
   keeps the strongest and names the rest ("lopit map (and 2 more that agree)"), so the agreement is reported
   as what it is -- corroboration -- rather than as five findings.
2. "How much of the list is measured" was fifteen paragraphs in a pane that showed two of them. As a sorted
   table with one summary line above it, the same height shows six and they can be read down the column.
3. The column of EIGHT next steps shrinks each card to an 84 px square, and "Expression" was drawn as
   "xpression" -- clipped at both ends, so neither end said which module it was. `SuggestionCard.fit` now
   elides. The details panel under it printed "Opens with TGME49_248510,TGME49_304720,…"; it says "Opens with
   28 genes" now (`navigator.subject_words`, and `next_steps.Offer.shown` for the row on the card).

`FLOW_FIRST` is now **None -- all of them**, because Einar asked for all of them and a number here was wrong
twice in one day: six was right for a gene until a list wanted eight, and eight was right for a list until
WP10 landed Predict and it wanted nine. The column shrinks its cards to fit and elides their names rather
than dropping any, so there is nothing for a limit to protect. The first column after a clicked card is also
computed from the SUBJECT IN HAND (`NavigatorPage.following`) rather than from the first N of everything the
module hands on -- which only happened to be right for a gene because there were exactly six of those and the
column showed six.

**Tests**: `tests/test_list_reading.py` (27, new), `tests/test_list_card_in_window.py` (15, new, on the real
`app.Window`), the typed list and the reveal in `tests/test_navigator.py` (10, new), the offers in
`tests/test_next_steps.py` (5, new), and `Port.one_of` in `tests/test_ports.py` (4, new). Verified by
rendering the real `app.Window` offscreen at 1920x1080 with `XDG_CONFIG_HOME` in a scratch directory, typing
a list of 28 dense-granule genes, opening the card, clicking it, and looking at every image.

## WP10 as built, 2026-11-12

The words of Einar, which are the acceptance criteria: *"now when i put in a gene and pressed the modual that came
up for the gene, i was brought to the old version of parascape. this also needs to be overhauled."* And:
*"the new interface is really close to what i wanted in terms of the navigator screen. but everything behind
it seems to be unchanged."*

WP8 stopped the GENE route reaching the old window. WP10 takes the window away, and gives every view behind
the navigator a name and a route.

### The mechanism: a view is declared on the module it belongs to

Twenty-nine more tiles on Home would have answered "everything behind the navigator is unchanged" with a
wall. Deleting twenty-nine panels would have dropped measurements. The middle answer is `modules.View`: a
name, one plain line saying what data is on screen, and one destination (`dock`, `dock:tab` or
`action:method`), declared on the catalogue row of the module whose panel already holds it. There is no
second table -- the search index (`help_index.panel_entries`) is computed from these rows, and
`landing_install.open_view` is the one route to them. A view of a dock the module borrows opens **inside the
module shell** with its tab brought forward; a view somewhere else leaves the landing and opens there,
honestly, as a panel. `tests/test_views.py` holds both directions: every declared view exists in a real
window, and every tab in a real window is a view, a module, or a signpost left where a module took the
controls.

### 1. `workflows.py` -- RETIRED, after moving what only it could do

The three tabs were checked by opening them, not by their names.

| tab | what it alone could do | decision |
|---|---|---|
| Explore a gene | gene search, then the evidence rows of that gene with sources | **(c) deleted** — the Gene card (WP8) does all of it and more, and this was the screen Einar landed on |
| Predict a trait | `prediction.run`: eight methods, holdout by family, folds, abstention, the full run export. **Nothing else in ParaScape calls it** | **(b) its own catalogue row** — the new `predict` module, controls MOVED into `parascape/predict_module.py` |
| Compare a screen | `prioritization.compare_screen`: your CSV joined to the gene table with the source hash kept. **Nothing else in ParaScape calls it** | **(a) routed to the module that claims to do it** — `screen`, rebuilt as `parascape/screen_module.py`, controls MOVED |

`parascape/workflows.py` and `tests/test_workflows.py` are deleted, with the three Tools menu items,
`Window._open_workflows`, the `workflow` destination kind of `nodes.py`, the branch in `home_panel`, the
branch in `help_search.open_panel` and the three rows of `help_index`. `tests/test_module_shell.py` holds
that the file is gone from the package and that both replacements were built.

**The Screen tile was wrong before this, which is the complaint in one click.** It said "join your screen to
everything already recorded" and opened the `compare` dock -- somebody else published RNA-seq. It reads a
screen now, and its declared ports (`SCREEN` in, `GENE_LIST` and `TRAIT` out) finally describe what it does.

### 2. `star map` -- KEPT, and the deletion is Einar to make

WP1 recorded it as dead weight the Network module replaced. **Opening both says otherwise.** The Network
module reads `star_edges.measured_edges` alone. The star map also draws `shipped_edges` -- the inferred
links of the 39 strategies -- and `UserEdgeStore`, the links from runs of your own made in the Strategies
tab, kept on disk between sessions, with saved link definitions (`DefinitionStore`). Deleting it would lose
the only view of every inferred link in ParaScape, which is the one thing WP10 is not allowed to do.

So: **(c) candidate for deletion, needs Einar**, and the code stays. It is now a declared view of Network --
"Links with their provenance" -- so it has a name and a route instead of being a dock nobody could ask for.
The two honest options are for Einar to choose between: either the Network module absorbs the inferred and
your-run layers (a work package of its own, and the right end state), or the star map stays as the one place
links are drawn with their provenance.

### 3. Every orphan tab and window, with its decision

Decisions: **(a)** a view of a module that already exists, routed there; **(b)** a catalogue row of its own;
**(c)** dead. Each was opened and read; none was judged from its name.

#### The `analysis` dock, 8 tabs

| tab | what it shows | decision | reason |
|---|---|---|---|
| 1 · Data | the 96 measurement blocks as a tree in three hierarchies, each with its coverage, plus the scaling and missing-value rules | **(a)** view of **Map**, "Choose individual measurements" | the Map module "What makes genes similar" is the coarse version of this, and its Advanced already had a button here; the fine choice had no name |
| 2 · Map | *(moved into the Map module in instruction 86; the tab is a signpost)* | — | already done |
| 3 · Clusters | *(moved into the Map module; the tab is a signpost)* | — | already done |
| 4 · Inference | the held-out battery: one row per measurement the map never saw, scored against the groups it found, with the best single group per value | **(a)** view of **Infer**, "What the groups mean, held out" | the only place the battery runs, and the Map module already sent people here by a button with no name |
| 5 · Search | walks combinations of data and settings, scoring how well a held-out label comes back | **(a)** view of **Map**, "Find the best map for a label" | the Map module has a folded section fed by this walk; the target choice and the per-run table live only here |
| 6 · Validation | hides a share of the genes that carry a label, picks the group holding most of the rest, scores against the hidden ones, per category | **(a)** view of **Scorecard**, "Hide labels and see if they come back" | it measures how often a call of that kind is right, which is what Scorecard is for |
| 7 · Discover | hill-climbs maps and groupings for structure that PREDICTS, read back as ranked claims with their caveats | **(a)** view of **Infer**, "Hill-climb for claims" | the only place the discovery climb runs, and its output is claims, which is what Infer produces |
| 8 · Questions | *(moved into the Infer module; the tab is a signpost)* | — | already done |

#### The `compare` dock, 11 tabs

The Expression module borrows this dock, so its ten other tabs are views of Expression. Ownership here is
not a judgement about subject matter: it is the module whose main section literally contains them, which is
the one fact a test can check.

| tab | what it shows | decision |
|---|---|---|
| Gene expression | where and when the RNA and protein of one gene are measured | — it IS the Expression module (`opens`) |
| Studies & samples | each quantified study, its samples, and how they were grouped | **(a)** view of Expression |
| Coexpression | the genes whose expression rises and falls with the one chosen, per source | **(a)** view of Expression |
| Published expression | the values the original papers published, as published | **(a)** view of Expression |
| Deposited matrices | whole deposited measurement matrices, read from the source file | **(a)** view of Expression |
| Other measurements | deposited values that are not expression, with source and units | **(a)** view of Expression |
| Sparse measurements | measurements recorded for a few genes, kept sparse rather than filled with zeros | **(a)** view of Expression |
| Original workbooks | the supplementary workbooks as the authors shipped them, sheet by sheet | **(a)** view of Expression |
| Reference identifiers | which identifier in a source maps to which gene here, and what did not map | **(a)** view of Expression |
| Dataset coverage | how many genes each dataset actually covers | **(a)** view of Expression |
| Condition comparison | two groups of samples from one study: which genes differ, by how much, how sure | **(a)** view of Expression |

Nothing here was deleted. Every one of the ten reads a different file and shows numbers that appear nowhere
else; "dead" was not true of any of them.

#### The `discoveries` dock, 5 tabs — all views of **Scorecard**, which borrows it

| tab | what it shows | decision |
|---|---|---|
| Labels and functions | every label, how many genes carry it, where that membership was recorded | **(a)** |
| Inferred claims | the claims already computed, each with its certainty and whether anything independent agreed | **(a)** |
| Functional tests | the frozen functional benchmarks and what each method scored on them | **(a)** |
| Strategy coverage | which strategies have been run on which labels, and which pairs have not | **(a)** |
| Host source status | each host-facing source, what it covers, whether it is loaded | **(a)** |

#### The `information map` dock, 2 tabs — views of **Information**, which borrows it

| tab | what it shows | decision |
|---|---|---|
| One gene | one gene trait by trait, each cell marked measured, a measured no, or not measured at all | **(a)** |
| One trait | one trait across every gene: coverage by source tier and the values it takes | **(a)** |

#### The `strategies` dock, 3 tabs — views of **Strategies**, which borrows it

| tab | what it shows | decision |
|---|---|---|
| Guide | what the chosen strategy does and what it may use | **(a)**; `open_view` calls `show_details()` first, because these sit behind the card |
| Settings | every setting of the chosen strategy, with the range each is allowed | **(a)** |
| Results | what it returned, its verdict and its tables | **(a)** |

#### The two docks that had no module at all

| dock | what it shows | decision | reason |
|---|---|---|---|
| `maps` | the maps that ship with ParaScape, grouped by what each was built from, with how well every label maps onto every one of them | **(a)** view of **Gallery**, "Maps shipped ready-made" | Gallery is "every map built so far, to compare rather than rebuild", which is this; the label-by-map score table exists nowhere else |
| `start here` | one question at a time, ending at the strategies worth running | **(a)** view of **Strategies**, "Start here" | it is a way INTO the strategies and ends on their cards; a tile of its own would be a second front door beside the navigator, which is the thing being removed |

#### The further windows

| window | what it shows | decision | reason |
|---|---|---|---|
| Dataset browser | all sources, their coverage and provenance | — it IS the Datasets module | already reachable |
| Slot tree (`slot_tree.py`) | every information category in three hierarchies, what fills it, the empty ones coloured | **(a)** view of **Datasets**, "Slot tree" | it is the catalogue of what is measured, which is what Datasets is about; the Tools toggle stays |
| Gene evidence browser (`gene_evidence_browser.py`) | one gene, every record behind it, missing values included | **(a)** view of **Gene**, "Every record behind one gene" | it was reachable only by clicking a link in the evidence dock. `Window.open_gene_evidence` now says in the status bar when no gene is chosen, instead of returning None silently |
| Label space browser (`label_space_browser.py`) | label membership, class balance, source records, recorded scorecards | **(a)** view of **Scorecard**, "Labels and classes" | it is memberships and the held-out scorecards against them |
| Question scoreboard (`question_scoreboard.py`) | every registered question scored held out, and whether it may fill the map | **(a)** view of **Questions**, "Question scoreboard" | Questions is "the catalogue of what can be asked, and of what cannot"; this is how each one has done |
| Workflows (`workflows.py`) | three tabs | **(c) DELETED**, after moving two of them into modules — see above |
| Preferences | Appearance, Display, Window | — left in the menus | settings are not a module, and the search field still finds each page |
| Download dialogs | what can be fetched and what is already here | — left in the menus | the same |

### What else changed, and why

* **`opens` is empty for Map, Network, Predict and Screen.** It is where a module GOES, and these four go to
  a shell of their own. It used to name a dock each was loosely associated with, which greyed Network out
  whenever the star map dock failed to build -- a dock Network does not use. A module with a `shell:` builder
  is no longer greyed for want of a view (`modules.unavailable_reason`), which is also what lets Predict and
  Screen be available at all: they come from a window that no longer exists.
* **Two new rows in the catalogue**, so 23 modules rather than 22, and 32 declared views.
* **`ports.py` gained `predict`**: it takes a trait and optionally your genes, and hands on a claim, a
  scorecard, a saved run and a job. It is the free-hand counterpart of Infer, where the question arrives with
  its holdout already chosen and checked; the module says so on screen rather than leaving the two looking
  interchangeable.

### Verified by looking

The real `app.Window` was rendered offscreen at 1920x1080 with `XDG_CONFIG_HOME` in a scratch directory, and
walked: the navigator, both new modules, and every declared view. **No view reported "not available", and no
menu item pointing at the retired window was left anywhere in the menu bar.** The Coexpression view opens
inside the Expression shell with its tab forward; the star map view leaves the landing and raises its dock.

### Tests

`tests/test_views.py` (44, new). Changed: `test_module_catalogue.py` (the six shell builders),
`test_ports.py` (Predict produces a job), `test_module_shell.py` (the retired window is gone from the
package), `test_home_panel.py` and `test_glass.py` (their `workflow` branches).

## WP7 and WP4 as built, 2026-11-12

The words of Einar, which are the acceptance criteria: *"there should be a tree of things that parascape can
do, and this progression should guide the user to the modual they want. just like in spacr. the one in spacr
works really well, so study it and implement conceptually the same thing in parascape."* And: *"the navigator
should have a search bar like spacr."*

### The tree (WP7)

**Where it is.** The rail down the left, on the navigator and beside every module shell alike. It was a flat
list -- six headings, and under each the modules of that band -- which answered "what is there" and not "what
do I do next", and could not name the 29 tabs and four windows behind the navigator at all.

**How it is computed.** `parascape/tree.py`, no Qt, from the one catalogue and the declared subjects and
nothing else. `parascape/module_tree.py` turns its rows into widgets and decides nothing. Four levels:

```
Start                        the progression. Open when the rail appears; the other five shut
  Home                       back to the Navigator and the module tiles
  One gene  ▾                a subject you can bring -- from ports.SUBJECT_SOURCES
    → Gene  → Partners  → Expression  → Network  → Cell  → Information  → Structure (greyed)
  A list of genes · A trait · A screen of your own · A question · A job
Look up  ▾                   the five bands of the catalogue, in catalogue order
  Gene ▸ · Structure · Partners ▸ · Expression ▾
    Studies & samples … Condition comparison      the 10 views behind Expression
    → Map  → Cell  → Infer  → Predict  → Strategies  → Scorecard  → Export  → Information
Explore · Ask · Bring your own · Data
```

* A **heading** is one of the six. Click to open or shut; shutting one forgets which module under it was
  unfolded, because reopening it to a module unfolded minutes ago is not what "show me Look up" means.
* A **subject** is the progression, and it is the one real difference from spaCR. spaCR progresses because its
  sections are stages of one imaging run; ParaScape has no single run. What it has is the subject in hand, so
  the tree answers "what becomes available once something is in hand" with
  `modules.downstream_of_subject(kind)` -- the same function the gene card and the gene list card offer, ready
  ones first and a blocked one kept with its reason.
* A **module** is filed under the heading its catalogue row declares, with what it TAKES and what it HANDS ON
  in its tooltip, from `modules.takes_words` and `modules.hands_on_words` (which read the ports). Clicking it
  opens it **and** unfolds it; only one module per heading stays unfolded, because nine modules each unfolding
  ten views is a rail past the length of any screen.
* Under a module: **its declared views** (all 32 are in the tree, which is the point -- a view not in it is a
  panel behind the navigator that nothing can name), then **what follows it**, `modules.downstream(key)`.

**Does any module need two headings? No, and the tree is built so that the question has one answer.** Each of
the 23 is filed once, under its catalogue band. The rows under *Start* are REFERENCES
(`tree.Node.reference`), marked as such, and `tree.find` and `tree.section_of` never return one, so "where
does this module live" has one answer. `test_module_catalogue.py::test_no_module_is_in_two_sections` is
**unchanged and still passes**, because it is about the catalogue and the catalogue did not change; the new
`tests/test_tree.py` holds the tree side of the same rule. Offering a module under a subject is not filing it
twice -- it is naming it where it is useful, which is the whole of the progression.

**Start is a heading and not a band**, and deliberately. It holds no module: Home is not one and a subject is
not one, and `modules.BANDS` is the declaration that a section HAS modules (`test_no_section_is_empty` says so
by name). spaCR keeps `SECTION_HELP` out of its own `SECTION_ORDER` for exactly this reason and records that
putting it in was tried first and reverted.

### The search (WP4)

**What it ranks on**, in order, and the order IS the ranking. The bands are spaCR's
(`spacr.qt.help_index._term_score`) with one added, declared in `modules.py` beside the reason:

| band | | matched against |
|---|---|---|
| 100 | typed in full | the name of the module or of the view |
| 80 | a whole word of the name | |
| 60 | the name begins with it | plus a little for how much of the name was typed |
| 50 | a word of the name begins with it | |
| 40 | somewhere inside the name | |
| **30** | **one of the OTHER WORDS declared for it** | `Module.words` -- "interactome", "LOPIT", "CRISPR" |
| 20 | the letters in order and close together | the fuzzy spelling |
| 10 | what it DOES | the note, what to do first, the section, **every subject it reads and hands on**, and what a view says is on screen |

The words of a query are ANDed, so a second word narrows. Modules and views are ranked **together**
(`modules.search_all`): the best answer to "coexpression" is a tab of the Expression panel and the best answer
to "expression" is the module, and two lists stapled together could not say which. A tie puts the module
first.

**What was added and why each.** Module names and what a module does were already matched, against the name
and the one-line note only. New: the **32 view names**, so WP10's naming of everything behind the navigator is
reachable by typing it (before this, a view could be found only through the help index of a LIVE window, and
then only as a row called "panel"); the **subjects from the ports**, so "trait" finds every module that reads
one without the word being written on any row; and the **words a parasitologist uses** -- "interactome",
"pulldown", "LOPIT", "rhoptry", "CRISPR", "UMAP", "AlphaFold", "RNA-seq", "gene set", "enrichment" -- none of
which is the name of anything in ParaScape, so each of them found nothing at all before. They are declared on
the catalogue row rather than in a synonym table, because a second table about the same 23 modules is the
failure this instruction exists to remove.

**Not regressed, and tested:** a typed gene still offers the gene card first with the gene in it, and a typed
list still offers the gene list card. **The AI symbol is unchanged**: a search WITH matches shows cards and
never asks the assistant; only a search with no match does. That differs from spaCR on purpose and
`test_a_search_with_matches_shows_cards_and_never_asks_the_assistant` now says so by name.

### Four things found by LOOKING that no test caught

1. **The rail sized itself on the text and ignored the indent.** A row three levels in starts 38 px from the
   edge, so "Published expression", "Other measurements" and "Sparse measurements" were elided in a rail wide
   enough for every one of them.
2. **The module a named view belongs to was dropped.** Typing "coexpression" offered the Coexpression tab and
   no way to the module that shows it, because only cards whose TITLE holds the typed word survive the filter
   that keeps "map" from returning two dozen cards -- and "Expression" does not hold "coexpression". The owner
   of a named view is now offered straight after it, and only the owner.
3. **A strategy is named by a whole sentence** and a 224 px card with a 96 px symbol in it has room for about
   three lines. "20 · Learn what makes your list special, from positives…" ran to four and the last was drawn
   over the edge of its own card. Names on cards are shortened to what the card can draw.
4. **One view name fits in the rail at no width the dock allows** ("Hide labels and see if they come back"), so
   every row tooltip opens with the name of the row.

### One race, turned into a reliable failure and then fixed

`test_the_dock_watcher_ignores_the_window_appearing` is the test this instruction records as "reported as a
possible defect; passes on committed code". It was a race: the platform delivers the dock visibility changes
of a window being shown asynchronously, and whether the landing survived depended on whether a 400 ms
`QTimer` had fired first. The tree made the window a tenth of a second slower to put on screen and the race
became a reliable failure. `landing_install._arm_later` now records a DEADLINE that every show and every
resize pushes forward, and `_watcher_awake` reads it when the signal arrives, so there is nothing left to
race.

### Tests

`tests/test_tree.py` (34, new: what the tree computes, and the rail built from it). Added to
`tests/test_navigator.py` (12: the bands, the parasitologist words, the views, the assistant rule). The three
Home tests left red by WP9 and WP10 were corrected on nightly while this was being built (04585995) and are
not touched here.

### Verified by looking

The real `app.Window` rendered offscreen at 1920x1080 with `XDG_CONFIG_HOME` in a scratch directory and the
animation at its shipped default: the rail with Start open, Look up open, Expression unfolded with its ten
views and eight next steps, Scorecard unfolded with its seven, the progression under "One gene", and the
field answering "interactome", "coexpression", "map", "LOPIT", "what is known", a gene accession and a query
that matches nothing.

### Left for Einar

1. **The star map.** Delete it only once the Network module draws the inferred and your-run links, or keep
   it. The code is untouched either way.
2. **The deck** (`docs/deck/`) still says "Open these guided workflows from the Tools menu" on slide 4 and
   shows a picture of that window. The slide images are rendered from a running application, so bringing them
   up to date is a pass of its own rather than an edit. `docs/workflows.md`, `docs/index.md`, `docs/guide.md`
   and `docs/API.md` are correct as of this change.

## WP11 as built, 2026-11-12

The words of Einar, which are the acceptance criteria: *"every module should be thematically consistent. every
analasys should be thematically consistent and intuitive to use and get to."*

Ten work packages built that consistency. WP11 makes it enforceable, because a contract nobody checks lasts
until the next module is added.

### What was already held, and where

Most of the contract was already asserted somewhere, and **nothing of it is asserted a second time**. Two
tests of one rule are two things to keep in step; they drift, and the one that drifts is the one nobody reads.

| clause | already held by |
|---|---|
| declares its ports | `test_module_catalogue::test_every_module_in_the_catalogue_has_ports_declared`, `::test_the_two_tables_are_keyed_the_same` |
| has an icon that exists and renders | `test_module_catalogue::test_every_module_has_an_icon_name`, `::test_the_symbol_of_a_module_is_declared_once`, `test_symbols_and_no_cross::test_every_module_has_a_shipped_spacr_symbol`, `::test_every_symbol_renders_visible_pixels_in_the_ink_of_the_theme` |
| appears in exactly one section | `test_module_catalogue::test_no_section_is_empty`, `::test_no_module_is_in_two_sections`, `::test_the_grouping_covers_every_module_exactly_once` |
| is reachable in the interface | `test_module_catalogue::test_every_catalogue_module_is_reachable_in_the_interface` and its converse, `test_modules_in_window::test_every_destination_exists_in_the_window` |
| is reachable from the tree | `test_tree::test_every_module_is_filed_exactly_once`, `::test_the_rail_lists_home_and_every_module_once` |
| its views exist and open | `test_views` (four tests, both directions) and `test_tree::test_every_declared_view_is_in_the_tree` |
| carries a docstring | `test_docstrings` (three tests over every file of the package) |
| says what to do first | `test_module_catalogue::test_every_module_says_what_to_do_first`, `test_home::test_every_module_says_what_you_get_not_what_it_is` |
| a tooltip on every control | `test_app_menus::test_no_control_anywhere_is_left_without_an_explanation` — it walks the real window, and every module shell is built with the window, so it reaches them all |
| shows no "×" | `test_symbols_and_no_cross::test_the_real_window_has_no_small_cross_anywhere`, `::test_no_cross_survives_opening_every_module` |

**Six clauses had no test anywhere**, and each of the six is a failure that has actually happened here rather
than one imagined for a test:

1. **Reachable by typing its own name.** `test_navigator` checked one module. It is a property of every row.
2. **A recipe surface wherever there are settings.** Map and Network had one; the other six shells did not.
3. **A control only for a subject the module reads.** Every tool module was given a "Your gene" field.
4. **Somewhere to put what it requires**, so an offer carrying a gene does not open the module empty.
5. **The same words on the tile and on the masthead.** Five of the seven shells had already drifted.
6. **Greyed only for a reason its own row declares.** This is the Network-and-the-star-map failure.

### The one place the contract lives: `tests/test_module_contract.py`

`CONTRACT` is a table with one line per clause. Each line either carries a `check` that runs against every
catalogue row, or names the tests that already hold it. The table-driven test is parametrised by module AND by
clause, so a failure reads

```
test_the_contract[jobs-offers a control only for a subject it reads]
```

which names the module and the missing part without anybody reading the assertion. Three tests keep the table
itself honest: a clause with neither a check nor a holder fails; a clause with BOTH fails, because that is the
duplication this file exists to avoid; and every test a clause leans on is looked up by name in the file it
names, so a rename over there is a failure here rather than a clause that quietly stops being checked.

**The guard against a fourth #16i.** `test_nothing_keeps_a_second_table_about_the_modules` walks the package
for any module-level dict, tuple, list or set whose keys are four or more module keys, outside `modules.py`
and `ports.py`. `iconset.MODULE_ICON` and `modules.DOWNSTREAM` were each exactly that, and each drifted; both
are closed by name elsewhere. This finds the next one by SHAPE, because the next one will have a different
name. Four keys is the bar: two or three is a list of special cases with a comment beside it.

### What the contract found, and what was changed to keep it

* **A third pair of tables, found by WP11.** Each of the seven modules built on a shell of its own wrote its
  own name, note and first instruction beside the same three facts on its catalogue row, and five of the seven
  had drifted: the tile said "Gene" and the masthead "Gene card"; the tile said Network was "a gene and its
  neighbourhood, edges you define" and the masthead "every gene and the links between them, yours in focus".
  `ModuleShell.for_module(key)` now builds the masthead from the row, as `dock_module` has since WP3. Nothing
  was lost: where the masthead had the plainer words they were moved INTO the catalogue, so the tiles, the
  tree and the search read better too. **The `gene` row is now called "Gene card"** -- the phrase of Einar himself, and
  it pairs with "Gene list".
* **Six modules asked for a gene they do not read.** `dock_module` gave every tool module a "Your gene" field
  and a "Read it in Gene" action, so Jobs asked for "a gene id or a word of its product" above the list of
  what is running, and Console, Gallery, Drop, Strategies and Scorecard did the same. The field is built from
  `Module.consumes` now, so Expression, Cell and Information keep it and the other six do not.
* **Five modules had settings and no way to save them.** Gene card, Gene list, Infer, Predict and Screen each
  gained a Recipes group, the same one spaCR gives every screen. **The rule of Einar is encoded rather than
  guessed at**: `ModuleShell.add_setting(..., subject=...)` records which groups of the Settings column hold
  THE SUBJECT -- the gene field, the question chooser, the file opener -- declared where the group is built
  and by the only person who knows. A module whose Settings column holds nothing but subjects needs no recipe,
  which is why the nine tool modules have none and fail nothing.
* **The word "recipe" now appears in the Infer module**, once, meaning the settled thing. This instruction
  recorded the risk by name; the rename to "question template" is what makes it safe, and
  `test_the_infer_module_offers_question_templates_and_never_calls_one_a_recipe` was narrowed to say exactly
  that -- nothing about a QUESTION is called a recipe.
* **Infer offers what follows an answer.** Predict had the row and Infer did not, although the two are the two
  halves of one thing and both hand on a claim. Both now read `next_steps.offers`, so the two cannot come to
  differ about what follows a claim.

### Four more things found by LOOKING that no test caught

The real `app.Window` was rendered offscreen at 1920x1080 with `XDG_CONFIG_HOME` in a scratch directory, and
every module and every declared view walked and looked at, four times over as the fixes landed.

1. **An action row cut its buttons at both ends.** Eight actions in the Network module left "Open a saved
   bundle…" drawn as "pen a saved bundle" and "Save everything as a bundle…" as "everything as a bur" -- so
   neither end said which button it was and both read as a different, wrong action. A `QPushButton` clips;
   `shell.ActionButton` shortens with an ellipsis, which SAYS it has been shortened, and the whole label comes
   first in the tooltip. `text()` is unchanged, because that is what the module keyed the button by.
2. **An empty Settings column took a third of the screen.** Once the six tool modules stopped offering a gene
   field, their Settings column was an empty bordered box with the panel squeezed into what was left. A shell
   with no settings hides the pane and gives the room to what it shows -- which is how Jobs and Strategies came
   to be readable across the full width.
3. **A table cut the names of its own columns.** The Drop module said "NTIFIER IN THE F" for "identifier in
   the file", the Scorecard said "NFERRED CLAIMS", and the Information module cut "evidence (reference ·
   dataset)" in half, because the columns were sized only once something had been loaded.
   `parascape/panes/tables.show_headings` is the one fix for the one mistake.
4. **The Gallery pointed at a tab that no longer exists.** "Run *walk hyperparameters* on the Map tab" -- the
   Map tab became the Map module in instruction 86, and "walk hyperparameters" is the name of a button rather
   than of anything a person wants. It says "Open the **Map** module and use **Find the best map for a
   label**" now.

Each of the four is now a test: `test_a_squeezed_action_is_shortened_with_an_ellipsis_rather_than_cut`,
`test_a_module_with_nothing_to_set_gives_the_room_to_what_it_shows`,
`test_every_table_says_what_its_columns_are_before_anything_is_loaded`, and the Gallery wording is held by the
clause that every module says one thing in one place.

### What the contract is, and what it does not cover

**What it is.** Every row of `modules.MODULES` declares its ports, has an icon that exists and renders, sits in
exactly one section, is reachable from the dock, from the tree and by typing its own name, has views that exist
and open, has a recipe surface wherever it has settings of its own, carries a docstring, says what to do first,
has a tooltip on every control, shows no "×", offers a control only for a subject its ports declare, has
somewhere to put what it requires, says the same thing on its tile as on its masthead, and is greyed only for
what its own row declares -- the `needs` written on it, or the destination it names failing to build. Adding a
module and forgetting any of that fails the suite, by name.

**What it does not cover, said plainly so that nobody reads more into it than is there.**

1. **Whether a module is any good.** The contract is about shape. That the Map module reads the right columns
   is held by `test_module_shell`, `test_map_reading` and the rest; nothing here replaces them.
2. **Whether the words are the right words.** A tooltip that exists and says nothing passes, and so does a note
   that is accurate and dull. Plain wording is a standing rule and a human reading, not a test.
3. **What a module DOES with a subject it is handed.** The contract holds that a module declaring it reads a
   gene has somewhere to put one. That the Expression module then actually profiles that gene is a test of its
   own, per module.
4. **Anything outside a module.** The preferences, the download dialogs and the menu bar are not catalogue rows
   and are held by `test_app_controls` and `test_app_menus`.
5. **The layout at other sizes.** Everything here is judged at 1920x1080. `test_module_shell` reads the gene
   card at several widths; no other module is checked at more than one.
6. **Map and Network offer no next step.** Gene card, Gene list, Screen, Predict and now Infer do. The two that
   do not hand on a SPACE and a GENE PAIR, which nothing can pre-fill, so the row would be a row of offers all
   saying "none chosen". Left as it is rather than half-built; the rail under each module already lists what
   follows it.
7. **Tables outside a module shell.** `test_every_table_says_what_its_columns_are_before_anything_is_loaded`
   walks the tables inside the module shells. A table in a dialog or a further window is not reached by it.
