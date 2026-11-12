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
