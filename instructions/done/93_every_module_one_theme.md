# 93 · Every module in the same spaCR theme: fields, translucent containers, the gene card

**Status: built 2026-11-06, waiting on Einar** (every module now in the one shell theme -- see 86). Done: View ▸ Container opacity (30-100 %, default 60 % =
Einar's spaCR setting, key `prefs/pane_opacity`); every container in a module shell at that opacity with the
animation showing through; spaCR's input rules and its field fade (vendored `backdrops/field_fade.py`, scoped
to module shells by `parascapeFieldFade`). Gene module built 2026-11-06 with the
skill (`parascape/gene_module.py`, `parascape/gene_reading.py`): left -- your gene (field, list, drop,
stepping) and what to read (the information categories as plain questions); right -- **What is known** (one
sentence per kind of knowledge, each with linked sources: where it is, what it binds, how much the parasite
needs it, when its RNA is made...), **The evidence** grouped the same way, **What is not known, and what can
be inferred** (key traits missing, predicted cells kept apart with their certainty, and what its Map
neighbours suggest with how often that is right against chance). Tests: `tests/test_gene_reading.py`,
`tests/test_module_shell.py`. Remaining: the other modules through the skill (instruction 86).

Requested 2026-11-06. Einar: "the fields should look more like spaCR and the containers
should have opacity like spaCR's do ... the gene card is not what I want. All modules need to be thematically
the same."

## What spaCR does
- Inputs: `surface_alt` fill, 1 px `border`, small radius, accent border on focus, transparent spin buttons
  (`spacr/qt/theme.py`, "Inputs" block), and the field fade (`backdrops/field_fade.py`, already vendored).
- Containers: every card, settings section, console and preview surface painted at the user's **pane
  opacity** (`prefs/pane_opacity`; Einar's spaCR is at 60 %), so the animation shows through.

## The plan
1. ParaScape gets the same preference (View ▸ Container opacity, default 60 %) and the shell paints every
   container at it.
2. spaCR's input rules in the shell stylesheet.
3. Every module through the same shell -- the Gene module first (the gene card becomes the Gene module:
   gene search and settings on the left, its evidence sections on the right), then the remaining ones in the
   order of instruction 86.

## Gate
Map, Gene and every later module look like one application, and like spaCR beside it.
