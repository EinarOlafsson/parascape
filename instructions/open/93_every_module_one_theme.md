# 93 · Every module in the same spaCR theme: fields, translucent containers, the gene card

**Status: open (requested 2026-11-06).** Einar: "the fields should look more like spaCR and the containers
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
