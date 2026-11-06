# 85 · Retire the old ParaScape home screen

**Status: open (requested 2026-11-06).** Einar: "the home screen seems to still be the home screen for
ParaScape".

## What is there
Two home screens. The new landing (instruction 74 overhaul: Navigator and module tiles, ported from spaCR)
covers the window at launch. Behind it the OLD home -- the `home` dock of Gene / List / Label node cards,
built before Einar asked for spaCR's -- is still installed and is still the raised tab in the right-hand
dock group, so it is what appears the moment the landing is left.

## The plan
1. The landing is the only home. The old `home` dock is no longer raised at startup and no longer offered as
   a home.
2. Nothing it did is lost. Inventory what only the old dock provides -- the dropped-file target, the list
   view that answers a question for every gene in a dropped list, the gene-list workflows -- and give each a
   home in the new structure (the Drop and List modules) before the dock is removed.
3. Then remove the dock, its tests moved to the new owners of each behaviour.

## Gate
Leaving the landing never shows a second home screen, and every behaviour the old dock had is reachable
from the landing.
