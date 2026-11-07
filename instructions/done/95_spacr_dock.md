# 95 · spaCR's dock

**Status: built 2026-11-06.** spaCR's `dock.py` and `eliding.py` vendored into `parascape/panes/`
(imports pointed at `panes/_compat.py`, title "ParaScape"); `Landing._build_dock` fills it from
`modules.py` (a Start ▸ Home row, then each band), with spaCR's icons and each module's note as its
tooltip. It sits left of the landing, so it is beside Home and every module shell; the edge drags its
width (remembered) and clicking the edge folds it away. Test:
`test_modules_in_window.py::test_the_dock_lists_home_and_every_module_and_opens_them`.

Requested 2026-11-06. Einar: "there should also be a dock like there is in spaCR".

## What spaCR does (`spacr/qt/widgets/dock.py`)
A left navigation column: Home, then folding categories, each module a row of icon and name (always both);
hover only recolours; the column folds against its edge.

## The plan
Port `dock.py` (vendored with a `_compat` like the panes) and fill it from `parascape/modules.py`'s bands.
A row opens its module exactly as its tile does. It sits left of the landing and of every module shell.

## Gate
The dock is there on Home and in every module, folds away, and every row opens its module.
