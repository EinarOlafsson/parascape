# 95 · spaCR's dock

**Status: open (requested 2026-11-06).** Einar: "there should also be a dock like there is in spaCR".

## What spaCR does (`spacr/qt/widgets/dock.py`)
A left navigation column: Home, then folding categories, each module a row of icon and name (always both);
hover only recolours; the column folds against its edge.

## The plan
Port `dock.py` (vendored with a `_compat` like the panes) and fill it from `parascape/modules.py`'s bands.
A row opens its module exactly as its tile does. It sits left of the landing and of every module shell.

## Gate
The dock is there on Home and in every module, folds away, and every row opens its module.
