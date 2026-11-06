# 86 · Every module in spaCR's module shell

**Status: open (instruction 74, goals 6-7; requested 2026-11-05).** Einar: "spaCR has its settings in a
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
