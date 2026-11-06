"""spaCR's pane widgets: the collapsible splitter, sections and fold headers.

What makes a spaCR module screen give its room to whatever is being worked on:
every pane resizes by its edge, folds to a header at the limit, and remembers how
it was left. Vendored from ``spacr/qt/widgets/`` for :mod:`parascape.shell`, the
module shell of instruction 86; what they asked of the rest of spaCR is answered in
:mod:`parascape.panes._compat`.

Licensed under spaCR's terms -- LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt.
"""
