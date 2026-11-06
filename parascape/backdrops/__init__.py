"""Animated backdrops from spaCR: the spaCR field and its sibling themes.

ParaScape's home screen and Navigator sit over a live, slowly moving field, the way
spaCR's do, because Einar asked for ParaScape to look like spaCR rather than to be
reminded of it. The animation is spaCR's own code, vendored into
:mod:`parascape.backdrops.engines` nearly verbatim -- possible only because ParaScape
moved to PySide6, so it transfers rather than being rewritten. What ParaScape supplies
in place of the rest of spaCR (preferences, theme lookups, event hubs) is in
:mod:`parascape.backdrops._compat`.

The package is named ``backdrops`` and not ``ambient`` for a reason worth keeping:
ParaScape already has ``parascape/ambient.py``, its own blob backdrop behind the map,
and a package directory silently wins over a same-named module. The first version of
this package was called ``ambient``, shadowed that module, and stopped the main window
opening.

Licensed under spaCR's terms -- LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt.
"""
