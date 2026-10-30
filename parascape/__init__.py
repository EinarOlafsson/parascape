"""Explore gene evidence in Toxoplasma gondii and Plasmodium falciparum.

Use the analysis modules from Python, ``parascape-discover`` for batch searches,
or ``parascape`` for the desktop browser. The user guide and API examples are
available at https://einarolafsson.github.io/parascape/.

The main modules load on first use, so ``import parascape`` stays light and never imports Qt::

    import parascape
    parascape.strategies.overview("Tg")          # every strategy, its task, grade and skill
    parascape.strategies.test("stacking", target="compartment").card()
    parascape.scorecard.glossary("ranking")      # what each metric means
"""
__version__ = "0.0.5.0"

#: Submodules reachable as attributes of the package, imported when first touched.
_LAZY = ("strategies", "scorecard", "techniques", "calibration", "graphspace", "deposits",
         "datasets", "slots", "paths", "search", "leakage", "organisms", "packs", "spaces")


def __getattr__(name: str):
    if name in _LAZY:
        import importlib
        module = importlib.import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module 'parascape' has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(_LAZY))
