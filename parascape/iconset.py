"""The one place ParaScape gets a symbol from, and the symbols are spaCR's.

Einar asked that the symbols be "exactly the same as spaCR". The way to keep that true is to have
one provider and no other: every module card, dock row and result card asks this file for its
symbol by the name of the module or the kind of result, and nothing else in ParaScape builds a
``QIcon`` from a path. A symbol that was loaded some other way could drift from spaCR without
anyone noticing.

How spaCR does it (``spacr/qt/iconset.py``), and what is kept here:

* **The artwork is flat white line art on transparency.** The shape is the alpha channel and the
  colour is whatever the theme wants, so the artwork is *re-inked* for the theme instead of being
  drawn as it came. White on the dark theme is the artwork unchanged; on a light theme it would be
  invisible, so it is painted in the ink colour of the theme.
* **Re-inked pixmaps are cached** by (file, theme, size), because a card grid asks for the same
  handful of files over and over.
* **A missing symbol is never a blank square.** spaCR falls back to a puzzle piece. Here the
  fallback is the settings symbol, and the miss is logged, and in a test run (or with
  ``PARASCAPE_STRICT_ICONS=1``) it raises :class:`MissingIcon` instead, so a module that has no
  symbol fails the suite rather than being noticed on screen by a person.

No Qt is imported until a pixmap is asked for, so the tables can be read and checked without a
display.
"""
from __future__ import annotations

import logging
import os
import sys
from functools import lru_cache
from pathlib import Path

_log = logging.getLogger("parascape.iconset")

#: Where the spaCR artwork lives inside ParaScape. Copied from spaCR, unchanged, under its licence.
ICON_DIR = Path(__file__).with_name("data") / "icons" / "spacr"

#: What stands in for a symbol that cannot be found, outside strict mode.
FALLBACK = "settings"

#: Setting this to ``1`` makes a missing symbol raise instead of falling back.
ENV_STRICT = "PARASCAPE_STRICT_ICONS"

#: The spaCR artwork used for each ParaScape module, by module key. One symbol per module, and the
#: symbol is the one spaCR uses for the nearest thing it already has.
MODULE_ICON = {
    "home": "home",
    "gene": "profiler", "structure": "layer_viewer", "partners": "graph_builder",
    "expression": "pca", "map": "umap", "network": "graph_builder", "cell": "image_scatter",
    "gallery": "trellis", "questions": "experiment_design", "infer": "ml_analyze",
    "strategies": "model_zoo", "scorecard": "model_compare", "runs": "run_history",
    "fill": "regression", "drop": "download", "screen": "hit_list", "list": "tabulate",
    "export": "methods_export", "datasets": "data_manager", "information": "feature_dict",
    "jobs": "queue", "console": "report",
}

#: The spaCR artwork for a result that is not a module, by the kind of result.
KIND_ICON = {
    "gene": "profiler", "dataset": "db_browser", "question": "experiment_design",
    "strategy": "model_zoo", "guide": "report", "panel": "settings", "node": "settings",
}


class MissingIcon(LookupError):
    """A symbol was asked for that ParaScape does not have.

    Raised in strict mode so that a missing symbol is a test failure and not a quiet empty square.
    """


def strict() -> bool:
    """Whether a missing symbol should raise rather than fall back.

    True under pytest and when ``PARASCAPE_STRICT_ICONS`` is set to anything but ``0``. Strict by
    default in a test run is the point: the person who would otherwise meet the blank square is
    Einar, and the suite should meet it first.
    """
    flag = os.environ.get(ENV_STRICT)
    if flag is not None:
        return flag.strip() not in ("", "0", "false", "no")
    return "pytest" in sys.modules


def _missing(what: str) -> str:
    """Handle a symbol that cannot be found: raise in strict mode, else log and fall back.

    :param what: a sentence saying what was asked for, so the failure names it.
    :returns: the name of the fallback symbol.
    """
    if strict():
        raise MissingIcon(what)
    _log.warning("%s; using the %s symbol instead", what, FALLBACK)
    return FALLBACK


def module_icon_name(key: str) -> str:
    """The spaCR artwork name for a module, or the fallback (or an error) if it has none.

    :param key: the module key, as in the module catalogue.
    """
    name = MODULE_ICON.get(key)
    return name if name else _missing(f"module {key!r} has no spaCR symbol in MODULE_ICON")


def kind_icon_name(kind: str) -> str:
    """The spaCR artwork name for a kind of result (a gene, a dataset, a question and so on).

    :param kind: the kind of result.
    """
    name = KIND_ICON.get(kind)
    return name if name else _missing(f"result kind {kind!r} has no spaCR symbol in KIND_ICON")


def icon_file(name: str) -> str:
    """The path of one piece of spaCR artwork, checked to exist.

    :param name: the file name without ``.png``, as spaCR names it.
    :returns: the path; the fallback artwork (or an error in strict mode) when the file is absent.
    """
    path = ICON_DIR / f"{name}.png"
    if path.is_file():
        return str(path)
    fallback = _missing(f"the spaCR artwork {name}.png is not in {ICON_DIR}")
    return str(ICON_DIR / f"{fallback}.png")


def module_icon_file(key: str) -> str:
    """The artwork path for a module.

    :param key: the module key.
    """
    return icon_file(module_icon_name(key))


def kind_icon_file(kind: str) -> str:
    """The artwork path for a kind of result.

    :param kind: the kind of result.
    """
    return icon_file(kind_icon_name(kind))


def ink_for(theme: str = "dark") -> str:
    """The colour the artwork is painted in for a theme: the ink of the theme.

    :param theme: a theme name from :mod:`parascape.theme`.
    """
    from . import theme as TH
    return TH.palette_for(theme)["fg"]


@lru_cache(maxsize=256)
def _inked(path: str, ink: str, px: int):
    """The artwork at ``px`` square, with its shape filled in ``ink``, kept for next time.

    The shape of the artwork is its alpha channel, so painting the ink through that channel
    re-inks every symbol however it was exported. Cached because a card grid asks for the same few
    files at the same size again and again, and each decode of a 1024 pixel PNG is not free.

    :param path: the artwork file.
    :param ink: the colour to fill the shape with, as ``#rrggbb``.
    :param px: the edge of the square to return, in physical pixels.
    """
    from PySide6 import QtCore, QtGui
    source = QtGui.QImage(path)
    if source.isNull():
        return QtGui.QPixmap()
    source = source.convertToFormat(QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    source = source.scaled(px, px, QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                           QtCore.Qt.TransformationMode.SmoothTransformation)
    out = QtGui.QImage(source.size(), QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    out.fill(QtGui.QColor(ink))
    painter = QtGui.QPainter(out)
    painter.setCompositionMode(QtGui.QPainter.CompositionMode.CompositionMode_DestinationIn)
    painter.drawImage(0, 0, source)
    painter.end()
    return QtGui.QPixmap.fromImage(out)


def pixmap(path: str, px: int, theme: str = "dark", ratio: float = 1.0):
    """A re-inked pixmap of one piece of artwork, ready to put in a label or a card.

    :param path: the artwork file, from :func:`icon_file` and its relatives.
    :param px: the edge in device-independent pixels.
    :param theme: the theme to ink it for.
    :param ratio: the device pixel ratio of the screen, so the symbol stays sharp on a dense one.
    """
    pm = _inked(path, ink_for(theme), max(1, int(px * ratio))).copy()
    pm.setDevicePixelRatio(ratio)
    return pm


def icon(path: str, theme: str = "dark"):
    """A ``QIcon`` of one piece of artwork, inked for the theme, at the sizes a dock row asks for.

    :param path: the artwork file.
    :param theme: the theme to ink it for.
    """
    from PySide6 import QtGui
    result = QtGui.QIcon()
    for px in (16, 24, 32, 48, 64):
        result.addPixmap(pixmap(path, px, theme))
    return result


def module_icon(key: str, theme: str = "dark"):
    """The ``QIcon`` for a module. The one call every dock row and tile goes through.

    :param key: the module key.
    :param theme: the theme to ink it for.
    """
    return icon(module_icon_file(key), theme)


def kind_icon(kind: str, theme: str = "dark"):
    """The ``QIcon`` for a kind of result.

    :param kind: the kind of result.
    :param theme: the theme to ink it for.
    """
    return icon(kind_icon_file(kind), theme)


def problems(module_keys) -> list:
    """Every way the symbols are wrong for a list of modules, as plain sentences.

    The test of the module catalogue calls this with every module key, so a module added without
    a symbol fails the suite. Never raises: it returns what is wrong so one run reports all of it.

    :param module_keys: the keys of every module that must have a symbol.
    :returns: an empty list when every module and every kind has a symbol file.
    """
    wrong = []
    for key in module_keys:
        name = MODULE_ICON.get(key)
        if not name:
            wrong.append(f"module {key!r} has no spaCR symbol")
        elif not (ICON_DIR / f"{name}.png").is_file():
            wrong.append(f"module {key!r} names {name}.png, which is not shipped")
    for kind, name in KIND_ICON.items():
        if not (ICON_DIR / f"{name}.png").is_file():
            wrong.append(f"result kind {kind!r} names {name}.png, which is not shipped")
    if not (ICON_DIR / f"{FALLBACK}.png").is_file():
        wrong.append(f"the fallback {FALLBACK}.png is not shipped")
    return wrong
