"""Save a module view as a figure: SVG or PDF as vectors, PNG at print resolution (instruction 60, change 4).

A view that draws itself with a painter (the Map scatter, the Network view) is drawn again onto the
export device, so the SVG and PDF are true vector figures -- points and lines a journal can scale -- and
not a screenshot wrapped in a vector file. A view may offer `paint_figure(painter)` to draw a figure
differently from the screen (the Network view does: on screen its background is a cached image).
"""
from __future__ import annotations

import os

from PySide6 import QtCore, QtGui, QtWidgets

#: The formats offered, as the file dialog filter.
FILTER = "SVG figure (*.svg);;PDF figure (*.pdf);;PNG image, 3x (*.png)"
PNG_SCALE = 3


def _background(widget) -> QtGui.QColor:
    colour = widget.palette().color(QtGui.QPalette.ColorRole.Window)
    try:
        from . import theme as TH
        from .shell import theme_name
        colour = QtGui.QColor(TH.palette_for(theme_name()).get("page", colour.name()))
    except Exception:
        pass
    return colour


def _draw(widget, painter, size) -> None:
    painter.fillRect(QtCore.QRectF(0, 0, size.width(), size.height()), _background(widget))
    if hasattr(widget, "paint_figure"):
        widget.paint_figure(painter)
    else:
        widget.render(painter, QtCore.QPoint(), QtGui.QRegion(),
                      QtWidgets.QWidget.RenderFlag.DrawChildren)


def save(widget: QtWidgets.QWidget, path: str, title: str = "ParaScape figure") -> str:
    """Write `widget` to `path`; the suffix chooses SVG, PDF or PNG. Returns the path."""
    suffix = os.path.splitext(path)[1].lower()
    size = widget.size()
    if suffix == ".svg":
        from PySide6.QtSvg import QSvgGenerator
        generator = QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(size)
        generator.setViewBox(QtCore.QRect(0, 0, size.width(), size.height()))
        generator.setTitle(title)
        generator.setDescription("Exported from ParaScape")
        painter = QtGui.QPainter(generator)
        _draw(widget, painter, size)
        painter.end()
    elif suffix == ".pdf":
        writer = QtGui.QPdfWriter(path)
        writer.setTitle(title)
        writer.setCreator("ParaScape")
        writer.setResolution(72)
        writer.setPageSize(QtGui.QPageSize(QtCore.QSizeF(size.width(), size.height()),
                                           QtGui.QPageSize.Unit.Point))
        writer.setPageMargins(QtCore.QMarginsF(0, 0, 0, 0))
        painter = QtGui.QPainter(writer)
        _draw(widget, painter, size)
        painter.end()
    elif suffix == ".png":
        image = QtGui.QImage(size.width() * PNG_SCALE, size.height() * PNG_SCALE,
                             QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        image.setDevicePixelRatio(PNG_SCALE)
        painter = QtGui.QPainter(image)
        _draw(widget, painter, size)
        painter.end()
        if not image.save(path):
            raise OSError(f"could not write {path}")
    else:
        raise ValueError(f"choose .svg, .pdf or .png, not {suffix or 'no suffix'}")
    return path


def ask_and_save(widget: QtWidgets.QWidget, parent=None, name: str = "figure",
                 title: str = "ParaScape figure") -> str:
    """A save dialog, then `save`. Returns the path written, or "" when cancelled."""
    path, chosen = QtWidgets.QFileDialog.getSaveFileName(parent or widget, "Save the figure",
                                                         name + ".svg", FILTER)
    if not path:
        return ""
    if not os.path.splitext(path)[1]:
        path += {"PDF": ".pdf", "PNG": ".png"}.get(chosen.split()[0], ".svg")
    return save(widget, path, title)
