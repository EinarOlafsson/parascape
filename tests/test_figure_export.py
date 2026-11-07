"""Figures from the module views: vector SVG and PDF, and a 3x PNG (instruction 60, change 4)."""
from __future__ import annotations

import pytest
from PySide6 import QtCore, QtGui, QtWidgets

from parascape import figure_export as F


class _Dots(QtWidgets.QWidget):
    """A view that paints itself, the way the module views do."""

    def paintEvent(self, _event):                         # noqa: N802
        painter = QtGui.QPainter(self)
        painter.setBrush(QtGui.QColor("#43c6d8"))
        for i in range(20):
            painter.drawEllipse(QtCore.QPointF(10 + 9 * i, 40), 3, 3)
        painter.drawLine(0, 0, 200, 80)
        painter.end()


@pytest.fixture
def view(qtbot):
    w = _Dots()
    w.resize(200, 80)
    qtbot.addWidget(w)
    return w


def test_an_svg_figure_is_vector_not_a_picture(view, tmp_path):
    path = F.save(view, str(tmp_path / "f.svg"))
    text = open(path).read()
    assert "<svg" in text and ("<circle" in text or "<ellipse" in text or "<path" in text)
    assert "<image" not in text


def test_a_pdf_figure_is_written_at_the_size_of_the_view(view, tmp_path):
    path = F.save(view, str(tmp_path / "f.pdf"))
    data = open(path, "rb").read()
    assert data.startswith(b"%PDF")
    import re
    box = re.search(rb"/MediaBox \[([^\]]*)\]", data)
    assert box and [round(float(v)) for v in box.group(1).split()] == [0, 0, 200, 80]


def test_a_png_figure_is_three_times_the_screen(view, tmp_path):
    path = F.save(view, str(tmp_path / "f.png"))
    image = QtGui.QImage(path)
    assert (image.width(), image.height()) == (200 * F.PNG_SCALE, 80 * F.PNG_SCALE)


def test_an_unknown_format_is_refused(view, tmp_path):
    with pytest.raises(ValueError):
        F.save(view, str(tmp_path / "f.jpg"))


def test_a_view_that_offers_a_figure_painter_is_drawn_with_it(view, tmp_path):
    calls = []
    view.paint_figure = lambda painter: calls.append(painter)
    F.save(view, str(tmp_path / "g.svg"))
    assert len(calls) == 1
