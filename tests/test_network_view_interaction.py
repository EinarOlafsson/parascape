"""The network view: a shift+drag box selects genes, a click between genes picks the nearest
bright link, and the camera transform agrees with `to_screen` (Einar, 2026-11-24)."""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

QtWidgets = pytest.importorskip("PySide6.QtWidgets")
from PySide6 import QtCore  # noqa: E402


@pytest.fixture(scope="module")
def view():
    QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from parascape.network_module import NetworkView
    v = NetworkView()
    v.resize(400, 400)
    xy = np.array([[-0.5, 0.0], [0.5, 0.0], [0.0, 0.5], [0.9, 0.9]])
    links = pd.DataFrame({"a": [0, 0], "b": [1, 2], "source": ["xlms", "cofitness"],
                          "strength": [0.9, 0.5]})
    v.set_data(xy, [(255, 0, 0, 255)] * 4, links, np.array([0, 1, 1, np.inf]), [0])
    return v


def test_transform_matches_to_screen(view):
    p = view.to_screen(view.xy)
    for i, (x, y) in enumerate(view.xy):
        q = view.transform().map(QtCore.QPointF(x, y))
        assert abs(q.x() - p[i, 0]) < 1e-6 and abs(q.y() - p[i, 1]) < 1e-6


def test_a_box_selects_the_genes_inside_it(view):
    p = view.to_screen(view.xy)
    rect = QtCore.QRect(int(p[:3, 0].min()) - 2, int(p[:3, 1].min()) - 2,
                        int(np.ptp(p[:3, 0])) + 4, int(np.ptp(p[:3, 1])) + 4)
    assert view.genes_in(rect) == [0, 1, 2]
    assert view.genes_in(QtCore.QRect(0, 0, 1, 1)) == []


def test_a_click_on_a_link_picks_it_and_off_links_picks_nothing(view):
    p = view.to_screen(view.xy)
    mid = QtCore.QPointF(*(p[0] + p[1]) / 2)
    assert view.link_at(mid) == (0, 1)
    assert view.link_at(QtCore.QPointF(5, 5)) is None


def test_faint_genes_are_built_once_as_shapes(view):
    assert view._far_dots, "the genes off the focus are drawn as one shape per colour"
    assert set(view._line_paths) == {"xlms", "cofitness"}
