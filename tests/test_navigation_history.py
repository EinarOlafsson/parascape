"""Back and forward through ParaScape (Einar, 2026-11-24): the pure history, and the window."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from parascape.navigation import History  # noqa: E402


def test_back_and_forward_retrace_the_places_visited():
    h = History()
    for place in [("page", "Navigator"), ("module", "network"), ("module", "gene")]:
        assert h.record(place)
    assert h.step(-1) == ("module", "network")
    assert h.step(-1) == ("page", "Navigator")
    assert not h.can_go(-1) and h.step(-1) is None
    assert h.step(1) == ("module", "network")
    assert h.can_go(1)


def test_a_new_place_after_going_back_drops_the_places_ahead():
    h = History()
    for place in [("page", "Navigator"), ("module", "network"), ("module", "gene")]:
        h.record(place)
    h.step(-1)
    h.record(("module", "map"))
    assert h.places == [("page", "Navigator"), ("module", "network"), ("module", "map")]
    assert not h.can_go(1)


def test_the_same_place_twice_is_one_visit_and_replaying_records_nothing():
    h = History()
    h.record(("module", "gene"))
    assert not h.record(("module", "gene"))
    h.replaying = True
    assert not h.record(("module", "map"))
    assert h.places == [("module", "gene")]


@pytest.fixture(scope="module")
def win():
    from parascape import paths
    if not os.path.exists(os.path.join(paths.data_dir(), "graph.npz")):
        pytest.skip("no built cache")
    from PySide6 import QtWidgets
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    from parascape import app as A
    try:
        w = A.Window()
    except Exception as e:                                  # pragma: no cover
        pytest.skip(f"cannot create the window offscreen: {e}")
    w.resize(1400, 900)
    w.show()
    qapp.processEvents()
    yield w
    w.close()
    qapp.processEvents()


def test_back_and_forward_move_through_modules_and_the_main_window(win):
    from parascape import landing_install as LI, navigation as NAV
    LI.open_module(win, "network")
    LI.open_module(win, "gene")
    win.gene_module.module._full_panel()
    assert not LI.landing_up(win)
    assert win.rail_dock.isVisible(), "the tree must stay on screen in the main window"
    assert win.nav_back_action.isEnabled()
    NAV.go(win, -1)
    assert LI.landing_up(win) and NAV.landing_place(win.landing) == ("module", "gene")
    NAV.go(win, -1)
    assert NAV.landing_place(win.landing) == ("module", "network")
    assert win.nav_forward_action.isEnabled()
    NAV.go(win, 1)
    assert NAV.landing_place(win.landing) == ("module", "gene")
    assert win.landing.rail.parent() is win.landing and not win.rail_dock.isVisible()


def test_the_evidence_panel_closes_and_does_not_outlive_home(win):
    from PySide6 import QtWidgets
    from parascape import landing_install as LI
    win.gene_module.module._full_panel()
    assert win.right_dock.isVisible()
    win.right_dock.close()
    assert not win.right_dock.isVisible()
    win.open_gene_evidence(10)
    dialog = win._gene_evidence_dialog
    assert dialog.isVisible()
    LI.go_home(win)
    assert not dialog.isVisible()
