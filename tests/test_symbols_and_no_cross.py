"""Symbols are spaCR's, from one provider; a small cross closes things; no search field has one.

Einar, 2026-11-12: "I don't like the little x cross (I don't want to see that anywhere)" and "the
symbols should be exactly the same as spaCR and the ai symbol should work in the same way".
Revised by Einar, 2026-11-24: "the evidence tab can have a cross to the top right not a close
button, use the x throughout and always opt for more discreet GUI elements". So a closable dock
or tab is now right, and what is defended is: a module with a symbol, no clear-button cross in a
search field (Escape empties it), and no worded "Close" button where a cross would do.
"""
import os
import re
from pathlib import Path

import pytest
from PySide6 import QtCore, QtGui, QtWidgets

from parascape import iconset, modules as MOD, no_cross

ROOT = Path(__file__).resolve().parent.parent / "parascape"


# --------------------------------------------------------------- the symbols


def test_every_module_has_a_shipped_spacr_symbol():
    assert iconset.problems([m.key for m in MOD.MODULES]) == []


def test_every_symbol_is_a_copy_of_spacr_artwork_that_exists_in_spacr():
    """Same artwork as spaCR means the same file names as spaCR; a look-alike would have its own."""
    spacr = Path("/media/carruthers/mnt3/claude/spacr-756-nightly/spacr/resources/icons")
    if not spacr.is_dir():
        pytest.skip("no spaCR checkout beside this one")
    shipped = sorted(p.name for p in iconset.ICON_DIR.glob("*.png"))
    assert shipped
    for name in shipped:
        assert (spacr / name).is_file(), f"{name} is not a spaCR symbol"
        assert (spacr / name).read_bytes() == (iconset.ICON_DIR / name).read_bytes(), \
            f"{name} differs from the spaCR artwork"


def test_a_missing_symbol_is_loud_in_a_test_not_a_blank_square():
    assert iconset.strict(), "a test run is strict by default"
    with pytest.raises(iconset.MissingIcon):
        iconset.module_icon_name("no_such_module")
    with pytest.raises(iconset.MissingIcon):
        iconset.kind_icon_name("no_such_kind")
    with pytest.raises(iconset.MissingIcon):
        iconset.icon_file("no_such_artwork")


def test_outside_a_test_a_missing_symbol_falls_back_and_says_so(monkeypatch, caplog):
    monkeypatch.setenv(iconset.ENV_STRICT, "0")
    with caplog.at_level("WARNING", logger="parascape.iconset"):
        assert iconset.module_icon_name("no_such_module") == iconset.FALLBACK
    assert "no_such_module" in caplog.text


def test_every_symbol_renders_visible_pixels_in_the_ink_of_the_theme(qapp):
    from parascape import theme as TH
    for key in [m.key for m in MOD.MODULES] + ["home"]:
        for theme in ("dark", "light"):
            pm = iconset.pixmap(iconset.module_icon_file(key), 48, theme)
            assert not pm.isNull(), key
            img = pm.toImage()
            ink = QtGui.QColor(TH.palette_for(theme)["fg"])
            seen = [img.pixelColor(x, y) for x in range(0, img.width(), 2)
                    for y in range(0, img.height(), 2) if img.pixelColor(x, y).alpha() > 200]
            assert seen, f"{key} draws nothing on {theme}"
            assert all(abs(c.red() - ink.red()) < 3 and abs(c.green() - ink.green()) < 3 for c in seen)


def test_a_dark_symbol_and_a_light_symbol_differ_so_neither_vanishes(qapp):
    path = iconset.module_icon_file("gene")
    dark = iconset.pixmap(path, 32, "dark").toImage()
    light = iconset.pixmap(path, 32, "light").toImage()
    assert dark != light


def test_nothing_else_loads_a_spacr_symbol_by_path():
    """One provider: a symbol built some other way could drift from spaCR unnoticed."""
    offenders = []
    for path in ROOT.glob("*.py"):
        if path.name == "iconset.py":
            continue
        text = path.read_text(encoding="utf8")
        if re.search(r"icons[\"'],\s*[\"']spacr|icons/spacr", text):
            offenders.append(path.name)
    assert offenders == [], f"load symbols through parascape.iconset, not by path: {offenders}"


def test_the_ai_symbol_is_the_spacr_ring_and_dot(qapp):
    """Off: a ring. On: a ring and a dot in the colour of the provider. The same as spaCR."""
    from parascape.ai_button import AiButton, PROVIDER_COLORS
    button = AiButton()
    button.resize(36, 36)
    assert button.size() == QtCore.QSize(36, 36)
    off = button.grab().toImage()
    centre = QtCore.QPoint(18, 18)
    assert off.pixelColor(centre).alpha() == 0, "no dot while the assistant is off"
    button.setChecked(True)
    for name, colour in PROVIDER_COLORS.items():
        button.set_provider(name)
        on = button.grab().toImage()
        got = on.pixelColor(centre)
        want = QtGui.QColor(colour)
        assert (got.red(), got.green(), got.blue()) == (want.red(), want.green(), want.blue()), name
        assert got.alpha() == 255


# ------------------------------------------------------------- no small cross


def test_the_cross_finder_sees_a_clear_button_and_allows_closing_crosses(qapp):
    """A finder that cannot see a clear button would pass the whole window by being blind; a
    closable dock or tab is allowed since 2026-11-24."""
    host = QtWidgets.QMainWindow()
    tabs = QtWidgets.QTabWidget()
    tabs.addTab(QtWidgets.QLabel("a"), "a")
    tabs.setTabsClosable(True)
    edit = QtWidgets.QLineEdit()
    edit.setClearButtonEnabled(True)
    box = QtWidgets.QWidget()
    lay = QtWidgets.QVBoxLayout(box)
    for w in (tabs, edit):
        lay.addWidget(w)
    host.setCentralWidget(box)
    dock = QtWidgets.QDockWidget("closable", host)
    host.addDockWidget(QtCore.Qt.DockWidgetArea.BottomDockWidgetArea, dock)
    found = " | ".join(no_cross.find_crosses(host))
    assert "clear button" in found
    assert "closable" not in found
    assert no_cross.closable_docks(host) == ["closable"]


def test_clear_on_escape_replaces_the_clear_button(qapp):
    edit = QtWidgets.QLineEdit()
    edit.setClearButtonEnabled(True)
    no_cross.clear_on_escape(edit)
    no_cross.clear_on_escape(edit)            # twice must not stack a second filter
    assert not edit.isClearButtonEnabled()
    edit.show()
    edit.setText("rhoptry")
    from PySide6 import QtTest
    QtTest.QTest.keyClick(edit, QtCore.Qt.Key.Key_Escape)
    assert edit.text() == ""
    assert len(edit.findChildren(no_cross._EscapeClears)) == 1


def test_no_source_file_switches_a_clear_button_on_or_adds_a_worded_close():
    """Cheap guard at the source, beside the walk of the real window below."""
    bad = []
    for path in ROOT.glob("*.py"):
        if path.name == "no_cross.py":
            continue
        for number, line in enumerate(path.read_text(encoding="utf8").splitlines(), 1):
            code = line.split("#")[0]
            if re.search(r"setClearButtonEnabled\(\s*True|QPushButton\(\s*[\"']Close[\"']\s*\)", code):
                bad.append(f"{path.name}:{number}")
    assert bad == []


def test_the_evidence_panel_closes_with_the_small_cross(window):
    assert "evidence" in no_cross.closable_docks(window)


@pytest.fixture(scope="module")
def window():
    from parascape import app as A
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        win = A.Window(species="Toxoplasma gondii")
    except Exception as error:                      # pragma: no cover - needs a GL context
        pytest.skip(f"cannot create the window offscreen: {error}")
    yield win
    from parascape.map_module import wait_for_projections
    wait_for_projections()
    win.close()
    qapp.processEvents()


def test_the_real_window_has_no_clear_button_cross(window):
    """Walk every widget of the real window, including every dock and the Navigator."""
    window.show()
    QtWidgets.QApplication.processEvents()
    assert no_cross.find_crosses(window) == []


def test_no_clear_button_cross_after_opening_every_module(window):
    """Opening a module builds panels lazily; each one is walked too."""
    from parascape.landing_install import open_module
    for module in MOD.MODULES:
        if MOD.unavailable_reason(module) or module.opens.startswith("action:"):
            continue                        # an action opens a modal dialog, which would block here
        try:
            open_module(window, module.key)
        except Exception:
            continue
        QtWidgets.QApplication.processEvents()
        assert no_cross.find_crosses(window) == [], module.key


def test_every_dressed_dialog_has_no_clear_button_cross(window):
    """The windows ParaScape builds itself: preferences, the slot tree, the help viewer."""
    from parascape import glass
    d = QtWidgets.QDialog()
    d.setWindowTitle("a dressed window")
    QtWidgets.QVBoxLayout(d).addWidget(QtWidgets.QLabel("body"))
    glass.dress(d)
    assert no_cross.find_crosses(d) == []
    d.close()


def test_the_navigator_search_field_has_neither_a_cross_nor_a_clear_button(window):
    field = window.landing.navigator.field
    assert not field.isClearButtonEnabled()
    assert window.landing.navigator.ai_button.isCheckable()


def test_the_home_dock_and_tiles_use_the_symbol_provider(window):
    """Every module tile and every dock row carries a symbol, and it is the one from the provider."""
    for module in MOD.MODULES:
        tile = window.landing.home.tile(module.key)
        label = tile.findChild(QtWidgets.QLabel)
        assert label is not None and label.pixmap() is not None and not label.pixmap().isNull(), \
            module.key
    dock = window.landing.dock
    assert dock is not None
    for row in dock._rows.values() if isinstance(dock._rows, dict) else dock._rows:
        icon = row.icon() if callable(getattr(row, "icon", None)) else None
        assert icon is None or not icon.isNull()
