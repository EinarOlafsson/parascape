"""The landing view: the Navigator and the module tiles, over one animated field.

PORTED AND MODIFIED FROM spaCR
==============================
spaCR keeps its Navigator and its tile Home as two pages of one stack over an
"unchanged ambient field": switching page swaps what floats over the animation and
leaves the animation itself running. That is the structure here.

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

**Modifications made in porting.** spaCR switches pages from its sidebar; ParaScape
has no sidebar yet, so a two-way switch sits at the top right. The backdrop is
installed once, on this widget, and both pages go transparent over it -- installing
one per page would restart the animation on every switch and run two timers for
one picture.
"""
from __future__ import annotations

from PySide6 import QtCore, QtWidgets

from . import iconset
from . import theme as TH
from .home import HomeScreen
from .navigator import NavigatorPage
from .tiles import tile_qss


class Landing(QtWidgets.QWidget):
    """Navigator and Modules, one switch between them, one backdrop under both.

    `module_chosen(key)` fires for a module, from either page. `item_chosen(kind, key)`
    fires for anything else the Navigator offers -- a gene, a question, a dataset.
    """

    module_chosen = QtCore.Signal(str)
    item_chosen = QtCore.Signal(str, str)

    PAGES = ("Navigator", "Modules")

    def __init__(self, counts=None, reasons=None, genes=None, index=None,
                 animate: bool = True, start: str = "Navigator", parent=None):
        super().__init__(parent)
        self.setObjectName("Landing")
        palette = TH.palette_for("dark")
        self.setStyleSheet(tile_qss(palette) + _switch_qss(palette)
                           + f"QWidget#Landing {{ background-color: {palette['bg']}; }}"
                           + "QWidget#Landing[ambient=\"on\"] { background: transparent; }")

        # spaCR's dock down the left (instruction 95): Home, then every module under its
        # band, on Home and beside every module shell alike. Folds against its edge.
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(0)
        self.dock, self.dock_edge = self._build_dock()
        if self.dock is not None:
            row.addWidget(self.dock)
            row.addWidget(self.dock_edge)
        body = QtWidgets.QWidget()
        body.setObjectName("LandingBody")
        body.setStyleSheet("QWidget#LandingBody { background: transparent; }")
        row.addWidget(body, 1)
        col = QtWidgets.QVBoxLayout(body)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(0)

        bar = QtWidgets.QHBoxLayout()
        bar.setContentsMargins(TH.SPACING["xl"], TH.SPACING["md"], TH.SPACING["xl"], 0)
        bar.addStretch(1)
        self._switch = QtWidgets.QButtonGroup(self)
        self._switch.setExclusive(True)
        for i, name in enumerate(self.PAGES):
            button = QtWidgets.QPushButton(name)
            button.setObjectName("LandingSwitch")
            button.setCheckable(True)
            button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            button.setProperty("position", "first" if i == 0 else "last")
            button.setToolTip(
                "One field for anything: a gene, a product, a question, a module, a dataset "
                "or a method, answered as you type." if name == "Navigator" else
                "Every module as a tile, grouped by what you came to do, with the ones that "
                "cannot open yet greyed and the reason on them.")
            self._switch.addButton(button, i)
            bar.addWidget(button)
        col.addLayout(bar)

        self.stack = QtWidgets.QStackedWidget()
        self.stack.setObjectName("LandingStack")
        self.navigator = NavigatorPage(genes=genes, index=index)
        self.home = HomeScreen(counts=counts, reasons=reasons, animate=False)
        self.stack.addWidget(self.navigator)
        self.stack.addWidget(self.home)
        col.addWidget(self.stack, 1)

        self._switch.idClicked.connect(self.show_page_index)
        self.navigator.chosen.connect(self._from_navigator)
        self.home.module_chosen.connect(self.module_chosen)

        self._ambient = None
        if animate:
            self._install_ambient()
        if not self.navigator_enabled():
            for button in self._switch.buttons():
                button.hide()
            start = "Modules"
        self.show_page(start)

    # -------------------------------------------------------------------- dock
    def _build_dock(self):
        """spaCR's dock, filled from `modules.py`. (None, None) if it cannot be built -- the
        landing must open without it."""
        try:
            from . import modules as MOD
            from .panes.dock import Dock, DockEdge
        except Exception:
            return None, None
        rows = [("home", "Home", "the Navigator and the module tiles", "Start")]
        rows += [(m.key, m.name, m.note, m.band) for m in MOD.MODULES]
        # Every symbol comes from the one provider, and is asked for here, outside the guarded
        # block below: a module with no symbol must fail loudly in a test, not silently leave the
        # dock out.
        icons = {key: iconset.module_icon(key) for key, _n, _d, _b in rows}

        def icon_for(key):
            """The symbol of a dock row, from the one icon provider."""
            return icons.get(key)
        try:
            dock = Dock(rows, icon_for=icon_for)
            edge = DockEdge(dock)
        except Exception:
            return None, None
        # spaCR explains a row in a strip along the window's foot rather than a tooltip;
        # ParaScape has no strip, and every control here says what it does on hover.
        notes = {key: (name, note) for key, name, note, _band in rows}
        for row in dock._rows:
            name, note = notes.get(row.key, (row.key, ""))
            row.setToolTip(f"{name}: {note}." if note else name)
        dock.nav_selected.connect(self._from_dock)
        dock.setToolTip("Every module. Click one to open it; drag the edge to resize, click the "
                        "edge to hide the dock.")
        return dock, edge

    def _from_dock(self, key: str) -> None:
        if key == "home":
            self.show_page(self.home_page())
        else:
            self.module_chosen.emit(key)

    # ------------------------------------------------------------------- pages
    def navigator_enabled(self) -> bool:
        """Whether the Navigator is on. On by default; remembered once changed."""
        from .backdrops import _compat as AC
        return AC._as_bool(AC._settings().value("landing/navigator", True), True)

    def set_navigator_enabled(self, on: bool) -> None:
        """Turn the Navigator on or off, as spaCR's Help ▸ Navigator does.

        Off, the landing is the module tiles alone: the Navigator page and the two-way
        switch go, and the tiles are what opens. On, the Navigator is the first page again.
        """
        from .backdrops import _compat as AC
        AC._settings().setValue("landing/navigator", bool(on))
        for button in self._switch.buttons():
            button.setVisible(bool(on))
        self.show_page("Navigator" if on else "Modules")

    def show_page(self, name: str) -> None:
        """Bring `Navigator` or `Modules` to the front."""
        self.show_page_index(self.PAGES.index(name))

    def show_page_index(self, i: int) -> None:
        """Bring page `i` forward and check its switch; the Navigator takes focus so
        typing works the moment it appears."""
        self.stack.setCurrentIndex(i)
        self._switch.button(i).setChecked(True)
        on = self.navigator_enabled()
        for button in self._switch.buttons():
            button.setVisible(on)
        if self.stack.currentWidget() is self.navigator:
            self.navigator.reset()
            self.navigator.field.setFocus()

    def show_module(self, shell: QtWidgets.QWidget) -> None:
        """Open a module in spaCR's module shell, over the same backdrop (instruction 86).

        The shell is a page of this stack, as a spaCR module is a page of spaCR's window:
        the animation keeps running behind it, and Home comes back to the tiles. The
        Navigator/Modules switch is hidden while a module is up -- it switches home pages.
        """
        if self.stack.indexOf(shell) < 0:
            self.stack.addWidget(shell)
        self.stack.setCurrentWidget(shell)
        for button in self._switch.buttons():
            button.hide()

    @property
    def page(self) -> str:
        i = self.stack.currentIndex()
        return self.PAGES[i] if i < len(self.PAGES) else "Module"

    def home_page(self) -> str:
        """The page Home returns to: the Navigator when it is on, else the tiles."""
        return "Navigator" if self.navigator_enabled() else "Modules"

    def _from_navigator(self, kind: str, key: str) -> None:
        if kind == "module":
            self.module_chosen.emit(key)
        else:
            self.item_chosen.emit(kind, key)

    # --------------------------------------------------------------- animation
    def _install_ambient(self) -> None:
        """One backdrop for both pages, installed here so a page switch never restarts it."""
        try:
            from .backdrops import _compat as AC
            from .backdrops import engines as AE
        except Exception as exc:
            self._ambient_error = f"{type(exc).__name__}: {exc}"
            return
        if not AC.get_ambient_enabled():
            return
        try:
            self._ambient = AE.install_ambient(
                self, theme=AC.get_ambient_theme(), palette=AC.get_ambient_palette(),
                background=AC.active_page_colour())
        except Exception as exc:
            self._ambient = None
            self._ambient_error = f"{type(exc).__name__}: {exc}"
            return
        self.setProperty("ambient", "on")
        self.stack.setAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.stack.setStyleSheet("QStackedWidget#LandingStack { background: transparent; }")
        self.home._clear_page_surfaces()
        self.style().unpolish(self)
        self.style().polish(self)

    @property
    def animated(self) -> bool:
        return self._ambient is not None

    def set_theme(self, theme: str) -> None:
        """Switch the animation, live, and remember the choice."""
        from .backdrops import _compat as AC
        AC._settings().setValue("ambient/theme", theme)
        if self._ambient is not None:
            self._ambient.hide()
            self._ambient.setParent(None)
            self._ambient.deleteLater()
            self._ambient = None
        self._install_ambient()
        if self._ambient is None:
            # Off, or failed: the containers were transparent for the backdrop, and with
            # nothing behind them they would show through to whatever is underneath.
            self.setProperty("ambient", "off")
            self.stack.setAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground, False)
            self.stack.setStyleSheet("")
            self.home._restore_page_surfaces()
            self.style().unpolish(self)
            self.style().polish(self)


def _switch_qss(P: dict) -> str:
    """The two-way switch: a segmented pill, accent when on."""
    return f"""
QPushButton#LandingSwitch {{
    color: {P["fg_muted"]}; background: {P["surface"]};
    border: 1px solid {P["border_soft"]}; padding: 5px 16px; font-size: 12px;
    min-height: 0px;
}}
QPushButton#LandingSwitch[position="first"] {{
    border-top-left-radius: 13px; border-bottom-left-radius: 13px; border-right: none;
}}
QPushButton#LandingSwitch[position="last"] {{
    border-top-right-radius: 13px; border-bottom-right-radius: 13px;
}}
QPushButton#LandingSwitch:checked {{ color: {P["bg"]}; background: {P["accent"]};
    border-color: {P["accent"]}; }}
QPushButton#LandingSwitch:hover:!checked {{ color: {P["fg"]}; }}
"""
