"""The home screen: a masthead, bands of module tiles, and a column of numbers.

PORTED AND MODIFIED FROM spaCR
==============================
The layout is ``spacr/qt/widgets/home.py`` (``HomePage``): masthead across the
top, tabs under it, a flowing grid of tiles in bands on the left, a fixed-width
aside of panels on the right.

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

**Modifications made in porting**: PySide6 against ParaScape's palette and
tokens; spaCR's news, queue and stage-legend panels replaced by panels that
answer ParaScape's questions (how much is measured, how much is inferred, what
is running); the tabs are ParaScape's bands; and the tiles flow rather than sit
on a fixed grid, because ParaScape's bands are uneven in size where spaCR's
pipeline stages are not.

What is deliberately kept from spaCR: the tile proportion, the single tile size
across every band, the absence of a resting rim, the fixed aside width, and the
rule that a module which cannot be opened is **greyed with its reason** rather
than hidden.

The UMAP and the parasite drawing are **not** here. They were the old landing
view; they are now the `Map` and `Cell` modules, reached from the Explore band,
which is what Einar asked for -- they are the best thing in the application and
they are not a front door.
"""
from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from . import modules as MOD
from . import theme as TH
from .tiles import ASIDE_W, TILE_MAX_W, AppTile, HomePanel, tile_qss


class _FlowLayout(QtWidgets.QLayout):
    """Left-to-right wrapping layout for the tiles.

    Qt ships no flow layout. A `QGridLayout` would need a column count chosen in
    advance, and the right number depends on the width the page actually gets --
    which is how spaCR's bands ended up either stopping short of the right edge
    or drawing two tiles half a metre wide. This measures instead.
    """

    def __init__(self, parent=None, spacing: int = 12):
        super().__init__(parent)
        self._items: list = []
        self._space = spacing
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item):                       # noqa: N802
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, i):                           # noqa: N802
        return self._items[i] if 0 <= i < len(self._items) else None

    def takeAt(self, i):                           # noqa: N802
        return self._items.pop(i) if 0 <= i < len(self._items) else None

    def expandingDirections(self):                 # noqa: N802
        return QtCore.Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:           # noqa: N802
        return True

    def heightForWidth(self, width: int) -> int:   # noqa: N802
        return self._lay(QtCore.QRect(0, 0, width, 0), apply=False)

    def setGeometry(self, rect):                   # noqa: N802
        super().setGeometry(rect)
        self._lay(rect, apply=True)

    def sizeHint(self) -> QtCore.QSize:            # noqa: N802
        return self.minimumSize()

    def minimumSize(self) -> QtCore.QSize:         # noqa: N802
        size = QtCore.QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        return size

    def _lay(self, rect: QtCore.QRect, apply: bool) -> int:
        """Place the tiles, or measure what placing them would need.

        Tiles stretch up to `TILE_MAX_W` to reach the right edge, so a band with
        two modules in it does not leave the row looking unfinished -- but they
        are capped, because without a cap those two tiles grow absurd.
        """
        x, y, line = rect.x(), rect.y(), 0
        right = rect.right() + 1
        for item in self._items:
            hint = item.sizeHint()
            w, h = hint.width(), hint.height()
            if line and x + w > right:
                x, y, line = rect.x(), y + line + self._space, 0
            if apply:
                item.setGeometry(QtCore.QRect(QtCore.QPoint(x, y), QtCore.QSize(w, h)))
            x += w + self._space
            line = max(line, h)
        return y + line - rect.y()


class HomeScreen(QtWidgets.QWidget):
    """The landing view. `module_chosen(str key)` fires when a tile is pressed.

    :param counts: optional `{name: (value, tip)}` drawn in the Measured panel, so
        the widget does not have to know how to count anything itself.
    :param reasons: optional `{module key: reason}` overriding the computed
        availability, which is how the window reports what it actually built.
    """

    module_chosen = QtCore.Signal(str)

    def __init__(self, counts: dict | None = None, reasons: dict | None = None,
                 parent=None):
        super().__init__(parent)
        self.setObjectName("HomeScreen")
        self._reasons = dict(reasons or {})
        self._tiles: dict = {}
        palette = TH.palette_for(TH.current_theme() if hasattr(TH, "current_theme") else "dark")
        self.setStyleSheet(tile_qss(palette))

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(*(TH.SPACING["xl"],) * 2, TH.SPACING["xl"], TH.SPACING["lg"])
        outer.setSpacing(TH.SPACING["lg"])
        outer.addLayout(self._masthead())

        split = QtWidgets.QHBoxLayout()
        split.setSpacing(TH.SPACING["xl"])
        split.addWidget(self._bands(), 1)
        split.addWidget(self._aside(counts or {}), 0)
        outer.addLayout(split, 1)

    # ------------------------------------------------------------------ pieces
    def _masthead(self) -> QtWidgets.QVBoxLayout:
        """Name over one sentence. The sentence says what the program does, because
        a masthead that only repeats the name wastes the one place everybody looks."""
        box = QtWidgets.QVBoxLayout()
        box.setSpacing(2)
        title = QtWidgets.QLabel("ParaScape")
        title.setObjectName("HeroTitle")
        box.addWidget(title)
        sub = QtWidgets.QLabel(
            "Everything published about a Toxoplasma or Plasmodium gene, in one place "
            "— and what has not been measured, predicted with an uncertainty on it.")
        sub.setObjectName("HeroSub")
        sub.setWordWrap(True)
        box.addWidget(sub)
        return box

    def _bands(self) -> QtWidgets.QWidget:
        """The tiles, in bands, each band under its heading and its question."""
        page = QtWidgets.QWidget()
        col = QtWidgets.QVBoxLayout(page)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(TH.SPACING["lg"])
        for title, note, members in MOD.bands():
            head = QtWidgets.QLabel(title.upper())
            head.setObjectName("PanelTitle")
            col.addWidget(head)
            hint = QtWidgets.QLabel(note)
            hint.setObjectName("TileNote")
            col.addWidget(hint)
            row = QtWidgets.QWidget()
            flow = _FlowLayout(row, spacing=TH.SPACING["md"])
            for module in members:
                tile = self._tile(module)
                flow.addWidget(tile)
                self._tiles[module.key] = tile
            col.addWidget(row)
        col.addStretch(1)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        # Both the area AND its viewport, by object name. A QScrollArea paints its
        # viewport from the style's Base role, which is near-white on every theme
        # here; without this the band column came out white and the headings --
        # light text, correctly -- disappeared into it.
        scroll.setObjectName("HomeScroll")
        scroll.viewport().setObjectName("HomeScrollViewport")
        page.setObjectName("HomeBands")
        scroll.setWidget(page)
        return scroll

    def _tile(self, module) -> AppTile:
        """One tile, with its availability resolved. A reason supplied by the window
        wins over the computed one: the window knows which docks it actually built."""
        reason = self._reasons.get(module.key)
        if reason is None:
            reason = MOD.unavailable_reason(module)
        tile = AppTile(module.key, module.name, module.note, reason=reason,
                       width=TILE_MAX_W)
        if tile.available:
            tile.clicked.connect(lambda _=False, k=module.key: self.module_chosen.emit(k))
        return tile

    def _aside(self, counts: dict) -> QtWidgets.QWidget:
        """The right-hand column: what is measured, what is inferred, what is running.

        Three panels rather than spaCR's four, and none of them is news. The second
        one is the one that matters: it says how much of the map is predicted rather
        than observed, which is the number a person has to see before trusting
        anything else on screen.
        """
        side = QtWidgets.QWidget()
        side.setFixedWidth(ASIDE_W)
        col = QtWidgets.QVBoxLayout(side)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(TH.SPACING["md"])

        measured = HomePanel("Measured")
        if counts:
            for label, (value, tip) in counts.items():
                measured.add_row(label, value, tip)
        else:
            measured.add_note("Counts appear once the data release is installed.")
        col.addWidget(measured)

        inferred = HomePanel("Inferred")
        inferred.add_note(
            "Nothing yet. A prediction is written into the map only when it has been "
            "graded on genes it was not shown and was right at least 90% of the time. "
            "No question has cleared that bar, so every value you see is measured.")
        col.addWidget(inferred)

        running = HomePanel("Running")
        running.add_note("Nothing running. Jobs started here appear in this panel.")
        self._running_panel = running
        col.addWidget(running)
        col.addStretch(1)
        return side

    # ------------------------------------------------------------------- query
    def tile(self, key: str) -> AppTile:
        """One tile by module key, for a test or the navigator to act on."""
        return self._tiles[key]

    def unavailable(self) -> dict:
        """Every module that cannot be opened, and why.

        The honesty check: this is what the home screen is admitting to, and a test
        asserts it is never empty of reasons where a reason is due.
        """
        return {k: t.reason for k, t in self._tiles.items() if not t.available}

    def open_module(self, key: str) -> bool:
        """Act as though the tile were pressed. Returns whether it was openable."""
        tile = self._tiles.get(key)
        if tile is None or not tile.available:
            return False
        self.module_chosen.emit(key)
        return True
