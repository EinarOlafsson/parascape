"""The home screen's tile and panel widgets, ported from spaCR.

PORTED AND MODIFIED FROM spaCR
==============================
Derived from ``spacr/qt/widgets/home.py`` and ``spacr/qt/theme.py``
(https://github.com/EinarOlafsson/spacr).

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

**Modifications made in porting**: rewritten against ParaScape's own palette
(``parascape.theme.palette_for``) and geometry tokens rather than spaCR's; the
maturity "stage" vocabulary replaced by ParaScape's availability model, where a
module whose data or view is absent is drawn with a dashed rim and carries the
reason; ``ElidingLabel`` replaced with Qt's own elision since ParaScape has no
equivalent widget; and the description restored to the tile, because ParaScape's
modules are not guessable from a one-word name the way an image-analysis
pipeline's stages are.

The rest of the file follows spaCR's decisions deliberately, including the ones
that look odd:

* **The height floor lives in the QSS, not in ``setFixedSize``.** A blanket
  ``QPushButton { min-height: 22px }`` in an app stylesheet becomes a real
  ``setMinimumHeight(22)``, and Qt's ``qSmartMinSize`` lets an explicit minimum
  beat a size hint. Without the QSS floor, a page that does not quite fit
  collapses every tile to 22 px and paints the name over the icon.
* **No resting rim.** A rim at rest is decoration; a rim means hover, focus or
  unavailability, so it carries information when it appears.
* **``heightForWidth`` as well as ``sizeHint``**, because
  ``QWidgetItem::sizeHint`` consults the former whenever it is available, and a
  tile whose height does not follow its width lets the grid stretch one tile and
  not its neighbours.
"""
from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from . import theme as TH

#: Tile geometry. spaCR's numbers, because the proportion is what makes a grid of
#: them read as a launcher rather than a list: 172x124 at a 48 px icon.
TILE_W = 172
TILE_MAX_W = 260
TILE_H = 124
TILE_ICON_PX = 48

#: The right-hand column. Fixed, because it holds numbers and a column of numbers
#: that changes width as they change is unreadable.
ASIDE_W = 268


def tile_qss(palette: dict) -> str:
    """The stylesheet for the tiles and panels of the home screen.

    Separate from `theme.stylesheet` so the home screen can be restyled without a
    full repolish of thirteen docks, and so the height floor sits beside the
    widget that depends on it.

    :param palette: a `theme.palette_for` mapping.
    :returns: QSS for `#AppTile`, its states, and `#HomePanel`.
    """
    P, R = palette, TH.RADIUS
    accent = P["accent"]
    return f"""
/* The tile. min-height/min-width are the reason this rule exists separately:
   see the module docstring on the 22 px blanket minimum. Padding is zero
   because the tile's own layout owns its margins. */
QPushButton#AppTile {{
    background-color: {P["surface"]};
    color: {P["fg"]};
    border: none;
    border-radius: {R["tile"]}px;
    padding: 0px;
    min-height: {TILE_H}px;
    min-width: {TILE_W}px;
    text-align: center;
    font-family: {TH.FONT_FAMILY};
}}
QPushButton#AppTile:hover {{
    background-color: {P["surface_hi"]};
    border: 1px solid {accent};
}}
QPushButton#AppTile:pressed {{
    background-color: {P["accent_soft"]};
    border: 1px solid {accent};
}}
/* Unavailable: a dashed rim and muted text, never hidden. A module that has
   disappeared teaches the user it never existed; one that is visibly greyed
   with a reason teaches them what is missing. */
QPushButton#AppTile[available="false"] {{
    background-color: {P["surface_alt"]};
    color: {P["fg_muted"]};
    border: 1px dashed {P["border"]};
}}
QPushButton#AppTile[available="false"]:hover {{
    background-color: {P["surface_alt"]};
    border: 1px dashed {P["border_soft"]};
}}
/* After the state rules, so focus stays visible while a tile is hovered. */
QPushButton#AppTile:focus {{
    border: 1px solid {accent};
}}
QLabel#TileName {{
    color: {P["fg"]};
    font-family: {TH.FONT_FAMILY};
    font-size: 14px;
    font-weight: 500;
    background: transparent;
}}
QLabel#TileNote {{
    color: {P["fg_muted"]};
    font-family: {TH.FONT_FAMILY};
    font-size: 11px;
    font-weight: 300;
    background: transparent;
}}
/* The page the tiles sit on, and the scroll area holding it. A QScrollArea
   paints its viewport from the Base role, which is near-white; unstyled, the
   band column renders white and every heading on it -- light text -- becomes
   invisible. Caught by rendering the page to an image rather than by reading it. */
QWidget#HomeScreen, QScrollArea#HomeScroll, QWidget#HomeScrollViewport,
QWidget#HomeBands {{
    background-color: {P["bg"]};
    border: none;
}}
/* With the animated backdrop on, the same containers go transparent: a backdrop
   is only as visible as what sits above it, and these cover the whole page. The
   opaque rule above is what applies when animation is off, so nothing reverts to
   the near-white Base colour. */
QWidget#HomeScreen[ambient="on"], QScrollArea#HomeScroll[ambient="on"],
QWidget#HomeScrollViewport[ambient="on"], QWidget#HomeBands[ambient="on"] {{
    background: transparent;
}}
/* Pure layout containers: one row of tiles, the right-hand column. They draw
   nothing of their own and must never cover the page. Standalone they were
   invisible; inside the real window the application stylesheet gives a plain
   QWidget a background, and each row and the aside became an opaque slab that hid
   the animation behind them. */
QWidget#HomeBandRow, QWidget#HomeAside {{
    background: transparent;
}}
/* The scrollbar, which otherwise renders as a bright slab against the page. */
QScrollArea#HomeScroll QScrollBar:vertical {{
    background: transparent;
    width: 10px;
    margin: 0px;
}}
QScrollArea#HomeScroll QScrollBar::handle:vertical {{
    background: {P["border"]};
    border-radius: 5px;
    min-height: 32px;
}}
QScrollArea#HomeScroll QScrollBar::handle:vertical:hover {{
    background: {P["fg_muted"]};
}}
QScrollArea#HomeScroll QScrollBar::add-line:vertical,
QScrollArea#HomeScroll QScrollBar::sub-line:vertical {{
    height: 0px;
}}
QScrollArea#HomeScroll QScrollBar::add-page:vertical,
QScrollArea#HomeScroll QScrollBar::sub-page:vertical {{
    background: transparent;
}}
QFrame#HomePanel {{
    background-color: {P["surface"]};
    border: 1px solid {P["border_soft"]};
    border-radius: {R["card"]}px;
}}
QLabel#PanelTitle {{
    color: {P["fg"]};
    font-family: {TH.FONT_FAMILY};
    font-size: 12px;
    font-weight: 600;
    background: transparent;
}}
QLabel#HeroTitle {{
    color: {P["fg"]};
    font-family: {TH.FONT_FAMILY};
    font-size: 26px;
    font-weight: 300;
    background: transparent;
}}
QLabel#HeroSub {{
    color: {P["fg_muted"]};
    font-family: {TH.FONT_FAMILY};
    font-size: 13px;
    font-weight: 300;
    background: transparent;
}}
"""


class _TileIcon(QtWidgets.QLabel):
    """The tile's mark. Its own label so the icon can be centred independently of
    the name, and so it stays transparent against the tile's own background."""

    def __init__(self, icon: QtGui.QIcon, px: int, parent=None):
        super().__init__(parent)
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.setFixedHeight(px)
        self.setStyleSheet("background: transparent;")
        if icon is not None and not icon.isNull():
            self.setPixmap(icon.pixmap(QtCore.QSize(px, px)))


class AppTile(QtWidgets.QPushButton):
    """One module, as a tile: icon over name over a one-line note.

    :param key: the module's identifier, which is what `clicked` is reported as.
    :param name: the module's name, drawn on the tile.
    :param note: one line under the name saying what the module is for.
    :param icon: the module's mark, or None.
    :param reason: why the module cannot be opened, or empty when it can. A
        non-empty reason makes the tile unavailable: dashed rim, muted text,
        disabled, and the reason on its tooltip.
    """

    def __init__(self, key: str, name: str, note: str = "", icon=None,
                 reason: str = "", *, width: int = TILE_W, height: int = TILE_H,
                 icon_px: int = TILE_ICON_PX, parent=None):
        super().__init__(parent)
        self.key = key
        self._name = name
        self._reason = str(reason or "")
        self._size = QtCore.QSize(int(width), int(height))
        self.setObjectName("AppTile")
        # Set BEFORE first polish, or the [available="false"] rule does not apply
        # until something else forces a repolish.
        self.setProperty("available", "false" if self._reason else "true")
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor if not self._reason
                       else QtCore.Qt.CursorShape.ForbiddenCursor)
        self.setAccessibleName(name)
        self.setAccessibleDescription(self._reason or note)
        self.setToolTip(f"{name} — {self._reason}" if self._reason else f"{name} — {note}")
        self.setEnabled(not self._reason)

        col = QtWidgets.QVBoxLayout(self)
        col.setContentsMargins(10, 10, 10, 10)
        col.setSpacing(6)
        col.addStretch(1)
        if icon is not None:
            col.addWidget(_TileIcon(icon, icon_px), 0,
                          QtCore.Qt.AlignmentFlag.AlignHCenter)
        label = QtWidgets.QLabel(name)
        label.setObjectName("TileName")
        label.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
        label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                            QtWidgets.QSizePolicy.Policy.Fixed)
        col.addWidget(label)
        self._name_label = label

        shown = self._reason or note
        if shown:
            sub = QtWidgets.QLabel(shown)
            sub.setObjectName("TileNote")
            sub.setAlignment(QtCore.Qt.AlignmentFlag.AlignHCenter)
            sub.setWordWrap(True)
            sub.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                              QtWidgets.QSizePolicy.Policy.Preferred)
            col.addWidget(sub)
            self._note_label = sub
        else:
            self._note_label = None
        col.addStretch(1)

    @property
    def available(self) -> bool:
        """Whether this module can be opened. False when it carries a reason."""
        return not self._reason

    @property
    def reason(self) -> str:
        """Why it cannot be opened, or empty. Shown on the tile and its tooltip."""
        return self._reason

    @property
    def name_label(self) -> QtWidgets.QLabel:
        """The label carrying the module name, so a text-fits sweep can measure it."""
        return self._name_label

    def heightForWidth(self, width: int) -> int:        # noqa: N802
        """At least the tile height, more if a child needs it.

        Overridden as well as `sizeHint` because `QWidgetItem::sizeHint` reads
        this in preference whenever it is available.
        """
        return max(self._size.height(), super().heightForWidth(width))

    def sizeHint(self) -> QtCore.QSize:                 # noqa: N802
        """Preferred size, height derived from width so the grid cannot stretch
        one tile and not its neighbours."""
        return QtCore.QSize(self._size.width(), self.heightForWidth(self._size.width()))

    def minimumSizeHint(self) -> QtCore.QSize:          # noqa: N802
        """The same as `sizeHint`: a tile does not shrink below its shape."""
        return self.sizeHint()


class HomePanel(QtWidgets.QFrame):
    """A titled box in the right-hand column: a heading over rows of label/value.

    :param title: the heading, drawn in small caps-weight above the rows.
    """

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setObjectName("HomePanel")
        self._rows = QtWidgets.QVBoxLayout(self)
        self._rows.setContentsMargins(12, 10, 12, 10)
        self._rows.setSpacing(4)
        head = QtWidgets.QLabel(title.upper())
        head.setObjectName("PanelTitle")
        self._rows.addWidget(head)
        self._title = head

    def add_row(self, label: str, value: str, tip: str = "") -> None:
        """One label/value line. Values are right-aligned so a column of numbers
        reads down rather than ragged."""
        row = QtWidgets.QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        left = QtWidgets.QLabel(label)
        left.setObjectName("TileNote")
        right = QtWidgets.QLabel(value)
        right.setObjectName("TileName")
        right.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                           | QtCore.Qt.AlignmentFlag.AlignVCenter)
        if tip:
            left.setToolTip(tip)
            right.setToolTip(tip)
        row.addWidget(left)
        row.addStretch(1)
        row.addWidget(right)
        self._rows.addLayout(row)

    def add_note(self, text: str) -> None:
        """A wrapped sentence rather than a row, for a panel that must explain
        itself instead of only counting."""
        note = QtWidgets.QLabel(text)
        note.setObjectName("TileNote")
        note.setWordWrap(True)
        self._rows.addWidget(note)
