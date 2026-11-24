"""Back and forward through ParaScape, as in a browser.

Einar, 2026-11-24: "when navigating through parascape the user should always be able to go
back and forward" and "the dock should always be there for the user to navigate where they
want".

A *place* is where the person is: a page of the landing (``("page", "Navigator")``), a
module (``("module", "network")``), or the main window with its panels
(``("window", "")``). Every route that moves the person -- a tile, the tree, the Navigator,
a next step, Home, a panel brought forward -- passes through `Landing.show_page_index`,
`Landing.show_module`, `enter_landing` or `leave_landing`, and those record the place here.
Going back or forward replays a place through the same routes, with recording paused, so a
replay never adds a place of its own.

`History` is pure and testable; `install` gives a window the history, the Back and Forward
actions (Alt+Left and Alt+Right, and the mouse side buttons) and the buttons at the top of
the tree.
"""
from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets


class History:
    """Places visited, and where in them the person is now.

    Recording a place after going back drops the places that were ahead, as a browser does.
    The same place twice in a row is one visit.
    """

    #: The most places kept.
    LIMIT = 200

    def __init__(self):
        self.places: list = []
        self.index = -1
        self.replaying = False

    @property
    def current(self):
        """The place the person is at, or None before any."""
        return self.places[self.index] if 0 <= self.index < len(self.places) else None

    def record(self, place) -> bool:
        """Note that the person is now at `place`. Returns whether it was a new visit."""
        if self.replaying or place is None or place == self.current:
            return False
        del self.places[self.index + 1:]
        self.places.append(place)
        if len(self.places) > self.LIMIT:
            del self.places[0]
        self.index = len(self.places) - 1
        return True

    def can_go(self, step: int) -> bool:
        """Whether there is a place `step` away (-1 back, +1 forward)."""
        return 0 <= self.index + step < len(self.places)

    def step(self, step: int):
        """Move `step` places and return the place now current, or None if there is none."""
        if not self.can_go(step):
            return None
        self.index += step
        return self.places[self.index]


def history(window) -> History:
    """The window's history, made on first use."""
    found = getattr(window, "nav_history", None)
    if found is None:
        found = History()
        window.nav_history = found
    return found


def landing_place(landing):
    """The place the landing is showing: its page, or the module in its stack."""
    if landing is None:
        return None
    page = landing.page
    if page != "Module":
        return ("page", page)
    shell = landing.stack.currentWidget()
    key = getattr(shell, "key", "") or (shell.objectName() if shell is not None else "")
    return ("module", str(key)) if key else None


def record(window, place) -> None:
    """Record a place for `window` and refresh the Back and Forward controls."""
    if window is None or place is None:
        return
    history(window).record(place)
    refresh(window)


def go(window, step: int) -> str:
    """Go `step` places back (-1) or forward (+1). Returns what the status bar was told."""
    from . import landing_install as LI
    h = history(window)
    place = h.step(step)
    if place is None:
        return ""
    kind, name = place
    h.replaying = True
    try:
        if kind == "window":
            LI.leave_landing(window)
            message = "the main window"
        elif kind == "page":
            window.landing.show_page(name)
            LI.enter_landing(window)
            message = name
        else:
            message = LI.open_module(window, name)
    finally:
        h.replaying = False
    refresh(window)
    try:
        window.statusBar().showMessage(f"{'Back' if step < 0 else 'Forward'}: {message}", 4000)
    except Exception:                       # a window without a bar must still navigate
        pass
    return message


def refresh(window) -> None:
    """Enable Back and Forward as the history allows."""
    h = history(window)
    for attr, step in (("nav_back_action", -1), ("nav_forward_action", 1)):
        action = getattr(window, attr, None)
        if action is not None:
            action.setEnabled(h.can_go(step))


class _SideButtons(QtCore.QObject):
    """The mouse side buttons go back and forward anywhere in the window."""

    def __init__(self, window):
        super().__init__(window)
        self._window = window

    def eventFilter(self, obj, event):                    # noqa: N802
        """Back/forward on the mouse side buttons; every other event passes through."""
        if event.type() == QtCore.QEvent.Type.MouseButtonPress and isinstance(obj, QtWidgets.QWidget) \
                and obj.window() is self._window:
            button = event.button()
            if button == QtCore.Qt.MouseButton.BackButton:
                go(self._window, -1)
                return True
            if button == QtCore.Qt.MouseButton.ForwardButton:
                go(self._window, 1)
                return True
        return False


def nav_bar(window) -> QtWidgets.QWidget:
    """A row with Back and Forward, for the top of the tree."""
    bar = QtWidgets.QWidget()
    bar.setObjectName("NavBar")
    row = QtWidgets.QHBoxLayout(bar)
    row.setContentsMargins(6, 4, 6, 2)
    row.setSpacing(4)
    for attr in ("nav_back_action", "nav_forward_action"):
        button = QtWidgets.QToolButton()
        button.setDefaultAction(getattr(window, attr))
        button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextOnly)
        button.setAutoRaise(True)
        button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        row.addWidget(button)
    row.addStretch(1)
    return bar


def install(window) -> None:
    """Give `window` the history, Back and Forward (Alt+Left, Alt+Right and the mouse side
    buttons), and record where it is now."""
    history(window)
    for attr, text, key, _icon, tip, step in (
            ("nav_back_action", "◀  Back", "Alt+Left", None,
             "Go back to where you were before (Alt+Left, or the mouse back button).", -1),
            ("nav_forward_action", "Forward  ▶", "Alt+Right", None,
             "Go forward again after going back (Alt+Right, or the mouse forward button).", 1)):
        action = QtGui.QAction(text, window)
        action.setShortcut(QtGui.QKeySequence(key))
        action.setToolTip(tip)
        action.setStatusTip(tip)
        action.triggered.connect(lambda _=False, s=step: go(window, s))
        window.addAction(action)
        setattr(window, attr, action)
    app = QtWidgets.QApplication.instance()
    if app is not None:
        window._nav_side_buttons = _SideButtons(window)
        app.installEventFilter(window._nav_side_buttons)
    refresh(window)
