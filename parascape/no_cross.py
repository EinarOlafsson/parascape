"""Crosses: a small x CLOSES things; a search field never gets a clear-button cross.

Einar, 2026-11-12: "I don't like the little x cross (I don't want to see that anywhere)".
REVISED by Einar, 2026-11-24: "the evidence tab can have a cross to the top right not a close
button, use the x throughout and always opt for more discreet GUI elements for things like
this." So a panel, tab or window that can be closed is closed with the small cross Qt draws (a
closable dock or tab, a title-bar cross), never with a worded "Close" button.

What stays from the first rule is the clear button inside a search field: that cross empties a
field rather than closing anything, and Escape does it instead, as in spaCR. This module holds:

* :func:`clear_on_escape`, what a search field does instead of a clear button: Escape empties it.
* :func:`find_crosses`, which walks a window and reports every search field that would draw a
  clear-button cross. The test of the whole application calls it on the real window.

It looks at what Qt would draw, not at what the code says. A grep for ``setClearButtonEnabled``
misses a line edit that got its button from somewhere else.
"""
from __future__ import annotations

from PySide6 import QtCore, QtWidgets

#: Button object names and glyphs that mean "close this" when found as the whole text of a button.
CROSS_TEXT = frozenset({"x", "X", "×", "✕", "✖", "✗", "╳", "⨯", "❌", "✘"})

#: The Qt feature that draws the cross in the title of a dock.
_CLOSABLE = QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetClosable


class _EscapeClears(QtCore.QObject):
    """Empties a line edit when Escape is pressed in it and it has text."""

    def eventFilter(self, obj, event):                    # noqa: N802
        """Empty the field on Escape, and say the event was used only when there was text."""
        if (event.type() == QtCore.QEvent.Type.KeyPress
                and event.key() == QtCore.Qt.Key.Key_Escape and obj.text()):
            obj.clear()
            obj.textEdited.emit("")
            return True
        return False


def clear_on_escape(edit: QtWidgets.QLineEdit) -> None:
    """Give a search field the way to empty it that is not a small cross.

    Turns the built-in clear button off, and makes Escape empty the field. Escape with nothing in
    the field is left alone, so it still reaches whatever else listens for it.

    :param edit: the field.
    """
    edit.setClearButtonEnabled(False)
    if edit.property("parascapeEscapeClears"):
        return
    edit.setProperty("parascapeEscapeClears", True)
    helper = _EscapeClears(edit)
    edit.installEventFilter(helper)


def find_crosses(root: QtWidgets.QWidget) -> list:
    """Every search field under ``root`` that would draw a clear-button cross, as sentences.

    Closing crosses (closable docks and tabs) are allowed since 2026-11-24 and are not reported.

    :param root: usually the main window.
    :returns: one sentence per offender, naming the widget; empty when there are none.
    """
    found = []

    def name(widget) -> str:
        """The class and object name of a widget, so the failure says which one it is."""
        return f"{type(widget).__name__}#{widget.objectName() or '?'}"

    for edit in root.findChildren(QtWidgets.QLineEdit):
        if edit.isClearButtonEnabled():
            found.append(f"{name(edit)} has a clear button")
    return found


def closable_docks(root: QtWidgets.QWidget) -> list:
    """The titles of the docks under ``root`` that close with the small cross."""
    return [dock.windowTitle() for dock in root.findChildren(QtWidgets.QDockWidget)
            if dock.features() & _CLOSABLE]
