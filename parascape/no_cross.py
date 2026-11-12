"""No small cross, anywhere. Einar: "I don't want to see that anywhere."

Qt draws a small cross for you in four places that are easy to switch on without meaning to: a
closable tab, a closable dock, a clear button inside a line edit, and the close mark of a window
without a title bar. This module holds the two things that keep them out:

* :func:`clear_on_escape`, what a search field does instead of a clear button, which is also what
  spaCR does: Escape empties the field.
* :func:`find_crosses`, which walks a window and reports every widget that would draw one. The
  test of the whole application calls it on the real window, so a cross that creeps back fails the
  suite instead of being found by a person.

It deliberately looks at what Qt would draw, not at what the code says. A grep for
``setClearButtonEnabled`` misses a line edit that got its button from somewhere else.
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
    """Every widget under ``root`` that would draw a small cross, as plain sentences.

    Looks for closable tabs, closable docks, line edits with a clear button (the cross inside a
    search field), the close buttons Qt creates inside tab bars and dock titles, and any button or
    label whose whole text is a cross glyph.

    :param root: usually the main window.
    :returns: one sentence per offender, naming the widget; empty when there are none.
    """
    found = []

    def name(widget) -> str:
        """The class and object name of a widget, so the failure says which one it is."""
        return f"{type(widget).__name__}#{widget.objectName() or '?'}"

    for bar in root.findChildren(QtWidgets.QTabBar):
        if bar.tabsClosable():
            found.append(f"{name(bar)} has closable tabs")
    for tabs in root.findChildren(QtWidgets.QTabWidget):
        if tabs.tabsClosable():
            found.append(f"{name(tabs)} has closable tabs")
    for dock in root.findChildren(QtWidgets.QDockWidget):
        if dock.features() & _CLOSABLE:
            found.append(f"dock {dock.windowTitle()!r} can be closed with a cross")
    for edit in root.findChildren(QtWidgets.QLineEdit):
        if edit.isClearButtonEnabled():
            found.append(f"{name(edit)} has a clear button")
    for button in root.findChildren(QtWidgets.QAbstractButton):
        if button.text().strip() in CROSS_TEXT:
            found.append(f"{name(button)} is a cross ({button.text()!r})")
        # Qt keeps a hidden close button inside every dock title even when the dock is not
        # closable, so only one that is not hidden counts.
        if ("close" in button.objectName().lower() and not button.text().strip()
                and not button.isHidden()):
            found.append(f"{name(button)} looks like an icon-only close button")
    for label in root.findChildren(QtWidgets.QLabel):
        if label.text().strip() in CROSS_TEXT:
            found.append(f"{name(label)} is a cross ({label.text()!r})")
    return found
