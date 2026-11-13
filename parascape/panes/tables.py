"""A table that shows the names of its own columns, before anything is in it.

WHY THIS EXISTS (instruction 98, WP11, found by rendering the modules and looking).

A `QTableWidget` gives every column the same default width and only sizes to the content when
it is asked to -- which the panels here do after loading something. So the table a person
meets FIRST, with nothing in it, draws each heading centred in a box too narrow for it and
CUTS it at both ends. The Drop module said "NTIFIER IN THE F" where it meant "identifier in
the file", which says neither which column it is nor that anything has been cut, and the
Scorecard said "NFERRED CLAIMS".

The fix belongs in one place because the mistake is the same every time: size the columns to
whatever is widest, headings included, and let the last one take the slack so the row of
headings reaches the right edge. `tests/test_module_contract.py` holds it for every table
inside a module shell.
"""
from __future__ import annotations

from PySide6 import QtCore, QtWidgets

#: Room kept beside a heading so that the first and last letters are not drawn on the edge of
#: the section. Qt sizes a section to the exact width of its text, and the outermost pixel of
#: an antialiased letter then falls outside it.
HEADING_MARGIN = 8


def _widen(table: QtWidgets.QTableWidget) -> None:
    """Widen any column narrower than its own heading. Never narrows one."""
    header = table.horizontalHeader()
    metrics = header.fontMetrics()
    for column in range(table.columnCount()):
        item = table.horizontalHeaderItem(column)
        if item is None:
            continue
        wanted = metrics.horizontalAdvance(item.text()) + HEADING_MARGIN
        if header.sectionSize(column) < wanted:
            header.resizeSection(column, wanted)


class _KeepHeadings(QtCore.QObject):
    """Re-measure when the font the headings are drawn in can have changed.

    Measuring once where the table is built is not enough, and the reason is worth recording:
    a table is built before the module shell has put its stylesheet on it, so the font it is
    measured in is the default one and the font it is DRAWN in is larger. The Information
    module was out by 37 px that way, which is most of a word.
    """

    WATCHED = (QtCore.QEvent.Type.Show, QtCore.QEvent.Type.FontChange,
               QtCore.QEvent.Type.StyleChange, QtCore.QEvent.Type.ApplicationFontChange)

    def eventFilter(self, obj, event):             # noqa: N802 - Qt's name
        if event.type() in self.WATCHED and isinstance(obj, QtWidgets.QTableWidget):
            _widen(obj)
        return False


def show_headings(table: QtWidgets.QTableWidget, stretch_last: bool = True) -> None:
    """Give every column of `table` the room its heading needs, now and from now on.

    Called once, where the table is built, before it has any rows. Panels that resize to the
    content again after filling lose nothing: this only ever widens a column, and a heading
    narrower than its data was never the problem.

    :param table: the table.
    :param stretch_last: whether the last column takes the leftover width, so that the row of
        headings reaches the right edge rather than stopping in the middle of the table.
    """
    table.resizeColumnsToContents()
    _widen(table)
    table.horizontalHeader().setStretchLastSection(bool(stretch_last))
    if not table.findChildren(_KeepHeadings):
        table.installEventFilter(_KeepHeadings(table))
