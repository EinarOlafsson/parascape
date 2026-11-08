"""The question scoreboard in the application (instructions 70.06 and 73): how well each question is answered.

One row per yes/no question scored through the harness (`scripts/score_question_catalogue.py`): its category, the
label it reads, how many genes are labelled, the frozen-test AUC against chance (0.5), average precision against the
share of yes, and whether it may write cells into the map (frozen-test lower bound >= 0.90). Read from the shipped
`data/question_scoreboard.tsv`; nothing is computed here.
"""
from __future__ import annotations

import os

import pandas as pd

COLUMNS = (("question", "question"), ("category", "kind of knowledge"), ("labelled", "genes labelled"),
           ("frozen_auc", "held-out AUC (0.5 = chance)"), ("frozen_ap", "average precision"),
           ("prevalence", "share of yes"), ("may_fill", "may fill the map"), ("top_families", "evidence that carries it"))


def load(path: str | None = None) -> pd.DataFrame:
    """The scoreboard table, or an empty frame when it is not shipped."""
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "question_scoreboard.tsv")
    if not os.path.exists(path):
        return pd.DataFrame(columns=[c for c, _ in COLUMNS])
    return pd.read_csv(path, sep="\t")


def summary(frame: pd.DataFrame) -> str:
    """One sentence for the top of the window."""
    if not len(frame):
        return "No question scoreboard is installed."
    return (f"{len(frame)} questions scored on held-out genes: median AUC {frame['frozen_auc'].median():.2f} "
            f"(0.5 is chance); {int((frame['frozen_auc'] >= 0.7).sum())} at 0.70 or better; "
            f"{int(frame['may_fill'].sum())} trusted enough to fill the map.")


def dialog(parent=None, frame: pd.DataFrame | None = None):
    """A window with the summary sentence and a sortable table of every question."""
    from PySide6 import QtCore, QtWidgets
    frame = load() if frame is None else frame
    d = QtWidgets.QDialog(parent)
    d.setWindowTitle("Question scoreboard")
    d.resize(1100, 640)
    lay = QtWidgets.QVBoxLayout(d)
    note = QtWidgets.QLabel(summary(frame))
    note.setWordWrap(True)
    lay.addWidget(note)
    table = QtWidgets.QTableWidget(len(frame), len(COLUMNS))
    table.setObjectName("QuestionScoreboardTable")
    table.setToolTip("Every question scored through the harness, one row each; click a heading to sort.")
    table.setHorizontalHeaderLabels([label for _c, label in COLUMNS])
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    for i, row in enumerate(frame.itertuples(index=False)):
        for j, (col, _label) in enumerate(COLUMNS):
            value = getattr(row, col, "")
            item = QtWidgets.QTableWidgetItem()
            if isinstance(value, float) and col != "prevalence":
                item.setData(QtCore.Qt.ItemDataRole.DisplayRole, round(value, 3))
            elif col == "prevalence":
                item.setData(QtCore.Qt.ItemDataRole.DisplayRole, f"{float(value):.0%}")
            else:
                item.setData(QtCore.Qt.ItemDataRole.DisplayRole, str(value))
            table.setItem(i, j, item)
    table.setSortingEnabled(True)
    table.resizeColumnsToContents()
    lay.addWidget(table)
    d.table = table
    return d
