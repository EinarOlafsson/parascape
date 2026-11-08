"""The question scoreboard (instructions 70.06, 73): shipped, summarised, shown."""
import pandas as pd

from parascape import question_scoreboard as QS


def test_the_shipped_scoreboard_has_every_column_the_window_shows():
    frame = QS.load()
    assert len(frame) > 100
    assert {c for c, _ in QS.COLUMNS} <= set(frame.columns)
    assert frame["frozen_auc"].dropna().between(0, 1).all()        # a question without a frozen test has none


def test_the_summary_says_chance_and_counts_fills():
    frame = pd.DataFrame({"frozen_auc": [0.5, 0.8, 0.95], "may_fill": [False, False, True]})
    text = QS.summary(frame)
    assert "3 questions" in text and "0.5 is chance" in text and "1 trusted" in text
    assert QS.summary(pd.DataFrame()) == "No question scoreboard is installed."


def test_the_window_lists_every_question(qtbot):
    d = QS.dialog()
    qtbot.addWidget(d)
    assert d.table.rowCount() == len(QS.load())
    assert d.table.isSortingEnabled()
