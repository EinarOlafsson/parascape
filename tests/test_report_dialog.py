"""Error reports from the window (instruction 97): asked once, shown in full, sent on a click.
The relay is replaced by a fake; nothing leaves the machine."""
import sys

import pytest
from PySide6 import QtCore, QtWidgets

from parascape import relay as R
from parascape import report_dialog as RD


@pytest.fixture
def window(qtbot, tmp_path, monkeypatch):
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "relay.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))
    sent = []
    monkeypatch.setattr(R, "send_report", lambda report: (sent.append(report) or
                        {"issue_url": "https://github.com/EinarOlafsson/parascape/issues/9"}))
    win = QtWidgets.QMainWindow()
    win.menuBar().addMenu("&Help")
    qtbot.addWidget(win)
    reporter = RD.install(win)
    yield win, reporter, sent
    reporter.uninstall()


def _report():
    try:
        raise ValueError("broken")
    except ValueError:
        return R.report_from_exception(*sys.exc_info())


def _dialogs():
    return [w for w in QtWidgets.QApplication.topLevelWidgets()
            if isinstance(w, RD.ErrorReportDialog) and w.isVisible()]


def test_help_offers_report_a_problem_and_the_three_choices(window):
    win, _reporter, _sent = window
    help_menu = win.menuBar().actions()[0].menu()
    texts = [a.text() for a in help_menu.actions()]
    assert "Report a problem…" in texts
    assert [a.text() for a in win.error_reports_menu.actions()] == [
        "Ask before sending", "Send without asking", "Never send"]


def test_an_error_asks_once_and_sends_only_on_send(window, qtbot):
    _win, reporter, sent = window
    for d in _dialogs():
        d.close()
    reporter._handle(_report())
    reporter._handle(_report())
    dialogs = _dialogs()
    assert len(dialogs) == 1, "the same error is asked about once"
    assert sent == []
    dialog = dialogs[0]
    assert "ValueError" in dialog.preview.toPlainText()
    dialog.note.setPlainText("I clicked Find groups")
    dialog.send()
    qtbot.waitUntil(lambda: bool(sent), timeout=5000)
    assert sent[0]["user_note"] == "I clicked Find groups"
    qtbot.waitUntil(lambda: "issues/9" in dialog.result.text(), timeout=5000)
    dialog.close()


def test_never_send_means_no_dialog_and_nothing_sent(window):
    _win, reporter, sent = window
    for d in _dialogs():
        d.close()
    R.set_consent("never")
    reporter._handle(_report())
    assert _dialogs() == [] and sent == []


def test_always_send_sends_without_asking(window, qtbot):
    _win, reporter, sent = window
    for d in _dialogs():
        d.close()
    R.set_consent("always")
    reporter._handle(_report())
    qtbot.waitUntil(lambda: bool(sent), timeout=5000)
    assert _dialogs() == []
