"""Error reports from the window: ask, show exactly what would be sent, send (instruction 97).

An uncaught error opens `ErrorReportDialog`: what happened in a sentence, the anonymised report
word for word, an optional note, and Send / Don't send. "Send future reports without asking"
turns asking off; Help ▸ Error reports turns it back on, or off altogether. The same error is
asked about once per session. Help ▸ Report a problem sends a report written by hand.

Sending happens on a thread -- the relay can take seconds -- and its answer (the issue link, or
why it was refused) is shown in the dialog. Nothing here raises: an error reporter that can
itself crash the window is worse than none.
"""
from __future__ import annotations

import sys
import threading

from PySide6 import QtCore, QtGui, QtWidgets

from . import relay as R


class _Sender(QtCore.QThread):
    done = QtCore.Signal(object, str)                   # answer dict or None, error text

    def __init__(self, report: dict, parent=None):
        super().__init__(parent)
        self.report = report

    def run(self):
        try:
            self.done.emit(R.send_report(self.report), "")
        except Exception as error:                      # shown, never raised
            self.done.emit(None, str(error))


#: Senders still running, owned here so a dialog closed mid-send does not destroy one.
_SENDING: set = set()


def _send(report: dict, on_done) -> None:
    sender = _Sender(report)
    _SENDING.add(sender)
    sender.done.connect(on_done)
    sender.finished.connect(lambda s=sender: (_SENDING.discard(s), s.deleteLater()))
    sender.start()


def _answer_text(answer, error: str) -> str:
    if error:
        return f"Not sent: {error}"
    url = (answer or {}).get("issue_url", "")
    seen = " (this error was already reported; your report was added to it)" \
        if (answer or {}).get("deduplicated") else ""
    return f"Sent{seen}. <a href='{url}'>{url}</a>" if url else "Sent."


class ErrorReportDialog(QtWidgets.QDialog):
    """Asks before sending one report, and shows it in full."""

    def __init__(self, report: dict, parent=None, title: str = "ParaScape ran into an error"):
        super().__init__(parent)
        self.report = report
        self.setWindowTitle(title)
        self.setMinimumWidth(620)
        v = QtWidgets.QVBoxLayout(self)
        what = report.get("exception_type", "")
        headline = QtWidgets.QLabel(
            f"<b>{title}.</b><br>{what}: {report.get('message', '')}" if what != "ProblemReport"
            else "<b>Report a problem.</b> Describe what went wrong; the report below is all that "
                 "is sent.")
        headline.setWordWrap(True)
        v.addWidget(headline)
        v.addWidget(QtWidgets.QLabel(
            "Sending files it as an issue on github.com/EinarOlafsson/parascape so it gets fixed. "
            "Paths, names and addresses were removed; this is exactly what would be sent:"))
        self.preview = QtWidgets.QPlainTextEdit(R.preview(report))
        self.preview.setReadOnly(True)
        self.preview.setMinimumHeight(220)
        v.addWidget(self.preview, 1)
        self.note = QtWidgets.QPlainTextEdit()
        self.note.setPlaceholderText("What were you doing? (optional; it is sent with the report)")
        self.note.setMaximumHeight(80)
        v.addWidget(self.note)
        self.always = QtWidgets.QCheckBox("Send future error reports without asking")
        self.always.setToolTip("Change this later in Help ▸ Error reports.")
        if what != "ProblemReport":
            v.addWidget(self.always)
        self.result = QtWidgets.QLabel("")
        self.result.setOpenExternalLinks(True)
        self.result.setWordWrap(True)
        v.addWidget(self.result)
        buttons = QtWidgets.QDialogButtonBox()
        self.send_button = buttons.addButton("Send report", QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        self.send_button.setToolTip("File this report as a ParaScape issue.")
        self.cancel_button = buttons.addButton("Don't send", QtWidgets.QDialogButtonBox.ButtonRole.RejectRole)
        self.cancel_button.setToolTip("Close without sending anything.")
        self.send_button.clicked.connect(self.send)
        self.cancel_button.clicked.connect(self.reject)
        v.addWidget(buttons)

    def send(self) -> None:
        """Send the report (with the note) and show what the relay answered."""
        note = self.note.toPlainText().strip()
        if note:
            if self.report.get("exception_type") == "ProblemReport":
                self.report = R.problem_report(note, self.report.get("screen", ""))
            else:
                self.report["user_note"] = R.scrub(note, R.MAX_NOTE)
        if self.always.isChecked():
            R.set_consent("always")
        self.send_button.setEnabled(False)
        self.result.setText("Sending…")
        _send(self.report, self._sent)

    def _sent(self, answer, error: str) -> None:
        self.result.setText(_answer_text(answer, error))
        self.cancel_button.setText("Close")
        if error:
            self.send_button.setEnabled(True)


class _Reporter(QtCore.QObject):
    """Receives exceptions from any thread and handles them on the GUI thread."""

    caught = QtCore.Signal(object)

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.asked: set = set()
        self.caught.connect(self._handle)
        self._previous = sys.excepthook
        self._previous_thread = threading.excepthook
        sys.excepthook = self._hook
        threading.excepthook = self._thread_hook

    def _hook(self, exc_type, exc, tb):
        try:
            self._previous(exc_type, exc, tb)           # the console still shows it
        finally:
            self._queue(exc_type, exc, tb)

    def _thread_hook(self, args):
        try:
            self._previous_thread(args)
        finally:
            self._queue(args.exc_type, args.exc_value, args.exc_traceback)

    def _queue(self, exc_type, exc, tb):
        if exc_type is None or issubclass(exc_type, (KeyboardInterrupt, SystemExit)):
            return
        try:
            self.caught.emit(R.report_from_exception(exc_type, exc, tb))
        except Exception:                               # never let reporting raise
            pass

    def _handle(self, report: dict) -> None:
        mode = R.consent()
        key = R.signature(report)
        if mode == "never" or key in self.asked:
            return
        self.asked.add(key)
        if mode == "always":
            _send(report, lambda answer, error: self._status(_answer_text(answer, error)))
            return
        dialog = ErrorReportDialog(report, self.window)
        dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.show()

    def _status(self, text: str) -> None:
        try:
            self.window.statusBar().showMessage(f"Error report: {text}", 10000)
        except Exception:
            pass

    def uninstall(self) -> None:
        sys.excepthook = self._previous
        threading.excepthook = self._previous_thread


def report_problem(window) -> ErrorReportDialog:
    """Help ▸ Report a problem: a report written by hand, previewed like any other."""
    dialog = ErrorReportDialog(R.problem_report("(describe the problem below)"), window,
                               title="Report a problem")
    dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)
    dialog.show()
    return dialog


def install(window) -> _Reporter:
    """Catch errors in `window`'s process, and add Help ▸ Report a problem and Help ▸ Error
    reports (ask / always / never)."""
    reporter = _Reporter(window)
    window.error_reporter = reporter
    help_menu = next((a.menu() for a in window.menuBar().actions()
                      if a.text().replace("&", "") == "Help" and a.menu() is not None), None)
    if help_menu is None:
        return reporter
    help_menu.addSeparator()
    act = help_menu.addAction("Report a problem…")
    act.setToolTip("Send a description of a problem as a ParaScape issue. You see exactly what "
                   "is sent before it goes.")
    act.setStatusTip("Report a problem to the ParaScape maintainers")
    act.triggered.connect(lambda: report_problem(window))
    sub = help_menu.addMenu("Error reports")
    sub.setToolTipsVisible(True)
    group = QtGui.QActionGroup(sub)
    group.setExclusive(True)
    for value, label, tip in (
            ("ask", "Ask before sending", "Show each error report and send it only if you click "
                                          "Send."),
            ("always", "Send without asking", "Send anonymised error reports as they happen."),
            ("never", "Never send", "Do not send error reports, and do not ask.")):
        a = QtGui.QAction(label, sub, checkable=True)
        a.setChecked(R.consent() == value)
        a.setToolTip(tip)
        a.setStatusTip(tip)
        a.triggered.connect(lambda _=False, v=value: R.set_consent(v))
        group.addAction(a)
        sub.addAction(a)
    window.error_reports_menu = sub
    return reporter
