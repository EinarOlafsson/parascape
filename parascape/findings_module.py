"""The Findings module: every analysis behind the paper, runnable here.

Instruction 101, items B2 and B3. The analyses are corpus-level -- they are about what the
published literature can and cannot say -- so this sits in the Data band beside Datasets and
Information rather than in a module that takes a gene.

WHAT IT IS FOR. A reader of the paper should be able to open the program, pick a figure, run the
analysis that produced it, and see both the answer and the record of what the answer was made
from. A result with no record of its inputs is a number somebody has to take on trust.

WHAT IT REFUSES TO DO. It never shows a cached number as though it had just been computed: a
finding that has not been run in this installation says so, with the date of the result on disk.
The claim shown beside each finding is the sentence about the parasite that the analysis ends on,
not a description of the method.
"""
from __future__ import annotations

import datetime as dt
import json
import os

from PySide6 import QtCore, QtWidgets

from . import paper_findings as PF
from . import ports as PORTS
from . import theme as TH
from .recipes_panel import widget_setting
from .shell import ModuleShell, solid_dropdown


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


class FindingsModule(QtCore.QObject):
    """The Findings module for one window."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.process: QtCore.QProcess | None = None
        self.shell = shell = ModuleShell.for_module("findings")
        shell.module = self
        self._build_settings()
        self._build_runtime()
        self._describe()
        shell.home_requested.connect(lambda: _go_home(window))

    # ------------------------------------------------------------------ the left column
    def _build_settings(self) -> None:
        shell = self.shell
        which = QtWidgets.QWidget(shell)
        wl = QtWidgets.QVBoxLayout(which)
        wl.setContentsMargins(6, 2, 6, 6)
        self.finding = QtWidgets.QComboBox()
        self.finding.setToolTip(TH.tip(
            "Each entry is one analysis of the paper. Choosing one shows the claim it ends on "
            "and when it was last run here."))
        solid_dropdown(self.finding)
        for f in PF.FINDINGS:
            self.finding.addItem(f"Figure {f.figure} — {f.name}", f.key)
        self.finding.currentIndexChanged.connect(self._describe)
        wl.addWidget(self.finding)
        self.claim = _muted("")
        wl.addWidget(self.claim)
        shell.add_setting("Which finding", which, subject=PORTS.QUESTION)

        where = QtWidgets.QWidget(shell)
        el = QtWidgets.QVBoxLayout(where)
        el.setContentsMargins(6, 2, 6, 6)
        self.cache = QtWidgets.QLineEdit(os.environ.get("PARASCAPE_CACHE", ""))
        self.cache.setToolTip(TH.tip("The data release the analysis reads."))
        self.omics = QtWidgets.QLineEdit(os.environ.get("PARASCAPE_OMICS", ""))
        self.omics.setToolTip(TH.tip(
            "The omics store, for the findings that read per-study matrices rather than the "
            "release alone. Left empty, those findings use the release."))
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.addRow("Data release", self.cache)
        form.addRow("Omics store", self.omics)
        holder = QtWidgets.QWidget(where)
        holder.setLayout(form)
        el.addWidget(holder)
        shell.add_setting("Where the data is", where, expanded=False)
        shell.add_recipes(self.recipe_settings())

    def recipe_settings(self) -> list:
        """What a recipe keeps. The finding itself is the SUBJECT and is not one of them."""
        return [widget_setting("cache", self.cache, "data release"),
                widget_setting("omics", self.omics, "omics store")]

    # ----------------------------------------------------------------- the right column
    def _build_runtime(self) -> None:
        shell = self.shell
        self.reading = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.reading)
        shell.add_section("What this finding says", self.reading, extent=230, stretch=2,
                          focus=True)

        self.made_from = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.made_from)
        shell.add_section("What it was made from", self.made_from, extent=220, stretch=2)

        self.output = QtWidgets.QPlainTextEdit(shell)
        self.output.setReadOnly(True)
        self.output.setPlaceholderText("Run it and what the analysis prints appears here.")
        shell.add_section("While it runs", self.output, extent=140, stretch=1)

        self.run_button = shell.add_action(
            "Run this analysis", self.run, primary=True,
            tip="Run the script that produced this figure, in its own process, and show what it "
                "printed. Nothing already on disk is overwritten until it finishes.")
        shell.add_action("Stop", self.stop, tip="Stop a running analysis.")
        shell.borrow("Console and AI chat",
                     [("Console", getattr(self.window, "console_dock", None)),
                      ("AI chat", getattr(self.window, "chat_dock", None))],
                     folded=True)

    # ------------------------------------------------------------------------ what it says
    def chosen(self) -> PF.Finding:
        return PF.get(self.finding.currentData())

    def _describe(self) -> None:
        """The claim, and whether this installation has ever produced it."""
        finding = self.chosen()
        self.claim.setText(finding.claim)
        record = PF.record(finding)
        directory = record.get("result_directory", "")
        if not directory:
            self.reading.setHtml(
                f"<p><b>{finding.claim}</b></p><p>This analysis has not been run in this "
                f"installation. The claim above is the one the paper reports; run it to see what "
                f"it gives here.</p>")
            self.made_from.setHtml("<p>Nothing on disk yet.</p>")
            return
        when = dt.datetime.fromtimestamp(os.path.getmtime(directory)).strftime("%-d %B %Y")
        self.reading.setHtml(
            f"<p><b>{finding.claim}</b></p>"
            f"<p>Produced by <code>scripts/{finding.script}</code>, last run here on {when}.</p>")
        rows = "".join(f"<li><b>{k}</b>: {json.dumps(v)[:160]}</li>"
                       for k, v in record.items() if k != "result_directory")
        self.made_from.setHtml(
            f"<p>From <code>{os.path.relpath(directory, PF.ROOT)}</code></p><ul>{rows}</ul>")

    # ---------------------------------------------------------------------------- running
    def run(self) -> None:
        """Run the chosen analysis in its own process, streaming what it prints."""
        if self.process is not None:
            self.shell.say("an analysis is already running; stop it first")
            return
        finding = self.chosen()
        self.output.clear()
        self.run_button.setEnabled(False)
        environment = QtCore.QProcessEnvironment.systemEnvironment()
        if self.cache.text().strip():
            environment.insert("PARASCAPE_CACHE", self.cache.text().strip())
        if self.omics.text().strip():
            environment.insert("PARASCAPE_OMICS", self.omics.text().strip())
        self.process = QtCore.QProcess(self)
        self.process.setProcessEnvironment(environment)
        self.process.setWorkingDirectory(PF.ROOT)
        self.process.setProcessChannelMode(QtCore.QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read)
        self.process.finished.connect(self._finished)
        import sys
        self.process.start(sys.executable, [PF.script_path(finding)])
        self.shell.say(f"running {finding.script}")

    def _read(self) -> None:
        if self.process is None:
            return
        text = bytes(self.process.readAllStandardOutput()).decode("utf-8", "replace")
        self.output.appendPlainText(text.rstrip("\n"))

    def _finished(self, code: int, _status) -> None:
        self.run_button.setEnabled(True)
        self.process = None
        if code == 0:
            self.shell.say("finished; the result and what it was made from are shown above")
            self._describe()
        else:
            self.shell.say(f"the analysis stopped with code {code}; see what it printed")

    def stop(self) -> None:
        if self.process is None:
            self.shell.say("nothing is running")
            return
        self.process.kill()
        self.shell.say("stopped")


def build(window) -> ModuleShell | None:
    """The Findings module for `window`."""
    return FindingsModule(window).shell


def open_in(window) -> str:
    """Show the Findings module over the window."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "findings_module", None)
    if landing is None or shell is None:
        return "the Findings module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return "Findings"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
