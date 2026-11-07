"""The tool modules in the module shell: every remaining tile looks and behaves alike.

Einar: "all modules need to be thematically the same." The gene-reading modules (Map, Gene,
Network, Infer) are built one by one with the parascape-module mould. The rest are TOOLS --
the strategy table, the scorecard, the comparisons, the jobs list -- whose panels already exist
as docks. Each is shown here inside the same shell: the masthead and Home, "Your gene" (which
selects the gene for the panel, as clicking it on the map would), the panel itself as the main
section, then actions and the folded console and AI chat. The panel is LENT to the module while
it is on screen and handed back to its dock after (`shell._Borrowed`), so nothing is copied and
the dock works as before.

`SPECS` is the one list of which tile is which panel; a module rebuilt with the mould later
simply leaves it.
"""
from __future__ import annotations

from dataclasses import dataclass

from PySide6 import QtCore, QtWidgets

from . import theme as TH
from .shell import ModuleShell, _Borrowed


@dataclass(frozen=True)
class Spec:
    """Which panel a tool tile shows, and how the module introduces it."""

    key: str                  # the module key (modules.py)
    title: str
    description: str          # what it is
    instruction: str          # what to do first
    dock: str                 # the dock title whose panel is the main section
    tab: str = ""             # a tab of that panel to bring forward
    console: bool = True      # False for the console itself (it cannot be lent twice)


SPECS = (
    Spec("expression", "Expression", "where and when a gene's RNA and protein are measured",
         "Choose a gene; compare its expression across stages, strains and conditions.",
         "compare", "Gene expression"),
    Spec("screen", "Screen", "your screen joined to everything already recorded",
         "Load your screen table and compare it with the recorded conditions.", "compare"),
    Spec("cell", "Cell", "the parasite drawn, with genes placed where they are found",
         "Choose a compartment or a gene to see where it sits in the cell.", "find"),
    Spec("gallery", "Gallery", "every map built so far, side by side",
         "Click a map to show it; compare maps rather than rebuild them.", "gallery"),
    Spec("strategies", "Strategies", "the ways of turning the data into a claim, each graded",
         "Choose a strategy to see what it does, what it may use and how well it did.",
         "strategies"),
    Spec("scorecard", "Scorecard", "how well each answer did on genes it was not shown",
         "Choose a label to see which methods recover it, and how often they are right.",
         "discoveries"),
    Spec("information", "Information", "the gene x trait map, and how sure each cell is",
         "Search a gene or choose a trait; every cell shows its source and certainty.",
         "information map"),
    Spec("drop", "Drop", "a file of genes, read and resolved",
         "Drop a CSV, TSV or Excel file of genes here, or choose one.", "gene list"),
    Spec("jobs", "Jobs", "what is running, what it is costing, and how to stop it",
         "Select a job to stop it; finished jobs keep their outcome here.", "jobs"),
    Spec("console", "Console", "what the program is doing, in its own words",
         "Filter by level or text; copy what you need for a report.", "console", console=False),
)
BY_KEY = {spec.key: spec for spec in SPECS}


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _dock(window, title: str):
    return next((d for d in window.findChildren(QtWidgets.QDockWidget)
                 if d.windowTitle() == title), None)


class _MainSection(_Borrowed):
    """The lent panel, with its tab brought forward each time it is taken."""

    def __init__(self, dock, where: str, tab: str = "", parent=None):
        super().__init__(dock, where, parent)
        self.tab = tab

    def take(self) -> None:
        super().take()
        if self.tab and self.panel is not None:
            for tabs in [self.panel] + self.panel.findChildren(QtWidgets.QTabWidget):
                if isinstance(tabs, QtWidgets.QTabWidget):
                    for i in range(tabs.count()):
                        if tabs.tabText(i) == self.tab:
                            tabs.setCurrentIndex(i)
                            return


class ToolModule(QtCore.QObject):
    """One tool module: a dock panel in the shell."""

    def __init__(self, window, spec: Spec):
        super().__init__(window)
        self.window, self.spec = window, spec
        dock = _dock(window, spec.dock)
        self.shell = shell = ModuleShell(spec.key, spec.title, spec.description,
                                         instruction=spec.instruction)
        shell.module = self

        gene = QtWidgets.QWidget(shell)
        gl = QtWidgets.QVBoxLayout(gene)
        gl.setContentsMargins(6, 2, 6, 6)
        self.field = QtWidgets.QLineEdit()
        self.field.setPlaceholderText("A gene id or a word of its product, then Return")
        self.field.setToolTip(TH.tip("Select a gene for this panel, as clicking it on the map "
                                     "would."))
        self.field.returnPressed.connect(self._select)
        gl.addWidget(self.field)
        gl.addWidget(_muted(f"{spec.title}: {spec.description}."))
        shell.add_setting("Your gene", gene)

        self.main = _MainSection(dock, f"the {spec.title} module", spec.tab, shell)
        shell.add_section(spec.title, self.main, extent=620, stretch=4)
        shell._borrowed.append(self.main)

        shell.add_action("Read it in Gene", self._to_gene, primary=True,
                         tip="Open the selected gene in the Gene module: what is known, with sources.")
        shell.add_action("Open as a panel", self._as_panel,
                         tip=f"Leave the module and use {spec.title} as a panel beside the map.")
        if spec.console:
            shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                                 ("AI chat", getattr(window, "chat_dock", None))],
                         folded=True)
        else:
            shell.borrow("AI chat", [("AI chat", getattr(window, "chat_dock", None))], folded=True)
        shell.home_requested.connect(lambda: _go_home(window))

    def _select(self) -> None:
        text = self.field.text().strip()
        if text and hasattr(self.window, "_workflow_gene"):
            try:
                self.window._workflow_gene(text)
                self.shell.say(f"{text} selected")
            except Exception as error:
                self.shell.say(f"could not select {text}: {error}")

    def _to_gene(self) -> None:
        text = self.field.text().strip()
        if not text:
            sel = getattr(self.window, "sel", None)
            if sel is not None:
                text = str(self.window.nodes.gene_id.iloc[int(sel)])
        from .gene_module import open_in
        open_in(self.window, text)

    def _as_panel(self) -> None:
        from .landing_install import leave_landing, raise_dock
        leave_landing(self.window)
        raise_dock(self.window, self.spec.dock, self.spec.tab)


def build_all(window, landing) -> dict:
    """Build every tool module whose dock exists, add each to the landing, and return
    {key: shell}."""
    out = {}
    for spec in SPECS:
        if _dock(window, spec.dock) is None:
            continue
        try:
            shell = ToolModule(window, spec).shell
        except Exception as error:                   # one tool must not stop the others
            window.statusBar().showMessage(f"{spec.title} module unavailable: {error}", 8000)
            continue
        landing.stack.addWidget(shell)
        out[spec.key] = shell
    return out


def open_in(window, key: str) -> str:
    """Show a tool module over the window."""
    from .landing_install import enter_landing
    shell = getattr(window, "tool_modules", {}).get(key)
    landing = getattr(window, "landing", None)
    if shell is None or landing is None:
        return f"the {key} module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return BY_KEY[key].title


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
