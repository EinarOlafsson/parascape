"""The tool modules in the module shell: every remaining tile looks and behaves alike.

Einar: "all modules need to be thematically the same." The gene-reading modules (Map, Gene,
Network, Infer) are built one by one with the parascape-module mould. The rest are TOOLS --
the strategy table, the scorecard, the comparisons, the jobs list -- whose panels already exist
as docks. Each is shown here inside the same shell: the masthead and Home, "Your gene" for a
module that reads one (which selects the gene for the panel, as clicking it on the map would),
the panel itself as the main section, then actions and the folded console and AI chat. The panel
is LENT to the module while it is on screen and handed back to its dock after
(`shell._Borrowed`), so nothing is copied and the dock works as before.

THE SUBJECT FIELD IS NOT GIVEN TO EVERY MODULE (instruction 98, WP11)
=====================================================================
It was, which meant Jobs asked for "a gene id or a word of its product" above the list of what
is running, and Console, Gallery, Drop, Strategies and Scorecard did the same. None of the six
declares that it reads a gene. The field is now built from `Module.consumes`, so a tool module
offers a control for exactly the subject its ports say it takes.

`SPECS` is the one list of which tile is which panel; a module rebuilt with the mould later
simply leaves it.

COMPUTED FROM THE CATALOGUE (instruction 98, WP3)
=================================================
`SPECS` used to be written out here: ten keys, each with its own title, its own description
and its own dock name, beside the same ten rows in `modules.py` carrying the same three
things in slightly different words. Two tables saying one thing is how a module ends up
named one way on its tile and another on its own masthead. The rows are now derived from
the catalogue -- a module joins this list by setting `builder="tool"` there, and leaves it
by setting a shell builder instead.
"""
from __future__ import annotations

from dataclasses import dataclass

from PySide6 import QtCore, QtWidgets

from . import modules as MOD
from . import ports as PORTS
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


def _spec(module) -> Spec:
    """Turn one catalogue row into the spec this file builds a module from.

    The dock and the tab come from `opens`, which is the destination the window
    was checked against (`tests/test_modules_in_window.py`); naming them again
    here is what let six of them drift to views that did not exist.

    :param module: a `modules.Module` whose builder is ``"tool"``.
    """
    dock, _, tab = module.opens.partition(":")
    # The console cannot lend itself its own console: `shell.borrow` would take
    # the panel the module IS showing and put it in the folded section under it.
    return Spec(module.key, module.name, module.note, module.instruction,
                dock, tab, console=module.key != "console")


#: Which tile is which panel, derived from `modules.MODULES`.
SPECS = tuple(_spec(module) for module in MOD.tool_modules())
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

        # The gene field only where the module SAYS IT READS A GENE (instruction 98, WP11).
        # Every tool module used to get one, so Jobs offered "Your gene: a gene id or a word
        # of its product" above a list of running jobs, and Console, Gallery, Drop, Strategies
        # and Scorecard did the same -- six modules with a control for a subject their ports do
        # not declare, which is the "tile that says one thing and does another" this
        # instruction exists to remove. The answer comes from the declared ports rather than
        # from a list of which six, so a module that starts reading a gene gets the field by
        # declaring it. `tests/test_module_contract.py` holds both directions.
        self.field = None
        if PORTS.GENE in MOD.get(spec.key).consumes:
            gene = QtWidgets.QWidget(shell)
            gl = QtWidgets.QVBoxLayout(gene)
            gl.setContentsMargins(6, 2, 6, 6)
            self.field = QtWidgets.QLineEdit()
            self.field.setPlaceholderText("A gene id or a word of its product, then Return")
            self.field.setToolTip(TH.tip("Select a gene for this panel, as clicking it on the "
                                         "map would."))
            self.field.returnPressed.connect(self._select)
            gl.addWidget(self.field)
            gl.addWidget(_muted(f"{spec.title}: {spec.description}."))
            shell.add_setting("Your gene", gene, subject=PORTS.GENE)

        self.main = _MainSection(dock, f"the {spec.title} module", spec.tab, shell)
        shell.add_section(spec.title, self.main, extent=620, stretch=4)
        shell._borrowed.append(self.main)

        if self.field is not None:
            shell.add_action("Read it in Gene", self._to_gene, primary=True,
                             tip="Open the selected gene in the Gene module: what is known, "
                                 "with sources.")
        shell.add_action("Open as a panel", self._as_panel,
                         primary=self.field is None,
                         tip=f"Leave the module and use {spec.title} as a panel beside the map.")
        if spec.console:
            shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                                 ("AI chat", getattr(window, "chat_dock", None))],
                         folded=True)
        else:
            shell.borrow("AI chat", [("AI chat", getattr(window, "chat_dock", None))], folded=True)
        shell.home_requested.connect(lambda: _go_home(window))

    def _select(self) -> None:
        text = self.field.text().strip() if self.field is not None else ""
        if text and hasattr(self.window, "_workflow_gene"):
            try:
                self.window._workflow_gene(text)
                self.shell.say(f"{text} selected")
            except Exception as error:
                self.shell.say(f"could not select {text}: {error}")

    def _to_gene(self) -> None:
        text = self.field.text().strip() if self.field is not None else ""
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
