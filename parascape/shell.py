"""spaCR's module shell, for ParaScape's modules (instruction 86).

Einar: "spaCR has its settings in a container on the left, a console, AI chatbox, and an
actions section. There is a live section to live test parameters and a hyperparameter
search ... spaCR containers are also collapsible so after setting the settings you can
give more space to the content."

The layout is spaCR's `AppScreen` (`spacr/qt/screens/app_screen.py`), rebuilt rather than
copied -- that file is 14,000 lines, nearly all of it spaCR's own modules -- on spaCR's
own pane widgets, which ARE copied (:mod:`parascape.panes`):

* a masthead: the module name large and light, what it is for, what to do first, and a
  way back to the home screen;
* a horizontal splitter: **Settings** on the left, an EDGE pane that folds to nothing
  against its handle, and the **runtime column** on the right;
* the runtime column, a vertical splitter of sections that each fold to their heading:
  whatever the module shows (a live preview, a table), the hyperparameter search,
  **Actions**, **Console** and **AI chat**.

Console and AI chat are the window's own panels, borrowed while the module is on screen
and handed back to their docks when it is not -- one console and one conversation, not
one per module, which is what the docks already were.
"""
from __future__ import annotations

from PySide6 import QtCore, QtWidgets

from . import theme as TH
from .panes.collapsible_section import CollapsibleSection
from .panes.collapsible_splitter import EDGE, CollapsibleSplitter

#: spaCR's settings column opens at 400 px beside an 800 px runtime column.
SETTINGS_EXTENT, RUNTIME_EXTENT = 400, 800


def shell_qss(palette: dict) -> str:
    """The shell's stylesheet: spaCR's masthead type and translucent section surfaces."""
    p = palette
    return f"""
    QWidget#ModuleShell {{ background: transparent; }}
    QLabel#ShellTitle {{ font-size: 30px; font-weight: 300; color: {p['fg']};
                         background: transparent; }}
    QLabel#ShellMuted {{ color: {p['fg_dim']}; background: transparent; }}
    QWidget#ShellSettings, QWidget#ShellRuntimePage {{ background: transparent; }}
    QScrollArea#ShellSettingsScroll, QScrollArea#ShellSettingsScroll > QWidget > QWidget {{
        background: {TH.rgba(p['surface'], 0.82)}; border: 1px solid {p['border_soft']};
        border-radius: 12px; }}
    QWidget#ShellSection {{ background: {TH.rgba(p['surface'], 0.82)};
                            border: 1px solid {p['border_soft']}; border-radius: 12px; }}
    QWidget#ModuleShell QLineEdit, QWidget#ModuleShell QSpinBox,
    QWidget#ModuleShell QDoubleSpinBox {{ min-height: 20px; }}
    QLabel#ShellStatus {{ color: {p['fg_dim']}; background: transparent; }}
    QPushButton#ShellHome {{ background: transparent; border: 1px solid {p['border']};
                             border-radius: 14px; padding: 4px 14px; color: {p['fg']}; }}
    QPushButton#ShellHome:hover {{ border-color: {p['accent']}; color: {p['accent']}; }}
    """


class _Borrowed(QtWidgets.QWidget):
    """A place in the shell for a dock's panel, which it holds only while it is shown.

    `QDockWidget.setWidget` leaves the old widget parented to the dock, so the panel is
    moved explicitly both ways, and the dock shows a one-line note while it is away
    rather than an empty frame.
    """

    def __init__(self, dock: QtWidgets.QDockWidget, where: str, parent=None):
        super().__init__(parent)
        self.dock = dock
        self.panel = None
        self._layout = QtWidgets.QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._note = QtWidgets.QLabel(f"Open in {where}. It comes back here when you leave it.")
        self._note.setObjectName("ShellMuted")
        self._note.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)

    def take(self) -> None:
        """Move the dock's panel here. Does nothing if it is already here."""
        if self.panel is not None or self.dock is None:
            return
        panel = self.dock.widget()
        if panel is None or panel is self._note:
            return
        self.dock.setWidget(self._note)
        self._layout.addWidget(panel)
        panel.show()
        self.panel = panel

    def give_back(self) -> None:
        """Return the panel to its dock."""
        panel, self.panel = self.panel, None
        if panel is None or self.dock is None:
            return
        self._layout.removeWidget(panel)
        # Out of the dock before the panel goes back in, and parked here hidden rather than
        # unparented: an unparented label is a stray top-level window, and deleting those at
        # teardown alongside the dock it was taken from crashed the process.
        self._note.setParent(self)
        self._note.hide()
        self.dock.setWidget(panel)


class ModuleShell(QtWidgets.QWidget):
    """One module in spaCR's shell. Modules fill it through `add_setting`, `add_section`
    and `add_action`; nothing is drawn for a part a module does not use."""

    home_requested = QtCore.Signal()

    def __init__(self, key: str, title: str, description: str = "", instruction: str = "",
                 parent=None):
        super().__init__(parent)
        self.key = key
        self.setObjectName("ModuleShell")
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self._borrowed: list[_Borrowed] = []
        self.sections: dict[str, QtWidgets.QWidget] = {}
        self.actions: dict[str, QtWidgets.QPushButton] = {}

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(24, 20, 24, 16)
        outer.setSpacing(12)

        head = QtWidgets.QHBoxLayout()
        head.setSpacing(16)
        names = QtWidgets.QVBoxLayout()
        names.setSpacing(2)
        self.title_label = QtWidgets.QLabel(title)
        self.title_label.setObjectName("ShellTitle")
        names.addWidget(self.title_label)
        self.instruction_label = QtWidgets.QLabel(instruction)
        self.instruction_label.setObjectName("ShellMuted")
        self.instruction_label.setVisible(bool(instruction))
        names.addWidget(self.instruction_label)
        head.addLayout(names)
        self.description_label = QtWidgets.QLabel(description)
        self.description_label.setObjectName("ShellMuted")
        if description:
            self.description_label.setToolTip(TH.tip(description))
        self.description_label.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum,
                                             QtWidgets.QSizePolicy.Policy.Preferred)
        self.description_label.setMinimumWidth(0)
        head.addWidget(self.description_label)
        head.addStretch(1)
        self.trailing = QtWidgets.QHBoxLayout()
        head.addLayout(self.trailing)
        self.home_button = QtWidgets.QPushButton("Home")
        self.home_button.setObjectName("ShellHome")
        self.home_button.setToolTip("Back to the home screen (Alt+Home).")
        self.home_button.clicked.connect(self.home_requested)
        head.addWidget(self.home_button)
        outer.addLayout(head)

        self.body = CollapsibleSplitter(QtCore.Qt.Orientation.Horizontal,
                                        persist_key=f"{key}::body")
        self.body.setChildrenCollapsible(False)

        self.settings = QtWidgets.QWidget()
        self.settings.setObjectName("ShellSettings")
        self._settings_layout = QtWidgets.QVBoxLayout(self.settings)
        self._settings_layout.setContentsMargins(10, 10, 10, 10)
        self._settings_layout.setSpacing(6)
        self._settings_layout.addStretch(1)
        scroll = QtWidgets.QScrollArea()
        scroll.setObjectName("ShellSettingsScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setWidget(self.settings)
        self.settings_scroll = scroll
        self.body.add_pane(scroll, "Settings", mode=EDGE, fold_key=f"{key}/Settings",
                           stretch=1, extent=SETTINGS_EXTENT)

        self.runtime = CollapsibleSplitter(QtCore.Qt.Orientation.Vertical,
                                           persist_key=f"{key}::runtime")
        self.runtime.setObjectName("ShellRuntimePage")
        self.body.add_pane(self.runtime, "Runtime", stretch=2, extent=RUNTIME_EXTENT)
        self.body.setSizes([SETTINGS_EXTENT, RUNTIME_EXTENT])
        outer.addWidget(self.body, 1)

        self._actions_page = QtWidgets.QWidget()
        actions = QtWidgets.QVBoxLayout(self._actions_page)
        actions.setContentsMargins(8, 4, 8, 8)
        self._action_row = QtWidgets.QHBoxLayout()
        self._action_row.setSpacing(8)
        self._action_row.addStretch(1)
        actions.addLayout(self._action_row)
        self.status_line = QtWidgets.QLabel("")
        self.status_line.setObjectName("ShellStatus")
        self.status_line.setWordWrap(True)
        actions.addWidget(self.status_line)
        self._actions_added = False

        self._palette = TH.palette_for("dark")
        self.set_palette(self._palette)

    # ------------------------------------------------------------------ filling it
    def set_palette(self, palette: dict) -> None:
        """Restyle the shell for a theme's palette."""
        self._palette = palette
        self.setStyleSheet(shell_qss(palette))

    def add_setting(self, name: str, widget: QtWidgets.QWidget,
                    expanded: bool = True) -> CollapsibleSection:
        """A foldable group in the Settings column."""
        section = CollapsibleSection(name, widget, expanded=expanded)
        # spaCR's header colours itself `palette(text)`, which over the landing's backdrop is the
        # dark default palette and left the headings near-invisible; the theme's own foreground.
        header = section._header
        header.setStyleSheet(header.styleSheet().replace("palette(text)", self._palette["fg"]))
        self._settings_layout.insertWidget(self._settings_layout.count() - 1, section)
        self.sections[name] = section
        return section

    def add_section(self, name: str, widget: QtWidgets.QWidget, extent: int = 240,
                    stretch: int = 1, focus: bool = False):
        """A section of the runtime column that folds to its heading."""
        widget.setObjectName(widget.objectName() or "ShellSection")
        section = self.runtime.add_section(widget, name, persist_key=f"{self.key}/{name}",
                                           extent=extent, stretch=stretch, focus=focus)
        self.sections[name] = section
        return section

    def add_action(self, label: str, slot, tip: str = "",
                   primary: bool = False) -> QtWidgets.QPushButton:
        """A button in the Actions section, which appears with the first one."""
        button = QtWidgets.QPushButton(label)
        if primary:
            button.setProperty("primary", True)
        if tip:
            button.setToolTip(TH.tip(tip))
        button.clicked.connect(slot)
        self._action_row.insertWidget(self._action_row.count() - 1, button)
        self.actions[label] = button
        if not self._actions_added:
            self._actions_added = True
            self.add_section("Actions", self._actions_page, extent=90, stretch=0)
        return button

    def add_action_widget(self, widget: QtWidgets.QWidget) -> None:
        """A ready-made control in the Actions row (a Stop button that already exists)."""
        self._action_row.insertWidget(self._action_row.count() - 1, widget)

    def borrow(self, name: str, docks, extent: int = 200) -> list:
        """Show docks' panels (the console, the assistant) side by side in one section while
        the module is on screen -- spaCR's console and chat share a row the same way, so the
        two do not take the column's height twice."""
        docks = [(title, d) for title, d in docks if d is not None]
        if not docks:
            return []
        row = CollapsibleSplitter(QtCore.Qt.Orientation.Horizontal,
                                  persist_key=f"{self.key}::{name}")
        holders = []
        for title, dock in docks:
            holder = _Borrowed(dock, f"the {self.title_label.text()} module")
            row.add_pane(holder, title, stretch=1, extent=400)
            holders.append(holder)
        self.add_section(name, row, extent=extent, stretch=1)
        self._borrowed.extend(holders)
        return holders

    def say(self, text: str) -> None:
        """Put a line under the actions."""
        self.status_line.setText(str(text))

    # ------------------------------------------------------------------ on screen
    def showEvent(self, event):
        for holder in self._borrowed:
            holder.take()
        super().showEvent(event)

    def hideEvent(self, event):
        for holder in self._borrowed:
            holder.give_back()
        super().hideEvent(event)
