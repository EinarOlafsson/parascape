"""The Map module in spaCR's module shell -- the first module moved onto it (instruction 86).

What it holds was the analysis panel's "2 · Map" tab: the projection settings, the
hyperparameter walk and its table, building and saving a map. Those controls are the
panel's own widgets, MOVED here rather than copied, so `AnalysisPanel.spec()`,
`run_embed` and `run_umap_walk` read exactly what they read before and nothing can
drift between two copies. The tab keeps a button that opens this module.

New, because the shell has a place for it and spaCR has one: a **live preview**. spaCR
re-segments a sample image as its settings change; here a sample of genes is projected
in 2D with the n_neighbors, min_dist and seed being set, so the effect of a setting is
seen before a full 3D map is built. It is a preview -- a sample, two dimensions, drawn
flat -- and says so on it.
"""
from __future__ import annotations

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from . import theme as TH
from .shell import ModuleShell

#: Genes in the live preview. Enough for the structure to show, few enough that a
#: projection takes a second or two once umap has compiled.
PREVIEW_GENES = 1500
#: How long the settings must rest before the preview is recomputed, in ms.
PREVIEW_DEBOUNCE_MS = 600
#: Where the Map tab used to be.
MAP_TAB = "2 · Map"


#: Projections still running. A QThread destroyed while it runs aborts the process, and a
#: preview can be closed (or the window quit) in the middle of one, so a projection is
#: owned here, not by its widget, until it has finished.
_RUNNING: set = set()


def wait_for_projections(timeout_ms: int = 60000) -> None:
    """Block until every running projection has finished (on quit, and in tests)."""
    for worker in list(_RUNNING):
        worker.wait(timeout_ms)


class _Projection(QtCore.QThread):
    """One 2D projection of a sample, off the GUI thread."""

    ready = QtCore.Signal(int, object, object)        # generation, xy, colour codes
    failed = QtCore.Signal(int, str)

    def __init__(self, generation, nodes, spec, size, parent=None):
        super().__init__(parent)
        self.generation, self.nodes, self.spec, self.size = generation, nodes, spec, size

    def run(self):
        try:
            from .embedding import build_matrix
            X, _names, rows = build_matrix(self.nodes, self.spec, log=lambda *_a, **_k: None)
            rng = np.random.default_rng(self.spec.random_state)
            pick = np.sort(rng.choice(len(X), size=min(self.size, len(X)), replace=False))
            sample = X[pick]
            import umap
            reducer = umap.UMAP(n_components=2,
                                n_neighbors=int(min(self.spec.n_neighbors, len(sample) - 1)),
                                min_dist=float(self.spec.min_dist),
                                random_state=int(self.spec.random_state))
            xy = reducer.fit_transform(sample)
            codes = None
            if "compartment" in self.nodes.columns:
                kept = np.asarray(rows)
                kept = np.flatnonzero(kept) if kept.dtype == bool else kept
                values = self.nodes["compartment"].to_numpy()[kept[pick]]
                codes = np.unique(values.astype(str), return_inverse=True)[1]
            self.ready.emit(self.generation, xy, codes)
        except Exception as error:                    # a preview must never take the module
            self.failed.emit(self.generation, f"{type(error).__name__}: {error}")


class ScatterView(QtWidgets.QWidget):
    """Points in two dimensions, drawn to fit with their aspect kept.

    Not a pyqtgraph PlotWidget: deleting one of those while the application keeps running
    -- which is what closing a window does -- crashed PySide6 inside `QGraphicsScene`'s
    destructor. A preview of a couple of thousand points needs none of what it brings.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.xy = np.zeros((0, 2))
        self.colours: list = []
        self.setMinimumHeight(180)

    def set_points(self, xy, colours) -> None:
        """Show `xy` (n x 2), one colour per point."""
        self.xy = np.asarray(xy, dtype=float).reshape(-1, 2)
        self.colours = [QtGui.QColor(*_rgb255(c)) for c in colours]
        self.update()

    def paintEvent(self, _event):
        if not len(self.xy):
            return
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        lo, hi = self.xy.min(axis=0), self.xy.max(axis=0)
        span = np.maximum(hi - lo, 1e-9)
        margin = 10
        scale = min((self.width() - 2 * margin) / span[0], (self.height() - 2 * margin) / span[1])
        x0 = (self.width() - span[0] * scale) / 2
        y0 = (self.height() - span[1] * scale) / 2
        radius = 2.0
        for (x, y), colour in zip(self.xy, self.colours):
            painter.setBrush(colour)
            painter.drawEllipse(QtCore.QPointF(x0 + (x - lo[0]) * scale,
                                               y0 + (hi[1] - y) * scale), radius, radius)
        painter.end()


class LivePreview(QtWidgets.QWidget):
    """A flat projection of a sample of genes, redrawn as the map settings change."""

    def __init__(self, panel, parent=None):
        super().__init__(parent)
        self.panel = panel
        self.generation = 0
        self._worker = None
        self._pending = False
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 8)
        row = QtWidgets.QHBoxLayout()
        self.live = QtWidgets.QCheckBox("Live")
        self.live.setChecked(True)
        self.live.setToolTip(TH.tip(
            "Recompute the preview whenever n_neighbors, min_dist or the seed changes. Off, "
            "the preview stays as it is until you turn it on again."))
        self.size = QtWidgets.QSpinBox()
        self.size.setRange(200, 5000)
        self.size.setSingleStep(250)
        self.size.setValue(PREVIEW_GENES)
        self.size.setSuffix(" genes")
        self.size.setToolTip(TH.tip("How many genes the preview projects. More shows more "
                                    "structure and takes longer; the full map uses them all."))
        self.note = QtWidgets.QLabel("")
        self.note.setObjectName("ShellMuted")
        row.addWidget(self.live)
        row.addWidget(self.size)
        row.addWidget(self.note, 1)
        layout.addLayout(row)
        self.plot = ScatterView(self)
        self.plot.setToolTip(TH.tip(
            "A preview: a sample of genes in two dimensions with the settings on the left. "
            "Build this map draws all of them in the 3D view."))
        layout.addWidget(self.plot, 1)

        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(PREVIEW_DEBOUNCE_MS)
        self.timer.timeout.connect(self.refresh)
        for box in (panel.nn, panel.md, panel.seed):
            box.valueChanged.connect(self.changed)
        self.size.valueChanged.connect(self.changed)
        self.live.toggled.connect(lambda on: self.changed() if on else None)

    def changed(self, *_):
        """A setting moved: recompute once the settings have rested, if Live is on."""
        if self.live.isChecked() and self.isVisible():
            self.timer.start()

    def showEvent(self, event):
        super().showEvent(event)
        if self.generation == 0 and self.live.isChecked():
            self.timer.start()

    def hideEvent(self, event):
        # A projection asked for by a timer that fires after the module has gone -- or while
        # the window is being torn down -- is work nobody will see.
        self.timer.stop()
        self._pending = False
        super().hideEvent(event)

    def refresh(self) -> None:
        """Project a sample with the current settings. One projection at a time; a change
        that arrives during one is done as soon as it finishes."""
        if not self.isVisible():
            return
        if self._worker is not None:
            self._pending = True
            return
        try:
            spec = self.panel.spec()
        except Exception as error:
            self.note.setText(f"no preview: {error}")
            return
        self.generation += 1
        first = " -- the first one also compiles umap, up to a minute" \
            if self.generation == 1 else ""
        self.note.setText(f"projecting {self.size.value():,} genes "
                          f"(n_neighbors {spec.n_neighbors}, min_dist {spec.min_dist:g})…{first}")
        worker = _Projection(self.generation, self.panel.nodes, spec, self.size.value())
        _RUNNING.add(worker)
        worker.finished.connect(lambda w=worker: _RUNNING.discard(w))
        worker.ready.connect(self._drawn)
        worker.failed.connect(lambda _g, message: self.note.setText(f"no preview: {message}"))
        worker.finished.connect(self._done)
        self._worker = worker
        worker.start()

    def _drawn(self, generation: int, xy, codes) -> None:
        if generation != self.generation:
            return
        colours = TH.categorical_colors(int(codes.max()) + 1 if codes is not None else 1)
        self.plot.set_points(xy, [colours[c] for c in codes] if codes is not None
                             else [colours[0]] * len(xy))
        self.note.setText(f"preview · {len(xy):,} genes · coloured by compartment")

    def _done(self) -> None:
        worker, self._worker = self._worker, None
        if worker is not None:
            worker.deleteLater()
        if self._pending:
            self._pending = False
            self.refresh()

    def stop(self) -> None:
        """Stop asking for projections and wait for the one running, if any."""
        self.timer.stop()
        self._pending = False
        wait_for_projections()


def _rgb255(colour):
    """A categorical colour as 0-255 RGBA, whichever form the theme returns."""
    if isinstance(colour, str):
        from PySide6.QtGui import QColor
        q = QColor(colour)
        return q.red(), q.green(), q.blue(), 220
    values = list(colour)
    if values and max(values[:3]) <= 1.0:
        values = [int(round(v * 255)) for v in values]
    return (*values[:3], 220)


def _moved(parent_layout_widget: QtWidgets.QWidget, rows) -> QtWidgets.QWidget:
    """A form of (label, widget) rows, the widgets reparented into it."""
    page = QtWidgets.QWidget(parent_layout_widget)
    form = QtWidgets.QFormLayout(page)
    form.setContentsMargins(6, 2, 6, 6)
    for label, widget in rows:
        form.addRow(label, widget)
    return page


def build(window) -> ModuleShell | None:
    """The Map module for `window`, with the Map tab's controls moved into it.

    None when the window has no analysis panel (it needs scikit-learn and umap-learn),
    and the Map tile then keeps opening the 3D view directly.
    """
    panel = getattr(window, "panel", None)
    if panel is None or not hasattr(panel, "nn"):
        return None
    shell = ModuleShell(
        "map", "Map", "the genome as a space you choose the axes of",
        instruction="Set the projection on the left, watch the preview, then build the map.")

    shell.add_setting("Projection", _moved(shell, [
        ("n_neighbors", panel.nn), ("min_dist", panel.md), ("seed", panel.seed)]))
    data = QtWidgets.QWidget(shell)
    dl = QtWidgets.QVBoxLayout(data)
    dl.setContentsMargins(6, 2, 6, 6)
    which = QtWidgets.QLabel("Which measurements go into the map, how they are scaled and what "
                             "happens to missing values is set in analysis ▸ 1 · Data, because "
                             "clustering and inference read the same choice.")
    which.setWordWrap(True)
    which.setObjectName("ShellMuted")
    open_data = QtWidgets.QPushButton("Open 1 · Data")
    open_data.setToolTip("Leave this module and show the data settings.")
    open_data.clicked.connect(lambda: _open_data_tab(window))
    dl.addWidget(which)
    dl.addWidget(open_data, 0, QtCore.Qt.AlignmentFlag.AlignLeft)
    shell.add_setting("Data", data)
    save = QtWidgets.QWidget(shell)
    sl = QtWidgets.QHBoxLayout(save)
    sl.setContentsMargins(6, 2, 6, 6)
    save_button = QtWidgets.QPushButton("save")
    save_button.setToolTip(TH.tip("Store the current map with its full recipe, so it can be "
                                  "rebuilt exactly."))
    save_button.clicked.connect(panel.save_embedding)
    sl.addWidget(panel.emb_name, 1)
    sl.addWidget(save_button)
    shell.add_setting("Save", save)

    shell.preview = LivePreview(panel, shell)
    app = QtWidgets.QApplication.instance()
    if app is not None:
        app.aboutToQuit.connect(shell.preview.stop)
    shell.add_section("Live preview", shell.preview, extent=440, stretch=3)

    search = QtWidgets.QWidget(shell)
    sv = QtWidgets.QVBoxLayout(search)
    sv.setContentsMargins(8, 4, 8, 8)
    grid = _moved(search, [("sample size", panel.sample),
                           ("n_neighbors values", panel.nn_grid),
                           ("min_dist values", panel.md_grid),
                           ("min_cluster_size values", panel.mcs_grid)])
    # The splitter shrinks a section to its minimum, and a form's minimum lets its rows
    # overlap; the grid keeps its height and the table gives way instead.
    grid.setSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred, QtWidgets.QSizePolicy.Policy.Fixed)
    sv.addWidget(grid)
    panel.walk_table.setMinimumHeight(90)
    sv.addWidget(panel.walk_table, 1)
    hint = QtWidgets.QLabel("Click a row to build that map and show it in the 3D view. "
                            "Right-click to save the table as CSV.")
    hint.setObjectName("ShellMuted")
    hint.setWordWrap(True)
    sv.addWidget(hint)
    search.setMinimumHeight(grid.sizeHint().height() + 90 + 40)
    shell.add_section("Hyperparameter search", search, extent=290, stretch=1)

    shell.add_action("build this map", panel.run_embed, primary=True,
                     tip="Project every gene with these settings and show the map in the 3D view.")
    shell.add_action("walk hyperparameters", panel.run_umap_walk,
                     tip="Build a map for every combination in Hyperparameter search and score "
                         "each; the table fills a row at a time.")
    shell.add_action_widget(panel.stop_button())
    shell.add_action("show the 3D map", lambda: _show_map(window),
                     tip="Leave the module and look at the map itself.")
    panel.status.connect(shell.say)

    shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                         ("AI chat", getattr(window, "chat_dock", None))])
    shell.home_requested.connect(lambda: _go_home(window))

    _replace_tab(panel, window)
    return shell


def _replace_tab(panel, window) -> None:
    """Put a way to the module where the Map tab's controls were."""
    tabs = panel.findChild(QtWidgets.QTabWidget)
    if tabs is None:
        return
    for i in range(tabs.count()):
        if tabs.tabText(i) == MAP_TAB:
            page = QtWidgets.QWidget()
            v = QtWidgets.QVBoxLayout(page)
            note = QtWidgets.QLabel("Building a map -- its settings, the hyperparameter walk and a "
                                    "live preview -- is now the Map module.")
            note.setWordWrap(True)
            go = QtWidgets.QPushButton("Open the Map module")
            go.setProperty("primary", True)
            go.setToolTip("Open the Map module: projection settings, a live preview, the "
                          "hyperparameter walk and building the map.")
            go.clicked.connect(lambda: open_in(window))
            v.addWidget(note)
            v.addWidget(go, 0, QtCore.Qt.AlignmentFlag.AlignLeft)
            v.addStretch(1)
            scroll = tabs.widget(i)
            if isinstance(scroll, QtWidgets.QScrollArea):
                scroll.setWidget(page)              # deletes the emptied old page
            return


def open_in(window) -> str:
    """Show the Map module over the window. Returns what the status bar was told."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "map_module", None)
    if landing is None or shell is None:
        return "the Map module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    return "Map"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)


def _show_map(window) -> None:
    from .landing_install import leave_landing
    leave_landing(window)


def _open_data_tab(window) -> None:
    from .landing_install import leave_landing, raise_dock
    leave_landing(window)
    raise_dock(window, "analysis", "1 · Data")
