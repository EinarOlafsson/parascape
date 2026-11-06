"""The Map module: where a gene or a list sits among all the others, and what that suggests.

Every gene is placed so that genes with similar measurements sit close together. The
module is built around what that is FOR (instruction 96, agreed with Einar 2026-11-06),
not around the method that draws it:

* **Synthesis** -- put your gene or list on the map and read its neighbourhood: "8 of
  10 genes around GRA16 are dense granule proteins".
* **Inference** -- a gene of unknown function among genes of known function is a
  candidate for it, and the module says how often that kind of call is right on this
  same map (`map_reading.agreement`), against chance.
* **Choosing the map** -- which measurements make genes "similar" is the one choice
  that changes what the map means, so it is on the left in plain words; the UMAP
  settings are under Advanced, and "Find the best map for a label" tries them for you.

The first version of this module was the analysis panel's "2 · Map" tab moved over --
n_neighbors, min_dist, seed and a grid of values to sweep -- and Einar found it hard to
understand: it showed the machinery and hid the map. The panel's own widgets are still
the ones used (MOVED here, not copied), so `AnalysisPanel.spec()`, `run_embed` and
`run_umap_walk` read exactly what they read before.

The map is drawn in software (`ScatterView`), never as a second GL view: the window's
GL map must not share or lose a context, and a pyqtgraph PlotWidget deleted mid-session
crashed PySide6 inside `QGraphicsScene`'s destructor.
"""
from __future__ import annotations

import os
from html import escape

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from . import map_reading as MR
from . import theme as TH
from .shell import ModuleShell

#: Genes in a preview of new settings: enough for the structure to show, few enough that
#: a projection takes a second or two once umap has compiled.
PREVIEW_GENES = 1500
#: How long the settings must rest before the preview is recomputed, in ms.
PREVIEW_DEBOUNCE_MS = 600
#: Where the Map tab used to be.
MAP_TAB = "2 · Map"
#: Neighbours listed by name under a single gene.
LISTED = 10


#: Projections still running. A QThread destroyed while it runs aborts the process, and a
#: preview can be closed (or the window quit) in the middle of one, so a projection is
#: owned here, not by its widget, until it has finished.
_RUNNING: set = set()


def wait_for_projections(timeout_ms: int = 60000) -> None:
    """Block until every running projection has finished (on quit, and in tests)."""
    for worker in list(_RUNNING):
        worker.wait(timeout_ms)


class _Projection(QtCore.QThread):
    """One 3D projection of a sample of genes, off the GUI thread."""

    ready = QtCore.Signal(int, object, object)        # generation, xyz, node-table positions
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
            import umap
            reducer = umap.UMAP(n_components=3,
                                n_neighbors=int(min(self.spec.n_neighbors, len(pick) - 1)),
                                min_dist=float(self.spec.min_dist),
                                random_state=int(self.spec.random_state))
            xyz = reducer.fit_transform(X[pick])
            kept = np.asarray(rows)
            kept = np.flatnonzero(kept) if kept.dtype == bool else kept.astype(int)
            self.ready.emit(self.generation, xyz, kept[pick])
        except Exception as error:                    # a preview must never take the module
            self.failed.emit(self.generation, f"{type(error).__name__}: {error}")


def _rgb255(colour):
    """A colour as 0-255 RGBA, whichever form it arrives in."""
    if isinstance(colour, QtGui.QColor):
        return colour.red(), colour.green(), colour.blue(), 220
    if isinstance(colour, str):
        q = QtGui.QColor(colour)
        return q.red(), q.green(), q.blue(), 220
    values = list(colour)
    if values and max(values[:3]) <= 1.0:
        values = [int(round(v * 255)) for v in values]
    return (*values[:3], 220)


class ScatterView(QtWidgets.QWidget):
    """Genes in three dimensions, drawn in software: drag to turn, wheel to zoom, click
    to pick. Nearer points are larger and brighter, so depth reads without a GL depth
    buffer. Highlighted genes are drawn last, larger, with a white rim."""

    picked = QtCore.Signal(int)                        # a node-table position
    FPS = 20
    SPIN = 0.010

    def __init__(self, parent=None):
        super().__init__(parent)
        self.xyz = np.zeros((0, 3))
        self.ids = np.zeros(0, int)
        self.colours: list = []
        self.highlight: set = set()
        self.yaw, self.pitch, self.zoom = 0.6, 0.35, 1.0
        self.spinning = False
        self._drag = None
        self._moved = False
        self.setMinimumHeight(220)
        self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)
        self.setMouseTracking(False)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(int(1000 / self.FPS))
        self.timer.timeout.connect(self._turn)

    def set_points(self, xyz, colours, ids=None) -> None:
        """Show `xyz` (n x 3), one colour per point; `ids` are their node-table positions."""
        xyz = np.asarray(xyz, dtype=float)
        if xyz.ndim == 2 and xyz.shape[1] == 2:
            xyz = np.column_stack([xyz, np.zeros(len(xyz))])
        xyz = xyz.reshape(-1, 3)
        if len(xyz):
            # Centred on the median and scaled to the 80th-percentile radius: scaled to the
            # farthest point (or even the 98th percentile), the small outlying clusters a
            # UMAP throws far out shrank the body of the map to a dot.
            xyz = xyz - np.median(xyz, axis=0)
            radius = float(np.percentile(np.linalg.norm(xyz, axis=1), 80))
            xyz = xyz / max(1e-9, radius)
        self.xyz = xyz
        self.ids = np.arange(len(xyz)) if ids is None else np.asarray(ids, int)
        self.colours = [QtGui.QColor(*_rgb255(c)) for c in colours]
        self.update()

    def set_highlight(self, ids) -> None:
        """Mark these node-table positions."""
        self.highlight = {int(i) for i in ids}
        self.update()

    def set_spinning(self, on: bool) -> None:
        """Turn the map slowly on its own, or stop."""
        self.spinning = bool(on)
        if on and self.isVisible():
            self.timer.start()
        else:
            self.timer.stop()

    def showEvent(self, event):
        super().showEvent(event)
        if self.spinning:
            self.timer.start()

    def hideEvent(self, event):
        self.timer.stop()
        super().hideEvent(event)

    def _turn(self):
        if self._drag is None:
            self.yaw += self.SPIN
            self.update()

    def mousePressEvent(self, event):
        self._drag = event.position()
        self._moved = False
        self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)

    def mouseMoveEvent(self, event):
        if self._drag is None:
            return
        delta = event.position() - self._drag
        if abs(delta.x()) + abs(delta.y()) > 1:
            self._moved = True
        self._drag = event.position()
        self.yaw += delta.x() * 0.01
        self.pitch = float(np.clip(self.pitch + delta.y() * 0.01, -1.4, 1.4))
        self.update()

    def mouseReleaseEvent(self, event):
        if not self._moved:
            hit = self.point_at(event.position())
            if hit is not None:
                self.picked.emit(hit)
        self._drag = None
        self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)

    def wheelEvent(self, event):
        steps = event.angleDelta().y() / 120.0
        self.zoom = float(np.clip(self.zoom * (1.15 ** steps), 0.3, 12.0))
        self.update()
        event.accept()

    def projected(self):
        """Screen x, y, a 0-1 nearness per point, and the back-to-front drawing order."""
        cy, sy, cp, sp = np.cos(self.yaw), np.sin(self.yaw), np.cos(self.pitch), np.sin(self.pitch)
        x, y, z = self.xyz[:, 0], self.xyz[:, 1], self.xyz[:, 2]
        x1, z1 = cy * x + sy * z, -sy * x + cy * z
        y1, z2 = cp * y - sp * z1, sp * y + cp * z1
        camera = 3.2
        scale = camera / np.maximum(0.2, camera - z2)
        side = 0.40 * min(self.width(), 1.2 * self.height()) * self.zoom
        sx = self.width() / 2 + x1 * scale * side
        sy_ = self.height() / 2 - y1 * scale * side
        near = np.clip((z2 + 1) / 2, 0, 1)
        return sx, sy_, near, np.argsort(z2)

    def point_at(self, pos, radius: float = 8.0):
        """The node-table position of the nearest point to `pos` within `radius` px."""
        if not len(self.xyz):
            return None
        sx, sy, _near, _order = self.projected()
        d = (sx - pos.x()) ** 2 + (sy - pos.y()) ** 2
        i = int(np.argmin(d))
        return int(self.ids[i]) if d[i] <= radius ** 2 else None

    def paintEvent(self, _event):
        if not len(self.xyz):
            return
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        sx, sy, near, order = self.projected()
        marked = []
        for i in order:
            if int(self.ids[i]) in self.highlight:
                marked.append(i)
                continue
            colour = QtGui.QColor(self.colours[i] if i < len(self.colours)
                                  else QtGui.QColor("white"))
            colour.setAlphaF(0.30 + 0.6 * float(near[i]))
            painter.setBrush(colour)
            r = 1.1 + 2.2 * float(near[i])
            painter.drawEllipse(QtCore.QPointF(float(sx[i]), float(sy[i])), r, r)
        rim = QtGui.QPen(QtGui.QColor(255, 255, 255, 235), 1.6)
        for i in marked:
            painter.setPen(rim)
            colour = QtGui.QColor(self.colours[i] if i < len(self.colours)
                                  else QtGui.QColor("white"))
            painter.setBrush(colour)
            painter.drawEllipse(QtCore.QPointF(float(sx[i]), float(sy[i])), 6.0, 6.0)
        painter.end()


def _moved(parent: QtWidgets.QWidget, rows) -> QtWidgets.QWidget:
    """A form of (label, widget) rows, the widgets reparented into it."""
    page = QtWidgets.QWidget(parent)
    form = QtWidgets.QFormLayout(page)
    form.setContentsMargins(6, 2, 6, 6)
    for label, widget in rows:
        form.addRow(label, widget)
    return page


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


class MapModule(QtCore.QObject):
    """The Map module for one window: its shell, its state, and what each control does."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.panel = window.panel
        self.selected: list = []                       # node-table positions
        self.selection_note = ""
        self.mode = "map"                              # "map", or "preview" of new settings
        self.generation = 0
        self._worker = None
        self._pending = False
        self._agreement_key = None
        self._agreement = None

        self.shell = shell = ModuleShell(
            "map", "Map", "every gene, placed by how alike its measurements are",
            instruction="Put your genes on the map and read what surrounds them.")
        shell.module = self
        # A list dropped on the module is read here, onto the map, rather than going to the
        # gene-list panel behind the module as a drop anywhere else on the window does.
        shell.setAcceptDrops(True)
        shell.installEventFilter(self)

        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))

        self.preview_timer = QtCore.QTimer(self)
        self.preview_timer.setSingleShot(True)
        self.preview_timer.setInterval(PREVIEW_DEBOUNCE_MS)
        self.preview_timer.timeout.connect(self.preview)
        for box in (self.panel.nn, self.panel.md, self.panel.seed):
            box.valueChanged.connect(self._settings_changed)
        self.panel.blocks.itemChanged.connect(lambda *_: self._sync_similarity())
        self.panel.embedding_ready.connect(lambda *_: self.show_current())
        self.panel.status.connect(shell.say)
        self.panel.walk_started.connect(self._walk_started)
        self.panel.walk_step.connect(self._walk_step)
        window.category_box.currentTextChanged.connect(self._colour_changed_elsewhere)
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self.stop)
        _replace_tab(self.panel, window)
        self.show_current()

    # ------------------------------------------------------------------ the left column
    def _build_settings(self) -> None:
        shell, panel, window = self.shell, self.panel, self.window

        mine = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(mine)
        v.setContentsMargins(6, 2, 6, 6)
        self.gene_field = QtWidgets.QLineEdit()
        self.gene_field.setPlaceholderText("A gene id or product, or several, then Return")
        self.gene_field.setToolTip(TH.tip(
            "Type a gene id (TGME49_208830), a name in its product (GRA16), or several "
            "separated by commas or spaces. Return puts them on the map."))
        self.gene_field.returnPressed.connect(self._genes_typed)
        v.addWidget(self.gene_field)
        row = QtWidgets.QHBoxLayout()
        choose = QtWidgets.QPushButton("Choose a list…")
        choose.setToolTip(TH.tip("A CSV, TSV, XLS or XLSX with gene ids in a column. You can "
                                 "also drop the file anywhere on this module."))
        choose.clicked.connect(self._choose_list)
        use = QtWidgets.QPushButton("Use the gene list")
        use.setToolTip(TH.tip("The list last dropped on ParaScape (the gene list panel)."))
        use.clicked.connect(self._use_gene_list)
        clear = QtWidgets.QPushButton("Clear")
        clear.setToolTip("Take your genes off the map.")
        clear.clicked.connect(lambda: self.select([], ""))
        for b in (choose, use, clear):
            row.addWidget(b)
        v.addLayout(row)
        self.selection_label = _muted("Nothing selected. Click a point on the map, too.")
        v.addWidget(self.selection_label)
        shell.add_setting("Your genes", mine)

        colour = QtWidgets.QWidget(shell)
        cl = QtWidgets.QVBoxLayout(colour)
        cl.setContentsMargins(6, 2, 6, 6)
        self.colour_box = QtWidgets.QComboBox()
        self.colour_box.setToolTip(TH.tip(
            "What the colours mean, and which label your genes' neighbours are read by. The "
            "same choice as 'color by' in the find panel."))
        self.colour_box.addItems(window.color_sources())
        self.colour_box.setCurrentText(window.category_box.currentText())
        self.colour_box.currentTextChanged.connect(self._colour_chosen)
        cl.addWidget(self.colour_box)
        shell.add_setting("Colour by", colour)

        similar = QtWidgets.QWidget(shell)
        sl = QtWidgets.QVBoxLayout(similar)
        sl.setContentsMargins(6, 2, 6, 6)
        sl.addWidget(_muted("Genes are placed near each other when these are alike:"))
        self.similarity: dict = {}
        tree = panel.blocks
        for i in range(tree.topLevelItemCount()):
            item = tree.topLevelItem(i)
            usable = [leaf for leaf in tree.leaves(item)
                      if leaf.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled]
            if not usable:
                continue
            box = QtWidgets.QCheckBox(item.text(0).split("  (")[0])
            box.setTristate(False)
            box.setToolTip(TH.tip(
                f"{len(usable)} kinds of measurement about {item.text(0).split('  (')[0]}. "
                "Ticked, they count towards which genes are placed together. The exact "
                "measurements are in Advanced ▸ Choose individual measurements."))
            box.clicked.connect(lambda on, it=item: self._similarity_clicked(it, on))
            sl.addWidget(box)
            self.similarity[box] = item
        self.circular_note = _muted("")
        self.circular_note.hide()
        sl.addWidget(self.circular_note)
        shell.add_setting("What makes genes similar", similar)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        al.addWidget(_muted("How the map is drawn. The defaults suit most uses; Find the best "
                            "map for a label tries these for you."))
        al.addWidget(_moved(advanced, [("neighbours considered (n_neighbors)", panel.nn),
                                       ("spacing (min_dist)", panel.md), ("seed", panel.seed)]))
        self.live = QtWidgets.QCheckBox("Preview a change before building it")
        self.live.setChecked(True)
        self.live.setToolTip(TH.tip("When a setting changes, project a sample of genes with "
                                    "it and show that, so you see the effect before building."))
        al.addWidget(self.live)
        al.addWidget(_muted("Values tried by Find the best map:"))
        al.addWidget(_moved(advanced, [("genes per try", panel.sample),
                                       ("neighbours considered", panel.nn_grid),
                                       ("spacing", panel.md_grid),
                                       ("cluster sizes", panel.mcs_grid)]))
        tree_button = QtWidgets.QPushButton("Choose individual measurements…")
        tree_button.setToolTip("The full measurement tree, block by block (analysis ▸ 1 · Data).")
        tree_button.clicked.connect(lambda: _open_data_tab(window))
        al.addWidget(tree_button, 0, QtCore.Qt.AlignmentFlag.AlignLeft)
        shell.add_setting("Advanced", advanced, expanded=False)

        save = QtWidgets.QWidget(shell)
        sv = QtWidgets.QHBoxLayout(save)
        sv.setContentsMargins(6, 2, 6, 6)
        save_button = QtWidgets.QPushButton("Save")
        save_button.setToolTip(TH.tip("Keep this map with its full recipe, so it can be "
                                      "rebuilt exactly."))
        save_button.clicked.connect(panel.save_embedding)
        sv.addWidget(panel.emb_name, 1)
        sv.addWidget(save_button)
        shell.add_setting("Save this map", save, expanded=False)
        self._sync_similarity()

    # ------------------------------------------------------------------ the right column
    def _build_runtime(self) -> None:
        shell, panel, window = self.shell, self.panel, self.window

        mapbox = QtWidgets.QWidget(shell)
        ml = QtWidgets.QVBoxLayout(mapbox)
        ml.setContentsMargins(8, 4, 8, 8)
        top = QtWidgets.QHBoxLayout()
        self.map_note = _muted("")
        self.back_to_map = QtWidgets.QPushButton("Back to the current map")
        self.back_to_map.setToolTip("Stop previewing and show the map as it is built now.")
        self.back_to_map.clicked.connect(self.show_current)
        self.back_to_map.hide()
        self.turn = QtWidgets.QCheckBox("Turn")
        self.turn.setToolTip("Turn the map slowly. Drag to turn it yourself, wheel to zoom.")
        self.turn.toggled.connect(lambda on: self.view.set_spinning(on))
        top.addWidget(self.map_note, 1)
        top.addWidget(self.back_to_map)
        top.addWidget(self.turn)
        ml.addLayout(top)
        self.view = ScatterView(mapbox)
        self.view.setToolTip(TH.tip("Drag to turn, wheel to zoom, click a gene to read its "
                                    "neighbours."))
        self.view.picked.connect(lambda i: self.select([i], ""))
        ml.addWidget(self.view, 1)
        self.legend = QtWidgets.QLabel("")
        self.legend.setObjectName("ShellMuted")
        self.legend.setWordWrap(True)
        self.legend.setTextFormat(QtCore.Qt.TextFormat.RichText)
        ml.addWidget(self.legend)
        shell.add_section("Map", mapbox, extent=470, stretch=3)

        self.reading = QtWidgets.QTextBrowser(shell)
        self.reading.setOpenLinks(False)
        self.reading.anchorClicked.connect(self._link)
        shell.add_section("Your genes and their neighbours", self.reading, extent=210, stretch=2)

        best = QtWidgets.QWidget(shell)
        bl = QtWidgets.QVBoxLayout(best)
        bl.setContentsMargins(8, 4, 8, 8)
        self.best_note = _muted("")
        bl.addWidget(self.best_note)
        self.tried = QtWidgets.QTableWidget(0, 4)
        self.tried.setHorizontalHeaderLabels(["neighbours considered", "spacing",
                                              "neighbours agree", "by chance"])
        self.tried.horizontalHeader().setStretchLastSection(True)
        self.tried.verticalHeader().hide()
        self.tried.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tried.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.tried.setToolTip(TH.tip("Each map tried, best first. Click one to use its settings "
                                     "and build it."))
        self.tried.cellClicked.connect(self._tried_clicked)
        self.tried.setMinimumHeight(90)
        bl.addWidget(self.tried, 1)
        shell.add_section("Find the best map for a label", best, extent=200, stretch=1,
                          folded=True)

        shell.add_action("Build this map", panel.run_embed, primary=True,
                         tip="Place every gene with the settings on the left. Takes a minute.")
        shell.add_action("Predict from neighbours", lambda: _open_inference(window),
                         tip="Call labels for unlabelled genes from their neighbours, with how "
                             "often such calls are right (Infer).")
        shell.add_action("Find the best map", self._find_best,
                         tip="Try the Advanced settings in turn and score each by how well "
                             "neighbours agree on the label you are colouring by.")
        shell.add_action_widget(panel.stop_button())
        shell.add_action("Full 3D view", lambda: _show_map(window),
                         tip="Leave the module for the map in the main window, with every tool "
                             "it has.")
        shell.add_action("Export", lambda: _export(window),
                         tip="Save the genes shown, with their values, as a table.")
        shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                             ("AI chat", getattr(window, "chat_dock", None))],
                     folded=True)
        self._update_best_note()

    # ------------------------------------------------------------------ the map shown
    def labels(self):
        """The current colour-by label of every gene, as strings."""
        return self.window.category_values().astype(str).to_numpy()

    def show_current(self) -> None:
        """Show the map as built now, coloured as the window colours it."""
        w = self.window
        self.mode = "map"
        self.back_to_map.hide()
        xyz = getattr(w, "xyz", None)
        if xyz is None or not len(xyz):
            self.map_note.setText("No map yet: Build this map.")
            return
        vis = w.visible_mask()
        ids = np.flatnonzero(vis)
        colours = w.colors(vis)[ids]
        self.view.set_points(np.asarray(xyz)[ids], colours, ids)
        self.view.set_highlight(self.selected)
        self.map_note.setText(f"The current map · {len(ids):,} genes · coloured by "
                              f"{self.colour_box.currentText()}")
        self._legend()
        self._agreement = None
        self.read()

    def _legend(self) -> None:
        w = self.window
        values = w.category_values().astype(str)
        known = ~values.str.strip().str.lower().isin(MR.UNLABELLED)
        counts = values[known].value_counts()
        colour_of = getattr(w, "color_of", {}) or {}
        parts = []
        for label, n in counts.head(14).items():
            rgb = colour_of.get(label)
            if rgb is None:
                continue
            r, g, b = [int(round(c * 255)) if max(rgb[:3]) <= 1 else int(c) for c in rgb[:3]]
            parts.append(f"<span style='color: rgb({r},{g},{b})'>●</span> {escape(label)} "
                         f"<span style='color:#888'>{n}</span>")
        more = len(counts) - len(parts)
        self.legend.setText("&nbsp;&nbsp; ".join(parts)
                            + (f"&nbsp;&nbsp; <span style='color:#888'>+{more} more</span>"
                               if more > 0 else "")
                            + "&nbsp;&nbsp; <span style='color:#888'>grey: unknown</span>")

    def _colour_chosen(self, text: str) -> None:
        if self.window.category_box.currentText() != text:
            self.window.category_box.setCurrentText(text)
        self._update_best_note()
        self._sync_similarity()
        if self.mode == "map":
            self.show_current()

    def _colour_changed_elsewhere(self, text: str) -> None:
        if self.colour_box.currentText() != text:
            self.colour_box.blockSignals(True)
            self.colour_box.setCurrentText(text)
            self.colour_box.blockSignals(False)
            self._update_best_note()
            if self.mode == "map":
                self.show_current()

    # ------------------------------------------------------------------ your genes
    def select(self, ids, note: str) -> None:
        """Put these genes on the map and read their neighbourhood."""
        self.selected = [int(i) for i in ids]
        self.selection_note = note
        self.view.set_highlight(self.selected)
        if not self.selected:
            self.selection_label.setText("Nothing selected. Click a point on the map, too.")
        elif len(self.selected) == 1:
            row = self.window.nodes.iloc[self.selected[0]]
            self.selection_label.setText(f"{row.get('gene_id', '')} — "
                                         f"{str(row.get('product', ''))[:80]}")
        else:
            self.selection_label.setText(f"{len(self.selected)} genes" + (f" — {note}" if note else ""))
        self.read()

    def resolve(self, text: str) -> tuple:
        """Gene positions for typed text, and the tokens that matched nothing."""
        nodes = self.window.nodes
        ids = nodes["gene_id"].astype(str)
        lower = {g.lower(): i for i, g in enumerate(ids)}
        tokens = [t for t in text.replace(",", " ").replace(";", " ").split() if t]
        found, missing = [], []
        product = nodes["product"].astype(str).str.lower() if "product" in nodes else None
        for token in tokens:
            i = lower.get(token.lower())
            if i is None and product is not None:
                hit = np.flatnonzero(product.str.contains(token.lower(), regex=False).to_numpy())
                i = int(hit[0]) if len(hit) == 1 else None
                if len(hit) > 1:
                    found.extend(int(h) for h in hit[:50])
                    continue
            if i is None:
                missing.append(token)
            else:
                found.append(int(i))
        return list(dict.fromkeys(found)), missing

    def _genes_typed(self) -> None:
        found, missing = self.resolve(self.gene_field.text())
        note = f"{len(missing)} not found: {', '.join(missing[:5])}" if missing else ""
        self.select(found, note)

    def load_list(self, path: str) -> None:
        """Read a file of genes and put them on the map."""
        from . import gene_drop as GD
        known = self.window.nodes["gene_id"].astype(str).tolist()
        try:
            drop = GD.read_drop(path, known_ids=known)
        except Exception as error:
            self.shell.say(f"{os.path.basename(path)} could not be read: {error}")
            return
        index = {g: i for i, g in enumerate(known)}
        ids = [index[g] for g in drop.genes if g in index]
        self.select(ids, os.path.basename(path))

    def _choose_list(self) -> None:
        from . import gene_drop as GD
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.shell, "Choose a list of genes", "",
            "Tables (" + " ".join(f"*{e}" for e in GD.READABLE) + ")")
        if path:
            self.load_list(path)

    def _use_gene_list(self) -> None:
        panel = getattr(self.window, "home", None)
        genes = panel.genes() if panel is not None else []
        if not genes:
            self.shell.say("No gene list loaded: drop a file on this module or choose one.")
            return
        index = {g: i for i, g in enumerate(self.window.nodes["gene_id"].astype(str))}
        self.select([index[g] for g in genes if g in index], "the gene list")

    def eventFilter(self, obj, event):                    # noqa: N802 - Qt name
        """Accept a file dragged onto the module and read it as a list of genes."""
        T = QtCore.QEvent.Type
        if obj is self.shell and event.type() in (T.DragEnter, T.DragMove, T.Drop):
            from .home_panel import HomePanel
            paths = HomePanel.dropped_paths(event)
            if not paths:
                return False
            event.acceptProposedAction()
            if event.type() == T.Drop:
                self.load_list(paths[0])
            return True
        return False

    # ------------------------------------------------------------------ reading the map
    def agreement(self):
        """How often neighbours name a gene's label on the current map (cached)."""
        w = self.window
        key = (id(getattr(w, "xyz", None)), self.colour_box.currentText())
        if key != self._agreement_key or self._agreement is None:
            xyz = getattr(w, "xyz", None)
            self._agreement = (None if xyz is None
                               else MR.agreement(xyz, self.labels(), getattr(w, "placed", None)))
            self._agreement_key = key
        return self._agreement

    def read(self) -> None:
        """Write the neighbourhood of the selection, and how far to trust it."""
        w = self.window
        label_name = self.colour_box.currentText()
        if not self.selected:
            agree = self.agreement()
            text = ("<p>Put a gene or a list on the map -- type it, choose a file, drop one "
                    "here, or click a point -- to read what surrounds it.</p>")
            if agree is not None:
                text += (f"<p>On this map, a gene's 10 nearest neighbours name its own "
                         f"{escape(label_name)} {agree.accuracy:.0%} of the time "
                         f"({agree.n:,} genes with a known {escape(label_name)}; "
                         f"{agree.chance:.0%} by chance).</p>")
            self.reading.setHtml(text + self._circular_html())
            return
        labels = self.labels()
        hood = MR.neighbours(w.xyz, labels, self.selected, getattr(w, "placed", None))
        if not hood.members:
            self.reading.setHtml("<p>None of these genes has a place on the current map.</p>")
            return
        nodes = w.nodes
        html = []
        if len(hood.members) == 1:
            g = hood.members[0]
            own = labels[g].strip()
            gid = escape(str(nodes.gene_id.iloc[g]))
            html.append(f"<h3>{gid}</h3><p>{escape(str(nodes.iloc[g].get('product', '')))}</p>")
            html.append(self._summary(hood, label_name))
            if own and own.lower() not in MR.UNLABELLED:
                html.append(f"<p>Its own {escape(label_name)}: <b>{escape(own)}</b>.</p>")
            else:
                html.append(f"<p>Its own {escape(label_name)} is <b>not known</b>"
                            + (f" -- its neighbours suggest <b>{escape(hood.top[0])}</b>."
                               if hood.top else ".") + "</p>")
            html.append(self._trust(hood, label_name))
            rows = "".join(
                f"<tr><td><a href='gene:{j}'>{escape(str(nodes.gene_id.iloc[j]))}</a></td>"
                f"<td>{escape(str(nodes.iloc[j].get('product', ''))[:60])}</td>"
                f"<td>{escape(labels[j]) or '<i>unknown</i>'}</td></tr>"
                for j in hood.neighbours[:LISTED])
            html.append(f"<table cellspacing='4'><tr><th align='left'>neighbour</th>"
                        f"<th align='left'>product</th><th align='left'>{escape(label_name)}"
                        f"</th></tr>{rows}</table>")
        else:
            html.append(f"<h3>{len(hood.members)} genes"
                        + (f" — {escape(self.selection_note)}" if self.selection_note else "")
                        + "</h3>")
            own = [labels[g].strip() for g in hood.members]
            known = [x for x in own if x and x.lower() not in MR.UNLABELLED]
            if known:
                from collections import Counter
                top = Counter(known).most_common(3)
                html.append(f"<p>Their own {escape(label_name)}: " + ", ".join(
                    f"{escape(k)} ({n})" for k, n in top)
                    + f"; {len(own) - len(known)} unknown.</p>")
            html.append(self._summary(hood, label_name, around=True))
            html.append(self._trust(hood, label_name))
        self.reading.setHtml("".join(html) + self._circular_html())

    @staticmethod
    def _summary(hood, label_name: str, around: bool = False) -> str:
        if not hood.counts:
            return (f"<p>None of the {len(hood.neighbours)} genes around it has a known "
                    f"{escape(label_name)}.</p>")
        where = "Around them" if around else f"Its {len(hood.neighbours)} nearest genes"
        parts = ", ".join(f"<b>{n}</b> {escape(k)}" for k, n in hood.counts[:4])
        rest = sum(n for _, n in hood.counts[4:])
        return (f"<p>{where}: {parts}" + (f", {rest} other" if rest else "")
                + (f", {hood.unlabelled} unknown" if hood.unlabelled else "") + ".</p>")

    def _trust(self, hood, label_name: str) -> str:
        agree = self.agreement()
        if agree is None or not hood.top:
            return ""
        top = hood.top[0]
        call = agree.for_label(top)
        if call is None:
            return (f"<p style='color:#888'>On this map neighbours never call "
                    f"{escape(top)} on their own, so there is no record of how often that call "
                    f"is right.</p>")
        precision, n = call
        return (f"<p style='color:#888'>On this map, when a gene's neighbours mostly say "
                f"{escape(top)}, they are right <b>{precision:.0%}</b> of the time "
                f"({n:,} such calls among genes whose {escape(label_name)} is known). "
                f"Guessing the commonest label would be right {agree.chance:.0%}.</p>")

    def _circular_html(self) -> str:
        text = self.circular_note.text()
        return f"<p style='color:#d9a441'>{escape(text)}</p>" if text else ""

    def _link(self, url) -> None:
        text = url.toString()
        if text.startswith("gene:"):
            self.select([int(text.split(":", 1)[1])], "")

    # ------------------------------------------------------------------ what makes genes similar
    def _similarity_clicked(self, item, on: bool) -> None:
        tree = self.panel.blocks
        keys = [leaf.data(0, QtCore.Qt.ItemDataRole.UserRole) for leaf in tree.leaves(item)
                if leaf.flags() & QtCore.Qt.ItemFlag.ItemIsEnabled]
        current = set(tree.checked())
        current = (current | set(keys)) if on else (current - set(keys))
        tree.set_checked(sorted(current))
        self._sync_similarity()
        self._settings_changed()

    def _sync_similarity(self) -> None:
        """Tick each plain-language box when its measurements are in the map; warn when the
        colouring label is itself one of those measurements."""
        tree = self.panel.blocks
        checked = set(tree.checked())
        for box, item in self.similarity.items():
            keys = {leaf.data(0, QtCore.Qt.ItemDataRole.UserRole) for leaf in tree.leaves(item)}
            box.blockSignals(True)
            box.setChecked(bool(keys & checked))
            box.blockSignals(False)
        column = self.colour_box.currentText()
        used = set()
        try:
            from .embedding import columns_for
            for cols in columns_for(self.window.nodes, self.panel.spec()).values():
                used.update(cols)
        except Exception:
            pass
        if column in used:
            self.circular_note.setText(
                f"{column} is one of the measurements this map is built from, so neighbours "
                f"agreeing on it is partly built in -- not evidence.")
            self.circular_note.show()
        else:
            self.circular_note.setText("")
            self.circular_note.hide()

    # ------------------------------------------------------------------ previews
    def _settings_changed(self, *_):
        if self.live.isChecked() and self.shell.isVisible():
            self.preview_timer.start()

    def preview(self) -> None:
        """Project a sample with the settings on the left and show it in place of the map."""
        if not self.shell.isVisible():
            return
        if self._worker is not None:
            self._pending = True
            return
        try:
            spec = self.panel.spec()
        except Exception as error:
            self.map_note.setText(f"No preview: {error}")
            return
        self.generation += 1
        first = " The first one also compiles umap, up to a minute." if self.generation == 1 else ""
        self.map_note.setText(f"Previewing the new settings on {PREVIEW_GENES:,} genes…{first}")
        worker = _Projection(self.generation, self.panel.nodes, spec, PREVIEW_GENES)
        _RUNNING.add(worker)
        worker.finished.connect(lambda w=worker: _RUNNING.discard(w))
        worker.ready.connect(self._previewed)
        worker.failed.connect(lambda _g, message: self.map_note.setText(f"No preview: {message}"))
        worker.finished.connect(self._preview_done)
        self._worker = worker
        worker.start()

    def _previewed(self, generation: int, xyz, ids) -> None:
        if generation != self.generation:
            return
        colours = self.window.colors(np.ones(self.window.n, bool))[np.asarray(ids, int)]
        self.mode = "preview"
        self.view.set_points(xyz, colours, ids)
        self.view.set_highlight(self.selected)
        self.back_to_map.show()
        self.map_note.setText(f"Preview of the new settings · {len(ids):,} genes · "
                              "Build this map to apply them")

    def _preview_done(self) -> None:
        worker, self._worker = self._worker, None
        if worker is not None:
            worker.deleteLater()
        if self._pending:
            self._pending = False
            self.preview()

    def stop(self) -> None:
        """Stop asking for previews and wait for the one running, if any."""
        self.preview_timer.stop()
        self._pending = False
        wait_for_projections()

    # ------------------------------------------------------------------ find the best map
    def _update_best_note(self) -> None:
        self.best_note.setText(
            f"Tries the Advanced settings in turn on a sample of genes and scores each map by "
            f"how often neighbours agree on {self.colour_box.currentText()} -- the label you "
            f"are colouring by. Click a row to use those settings.")

    def _find_best(self) -> None:
        section = self.shell.sections.get("Find the best map for a label")
        if section is not None and section.folder.shut:
            section.folder.set_shut(False, by_user=False)
        self.panel.run_umap_walk()

    def _walk_started(self) -> None:
        self.tried.setRowCount(0)
        self._tried_rows = []

    def _walk_step(self, step) -> None:
        labels = self.labels()[np.asarray(step.genes, bool)]
        acc, chance = MR.label_agreement(step.coords, labels)
        row = dict(nn=int(step.row.get("n_neighbors", 0)), md=float(step.row.get("min_dist", 0)),
                   acc=acc, chance=chance)
        rows = getattr(self, "_tried_rows", [])
        rows.append(row)
        rows.sort(key=lambda r: -(r["acc"] if r["acc"] == r["acc"] else -1))
        self._tried_rows = rows
        self.tried.setRowCount(len(rows))
        for r, entry in enumerate(rows):
            cells = (str(entry["nn"]), f"{entry['md']:g}",
                     "" if entry["acc"] != entry["acc"] else f"{entry['acc']:.0%}",
                     "" if entry["chance"] != entry["chance"] else f"{entry['chance']:.0%}")
            for c, value in enumerate(cells):
                self.tried.setItem(r, c, QtWidgets.QTableWidgetItem(value))

    def _tried_clicked(self, row: int, _col: int) -> None:
        rows = getattr(self, "_tried_rows", [])
        if 0 <= row < len(rows):
            self.panel.nn.setValue(rows[row]["nn"])
            self.panel.md.setValue(rows[row]["md"])
            self.panel.run_embed()


def build(window) -> ModuleShell | None:
    """The Map module for `window`. None when the window has no analysis panel (it needs
    scikit-learn and umap-learn); the Map tile then opens the 3D view directly."""
    panel = getattr(window, "panel", None)
    if panel is None or not hasattr(panel, "nn"):
        return None
    return MapModule(window).shell


def _replace_tab(panel, window) -> None:
    """Put a way to the module where the Map tab's controls were."""
    tabs = panel.findChild(QtWidgets.QTabWidget)
    if tabs is None:
        return
    for i in range(tabs.count()):
        if tabs.tabText(i) == MAP_TAB:
            page = QtWidgets.QWidget()
            v = QtWidgets.QVBoxLayout(page)
            note = QtWidgets.QLabel("Building and reading the map is now the Map module.")
            note.setWordWrap(True)
            go = QtWidgets.QPushButton("Open the Map module")
            go.setProperty("primary", True)
            go.setToolTip("Open the Map module: put your genes on the map, read their "
                          "neighbours, and choose what the map is built from.")
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


def _open_inference(window) -> None:
    from .landing_install import leave_landing, raise_dock
    leave_landing(window)
    raise_dock(window, "analysis", "4 · Inference")


def _export(window) -> None:
    try:
        window.export_visible()
    except Exception as error:                        # an export must not take the module
        window.statusBar().showMessage(f"export failed: {error}", 8000)
