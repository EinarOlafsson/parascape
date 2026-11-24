"""The Network module: the whole network, your gene in focus, the links you choose.

Einar (instruction 87, restated 2026-11-06): "the network still looks nothing like I want. Now
you just get a gene with its connected genes and no context. I want to see the entire network
with the gene's region zoomed in, with the ability to zoom the network and change the edge
definitions of the network."

Built with the parascape-module mould. Every gene is drawn, placed as on the map (or arranged
by the links themselves); the links of the layers you tick, above the strength you set, are
drawn faintly everywhere; around your gene the network is sharp -- its links bright, the genes
one and two link steps away larger and named -- and fades with distance, the depth of field of
the "universe" plan. Wheel zooms about the pointer, drag pans, click a gene to travel to it.

Drawn in software (`NetworkView`), never as a second GL view or a pyqtgraph plot (see the
skill): the faint whole-network layer is rendered once into an image and re-rendered after
the view settles, so zooming and panning stay smooth over 300,000 links.
"""
from __future__ import annotations

import os
from html import escape

import numpy as np
import pandas as pd
from PySide6 import QtCore, QtGui, QtWidgets

from . import map_reading as MR
from . import network_reading as NR
from . import ports as PORTS
from . import theme as TH
from .recipes_panel import Setting, widget_setting
from .shell import ModuleShell, solid_dropdown
from . import organisms

#: A colour per measured layer, for links and the legend.
LAYER_COLOURS = {
    "xlms": "#ff6b6b", "ip_ms": "#ff9f43", "unwritten_interaction": "#feca57",
    "cofitness": "#48dbfb", "coexpression": "#1dd1a1", "cotranslation": "#10ac84",
    "compartment": "#a29bfe", "domain": "#fd79a8", "struct": "#e17055",
    "orthogroup": "#81ecec", "comention": "#dfe6e9", "comention_ft": "#b2bec3",
    "structural_hole": "#fab1a0",
}
#: Layers ticked when the module opens: the measured physical and functional ones. Shared
#: compartment is off -- it links every gene of a compartment to every other and would make
#: "where it is" look like "what it binds".
DEFAULT_LAYERS = ("xlms", "ip_ms", "cofitness", "coexpression", "cotranslation")
#: The most faint links drawn at once over the whole network, strongest first.
MAX_BACKGROUND_LINKS = 30000


class NetworkView(QtWidgets.QWidget):
    """Every gene and the chosen links, in software, with a camera: wheel zooms about the
    pointer, drag pans, a click picks the nearest gene or, between genes, the nearest bright
    link; shift+drag draws a box and selects every gene inside it.

    Speed (Einar 2026-11-24: "takes a while to load"): the faint links and the faint genes
    are built ONCE per choice of links / focus as shapes in data coordinates (`QPainterPath`
    per layer, `QPolygonF` per colour), and drawing them is one transformed call each with
    a cosmetic pen. Building 30,000 `QLineF` and 8,000 ellipses in Python on every redraw
    took 0.6-0.7 s. The genes are drawn live on every paint and never from the cached
    image, so they stay round and sharp while the image of the links is stretched during
    a zoom ("the nodes look pixelated").
    """

    picked = QtCore.Signal(int)
    #: Genes inside a shift+drag box, as positions in the gene table.
    selected = QtCore.Signal(list)
    #: A bright link clicked: the two genes it joins.
    link_picked = QtCore.Signal(int, int)

    #: Radius in pixels of a faint (far) gene.
    FAR_RADIUS = 1.7

    def __init__(self, parent=None):
        super().__init__(parent)
        self.xy = np.zeros((0, 2))
        self.colours: list = []
        self.links = pd.DataFrame(columns=["a", "b", "source", "strength"])
        self.dist = np.zeros(0)
        self.focus: list = []
        self.names: dict = {}
        self.show_names = True
        self.sharp = NR.SHARP_HOPS
        self.center = np.zeros(2)
        self.zoom = 1.0                                 # pixels per data unit / base
        self.chosen_link = None                         # (a, b) of a clicked link, drawn marked
        self._drag = None
        self._moved = False
        self._band = None                               # (origin, QRubberBand) while shift-dragging
        self._cache = None                              # (image, center, zoom, size)
        self._line_paths: dict = {}                     # layer -> QPainterPath in data coordinates
        self._line_key = None
        self._far_dots: list = []                       # [(QColor, QPolygonF)] in data coordinates
        self._animation = None
        self.redraw_timer = QtCore.QTimer(self)
        self.redraw_timer.setSingleShot(True)
        self.redraw_timer.setInterval(90)
        self.redraw_timer.timeout.connect(self._render_background)
        self.setMinimumHeight(260)
        self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)
        self.setMouseTracking(False)

    # ------------------------------------------------------------------ data
    def set_data(self, xy, colours, links, dist, focus, names=None, links_key=None) -> None:
        """Positions, gene colours, chosen links, link steps from the focus, and the focus.

        `links_key` names the choice the links came from (kinds, strength, arrangement); while
        it is unchanged the faint link shapes are reused, so moving the focus is quick.
        """
        self._links_key = links_key
        self.xy = np.asarray(xy, dtype=float)
        self.colours = [QtGui.QColor(*c) for c in colours]
        self.links = links
        self.dist = np.asarray(dist, dtype=float)
        self.focus = list(focus)
        self.names = names or {}
        self.chosen_link = None
        self._build_shapes()
        self._cache = None
        self._render_background()

    def _build_shapes(self) -> None:
        """The faint links (per layer) and the faint genes (per colour) as shapes in data
        coordinates. The links are rebuilt only when the links or the positions change."""
        links = self.links
        key = ((getattr(self, "_links_key", None),) if getattr(self, "_links_key", None) is not None
               else (id(links),)) + (len(links) if links is not None else 0, len(self.xy),
                                     float(self.xy[:, 0].sum()) if len(self.xy) else 0.0)
        if key != self._line_key:
            self._line_paths = {}
            if links is not None and len(links) and len(self.xy):
                if len(links) > MAX_BACKGROUND_LINKS:
                    links = links.nlargest(MAX_BACKGROUND_LINKS, "strength")
                a, b = links["a"].to_numpy(int), links["b"].to_numpy(int)
                sources = links["source"].astype(str).to_numpy()
                for layer in pd.unique(sources):
                    sel = sources == layer
                    path = QtGui.QPainterPath()
                    for (x1, y1), (x2, y2) in zip(self.xy[a[sel]].tolist(), self.xy[b[sel]].tolist()):
                        path.moveTo(x1, y1)
                        path.lineTo(x2, y2)
                    self._line_paths[layer] = path
            self._line_key = key
        far = (~(self.dist <= self.sharp) if len(self.dist) == len(self.xy)
               else np.ones(len(self.xy), bool))
        groups: dict = {}
        for i in np.flatnonzero(far).tolist():
            colour = self.colours[i] if i < len(self.colours) else QtGui.QColor("grey")
            groups.setdefault(colour.rgb(), (colour, []))[1].append(i)
        self._far_dots = []
        for colour, members in groups.values():
            faint = QtGui.QColor(colour)
            faint.setAlphaF(0.4)
            points = self.xy[members]
            self._far_dots.append((faint, QtGui.QPolygonF(
                [QtCore.QPointF(x, y) for x, y in points.tolist()])))

    # ------------------------------------------------------------------ camera
    def _base(self) -> float:
        return 0.45 * min(self.width(), self.height())

    def to_screen(self, xy) -> np.ndarray:
        """Data positions to widget pixels under the current camera."""
        s = self._base() * self.zoom
        xy = np.asarray(xy, dtype=float)
        return np.column_stack([self.width() / 2 + (xy[:, 0] - self.center[0]) * s,
                                self.height() / 2 - (xy[:, 1] - self.center[1]) * s])

    def transform(self) -> QtGui.QTransform:
        """The camera as a painter transform from data coordinates to widget pixels (the same
        mapping as `to_screen`)."""
        s = self._base() * self.zoom
        return QtGui.QTransform(s, 0.0, 0.0, -s, self.width() / 2 - self.center[0] * s,
                                self.height() / 2 + self.center[1] * s)

    def fit(self) -> None:
        """Show the whole network."""
        self.center = np.zeros(2)
        self.zoom = 1.0
        self._changed()

    def travel_to(self, index: int, zoom: float = 4.0) -> None:
        """Glide the camera to a gene and zoom into its region."""
        if not (0 <= index < len(self.xy)):
            return
        start_c, start_z = self.center.copy(), self.zoom
        end_c, end_z = self.xy[index].copy(), max(zoom, 1.0)
        if self._animation is not None:
            self._animation.stop()
        animation = QtCore.QVariantAnimation(self)
        animation.setDuration(450)
        animation.setStartValue(0.0)
        animation.setEndValue(1.0)
        animation.setEasingCurve(QtCore.QEasingCurve.Type.InOutCubic)

        def step(t):
            t = float(t)
            self.center = start_c + (end_c - start_c) * t
            self.zoom = start_z * (end_z / start_z) ** t
            self.update()
        animation.valueChanged.connect(step)
        animation.finished.connect(self._changed)
        self._animation = animation
        animation.start()

    def _changed(self) -> None:
        self.update()
        self.redraw_timer.start()

    def wheelEvent(self, event):                          # noqa: N802
        steps = event.angleDelta().y() / 120.0
        factor = 1.2 ** steps
        pos = event.position()
        s = self._base() * self.zoom
        at = np.array([self.center[0] + (pos.x() - self.width() / 2) / s,
                       self.center[1] - (pos.y() - self.height() / 2) / s])
        self.zoom = float(np.clip(self.zoom * factor, 0.3, 200.0))
        self.center = at + (self.center - at) / factor
        self._changed()
        event.accept()

    def mousePressEvent(self, event):                     # noqa: N802
        if event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier:
            origin = event.position().toPoint()
            band = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Shape.Rectangle, self)
            band.setGeometry(QtCore.QRect(origin, QtCore.QSize()))
            band.show()
            self._band = (origin, band)
            self.setCursor(QtCore.Qt.CursorShape.CrossCursor)
            return
        self._drag = event.position()
        self._moved = False
        self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)

    def mouseMoveEvent(self, event):                      # noqa: N802
        if self._band is not None:
            origin, band = self._band
            band.setGeometry(QtCore.QRect(origin, event.position().toPoint()).normalized())
            return
        if self._drag is None:
            return
        delta = event.position() - self._drag
        if abs(delta.x()) + abs(delta.y()) > 2:
            self._moved = True
        self._drag = event.position()
        s = self._base() * self.zoom
        self.center = self.center - np.array([delta.x() / s, -delta.y() / s])
        self.update()

    def mouseReleaseEvent(self, event):                   # noqa: N802
        if self._band is not None:
            origin, band = self._band
            rect = QtCore.QRect(origin, event.position().toPoint()).normalized()
            band.hide()
            band.deleteLater()
            self._band = None
            self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)
            genes = self.genes_in(rect)
            if genes:
                self.selected.emit(genes)
            return
        if self._moved:
            self.redraw_timer.start()
        else:
            hit = self.gene_at(event.position())
            if hit is not None:
                self.picked.emit(hit)
            else:
                link = self.link_at(event.position())
                if link is not None:
                    self.chosen_link = link
                    self.update()
                    self.link_picked.emit(*link)
        self._drag = None
        self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)

    def genes_in(self, rect) -> list:
        """The genes drawn inside a rectangle of the widget, as positions in the gene table."""
        if not len(self.xy) or rect.width() < 2 or rect.height() < 2:
            return []
        p = self.to_screen(self.xy)
        inside = ((p[:, 0] >= rect.left()) & (p[:, 0] <= rect.right())
                  & (p[:, 1] >= rect.top()) & (p[:, 1] <= rect.bottom()))
        return np.flatnonzero(inside).astype(int).tolist()

    def drawn_links(self) -> pd.DataFrame:
        """The links drawn brightly: those of the genes in focus, and those among their direct
        partners -- the ones a click can pick."""
        links = self.links
        if links is None or not len(links) or not len(self.dist) or not self.focus:
            return pd.DataFrame(columns=["a", "b", "source", "strength"])
        a, b = links["a"].to_numpy(int), links["b"].to_numpy(int)
        da, db = self.dist[a], self.dist[b]
        return links[((da == 0) | (db == 0)) | ((da == 1) & (db == 1))]

    def link_at(self, pos, radius: float = 5.0):
        """The bright link nearest a point on screen, within `radius` pixels, as (a, b)."""
        links = self.drawn_links()
        if not len(links):
            return None
        p = self.to_screen(self.xy)
        a, b = links["a"].to_numpy(int), links["b"].to_numpy(int)
        pa, pb = p[a], p[b]
        q = np.array([pos.x(), pos.y()])
        seg = pb - pa
        length2 = np.maximum((seg ** 2).sum(axis=1), 1e-9)
        t = np.clip(((q - pa) * seg).sum(axis=1) / length2, 0.0, 1.0)
        d2 = ((pa + seg * t[:, None] - q) ** 2).sum(axis=1)
        i = int(np.argmin(d2))
        return (int(a[i]), int(b[i])) if d2[i] <= radius ** 2 else None

    def resizeEvent(self, event):                         # noqa: N802
        super().resizeEvent(event)
        self.redraw_timer.start()

    def gene_at(self, pos, radius: float = 9.0):
        """The gene nearest a point on screen, within `radius` pixels."""
        if not len(self.xy):
            return None
        p = self.to_screen(self.xy)
        d = (p[:, 0] - pos.x()) ** 2 + (p[:, 1] - pos.y()) ** 2
        i = int(np.argmin(d))
        return i if d[i] <= radius ** 2 else None

    # ------------------------------------------------------------------ drawing
    def _render_background(self) -> None:
        """The faint links of the whole network into an image at the current camera; painting
        moves and scales it until the next render. Genes are not in it: they are drawn live."""
        if self.width() < 10 or self.height() < 10 or not len(self.xy):
            return
        ratio = self.devicePixelRatioF() or 1.0
        image = QtGui.QImage(int(self.width() * ratio), int(self.height() * ratio),
                             QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        image.setDevicePixelRatio(ratio)
        image.fill(QtCore.Qt.GlobalColor.transparent)
        painter = QtGui.QPainter(image)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self._draw_links(painter)
        painter.end()
        self._cache = (image, self.center.copy(), self.zoom, (self.width(), self.height()))
        self.update()

    def _draw_links(self, painter) -> None:
        """The strongest chosen links, faint, one transformed path per layer."""
        if not self._line_paths:
            return
        painter.save()
        painter.setTransform(self.transform(), True)
        for layer, path in self._line_paths.items():
            colour = QtGui.QColor(LAYER_COLOURS.get(str(layer), "#888888"))
            colour.setAlphaF(0.045)
            pen = QtGui.QPen(colour, 0.6)
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            painter.drawPath(path)
        painter.restore()

    def _draw_far_genes(self, painter) -> None:
        """Every gene off the focus, faint: round points of a fixed size in pixels, one
        transformed call per colour."""
        if not self._far_dots:
            return
        # The points are mapped through the camera in C++ and drawn untransformed: a wide
        # round point drawn UNDER a painter transform comes out as a short horizontal dash.
        camera = self.transform()
        painter.save()
        for colour, points in self._far_dots:
            pen = QtGui.QPen(colour, 2 * self.FAR_RADIUS)
            pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.drawPoints(camera.map(points))
        painter.restore()

    def _draw_background(self, painter) -> None:
        """The whole network, faint: the strongest chosen links and every gene off the focus."""
        self._draw_links(painter)
        self._draw_far_genes(painter)

    def paint_figure(self, painter) -> None:
        """The view as a figure: background and focus drawn as vectors, never the cached image."""
        if not len(self.xy):
            return
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self._draw_background(painter)
        self._paint_focus(painter)

    def paintEvent(self, _event):                         # noqa: N802
        if not len(self.xy):
            return
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        if self._cache is not None:
            image, c0, z0, size0 = self._cache
            # The cached links, moved and scaled to the current camera until redrawn.
            k = self.zoom / z0
            s = self._base() * self.zoom
            dx = (c0[0] - self.center[0]) * s + self.width() / 2 * (1 - k)
            dy = -(c0[1] - self.center[1]) * s + self.height() / 2 * (1 - k)
            painter.save()
            painter.setRenderHint(QtGui.QPainter.RenderHint.SmoothPixmapTransform, True)
            painter.translate(dx, dy)
            painter.scale(k, k)
            painter.drawImage(QtCore.QPointF(0, 0), image)
            painter.restore()
        self._draw_far_genes(painter)
        self._paint_focus(painter)
        painter.end()

    def _paint_focus(self, painter) -> None:
        """The sharp part around the focus: its own links bright, links among its direct
        partners faint, genes up to `sharp` steps away as dots that grow the nearer they are,
        and names for the focus and its strongest partners only -- two full link steps of
        a well-studied gene are thousands of genes, and drawing all their links is a hairball."""
        if not len(self.dist) or not self.focus:
            return
        p = self.to_screen(self.xy)
        links = self.links
        if links is not None and len(links):
            a, b = links["a"].to_numpy(int), links["b"].to_numpy(int)
            da, db = self.dist[a], self.dist[b]
            own = (da == 0) | (db == 0)
            among = (da == 1) & (db == 1)
            if among.sum() > 2500:                       # keep the strongest of a dense core
                idx = np.flatnonzero(among)
                strongest = idx[np.argsort(-links["strength"].to_numpy()[idx])[:2500]]
                among = np.zeros(len(a), bool)
                among[strongest] = True
            for mask, alpha, width in ((among, 0.16, 0.8), (own, 0.9, 1.8)):
                sub = links[mask]
                for layer in pd.unique(sub["source"]):
                    part = sub[sub["source"] == layer]
                    colour = QtGui.QColor(LAYER_COLOURS.get(str(layer), "#cccccc"))
                    colour.setAlphaF(alpha)
                    painter.setPen(QtGui.QPen(colour, width))
                    pa = p[part["a"].to_numpy(int)]
                    pb = p[part["b"].to_numpy(int)]
                    painter.drawLines([QtCore.QLineF(x1, y1, x2, y2)
                                       for (x1, y1), (x2, y2) in zip(pa, pb)])
        if self.chosen_link is not None:
            a0, b0 = self.chosen_link
            if 0 <= a0 < len(p) and 0 <= b0 < len(p):
                painter.setPen(QtGui.QPen(QtGui.QColor(255, 255, 255, 240), 3.2))
                painter.drawLine(QtCore.QPointF(*p[a0]), QtCore.QPointF(*p[b0]))
        near = np.flatnonzero(self.dist <= self.sharp)
        for i in sorted(near.tolist(), key=lambda j: -self.dist[j]):
            d = self.dist[i]
            colour = QtGui.QColor(self.colours[i] if i < len(self.colours) else QtGui.QColor("white"))
            colour.setAlphaF(1.0 if d == 0 else 0.9 if d == 1 else 0.45)
            r = 8.0 if d == 0 else 4.5 if d == 1 else 2.4
            painter.setPen(QtGui.QPen(QtGui.QColor(255, 255, 255, 230), 2.0) if d == 0
                           else QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(colour)
            painter.drawEllipse(QtCore.QPointF(*p[i]), r, r)
        if self.show_names:
            for i, name in self.names.items():
                if 0 <= i < len(p):
                    painter.setPen(QtGui.QColor(240, 240, 240, 235 if self.dist[i] == 0 else 185))
                    painter.drawText(QtCore.QPointF(p[i, 0] + 9, p[i, 1] + 4), name)


#: Layouts still being computed, owned here so a closed module does not destroy one running.
_RUNNING: set = set()


class _Arrange(QtCore.QThread):
    done = QtCore.Signal(object, object)                 # key, xy

    def __init__(self, key, adj):
        super().__init__(None)
        self.key, self.adj = key, adj

    def run(self):
        try:
            self.done.emit(self.key, NR.spectral_layout(self.adj))
        except Exception:
            self.done.emit(self.key, None)


def wait_for_layouts(timeout_ms: int = 120000) -> None:
    """Block until every layout being computed has finished (on quit, and in tests)."""
    for worker in list(_RUNNING):
        worker.wait(timeout_ms)


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _rgb255(c) -> tuple:
    values = list(c)
    if values and max(values[:3]) <= 1.0:
        values = [int(round(v * 255)) for v in values[:3]]
    return (*values[:3], 255)


class NetworkModule(QtCore.QObject):
    """The Network module for one window."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.edges = self._load_edges()
        self.layout_xy = None
        self.layout_kind = "map"
        self._layouts: dict = {}
        self._arranging: set = set()
        self.focus: list = []
        self._links = None
        self._adj = None
        self._agreement = None
        self.shell = shell = ModuleShell.for_module("network")
        shell.module = self
        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))
        window.category_box.currentTextChanged.connect(lambda *_: self.refresh())

    def _load_edges(self) -> pd.DataFrame:
        from . import organisms
        from . import star_edges as SE
        try:
            code = organisms.by_species(self.window.species).code
            return SE.measured_edges(SE.load_graph(code, self.window.nodes.gene_id.astype(str).tolist()))
        except Exception:
            return SE.empty()

    # ------------------------------------------------------------------ left
    def _build_settings(self) -> None:
        shell, window = self.shell, self.window
        mine = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(mine)
        v.setContentsMargins(6, 2, 6, 6)
        self.field = QtWidgets.QLineEdit()
        self.field.setPlaceholderText("A gene id or a word of its product, then Return")
        self.field.setToolTip(TH.tip("The gene to put in focus. Several separated by commas "
                                     "focus all of them -- a list's own network."))
        self.field.returnPressed.connect(self._typed)
        v.addWidget(self.field)
        self.focus_note = _muted("No gene in focus: the whole network is shown.")
        v.addWidget(self.focus_note)
        shell.add_setting("Your gene", mine, subject=PORTS.GENE)

        colour = QtWidgets.QWidget(shell)
        cl = QtWidgets.QVBoxLayout(colour)
        cl.setContentsMargins(6, 2, 6, 6)
        self.colour_box = QtWidgets.QComboBox()
        self.colour_box.addItems(window.color_sources())
        solid_dropdown(self.colour_box)
        self.colour_box.setCurrentText(window.category_box.currentText())
        self.colour_box.setToolTip(TH.tip("What the gene colours mean; the same choice as the "
                                          "map's. Linked genes are read by it."))
        self.colour_box.currentTextChanged.connect(
            lambda text: window.category_box.setCurrentText(text))
        cl.addWidget(self.colour_box)
        shell.add_setting("Colour by", colour)

        kinds = QtWidgets.QWidget(shell)
        kl = QtWidgets.QVBoxLayout(kinds)
        kl.setContentsMargins(6, 2, 6, 6)
        kl.addWidget(_muted("Two genes are linked when they share:"))
        from . import star_edges as SE
        self.layer_boxes: dict = {}
        counts = self.edges["source"].value_counts() if len(self.edges) else pd.Series(dtype=int)
        for layer, n in counts.items():
            box = QtWidgets.QCheckBox(f"{SE.layer_text(layer)}  ({n:,})")
            box.setChecked(layer in DEFAULT_LAYERS)
            box.setToolTip(TH.tip(f"{SE.layer_text(layer)}: {n:,} links."))
            box.setStyleSheet(f"QCheckBox {{ color: {LAYER_COLOURS.get(layer, '#ddd')}; }}")
            box.toggled.connect(lambda *_: self.refresh())
            kl.addWidget(box)
            self.layer_boxes[layer] = box
        kl.addWidget(_muted("How strong a link must be, within its kind:"))
        self.strength = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.strength.setRange(0, 95)
        self.strength.setValue(50)
        self.strength.setToolTip(TH.tip("Keep links at or above this percentile of their own "
                                        "kind. 50 keeps the stronger half of every kind."))
        self.strength_note = _muted("")
        self.strength.valueChanged.connect(lambda v: self.strength_note.setText(
            f"the strongest {100 - v}% of each kind"))
        self.strength.sliderReleased.connect(self.refresh)
        self.strength_note.setText("the strongest 50% of each kind")
        kl.addWidget(self.strength)
        kl.addWidget(self.strength_note)
        shell.add_setting("Which links count", kinds)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QFormLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        self.arrange = QtWidgets.QComboBox()
        self.arrange.addItems(["As on the map", "By the links"])
        solid_dropdown(self.arrange)
        self.arrange.setToolTip(TH.tip("Where genes sit: their places on the current map, or "
                                       "arranged so that linked genes sit together (spectral "
                                       "layout of the links you chose)."))
        self.arrange.currentTextChanged.connect(lambda *_: self._relayout())
        al.addRow("arrange genes", self.arrange)
        self.sharp = QtWidgets.QSpinBox()
        self.sharp.setRange(1, 4)
        self.sharp.setValue(NR.SHARP_HOPS)
        self.sharp.setToolTip(TH.tip("How many link steps around your gene are drawn sharply; "
                                     "the rest of the network fades."))
        self.sharp.valueChanged.connect(lambda *_: self.refresh())
        al.addRow("sharp link steps", self.sharp)
        self.names = QtWidgets.QCheckBox("name the genes linked to it")
        self.names.setChecked(True)
        self.names.setToolTip("Write gene names beside your gene and its direct links.")
        self.names.toggled.connect(lambda on: (setattr(self.view, "show_names", on), self.view.update()))
        al.addRow(self.names)
        shell.add_setting("Advanced", advanced, expanded=False)
        shell.add_recipes(self.recipe_settings(), after_apply=self.refresh)

    def recipe_settings(self) -> list:
        """The settings of the Network module that a recipe keeps (instruction 98, work package 6).

        The gene in focus is not kept: a recipe is how the module is set up, not which gene it is
        looking at. The kinds of link are kept by their layer name, one setting for each, so a
        layer that a data release drops is reported by name when an old recipe is applied.
        """
        out = [widget_setting("colour_by", self.colour_box, "colour by"),
               widget_setting("link_strength", self.strength, "how strong a link must be"),
               widget_setting("arrange_genes", self.arrange, "arrange genes"),
               widget_setting("sharp_link_steps", self.sharp, "sharp link steps"),
               widget_setting("name_linked_genes", self.names, "name the genes linked to it")]
        for layer, box in self.layer_boxes.items():
            out.append(widget_setting(f"links_{layer}", box, f"links: {layer}"))
        return out

    # ------------------------------------------------------------------ right
    def _build_runtime(self) -> None:
        shell, window = self.shell, self.window
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(wait_for_layouts)
        box = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(box)
        v.setContentsMargins(8, 4, 8, 8)
        self.note = _muted("")
        v.addWidget(self.note)
        self.view = NetworkView(box)
        self.view.setToolTip(TH.tip("Wheel to zoom, drag to move, click a gene to travel to it "
                                    "and read its links. Click a bright link to read how its two "
                                    "genes are linked. Hold Shift and drag a box to select "
                                    "every gene inside it."))
        self.view.picked.connect(lambda i: self.set_focus([i], travel=True))
        self.view.selected.connect(self._boxed)
        self.view.link_picked.connect(self.read_link)
        v.addWidget(self.view, 1)
        self.legend = QtWidgets.QLabel("")
        self.legend.setObjectName("ShellMuted")
        self.legend.setWordWrap(True)
        self.legend.setTextFormat(QtCore.Qt.TextFormat.RichText)
        v.addWidget(self.legend)
        shell.add_section("The network", box, extent=520, stretch=3)

        self.reading = QtWidgets.QTextBrowser(shell)
        shell.link_style(self.reading)
        self.reading.setOpenLinks(False)
        self.reading.anchorClicked.connect(self._link)
        shell.add_section("Your gene and its links", self.reading, extent=220, stretch=2)

        shell.add_action("Read it in Gene", self._to_gene, primary=True,
                         tip="Open the gene in focus in the Gene module: what is known, with sources.")
        shell.add_action("Put it on the map", self._to_map,
                         tip="Open the Map module with the gene in focus marked.")
        shell.add_action("Predict from links", self._to_infer,
                         tip="Call labels for unlabelled genes from their links, with how often such "
                             "calls are right (Infer). The genes in focus are filled in.")
        shell.add_action("Whole network", lambda: self.view.fit(),
                         tip="Zoom out to the whole network.")
        shell.add_action("Open a saved bundle…", lambda: self.open_bundle(),
                         tip="Read a bundle saved here back: its files are checked against their hashes "
                             "and its view is restored.")
        shell.add_action("Save everything as a bundle…", lambda: self.save_bundle(),
                         tip="One folder: the tables (each marked measured, inferred or computed), "
                             "the view as a vector figure, the methods text with its references, and "
                             "a manifest whose hashes show if anything was edited later.")
        shell.add_action("Save the figure…", lambda: self._save_figure(),
                         tip="The view as a vector figure (SVG or PDF) or a 3x PNG, drawn again for the "
                             "file rather than captured from the screen.")
        shell.add_action("Export its links", self._export,
                         tip="Save the links of the gene in focus as a table.")
        shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                             ("AI chat", getattr(window, "chat_dock", None))],
                     folded=True)

    # ------------------------------------------------------------------ state
    def layers(self) -> list:
        """The kinds of link ticked."""
        return [layer for layer, box in self.layer_boxes.items() if box.isChecked()]

    def labels(self) -> np.ndarray:
        """The colouring label of every gene, as text."""
        return self.window.category_values().astype(str).to_numpy()

    def _relayout(self) -> None:
        self.refresh()

    def _map_positions(self) -> np.ndarray:
        xyz = getattr(self.window, "xyz", None)
        return (NR.map_layout(xyz) if xyz is not None
                else np.random.default_rng(0).normal(size=(self.window.n, 2)))

    def _positions(self) -> np.ndarray:
        """Where genes sit. Arranged by the links, the layout is computed in the background
        (several seconds over the whole network) and kept per choice of links; until it
        arrives the map places are used."""
        if self.arrange.currentText() != "By the links":
            if self.layout_xy is None or self.layout_kind != "map":
                self.layout_xy, self.layout_kind = self._map_positions(), "map"
            return self.layout_xy
        key = (tuple(sorted(self.layers())), self.strength.value())
        cached = self._layouts.get(key)
        if cached is not None:
            self.layout_xy, self.layout_kind = cached, "links"
            return cached
        if key not in self._arranging:
            self._arranging.add(key)
            worker = _Arrange(key, self._adj)
            _RUNNING.add(worker)
            worker.done.connect(self._arranged)
            worker.finished.connect(lambda w=worker: (_RUNNING.discard(w), w.deleteLater()))
            worker.start()
            self.shell.say("Arranging the genes by their links… (several seconds)")
        if self.layout_xy is None:
            self.layout_xy, self.layout_kind = self._map_positions(), "map"
        return self.layout_xy

    def _arranged(self, key, xy) -> None:
        self._arranging.discard(key)
        if xy is None:
            self.shell.say("Could not arrange by the links; genes stay where the map puts them.")
            return
        self._layouts[key] = xy
        if self.arrange.currentText() == "By the links":
            self.shell.say("Arranged by the links.")
            self.refresh()
            if self.focus:
                self.view.travel_to(self.focus[0], zoom=5.0)

    def refresh(self) -> None:
        """Recompute the chosen links and redraw everything around the focus."""
        n = self.window.n
        self._links = NR.chosen(self.edges, self.layers(), self.strength.value() / 100.0)
        self._adj = NR.adjacency(self._links, n)
        self._agreement = None
        xy = self._positions()
        w = self.window
        colours = [_rgb255(c) for c in w.colors(np.ones(n, bool))]
        dist = NR.hops(self._adj, self.focus, max_hops=6) if self.focus else np.full(n, np.inf)
        ids = w.nodes["gene_id"].astype(str).to_numpy()
        names = {}
        for g in self.focus[:3]:
            names[int(g)] = ids[g]
            for partner in NR.links_of(self._links, int(g))["partner"].head(12):
                names[int(partner)] = ids[int(partner)]
        self.view.sharp = self.sharp.value()
        self.view.set_data(xy, colours, self._links, dist, self.focus, names,
                           links_key=(tuple(sorted(self.layers())), self.strength.value(),
                                      self.layout_kind))
        self.note.setText(f"{n:,} genes · {len(self._links):,} links of "
                          f"{len(self.layers())} kinds · coloured by "
                          f"{self.colour_box.currentText()}")
        self.legend.setText("&nbsp;&nbsp; ".join(
            f"<span style='color:{LAYER_COLOURS.get(layer, '#ddd')}'>━</span> {escape(layer)}"
            for layer in self.layers()))
        self.read()

    def set_focus(self, genes, travel: bool = True) -> None:
        """Put genes in focus: the network sharpens around them and the camera goes there."""
        self.focus = [int(g) for g in genes]
        ids = self.window.nodes["gene_id"].astype(str).to_numpy()
        self.focus_note.setText(", ".join(ids[g] for g in self.focus[:5])
                                + (f" and {len(self.focus) - 5} more" if len(self.focus) > 5 else "")
                                if self.focus else "No gene in focus: the whole network is shown.")
        self.refresh()
        if travel and self.focus:
            self.view.travel_to(self.focus[0], zoom=5.0)

    def focus_ids(self) -> list:
        """The gene identifiers in focus, in the order they were chosen."""
        ids = self.window.nodes["gene_id"].astype(str).to_numpy()
        return [str(ids[g]) for g in self.focus]

    def _boxed(self, genes: list) -> None:
        """Genes selected with a shift+drag box: all in focus, and named in the field so the
        selection can be read, edited and carried on."""
        self.set_focus(genes, travel=False)
        self.field.setText(", ".join(self.focus_ids()))
        self.shell.say(f"{len(genes):,} genes selected with the box.")

    def _to_infer(self) -> None:
        """Open Infer (predict from links) with the genes in focus already filled in."""
        from . import ports as PORTS
        from .landing_install import open_module
        genes = self.focus_ids()
        if not genes:
            open_module(self.window, "infer")
            return
        subject = PORTS.GENE if len(genes) == 1 else PORTS.GENE_LIST
        open_module(self.window, "infer", subject, ",".join(genes))

    def read_link(self, a: int, b: int) -> None:
        """How two genes are linked: every kind of link between them in the whole network --
        chosen or not -- with its strength within its kind, and what each gene is."""
        from . import star_edges as SE
        w = self.window
        ids = w.nodes["gene_id"].astype(str).to_numpy()
        products = w.nodes["product"].astype(str).to_numpy() if "product" in w.nodes else ids
        labels = self.labels()
        name = self.colour_box.currentText()
        edges = self.edges
        if len(edges):
            ea, eb = edges["a"].to_numpy(int), edges["b"].to_numpy(int)
            pair = edges[((ea == a) & (eb == b)) | ((ea == b) & (eb == a))]
        else:
            pair = edges
        chosen = set(self.layers())
        html = [f"<h3><a href='gene:{a}'>{escape(ids[a])}</a> — "
                f"<a href='gene:{b}'>{escape(ids[b])}</a></h3>",
                f"<p style='color:#999'>{escape(products[a])}<br>{escape(products[b])}</p>"]
        if not len(pair):
            html.append("<p>No recorded link joins these two genes.</p>")
        else:
            rows = []
            for r in pair.sort_values("strength", ascending=False).itertuples(index=False):
                layer = str(r.source)
                shown = "" if layer in chosen else " <i>(not ticked)</i>"
                rows.append(
                    f"<tr><td style='color:{LAYER_COLOURS.get(layer, '#ddd')}'>"
                    f"{escape(SE.layer_text(layer))}{shown}</td>"
                    f"<td>stronger than {float(r.strength):.0%} of links of this kind</td>"
                    f"<td style='color:#999'>score {float(r.score):.3g}</td></tr>")
            html.append(f"<p>Linked by <b>{len(pair)}</b> kind{'s' if len(pair) != 1 else ''} "
                        f"of evidence:</p><table cellspacing='4'>{''.join(rows)}</table>")
        la, lb = labels[a], labels[b]
        known_a = la.strip().lower() not in MR.UNLABELLED
        known_b = lb.strip().lower() not in MR.UNLABELLED
        if known_a and known_b:
            same = "the same" if la == lb else "different"
            html.append(f"<p>By {escape(name)} they are {same}: {escape(la)} and {escape(lb)}.</p>")
        else:
            html.append(f"<p>By {escape(name)}: {escape(la) if known_a else '<i>unknown</i>'} and "
                        f"{escape(lb) if known_b else '<i>unknown</i>'}.</p>")
        html.append("<p style='color:#888'>Click either gene to put it in focus.</p>")
        self.reading.setHtml("".join(html))

    #: The kinds of link that are evidence of binding, which the Partners tile reads.
    BINDING_LAYERS = ("xlms", "ip_ms", "unwritten_interaction")

    def show_partners(self) -> None:
        """The Partners tile: the network read by binding evidence alone."""
        for layer, box in self.layer_boxes.items():
            box.blockSignals(True)
            box.setChecked(layer in self.BINDING_LAYERS)
            box.blockSignals(False)
        self.strength.setValue(0)
        self.refresh()
        self.shell.say("Partners: only links that are evidence of binding are ticked.")

    def _typed(self) -> None:
        nodes = self.window.nodes
        ids = nodes["gene_id"].astype(str)
        tokens = [t for t in self.field.text().replace(",", " ").split() if t]
        found = []
        for token in tokens:
            hit = np.flatnonzero(ids.str.lower().to_numpy() == token.lower())
            if not len(hit) and "product" in nodes:
                hit = np.flatnonzero(nodes["product"].astype(str).str.lower()
                                     .str.contains(token.lower(), regex=False).to_numpy())[:1]
            found.extend(int(h) for h in hit)
        self.set_focus(found)

    # ------------------------------------------------------------------ reading
    def read(self) -> None:
        """What the gene in focus is linked to, by kind, and what its links suggest."""
        from . import star_edges as SE
        if not self.focus:
            self.reading.setHtml("<p>Choose a gene, or click one in the network, to read its "
                                 "links.</p>")
            return
        w = self.window
        ids = w.nodes["gene_id"].astype(str).to_numpy()
        products = w.nodes["product"].astype(str).to_numpy() if "product" in w.nodes else ids
        labels = self.labels()
        g = self.focus[0]
        link = NR.linkage(self._links, g)
        name = self.colour_box.currentText()
        html = [f"<h3>{escape(ids[g])}</h3><p style='color:#999'>{escape(products[g])}</p>"]
        if not link.per_layer:
            html.append("<p>No links of the kinds chosen at this strength. Tick more kinds of "
                        "link, or lower the strength.</p>")
        else:
            html.append("<p>" + ", ".join(f"<b>{n}</b> by {escape(SE.layer_text(layer))}"
                                          for layer, n in link.per_layer) + ".</p>")
            rows = "".join(
                f"<tr><td><a href='gene:{int(r.partner)}'>{escape(ids[int(r.partner)])}</a></td>"
                f"<td>{escape(products[int(r.partner)][:50])}</td>"
                f"<td>{escape(labels[int(r.partner)]) or '<i>unknown</i>'}</td>"
                f"<td style='color:{LAYER_COLOURS.get(r.source, '#ddd')}'>{escape(r.source)}</td></tr>"
                for r in link.partners.head(12).itertuples(index=False))
            html.append(f"<table cellspacing='4'><tr><th align='left'>linked to</th>"
                        f"<th align='left'>product</th><th align='left'>{escape(name)}</th>"
                        f"<th align='left'>by</th></tr>{rows}</table>")
            known = [labels[int(p)] for p in link.partners["partner"]
                     if labels[int(p)].strip().lower() not in MR.UNLABELLED]
            if known:
                from collections import Counter
                top, k = Counter(known).most_common(1)[0]
                html.append(f"<p>{k} of its {len(known)} linked genes of known {escape(name)} "
                            f"are <b>{escape(top)}</b>.</p>")
                if self._agreement is None:
                    self._agreement = NR.agreement(self._adj, labels)
                call = self._agreement.for_label(top) if self._agreement is not None else None
                if call is not None:
                    html.append(f"<p style='color:#888'>In this network, when a gene's links "
                                f"mostly say {escape(top)} they are right <b>{call[0]:.0%}</b> of "
                                f"the time ({call[1]:,} such calls); the commonest label alone "
                                f"would be right {self._agreement.chance:.0%}.</p>")
        self.reading.setHtml("".join(html))

    def _link(self, url) -> None:
        text = url.toString()
        if text.startswith("gene:"):
            self.set_focus([int(text.split(":", 1)[1])])

    # ------------------------------------------------------------------ actions
    def _focus_id(self) -> str:
        return str(self.window.nodes["gene_id"].iloc[self.focus[0]]) if self.focus else ""

    def _to_gene(self) -> None:
        if not self.focus:
            return
        from .gene_module import open_in
        open_in(self.window, self._focus_id())

    def _to_map(self) -> None:
        if not self.focus or getattr(self.window, "map_module", None) is None:
            return
        from .map_module import open_in
        open_in(self.window)
        self.window.map_module.module.select(list(self.focus), "")

    def open_bundle(self, folder: str = "") -> dict:
        """Read a bundle saved by this module back and restore its view (64.39)."""
        from . import export_bundle
        return export_bundle.reopen(self, folder)

    def save_bundle(self, folder: str = "") -> str:
        """Write the focus gene's links, as measured, with the network figure; returns the folder."""
        from . import export_bundle
        if not self.focus or self._links is None:
            self.shell.say("choose a gene first; the bundle holds its links")
            return ""
        folder = folder or export_bundle.ask_folder(self.shell)
        if not folder:
            return ""
        ids = self.window.nodes["gene_id"].astype(str).to_numpy()
        table = NR.links_of(self._links, self.focus[0])
        table.insert(0, "gene", ids[self.focus[0]])
        table["partner"] = ids[table["partner"].to_numpy(int)]
        layers = sorted(map(str, pd.unique(self._links["source"]))) if "source" in self._links else []
        export_bundle.write(folder, organism=getattr(self.window.panel, "organism", organisms.TOXOPLASMA),
                            query={"module": "network", "gene": str(ids[self.focus[0]]), "layers": layers},
                            tables={"links": table}, kinds={"links": "measured"}, columns=[],
                            figures={"network": self.view}, nodes=self.window.nodes,
                            strategy="the recorded links of the chosen layers")
        self.shell.say(f"bundle written: {folder}")
        return folder

    def _save_figure(self) -> str:
        """Save the view as a figure (SVG, PDF or a 3x PNG), and say where it went."""
        from .figure_export import ask_and_save
        path = ask_and_save(self.view, self.shell, name=self.shell.key, title=self.figure_title())
        if path:
            self.shell.say(f"figure saved: {path}")
        return path

    def figure_title(self) -> str:
        """The figure title: the module and what the view shows."""
        return f"ParaScape {self.shell.key}"

    def _export(self) -> None:
        if not self.focus or self._links is None:
            return
        ids = self.window.nodes["gene_id"].astype(str).to_numpy()
        table = NR.links_of(self._links, self.focus[0])
        table.insert(0, "gene", ids[self.focus[0]])
        table["partner"] = ids[table["partner"].to_numpy(int)]
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self.shell, "Save the links", f"{ids[self.focus[0]]}_links.csv", "CSV (*.csv)")
        if path:
            table.to_csv(path, index=False)
            self.shell.say(f"Saved {len(table)} links to {os.path.basename(path)}")


def build(window) -> ModuleShell | None:
    """The Network module for `window`."""
    if getattr(window, "nodes", None) is None:
        return None
    module = NetworkModule(window)
    module.refresh()
    return module.shell


def open_in(window, gene: str = "") -> str:
    """Show the Network module over the window, focused on `gene` when given."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "network_module", None)
    if landing is None or shell is None:
        return "the Network module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    if gene:
        shell.module.field.setText(gene)
        shell.module._typed()
    return "Network"


def _open(window, key: str) -> None:
    from .landing_install import open_module
    open_module(window, key)


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
