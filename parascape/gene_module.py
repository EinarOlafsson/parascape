"""The Gene module: everything recorded about one gene, synthesised, with its sources.

Built with the parascape-module mould (instruction 96); Einar chose "synthesis first"
(2026-11-06). The module lets you take a gene -- or step through a list -- and read what is
known about it in sentences (where it is, what it binds, how much the parasite needs it, when
it is made...) each with its sources, then the evidence itself grouped the same way, then what
is NOT known and what can be inferred, with how often such a call is right.

The facts come from the published information map (`information_map.PublishedMap`, loaded by
the information-map panel) and the window's own gene table; `gene_reading` turns them into the
sentences, so every sentence is testable without a window.
"""
from __future__ import annotations

import os
from html import escape

import numpy as np
import pandas as pd
from PySide6 import QtCore, QtWidgets

from . import gene_reading as GR
from . import theme as TH
from .shell import ModuleShell


def _muted(text: str) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setObjectName("ShellMuted")
    label.setWordWrap(True)
    return label


def _browser() -> QtWidgets.QTextBrowser:
    view = QtWidgets.QTextBrowser()
    view.setOpenExternalLinks(True)
    return view


class GeneModule(QtCore.QObject):
    """The Gene module for one window."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.gene = ""
        self.genes: list = []                          # a list being stepped through
        self.reading = None
        self.shell = shell = ModuleShell(
            "gene", "Gene", "everything recorded about one gene, with its sources",
            instruction="Type a gene, or choose a list and step through it.")
        shell.module = self
        shell.setAcceptDrops(True)
        shell.installEventFilter(self)
        self._build_settings()
        self._build_runtime()
        shell.home_requested.connect(lambda: _go_home(window))
        self._waiting = QtCore.QTimer(self)
        self._waiting.setInterval(300)
        self._waiting.timeout.connect(self._map_ready)

    # ------------------------------------------------------------------ left
    def _build_settings(self) -> None:
        shell = self.shell
        mine = QtWidgets.QWidget(shell)
        v = QtWidgets.QVBoxLayout(mine)
        v.setContentsMargins(6, 2, 6, 6)
        self.field = QtWidgets.QLineEdit()
        self.field.setPlaceholderText("A gene id or a word of its product, then Return")
        self.field.setToolTip(TH.tip("TGME49_208830, or GRA16, or 'rhoptry neck'. Several "
                                     "matches are listed below to pick from."))
        self.field.returnPressed.connect(self._typed)
        v.addWidget(self.field)
        row = QtWidgets.QHBoxLayout()
        choose = QtWidgets.QPushButton("Choose a list…")
        choose.setToolTip(TH.tip("A table of gene ids; step through it with the list below. "
                                 "You can also drop the file on this module."))
        choose.clicked.connect(self._choose_list)
        use = QtWidgets.QPushButton("Use the gene list")
        use.setToolTip("The list last dropped on ParaScape (the gene list panel).")
        use.clicked.connect(self._use_gene_list)
        row.addWidget(choose)
        row.addWidget(use)
        v.addLayout(row)
        self.matches = QtWidgets.QListWidget()
        self.matches.setToolTip("Click a gene to read it. Up and Down step through the list.")
        self.matches.setMaximumHeight(170)
        self.matches.currentItemChanged.connect(
            lambda item, _old: self.show(item.data(QtCore.Qt.ItemDataRole.UserRole)) if item else None)
        v.addWidget(self.matches)
        self.match_note = _muted("")
        v.addWidget(self.match_note)
        shell.add_setting("Your gene", mine)

        kinds = QtWidgets.QWidget(shell)
        kl = QtWidgets.QVBoxLayout(kinds)
        kl.setContentsMargins(6, 2, 6, 6)
        kl.addWidget(_muted("What to read about it:"))
        self.include: dict = {}
        for category, title in GR.ORDER:
            box = QtWidgets.QCheckBox(title)
            box.setChecked(category != "sequence")
            box.setToolTip(TH.tip(f"Include {title.lower()} in the summary and the evidence."))
            box.toggled.connect(lambda *_: self.render())
            kl.addWidget(box)
            self.include[category] = box
        shell.add_setting("What to read", kinds)

        advanced = QtWidgets.QWidget(shell)
        al = QtWidgets.QVBoxLayout(advanced)
        al.setContentsMargins(6, 2, 6, 6)
        self.show_predicted = QtWidgets.QCheckBox("Show predicted values (never mixed with facts)")
        self.show_predicted.setChecked(True)
        self.show_predicted.setToolTip(TH.tip("Predictions appear only under 'What can be "
                                              "inferred', with their certainty."))
        self.show_predicted.toggled.connect(lambda *_: self.render())
        al.addWidget(self.show_predicted)
        shell.add_setting("Advanced", advanced, expanded=False)

    # ------------------------------------------------------------------ right
    def _build_runtime(self) -> None:
        shell, window = self.shell, self.window
        self.known = _browser()
        shell.link_style(self.known)
        self.known.setToolTip("What is known about the gene, one sentence per kind of "
                              "knowledge, each with its sources.")
        shell.add_section("What is known", self.known, extent=300, stretch=3)

        self.unusual = _browser()
        shell.link_style(self.unusual)
        self.unusual.setToolTip("Of everything measured on this gene, what is furthest from other "
                                "genes: each as a percentile among every gene measured, with its "
                                "source. Study effort and identifiers are left out.")
        shell.add_section("What is unusual about it", self.unusual, extent=170, stretch=1)

        holder = QtWidgets.QWidget(shell)
        hl = QtWidgets.QVBoxLayout(holder)
        hl.setContentsMargins(0, 0, 0, 0)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.evidence_page = QtWidgets.QWidget()
        self.evidence_layout = QtWidgets.QVBoxLayout(self.evidence_page)
        self.evidence_layout.setContentsMargins(8, 4, 8, 8)
        self.evidence_layout.addStretch(1)
        scroll.setWidget(self.evidence_page)
        hl.addWidget(scroll)
        shell.add_section("The evidence", holder, extent=260, stretch=2)

        self.inferred = _browser()
        shell.link_style(self.inferred)
        self.inferred.setOpenLinks(False)
        self.inferred.anchorClicked.connect(self._inferred_link)
        self.inferred.setToolTip("What is not known about the gene, and what the data suggest, "
                                 "each with how often such a call is right.")
        shell.add_section("What is not known, and what can be inferred", self.inferred,
                          extent=170, stretch=1)

        self.methods = _browser()
        shell.link_style(self.methods)
        self.methods.setOpenLinks(False)
        self.methods.anchorClicked.connect(lambda url: self.open_strategy(url.toString().split(":", 1)[-1]))
        self.methods.setToolTip("Every method in the catalogue: what it said about this gene when the "
                                "gene's own answer was hidden, and how its calibration graded it.")
        shell.add_section("What each method said about it, held out", self.methods, extent=200,
                          stretch=1, folded=True)

        shell.add_action("Put it on the map", self._on_map, primary=True,
                         tip="Open the Map module with this gene marked and its neighbours read.")
        shell.add_action("Its partners in the network", lambda: self._open("network"),
                         tip="See the gene among its recorded partners (Network).")
        shell.add_action("Predict what is missing", lambda: self._open("infer"),
                         tip="Answer a question about this gene from the evidence, with how often "
                             "such answers are right (Infer).")
        shell.add_action("Open a saved bundle…", lambda: self.open_bundle(),
                         tip="Read a bundle saved here back: its files are checked against their hashes "
                             "and its view is restored.")
        shell.add_action("Save everything as a bundle…", lambda: self.save_bundle(),
                         tip="One folder: what is measured about the gene, what is unusual, your kept "
                             "claims (marked inferred), the methods text with references, and a "
                             "manifest whose hashes show if anything was edited later.")
        shell.add_action("Full evidence panel", self._full_panel,
                         tip="The detailed evidence panel of the main window for this gene.")
        shell.add_action("Previous", lambda: self._step(-1), tip="The previous gene of the list.")
        shell.add_action("Next", lambda: self._step(1), tip="The next gene of the list.")
        shell.borrow("Console and AI chat", [("Console", getattr(window, "console_dock", None)),
                                             ("AI chat", getattr(window, "chat_dock", None))],
                     folded=True)
        self.known.setHtml("<p>Type a gene on the left to read what is known about it.</p>")

    # ------------------------------------------------------------------ choosing a gene
    def _nodes(self):
        return getattr(self.window, "nodes", None)

    def _product(self, gene: str) -> str:
        nodes = self._nodes()
        if nodes is None or "product" not in nodes:
            return ""
        hit = nodes[nodes["gene_id"].astype(str) == gene]
        return str(hit["product"].iloc[0]) if len(hit) else ""

    def matching(self, text: str) -> list:
        """Gene ids matching typed text: an exact id first, then ids and products containing it."""
        nodes = self._nodes()
        text = (text or "").strip()
        if nodes is None or not text:
            return []
        ids = nodes["gene_id"].astype(str)
        low = text.lower()
        exact = ids[ids.str.lower() == low].tolist()
        contains = ids[ids.str.lower().str.contains(low, regex=False)].tolist()
        if "product" in nodes:
            prod = nodes["product"].astype(str).str.lower()
            contains += ids[prod.str.contains(low, regex=False)].tolist()
        return list(dict.fromkeys(exact + contains))[:200]

    def _fill_matches(self, genes: list, note: str) -> None:
        self.matches.blockSignals(True)
        self.matches.clear()
        for gene in genes:
            item = QtWidgets.QListWidgetItem(f"{gene} — {self._product(gene)[:60]}")
            item.setData(QtCore.Qt.ItemDataRole.UserRole, gene)
            self.matches.addItem(item)
        self.matches.blockSignals(False)
        self.match_note.setText(note)
        if genes:
            self.matches.setCurrentRow(0)

    def _typed(self) -> None:
        genes = self.matching(self.field.text())
        self.genes = genes
        self._fill_matches(genes, f"{len(genes)} matching gene{'s' if len(genes) != 1 else ''}"
                           if genes else "No gene matches that.")

    def load_list(self, path: str) -> None:
        """Read a file of genes and step through them."""
        from . import gene_drop as GD
        nodes = self._nodes()
        known = nodes["gene_id"].astype(str).tolist() if nodes is not None else None
        try:
            drop = GD.read_drop(path, known_ids=known)
        except Exception as error:
            self.shell.say(f"{os.path.basename(path)} could not be read: {error}")
            return
        genes = [g for g in drop.genes if known is None or g in set(known)]
        self.genes = genes
        self._fill_matches(genes, f"{len(genes)} genes from {os.path.basename(path)}")

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
        self.genes = list(genes)
        self._fill_matches(self.genes, f"{len(genes)} genes from the gene list")

    def _step(self, by: int) -> None:
        row = self.matches.currentRow()
        if self.matches.count():
            self.matches.setCurrentRow(max(0, min(self.matches.count() - 1, row + by)))

    def eventFilter(self, obj, event):                    # noqa: N802 - Qt name
        """Accept a file dragged onto the module and step through its genes."""
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

    # ------------------------------------------------------------------ reading it
    def _published(self):
        panel = getattr(self.window, "information_map", None)
        if panel is None:
            return None
        if panel.map is None:
            panel.load()
        return panel.map

    def _map_ready(self) -> None:
        if self._published() is not None:
            self._waiting.stop()
            if self.gene:
                self.show(self.gene)

    def show(self, gene: str) -> None:
        """Read one gene and show it."""
        if not gene:
            return
        self.gene = gene
        published = self._published()
        if published is None:
            self._waiting.start()
            self.known.setHtml(f"<p>Reading the information map for {escape(gene)}…</p>")
        nodes = self._nodes()
        row = None
        if nodes is not None:
            hit = nodes[nodes["gene_id"].astype(str) == gene]
            row = hit.iloc[0] if len(hit) else None
        self.reading = GR.read(published, gene, self._product(gene), node_row=row)
        self.render()

    def _wanted(self, category: str) -> bool:
        box = self.include.get(category)
        return box is None or box.isChecked()

    def render(self) -> None:
        """Write the synthesis, the evidence sections and the inferences for the current gene."""
        r = self.reading
        if r is None:
            return
        sentences = [(c, s) for c, s in GR.synthesis(r) if self._wanted(c)]
        head = (f"<h2>{escape(r.gene)}</h2><p style='color:#999'>{escape(r.product)}</p>"
                if r.product else f"<h2>{escape(r.gene)}</h2>")
        body = "".join(f"<p>{s}</p>" for _c, s in sentences) or \
            "<p>Nothing recorded about this gene in the kinds of knowledge chosen.</p>"
        self.known.setHtml(head + body)
        self._render_unusual(r)
        self._render_methods(r)
        self._render_evidence(r)
        self._render_inferences(r)

    def open_bundle(self, folder: str = "") -> dict:
        """Read a bundle saved by this module back and restore its view (64.39)."""
        from . import export_bundle
        return export_bundle.reopen(self, folder)

    def save_bundle(self, folder: str = "") -> str:
        """Write the current gene as an export bundle; returns the folder."""
        from . import export_bundle
        from . import kept_claims
        r = self.reading
        if r is None:
            self.shell.say("choose a gene first; the bundle holds its reading")
            return ""
        folder = folder or export_bundle.ask_folder(self.shell)
        if not folder:
            return ""
        frames = [f.assign(category=c) for c, f in (r.facts or {}).items() if f is not None and len(f)]
        facts = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        tables, kinds = {"what is measured": facts}, {"what is measured": "measured"}
        unusual = GR.unusual(self._nodes(), r.gene)
        if len(unusual):
            tables["what is unusual"], kinds["what is unusual"] = unusual, "computed"
        if r.predicted is not None and len(r.predicted):
            tables["predicted cells"], kinds["predicted cells"] = r.predicted, "inferred"
        kept = kept_claims.claims(r.gene)
        if len(kept):
            tables["your kept claims"], kinds["your kept claims"] = kept, "inferred"
        columns = [str(c) for c in facts.get("column", pd.Series(dtype=str)).dropna().unique()
                   if c in self._nodes().columns] if len(facts) else []
        export_bundle.write(folder, organism=getattr(self.window.panel, "organism", "Tg"),
                            query={"module": "gene", "gene": r.gene, "product": r.product},
                            tables=tables, kinds=kinds, columns=columns, nodes=self._nodes())
        self.shell.say(f"bundle written: {folder}")
        return folder

    def _render_methods(self, r) -> None:
        """Every catalogue method and its held-out calls on this gene (instruction 64.30)."""
        try:
            profile = GR.strategy_profile(r.gene)
        except Exception as error:
            self.methods.setHtml(f"<p>The hold-out record could not be read: {escape(str(error))}</p>")
            return
        spoke = profile[profile["status"] == "held-out calls"]
        head = (f"<p><b>{len(spoke)}</b> of the {len(profile)} methods were tested on this gene with its "
                "answer hidden; the others are listed with why.</p>")
        rows = "".join(
            f"<tr><td><a href='strategy:{escape(p.key)}'>{escape(p.strategy)}</a></td><td>{escape(p.grade)}</td>"
            f"<td>{escape(p.calls or p.status)}</td></tr>" for p in profile.itertuples(index=False))
        self.methods.setHtml(head + "<table cellspacing=4><tr><th align=left>method</th>"
                             "<th align=left>graded</th><th align=left>held out, it said</th></tr>"
                             + rows + "</table>")

    def open_strategy(self, key: str) -> bool:
        """One click from a method to its card and scorecard: the Strategies module, that strategy shown."""
        panel = getattr(self.window, "strategy_panel", None)
        if panel is None or not key:
            return False
        try:
            panel.select(key)
        except KeyError:
            return False
        from .dock_module import open_in
        open_in(self.window, "strategies")
        return True

    def _render_unusual(self, r) -> None:
        """The gene's most extreme measurements, among the kinds of knowledge chosen."""
        frame = GR.unusual(self._nodes(), r.gene, top=30)
        frame = frame[[self._wanted(c) for c in frame["category"]]].head(10) if len(frame) else frame
        lines = GR.unusual_sentences(frame)
        self.unusual.setHtml("".join(f"<p>{s}</p>" for s in lines) or
                             "<p>Nothing measured on this gene stands out from other genes, or too "
                             "little is measured on it to say.</p>")

    def _render_evidence(self, r) -> None:
        while self.evidence_layout.count() > 1:
            item = self.evidence_layout.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
        for category, title in GR.ORDER:
            rows = r.facts.get(category)
            if rows is None or not len(rows) or not self._wanted(category):
                continue
            rows = GR.informative(rows) if category not in ("place",) else rows
            if not len(rows):
                continue
            table = QtWidgets.QTableWidget(len(rows), 4)
            table.setHorizontalHeaderLabels(["measurement", "value", "source", "how sure"])
            table.verticalHeader().hide()
            table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
            for i, row in enumerate(rows.itertuples(index=False)):
                source = GR._text(getattr(row, "reference", "")) or \
                    "; ".join(r.sources_for(str(row.column), 2)) or GR._text(getattr(row, "dataset", "")) \
                    or GR._text(row.source_tier)
                p = getattr(row, "probability", float("nan"))
                dataset = GR._text(getattr(row, "dataset", ""))
                what = (f"{dataset} · {str(row.column).replace('_', ' ')}" if dataset
                        else GR._text(row.description) or str(row.column))
                cells = (what,
                         GR.value_text(str(row.column), row.value), source,
                         "" if p != p else f"{p:.0%}")
                for j, text in enumerate(cells):
                    table.setItem(i, j, QtWidgets.QTableWidgetItem(text))
            table.resizeColumnsToContents()
            table.horizontalHeader().setStretchLastSection(True)
            table.setMinimumHeight(min(260, 30 + 24 * len(rows)))
            section = self.shell.section(f"{title} ({len(rows)})", table,
                                         expanded=category in ("place", "interaction", "fitness"))
            self.evidence_layout.insertWidget(self.evidence_layout.count() - 1, section)

    def _render_inferences(self, r) -> None:
        lines = [s for s in GR.inferences(r)
                 if self.show_predicted.isChecked() or not s.startswith("Predicted")]
        suggestion = self._map_suggestion(r.gene)
        if suggestion:
            lines.append(suggestion + " <a href='keep:map'>keep this</a>")
        lines.extend(self._kept_lines(r.gene))
        self.inferred.setHtml("".join(f"<p>{s}</p>" for s in lines)
                              or "<p>Every key trait of this gene is measured.</p>")

    def _kept_lines(self, gene: str) -> list:
        """The user's own kept claims for this gene -- shown as kept, never as facts -- and any
        that a release now measures, said loudly (instruction 78)."""
        from . import kept_claims as KC
        kept = KC.claims(gene)
        if not len(kept):
            return []
        superseded = {}
        published = self._published()
        if published is not None:
            checked = KC.check_against(published)
            superseded = {row.claim_id: row for row in checked.itertuples(index=False)}
        out = []
        for c in kept.itertuples(index=False):
            p = "" if c.probability != c.probability else f", {float(c.probability):.0%}"
            line = (f"<b>Kept by you</b>: {escape(str(c.trait))} = <b>{escape(str(c.value))}</b> "
                    f"({escape(str(c.method))}{p}; release {escape(str(c.release))}, "
                    f"{escape(str(c.kept_at)[:10])}) <a href='forget:{c.claim_id}'>forget</a>")
            hit = superseded.get(c.claim_id)
            if hit is not None:
                line += (f"<br><b style='color:#d9a441'>Now measured</b>: {escape(hit.measured_value)}"
                         f" -- the kept claim {'agrees' if hit.agrees else 'DISAGREES'}; the "
                         f"measurement replaces it.")
            out.append(line)
        return out

    def _inferred_link(self, url) -> None:
        from PySide6 import QtGui
        text = url.toString()
        from . import kept_claims as KC
        if text.startswith("forget:"):
            KC.revert(text.split(":", 1)[1])
            self.render()
        elif text == "keep:map":
            self.keep_map_suggestion()
        else:
            QtGui.QDesktopServices.openUrl(url)

    def keep_map_suggestion(self):
        """Keep what the gene's map neighbours suggest, with its hit rate, in your own layer."""
        from . import kept_claims as KC
        from . import map_reading as MR
        map_module = getattr(getattr(self.window, "map_module", None), "module", None)
        nodes = self._nodes()
        if map_module is None or nodes is None or not self.gene:
            return None
        index = np.flatnonzero(nodes["gene_id"].astype(str).to_numpy() == self.gene)
        if not len(index):
            return None
        hood = MR.neighbours(self.window.xyz, map_module.labels(), [int(index[0])],
                             getattr(self.window, "placed", None))
        if not hood.top:
            return None
        label, n, of = hood.top
        agree = map_module.agreement()
        call = agree.for_label(label) if agree is not None else None
        claim = KC.keep(KC.KeptClaim(
            gene=self.gene, trait=map_module.label_name(), value=label,
            probability=call[0] if call else float("nan"),
            uncertainty=(f"right {call[0]:.0%} of {call[1]} such calls on this map; chance "
                         f"{agree.chance:.0%}") if call else "no record on this map",
            method="map neighbours", method_version=f"k={MR.K}",
            question=f"{map_module.label_name()} from its {of} nearest labelled genes on the map",
            accepted_score=call[0] if call else float("nan")))
        self.shell.say(f"Kept: {self.gene} {claim.trait} = {claim.value}")
        self.render()
        return claim

    def _map_suggestion(self, gene: str) -> str:
        """What the gene's neighbours on the current map suggest, with how often that is right."""
        map_module = getattr(getattr(self.window, "map_module", None), "module", None)
        nodes = self._nodes()
        if map_module is None or nodes is None or getattr(self.window, "xyz", None) is None:
            return ""
        index = np.flatnonzero(nodes["gene_id"].astype(str).to_numpy() == gene)
        if not len(index):
            return ""
        from . import map_reading as MR
        labels = map_module.labels()
        hood = MR.neighbours(self.window.xyz, labels, [int(index[0])],
                             getattr(self.window, "placed", None))
        if not hood.top:
            return ""
        label, n, of = hood.top
        name = map_module.label_name()
        text = (f"On the current map, {n} of its {of} nearest genes of known {escape(name)} are "
                f"<b>{escape(label)}</b>")
        agree = map_module.agreement()
        call = agree.for_label(label) if agree is not None else None
        if call is not None:
            text += (f"; when neighbours say {escape(label)} they are right {call[0]:.0%} of the "
                     f"time on this map ({agree.chance:.0%} by chance)")
        return text + "."

    # ------------------------------------------------------------------ actions
    def _on_map(self) -> None:
        nodes = self._nodes()
        map_shell = getattr(self.window, "map_module", None)
        if map_shell is None or nodes is None or not self.gene:
            return
        index = np.flatnonzero(nodes["gene_id"].astype(str).to_numpy() == self.gene)
        from .map_module import open_in
        open_in(self.window)
        if len(index):
            map_shell.module.select([int(index[0])], "")

    def _open(self, key: str) -> None:
        from .landing_install import open_module
        if self.gene and hasattr(self.window, "_workflow_gene"):
            try:
                self.window._workflow_gene(self.gene)
            except Exception:
                pass
        open_module(self.window, key)

    def _full_panel(self) -> None:
        from .landing_install import leave_landing
        leave_landing(self.window)
        if self.gene and hasattr(self.window, "_workflow_gene"):
            self.window._workflow_gene(self.gene)
        dock = getattr(self.window, "right_dock", None)
        if dock is not None:
            dock.show()
            dock.raise_()


def build(window) -> ModuleShell | None:
    """The Gene module for `window`."""
    if getattr(window, "nodes", None) is None:
        return None
    return GeneModule(window).shell


def open_in(window, gene: str = "") -> str:
    """Show the Gene module over the window, reading `gene` when given."""
    from .landing_install import enter_landing
    landing = getattr(window, "landing", None)
    shell = getattr(window, "gene_module", None)
    if landing is None or shell is None:
        return "the Gene module is not available in this window"
    landing.show_module(shell)
    enter_landing(window)
    if gene:
        shell.module.field.setText(gene)
        shell.module._typed()
    return "Gene"


def _go_home(window) -> None:
    from .landing_install import go_home
    go_home(window)
