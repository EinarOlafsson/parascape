"""Browse what is known about one gene, and how well it is known (instruction 72, user-facing side).

The information map is a gene x trait table with an aligned probability and an aligned source, and the whole
point of it is the distinction an ordinary table destroys: **an empty cell is not a "no"**. A gene with no
recorded host-facing experiment is untested, not internal; a domain flag that is False was looked for and not
found. So every cell here says which of the three it is, in words -- ``not measured``, a measured negative, or
an observed fact with its probability and the tier of the source that produced it -- and the rows that carry no
measurement are shown in grey and only when asked for, so they can never be mistaken for evidence of absence.

Two views, because Einar asked the map two different questions:

* **Gene view** -- type or search a gene; its tier-1 conclusions (host-facing, compartment, function, partners)
  come first, then the tier-2 measurements they would be learned from, each row with the reference and dataset
  behind it from ``evidence.parquet`` and the gene's recorded partners from ``pairs.parquet``.
* **Trait view** -- pick a column; how many genes carry an observed fact for it, split by source tier, and the
  distribution of the values it takes. Coverage is of cells, so a gene counts once under the strongest source
  that spoke about it.

A banner states what kind of cells these are. As of the 2026-11-04 release **no predicted value exists**: no
strategy has reached the 90 % held-out precision bar, so nothing has been filled in and there is deliberately
no predictions display to look at. The banner is computed from the loaded release
(``PublishedMap.predicted_cells``) rather than written here, so it stops saying so the moment a filled release
is published.

Loading reads several parquet files, one of them 34 MB, so it happens in a worker thread; the panel is usable
(and testable) with a map handed to it directly.
"""
from __future__ import annotations

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from . import information_map as IM
from . import theme as TH

#: Every control's explanation. A release once failed its checks for a missing tooltip, and the costly
#: misreading here is not a mis-set option but an empty cell read as a negative result, so the tips say so.
TIPS = {
    'search': 'Type part of a gene identifier or product name and press Find. Matching genes are listed '
              'below; choosing one shows every trait recorded for it.',
    'find': 'Search the loaded map for genes whose identifier or product matches the text. The map holds the '
            '8,140 ME49 genes; a gene that is absent from the map is not absent from biology.',
    'matches': 'Genes matching the search. Click one to show its traits; double-click to open it in the rest '
               'of the application.',
    'open': 'Show the selected gene in the main gene view, alongside the map.',
    'filter': 'Narrow the trait rows below to those whose column name or description contains this text.',
    'absent': 'Also list the traits with NO measurement for this gene. They are shown in grey and read '
              '"not measured", never "no": nothing was recorded, which is not evidence that the answer is '
              'negative. Off by default so that absence never looks like data.',
    'facts': 'One row per trait. Tier 1 are the conclusions the map would predict; tier 2 are the '
             'measurements they would be learned from. Probability is the calibrated probability that the '
             'fact is true, provisional in this release (tier defaults, not yet measured agreement). '
             'A False value is a measured negative; an empty one is not measured.',
    'partners': 'Recorded binding partners of this gene from the pairs table, strongest first, in both '
                'directions. Crosslink probabilities are 1 - local FDR; pull-down members are provisional '
                'until a background frequency background exists. A missing partner is an unrecorded one.',
    'trait': 'Choose a column of the map to see how much of the proteome carries an observed fact for it and '
             'which values that fact takes.',
    'coverage': 'Genes with an observed fact for this column, split by the source tier of the cell. The tiers '
                'that lost a cell to a stronger source remain in the evidence table, so these are counts of '
                'cells, not of sources.',
    'distribution': 'The values this column takes among genes that have an observed fact, commonest first. '
                    'Genes with no measurement are counted separately above and are not a category here.',
}
FACT_COLUMNS = ('trait', 'tier', 'value', 'probability', 'source tier', 'evidence (reference · dataset)')
PARTNER_COLUMNS = ('partner', 'side', 'probability', 'source tier', 'reference')
#: Grey for a cell with no measurement, amber for a measured negative: the two must not look alike.
NOT_MEASURED_COLOUR = QtGui.QColor('#8b8f98')
NEGATIVE_COLOUR = QtGui.QColor('#c08a2e')
MAX_MATCHES = 300
MAX_PARTNERS = 300


#: Releases already read, by root. The sequence table alone is 34 MB and the contents never change while the
#: application runs, so a second panel (or a reopened dock) costs nothing.
_CACHE: dict = {}
#: Loaders still running, held here rather than only on the panel. A QThread destroyed while running aborts
#: the process, and a panel can be closed before its first load returns.
_RUNNING: set = set()


class _Loader(QtCore.QThread):
    ready = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)

    def __init__(self, root):
        super().__init__(None)
        self.root = root

    def run(self):
        try:
            self.ready.emit(IM.load_published(self.root))
        except (OSError, ValueError, KeyError, TypeError) as error:
            self.failed.emit(str(error))


def value_text(row) -> str:
    """How one fact reads in the table: an observed value, a measured negative, or ``not measured``.

    A separate function because this wording is the panel's entire reason to exist and it is tested directly.
    """
    if not bool(row['observed']):
        return IM.NOT_MEASURED
    value = row['value']
    if isinstance(value, bool) or str(value) in ('True', 'False'):
        return 'yes (observed)' if str(value) == 'True' else 'no — measured negative'
    if isinstance(value, float):
        return f'{value:.4g}'                   # a probability summary is not 16 digits of information
    text = str(value)
    return text if len(text) <= 120 else text[:117] + '…'


class InformationMapPanel(QtWidgets.QWidget):
    """The information map as two views: everything known about one gene, and the coverage of one trait."""

    gene_chosen = QtCore.pyqtSignal(str)        # a gene id, for the rest of the application

    def __init__(self, published=None, root=None, products=None, parent=None):
        super().__init__(parent)
        self.map = published
        self.root = root
        self.products = {str(k): str(v) for k, v in (products or {}).items()}
        self.loader = None
        self.gene = ''
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        self.banner = QtWidgets.QLabel('')
        self.banner.setWordWrap(True)
        self.banner.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.banner.setToolTip(TH.tip(
            'What kind of cells this release holds. Computed from the loaded source table, not written into '
            'the application: it reports predicted cells as soon as a release contains any.'))
        layout.addWidget(self.banner)
        self.status = QtWidgets.QLabel('Loading the published information map…')
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self._gene_page(), 'One gene')
        self.tabs.addTab(self._trait_page(), 'One trait')
        layout.addWidget(self.tabs, 1)
        if self.map is None:
            self.load()
        else:
            self._loaded(self.map)

    # ------------------------------------------------------------------ construction
    def _gene_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(page)
        row = QtWidgets.QHBoxLayout()
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText('Gene identifier or product, e.g. TGME49_208830 or GRA16')
        self.search.setToolTip(TH.tip(TIPS['search']))
        self.find = QtWidgets.QPushButton('Find')
        self.find.setToolTip(TH.tip(TIPS['find']))
        self.open_button = QtWidgets.QPushButton('Show in gene view')
        self.open_button.setToolTip(TH.tip(TIPS['open']))
        row.addWidget(self.search, 1)
        row.addWidget(self.find)
        row.addWidget(self.open_button)
        lay.addLayout(row)
        self.matches = QtWidgets.QListWidget()
        self.matches.setToolTip(TH.tip(TIPS['matches']))
        self.matches.setMaximumHeight(110)
        lay.addWidget(self.matches)
        self.match_note = QtWidgets.QLabel('')
        self.match_note.setWordWrap(True)
        lay.addWidget(self.match_note)
        self.gene_note = QtWidgets.QLabel('Search for a gene to see what is known about it.')
        self.gene_note.setWordWrap(True)
        lay.addWidget(self.gene_note)
        row = QtWidgets.QHBoxLayout()
        self.filter = QtWidgets.QLineEdit()
        self.filter.setPlaceholderText('Filter traits by name or description')
        self.filter.setToolTip(TH.tip(TIPS['filter']))
        self.absent = QtWidgets.QCheckBox('Also show traits with no measurement')
        self.absent.setToolTip(TH.tip(TIPS['absent']))
        row.addWidget(self.filter, 1)
        row.addWidget(self.absent)
        lay.addLayout(row)
        self.facts = QtWidgets.QTableWidget(0, len(FACT_COLUMNS))
        self.facts.setHorizontalHeaderLabels(list(FACT_COLUMNS))
        self.facts.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.facts.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.facts.verticalHeader().setVisible(False)
        self.facts.horizontalHeader().setStretchLastSection(True)
        self.facts.setToolTip(TH.tip(TIPS['facts']))
        lay.addWidget(self.facts, 2)
        self.partner_note = QtWidgets.QLabel('')
        self.partner_note.setWordWrap(True)
        lay.addWidget(self.partner_note)
        self.partners = QtWidgets.QTableWidget(0, len(PARTNER_COLUMNS))
        self.partners.setHorizontalHeaderLabels(list(PARTNER_COLUMNS))
        self.partners.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.partners.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.partners.verticalHeader().setVisible(False)
        self.partners.horizontalHeader().setStretchLastSection(True)
        self.partners.setToolTip(TH.tip(TIPS['partners']))
        lay.addWidget(self.partners, 1)
        self.search.returnPressed.connect(self.find_genes)
        self.find.clicked.connect(self.find_genes)
        self.matches.currentRowChanged.connect(self._match_chosen)
        self.matches.itemDoubleClicked.connect(lambda item: self.gene_chosen.emit(item.data(256) or item.text()))
        self.open_button.clicked.connect(lambda: self.gene and self.gene_chosen.emit(self.gene))
        self.filter.textChanged.connect(self._fill_facts)
        self.absent.toggled.connect(self._fill_facts)
        return page

    def _trait_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(page)
        self.trait = QtWidgets.QComboBox()
        self.trait.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.trait.setMinimumContentsLength(24)
        self.trait.setToolTip(TH.tip(TIPS['trait']))
        lay.addWidget(self.trait)
        self.trait_note = QtWidgets.QLabel('')
        self.trait_note.setWordWrap(True)
        lay.addWidget(self.trait_note)
        lay.addWidget(QtWidgets.QLabel('Coverage: genes with an observed fact, by the source tier of the cell'))
        self.coverage = QtWidgets.QTableWidget(0, 2)
        self.coverage.setHorizontalHeaderLabels(['source tier', 'genes'])
        self.coverage.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.coverage.verticalHeader().setVisible(False)
        self.coverage.horizontalHeader().setStretchLastSection(True)
        self.coverage.setToolTip(TH.tip(TIPS['coverage']))
        lay.addWidget(self.coverage, 1)
        lay.addWidget(QtWidgets.QLabel('Values this trait takes among the genes that have an observed fact'))
        self.distribution = QtWidgets.QTableWidget(0, 2)
        self.distribution.setHorizontalHeaderLabels(['value', 'genes'])
        self.distribution.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.distribution.verticalHeader().setVisible(False)
        self.distribution.horizontalHeader().setStretchLastSection(True)
        self.distribution.setToolTip(TH.tip(TIPS['distribution']))
        lay.addWidget(self.distribution, 1)
        self.trait.currentIndexChanged.connect(self._fill_trait)
        return page

    # ------------------------------------------------------------------ loading
    def load(self):
        """Read the published release in a worker thread; parquet reading must not block the window."""
        if self.loader is not None:
            return
        cached = _CACHE.get(str(self.root))
        if cached is not None:
            self._loaded(cached)
            return
        loader = _Loader(self.root)
        self.loader = loader
        _RUNNING.add(loader)
        loader.ready.connect(self._loaded)
        loader.failed.connect(self._failed)
        loader.finished.connect(lambda: _RUNNING.discard(loader))
        loader.finished.connect(self._finished)
        loader.start()

    def _finished(self):
        self.loader = None

    def _failed(self, message):
        self.status.setText(f'The information map could not be read: {message}')
        self.banner.setText('No map is loaded, so nothing is shown. Nothing is inferred either.')

    def _loaded(self, published):
        self.map = published
        _CACHE.setdefault(str(self.root), published)
        self.banner.setText(published.banner())
        notes = [published.summary()]
        if published.missing:
            notes.append('Not loaded: ' + ' '.join(published.missing))
        self.status.setText(' '.join(notes))
        self.trait.blockSignals(True)
        self.trait.clear()
        for record in published.catalogue.itertuples(index=False):
            label = f'tier {record.tier} · {record.column}'
            self.trait.addItem(label + (f' — {record.description[:70]}' if record.description else ''), record.column)
        self.trait.blockSignals(False)
        if self.trait.count():
            self._fill_trait()
        if self.gene:
            self.show_gene(self.gene)

    # ------------------------------------------------------------------ gene view
    def find_genes(self):
        """List the genes whose identifier or product matches the search text, and show the first."""
        if self.map is None:
            return
        query = self.search.text().strip().casefold()
        self.matches.clear()
        if not query:
            self.match_note.setText('Type part of a gene identifier or product name, then press Find.')
            return
        hits = [gene for gene in self.map.genes.astype(str)
                if query in gene.casefold() or query in self.products.get(gene, '').casefold()]
        exact = [gene for gene in hits if gene.casefold() == query]
        hits = exact + [gene for gene in hits if gene not in exact]
        for gene in hits[:MAX_MATCHES]:
            item = QtWidgets.QListWidgetItem(gene + (f' — {self.products[gene]}' if self.products.get(gene) else ''))
            item.setData(256, gene)
            self.matches.addItem(item)
        if not hits:
            self.match_note.setText(f'No gene of the map matches "{self.search.text().strip()}". '
                                    'A gene absent from the map has no row here; that is not a result about it.')
            self._clear_gene()
            return
        self.match_note.setText(f'{len(hits):,} matching genes'
                                + (f', showing the first {MAX_MATCHES}' if len(hits) > MAX_MATCHES else ''))
        self.matches.setCurrentRow(0)

    def _match_chosen(self, row):
        item = self.matches.item(row) if row is not None and row >= 0 else None
        if item is not None:
            self.show_gene(item.data(256) or item.text())

    def _clear_gene(self):
        self.gene = ''
        self.facts.setRowCount(0)
        self.partners.setRowCount(0)
        self.partner_note.setText('')

    def show_gene(self, gene_id: str):
        """Show every trait recorded for one gene, with its evidence and its recorded partners."""
        if self.map is None:
            return
        if not self.map.has_gene(gene_id):
            self._clear_gene()
            self.gene_note.setText(f'{gene_id} is not a row of this map, so it has no cells here. '
                                   'That is a gap in the map, not a statement about the gene.')
            return
        self.gene = str(gene_id)
        self._fill_facts()
        self._fill_partners()

    def gene_facts(self) -> pd.DataFrame:
        """The rows currently behind the gene table, including the not-measured ones when asked for."""
        if self.map is None or not self.gene:
            return pd.DataFrame()
        facts = self.map.gene_facts(self.gene, include_absent=self.absent.isChecked())
        query = self.filter.text().strip().casefold()
        if query and not facts.empty:
            keep = (facts.column.str.casefold().str.contains(query, regex=False)
                    | facts.description.fillna('').str.casefold().str.contains(query, regex=False))
            facts = facts[keep]
        # Tier 1 (the conclusions) first, then tier 2; within a tier the published column order is kept, so
        # host_facing leads rather than whatever happens to sort first.
        return facts.sort_values('tier', kind='stable')

    def _fill_facts(self):
        if self.map is None or not self.gene:
            return
        facts = self.gene_facts()
        evidence = self.map.gene_evidence(self.gene)
        self.facts.setRowCount(len(facts))
        for i, row in enumerate(facts.to_dict('records')):
            observed = bool(row['observed'])
            text = value_text(row)
            cells = [row['column'], f"tier {row['tier']} · " + ('conclusion' if row['tier'] == 1 else 'measurement'),
                     text,
                     '' if pd.isna(row['probability']) else f"{row['probability']:.2f}",
                     row['source_tier'] or '',
                     self._evidence_text(evidence, row)]
            for j, value in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(str(value))
                if not observed:
                    item.setForeground(NOT_MEASURED_COLOUR)
                    font = item.font(); font.setItalic(True); item.setFont(font)
                    item.setToolTip('Nothing has been measured for this trait on this gene. That is not a '
                                    'negative result: the cell is empty, not "no".')
                elif j == 2 and text.startswith('no —'):
                    item.setForeground(NEGATIVE_COLOUR)
                    item.setToolTip('A source looked and recorded the negative. Unlike an empty cell, this is '
                                    'a measurement.')
                elif j == 5:
                    item.setToolTip(str(value) or 'No long-form evidence row is recorded for this cell.')
                elif j == 0:
                    item.setToolTip(str(row['description'] or '') + f"\n\nfrom: {row['table']}")
                elif j == 2:
                    item.setToolTip(str(row['value']))
                self.facts.setItem(i, j, item)
        self.facts.resizeColumnsToContents()
        total = len(self.map.catalogue)
        observed_rows = int(facts.observed.sum()) if not facts.empty else 0
        self.gene_note.setText(
            f'{self.gene}: {observed_rows:,} of {total:,} traits carry an observed fact; the rest are not '
            'measured (empty, which is never a "no"). '
            + ('Not-measured traits are listed in grey.' if self.absent.isChecked()
               else 'Tick the box to list them as well.'))

    def _evidence_text(self, evidence, row) -> str:
        if row['table'] != 'information map':
            return f"measured in the {row['table']} table"
        if evidence is None or evidence.empty:
            return ''
        base = row['column'][:-2] if str(row['column']).endswith('_n') else row['column']
        rows = evidence[evidence.trait.isin({row['column'], base})]
        parts = [f"{str(r.reference or '').strip()} · {str(r.dataset or '').strip()}" for r in rows.itertuples(index=False)]
        return '; '.join(dict.fromkeys(p for p in parts if p.strip(' ·')))

    def _fill_partners(self):
        partners = self.map.gene_partners(self.gene, MAX_PARTNERS) if self.map is not None else pd.DataFrame()
        self.partners.setRowCount(len(partners))
        for i, row in enumerate(partners.to_dict('records')):
            for j, key in enumerate(('partner', 'side', 'probability', 'source_tier', 'reference')):
                value = row[key]
                text = f'{value:.3f}' if key == 'probability' and not pd.isna(value) else str(value)
                item = QtWidgets.QTableWidgetItem(text)
                item.setToolTip(str(row['reference']))
                self.partners.setItem(i, j, item)
        self.partners.resizeColumnsToContents()
        total = len(self.map.gene_partners(self.gene)) if self.map is not None else 0
        if total:
            shown = f', showing the strongest {MAX_PARTNERS}' if total > MAX_PARTNERS else ''
            note = (f'{total:,} recorded partners{shown}. A partner that is not listed is one that has not '
                    'been measured, not one that has been ruled out.')
        else:
            note = ('No binding partner is recorded for this gene: none has been measured, which is not the '
                    'same as having none.')
        self.partner_note.setText(note)

    # ------------------------------------------------------------------ trait view
    def _fill_trait(self):
        if self.map is None or not self.trait.count():
            return
        column = self.trait.currentData()
        record = self.map.catalogue[self.map.catalogue.column == column]
        coverage = self.map.coverage(column)
        distribution = self.map.distribution(column)
        genes = len(self.map.genes)
        observed = int(coverage.genes.max()) if not coverage.empty else 0
        description = record.description.iloc[0] if not record.empty else ''
        tier = int(record.tier.iloc[0]) if not record.empty else 0
        self.trait_note.setText(
            f'{column} — tier {tier} ' + ('conclusion' if tier == 1 else 'measurement')
            + (f'. {description}' if description else '')
            + f'. {observed:,} of {genes:,} genes carry an observed fact; the other {genes - observed:,} are '
              'not measured, which is not a negative value and is not counted as one below.')
        for table, frame, keys in ((self.coverage, coverage, ('source_tier', 'genes')),
                                   (self.distribution, distribution, ('value', 'genes'))):
            table.setRowCount(len(frame))
            for i, row in enumerate(frame.to_dict('records')):
                for j, key in enumerate(keys):
                    text = str(row[key])
                    # A sequence column's values are whole proteins; the cell shows the start and the
                    # tooltip the rest, rather than a row the width of the screen.
                    item = QtWidgets.QTableWidgetItem(text if len(text) <= 120 else text[:117] + '…')
                    item.setToolTip(text[:4000])
                    table.setItem(i, j, item)
            table.resizeColumnsToContents()

    def show_trait(self, column: str) -> bool:
        """Open the trait view on one column; False when the loaded map has no such column."""
        index = self.trait.findData(column)
        if index < 0:
            return False
        self.tabs.setCurrentIndex(1)
        self.trait.setCurrentIndex(index)
        return True

    def closeEvent(self, event):
        if self.loader is not None:
            self.loader.wait(2000)
        super().closeEvent(event)


def install(window) -> InformationMapPanel | None:
    """Add the information-map dock beside Discoveries and wire its gene selection to the gene view."""
    products = {}
    nodes = getattr(window, 'nodes', None)
    if nodes is not None and 'gene_id' in nodes:
        products = dict(zip(nodes['gene_id'].astype(str),
                            nodes['product'].astype(str) if 'product' in nodes else [''] * len(nodes)))
    panel = InformationMapPanel(products=products)
    if hasattr(window, '_workflow_gene'):
        panel.gene_chosen.connect(window._workflow_gene)
    dock = QtWidgets.QDockWidget('information map', window)
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    scroll = QtWidgets.QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
    scroll.setMinimumWidth(120)
    scroll.setWidget(panel)
    dock.setWidget(scroll)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    anchor = getattr(window, 'discoveries_dock', None) or getattr(window, 'strategies_dock', None)
    if anchor is not None:
        window.tabifyDockWidget(anchor, dock)
    if getattr(window, 'right_dock', None) is not None:
        window.right_dock.raise_()
    window.information_map, window.information_map_dock = panel, dock
    return panel
