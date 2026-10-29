"""Browse studies, source-declared conditions, submitted samples and linked runs.

Study-level summaries can hide missing controls and repeated sequencing of one sample. This hierarchy keeps each run attached to its submitted sample and source annotations, marks pending measurements explicitly, and leaves suggested condition groups and unverified pairing available for review before a comparison.
"""
from __future__ import annotations

import html

from PyQt6 import QtCore, QtWidgets

from . import expression_contrast as X
from . import rna_design as D


class StudyDesignPanel(QtWidgets.QWidget):
    """A searchable study → condition → sample → run tree with source evidence."""

    def __init__(self, root: str = X.NAS, parent=None):
        super().__init__(parent)
        self.root = root
        self.design = D.samples(root)
        layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel('Browse study → condition → submitted sample → sequencing run. '
                                 'Control labels and replicate labels retain their source. '
                                 'Suggested groups need review; a repeated run is not an independent biological replicate.')
        label.setWordWrap(True)
        layout.addWidget(label)
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText('Find a study, sample, control, replicate or run')
        layout.addWidget(self.search)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderLabels(['Study / condition / sample / run', 'Role or status'])
        self.details = QtWidgets.QTextBrowser()
        self.details.setOpenExternalLinks(True)
        split.addWidget(self.tree)
        split.addWidget(self.details)
        layout.addWidget(split, 1)
        self.export = QtWidgets.QPushButton('Export sample design…')
        self.export.setToolTip('Export all planned samples, processing status, source control and replicate labels, and linked run relationships to CSV.')
        self.export.clicked.connect(lambda: self.save())
        layout.addWidget(self.export)
        self._populate()
        self.tree.currentItemChanged.connect(self._selected)
        self.search.textChanged.connect(self._filter)

    def _populate(self):
        if self.design.empty:
            self.details.setPlainText('No sample metadata is available for this pack.')
            return
        from .compare_panel import catalogue_titles
        titles = catalogue_titles()
        from .rna_profiles import snapshot
        priority = snapshot(self.root).get('sample_design', {}).get('pending_priority_study')
        design = self.design.assign(_priority=self.design.study != priority).sort_values(['_priority', 'study'])
        for study, rows in design.groupby('study', sort=False):
            complete = int((rows.get('processing_status', 'verified') == 'verified').sum()) if 'processing_status' in rows else len(rows)
            title = titles.get(study, '')
            top = QtWidgets.QTreeWidgetItem(self.tree, [f'{study} — {title}', ('Priority · ' if study == priority else '') + f'{complete}/{len(rows)} runs verified'])
            top.setData(0, QtCore.Qt.ItemDataRole.UserRole, {'study': study, 'coverage': f'{complete}/{len(rows)} planned runs verified', 'note': 'Select a sample or run to inspect source metadata and relationships.'})
            for condition, group in rows.groupby('condition_group', sort=True):
                item = QtWidgets.QTreeWidgetItem(top, [str(condition), f'{group.sample_key.nunique()} submitted samples'])
                item.setData(0, QtCore.Qt.ItemDataRole.UserRole, {'condition': condition, 'group_evidence': '; '.join(sorted(set(group.group_evidence))), 'note': 'Control pairing is not automatically inferred from a matching replicate number.'})
                for sample, members in group.groupby('sample_key', sort=True):
                    first = members.iloc[0]
                    node = QtWidgets.QTreeWidgetItem(item, [f"{sample} — {first.get('sample_title', '')}", str(first.control_role)])
                    node.setData(0, QtCore.Qt.ItemDataRole.UserRole, first.to_dict())
                    for row in members.to_dict('records'):
                        leaf = QtWidgets.QTreeWidgetItem(node, [row['run_accession'], str(row.get('processing_status', 'verified'))])
                        leaf.setData(0, QtCore.Qt.ItemDataRole.UserRole, row)
        self.tree.header().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.tree.header().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)

    def _selected(self, item, previous):
        if not item:
            return
        row = item.data(0, QtCore.Qt.ItemDataRole.UserRole) or {}
        fields = ['study', 'coverage', 'sample_title', 'sample_accession', 'sample_alias',
                  'condition_group', 'group_evidence', 'control_role', 'control_link',
                  'replicate_label', 'replicate_type', 'replicate_evidence', 'linked_runs',
                  'run_relationship', 'shared_across_studies', 'processing_status',
                  'metadata_status', 'attributes_json', 'note']
        fields += ['library_strategy', 'library_layout', 'instrument_model', 'read_count']
        body = ''.join(f'<p><b>{html.escape(k.replace("_", " "))}:</b> {html.escape(str(row[k]))}</p>' for k in fields if k in row)
        url = str(row.get('metadata_url', ''))
        if url.startswith('https://www.ebi.ac.uk/ena/'):
            body += f'<p><a href="{html.escape(url, quote=True)}">Open source sample record</a></p>'
        self.details.setHtml(body)

    def _filter(self, text):
        query = text.lower().strip()

        def visit(item, inherited=False):
            own = inherited or query in (item.text(0) + ' ' + item.text(1)).lower()
            child_matches = [visit(item.child(i), own) for i in range(item.childCount())]
            show = own or any(child_matches)
            item.setHidden(not show)
            if query and show:
                item.setExpanded(True)
            return show

        for i in range(self.tree.topLevelItemCount()):
            visit(self.tree.topLevelItem(i))

    def save(self, path: str = '') -> int:
        """Export source annotations, suggested groups and run/sample relationships to CSV."""
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Export sample design', 'sample_design.csv', 'CSV (*.csv)')
        if not path:
            return 0
        self.design.to_csv(path, index=False)
        return len(self.design)
