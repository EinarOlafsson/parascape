"""Replay the coexpression UI on actual checked RNA and protein source tables.

This is a software/reader test, not biological validation or admission of a
cross-condition cohort. Every selected source column and operation is exported.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PyQt6 import QtWidgets
from scipy.stats import spearmanr

from parascape import coexpression as C
from parascape.coexpression_panel import CoexpressionPanel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rna-root', type=Path, required=True)
    parser.add_argument('--deposited-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    panel = CoexpressionPanel(args.rna_root, args.deposited_root)
    panel.resize(1200, 1400)
    proofs = []
    for level in ('RNA', 'Protein'):
        tested = False
        for source in panel.sources:
            if source['level'] != level:
                continue
            matrix = C.load(source, args.rna_root)
            columns, samples = [], set()
            for row in matrix.samples:
                identity = row.get('sample_accession') or row['column']
                if identity not in samples:
                    samples.add(identity); columns.append(row['column'])
            columns = columns[:12]
            if len(columns) < 6:
                continue
            for key in matrix.values.index:
                vector = matrix.values.loc[key, columns].to_numpy(dtype=float)
                if np.isfinite(vector).sum() >= 6 and np.unique(vector[np.isfinite(vector)]).size >= 3:
                    seed = matrix.labels.loc[key]
                    table, record = C.correlate(matrix, seed, columns, seed_row=key)
                    valid = table[table.correlation.notna()]
                    if valid.empty:
                        continue
                    for row in valid.head(20).itertuples(index=False):
                        x = matrix.values.loc[key, columns].to_numpy(dtype=float)
                        y = matrix.values.loc[row.source_row, columns].to_numpy(dtype=float)
                        finite = np.isfinite(x) & np.isfinite(y)
                        expected = spearmanr(x[finite], y[finite]).statistic
                        assert abs(expected - row.correlation) < 1e-12
                    panel.level.setCurrentText(level)
                    for source_index in range(panel.source.count()):
                        if panel.source.itemData(source_index) == source:
                            panel.source.setCurrentIndex(source_index)
                            break
                    else:
                        raise ValueError('Smoke source is absent from actual UI selector')
                    panel._loaded(matrix)
                    panel.seed.setText(str(seed))
                    panel.seed_row.setText(str(key))
                    panel.samples.blockSignals(True)
                    from PyQt6 import QtCore
                    for sample_index, sample in enumerate(matrix.samples):
                        panel.samples.item(sample_index, 0).setCheckState(
                            QtCore.Qt.CheckState.Checked if sample['column'] in columns else QtCore.Qt.CheckState.Unchecked)
                    panel.samples.blockSignals(False)
                    panel._calculated((table, record))
                    panel.show(); app.processEvents()
                    csv = args.output / (level.lower() + '.csv')
                    assert panel.save(str(csv)) == len(table)
                    panel.canvas.print_png(args.output / (level.lower() + '-profiles.png'))
                    assert panel.grab().save(str(args.output / (level.lower() + '-panel.png')))
                    proofs.append(dict(level=level, source=source, rows=len(table),
                        valid_correlations=len(valid), selected_columns=columns,
                        independent_scipy_checks=min(20, len(valid)),
                        csv_sha256=hashlib.sha256(csv.read_bytes()).hexdigest(),
                        biological_cohort_validation=False))
                    tested = True
                    break
            if tested:
                break
        if not tested:
            raise ValueError('No eligible real source for software smoke: ' + level)
    (args.output / 'verification.json').write_text(json.dumps(dict(
        status='actual_source_reader_and_UI_smoke_passed', proofs=proofs,
        scope='Software test only; no validated biological coexpression claim'), indent=2) + '\n')
    panel.close()
    print(json.dumps(proofs))


if __name__ == '__main__':
    main()
