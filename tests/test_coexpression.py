"""Coexpression numerical, source identity and Qt regression checks."""
import hashlib
import json

import numpy as np
import pandas as pd
import pytest
from scipy.stats import pearsonr, spearmanr

from parascape import coexpression as C


@pytest.fixture
def matrix():
    values = pd.DataFrame([
        [0, 1, 2, 3, 4, 5], [0, 2, 4, 6, 8, 10], [5, 4, 3, 2, 1, 0],
        [1, np.nan, 5, 2, np.inf, 3], [2, 2, 2, 2, 2, 2],
        [np.nan, np.nan, 1, 2, np.nan, np.nan]],
        index=[10, 11, 12, 13, 14, 15], columns=list('abcdef'), dtype=float)
    return C.ExpressionMatrix(values, pd.Series(['seed', 'same', 'opposite', 'missing', 'constant', 'short'],
        index=values.index), [dict(column=c, source_sample='original ' + c) for c in values],
        dict(quantity='Original LFQ', units='Original units'))


@pytest.mark.parametrize('method', ['spearman', 'pearson'])
def test_correlations_agree_with_independent_scipy_on_each_overlap(matrix, method):
    table, record = C.correlate(matrix, 'seed', list('abcdef'), method=method, min_pairs=4)
    table = table.set_index('feature_label')
    assert table.loc['same', 'correlation'] == pytest.approx(1)
    assert table.loc['opposite', 'correlation'] == pytest.approx(-1)
    expected = (spearmanr if method == 'spearman' else pearsonr)([0, 2, 3, 5], [1, 5, 2, 3]).statistic
    assert table.loc['missing', 'correlation'] == pytest.approx(expected)
    assert table.loc['missing', 'paired_columns'] == 4
    assert table.loc['constant', 'status'] == 'constant_on_overlap'
    assert table.loc['short', 'status'] == 'insufficient_pairs'
    assert record['input_matrix']['shape'] == [6, 6]
    assert record['selected_sample_metadata'] == matrix.samples
    assert record['biological_quality_score'] is None and not record['p_values_computed']
    assert not record['cross_study_pooling'] and not record['gene_aggregation']
    assert record['source_record']['units'] == 'Original units'


def test_duplicate_labels_need_original_row_and_repeated_samples_are_refused(matrix):
    matrix.labels.loc[11] = 'seed'
    with pytest.raises(ValueError, match='2 rows'):
        C.correlate(matrix, 'seed', list('abcdef'))
    result, record = C.correlate(matrix, 'seed', list('abcdef'), seed_row=10)
    assert 11 in set(result.source_row) and record['seed_row'] == 10
    matrix.samples[0]['sample_accession'] = matrix.samples[1]['sample_accession'] = 'SAME-SAMPLE'
    with pytest.raises(ValueError, match='repeated sequencing'):
        C.correlate(matrix, 'seed', list('abcdef'), seed_row=10)


def test_selection_transform_cancellation_and_input_immutability(matrix):
    original = matrix.values.copy(deep=True)
    with pytest.raises(ValueError, match='Select at least'):
        C.correlate(matrix, 'seed', list('abc'))
    with pytest.raises(ValueError, match='distinct columns'):
        C.correlate(matrix, 'seed', list('aabcde'))
    with pytest.raises(InterruptedError):
        C.correlate(matrix, 'seed', list('abcdef'), cancelled=lambda: True)
    table, record = C.correlate(matrix, 'seed', list('abcdef'), transform='log2p1', method='pearson')
    expected = pearsonr(np.log2(np.arange(6) + 1), np.log2(np.arange(6) * 2 + 1)).statistic
    assert table.set_index('feature_label').loc['same', 'correlation'] == pytest.approx(expected)
    assert record['transform'] == 'log2p1'
    pd.testing.assert_frame_equal(original, matrix.values)
    matrix.values.iloc[2, 2] = -1
    with pytest.raises(ValueError, match='nonnegative'):
        C.correlate(matrix, 'seed', list('abcdef'), transform='log2p1')


@pytest.mark.parametrize('scale', [1e-250, 1e300])
def test_pearson_does_not_lose_extreme_but_finite_source_scales(matrix, scale):
    matrix.values *= scale
    result, _ = C.correlate(matrix, 'seed', list('abcdef'), method='pearson')
    assert result.set_index('feature_label').loc['same', 'correlation'] == pytest.approx(1.)
    assert result.set_index('feature_label').loc['opposite', 'correlation'] == pytest.approx(-1.)


@pytest.fixture
def sources(tmp_path):
    """One on-disk RNA matrix and one independently bound original protein table."""
    columns = ['run' + str(i) for i in range(6)]
    values = pd.DataFrame([np.arange(6), np.arange(6)*2, np.arange(5, -1, -1)],
                          index=['geneA', 'geneB', 'geneC'], columns=columns, dtype=float)
    rna = tmp_path / 'rna'
    target = rna / 'matrices' / 'STUDY'
    target.mkdir(parents=True)
    pd.concat({'count': values, 'tpm': values}, axis=1).to_parquet(target / 'salmon_ME49_host.parquet')
    pd.DataFrame(dict(run_accession=columns, sample_accession=['sample' + str(i) for i in range(6)],
                      condition=['WT']*6)).to_csv(target / 'salmon_ME49_host_samples.tsv', sep='\t', index=False)
    root = tmp_path / 'deposited'
    pack = root / 'protein'
    pack.mkdir(parents=True)
    frame = values.reset_index(names='feature_label')
    frame.insert(0, 'source_row', [0, 1, 2])
    frame.to_parquet(pack / 'values.parquet', index=False)
    tokens = frame.copy()
    for c in columns:
        tokens[c] = tokens[c].astype(str)
    tokens.to_parquet(pack / 'source_tokens.parquet', index=False)
    (pack / 'samples.json').write_text(json.dumps(dict(study='PROTEIN', source_analysis='LFQ table',
        samples=[dict(column=c, source_sample='Submitted ' + c, host_cell='HFF') for c in columns])))
    (pack / 'source_metadata.json').write_text(json.dumps(dict(original='source context')))
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    info = dict(schema_version=1, study='PROTEIN', source_analysis='LFQ table',
        assay='deposited protein-group quantities', quantity='Original submitted LFQ intensity',
        units='Source LFQ intensity scale; absolute units unreviewed', quantity_domain='nonnegative',
        features=3, sample_columns=columns, source_samples=['Submitted ' + c for c in columns],
        source_url='https://example.org/source', input_sha256='b'*64,
        comparison_scope='Original source only', outputs={p.name: digest(p) for p in pack.iterdir()})
    (pack / 'manifest.json').write_text(json.dumps(info))
    (pack / 'validation.json').write_text(json.dumps(dict(status='source_cells_reconciled',
        manifest_sha256=digest(pack / 'manifest.json'), source_sha256='b'*64,
        verified_rows=3, verified_cells=18, literal_tokens_verified=True)))
    (root / 'catalogue.json').write_text(json.dumps(dict(packs=['protein'])))
    return rna, root


def test_sources_keep_assay_units_identifiers_metadata_and_corrupt_output_refusal(sources):
    rna, root = sources
    catalogue, errors = C.catalogue(rna, root)
    assert not errors and {s['level'] for s in catalogue} == {'RNA', 'Protein'}
    for source in catalogue:
        loaded = C.load(source, rna)
        table, record = C.correlate(loaded, 'geneA', list(loaded.values.columns))
        assert set(table.feature_label) == {'geneB', 'geneC'}
        assert record['selected_sample_metadata'][0]['column'] == 'run0'
        if source['level'] == 'Protein':
            assert record['source_record']['manifest']['input_sha256'] == 'b'*64
            assert 'LFQ' in record['source_record']['units']
            assert loaded.samples[0]['host_cell'] == 'HFF'
    (root / 'protein' / 'values.parquet').write_bytes(b'corrupted')
    with pytest.raises(ValueError, match='checksum'):
        C.load(next(s for s in catalogue if s['level'] == 'Protein'))


@pytest.mark.parametrize('level', ['RNA', 'Protein'])
def test_qt_load_calculate_plot_export_and_invalidate(qtbot, sources, level, tmp_path):
    from parascape.coexpression_panel import CoexpressionPanel
    rna, root = sources
    panel = CoexpressionPanel(rna, root)
    qtbot.addWidget(panel)
    panel.level.setCurrentText(level)
    panel.load_selected()
    qtbot.waitUntil(lambda: panel.matrix is not None and panel.worker is None, timeout=10000)
    assert not panel.selected_columns()  # Explicit review, no inferred sample roles.
    panel._select_all(True)
    panel.seed.setText('geneA')
    panel.calculate()
    qtbot.waitUntil(lambda: panel.result is not None and panel.worker is None, timeout=10000)
    assert panel.table.rowCount() == 2 and len(panel.figure.axes) == 3
    path = tmp_path / (level + '.csv')
    assert panel.save(str(path)) == 2
    record = json.loads(path.with_suffix('.csv.json').read_text())
    assert record['method'] == 'spearman' and record['selected_columns'] == ['run' + str(i) for i in range(6)]
    assert record['export_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record['biological_quality_score'] is None
    panel.seed.setText('geneB')
    assert panel.result is None and not panel.export_button.isEnabled()
    assert panel.save(str(tmp_path / 'stale.csv')) == 0
    panel.close()


def test_qt_late_background_result_cannot_replace_a_changed_source(qtbot, sources, matrix, monkeypatch):
    import threading
    from parascape.coexpression_panel import CoexpressionPanel
    entered, release = threading.Event(), threading.Event()

    def delayed(*args):
        entered.set()
        if not release.wait(5):
            raise ValueError('Test operation was not released')
        return matrix

    monkeypatch.setattr(C, 'load', delayed)
    panel = CoexpressionPanel(*sources)
    qtbot.addWidget(panel)
    panel.load_selected()
    try:
        qtbot.waitUntil(entered.is_set)
        panel.level.setCurrentText('Protein')
    finally:
        release.set()
    qtbot.waitUntil(lambda: panel.worker is None)
    assert panel.matrix is None and panel.result is None and not panel.export_button.isEnabled()
    panel.close()
