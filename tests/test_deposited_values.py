"""Original deposited values retain source rows, literal precision and sample context."""
import hashlib
import json
import shutil

import numpy as np
import pandas as pd
import pytest

from parascape import deposited_values as D


def digest(path):
    """Pin synthetic original-value fixture bytes independently of the runtime reader."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def pack(tmp_path):
    """Repeated labels describe distinct source rows, with signed zero and missingness."""
    folder = tmp_path / 'original-source'
    folder.mkdir()
    axes = dict(source_row=[0, 1, 2], feature_label=['001', '001', 'CaseID'])
    pd.DataFrame(dict(**axes, sample_0=[-0., 1.25, np.nan], sample_1=[0., -2.5, 2.5])).to_parquet(folder / 'values.parquet', index=False)
    pd.DataFrame(dict(**axes, sample_0=['-0.000', '1.2500', ''], sample_1=['0', '-2.500', '2.5'])).to_parquet(folder / 'source_tokens.parquet', index=False)
    (folder / 'samples.json').write_text(json.dumps(dict(study='SOURCE1', source_analysis='A1', samples=[
        dict(column='sample_0', source_sample='Original 01', submitted_factor_rows=[{'batch': 'A'}, {'batch': 'B'}]),
        dict(column='sample_1', source_sample='Original 02', submitted_factor_rows=[{'batch': 'C'}])])))
    (folder / 'source_metadata.json').write_text(json.dumps(dict(source_units='Original units')))
    info = dict(schema_version=1, study='SOURCE1', source_analysis='A1', assay='synthetic assay',
        quantity='Original quantity', units='Original units', quantity_domain='unspecified_real',
        features=3, sample_columns=['sample_0', 'sample_1'], source_samples=['Original 01', 'Original 02'],
        source_url='https://example.org/original', input_sha256='a' * 64,
        comparison_scope='Original source only; no row collapse or inferred independence',
        outputs={p.name: digest(p) for p in folder.iterdir()})
    (folder / 'manifest.json').write_text(json.dumps(info))
    (folder / 'validation.json').write_text(json.dumps(dict(status='source_cells_reconciled',
        manifest_sha256=digest(folder / 'manifest.json'), source_sha256='a'*64,
        verified_rows=3, verified_cells=6, literal_tokens_verified=True)))
    (tmp_path / 'catalogue.json').write_text(json.dumps(dict(packs=[folder.name])))
    return folder


def test_original_rows_precision_signed_zero_and_ambiguous_labels_are_not_collapsed(pack):
    frame = D.profiles(pack, '001')
    assert len(frame) == 4 and set(frame.source_row) == {0, 1}
    assert frame.source_token.tolist() == ['-0.000', '1.2500', '0', '-2.500']
    assert np.signbit(frame.iloc[0]['value'])
    assert len(json.loads(frame.iloc[0].sample_metadata)['submitted_factor_rows']) == 2
    assert D.profiles(pack, 'caseid').empty
    assert D.profiles(pack, 'CaseID')['value'].isna().sum() == 1
    assert frame.units.eq('Original units').all() and frame.matrix_sha256.eq(digest(pack / 'values.parquet')).all()


def test_full_matrix_qc_keeps_original_row_denominator_missingness_and_signed_values(pack):
    before = digest(pack / 'values.parquet')
    frame = D.technical_qc(pack)
    assert frame.features.tolist() == [3, 3]
    assert frame.missing.tolist() == [1, 0]
    assert frame.measured_zero.tolist() == [1, 1]
    assert frame.negative.tolist() == [0, 1]
    assert frame.numeric_validity_observed_pct.tolist() == [100., 100.]
    assert not frame.negative_values_invalid.any() and frame.quality_score.isna().all()
    assert digest(pack / 'values.parquet') == before


def test_changed_outputs_and_incomplete_independent_cell_coverage_are_rejected(pack):
    receipt = json.loads((pack / 'validation.json').read_text())
    receipt['verified_cells'] = 5
    (pack / 'validation.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match='coverage differs'):
        D.profiles(pack, '001')
    receipt['verified_cells'] = 6
    (pack / 'validation.json').write_text(json.dumps(receipt))
    (pack / 'source_tokens.parquet').write_bytes(b'changed source tokens')
    with pytest.raises(ValueError, match='checksum'):
        D.technical_qc(pack)


def test_catalogue_refuses_paths_outside_the_installed_collection(pack):
    path = pack.parent / 'catalogue.json'
    path.write_text(json.dumps(dict(packs=['../outside'])))
    packs, errors = D.catalogue(pack.parent)
    assert not packs and errors


def test_panel_retains_duplicate_source_rows_and_exports_qc_without_a_query(qtbot, pack, tmp_path):
    from parascape.deposited_values_panel import DepositedValuesPanel
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.features.setPlainText('001')
    panel.show_values()
    assert panel.table.rowCount() == 4
    assert panel.table.item(0, 4).text() == '-0.000'
    path = tmp_path / 'values.csv'
    assert panel.save(str(path)) == 4
    frame = pd.read_csv(path, dtype={'source_token': str}, keep_default_na=False)
    assert frame.source_token.tolist() == ['-0.000', '1.2500', '0', '-2.500']
    panel.refresh()
    assert panel.save(str(path)) == 0 and panel.table.rowCount() == 0
    panel.features.setPlainText('absent label')
    qc = tmp_path / 'qc.csv'
    assert panel.save_qc(str(qc)) == 2
    assert pd.read_csv(qc).features.tolist() == [3, 3]


def test_independent_source_namespaces_coexist_without_replacing_the_root_catalogue(pack):
    namespace = pack.parent / 'protein_groups'
    namespace.mkdir()
    copied = namespace / 'protein-source'
    shutil.copytree(pack, copied)
    (namespace / 'catalogue.json').write_text(json.dumps(dict(packs=[copied.name])))
    packs, errors = D.catalogue(pack.parent)
    assert len(packs) == 2 and not errors
    assert {row['path'] for row in packs} == {str(pack), str(copied)}
    assert len(D.profiles(copied, '001')) == 4


def test_original_row_quality_flags_are_exported_without_filtering_or_group_collapse(pack):
    path = pack / 'source_metadata.json'
    metadata = json.loads(path.read_text())
    metadata['row_annotations'] = [dict(source_row=i, feature_label=label,
        source_fields={'Reverse': '+' if i == 0 else '', 'Potential contaminant': '+' if i == 1 else ''})
        for i, label in enumerate(['001', '001', 'CaseID'])]
    path.write_text(json.dumps(metadata))
    manifest = json.loads((pack / 'manifest.json').read_text())
    manifest['outputs']['source_metadata.json'] = digest(path)
    (pack / 'manifest.json').write_text(json.dumps(manifest))
    validation = json.loads((pack / 'validation.json').read_text())
    validation['manifest_sha256'] = digest(pack / 'manifest.json')
    (pack / 'validation.json').write_text(json.dumps(validation))
    frame = D.profiles(pack, '001')
    assert len(frame) == 4
    flags = {row.source_row: json.loads(row.source_row_metadata)['source_fields'] for row in frame.itertuples()}
    assert flags[0]['Reverse'] == '+' and flags[1]['Potential contaminant'] == '+'
    assert frame.units.eq('Original units').all()
