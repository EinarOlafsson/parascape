"""Original-value packs preserve source probes, signs, missingness and sample identities."""
import hashlib
import json

import pandas as pd
import pytest

from parascape import source_matrices as S


def digest(path):
    """Compute fixture checksums independently of the source-matrix reader."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def pack(tmp_path):
    source = tmp_path / 'study-original'
    source.mkdir()
    pd.DataFrame({'probe_id': ['001', 'CaseID'], 'GSM1': [-1.5, 0.], 'GSM2': [float('nan'), 2.5]}).to_parquet(source / 'values.parquet', index=False)
    (source / 'samples.json').write_text(json.dumps({'study': 'GSE1', 'samples': [
        {'geo_accession': 'GSM1', 'title': 'Original sample one', 'channels': [{}, {}]},
        {'geo_accession': 'GSM2', 'title': 'Original sample two', 'channels': [{}]}]}))
    (source / 'original_scale_and_platform.json').write_text(json.dumps({'!Series_platform_id': [['GPL1']]}))
    info = dict(schema_version=1, study='GSE1', catalogue_readout='synthetic assay', sample_order=['GSM1', 'GSM2'],
                features=2, quantity='Original submitted VALUE', units='source scale not inferred',
                source_url='https://example.org/original.gz', input_sha256='a' * 64,
                feature_identity='exact source ID; no remapping', comparison_scope='descriptive source values only',
                outputs={p.name: digest(p) for p in source.iterdir()})
    (source / 'manifest.json').write_text(json.dumps(info))
    (source / 'validation.json').write_text(json.dumps(dict(manifest_sha256=digest(source / 'manifest.json'),
        original_source_sha256=info['input_sha256'], verified_rows=2, verified_numeric_or_missing_cells=4)))
    return source


def test_original_probe_ids_signed_values_missingness_and_source_metadata(pack):
    frame = S.profiles(pack, ['001', 'CaseID'])
    assert len(frame) == 4 and frame['value'].isna().sum() == 1
    assert set(frame['value'].dropna()) == {-1.5, 0., 2.5}
    assert set(frame.probe_id) == {'001', 'CaseID'}
    assert S.profiles(pack, 'caseid').empty
    assert len(json.loads(frame.iloc[0].sample_metadata)['channels']) == 2
    assert frame.matrix_sha256.notna().all() and frame.source_url.notna().all()


def test_changed_outputs_and_incomplete_source_validation_are_rejected(pack):
    (pack / 'values.parquet').write_bytes(b'changed')
    with pytest.raises(ValueError, match='checksum'):
        S.profiles(pack, '001')
    validation = json.loads((pack / 'validation.json').read_text())
    validation['verified_numeric_or_missing_cells'] = 3
    (pack / 'validation.json').write_text(json.dumps(validation))
    with pytest.raises(ValueError, match='coverage differs'):
        S.manifest(pack)


def test_panel_exports_and_clears_old_selection(qtbot, pack, tmp_path):
    from parascape.source_matrices_panel import SourceMatricesPanel
    panel = SourceMatricesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.probes.setPlainText('001\nabsent')
    panel.compare()
    assert panel.table.rowCount() == 2 and 'absent' in panel.summary.text()
    path = tmp_path / 'values.csv'
    assert panel.save(str(path)) == 2
    assert 'comparison_scope' in pd.read_csv(path)
    panel.refresh()
    assert panel.save(str(path)) == 0 and panel.table.rowCount() == 0


@pytest.fixture
def identifier_index(pack, tmp_path, monkeypatch):
    root = tmp_path / 'platform_identifiers'
    folder = root / 'GPL1-fixture'
    folder.mkdir(parents=True)
    pd.DataFrame([dict(probe_id=probe, identifier=identifier, source_column='Gene Symbol',
        source_cell='GeneA /// geneB' if probe == 'CaseID' else 'GeneA',
        multivalue_source_field=probe == 'CaseID') for probe, identifier in
        [('001', 'GeneA'), ('CaseID', 'GeneA'), ('CaseID', 'geneB')]]).to_parquet(folder / 'links.parquet', index=False)
    info = dict(platform='GPL1', source_url='https://example.org/annotation', source_sha256='c'*64,
                identifier_link_records=3, observed_matrix_pack_manifests={pack.name: digest(pack / 'manifest.json')},
                outputs={'links.parquet': digest(folder / 'links.parquet')})
    (folder / 'manifest.json').write_text(json.dumps(info))
    receipt = dict(manifest_sha256=digest(folder / 'manifest.json'), verified_links=3, source_sha256='c'*64)
    (folder / 'independent_validation.json').write_text(json.dumps(receipt))
    (root / 'independent_validation_coverage.json').write_text(json.dumps([dict(platform='GPL1',
        status='source_identifier_links_reconciled', path='/original/NAS/GPL1-fixture',
        independent_validation_sha256=digest(folder / 'independent_validation.json'))]))
    monkeypatch.setenv('PARASCAPE_PLATFORM_IDENTIFIERS', str(root))
    return root


def test_annotation_queries_preserve_probes_signed_values_missingness_and_ambiguity(pack, identifier_index):
    frame = S.annotated_profiles(pack, 'GeneA', identifier_index)
    assert len(frame) == 4 and frame.probe_id.nunique() == 2
    assert frame['value'].isna().sum() == 1 and set(frame['value'].dropna()) == {-1.5, 0., 2.5}
    assert set(frame.query_identifier) == {'GeneA'}
    assert any(json.loads(value)[0]['multivalue_source_field'] for value in frame.annotation_evidence)
    assert S.annotated_profiles(pack, 'GENEA', identifier_index).empty
    path = pack / 'manifest.json'
    info = json.loads(path.read_text()); info['units'] = 'changed context'; path.write_text(json.dumps(info))
    receipt = json.loads((pack / 'validation.json').read_text()); receipt['manifest_sha256'] = digest(path)
    (pack / 'validation.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match='not prepared for this matrix'):
        S.annotated_profiles(pack, 'GeneA', identifier_index)


def test_panel_annotation_lookup_retains_original_probe_rows(qtbot, pack, identifier_index):
    from parascape.source_matrices_panel import SourceMatricesPanel
    panel = SourceMatricesPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.probes.setPlainText('GeneA')
    panel.annotation.setChecked(True)
    panel.compare()
    assert panel.table.rowCount() == 4 and 'not matched' not in panel.summary.text()
    assert {panel.table.item(row, 1).text() for row in range(4)} == {'001', 'CaseID'}
