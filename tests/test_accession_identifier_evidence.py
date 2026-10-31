"""Accession components retain whole source rows, distinct evidence classes and explicit unavailable mappings."""
import hashlib
import json

import pandas as pd
import pytest

from parascape import deposited_values as D
from parascape import protein_identifier_evidence as I
from parascape.deposited_values_panel import DepositedValuesPanel
from test_deposited_values import pack  # noqa: F401


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def accession_evidence(pack, tmp_path, monkeypatch):
    labels = ['001;001', '001;001', 'CaseID']
    for name in ('values.parquet', 'source_tokens.parquet'):
        frame = pd.read_parquet(pack / name)
        frame['feature_label'] = labels
        frame.to_parquet(pack / name, index=False)
    metadata = dict(original_feature_column='PG.ProteinAccessions', row_layout=dict(layout='original_rowwise_quantities'),
        row_annotations=[dict(source_row=i, feature_label=label, source_fields={'PG.ProteinAccessions': label,
            'R.FileName': 'Sample1', 'R.Condition': 'Condition ' + str(i), 'R.Replicate': '001'}) for i, label in enumerate(labels)])
    (pack / 'source_metadata.json').write_text(json.dumps(metadata))
    info = json.loads((pack / 'manifest.json').read_text())
    info['assay'] = 'deposited protein-accession quantities'
    info['outputs'] = {name: digest(pack / name) for name in info['outputs']}
    (pack / 'manifest.json').write_text(json.dumps(info))
    receipt = json.loads((pack / 'validation.json').read_text())
    receipt['manifest_sha256'] = digest(pack / 'manifest.json')
    (pack / 'validation.json').write_text(json.dumps(receipt))
    root = tmp_path / 'audit/accession_identifier_evidence'
    folder = root / ('2' * 64)
    folder.mkdir(parents=True)
    rows = [dict(source_row=i, feature_label=label, query_namespace='original_accession_cell_component',
        query_literal=component, protein_component_index=position, protein_component_literal=component,
        source_field='PG.ProteinAccessions', source_field_literal=label) for i, label in enumerate(labels)
        for position, component in enumerate(label.split(';'))]
    name = digest(pack / 'manifest.json') + '.parquet'
    pd.DataFrame(rows).to_parquet(folder / name, index=False)
    (folder / 'index.json').write_text(json.dumps(dict(schema_version=1, scope='Whole original accession rows; gene equivalence unreviewed',
        measurements=[dict(measurement_manifest_sha256=digest(pack / 'manifest.json'), input_sha256=info['input_sha256'],
            source_url=info['source_url'], features=3, assay=info['assay'], source_feature_field='PG.ProteinAccessions',
            file=name, namespaces={'original_accession_cell_component': 5})])))
    (folder / 'validation.json').write_text(json.dumps(dict(status='all_identifier_evidence_reconciled', measurement_packs=1,
        outputs={path.name: digest(path) for path in folder.iterdir()})))
    (root / 'collection.json').write_text(json.dumps(dict(schema_version=1, directory=folder.name,
        validation_sha256=digest(folder / 'validation.json'))))
    monkeypatch.setenv('PARASCAPE_AUDIT', str(tmp_path / 'audit'))
    monkeypatch.setenv('PARASCAPE_ACCESSION_IDENTIFIERS', str(root))
    return pack, root, folder / name


def test_repeated_cell_components_do_not_duplicate_or_merge_original_measurements(accession_evidence):
    pack, root, _ = accession_evidence
    evidence = I.records(pack, '001', namespace='original_accession_cell_component', root=root)
    assert len(evidence) == 4 and set(evidence.protein_component_index) == {0, 1}
    result = I.profiles(pack, '001', namespace='original_accession_cell_component', root=root)
    raw = D.profiles(pack, '001;001')
    pd.testing.assert_frame_equal(result[raw.columns], raw, check_exact=True)
    assert len(result) == 4 and result.source_token.tolist() == ['-0.000', '1.2500', '0', '-2.500']
    assert all(len(json.loads(value)) == 2 for value in result.identifier_evidence)
    assert I.profiles(pack, '1', namespace='original_accession_cell_component').empty
    assert I.profiles(pack, 'caseid', namespace='original_accession_cell_component').empty
    with pytest.raises(ValueError, match='does not apply'):
        I.records(pack, '001', namespace='original_protein_group_component')


@pytest.mark.parametrize('change', ['component_index', 'component_literal', 'original_field', 'full_cell'])
def test_rehashed_component_rows_must_still_match_original_cell_components(accession_evidence, change):
    pack, root, path = accession_evidence
    frame = pd.read_parquet(path)
    key, value = {'component_index': ('protein_component_index', 8), 'component_literal': ('protein_component_literal', '1'),
                  'original_field': ('source_field', 'Gene names'), 'full_cell': ('source_field_literal', '001')}[change]
    frame.loc[0, key] = value
    frame.to_parquet(path, index=False)
    folder = path.parent
    receipt = json.loads((folder / 'validation.json').read_text())
    receipt['outputs'][path.name] = digest(path)
    (folder / 'validation.json').write_text(json.dumps(receipt))
    pointer = json.loads((root / 'collection.json').read_text())
    pointer['validation_sha256'] = digest(folder / 'validation.json')
    (root / 'collection.json').write_text(json.dumps(pointer))
    with pytest.raises(ValueError, match='component identity'):
        I.records(pack, '001', namespace='original_accession_cell_component')


def test_browser_assay_modes_empty_mapping_and_row_context_remain_explicit(qtbot, accession_evidence, tmp_path):
    pack, _, _ = accession_evidence
    panel = DepositedValuesPanel(pack.parent)
    qtbot.addWidget(panel)
    group = panel.identifier_kind.findData('original_protein_group_component')
    accession = panel.identifier_kind.findData('original_accession_cell_component')
    assert not panel.identifier_kind.model().item(group).isEnabled()
    assert panel.identifier_kind.model().item(accession).isEnabled()
    panel.identifier_kind.setCurrentIndex(accession)
    panel.features.setPlainText('001')
    panel.show_values()
    assert len(panel.result) == 4 and 'row_design_metadata' in panel.result
    assert 'Original accession-row quantities retained' in panel.summary.text()
    panel.table.setCurrentCell(0, 0)
    panel.show_measurement_context()
    assert panel.measurement_dialog is not None
    assert panel.save(str(tmp_path / 'accession.csv')) == 4
    panel.identifier_kind.setCurrentIndex(panel.identifier_kind.findData('reference_declared_gene_identifier'))
    assert panel.result is None and panel.measurement_dialog is None
    panel.features.setPlainText('G001')
    panel.show_values()
    assert panel.result is not None and panel.result.empty
    assert 'no indexed declarations' in panel.summary.text() and 'not evidence of absent expression' in panel.summary.text()
    assert panel.treatment_record['indexed_identifier_class_records'] == 0
    assert panel.save_treatment(str(tmp_path / 'unavailable-mapping.json')) == 1
    assert panel.save(str(tmp_path / 'absent.csv')) == 0


def test_assay_roots_and_allowed_namespaces_do_not_mix_accessions_and_groups(monkeypatch, tmp_path):
    monkeypatch.setenv('PARASCAPE_PROTEIN_IDENTIFIERS', str(tmp_path / 'groups'))
    monkeypatch.setenv('PARASCAPE_ACCESSION_IDENTIFIERS', str(tmp_path / 'accessions'))
    assert I.data_root() == tmp_path / 'groups'
    assert I.data_root('deposited protein-accession quantities') == tmp_path / 'accessions'
    assert 'original_accession_cell_component' not in I.namespaces_for_assay('deposited protein-group quantities')
    assert I.namespaces_for_assay('unreviewed assay') == set()
