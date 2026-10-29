"""Deposited values preserve source identity, units, zeroes and annotation separation."""
import hashlib
import json

import pandas as pd
import pytest

from parascape import published_expression as P


def digest(path):
    """Compute fixture integrity independently from the reader."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_manifest(pack, info):
    """Write a synthetic validation receipt for an explicitly synthetic test pack."""
    (pack / 'manifest.json').write_text(json.dumps(info))
    (pack / 'independent_validation.json').write_text(json.dumps({'manifest_sha256': digest(pack / 'manifest.json')}))


@pytest.fixture
def source_pack(tmp_path):
    """Two annotation variants share two samples and retain different exact gene IDs."""
    pack = tmp_path / 'source_pack'
    pack.mkdir()
    fields = ['expected_count', 'TPM', 'FPKM', 'effective_length', 'length']
    variants, outputs, inputs = {}, {}, []
    for ref, gene in [('annotation_one', 'exactGene'), ('annotation_two', 'MSTRG.1')]:
        columns = pd.MultiIndex.from_product([fields, ['GSM1', 'GSM2']])
        frame = pd.DataFrame([[0., 2.5, 0., 30.25, 0., 10.125, 10., 20., 30., 40.]], index=[gene], columns=columns)
        frame.to_parquet(pack / (ref + '.parquet'))
        outputs[ref + '.parquet'] = digest(pack / (ref + '.parquet'))
        variants[ref] = {'genes': 1, 'samples': ['GSM1', 'GSM2'], 'fields': fields}
        inputs.extend({'reference': ref, 'sample': s, 'url': 'https://example.org/' + s,
                       'sha256': 'synthetic-input-hash'} for s in ['GSM1', 'GSM2'])
    samples = {'primary_archive_sha256': 'synthetic-archive', 'samples': [
        {'geo_accession': s, 'title': title, 'relations': [], 'pairing_status': 'unknown'}
        for s, title in [('GSM1', 'original one'), ('GSM2', 'original two')]]}
    (pack / 'samples.json').write_text(json.dumps(samples))
    outputs['samples.json'] = digest(pack / 'samples.json')
    info = {'schema_version': 1, 'study': 'synthetic_study', 'software': {'RSEM': 'synthetic-version'},
            'variants': variants, 'inputs': inputs, 'outputs': outputs,
            'primary_archive_sha256': 'synthetic-archive', 'normalization': 'None; source units retained'}
    write_manifest(pack, info)
    return pack, info


def test_original_values_units_samples_and_exact_annotation_ids(source_pack):
    pack, _ = source_pack
    frame = P.profiles(pack, 'exactGene missing', 'annotation_one')
    assert frame.gene.tolist() == ['exactGene', 'exactGene']
    assert frame.TPM.tolist() == [0., 30.25]
    assert frame.FPKM.tolist() == [0., 10.125]
    assert frame.expected_count.tolist() == [0., 2.5]
    assert frame['sample'].nunique() == 2
    assert set(frame.quantifier) == {'RSEM'}
    assert 'unknown' in frame.sample_metadata.iloc[0]
    assert P.profiles(pack, 'exactgene', 'annotation_one').empty
    assert P.profiles(pack, 'exactGene', 'annotation_two').empty
    assert len(P.profiles(pack, 'MSTRG.1', 'annotation_two')) == 2


def test_manifest_and_matrix_tampering_are_rejected(source_pack):
    pack, info = source_pack
    (pack / 'annotation_one.parquet').write_bytes(b'corrupted')
    with pytest.raises(ValueError, match='checksum'):
        P.profiles(pack, 'exactGene', 'annotation_one')
    packs, errors = P.catalogue(pack.parent)
    assert not packs and errors
    info['study'] = 'changed'
    (pack / 'manifest.json').write_text(json.dumps(info))
    with pytest.raises(ValueError, match='validation receipt'):
        P.manifest(pack)


def test_annotations_are_not_counted_as_extra_samples(source_pack):
    pack, info = source_pack
    packs, errors = P.catalogue(pack.parent)
    assert not errors and packs[0]['samples'] == 2
    info['variants']['annotation_one']['samples'].append('GSM1')
    write_manifest(pack, info)
    with pytest.raises(ValueError, match='Duplicated samples'):
        P.manifest(pack)


def test_undeclared_output_and_path_escape_are_rejected(source_pack):
    pack, info = source_pack
    del info['outputs']['samples.json']
    write_manifest(pack, info)
    with pytest.raises(ValueError, match='does not declare'):
        P.manifest(pack)
    with pytest.raises(ValueError, match='outside'):
        P._file(pack, '../other-data')


def test_panel_and_csv_preserve_source_metadata(qtbot, source_pack, tmp_path):
    from parascape.published_expression_panel import PublishedExpressionPanel
    pack, _ = source_pack
    panel = PublishedExpressionPanel(pack.parent)
    qtbot.addWidget(panel)
    panel.genes.setText('exactGene missing')
    panel.compare()
    assert panel.table.rowCount() == 2
    assert 'missing' in panel.summary.text()
    output = tmp_path / 'export.csv'
    assert panel.save(str(output)) == 2
    values = pd.read_csv(output)
    assert set(values.reference) == {'annotation_one'}
    assert 'sample_metadata' in values and 'source_sha256' in values
    panel.reference.setCurrentText('annotation_two')
    assert panel.table.rowCount() == 0 and panel.save(str(output)) == 0
