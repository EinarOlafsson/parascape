"""Identifier matches retain original group values and ambiguity without creating gene-level measurements."""
import hashlib
import json

import pandas as pd
import pytest

from parascape import protein_identifier_evidence as I
from test_deposited_values import pack  # noqa: F401


def digest(path):
    """Pin synthetic identifier evidence independently of the production reader."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def evidence(pack, tmp_path, monkeypatch):
    """Only one of two identical group labels declares a queried gene name; evidence stays row-bound."""
    info=json.loads((pack/'manifest.json').read_text());info['assay']='deposited protein-group quantities'
    (pack/'manifest.json').write_text(json.dumps(info))
    old=json.loads((pack/'validation.json').read_text());old['manifest_sha256']=digest(pack/'manifest.json')
    (pack/'validation.json').write_text(json.dumps(old))
    root=tmp_path/'audit/protein_identifier_evidence';folder=root/('1'*64);folder.mkdir(parents=True)
    rows=[dict(source_row=0,feature_label='001',query_namespace='submitted_gene_label',query_literal='G001',
               source_field='Gene names',source_field_literal='G001;G001',label_component_index=position)
          for position in [0,1]]
    rows += [dict(source_row=position,feature_label='001',query_namespace='original_protein_group_component',query_literal='001',
                  source_field='',source_field_literal='',label_component_index=-1) for position in [0,1]]
    name=digest(pack/'manifest.json')+'.parquet'
    pd.DataFrame(rows).to_parquet(folder/name,index=False)
    (folder/'index.json').write_text(json.dumps(dict(schema_version=1,scope='Whole original groups; unique-gene abundance and current identity remain unreviewed',
        measurements=[dict(measurement_manifest_sha256=digest(pack/'manifest.json'),input_sha256=info['input_sha256'],
            source_url=info['source_url'],features=info['features'],file=name)])))
    receipt=folder/'validation.json';receipt.write_text(json.dumps(dict(status='all_identifier_evidence_reconciled',
        measurement_packs=1,outputs={p.name:digest(p) for p in folder.iterdir()})))
    (root/'collection.json').write_text(json.dumps(dict(schema_version=1,directory=folder.name,validation_sha256=digest(receipt))))
    monkeypatch.setenv('PARASCAPE_PROTEIN_IDENTIFIERS',str(root))
    monkeypatch.setenv('PARASCAPE_AUDIT',str(tmp_path/'audit'))
    return pack,root,folder/name


def test_repeated_evidence_does_not_duplicate_measurements_or_match_other_identical_group_rows(evidence):
    pack,root,_=evidence
    records=I.records(pack,'G001',namespace='submitted_gene_label',root=root)
    assert len(records)==2 and records.label_component_index.tolist()==[0,1]
    values=I.profiles(pack,'G001',namespace='submitted_gene_label',root=root)
    assert len(values)==2 and set(values.source_row)=={0}
    assert values.source_token.tolist()==['-0.000','0']
    assert all(len(json.loads(value))==2 for value in values.identifier_evidence)
    raw=I.profiles(pack,'001',namespace='original_protein_group_component',root=root)
    assert len(raw)==4 and set(raw.source_row)=={0,1}
    assert I.profiles(pack,'g001',namespace='submitted_gene_label',root=root).empty


def test_changed_evidence_and_stale_measurement_manifests_are_rejected(evidence):
    pack,root,path=evidence
    original=path.read_bytes();path.write_bytes(b'changed evidence')
    with pytest.raises(ValueError,match='output checksum'):
        I.records(pack,'G001',namespace='submitted_gene_label',root=root)
    path.write_bytes(original)
    info=json.loads((pack/'manifest.json').read_text());info['units']='Changed measurement context'
    (pack/'manifest.json').write_text(json.dumps(info))
    receipt=json.loads((pack/'validation.json').read_text());receipt['manifest_sha256']=digest(pack/'manifest.json')
    (pack/'validation.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='exact measurement manifest'):
        I.records(pack,'G001',namespace='submitted_gene_label',root=root)


def test_query_namespaces_are_distinct_and_results_clear_when_selection_changes(qtbot,evidence):
    from parascape.deposited_values_panel import DepositedValuesPanel
    pack,_,_=evidence
    panel=DepositedValuesPanel(pack.parent);qtbot.addWidget(panel)
    index=panel.identifier_kind.findData('submitted_gene_label');panel.identifier_kind.setCurrentIndex(index)
    panel.features.setPlainText('G001');panel.show_values()
    assert panel.table.rowCount()==2
    assert 'Labels not matched' not in panel.summary.text()
    assert 'unique-gene abundance' in panel.summary.text()
    panel.identifier_kind.setCurrentIndex(0)
    assert panel.result is None and panel.table.rowCount()==0
    panel.show_values()
    assert panel.table.rowCount()==0
