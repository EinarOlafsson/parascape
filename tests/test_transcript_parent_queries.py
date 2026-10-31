"""Reported parent identifiers select original transcript rows without aggregation."""
import json

import numpy as np
import pandas as pd
import pytest

from parascape import deposited_values as D
from parascape.deposited_values_panel import DepositedValuesPanel
from test_deposited_values import pack, digest


def rebind(folder, metadata=None, info=None):
    if metadata is not None:
        (folder/'source_metadata.json').write_text(json.dumps(metadata))
    info = info or json.loads((folder/'manifest.json').read_text())
    info['outputs']['source_metadata.json'] = digest(folder/'source_metadata.json')
    (folder/'manifest.json').write_text(json.dumps(info))
    receipt = json.loads((folder/'validation.json').read_text())
    receipt['manifest_sha256'] = digest(folder/'manifest.json')
    (folder/'validation.json').write_text(json.dumps(receipt))


@pytest.fixture
def transcript_pack(pack):
    metadata = dict(original_feature_column='transcript_id', original_annotation_columns=['gene_id','length'],
        row_annotations=[dict(source_row=i, feature_label=label, source_fields=dict(gene_id=parent,length='001.00'))
                         for i,(label,parent) in enumerate([('001','GeneA'),('001','GeneB'),('CaseID','GeneA')])])
    info = json.loads((pack/'manifest.json').read_text()); info['assay']='deposited original transcript-table quantities'
    rebind(pack,metadata,info)
    return pack


def test_parent_matches_every_original_row_and_keeps_shared_transcript_labels_distinct(transcript_pack):
    source = transcript_pack
    before = {name:digest(source/name) for name in D.manifest(source)['outputs']}
    result = D.profiles(source,['GeneA','GeneA','missing'],match='reported_gene_id')
    assert result.source_row.tolist()==[0,2,0,2]
    assert result.feature_label.tolist()==['001','CaseID','001','CaseID']
    assert result.source_token.tolist()==['-0.000','','0','2.5']
    assert np.signbit(result.iloc[0].value) and pd.isna(result.iloc[1].value)
    assert not (result.source_row==1).any()
    baseline=D.profiles(source,['001','CaseID']);baseline=baseline[baseline.source_row.isin([0,2])]
    assert result.to_csv(index=False)==baseline.to_csv(index=False)
    record=result.attrs['identifier_query']
    assert record['selected_source_rows']==[0,2] and record['unmatched_literals']==['missing']
    assert not record['gene_aggregation_performed'] and not record['canonical_gene_equivalence_established']
    assert before=={name:digest(source/name) for name in before}
    assert D.profiles(source,'genea',match='reported_gene_id').empty
    assert D.profiles(source,'GeneA-t26_1',match='reported_gene_id').empty
    assert 'identifier_query' not in D.profiles(source,'001').attrs


@pytest.mark.parametrize('damage',['missing_annotations','missing_parent','wrong_label','wrong_row','wrong_feature_column','nonstring_parent'])
def test_internally_rebound_but_structurally_wrong_source_relationships_are_refused(transcript_pack,damage):
    source=transcript_pack;metadata=json.loads((source/'source_metadata.json').read_text())
    if damage=='missing_annotations':metadata['row_annotations'].pop()
    elif damage=='missing_parent':metadata['row_annotations'][1]['source_fields'].pop('gene_id')
    elif damage=='wrong_label':metadata['row_annotations'][1]['feature_label']='invented'
    elif damage=='wrong_row':metadata['row_annotations'][1]['source_row']=True
    elif damage=='wrong_feature_column':metadata['original_feature_column']='gene_id'
    else:metadata['row_annotations'][1]['source_fields']['gene_id']=123
    rebind(source,metadata)
    with pytest.raises(ValueError):D.profiles(source,'GeneA',match='reported_gene_id')


def test_parent_query_requires_supported_transcript_source(pack):
    with pytest.raises(ValueError,match='transcript-table'):D.profiles(pack,'001',match='reported_gene_id')
    with pytest.raises(ValueError,match='matching mode'):D.profiles(pack,'001',match='guessed_gene')


def test_parent_selector_export_context_and_stale_state(qtbot,transcript_pack,tmp_path):
    panel=DepositedValuesPanel(transcript_pack.parent);qtbot.addWidget(panel)
    index=panel.identifier_kind.findData('reported_gene_id')
    assert panel.identifier_kind.model().item(index).isEnabled()
    panel.identifier_kind.setCurrentIndex(index);panel.features.setPlainText('GeneA');panel.show_values()
    assert panel.result is not None and panel.result.source_row.tolist()==[0,2,0,2]
    assert panel.treatment_record['reported_parent_gene_query']['selected_source_rows']==[0,2]
    assert 'No gene-level sum' in panel.summary.text()
    output=tmp_path/'parent.csv';assert panel.save(str(output))==4
    records=[json.loads(p.read_text()) for p in tmp_path.glob('parent.csv.treatment-*.json')]
    assert any(r['reported_parent_gene_query']['original_annotation_column']=='gene_id' for r in records)
    panel.features.setPlainText('GeneB');assert panel.result is None and not panel.export_button.isEnabled()
    panel.show_values();assert panel.result.source_row.tolist()==[1,1]
    panel.identifier_kind.setCurrentIndex(0);assert panel.result is None
    panel.show_values();assert panel.result.empty
