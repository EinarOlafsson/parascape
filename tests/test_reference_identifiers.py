"""Original reference header lookup preserves exact declarations and refuses changed provenance."""
import hashlib
import json

import pandas as pd
import pytest

from parascape import reference_identifiers as R


def digest(path):
    """Pin synthetic reference fixture bytes without using the production reader."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def references(tmp_path):
    """Repeated primary identifiers and gene declarations remain distinct source header rows."""
    folder=tmp_path/('1'*64)
    folder.mkdir()
    fields=[
        [dict(key='gene',separator='=',literal_value='G001'),dict(key='transcript',separator='=',literal_value='T001')],
        [dict(key='gene',separator='=',literal_value='G001'),dict(key='transcript',separator='=',literal_value='T002')],
        [dict(key='gene_id',separator=':',literal_value='CaseGene.2')],[]]
    pd.DataFrame(dict(header_row=[0,1,2,3],primary_identifier=['001','001','CaseID','SymbolOnly'],
        original_header=['001 | gene=G001 | transcript=T001  ','001 | gene=G001 | transcript=T002',
                         'CaseID gene_id:CaseGene.2','SymbolOnly GN=SYMBOL'],
        declared_identifier_fields=[json.dumps(value) for value in fields])).to_parquet(folder/'headers.parquet',index=False)
    info=dict(schema_version=1,input_sha256='a'*64,header_rows=4,
        fields=['header_row','primary_identifier','original_header','declared_identifier_fields'],
        scope='Original source declarations only; reference versions and database use are unreviewed',
        outputs={'headers.parquet':digest(folder/'headers.parquet')})
    (folder/'manifest.json').write_text(json.dumps(info))
    (folder/'validation.json').write_text(json.dumps(dict(status='all_header_cells_reconciled',
        bindings={'input_sha256':'a'*64},counts={'headers':4},outputs={p.name:digest(p) for p in folder.iterdir()})))
    source=dict(source_url='https://example.org/source.fasta',declared_projects=['PXD000001'],filename='source.fasta',
        status='all_header_cells_reconciled',pack=folder.name,input_sha256='a'*64,
        validation_sha256=digest(folder/'validation.json'),counts={'headers':4,'headers_with_explicit_gene_declarations':3})
    (tmp_path/'catalogue.json').write_text(json.dumps(dict(schema_version=1,sources=[source,
        dict(source_url='https://example.org/pending.fasta',status='awaiting_acquisition')])) )
    (tmp_path/'status.json').write_text(json.dumps(dict(planned=2,catalogue_sha256=digest(tmp_path/'catalogue.json'))))
    return tmp_path,source


def test_exact_source_identifiers_and_gene_declarations_do_not_collapse_or_change_versions(references):
    root,source=references
    original=R.lookup(source['source_url'],'001',root=root)
    assert len(original)==2 and original.header_row.tolist()==[0,1]
    assert original.original_header.iloc[0].endswith('  ')
    genes=R.lookup(source['source_url'],'G001',kind='gene',root=root)
    assert len(genes)==2 and genes.primary_identifier.tolist()==['001','001']
    assert json.loads(genes.matching_gene_declarations.iloc[0])[0]['literal_value']=='G001'
    assert R.lookup(source['source_url'],'caseid',root=root).empty
    assert R.lookup(source['source_url'],'CaseGene',kind='gene',root=root).empty
    assert len(R.lookup(source['source_url'],'CaseGene.2',kind='gene',root=root))==1
    assert R.lookup(source['source_url'],'SYMBOL',kind='gene',root=root).empty
    assert original.input_sha256.eq('a'*64).all()


def test_header_output_tampering_pending_sources_and_catalogue_changes_are_refused(references):
    root,source=references
    with pytest.raises(ValueError,match='not independently reconciled'):
        R.lookup('https://example.org/pending.fasta','001',root=root)
    (root/source['pack']/'headers.parquet').write_bytes(b'changed original headers')
    with pytest.raises(ValueError,match='output checksum'):
        R.lookup(source['source_url'],'001',root=root)
    (root/'catalogue.json').write_text('{}')
    with pytest.raises(ValueError,match='catalogue is updating'):
        R.catalogue(root)


def test_reference_pack_cannot_escape_its_collection(references):
    root,source=references
    catalogue=json.loads((root/'catalogue.json').read_text())
    catalogue['sources'][0]['pack']='../outside'
    (root/'catalogue.json').write_text(json.dumps(catalogue))
    status=json.loads((root/'status.json').read_text())
    status['catalogue_sha256']=digest(root/'catalogue.json')
    (root/'status.json').write_text(json.dumps(status))
    with pytest.raises(ValueError,match='not been reconciled'):
        R.lookup(source['source_url'],'001',root=root)


def test_reference_browser_exports_exact_rows_and_clears_stale_queries(qtbot,references,tmp_path):
    from parascape.reference_identifiers_panel import ReferenceIdentifiersPanel
    root,_=references
    panel=ReferenceIdentifiersPanel(root)
    qtbot.addWidget(panel)
    assert panel.source.count()==1 and '1 awaiting' in panel.summary.text()
    panel.query.setText('001');panel.show_records()
    assert panel.table.rowCount()==2 and panel.table.item(0,1).text()=='001'
    path=tmp_path/'original-identifiers.csv'
    assert panel.save(str(path))==2
    exported=pd.read_csv(path,dtype=str,keep_default_na=False)
    assert exported.primary_identifier.tolist()==['001','001']
    panel.kind.setCurrentIndex(1)
    assert panel.result is None and panel.table.rowCount()==0
    panel.query.setText('G001');panel.show_records()
    assert panel.table.rowCount()==2
    panel.query.setText('unavailable literal')
    assert panel.result is None and panel.table.rowCount()==0 and panel.save(str(path))==0
