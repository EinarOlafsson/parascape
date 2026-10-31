"""Real quantity packs install on demand without changing source identities or values."""
import copy
import json
from pathlib import Path
import shutil

import pandas as pd
import pytest

from parascape import data_components as C, deposited_values as D, optional_deposited_data as O
from scripts import build_data_components as B
from test_deposited_values import pack
from test_data_components import transport


@pytest.fixture
def optional(pack,tmp_path,monkeypatch):
    (pack/'preparation.json').write_text(json.dumps(dict(manifest_sha256=C._hash(pack/'manifest.json'))))
    D.technical_qc(pack).to_csv(pack/'technical_qc.csv',index=False)
    receipt=json.loads((pack/'validation.json').read_text());receipt['technical_qc_sha256']=C._hash(pack/'technical_qc.csv')
    (pack/'validation.json').write_text(json.dumps(receipt))
    files=B.deposited_files(pack);identity=B.deposited_identity(files)
    prepared=B.build_component(pack,'deposited_values/'+identity,tmp_path/'parts',tag='v0.0.9.0-quantities-data',part_limit=4096)
    core=tmp_path/'core';stub=core/'deposited_values/synthetic'/pack.name;stub.mkdir(parents=True)
    for name in ['manifest.json','validation.json']:shutil.copyfile(pack/name,stub/name)
    marker=dict(schema_version=1,component_id=identity,manifest_sha256=C._hash(stub/'manifest.json'),validation_sha256=C._hash(stub/'validation.json'))
    (stub/O.MARKER).write_text(json.dumps(marker))
    (stub.parent/'catalogue.json').write_text(json.dumps(dict(packs=[stub.name])))
    info=D.manifest(pack);component=prepared['component'];directory=component['directory']
    index=dict(schema_version=1,version='0.0.9.0',publication_status='published_verified',test_only_publication_marker=True,
        components={directory:component},deposited_values={identity:dict(components=[directory],
            source_url=info['source_url'],input_sha256=info['input_sha256'],manifest_sha256=marker['manifest_sha256'],
            validation_sha256=marker['validation_sha256'],metadata_directory=stub.relative_to(core).as_posix())})
    index_path=core/'data_components.json';index_path.write_text(json.dumps(index))
    folder=Path(prepared['path']);payloads={part['url']:(folder/part['url'].rsplit('/',1)[1]).read_bytes() for part in component['parts']}
    calls=transport(monkeypatch,payloads)
    cache=tmp_path/'cache';monkeypatch.setenv('PARASCAPE_COMPONENTS',str(cache))
    return dict(source=pack,stub=stub,core=core,component=component,index=index,index_path=index_path,calls=calls,cache=cache,prepared=prepared)


def test_real_installer_preserves_all_values_tokens_context_QC_and_operations(optional):
    x=optional;stub=x['stub'];before={p.name:p.read_bytes() for p in stub.iterdir()}
    rows,errors=D.catalogue(x['core']/'deposited_values');assert len(rows)==1 and not errors
    with pytest.raises(FileNotFoundError,match='Download'):D.profiles(stub,'001')
    assert not x['calls']
    plan=O.plan(stub);assert plan['missing']==[x['component']]
    installed=C.download(x['component']);assert C.status(x['component'],verify=True)['ready']
    for name in x['component']['files']:assert (installed/name).read_bytes()==(x['source']/name).read_bytes()
    pd.testing.assert_frame_equal(D.profiles(stub,'001'),D.profiles(x['source'],'001'))
    pd.testing.assert_frame_equal(D.technical_qc(stub),D.technical_qc(x['source']))
    assert D.column_context(stub)==D.column_context(x['source'])
    assert {p.name:p.read_bytes() for p in stub.iterdir()}==before
    x['calls'].clear();assert not O.plan(stub)['missing'];D.profiles(stub,'001');assert not x['calls']
    assert B.build_component(x['source'],x['component']['directory'],x['cache']/'staging',tag='v0.0.9.0-quantities-data',part_limit=4096)['component']==x['component']


@pytest.mark.parametrize('corruption',['unpublished','source','manifest','validation','location','dependency','missing_output','foreign_output'])
def test_wrong_registration_refused_before_network(optional,corruption):
    x=optional;index=copy.deepcopy(x['index']);component=next(iter(index['components'].values()));entry=next(iter(index['deposited_values'].values()))
    if corruption=='unpublished':index['publication_status']='provisional'
    elif corruption=='source':entry['input_sha256']='0'*64
    elif corruption=='manifest':entry['manifest_sha256']='0'*64
    elif corruption=='validation':entry['validation_sha256']='0'*64
    elif corruption=='location':entry['metadata_directory']='deposited_values/elsewhere'
    elif corruption=='dependency':entry['components']=[]
    elif corruption=='missing_output':component['files'].pop('source_tokens.parquet')
    else:component['files']['values.parquet']['sha256']='0'*64
    x['index_path'].write_text(json.dumps(index))
    with pytest.raises((ValueError,KeyError)):O.plan(x['stub'])
    assert not x['calls']


def test_changed_cached_payload_or_processing_records_never_silently_replace_source(optional):
    x=optional;installed=C.download(x['component'])
    for name in ['values.parquet','preparation.json']:
        original=(installed/name).read_bytes();(installed/name).write_bytes(b'X'+original[1:])
        with pytest.raises(ValueError):D.profiles(x['stub'],'001')
        (installed/name).write_bytes(original)
    assert not (x['stub']/'values.parquet').exists()


def test_unregistered_local_collection_is_not_redirected_to_cached_data(optional):
    x=optional;C.download(x['component']);(x['source']/'values.parquet').unlink()
    with pytest.raises(FileNotFoundError,match='local collection'):D.profiles(x['source'],'001')


def test_corrupt_QC_or_content_identity_is_refused_before_packaging(optional,tmp_path):
    x=optional
    with pytest.raises(ValueError,match='identity'):B.build_component(x['source'],'deposited_values/'+'0'*64,tmp_path/'wrong')
    (x['source']/'technical_qc.csv').write_text('changed')
    with pytest.raises(ValueError,match='QC checksum'):B.deposited_files(x['source'])


def test_core_metadata_contains_only_receipts_and_roundtrips_original_identity(optional,tmp_path):
    x=optional;core=tmp_path/'small_core';relative='deposited_values/ns/'+x['source'].name
    entry=B.write_deposited_metadata(x['source'],x['component'],core,relative)
    assert set(p.name for p in (core/relative).iterdir())=={'manifest.json','validation.json',O.MARKER}
    assert entry['manifest_sha256']==C._hash(x['source']/'manifest.json')
    assert B.write_deposited_metadata(x['source'],x['component'],core,relative)==entry
    with pytest.raises(ValueError,match='Unsafe|Invalid'):B.write_deposited_metadata(x['source'],x['component'],core,'../escape')
    (core/relative/'validation.json').write_text('changed')
    with pytest.raises(ValueError,match='Existing'):B.write_deposited_metadata(x['source'],x['component'],core,relative)


def test_UI_installs_selected_table_then_reopens_offline(qtbot,optional,monkeypatch,tmp_path):
    from PyQt6 import QtWidgets
    from parascape import deposited_values_panel as P
    x=optional;dialogs=[]
    def accept(dialog):
        dialogs.append(dialog.components)
        for component in dialog.components:C.download(component)
        return QtWidgets.QDialog.DialogCode.Accepted
    monkeypatch.setattr(P.ComponentDownloadDialog,'exec',accept)
    panel=P.DepositedValuesPanel(x['core']/'deposited_values');qtbot.addWidget(panel)
    panel.features.setPlainText('001');panel.show_values()
    assert len(dialogs)==1 and panel.result is not None and len(panel.result)==4
    assert panel.table.item(0,4).text()=='-0.000'
    assert panel.save(str(tmp_path/'values.csv'))==4
    x['calls'].clear();panel.show_values();assert len(dialogs)==1 and not x['calls']
    panel.show_column_context();assert panel.context_dialog is not None
    assert panel.save_qc(str(tmp_path/'qc.csv'))==2


@pytest.mark.parametrize('changed',[False,True])
def test_UI_cancel_or_query_change_during_download_leaves_no_stale_result(qtbot,optional,monkeypatch,changed):
    from PyQt6 import QtWidgets
    from parascape import deposited_values_panel as P
    x=optional;panel=P.DepositedValuesPanel(x['core']/'deposited_values');qtbot.addWidget(panel)
    def finish(dialog):
        if changed:
            C.download(x['component']);panel.features.setPlainText('CaseID')
            return QtWidgets.QDialog.DialogCode.Accepted
        return QtWidgets.QDialog.DialogCode.Rejected
    monkeypatch.setattr(P.ComponentDownloadDialog,'exec',finish)
    panel.features.setPlainText('001');panel.show_values()
    assert panel.result is None and panel.table.rowCount()==0 and not panel.export_button.isEnabled()
    assert ('changed' if changed else 'cancelled') in panel.summary.text()
