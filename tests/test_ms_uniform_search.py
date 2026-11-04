"""Uniform MS re-search: settings follow the MS2 analyser, and acquisition types are flagged, not guessed."""
import gzip
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('ms_uniform_search', ROOT / 'scripts' / 'ms_uniform_search.py')
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)


def mzml(tmp_path, filters):
    path = tmp_path / 'run.mzML.gz'
    with gzip.open(path, 'wt') as fh:
        for f in filters:
            fh.write(f'<cvParam cvRef="MS" accession="MS:1000512" value="{f}" name="filter string" />\n')
    return path


def test_fragment_tolerance_follows_the_ms2_analyser(tmp_path):
    orbitrap = mzml(tmp_path, ['FTMS + p NSI Full ms [350-1400]', 'FTMS + c NSI d Full ms2 584.3@hcd30.00 [120-2000]'])
    assert M.ms2_analyser(orbitrap) == 'FTMS'
    config = M.sage_config([orbitrap], 'FTMS', tmp_path)
    assert config['fragment_tol'] == {'ppm': [-20.0, 20.0]} and config['deisotope'] is True
    assert config['database']['static_mods'] == {'C': 57.021464} and config['quant'] == {'lfq': True}
    assert M.sage_config([orbitrap], 'ITMS', tmp_path)['fragment_tol'] == {'da': [-0.5, 0.5]}


def test_mixed_or_unknown_analysers_are_refused(tmp_path):
    mixed = mzml(tmp_path, ['FTMS + c NSI d Full ms2 1@hcd', 'ITMS + c NSI d Full ms2 2@cid'])
    try:
        M.ms2_analyser(mixed)
    except ValueError:
        pass
    else:
        raise AssertionError('mixed MS2 analysers must not be searched with one tolerance')


def test_acquisition_flags(tmp_path):
    assert M.acquisition(mzml(tmp_path, ['FTMS + c NSI d Full ms2 1@hcd']))['flag'] == 'DDA'
    assert M.acquisition(mzml(tmp_path, ['FTMS + c NSI d Full ms2 1@cid', 'FTMS + c NSI sps d Full ms3 1@hcd']))['flag'].startswith('possible TMT')
    assert M.acquisition(mzml(tmp_path, ['FTMS + c NSI Full ms2 500.0@hcd [145-1450]']))['flag'].startswith('possible DIA')


def test_file_kinds_cover_the_pride_raw_category():
    assert M.file_kind('run1.RAW') == 'raw' and M.file_kind('a.mzML.gz') == 'mzml'
    assert M.file_kind('batch.raw.zip') == 'archive' and M.file_kind('x.tar.gz') == 'archive'
    assert M.file_kind('timstof.d.zip') == 'bruker'                          # searched by Sage directly
    assert M.file_kind('r.rar') == 'archive' and M.file_kind('r.7z') == 'archive'   # unpacked with 7-Zip
    for name in ('s.wiff', 's.wiff.scan', 'GRA1-BirA-Raw.xml'):             # SCIEX: no Linux converter
        assert M.file_kind(name) == 'unsupported'


def test_author_mzml_is_kept_as_deposited_and_its_source_removed(tmp_path):
    spectra = tmp_path / 'spectra'; spectra.mkdir()
    source = tmp_path / 'author.mzML'
    source.write_text('<cvParam value="FTMS + c NSI d Full ms2 1@hcd" name="filter string" />\n')
    record = M.convert(source, spectra, dict(source='author.mzML'))
    assert record['mzml'] == 'author.mzML.gz' and record['ms2_analyser'] == 'FTMS' and not source.exists()
    assert 'as deposited' in record['conversion']


def test_unpack_refuses_paths_outside_the_folder(tmp_path):
    import zipfile
    bad = tmp_path / 'bad.zip'
    with zipfile.ZipFile(bad, 'w') as z:
        z.writestr('../escape.raw', b'x')
    try:
        M.unpack(bad, tmp_path / 'out')
    except ValueError:
        pass
    else:
        raise AssertionError('a member leaving the folder must be refused')


def test_timstof_runs_get_their_own_tolerances_without_ms1_quantities(tmp_path):
    config = M.sage_config([tmp_path / 'run.d'], 'TIMS', tmp_path)
    assert config['precursor_tol'] == {'ppm': [-15.0, 15.0]} and config['fragment_tol'] == {'ppm': [-20.0, 20.0]}
    assert config['quant'] == {'lfq': False} and config['deisotope'] is False


def test_bruker_folder_is_recorded_with_its_stored_zip(tmp_path):
    d = tmp_path / 'work' / 'run.d'; d.mkdir(parents=True); (d / 'analysis.tdf').write_bytes(b'tdf')
    spectra = tmp_path / 'spectra'; spectra.mkdir()
    record = M.bruker(d, spectra, dict(source='run.d.zip'))
    assert (spectra / record['stored']).exists() and record['ms2_analyser'] == 'TIMS'
    assert record['search_path'] == str(d) and len(record['stored_sha256']) == 64


def test_truncated_deposited_tar_keeps_members_before_the_cut(tmp_path):
    import io
    import tarfile
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode='w:gz') as t:
        for name, size in (('a.raw', 1000), ('b.raw', 200000)):
            info = tarfile.TarInfo(name); info.size = size
            t.addfile(info, io.BytesIO(bytes(range(256)) * (size // 256) + bytes(size % 256)))
    archive = tmp_path / 'cut.tar.gz'
    archive.write_bytes(buffer.getvalue()[:len(buffer.getvalue()) // 2])
    note = M.unpack(archive, tmp_path / 'out')
    assert 'truncated as deposited' in note
    assert (tmp_path / 'out' / 'a.raw').stat().st_size == 1000 and not (tmp_path / 'out' / 'b.raw').exists()
