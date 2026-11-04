"""The NAS raw tier has a hard byte cap and the quantifier's scratch margin is configurable.

Both limits exist because the storage is shared: the NAS tier once grew from 214 GB to
1.58 TB overnight when only a free-space reserve bounded it, and the scratch disk must
leave room for other work such as tutorial recordings. No network, NAS or Qt is used.
"""
import importlib.util
import os
from pathlib import Path
import sys
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    sys.path.insert(0, str(ROOT / 'scripts'))
    try:
        spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.remove(str(ROOT / 'scripts'))


def test_nas_cap_blocks_only_past_the_cap():
    N = load('nas_overflow_cache')
    gib = N.GIB
    assert not N.nas_cap_blocks(400 * gib, 66 * gib, 466 * gib)      # exactly at the cap is allowed
    assert N.nas_cap_blocks(400 * gib, 66 * gib + 1, 466 * gib)
    assert N.nas_cap_blocks(1583 * gib, 0 + 1, 466 * gib)             # over-full tier admits nothing new
    assert not N.nas_cap_blocks(10 ** 15, 10 ** 12, None)             # uncapped keeps the old behaviour


def test_scratch_margin_defaults_to_25_and_follows_the_environment():
    with mock.patch.dict(os.environ, {}, clear=False):
        os.environ.pop('PARASCAPE_SCRATCH_MARGIN_GB', None)
        assert load('omics_quant').SCRATCH_MARGIN_GB == 25
    with mock.patch.dict(os.environ, {'PARASCAPE_SCRATCH_MARGIN_GB': '125'}):
        assert load('omics_quant').SCRATCH_MARGIN_GB == 125


def test_queue_tools_exclude_single_cell_collections_and_unreferenced_hosts():
    sys.argv, saved = ['queue_unplanned_rna_studies.py'], sys.argv
    try:
        Q = load('queue_unplanned_rna_studies')
    finally:
        sys.argv = saved
    assert {'GSE145800', 'GSE145080'} <= set(Q.EXCLUDED)
    assert set(Q.REFERENCE_HOSTS) == {'Homo sapiens', 'Mus musculus'}
    hosts = load('omics_quant').HOST_TRANSCRIPTS
    assert all(any(name.replace(' ', '_').lower() in url.lower() for url in hosts) for name in Q.REFERENCE_HOSTS)


def test_every_host_species_has_its_own_competitive_reference(tmp_path, monkeypatch):
    """Dual RNA-seq from a host absent from the reference lets conserved host reads count as parasite."""
    import gzip
    Q = load('omics_quant')
    assert Q.HOST_REFERENCES['ME49_host'] == Q.HOST_TRANSCRIPTS                     # existing meaning unchanged
    for species, name in Q.HOST_SPECIES.items():
        urls = Q.HOST_REFERENCES[name]
        assert any(species.replace(' ', '_').lower() in url.lower() for url in urls)     # this host's own files
        assert any('/cdna/' in u for u in urls) and any('/ncrna/' in u for u in urls)
    me49 = tmp_path / 'references' / 'ME49'; me49.mkdir(parents=True)
    cdna_url, dna_url = Q.REFERENCES['ME49']
    with gzip.open(me49 / Path(cdna_url).name, 'wt') as fh:
        fh.write('>TGME49_t1 cdna gene:TGME49_g1\nACGT\n')
    with gzip.open(me49 / Path(dna_url).name, 'wt') as fh:
        fh.write('>chrIa\nACGTACGT\n')
    (me49 / 'decoys.txt').write_text('chrIa\n')
    commands = []

    def fetch(url, dest):
        with gzip.open(dest, 'wt') as fh:
            fh.write('>PIG_t1 cdna gene:PIG_g1\nTTTT\n')
        return 'f' * 64
    monkeypatch.setattr(Q, 'NAS', str(tmp_path)); monkeypatch.setattr(Q, 'SCRATCH', str(tmp_path / 'scratch'))
    monkeypatch.setattr(Q, 'fetch', fetch); monkeypatch.setattr(Q, 'version', lambda tool: tool + ' fixture')
    monkeypatch.setattr(Q, 'run', lambda cmd, **kw: commands.append(cmd))
    Q.host_reference('ME49_pig')
    out = tmp_path / 'references' / 'ME49_pig'
    rows = dict(line.split('\t') for line in (out / 'tx2gene.tsv').read_text().splitlines())
    assert rows['TGME49_t1'] == 'TGME49_g1' and rows['PIG_t1'] == 'HOST:PIG_g1'
    salmon = next(c for c in commands if c[1] == 'index' and 'salmon' in c[0])
    assert salmon[salmon.index('-d') + 1] == str(me49 / 'decoys.txt')                # ME49 genome decoys
    assert salmon[salmon.index('-i') + 1] == str(out / 'salmon_index')
    meta = __import__('json').loads((out / 'reference.json').read_text())
    assert meta['name'] == 'ME49_pig' and meta['host'] == list(Q.HOST_REFERENCES['ME49_pig'])
    (out / 'salmon_index').mkdir(); (out / 'salmon_index' / 'info.json').write_text('{}')  # salmon 2.x
    (out / 'kallisto.idx').write_text('')
    commands.clear(); Q.host_reference('ME49_pig')
    assert commands == []                                                             # built indexes are reused


def test_host_routing_names_one_reference_per_declared_host():
    R = load('route_host_references')
    assert R.host_of('Toxoplasma gondii; Sus scrofa') == 'ME49_pig'
    assert R.host_of('Toxoplasma gondii; Rattus norvegicus') == 'ME49_rat'
    assert R.host_of('Homo sapiens; Toxoplasma gondii') is None       # already served by ME49_host
    assert R.host_of('Toxoplasma gondii') is None
    try:
        R.host_of('Sus scrofa; Bos taurus')
    except ValueError:
        pass
    else:
        raise AssertionError('two foreign hosts must not be resolved silently')


def test_a_run_too_large_for_scratch_is_deferred_not_blocking_the_queue(tmp_path, monkeypatch):
    """One 328 GB run once stalled every other run behind it for hours."""
    Q = load('omics_quant')
    scratch = tmp_path / 'scratch'; scratch.mkdir()
    for name in ('SRR1_1.fastq', 'SRR2.sra', 'fasterq.tmp.x', 'keep.txt'):
        (scratch / name).write_text('left over')
    monkeypatch.setattr(Q, 'SCRATCH', str(scratch)); monkeypatch.setattr(Q, 'NAS', str(tmp_path / 'nas'))
    monkeypatch.setattr(Q, 'DOWNLOAD_CACHE', '')
    Q._clear_scratch(None)
    assert sorted(p.name for p in scratch.iterdir()) == ['keep.txt']                 # only run files removed
    done = []

    def reads(row):
        if row['run_accession'] == 'SRR9':
            Q._need_scratch('SRR9', 10 ** 9)                                         # never fits
        done.append(row['run_accession'])
        return []                                                                    # "no reads": skipped
    monkeypatch.setattr(Q, '_reads', reads)
    rows = [dict(run_accession=r, study='S', references=['ME49']) for r in ('SRR9', 'SRR3', 'SRR4')]
    deferred = Q._quantify_pass(rows, None, 0)
    assert done == ['SRR3', 'SRR4'] and [r['run_accession'] for r in deferred] == ['SRR9']


def test_a_run_too_large_locally_uses_nas_scratch_then_local_again(tmp_path, monkeypatch):
    Q = load('omics_quant')
    local, nas = tmp_path / 'local', tmp_path / 'nas_scratch'; local.mkdir()
    monkeypatch.setattr(Q, 'SCRATCH', str(local)); monkeypatch.setattr(Q, 'NAS', str(tmp_path / 'nas'))
    monkeypatch.setattr(Q, 'DOWNLOAD_CACHE', ''); monkeypatch.setattr(Q, 'NAS_SCRATCH', str(nas))
    monkeypatch.setattr(Q, 'SCRATCH_MARGIN_GB', 125.0); monkeypatch.setattr(Q, 'NAS_SCRATCH_RESERVE_GB', 150.0)
    free = {str(local): 160.0, str(nas): 2300.0}
    monkeypatch.setattr(Q, '_free_gb', lambda path: free[str(path)])
    used = []

    def reads(row):
        Q._need_scratch(row['run_accession'], 328.0 if row['run_accession'] == 'SRR1' else 130.0)
        used.append((row['run_accession'], Q.SCRATCH))
        return []
    monkeypatch.setattr(Q, '_reads', reads)
    rows = [dict(run_accession=r, study='S', references=['ME49']) for r in ('SRR1', 'SRR2')]
    assert Q._quantify_pass(rows, None, 0) == []
    assert used == [('SRR1', str(nas)), ('SRR2', str(local))] and Q.SCRATCH == str(local)
    monkeypatch.setattr(Q, 'NAS_SCRATCH', '')                                  # no NAS scratch: deferred as before
    used.clear()
    assert [r['run_accession'] for r in Q._quantify_pass(rows, None, 0)] == ['SRR1'] and used == [('SRR2', str(local))]


def test_reads_stream_gzip_locally_or_convert_on_nas_by_space(tmp_path, monkeypatch):
    """Stream gzip locally when 1.6x the archive fits; otherwise fastq-dump --gzip on NAS scratch."""
    Q = load('omics_quant')
    local, nas = tmp_path / 'local', tmp_path / 'nas'; local.mkdir(); nas.mkdir()
    monkeypatch.setattr(Q, 'SCRATCH', str(local)); monkeypatch.setattr(Q, 'NAS_SCRATCH', str(nas))
    monkeypatch.setattr(Q, 'DOWNLOAD_CACHE', ''); monkeypatch.setattr(Q, 'SCRATCH_MARGIN_GB', 125.0)
    monkeypatch.setattr(Q, '_head_size', lambda url: 25 * 10 ** 9)                  # 25 GB archive -> 325 GB needed
    monkeypatch.setattr(Q, 'fetch', lambda url, dest: Path(dest).write_bytes(b''))
    monkeypatch.setattr(Q, 'log', lambda message: None)
    commands = []

    def run(cmd, **kw):
        tool = Path(cmd[0]).name; commands.append(tool)
        out = Path(cmd[cmd.index('-O') + 1]) if '-O' in cmd else None
        if tool == 'fasterq-dump':
            (out / 'SRR1_1.fastq').write_text(''); (out / 'SRR1_2.fastq').write_text('')
        if tool == 'fastq-dump':
            (out / 'SRR1_1.fastq.gz').write_text(''); (out / 'SRR1_2.fastq.gz').write_text('')
    monkeypatch.setattr(Q, 'run', run)
    row = dict(run_accession='SRR1', fastq_ftp='ftp.sra.ebi.ac.uk/a/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/a/SRR1_2.fastq.gz',
               fastq_md5='', fastq_bytes=f'{19 * 10 ** 9};{19 * 10 ** 9}')
    streamed = []

    def stream(run_acc, sra):
        streamed.append(run_acc)
        return [str(local / f'{run_acc}_1.fastq.gz'), str(local / f'{run_acc}_2.fastq.gz')]
    monkeypatch.setattr(Q, '_stream_gzip', stream)
    free = {str(local): 200.0}                                                    # 25 GB x 1.6 + 125 fits locally
    monkeypatch.setattr(Q, '_free_gb', lambda path: free.get(str(path), 2300.0))
    assert [Path(p).name for p in Q._reads(row)] == ['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz'] and streamed == ['SRR1']
    assert 'fasterq-dump' not in commands
    free[str(local)] = 146.0; commands.clear()                                    # local too small: gzip on NAS
    files = Q._reads(row)
    assert [Path(p).name for p in files] == ['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz'] and 'fasterq-dump' not in commands
    assert all(Path(p).parent == nas / 'SRR1' for p in files) and 'fastq-dump' in commands
    monkeypatch.setattr(Q, 'SCRATCH', str(nas)); commands.clear()                # running on NAS scratch: same
    assert 'fasterq-dump' not in commands and Q._reads(row) and 'fastq-dump' in commands
    monkeypatch.setattr(Q, 'SCRATCH', str(local)); monkeypatch.setattr(Q, 'NAS_SCRATCH', '')
    free[str(local)] = 100.0                                                      # nothing fits and no NAS: deferred
    try:
        Q._reads(row)
    except Q.ScratchDeferred:
        pass
    else:
        raise AssertionError('expected deferral')


def test_runs_held_back_by_hand_stay_pending_without_blocking(tmp_path, monkeypatch):
    Q = load('omics_quant')
    monkeypatch.setattr(Q, 'SCRATCH', str(tmp_path)); monkeypatch.setattr(Q, 'NAS', str(tmp_path / 'nas'))
    monkeypatch.setattr(Q, 'DOWNLOAD_CACHE', ''); monkeypatch.setattr(Q, 'DEFER_RUNS', {'SRR9'})
    monkeypatch.setattr(Q, 'log', lambda message: None)
    seen = []
    monkeypatch.setattr(Q, '_reads', lambda row: seen.append(row['run_accession']) or [])
    rows = [dict(run_accession=r, study='S', references=['ME49']) for r in ('SRR9', 'SRR3')]
    assert Q._quantify_pass(rows, None, 0) == [] and seen == ['SRR3']


def test_preconverted_fastq_is_used_only_with_a_matching_marker(tmp_path, monkeypatch):
    import json
    Q = load('omics_quant')
    folder = tmp_path / 'SRR7'; folder.mkdir()
    for k in (1, 2):
        (folder / f'SRR7_{k}.fastq.gz').write_bytes(b'x' * (10 * k))
    monkeypatch.setattr(Q, 'PRECONVERTED', str(tmp_path)); monkeypatch.setattr(Q, 'log', lambda message: None)
    record = dict(run='SRR7', tool='fastq-dump --split-3 --gzip 3.4.1',
                  files={f'SRR7_{k}.fastq.gz': dict(bytes=10 * k, sha256='') for k in (1, 2)})
    (folder / 'done.json').write_text(json.dumps(record))
    assert [Path(p).name for p in Q._reads(dict(run_accession='SRR7'))] == ['SRR7_1.fastq.gz', 'SRR7_2.fastq.gz']
    (folder / 'SRR7_2.fastq.gz').write_bytes(b'truncated')
    try:
        Q._reads(dict(run_accession='SRR7'))
    except RuntimeError:
        pass
    else:
        raise AssertionError('a changed file must be refused')


def test_two_workers_process_each_run_once_with_their_own_scratch(tmp_path, monkeypatch):
    import os, time
    Q = load('omics_quant')
    base = tmp_path / 'scratch'; base.mkdir(); record = tmp_path / 'record.txt'
    monkeypatch.setattr(Q, 'SCRATCH', str(base)); monkeypatch.setattr(Q, 'NAS', str(tmp_path / 'nas'))
    monkeypatch.setattr(Q, 'DOWNLOAD_CACHE', ''); monkeypatch.setattr(Q, 'WORKERS', 2)
    monkeypatch.setattr(Q, 'THREADS', 16); monkeypatch.setattr(Q, 'log', lambda message: None)

    def reads(row):
        with open(record, 'a') as fh:
            fh.write(f"{os.getpid()} {row['run_accession']} {Q.SCRATCH} {Q.THREADS}\n")
        time.sleep(0.3)
        return []
    monkeypatch.setattr(Q, '_reads', reads)
    rows = [dict(run_accession=f'SRR{i}', study='S', references=['ME49']) for i in range(5)]
    Q._quantify_workers(rows, None)
    lines = [line.split() for line in record.read_text().splitlines()]
    assert sorted(run for _, run, _, _ in lines) == [f'SRR{i}' for i in range(5)]          # each run once
    assert len({pid for pid, *_ in lines}) == 2                                              # both workers worked
    assert {scratch for *_, scratch, _ in lines} == {str(base / 'w0'), str(base / 'w1')}
    assert {threads for *_, threads in lines} == {'8'}
    assert os.listdir(base / 'claims') == []                                                 # every claim released


def test_scratch_reserved_by_another_worker_counts_against_free_space(tmp_path, monkeypatch):
    Q = load('omics_quant')
    claims = tmp_path / 'claims'; claims.mkdir(); (claims / 'SRR1.claim').write_text('150.0')
    monkeypatch.setattr(Q, 'CLAIMS', str(claims)); monkeypatch.setattr(Q, 'SCRATCH', str(tmp_path))
    monkeypatch.setattr(Q, 'NAS_SCRATCH', ''); monkeypatch.setattr(Q, 'SCRATCH_MARGIN_GB', 125.0)
    monkeypatch.setattr(Q, '_free_gb', lambda path: 300.0)
    assert Q._claim('SRR2') and not Q._claim('SRR2')                                        # atomic claim
    Q._need_scratch('SRR2', 140.0)                                                          # 300-150 >= 140
    assert (claims / 'SRR2.claim').read_text() == '15.0'                                    # reserved minus margin
    try:
        Q._need_scratch('SRR3', 200.0)                                                      # 300-150-15 < 200
    except Q.ScratchDeferred:
        pass
    else:
        raise AssertionError('other workers\' reservations must be counted')



def test_stream_gzip_splits_pairs_and_refuses_unpaired_reads(tmp_path, monkeypatch):
    """The real pipeline (bash, awk, pigz) on a fake fastq-dump that prints interleaved reads."""
    import os, stat
    Q = load('omics_quant')
    bin_dir = tmp_path / 'bin'; bin_dir.mkdir()
    pigz = __import__('shutil').which('pigz') or str(Path.home() / 'anaconda3/envs/parascape-quant/bin/pigz')
    os.symlink(pigz, bin_dir / 'pigz')
    def fake(records):
        script = bin_dir / 'fastq-dump'
        script.write_text('#!/bin/bash\n' + ''.join(f"printf '@r.{i}/{k}\\nACGT\\n+\\nIIII\\n'\n" for i, k in records))
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(Q, 'BIN', str(bin_dir)); monkeypatch.setattr(Q, 'SCRATCH', str(tmp_path))
    monkeypatch.setattr(Q, 'THREADS', 4); monkeypatch.setattr(Q, 'log', lambda message: None)
    fake([(1, 1), (1, 2), (2, 1), (2, 2)])
    files = Q._stream_gzip('SRR5', 'x.sra')
    assert [Path(f).name for f in files] == ['SRR5_1.fastq.gz', 'SRR5_2.fastq.gz']
    import gzip
    assert gzip.open(files[0], 'rt').read().count('@r.') == 2 and '/2' in gzip.open(files[1], 'rt').read()
    fake([(1, 1), (2, 1)])                                                        # single-end
    assert [Path(f).name for f in Q._stream_gzip('SRR6', 'x.sra')] == ['SRR6_1.fastq.gz']
    fake([(1, 1), (1, 2), (2, 1)])                                                # one unpaired read
    assert Q._stream_gzip('SRR7', 'x.sra') is None and not (tmp_path / 'SRR7_1.fastq.gz').exists()


def test_a_redownloaded_run_already_on_the_nas_becomes_a_pointer(tmp_path):
    import hashlib
    import json
    C = load('nas_overflow_cache')
    local, nas = tmp_path / 'local', tmp_path / 'nas'
    local.mkdir(); nas.mkdir()
    payload = b'sra archive bytes'
    receipt = dict(url='https://example/SRR1', bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    for folder, finished in ((local, '2026-11-03'), (nas, '2026-10-31')):
        (folder / 'SRR1.sra').write_bytes(payload)
        (folder / 'SRR1.json').write_text(json.dumps(dict(receipt, finished=finished)))
    cache = C.Cache(local, nas, 0, 10 ** 9, tmp_path / 'events.jsonl')
    cache.spill('SRR1')
    assert (local / 'SRR1.sra').is_symlink() and (local / 'SRR1.sra').read_bytes() == payload
    (nas / 'SRR1.sra').write_bytes(b'different bytes!!')       # same size, other content: still refused
    (local / 'SRR1.sra').unlink(); (local / 'SRR1.sra').write_bytes(payload)
    try:
        cache.spill('SRR1')
    except ValueError:
        pass
    else:
        raise AssertionError('a NAS copy with other content must be refused')
