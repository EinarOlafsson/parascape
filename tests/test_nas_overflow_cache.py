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
