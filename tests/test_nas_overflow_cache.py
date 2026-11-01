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
