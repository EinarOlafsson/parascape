"""Resolve checksum-bound quantitative packs without network or source replacement."""
from __future__ import annotations

import json
from pathlib import Path

from . import data_components as C

MARKER = 'optional_data.json'


def plan(pack, *, index_path=None) -> dict:
    """Plan one registered metadata-only pack, retaining its original identity."""
    from . import deposited_values as D
    pack = Path(pack)
    marker_path = pack / MARKER
    if not marker_path.exists():
        return dict(registered=False, components=[], missing=[], runtime_pack=pack)
    if pack.is_symlink() or any((pack / name).is_symlink() for name in (MARKER, 'manifest.json', 'validation.json')):
        raise ValueError('Optional measurement metadata must be ordinary files')
    marker = json.loads(marker_path.read_text())
    if marker.get('schema_version') != 1:
        raise ValueError('Unsupported optional measurement metadata')
    identity = marker['component_id'];C._digest(identity)
    info = D.manifest(pack, verify_outputs=False)
    manifest_sha = C._hash(pack / 'manifest.json')
    validation_sha = C._hash(pack / 'validation.json')
    if (marker['manifest_sha256'] != manifest_sha
            or marker['validation_sha256'] != validation_sha):
        raise ValueError('Optional measurement metadata checksum differs')
    if index_path is None:
        candidates = [root / 'data_components.json' for root in
                      (pack.parent, pack.parent.parent, pack.parent.parent.parent)]
        candidates = [path for path in candidates if path.exists()]
        if len(candidates) != 1:
            raise ValueError('Optional measurement needs one installed component index')
        index_path = candidates[0]
    index_path = Path(index_path)
    if index_path.is_symlink():
        raise ValueError('Optional measurement index must be an ordinary file')
    index = C.read_index(index_path)
    declared = index['deposited_values'][identity]
    directory = 'deposited_values/' + identity
    component = index['components'][directory]
    expected = {**info['outputs'], 'manifest.json': manifest_sha, 'validation.json': validation_sha}
    extras = set(component['files']) - set(expected)
    if (component['directory'] != directory or extras - {'preparation.json', 'technical_qc.csv'}
            or any(component['files'].get(name, {}).get('sha256') != digest for name, digest in expected.items())
            or declared['components'] != [directory]
            or declared['source_url'] != info['source_url'] or declared['input_sha256'] != info['input_sha256']
            or declared['manifest_sha256'] != manifest_sha or declared['validation_sha256'] != validation_sha):
        raise ValueError('Optional measurement component and original source differ')
    try:
        relative = pack.resolve().relative_to(index_path.parent.resolve()).as_posix()
    except ValueError:
        raise ValueError('Optional measurement metadata lies outside its core collection') from None
    if relative != declared['metadata_directory']:
        raise ValueError('Optional measurement metadata location differs')
    if index.get('publication_status') != 'published_verified':
        raise ValueError('This optional measurement release is not yet available for download')
    status = C.status(component)
    return dict(registered=True, components=[component], missing=[] if status['ready'] else [component],
                runtime_pack=Path(status['directory']))


def runtime_pack(pack) -> Path:
    """Use complete local files or an already installed registered component; never fetch."""
    from . import deposited_values as D
    pack = Path(pack)
    info = D.manifest(pack, verify_outputs=False)
    if all((pack / name).is_file() for name in info['outputs']):
        return pack
    selected = plan(pack)
    if not selected['registered']:
        raise FileNotFoundError('Original measurement payload is missing from this local collection: ' + str(pack))
    if selected['missing']:
        raise FileNotFoundError('Download the selected original measurement table before querying it')
    runtime = selected['runtime_pack']
    # The normal reader rechecks all runtime output hashes; bind the original
    # manifest/receipt here as well before any identifiers or annotations are used.
    for name, entry in selected['components'][0]['files'].items():
        if name not in info['outputs'] and C._hash(runtime / name) != entry['sha256']:
            raise ValueError('Installed measurement identity or processing record differs')
    return runtime
