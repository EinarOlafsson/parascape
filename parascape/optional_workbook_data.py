"""Resolve hosted workbook cells while preserving explicit local collections."""
from __future__ import annotations

import os
from pathlib import Path

from . import data_components as C
from . import workbook_values as W


def plan(folder, *, index_path=None):
    """Bind one workbook's metadata to an immutable optional component."""
    folder = Path(folder)
    info = W.manifest(folder)
    outputs = info['outputs']
    if all((folder / name).is_file() for name in outputs):
        return dict(registered=False, folder=str(folder), components=[], missing=[])
    explicit = os.environ.get('PARASCAPE_WORKBOOK_VALUES')
    if explicit is not None and folder.parent.resolve() == Path(explicit).resolve():
        raise ValueError('Files are missing in the explicit PARASCAPE_WORKBOOK_VALUES collection. Restore that collection or remove the override to use hosted data.')
    C._digest(folder.name)
    index = C.read_index(index_path if index_path is not None else folder.parent.parent / 'data_components.json')
    identity = 'workbook_cells/' + folder.name
    if identity not in index['components']:
        return dict(registered=False, folder=str(folder), components=[], missing=[])
    if index.get('publication_status') != 'published_verified':
        raise ValueError('This workbook data release is not yet available for download')
    component = index['components'][identity]
    W.validation(folder, info)
    expected = {name: record['sha256'] for name, record in outputs.items()}
    expected.update({'manifest.json': W._sha(folder / 'manifest.json'),
                     'independent_validation.json': W._sha(folder / 'independent_validation.json')})
    if (set(component['files']) != set(expected) or any(
            component['files'][name]['sha256'] != digest for name, digest in expected.items())):
        raise ValueError('Hosted workbook component does not match its core source metadata')
    if any(component['files'][name]['size'] != outputs[name]['size'] for name in outputs):
        raise ValueError('Hosted workbook size binding differs')
    target = C.storage_root() / identity
    missing = [] if C.status(component)['ready'] else [component]
    return dict(registered=True, folder=str(target), components=[component], missing=missing)


def resolve(folder):
    """Return local or already installed cells without starting a network transfer."""
    result = plan(folder)
    if result['missing']:
        raise FileNotFoundError('Download the selected workbook before querying its source cells')
    return Path(result['folder'])
