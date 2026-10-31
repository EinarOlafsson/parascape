"""Resolve a sparse dataset's pinned optional dependencies without changing overrides."""
from __future__ import annotations

import json
import os
from pathlib import Path

from . import data_components as C


def plan(descriptor, *, coordinate_root=None, axis_root=None, index_path=None) -> dict:
    """Find missing shared components; explicit local roots are never replaced."""
    descriptor = Path(descriptor)
    C._digest(descriptor.stem)
    if descriptor.is_symlink() or C._hash(descriptor) != descriptor.stem:
        raise ValueError('Sparse dataset descriptor checksum differs')
    info = json.loads(descriptor.read_text())
    index = C.read_index(index_path if index_path is not None else descriptor.parent.parent/'data_components.json')
    if not index['components']:
        return dict(registered=False, components=[], missing=[], coordinate_root=coordinate_root, axis_root=axis_root)
    if index.get('publication_status') != 'published_verified':
        raise ValueError('This optional-data release is not yet available for download')
    references = [('sparse_coordinates', info['coordinate_artifact'])]
    references.extend(('sparse_axes', reference) for role in ('features', 'barcodes')
                      for reference in info['axis_roles'][role]['original_axis_artifacts'])
    required = {}
    for kind, reference in references:
        C._digest(reference['directory'])
        identity = kind + '/' + reference['directory']
        component = index['components'][identity]
        expected = {'manifest.json': reference['manifest_sha256'],
                    'validation.json': reference['validation_sha256'], reference['output_name']: reference['output_sha256']}
        if (component['directory'] != identity or set(component['files']) != set(expected)
                or any(component['files'][name]['sha256'] != digest for name, digest in expected.items())):
            raise ValueError('Optional component does not match the sparse dataset artifacts')
        required[identity] = component
    declared = index['datasets'][descriptor.stem]
    if (declared['source_url'] != info['source_url'] or len(declared['components']) != len(required)
            or set(declared['components']) != set(required)):
        raise ValueError('Optional dataset dependency mapping differs')

    def available(root, component):
        folder = root / component['directory'].split('/')[1]
        if folder.is_symlink():
            return False
        for name, entry in component['files'].items():
            file = folder / name
            if (not file.is_file() or file.is_symlink() or folder.resolve() not in file.resolve().parents
                    or file.stat().st_size != entry['size']):
                return False
            if name in ('manifest.json', 'validation.json') and C._hash(file) != entry['sha256']:
                return False
        return True

    roots, downloadable, missing = {}, [], []
    for kind, override, environment in (
            ('sparse_coordinates', coordinate_root, 'PARASCAPE_SPARSE_COORDINATES'),
            ('sparse_axes', axis_root, 'PARASCAPE_SPARSE_AXES')):
        components = [c for identity, c in required.items() if identity.startswith(kind + '/')]
        explicit = override if override is not None else os.environ.get(environment)
        if explicit is not None:
            root = Path(explicit)
            if any(not available(root, component) for component in components):
                raise ValueError(f'Files are missing or changed in the explicit {environment} location: {root}. Restore that collection or remove the override to use hosted data.')
        else:
            local = descriptor.parent.parent / kind
            if all(available(local, component) for component in components):
                root = local
            else:
                root = C.storage_root() / kind
                downloadable.extend(components)
                missing.extend(component for component in components if not C.status(component)['ready'])
        roots[kind] = root
    return dict(registered=True, components=downloadable, missing=missing,
                coordinate_root=roots['sparse_coordinates'], axis_root=roots['sparse_axes'])
