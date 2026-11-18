"""Build a deterministic data release archive and its small package descriptor.

The descriptor records where the archive is published: the Hugging Face dataset repository first,
then the GitHub release asset, which stays published and is the fallback. `parascape.data_manager`
tries them in that order and verifies the checksums here against whatever arrives.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parascape import __version__, data_manager


def _refuse_an_empty_primary_address(info: dict) -> None:
    """A release may not ship a descriptor whose primary address is empty or not a host we publish to.

    Instruction 99 step 3. The descriptor is what `data_manager` reads to find the archive, and an
    empty or malformed primary address is not noticed until somebody installs the release and the
    download fails. It is cheap to refuse here and expensive to discover there.
    """
    from urllib.parse import urlparse
    for field, required in (('url', True), ('fallback_url', False)):
        value = str(info.get(field) or '').strip()
        if not value:
            if required:
                raise ValueError(f'the release descriptor has no {field}; it cannot be published')
            continue
        parsed = urlparse(value)
        if parsed.scheme != 'https' or not parsed.netloc:
            raise ValueError(f'the release descriptor {field} is not an https address: {value!r}')
        if not parsed.path.rstrip('/').endswith(info['asset']):
            raise ValueError(f'the release descriptor {field} does not end in the asset name '
                             f'{info["asset"]!r}: {value!r}')
        if (parsed.hostname or '') not in data_manager.HOSTS:
            raise ValueError(f'the release descriptor {field} is on {parsed.hostname!r}, which is '
                             f'not one of the hosts ParaScape downloads from: {data_manager.HOSTS}')


def build(source: Path, output: Path, descriptor: Path):
    """Archive only application data; retain icons in the application wheel."""
    files = []
    suffixes = ('.npz', '.parquet', '.tsv', '.tsv.gz', '.json')
    files.extend(p for p in source.iterdir() if p.is_file() and p.name.endswith(suffixes)
                 and p.name != 'data_release.json')
    files.extend(p for p in (source / 'hf_release').glob('*') if p.suffix in ('.parquet', '.md'))
    files.extend((source / 'icons').glob('*.svg'))
    files.extend(p for p in (source / 'rnaseq').rglob('*') if p.is_file() and p.suffix in ('.parquet', '.tsv', '.json'))
    # The information map release: the three aligned tables, their evidence and pairs, and the tier-2
    # measurement tables the panel reads beside them. Without these the map panel ships with nothing to show.
    for folder in ('published_expression', 'processing_audit', 'source_matrices', 'platform_identifiers',
                   'deposited_values', 'sparse_datasets', 'information_map'):
        files.extend(p for p in (source / folder).rglob('*') if p.is_file()
                     and not any(part.startswith('.') for part in p.relative_to(source / folder).parts)
                      and p.suffix in ('.parquet', '.tsv', '.csv', '.json', '.jsonl', '.py', '.md'))
    # Workbook payloads are independent optional components. Only their discovery
    # and verification metadata belongs in the core download.
    files.extend(p for p in (source / 'workbook_cells').rglob('*.json') if p.is_file()
                 and not any(part.startswith('.') for part in p.relative_to(source / 'workbook_cells').parts))
    if not files or not (source / 'nodes.parquet').is_file():
        raise ValueError('Cannot publish an empty data release')
    output.mkdir(parents=True, exist_ok=True)
    name = f'parascape-data-{__version__}.tar.gz'
    asset = output / name
    entries = {}
    with asset.open('wb') as raw, gzip.GzipFile(fileobj=raw, filename='', mode='wb', mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as archive:
            for path in sorted(files):
                if path.is_symlink() or any(parent.is_symlink() for parent in path.parents
                                            if parent == source or source in parent.parents):
                    raise ValueError(f'Data release cannot contain symlinks: {path}')
                data = path.read_bytes()
                relative = path.relative_to(source).as_posix()
                entries[relative] = {'size': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
                info = tarfile.TarInfo(relative)
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                archive.addfile(info, io.BytesIO(data))
    digest = hashlib.sha256(asset.read_bytes()).hexdigest()
    info = {'schema_version': 1, 'version': __version__, 'asset': name,
            'url': f'{data_manager.HOST}/{__version__}/{name}',
            'fallback_url': f'{data_manager.GITHUB_RELEASES}/v{__version__}/{name}',
            'size': asset.stat().st_size, 'sha256': digest, 'files': entries}
    _refuse_an_empty_primary_address(info)
    descriptor.write_text(json.dumps(info, indent=2, sort_keys=True) + '\n')
    print(f'{asset}: {asset.stat().st_size / 1e6:.1f} MB, {len(entries)} verified files')
    return info


def main():
    """Build the asset deterministically so CI can verify the checked-in descriptor."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=Path('parascape/data'))
    p.add_argument('--output', type=Path, default=Path('data-dist'))
    p.add_argument('--descriptor', type=Path, default=Path('parascape/data_manifest.json'))
    a = p.parse_args()
    build(a.source, a.output, a.descriptor)


if __name__ == '__main__':
    main()
