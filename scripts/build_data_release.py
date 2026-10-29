"""Build a deterministic GitHub data asset and its small package descriptor."""
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
from parascape import __version__


def build(source: Path, output: Path, descriptor: Path):
    """Archive only application data; retain icons in the application wheel."""
    files = []
    suffixes = ('.npz', '.parquet', '.tsv', '.tsv.gz', '.json')
    files.extend(p for p in source.iterdir() if p.is_file() and p.name.endswith(suffixes))
    files.extend(p for p in (source / 'hf_release').glob('*') if p.suffix in ('.parquet', '.md'))
    files.extend((source / 'icons').glob('*.svg'))
    files.extend(p for p in (source / 'rnaseq').rglob('*') if p.is_file() and p.suffix in ('.parquet', '.tsv', '.json'))
    if not files or not (source / 'nodes.parquet').is_file():
        raise ValueError('Cannot publish an empty data release')
    output.mkdir(parents=True, exist_ok=True)
    name = f'parascape-data-{__version__}.tar.gz'
    asset = output / name
    entries = {}
    with asset.open('wb') as raw, gzip.GzipFile(fileobj=raw, filename='', mode='wb', mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as archive:
            for path in sorted(files):
                if path.is_symlink():
                    raise ValueError(f'Data release cannot contain symlinks: {path}')
                data = path.read_bytes()
                relative = path.relative_to(source).as_posix()
                entries[relative] = {'size': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
                info = tarfile.TarInfo(relative)
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                archive.addfile(info, io.BytesIO(data))
    digest = hashlib.sha256(asset.read_bytes()).hexdigest()
    info = {'schema_version': 1, 'version': __version__, 'asset': name,
            'url': f'https://github.com/EinarOlafsson/parascape/releases/download/v{__version__}/{name}',
            'size': asset.stat().st_size, 'sha256': digest, 'files': entries}
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
