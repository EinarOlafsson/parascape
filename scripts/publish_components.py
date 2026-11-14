"""Publish the optional component parts to the Hugging Face dataset repository.

One kind at a time, resumable. Every part is fetched from the public GitHub release asset the
component index names, checked against the SHA-256 in that index, uploaded to
`components/<kind>/<part name>`, checked again from what the repository stored, and only then
deleted from staging, so peak disk stays at about one batch. A part whose checksum differs is never
uploaded. Only assets that are already public are published; nothing is read from the NAS or from
`results/`.

Run it with an environment that has `huggingface_hub` and a write token. That library is a tool
dependency only: the application fetches parts over `urllib` with ranged requests and must keep no
dependency on it.

    python scripts/publish_components.py sparse_axes \\
        --index /path/to/data/data_components.json --staging /path/with/room

Run it again with the same arguments to continue; parts the repository already holds with the right
size and checksum are skipped.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

REPO = 'einarolafsson/parascape-data'
KINDS = ('sparse_axes', 'workbook_cells', 'deposited_values', 'sparse_coordinates')


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def parts_for(index_path, kind) -> list:
    """Every distinct part of one kind, as the published index records it."""
    index = json.loads(Path(index_path).read_text())
    seen, out = set(), []
    for component in index['components'].values():
        if component['directory'].split('/')[0] != kind:
            continue
        for part in component['parts']:
            name = part['url'].rsplit('/', 1)[1]
            if name not in seen:
                seen.add(name)
                out.append(dict(name=name, url=part['url'], size=part['size'], sha256=part['sha256']))
    return sorted(out, key=lambda part: part['name'])


def published(api, kind) -> dict:
    """What the repository already holds for this kind: part name -> (size, SHA-256 or None).

    A part small enough to be stored as an ordinary file rather than a large one has no recorded
    checksum, which is why `confirm` re-downloads those.
    """
    out = {}
    try:
        for entry in api.list_repo_tree(REPO, f'components/{kind}', repo_type='dataset',
                                        expand=True, recursive=False):
            if getattr(entry, 'size', None) is None:
                continue
            large = getattr(entry, 'lfs', None)
            out[entry.path.rsplit('/', 1)[1]] = (entry.size, large.sha256 if large else None)
    except Exception as error:
        if 'does not exist' not in str(error) and '404' not in str(error):
            raise
    return out


def fetch(part, directory, attempts=4) -> Path:
    """Fetch one part from its published GitHub asset and return it only when its checksum matches."""
    target = Path(directory) / part['name']
    for attempt in range(attempts):
        if target.exists() and target.stat().st_size == part['size'] and sha256(target) == part['sha256']:
            return target
        target.unlink(missing_ok=True)
        done = subprocess.run(['curl', '-sSL', '--retry', '5', '--retry-delay', '5',
                               '--max-time', '3600', '-o', str(target), part['url']],
                              capture_output=True, text=True)
        if done.returncode == 0 and target.exists() and target.stat().st_size == part['size']:
            if sha256(target) == part['sha256']:
                return target
            print(f"  checksum differs from the index: {part['name']}", flush=True)
        time.sleep(5 * (attempt + 1))
    raise SystemExit(f"could not fetch a verified copy of {part['name']}")


def confirm(api, kind, batch, staging):
    """Check what the repository stored for each part of this batch against the published index."""
    have = published(api, kind)
    for part in batch:
        size, digest = have.get(part['name'], (None, None))
        if size != part['size']:
            raise SystemExit(f"uploaded size differs for {part['name']}: {size} not {part['size']}")
        if digest is None:
            copy = api.hf_hub_download(REPO, f"components/{kind}/{part['name']}", repo_type='dataset',
                                       local_dir=str(Path(staging) / '.confirm'))
            digest = sha256(copy)
            os.unlink(copy)
        if digest != part['sha256']:
            raise SystemExit(f"uploaded checksum differs for {part['name']}")


def publish(kind, index_path, staging, batch_files=25, batch_bytes=3 * 1024 ** 3, workers=3):
    from huggingface_hub import HfApi
    from huggingface_hub.hf_api import CommitOperationAdd

    api = HfApi()
    parts = parts_for(index_path, kind)
    have = published(api, kind)
    todo = [p for p in parts if have.get(p['name'], (None, None)) != (p['size'], p['sha256'])]
    print(f"{kind}: {len(parts)} parts, {sum(p['size'] for p in parts) / 1e9:.1f} GB; "
          f"{len(parts) - len(todo)} already published, {len(todo)} to go", flush=True)
    directory = Path(staging) / kind
    directory.mkdir(parents=True, exist_ok=True)
    moved = 0
    while todo:
        batch, total = [], 0
        while todo and len(batch) < batch_files and (not batch or total + todo[0]['size'] <= batch_bytes):
            batch.append(todo.pop(0))
            total += batch[-1]['size']
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(lambda part: fetch(part, directory), batch))
        operations = [CommitOperationAdd(path_in_repo=f"components/{kind}/{part['name']}",
                                         path_or_fileobj=str(directory / part['name'])) for part in batch]
        for attempt in range(4):
            try:
                api.create_commit(REPO, operations, repo_type='dataset',
                                  commit_message=f'components/{kind}: {len(batch)} verified parts')
                break
            except Exception as error:
                print(f'  upload attempt {attempt + 1} failed: {error}', flush=True)
                if attempt == 3:
                    raise
                time.sleep(20 * (attempt + 1))
        confirm(api, kind, batch, staging)
        for part in batch:
            (directory / part['name']).unlink(missing_ok=True)
        moved += len(batch)
        print(f'  {kind}: {moved} published and verified, {len(todo)} left', flush=True)
    checksums = ''.join(f"{part['sha256']}  {part['name']}\n" for part in parts)
    api.upload_file(path_or_fileobj=checksums.encode(), path_in_repo=f'components/{kind}/SHA256SUMS',
                    repo_id=REPO, repo_type='dataset',
                    commit_message=f'components/{kind}: checksums for {len(parts)} parts')
    print(f'{kind}: all {len(parts)} parts published and verified', flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('kind', choices=KINDS)
    parser.add_argument('--index', required=True, help='data_components.json from verified core data')
    parser.add_argument('--staging', required=True, help='a directory with room for one batch')
    arguments = parser.parse_args(argv)
    publish(arguments.kind, arguments.index, arguments.staging)
    return 0


if __name__ == '__main__':
    sys.exit(main())
