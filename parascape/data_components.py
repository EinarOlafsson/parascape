"""Install checksum-pinned optional datasets from individually bounded GitHub assets.

Components have independent immutable directories. Shared axes download once;
large archives use ordered parts, and every part, archive and extracted file is
verified before a directory becomes visible. Core data and other components are
never replaced by an optional-data download.
"""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tempfile
from urllib.parse import urlsplit
import urllib.request
import uuid

from . import paths

PREFIX = '/EinarOlafsson/parascape/releases/download/'
KINDS = {'sparse_coordinates', 'sparse_axes', 'workbook_cells'}


def storage_root() -> Path:
    """Return the shared optional-data cache, honoring an explicit location."""
    return Path(os.environ.get('PARASCAPE_COMPONENTS', str(Path(paths.user_cache_dir()) / 'components')))


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _digest(value):
    if not isinstance(value, str) or not re.fullmatch('[a-f0-9]{64}', value):
        raise ValueError('Invalid optional-data checksum')


def _file_name(value):
    name = PurePosixPath(value)
    if (not value or '\\' in value or name.is_absolute() or '..' in name.parts
            or str(name) != value or name.parts[0] == '.component_release.json'):
        raise ValueError('Unsafe optional-data filename')


def validate(component: dict) -> dict:
    """Validate the pinned component contract before network or filesystem changes."""
    directory = component['directory']
    pieces = directory.split('/')
    if len(pieces) != 2 or pieces[0] not in KINDS:
        raise ValueError('Unsupported optional-data collection')
    _digest(pieces[1])
    if not component['files'] or not component['parts']:
        raise ValueError('Optional data must declare files and ordered archive parts')
    for name, item in component['files'].items():
        _file_name(name)
        _digest(item['sha256'])
        if type(item['size']) is not int or item['size'] < 0:
            raise ValueError('Invalid optional-data file size')
    if (type(component['archive_size']) is not int or component['archive_size'] <= 0
            or sum(item['size'] for item in component['parts']) != component['archive_size']):
        raise ValueError('Optional archive part lengths differ')
    _digest(component['archive_sha256'])
    for part in component['parts']:
        _digest(part['sha256'])
        url = urlsplit(part['url'])
        tail = url.path[len(PREFIX):] if url.path.startswith(PREFIX) else ''
        if (url.scheme != 'https' or url.netloc != 'github.com' or url.query or url.fragment
                or len(tail.split('/')) != 2 or any(not re.fullmatch('[A-Za-z0-9._-]+', token) or token in ('.', '..') for token in tail.split('/'))):
            raise ValueError('Unexpected optional-data download host or release path')
        if type(part['size']) is not int or not 0 < part['size'] <= 1024**3:
            raise ValueError('Optional archive parts must be at most one GiB')
    content = dict(directory=directory, files=component['files'])
    return dict(content_sha256=hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest(),
                download_bytes=component['archive_size'], extracted_bytes=sum(x['size'] for x in component['files'].values()))


def read_index(path=None) -> dict:
    """Read an optional index shipped inside the separately verified core data."""
    path = Path(path) if path is not None else Path(paths.data_dir()) / 'data_components.json'
    if not path.exists():
        return dict(schema_version=1, components={})
    index = json.loads(path.read_text())
    if index['schema_version'] != 1:
        raise ValueError('Unsupported optional-data index')
    for name, component in index['components'].items():
        validate(component)
        if name != component['directory']:
            raise ValueError('Optional component identity differs from its directory')
    return index


def _target(root, component):
    validate(component)
    kind, name = component['directory'].split('/')
    parent = root / kind
    if parent.is_symlink() or (parent.exists() and parent.resolve().parent != root.resolve()):
        raise ValueError('Optional collection escapes storage directory')
    target = parent / name
    if target.is_symlink():
        raise ValueError('Optional component cannot be a symbolic link')
    return target


def status(component, root=None, *, verify=False) -> dict:
    """Check a component without downloading; strict mode hashes every data file."""
    root = Path(root) if root is not None else storage_root()
    target = _target(root, component)
    metadata = validate(component)
    missing = []
    for name, item in component['files'].items():
        file = target.joinpath(*PurePosixPath(name).parts)
        if (not file.is_file() or file.is_symlink() or target.resolve() not in file.resolve().parents
                or file.stat().st_size != item['size'] or (verify and _hash(file) != item['sha256'])):
            missing.append(name)
    receipt = target / '.component_release.json'
    try:
        saved = json.loads(receipt.read_text())
        receipt_matches = (not receipt.is_symlink() and saved['content_sha256'] == metadata['content_sha256']
                           and saved['directory'] == component['directory'])
    except (OSError, ValueError, KeyError, TypeError):
        receipt_matches = False
    return dict(ready=not missing and receipt_matches, missing=missing, directory=str(target), **metadata)


def _checkpoint(cancelled):
    if cancelled and cancelled():
        raise InterruptedError('Optional-data download cancelled; completed parts retained')


@contextmanager
def _lock(root):
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.install.lock').open('a+b') as lock:
        try:
            if os.name == 'nt':
                import msvcrt
                lock.write(b'0'); lock.flush(); lock.seek(0)
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise RuntimeError('Another optional-data installation is running') from error
        yield


def _download_part(part, cache, progress, cancelled):
    target = cache / part['sha256']
    if target.is_symlink():
        raise ValueError('Unsafe cached optional-data file')
    if target.exists() and target.stat().st_size == part['size'] and _hash(target) == part['sha256']:
        return target
    target.unlink(missing_ok=True)
    partial = target.with_suffix('.part')
    if partial.is_symlink():
        raise ValueError('Unsafe partial optional-data file')
    offset = partial.stat().st_size if partial.exists() else 0
    if offset >= part['size']:
        if offset == part['size'] and _hash(partial) == part['sha256']:
            partial.replace(target)
            return target
        partial.unlink()
        offset = 0
    headers = {'User-Agent': 'ParaScape optional data installer'}
    if offset:
        headers['Range'] = f'bytes={offset}-'
    _checkpoint(cancelled)
    with urllib.request.urlopen(urllib.request.Request(part['url'], headers=headers), timeout=60) as response:
        if response.status == 206:
            if response.headers.get('Content-Range') != f'bytes {offset}-{part["size"] - 1}/{part["size"]}':
                raise ValueError('Unexpected optional-data range response')
        elif response.status == 200:
            offset = 0
        else:
            raise ValueError('Unexpected optional-data download response')
        with partial.open('ab' if offset else 'wb') as stream:
            while True:
                _checkpoint(cancelled)
                block = response.read(1024 * 1024)
                if not block:
                    break
                offset += len(block)
                if offset > part['size']:
                    raise ValueError('Optional-data part exceeds its declared length')
                stream.write(block)
                if progress:
                    progress('Downloading optional data', offset, part['size'])
    if partial.stat().st_size != part['size'] or _hash(partial) != part['sha256']:
        partial.unlink(missing_ok=True)
        raise ValueError('Optional-data part checksum or size mismatch')
    partial.replace(target)
    return target


def _install(component, root, progress, cancelled):
    metadata = validate(component)
    target = _target(root, component)
    if status(component, root, verify=True)['ready']:
        return target
    if target.exists():
        try:
            receipt = target / '.component_release.json'
            saved = json.loads(receipt.read_text())
            _digest(saved['content_sha256'])
            if receipt.is_symlink() or saved['directory'] != component['directory']:
                raise ValueError('Unmanaged optional-data target')
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise ValueError('Optional target contains unmanaged files; choose a separate component cache') from error
    cache = root / '.download-parts'
    if cache.is_symlink():
        raise ValueError('Unsafe optional-data part cache')
    cache.mkdir(exist_ok=True)
    required = 2 * component['archive_size'] + metadata['extracted_bytes']
    if shutil.disk_usage(root).free < required + 64 * 1024**2:
        raise OSError('Not enough free space for optional data and verified staging')
    parts = [_download_part(part, cache, progress, cancelled) for part in component['parts']]
    with tempfile.TemporaryDirectory(prefix='.component-stage-', dir=root) as temporary:
        staging = Path(temporary)
        archive = staging / 'archive'
        digest = hashlib.sha256()
        with archive.open('wb') as output:
            for part in parts:
                with part.open('rb') as source:
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        _checkpoint(cancelled)
                        output.write(block)
                        digest.update(block)
        if archive.stat().st_size != component['archive_size'] or digest.hexdigest() != component['archive_sha256']:
            raise ValueError('Combined optional-data archive checksum differs')
        prepared = staging / 'prepared'
        prepared.mkdir()
        expected, seen = component['files'], set()
        with tarfile.open(archive, 'r:*') as bundle:
            for member in bundle:
                _checkpoint(cancelled)
                if (not member.isfile() or member.name not in expected or member.name in seen
                        or member.size != expected[member.name]['size']):
                    raise ValueError('Unexpected optional-data archive member: ' + member.name)
                _file_name(member.name)
                output = prepared.joinpath(*PurePosixPath(member.name).parts)
                output.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(member)
                if source is None:
                    raise ValueError('Unreadable optional-data member')
                with source, output.open('wb') as stream:
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        _checkpoint(cancelled)
                        stream.write(block)
                if _hash(output) != expected[member.name]['sha256']:
                    raise ValueError('Optional-data extracted file checksum differs')
                seen.add(member.name)
                if progress:
                    progress('Verifying optional data', len(seen), len(expected))
        if seen != set(expected):
            raise ValueError('Optional-data archive is incomplete')
        (prepared / '.component_release.json').write_text(json.dumps(dict(content_sha256=metadata['content_sha256'],
            archive_sha256=component['archive_sha256'], directory=component['directory']), indent=2) + '\n')
        _checkpoint(cancelled)
        target.parent.mkdir(exist_ok=True)
        previous = None
        if target.exists():
            previous = target.with_name(target.name + '.previous-' + uuid.uuid4().hex[:8])
            target.rename(previous)
        try:
            prepared.replace(target)
        except Exception:
            if previous is not None:
                previous.rename(target)
            raise
    for part in parts:
        part.unlink(missing_ok=True)
    return target


def download(component, root=None, progress=None, cancelled=None) -> Path:
    """Install one component atomically; verified shared components are reused."""
    validate(component)
    root = Path(root) if root is not None else storage_root()
    with _lock(root):
        return _install(component, root, progress, cancelled)
