"""Install versioned data releases without placing datasets inside the wheel.

Archive and individual file checksums are pinned by the package. Installation is
staged and verified before becoming visible; failed downloads never replace data.

Keeping release data outside the wheel reduces installation size while checksum-pinned assets preserve reproducibility across machines and application upgrades.

The data is hosted on Hugging Face, as a dataset repository under Einar's account (decision of
2026-11-13: a GitHub release has a space limit the data would reach soon, and a dataset repository
does not). The GitHub release asset for each version stays where it is, and is tried second, so an
install of an earlier version keeps working and a Hugging Face outage is not a dead end. Which host
is used is a setting rather than a code change: `$PARASCAPE_DATA_HOST` replaces the base address and
`$PARASCAPE_DATA_URL` replaces the whole address of this version's archive.

Neither host is trusted with what it serves. The archive is accepted only when its size and SHA-256
match the descriptor in the wheel, and then only when every file inside it matches too, so changing
the host changes where the bytes come from and nothing about what is accepted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile
from urllib.parse import urlsplit
import urllib.error
import urllib.request
import uuid

from . import paths

#: Base address of the data host, and the per-version archive address, as settings. A move to
#: another host is then a configuration change for everybody who already has this version.
ENV_HOST = 'PARASCAPE_DATA_HOST'
ENV_URL = 'PARASCAPE_DATA_URL'

#: The Hugging Face dataset repository that now holds every release archive. Paths under it are
#: `<version>/<asset name>`, matching the layout in that repository.
HOST = 'https://huggingface.co/datasets/einarolafsson/parascape-data/resolve/main/releases'

#: Where the published archives were before the move, and where the ones for earlier versions
#: remain. Tried after the host above.
GITHUB_RELEASES = 'https://github.com/EinarOlafsson/parascape/releases/download'

#: Hosts a download may come from. The configured host is added to this set, so setting
#: `$PARASCAPE_DATA_HOST` is enough; an address that matches neither is refused before any request.
HOSTS = ('huggingface.co', 'github.com')


def host() -> str:
    """The configured base address for release archives, without its trailing slash."""
    return (os.environ.get(ENV_HOST) or HOST).rstrip('/')


def _checked_url(url: str) -> str:
    """Refuse an address that is not https on a known host before anything is requested."""
    parts = urlsplit(url)
    allowed = set(HOSTS) | {urlsplit(host()).netloc}
    if (parts.scheme != 'https' or parts.netloc not in allowed or parts.query or parts.fragment
            or '..' in PurePosixPath(parts.path).parts):
        raise ValueError('Unexpected data download host: ' + url)
    return url


def sources(description=None) -> list:
    """Every address for this version's archive, in the order they are tried.

    An explicit `$PARASCAPE_DATA_URL` is used alone, because somebody who names one address means
    that one. Otherwise the configured host comes first and the addresses recorded in the descriptor
    follow, which is what keeps a published GitHub asset working as a fallback.
    """
    description = description if description is not None else manifest()
    explicit = os.environ.get(ENV_URL)
    if explicit:
        return [_checked_url(explicit)]
    # A descriptor records its asset name; one written before it did is read from its address, and
    # one with neither -- a descriptor assembled in a test -- falls back to the name the builder gives.
    asset = (description.get('asset')
             or (PurePosixPath(urlsplit(description['url']).path).name if description.get('url') else '')
             or f"parascape-data-{description['version']}.tar.gz")
    candidates = [f"{host()}/{description['version']}/{asset}"]
    candidates += [description[key] for key in ('url', 'fallback_url') if description.get(key)]
    out = []
    for url in candidates:
        if url not in out:
            out.append(_checked_url(url))
    return out


def manifest() -> dict:
    """Read the release descriptor included with this application version."""
    return json.loads(Path(__file__).with_name('data_manifest.json').read_text())


#: What to call each host where a person reads it. An address nobody has named reads as its own
#: hostname, which is more use in a dialog than a word chosen for it here.
HOST_NAMES = {'huggingface.co': 'Hugging Face', 'github.com': 'GitHub'}


def host_name(description=None) -> str:
    """The first host a download would use, named the way it is written about."""
    netloc = urlsplit(sources(description)[0]).netloc
    return HOST_NAMES.get(netloc, netloc)


def status(destination=None, description=None) -> dict:
    """Report missing or wrong-sized data files, without any network access."""
    description = description if description is not None else manifest()
    root = Path(destination or paths.data_dir())
    missing = [name for name, entry in description['files'].items()
               if not (root / name).is_file() or (root / name).stat().st_size != entry['size']]
    return {'ready': not missing, 'directory': str(root), 'missing': missing,
            'version': description['version'], 'download_bytes': description['size'],
            'host': host_name(description), 'urls': sources(description)}


def _hash(path, progress=None, phase='Checking the downloaded archive'):
    """sha256 of a file, reporting progress. The report matters: this reads 898 MB, and with nothing
    emitted the download dialog sat at whatever it last showed and looked dead."""
    total = Path(path).stat().st_size
    done = 0
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
            done += len(block)
            if progress:
                progress(phase, done, total)
    return h.hexdigest()


class DataLocationError(ValueError):
    """The destination cannot receive a managed download, and no retry will change that.

    A `ValueError` subclass so that every caller which already caught `ValueError` from
    `install_archive` keeps working unchanged.

    Separate from every other failure because it is not retryable: the download dialog used to show
    this in a small label and then offer "Retry download", which cannot succeed however many times
    it is pressed. Raised by `check_destination` before any bytes move.
    """


def check_destination(destination=None) -> Path:
    """Refuse a destination a managed install must not touch, before downloading anything.

    `install_archive` has always refused to overwrite a populated folder that carries no
    `data_release.json`, which is correct -- it is somebody's own data, or, in a source checkout, the
    repository's committed `parascape/data`. The bug was when it refused: after the whole 898 MB
    archive had been fetched and hashed. Checking first turns a wasted download into an immediate,
    actionable message.
    """
    destination = Path(destination or paths.data_dir()).absolute()
    if destination.is_symlink():
        raise DataLocationError('Refusing to replace a symbolic-link data directory: ' + str(destination))
    if destination.exists() and any(destination.iterdir()) and not (destination / 'data_release.json').is_file():
        raise DataLocationError(
            f'{destination} already holds files that ParaScape did not install, so a download would '
            f'have to overwrite them. This is normal when running from a source checkout, where that '
            f'folder is part of the repository. Point PARASCAPE_CACHE at an empty folder and start '
            f'again, for example:\n\n'
            f'    PARASCAPE_CACHE={paths.download_dir()} parascape')
    return destination


def install_archive(archive, destination=None, description=None, progress=None, cancelled=None) -> Path:
    """Verify and install a local release archive, rejecting traversal, links and extra files."""
    description = description if description is not None else manifest()
    destination = Path(destination or paths.data_dir()).absolute()
    archive = Path(archive)
    if archive.stat().st_size != description['size'] or \
            _hash(archive, progress, 'Verifying the archive') != description['sha256']:
        raise ValueError('Data archive checksum or size mismatch')
    destination.parent.mkdir(parents=True, exist_ok=True)
    check_destination(destination)   # also checked before downloading; kept here as the last defence
    staging = Path(tempfile.mkdtemp(prefix='.parascape-data-', dir=destination.parent))
    try:
        seen = set()
        expected = description['files']
        with tarfile.open(archive, 'r:gz') as bundle:
            for member in bundle:
                if cancelled and cancelled():
                    raise InterruptedError('Data installation cancelled')
                name = PurePosixPath(member.name)
                if (not member.isfile() or name.is_absolute() or '..' in name.parts
                        or member.name not in expected or member.name in seen
                        or member.size != expected[member.name]['size']):
                    raise ValueError('Unexpected or unsafe archive member: ' + member.name)
                target = staging.joinpath(*name.parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(member)
                if source is None:
                    raise ValueError('Unreadable archive member')
                with source, target.open('wb') as out:
                    shutil.copyfileobj(source, out, 1024 * 1024)
                if _hash(target) != expected[member.name]['sha256']:
                    raise ValueError('Data file checksum mismatch: ' + member.name)
                seen.add(member.name)
                if progress:
                    progress('Verifying data files', len(seen), len(expected))
        if seen != set(expected):
            raise ValueError('Data archive is incomplete')
        (staging / 'data_release.json').write_text(json.dumps({
            'version': description['version'], 'archive_sha256': description['sha256']}, indent=2) + '\n')
        previous = None
        if destination.exists():
            previous = destination.with_name(destination.name + '.previous-' + uuid.uuid4().hex[:8])
            destination.rename(previous)
        try:
            os.replace(staging, destination)
        except Exception:
            if previous is not None:
                previous.rename(destination)
            raise
        # Keep any previous managed data for recovery; never delete user data here.
        return destination
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def _fetch(url, part, description, progress=None, cancelled=None):
    """Fetch the archive from one address into `part`, resuming whatever is already there."""
    offset = part.stat().st_size if part.exists() else 0
    if offset >= description['size'] and progress:
        # The download loop below is skipped when the part file is already complete, and it was
        # the only thing reporting progress -- so a resumed, finished download read 0 %.
        progress('Download already complete; checking it', description['size'], description['size'])
    if offset > description['size']:
        part.unlink(); offset = 0
    if offset == description['size']:
        return
    headers = {'User-Agent': 'ParaScape data installer'}
    if offset:
        headers['Range'] = f'bytes={offset}-'
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.status != 206:
            offset = 0
        elif response.headers.get('Content-Range', '').split('/')[0] != f'bytes {offset}-{description["size"] - 1}':
            raise ValueError('Unexpected download range response')
        with part.open('ab' if offset else 'wb') as out:
            while True:
                if cancelled and cancelled():
                    raise InterruptedError('Data download cancelled; partial download retained')
                block = response.read(1024 * 1024)
                if not block:
                    break
                offset += len(block)
                if offset > description['size']:
                    raise ValueError('Download exceeds the expected data size')
                out.write(block)
                if progress:
                    progress('Downloading data', offset, description['size'])
    if offset != description['size']:
        raise OSError(f'The data download ended early at {offset} of {description["size"]} bytes')


def _download_locked(destination=None, description=None, progress=None, cancelled=None) -> Path:
    description = description if description is not None else manifest()
    root = Path(destination or paths.data_dir()).absolute()
    root.parent.mkdir(parents=True, exist_ok=True)
    required = description['size'] + sum(x['size'] for x in description['files'].values())
    if shutil.disk_usage(root.parent).free < required + 64 * 1024**2:
        raise OSError('Not enough free space for the download and verified data')
    addresses = sources(description)
    # Before any bytes move. A destination that cannot receive the install is worth knowing now
    # rather than after 898 MB, which is exactly how this failed when Einar first ran it.
    check_destination(root)
    part = root.parent / ('.parascape-' + description['sha256'] + '.part')
    try:
        # Each address in turn, so an unreachable host is not a dead end while the same archive is
        # published on another. Whatever arrives still has to match the pinned checksum below, so
        # trying a second host cannot widen what gets installed. Cancellation stops the whole thing.
        failures = []
        for index, url in enumerate(addresses):
            try:
                _fetch(url, part, description, progress, cancelled)
                break
            except (OSError, ValueError) as error:
                if isinstance(error, InterruptedError):
                    raise
                failures.append(f'{urlsplit(url).netloc}: {error}')
                if index == len(addresses) - 1:
                    raise OSError('The data download failed on every host -- ' + '; '.join(failures)) from error
                if progress:
                    progress(f'{urlsplit(url).netloc} did not answer; trying another host', 0, description['size'])
        if _hash(part, progress) != description['sha256']:
            part.unlink()
            raise ValueError('Data download checksum mismatch; retry the download')
        result = install_archive(part, root, description, progress, cancelled)
        part.unlink(missing_ok=True)
        return result
    except InterruptedError:
        raise


def download(destination=None, description=None, progress=None, cancelled=None) -> Path:
    """Download with a cross-process lock, preventing simultaneous cache writers."""
    root = Path(destination or paths.data_dir()).absolute()
    root.parent.mkdir(parents=True, exist_ok=True)
    with (root.parent / '.parascape-data-install.lock').open('a+b') as lock:
        try:
            if os.name == 'nt':
                import msvcrt
                lock.seek(0); lock.write(b'0'); lock.flush(); lock.seek(0)
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError('Another ParaScape data download is running. Please retry when it finishes.') from exc
        return _download_locked(root, description, progress, cancelled)


def main(argv=None):
    """Inspect local data or explicitly download it without starting the desktop app."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--directory')
    args = parser.parse_args(argv)
    if args.download:
        print(download(args.directory))
    else:
        print(json.dumps(status(args.directory), indent=2))


if __name__ == '__main__':
    main()
