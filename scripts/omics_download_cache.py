"""Bounded, low-bandwidth lookahead downloads for the existing omics run plan.

Only compressed archives are cached. The existing worker still validates each archive,
processes one run at a time, and deletes raw inputs. No biological analysis is added here.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time
import urllib.request

GIB = 1024 ** 3
BLOCK = 1024 ** 2


def save(path: Path, value: dict) -> None:
    """Publish a small JSON record atomically."""
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def digest(path: Path) -> str:
    """Hash an archive without retaining it in memory."""
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(BLOCK), b''):
            h.update(block)
    return h.hexdigest()


def cached_archive(folder: Path, run: str, url: str, size: int) -> Path | None:
    """Return a complete, checksum-verified cache entry, never a partial download."""
    archive = folder / f'{run}.sra'
    with archive_lock(archive):
        return _cached_archive(archive, url, size)


@contextmanager
def archive_lock(archive: Path):
    """Serialize producer and consumer access, including interrupted downloads."""
    with archive.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def _cached_archive(archive: Path, url: str, size: int) -> Path | None:
    try:
        receipt = json.loads(archive.with_suffix('.json').read_text())
        if (receipt['url'] == url and receipt['bytes'] == size
                and archive.stat().st_size == size and digest(archive) == receipt['sha256']):
            return archive
    except (OSError, ValueError, KeyError):
        pass
    return None


def cache_bytes(folder: Path) -> int:
    """Count completed and partially downloaded raw archives against the same limit."""
    total = 0
    for path in folder.iterdir():
        if path.name.endswith(('.sra', '.part')):
            try:
                total += path.stat().st_size
            except FileNotFoundError:
                pass
    return total


def pressure_reason() -> str:
    """Yield under measured CPU/IO pressure or low available system memory."""
    for resource, threshold in (('cpu', 5.0), ('io', 8.0), ('memory', 0.5)):
        try:
            first = Path(f'/proc/pressure/{resource}').read_text().splitlines()[0]
            avg = float(dict(field.split('=') for field in first.split()[1:])['avg10'])
        except (OSError, ValueError, KeyError, IndexError):
            return f'{resource} pressure unavailable'
        if avg >= threshold:
            return f'{resource} pressure {avg:.2f}%'
    fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    if int(fields['MemAvailable'].split()[0]) * 1024 < 32 * GIB:
        return 'less than 32 GiB available RAM'
    return ''


def download(url: str, archive: Path, size: int, *, rate: int, reserve: int,
             budget: int, guard=pressure_reason, report=lambda **kw: None) -> dict:
    """Resume an HTTP archive with enforced transfer, capacity and free-space limits."""
    with archive_lock(archive):
        return _download(url, archive, size, rate=rate, reserve=reserve,
                         budget=budget, guard=guard, report=report)


def _download(url, archive, size, *, rate, reserve, budget, guard, report):
    part = archive.with_suffix('.part')
    offset = part.stat().st_size if part.exists() else 0
    if offset > size:
        raise ValueError('partial archive exceeds remote size')
    if cache_bytes(archive.parent) + size - offset > budget:
        raise RuntimeError('cache capacity reached')
    if shutil.disk_usage(archive.parent).free - (size - offset) < reserve:
        raise RuntimeError('disk reserve reached')
    if offset < size:
        headers = {'Range': f'bytes={offset}-'} if offset else {}
        request = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(request, timeout=120) as response:
            if offset and response.status == 206:
                content_range = response.headers.get('Content-Range', '')
                if not content_range.startswith(f'bytes {offset}-') or not content_range.endswith(f'/{size}'):
                    raise RuntimeError('unexpected resume range')
            elif response.status == 200:
                offset = 0  # Server ignored Range: replace, never append the entire file.
            else:
                raise RuntimeError(f'unexpected HTTP status {response.status}')
            expected = int(response.headers.get('Content-Length', size - offset))
            if expected != size - offset:
                raise RuntimeError('remote size changed')
            with part.open('ab' if offset else 'wb') as target:
                while offset < size:
                    reason = guard()
                    while reason:
                        report(phase='waiting', reason=reason, downloaded_bytes=offset)
                        time.sleep(5)
                        reason = guard()
                    chunk = min(BLOCK, size - offset)
                    if shutil.disk_usage(archive.parent).free - chunk < reserve:
                        raise RuntimeError('disk reserve reached during transfer')
                    if cache_bytes(archive.parent) + chunk > budget:
                        raise RuntimeError('cache capacity reached during transfer')
                    started = time.monotonic()
                    block = response.read(chunk)
                    if not block:
                        raise RuntimeError('truncated download; partial retained for resume')
                    target.write(block)
                    target.flush()
                    offset += len(block)
                    report(phase='downloading', reason='', downloaded_bytes=offset)
                    time.sleep(max(0, len(block) / rate - (time.monotonic() - started)))
    receipt = {'url': url, 'bytes': size, 'sha256': digest(part),
               'finished': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    part.replace(archive)
    save(archive.with_suffix('.json'), receipt)
    return receipt


def pending_rows(plan: list, outputs: Path, folder: Path) -> list:
    """Keep a shared archive until every planned study entry for that run is complete."""
    pending = [row for row in plan
               if not (outputs / row['study'] / row['run_accession'] / 'done.json').exists()]
    complete_runs = {row['run_accession'] for row in plan} - {row['run_accession'] for row in pending}
    for run in complete_runs:
        for suffix in ('.sra', '.part', '.json'):
            (folder / (run + suffix)).unlink(missing_ok=True)
    return pending


def main() -> None:
    """Maintain a small cache of pending inputs until the existing plan is complete."""
    try:
        from . import omics_quant as quant
    except ImportError:
        import omics_quant as quant
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--max-gib', type=float, default=40)
    parser.add_argument('--reserve-gib', type=float, default=90)
    parser.add_argument('--mib-per-second', type=float, default=4)
    args = parser.parse_args()
    if min(args.max_gib, args.reserve_gib, args.mib_per_second) <= 0:
        parser.error('capacity, reserve and transfer rate must be positive')
    folder = args.cache.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    lock = (folder / '.prefetch.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    plan = [json.loads(line) for line in Path(quant.PLAN).read_text().splitlines() if line]
    for row in plan:
        if not re.fullmatch(r'[SED]RR\d+', row['run_accession']) or not re.fullmatch(r'[A-Za-z0-9_-]+', row['study']):
            raise ValueError('invalid run/study identifier in existing plan')
    failures = {}
    last_report = 0.0
    status = {'pid': os.getpid(), 'max_bytes': int(args.max_gib * GIB),
              'reserve_bytes': int(args.reserve_gib * GIB),
              'rate_bytes_per_second': int(args.mib_per_second * BLOCK)}

    def report(**values):
        nonlocal last_report
        changed = values.get('phase') != status.get('phase') or values.get('reason') != status.get('reason')
        status.update(values)
        if changed or time.monotonic() - last_report >= 10:
            status['updated'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
            save(folder / 'status.json', status)
            last_report = time.monotonic()
            if changed:
                print(json.dumps(status), flush=True)

    while True:
        pending = pending_rows(plan, Path(quant.NAS) / 'rnaseq', folder)
        if not pending:
            report(phase='complete', reason='all planned runs complete')
            return
        try:
            active = json.loads((folder / 'active.json').read_text())['run']
        except (OSError, ValueError, KeyError):
            active = pending[0]['run_accession']
        attempted = False
        for row in pending:
            run = row['run_accession']
            archive = folder / f'{run}.sra'
            if run == active or archive.exists() or failures.get(run, 0) > time.monotonic():
                continue
            reason = pressure_reason()
            if reason:
                report(phase='waiting', reason=reason)
                break
            url = quant.SRA_ODP.format(run=run)
            size = quant._head_size(url)
            if not size:
                failures[run] = time.monotonic() + 3600
                continue
            part = archive.with_suffix('.part')
            remaining = size - (part.stat().st_size if part.exists() else 0)
            if (cache_bytes(folder) + remaining > status['max_bytes']
                    or shutil.disk_usage(folder).free - remaining < status['reserve_bytes']):
                report(phase='waiting', reason='cache full or disk reserve', run=run)
                break
            status.update(run=run, expected_bytes=size)
            try:
                download(url, archive, size, rate=status['rate_bytes_per_second'],
                         reserve=status['reserve_bytes'], budget=status['max_bytes'], report=report)
                report(phase='ready', reason='', run=run, downloaded_bytes=size)
            except Exception as error:
                failures[run] = time.monotonic() + 300
                report(phase='retry_later', reason=str(error), run=run)
            attempted = True
            break  # Re-read completion/active state between downloads.
        if not attempted:
            time.sleep(30)


if __name__ == '__main__':
    main()
