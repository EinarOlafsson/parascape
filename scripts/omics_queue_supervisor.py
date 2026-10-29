"""Recover the existing bounded queue independently of chat or model availability.

Run as a user-systemd timer. Completion means every planned receipt and expected
compressed output has been verified, including gzip integrity. Incomplete passes
are retried indefinitely with bounded backoff; unavailable inputs stay pending.
"""
from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time
import zlib

try:
    from .omics_download_cache import digest, save
except ImportError:
    from omics_download_cache import digest, save

SERVICES = ('parascape-downloads.service', 'parascape-resource-guard.service',
            'parascape-quant-cached.service')
TIMER = 'parascape-queue.timer'
HEADERS = {'salmon': b'Name\tLength\tEffectiveLength\tTPM\tNumReads',
           'kallisto': b'target_id\tlength\teff_length\test_counts\ttpm'}


def read_json(path: Path, default=None):
    """Read optional durable state; malformed state is rebuilt, never counted complete."""
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def load_plan(path: Path, expected_hash: str) -> list[dict]:
    """Validate the frozen, nonempty plan before controlling any worker."""
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_hash:
        raise ValueError('frozen plan checksum changed; operator review required')
    rows = [json.loads(line) for line in raw.splitlines() if line.strip()]
    seen = set()
    for row in rows:
        run, study, refs = row['run_accession'], row['study'], row['references']
        identity = (study, run)
        if (not re.fullmatch(r'[SED]RR\d+', run) or identity in seen
                or not re.fullmatch(r'[A-Za-z0-9_-]+', study) or not refs
                or len(refs) != len(set(refs))
                or any(not re.fullmatch(r'[A-Za-z0-9_-]+', ref) for ref in refs)):
            raise ValueError('invalid or duplicate entry in frozen plan')
        seen.add(identity)
    if not rows:
        raise ValueError('empty plan is not a completed queue')
    return rows


def signature(paths: list[Path]) -> list:
    """Detect receipt/output changes without rehashing unchanged verified results."""
    answer = []
    for path in paths:
        stat = path.stat()
        if not path.is_file() or not stat.st_size:
            raise ValueError(f'empty or non-file output: {path}')
        answer.append([str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns])
    return answer


def verify_table(path: Path, header: bytes) -> str:
    """Check the expected table, consume the gzip CRC, and hash compressed bytes."""
    with gzip.open(path, 'rb') as source:
        if source.readline(1024).rstrip(b'\r\n') != header:
            raise ValueError(f'wrong table header: {path}')
        row = source.readline(16384)
        if not row.strip() or len(row.rstrip(b'\r\n').split(b'\t')) != 5:
            raise ValueError(f'missing or malformed table rows: {path}')
        while source.read(1024 * 1024):
            pass
    return digest(path)


def inspect_run(row: dict, outputs: Path, cached: dict | None) -> dict:
    """Verify one run's identity and all requested tool/reference outputs."""
    folder = outputs / row['study'] / row['run_accession']
    receipt = folder / 'done.json'
    if not receipt.exists():
        return {'valid': False, 'reason': 'no completion receipt', 'repair_receipt': False}
    try:
        metadata = json.loads(receipt.read_text())
        if not isinstance(metadata, dict):
            raise ValueError('completion receipt must be a JSON object')
        if any(metadata.get(key) != value for key, value in row.items()):
            raise ValueError('completion receipt does not match the frozen plan')
        if not all(metadata.get(key) for key in ('finished', 'salmon', 'kallisto')):
            raise ValueError('completion receipt lacks completion/tool provenance')
        tables = [(folder / f'{tool}_{ref}' / filename, HEADERS[tool])
                  for ref in row['references']
                  for tool, filename in (('salmon', 'quant.sf.gz'), ('kallisto', 'abundance.tsv.gz'))]
        paths = [receipt] + [path for path, _ in tables]
        stamp = signature(paths)
        if isinstance(cached, dict) and cached.get('signature') == stamp:
            return {'valid': True, 'evidence': cached, 'reused': True}
        hashes = {str(path): verify_table(path, header) for path, header in tables}
        hashes[str(receipt)] = digest(receipt)
        if signature(paths) != stamp:
            return {'valid': False, 'reason': 'outputs changed during verification', 'repair_receipt': False}
        return {'valid': True, 'evidence': {'signature': stamp, 'sha256': hashes}, 'reused': False}
    except (OSError, EOFError, ValueError, TypeError, KeyError, zlib.error) as error:
        return {'valid': False, 'reason': str(error), 'repair_receipt': True}


def audit(rows: list[dict], outputs: Path, verified: dict, active: str | None,
          checkpoint=lambda verified, count: None) -> tuple[dict, list]:
    """Preserve good results and quarantine false receipts so existing retries can run."""
    accepted, pending = {}, []
    fresh = 0
    for row in rows:
        run = row['run_accession']
        key = f"{row['study']}/{run}"
        if active in (run, key):
            pending.append({'run': run, 'study': row['study'], 'reason': 'worker active; awaiting final cleanup'})
            continue
        result = inspect_run(row, outputs, verified.get(key))
        if result['valid']:
            accepted[key] = result['evidence']
            if not result['reused']:
                fresh += 1
                if fresh % 10 == 0:
                    checkpoint(accepted, len(pending))
        else:
            pending.append({'run': run, 'study': row['study'], 'reason': result['reason']})
            if result.get('repair_receipt'):
                receipt = outputs / row['study'] / run / 'done.json'
                try:
                    receipt.rename(receipt.with_name(f'done.invalid-{time.time_ns()}.json'))
                except FileNotFoundError:
                    pass
    return accepted, pending


class Services:
    """Control only the three owned workers and their queue timer."""
    def state(self, unit: str) -> str:
        result = subprocess.run(['systemctl', '--user', 'show', unit, '-p', 'ActiveState', '--value'],
                                capture_output=True, text=True, check=True, timeout=30)
        state = result.stdout.strip()
        if unit == SERVICES[2]:
            self.worker_unit = unit
        if unit == SERVICES[2] and state not in ('active', 'activating', 'reloading'):
            # Claude may have resumed the original transient worker after handoff.
            # Adopt its activity rather than starting a second output writer.
            legacy = subprocess.run(['systemctl', '--user', 'show', 'parascape-quant.service',
                                     '-p', 'ActiveState', '--value'], capture_output=True,
                                    text=True, check=True, timeout=30).stdout.strip()
            if legacy in ('active', 'activating', 'reloading'):
                self.worker_unit = 'parascape-quant.service'
                return legacy
        return state

    def command(self, action: str, *units: str) -> None:
        if action not in ('start', 'stop') or any(unit not in (*SERVICES, TIMER) for unit in units):
            raise ValueError('supervisor may control only its owned services')
        subprocess.run(['systemctl', '--user', action, *units], check=True, timeout=120)

    def worker_status(self) -> str:
        """Expose signal-paused processes without overriding another operator's pause."""
        unit = getattr(self, 'worker_unit', SERVICES[2])
        group = Path('/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice') / unit
        try:
            for pid in (group / 'cgroup.procs').read_text().split():
                text = (Path('/proc') / pid / 'stat').read_text()
                if text.rsplit(')', 1)[1].split()[0] in ('T', 't'):
                    return 'signal_paused'
        except (OSError, IndexError):
            pass
        return 'running'


def recover(services, previous: dict, completed: int, now: float) -> dict:
    """Restart missing workers; normal incomplete exits never mean queue completion."""
    state = dict(previous)
    if completed > state.get('last_completed', -1):
        state.update(attempt=0, retry_at=0)
    state['last_completed'] = completed
    for unit in SERVICES[:2]:
        if services.state(unit) not in ('active', 'activating', 'reloading'):
            services.command('start', unit)
    worker = SERVICES[2]
    if services.state(worker) in ('active', 'activating', 'reloading'):
        state['worker'] = services.worker_status() if hasattr(services, 'worker_status') else 'running'
    elif now >= state.get('retry_at', 0):
        services.command('start', worker)
        state['worker'] = 'restarted'
        state['attempt'] = state.get('attempt', 0) + 1
        state['retry_at'] = now + min(60 * 2 ** min(state['attempt'] - 1, 6), 3600)
    else:
        state['worker'] = 'retry_wait'
    return state


def tick(plan: Path, outputs: Path, folder: Path, services=None) -> dict:
    """Audit one recovery cycle and publish durable status or a verified completion receipt."""
    folder.mkdir(parents=True, exist_ok=True)
    status_path = folder / 'status.json'
    services = services or Services()
    if (folder / 'PAUSE').exists():
        status = {'phase': 'paused', 'reason': 'operator PAUSE marker present', 'updated_epoch': time.time()}
        save(status_path, status)
        return status
    config = read_json(folder / 'config.json', {})
    rows = load_plan(plan, config['plan_sha256'])
    active = read_json(folder.parent / 'active.json', {}) or {}
    # Ignore stale active metadata after a crash, so a final result can still be verified.
    try:
        command = (Path('/proc') / str(int(active.get('pid', 0))) / 'cmdline').read_bytes()
        process_alive = b'omics_quant.py' in command
    except (OSError, TypeError, ValueError):
        process_alive = False
    active_run = (f"{active['study']}/{active['run']}" if active.get('study') else active.get('run')) if process_alive else None
    if active_run is None:
        legacy = Path('/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice/parascape-quant.service/cgroup.procs')
        try:
            for pid in legacy.read_text().split():
                try:
                    args = (Path('/proc') / pid / 'cmdline').read_bytes().decode().split('\0')
                    if Path(args[0]).name not in ('salmon', 'kallisto') or '-o' not in args:
                        continue
                    parts = Path(args[args.index('-o') + 1]).relative_to(outputs).parts
                    if len(parts) >= 2:
                        active_run = '/'.join(parts[:2])
                        break
                except (OSError, ValueError, IndexError, UnicodeError):
                    continue
        except OSError:
            pass
    verified = read_json(folder / 'verified.json', {})
    previous = read_json(folder / 'recovery.json', {})
    verified = verified if isinstance(verified, dict) else {}
    previous = previous if isinstance(previous, dict) else {}
    # Recover crashed workers before a potentially long initial integrity audit.
    recovery = previous
    if not (folder / 'completion.json').exists():
        recovery = recover(services, previous, previous.get('last_completed', 0), time.time())
        save(folder / 'recovery.json', recovery)
    save(status_path, {'phase': 'verifying', 'planned': len(rows), 'updated_epoch': time.time()})

    def checkpoint(accepted, pending_count):
        # Preserve old cached evidence for rows not reached yet if this audit is interrupted.
        save(folder / 'verified.json', {**verified, **accepted})
        save(status_path, {'phase': 'verifying', 'planned': len(rows), 'verified_so_far': len(accepted),
                           'pending_so_far': pending_count, 'updated_epoch': time.time()})

    accepted, pending = audit(rows, outputs, verified, active_run, checkpoint)
    save(folder / 'verified.json', accepted)
    status = {'phase': 'running' if pending else 'complete', 'planned': len(rows),
              'verified_complete': len(accepted), 'remaining': len(pending),
              'pending': pending, 'updated_epoch': time.time(),
              'processing_service': getattr(services, 'worker_unit', SERVICES[2])}
    if pending:
        (folder / 'completion.json').unlink(missing_ok=True)
        recovery = recover(services, recovery, len(accepted), time.time())
        save(folder / 'recovery.json', recovery)
        status['recovery'] = recovery
    else:
        completion = {'phase': 'complete', 'planned': len(rows), 'verified_complete': len(accepted),
                      'plan_sha256': config['plan_sha256'],
                      'verified_manifest_sha256': digest(folder / 'verified.json'),
                      'finished': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
        save(folder / 'completion.json', completion)
        services.command('stop', *SERVICES)
        services.command('stop', TIMER)
    save(status_path, status)
    return status


def main() -> None:
    """Run one locked timer cycle; interrupted audits and missing storage are retried."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--outputs', type=Path, required=True)
    parser.add_argument('--state', type=Path, required=True)
    args = parser.parse_args()
    args.state.mkdir(parents=True, exist_ok=True)
    with (args.state / '.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        try:
            status = tick(args.plan, args.outputs, args.state)
        except Exception as error:
            status = {'phase': 'waiting', 'reason': str(error), 'updated_epoch': time.time()}
            save(args.state / 'status.json', status)
            print(json.dumps(status), flush=True)
            raise
        print(json.dumps({key: value for key, value in status.items() if key != 'pending'}), flush=True)


if __name__ == '__main__':
    main()
