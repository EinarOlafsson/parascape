"""Give only the background worker extra CPU capacity while the workstation is quiet.

Uses process/thread affinity, which is enforced even when user cgroups have no CPU
controller. Tutorial processes are observed but never modified.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time

try:
    from .omics_download_cache import pressure_reason, save
except ImportError:
    from omics_download_cache import pressure_reason, save

TUTORIAL_MARKERS = ('run_native_capture.py', 'gnome_capture_host.py',
                    'capture_clip_gate.py', 'spacr-tutorial-current', 't17-captures')


def tutorial_active(proc_root: Path = Path('/proc')) -> bool:
    """Detect capture/render workers without mistaking shell command text for a worker."""
    for process in proc_root.iterdir():
        if not process.name.isdigit():
            continue
        try:
            args = (process / 'cmdline').read_bytes().split(b'\0')
        except OSError:
            continue
        executable = Path(os.fsdecode(args[0])).name if args and args[0] else ''
        if not executable.startswith(('python', 'gst-launch', 'ffmpeg')):
            continue
        command = b' '.join(args).decode(errors='replace')
        if any(marker in command for marker in TUTORIAL_MARKERS):
            return True
    return False


def apply_affinity(cgroup: Path, cores: set[int], proc_root: Path = Path('/proc')) -> int:
    """Limit every existing thread in this worker's cgroup, including new child tools."""
    try:
        pids = (cgroup / 'cgroup.procs').read_text().split()
    except FileNotFoundError:
        return 0
    changed = 0
    for pid in pids:
        if not pid.isdigit():
            continue
        try:
            threads = list((proc_root / pid / 'task').iterdir())
        except FileNotFoundError:
            continue
        for thread in threads:
            try:
                tid = int(thread.name)
                if os.sched_getaffinity(tid) != cores:
                    os.sched_setaffinity(tid, cores)
                    changed += 1
            except ProcessLookupError:
                pass
    return changed


def select_mode(reason: str, now: float, quiet_since: float | None,
                quiet_seconds: float) -> tuple[str, float | None]:
    """Yield immediately and require sustained quiet before restoring extra capacity."""
    if reason:
        return 'protected', None
    if quiet_since is None:
        quiet_since = now
    return ('fast' if now - quiet_since >= quiet_seconds else 'protected'), quiet_since


def main() -> None:
    """Monitor resource pressure and tutorial workers; adapt only our processing tree."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cgroup', type=Path, required=True)
    parser.add_argument('--status', type=Path, required=True)
    parser.add_argument('--fast-cpus', default='4,5,6,7,8,9,10,11,12,13')
    parser.add_argument('--protected-cpus', default='10,11,12,13,14')
    parser.add_argument('--quiet-seconds', type=float, default=30)
    parser.add_argument('--restore', action='store_true')
    parser.add_argument('--full-resources', action='store_true')
    args = parser.parse_args()
    fast = {int(cpu) for cpu in args.fast_cpus.split(',')}
    protected = {int(cpu) for cpu in args.protected_cpus.split(',')}
    if args.full_resources:
        fast = protected = set(range(os.cpu_count() or 1))
    args.status.parent.mkdir(parents=True, exist_ok=True)
    if args.restore:
        apply_affinity(args.cgroup, protected)
        save(args.status, {'mode': 'protected', 'reason': 'guard stopped; conservative affinity restored',
                           'cpus': sorted(protected), 'updated': time.strftime('%Y-%m-%dT%H:%M:%S%z')})
        return
    quiet_since = None
    previous_mode = None
    while True:
        if args.full_resources:
            reason, mode, cores = 'explicit full-resource user setting', 'full', fast
        else:
            reason = 'tutorial capture or validation active' if tutorial_active() else pressure_reason()
            mode, quiet_since = select_mode(reason, time.monotonic(), quiet_since, args.quiet_seconds)
            cores = fast if mode == 'fast' else protected
        changed = apply_affinity(args.cgroup, cores)
        record = {'mode': mode, 'reason': reason or ('quiet' if mode == 'fast' else 'waiting for sustained quiet'),
                  'cpus': sorted(cores), 'threads_adjusted': changed,
                  'updated': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
        save(args.status, record)
        if previous_mode != mode or changed:
            print(json.dumps(record), flush=True)
        previous_mode = mode
        time.sleep(5)


if __name__ == '__main__':
    main()
