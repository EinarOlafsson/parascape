"""Give only the background worker extra CPU capacity while the workstation is quiet.

Einar 2026-11-03: ParaScape and autotrade share the machine; when autotrade (and everything else) is not
using the cores, ParaScape may. 'wide' mode adds the remaining CPUs and lets an extra quantifier worker
claim runs; it starts only after sustained idle and yields as soon as other work returns.

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


APP_SLICE = Path('/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice')


class OtherLoad:
    """Busy cores used by everything except ParaScape units, averaged over the last sample interval."""

    def __init__(self):
        self.last = None
        self.smooth = None

    @staticmethod
    def _sample():
        fields = [int(v) for v in Path('/proc/stat').read_text().splitlines()[0].split()[1:]]
        busy = (sum(fields[:8]) - fields[3] - fields[4]) / os.sysconf('SC_CLK_TCK')   # minus idle, iowait
        ours = 0.0
        for unit in APP_SLICE.glob('parascape-*.service'):
            try:
                stat = dict(line.split() for line in (unit / 'cpu.stat').read_text().splitlines())
                ours += int(stat['usage_usec']) / 1e6
            except (OSError, KeyError, ValueError):
                pass
        return time.monotonic(), busy, ours

    def cores(self) -> float | None:
        now = self._sample()
        previous, self.last = self.last, now
        if previous is None or now[0] - previous[0] <= 0:
            return None
        current = max(0.0, ((now[1] - previous[1]) - (now[2] - previous[2])) / (now[0] - previous[0]))
        self.smooth = current if self.smooth is None else 0.8 * self.smooth + 0.2 * current   # ~25 s at 5 s steps
        return self.smooth


def used_ram_gb() -> float:
    fields = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return (int(fields['MemTotal'].split()[0]) - int(fields['MemAvailable'].split()[0])) / 1024 ** 2


def select_wide(others: float | None, ram_gb: float, now: float, idle_since: float | None, wide: bool,
                enter_max: float, leave_max: float, ram_max: float, sustain: float) -> tuple[bool, float | None]:
    """Widen after `sustain` seconds with other work under enter_max cores; narrow at once above leave_max."""
    if others is None:
        return wide, idle_since
    if ram_gb >= ram_max or others > leave_max:
        return False, None
    if others > enter_max:
        return wide, None
    idle_since = now if idle_since is None else idle_since
    return (wide or now - idle_since >= sustain), idle_since


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
    parser.add_argument('--wide-cpus', default='', help='CPUs added when other work is idle (empty: never widen)')
    parser.add_argument('--wide-enter-cores', type=float, default=6.0)
    parser.add_argument('--wide-leave-cores', type=float, default=10.0)
    parser.add_argument('--wide-ram-gb', type=float, default=80.0)
    parser.add_argument('--wide-sustain-seconds', type=float, default=180.0)
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
    wide_cpus = {int(cpu) for cpu in args.wide_cpus.split(',') if cpu.strip()}
    load, wide, idle_since, others = OtherLoad(), False, None, None
    quiet_since = None
    previous_mode = None
    while True:
        if args.full_resources:
            reason, mode, cores = 'explicit full-resource user setting', 'full', fast
        else:
            reason = 'tutorial capture or validation active' if tutorial_active() else pressure_reason()
            mode, quiet_since = select_mode(reason, time.monotonic(), quiet_since, args.quiet_seconds)
            cores = fast if mode == 'fast' else protected
            if wide_cpus:
                others = load.cores()
                wide, idle_since = select_wide(others, used_ram_gb(), time.monotonic(), idle_since,
                                               wide and mode == 'fast', args.wide_enter_cores,
                                               args.wide_leave_cores, args.wide_ram_gb, args.wide_sustain_seconds)
                wide = wide and mode == 'fast'
                if wide:
                    mode, cores = 'wide', fast | wide_cpus
        changed = apply_affinity(args.cgroup, cores)
        record = {'mode': mode, 'extra_workers': 1 if mode == 'wide' else 0,
                  'other_cores': None if others is None else round(others, 1), 'reason': reason or ('quiet' if mode == 'fast' else 'waiting for sustained quiet'),
                  'cpus': sorted(cores), 'threads_adjusted': changed,
                  'updated': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
        save(args.status, record)
        if previous_mode != mode or changed:
            print(json.dumps(record), flush=True)
        previous_mode = mode
        time.sleep(5)


if __name__ == '__main__':
    main()
