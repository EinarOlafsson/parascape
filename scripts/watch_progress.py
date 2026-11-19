"""Watch every long-running job for a STALL, not just for failure.

Einar, 2026-11-19, after the re-quantification queue sat at 624 for nine hours: *"whenever you run
into blockers like this, learn from them, set up a rule that ensures this does not happen again,
you should routinely be checking the progress of each of your tasks."*

WHY A STALL NEEDS ITS OWN CHECK. A failed job announces itself: the unit goes inactive, the log
says what broke. A stalled job looks healthy from every angle that is usually looked at -- the
service is active, the process is running, the CPU is pegged. The re-quantification worker reported
`624 new runs still to quantify` once an hour for nine hours while three kallisto processes burned
eighteen cores having read their index and written nothing. Nothing was broken. Nothing was moving.

WHAT DISTINGUISHES THEM, AND IT IS NOT CPU. A working quantifier writes steadily; a deadlocked one
has read its index and written zero bytes. `write_bytes` in `/proc/<pid>/io` is the measurement that
tells a slow job from a hung one, and high CPU is no evidence either way -- a futex spin burns a
core just as hard as real work.

HOW IT DECIDES. Each run records a snapshot of every tracked counter. A counter that has not moved
since the previous snapshot is a SUSPECT, not yet a stall: a real job can finish nothing between
two checks when the box is loaded. A counter that has not moved across `STALL_CHECKS` consecutive
snapshots, while a process it owns has written nothing, is reported as STALLED with the evidence.
Load is printed beside it, because "waiting behind someone else's work" is the common innocent
explanation and the report should let a reader see it rather than guess.

It changes nothing. It measures, compares with last time, and says what it found.

    python scripts/watch_progress.py          # one check, appends to the log
    python scripts/watch_progress.py --show   # the last few checks, no new one
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import subprocess

STATE = os.path.expanduser("~/.cache/parascape-progress-watch.json")
LOG = os.path.expanduser("~/.cache/parascape-progress-watch.log")

#: Consecutive unchanged snapshots before a counter is called stalled rather than merely quiet.
STALL_CHECKS = 2

#: How long a counter must have been unchanged before silence means anything, in minutes. A
#: counter is only as fresh as the thing that writes it: the re-quantification worker logs hourly,
#: and the Globus monitor writes its free-space line once a DAY, so comparing either across two
#: checks a minute apart reports a stall that is nothing but the reporting interval. The first
#: version of this script did exactly that and flagged the download disk as stalled (2026-11-19).
QUIET_AFTER_MINUTES = {"requant_queue_remaining": 150, "download_disk_free_gb": 1800}

#: A snapshot taken sooner than this after the last one is not recorded, for the same reason.
MIN_SNAPSHOT_MINUTES = 20

#: Processes whose output is the evidence that work is happening.
WORKERS = ("kallisto", "salmon", "sage", "fastp", "fasterq-dump")


def shell(command: str) -> str:
    try:
        return subprocess.run(command, shell=True, capture_output=True, text=True,
                              timeout=120).stdout.strip()
    except subprocess.SubprocessError:
        return ""


def counters() -> dict:
    """Every number that should be going up, or down, if work is being done."""
    out = {}
    log = shell("systemctl --user status parascape-requant --no-pager -n 1 2>/dev/null")
    m = re.search(r"waits:\s*(\d+)", log)
    if m:
        out["requant_queue_remaining"] = int(m.group(1))
    globus = shell("tail -1 /nas_mnt/carruthers/public_datasets/_download/globus_monitor.log "
                   "2>/dev/null")
    m = re.search(r"(\d+)\s*G free", globus)
    if m:
        out["download_disk_free_gb"] = int(m.group(1))
    return out


def live_workers() -> list:
    """Each running worker process with what it has read and written so far."""
    rows = []
    for name in WORKERS:
        for pid in shell(f"pgrep -x {name}").split():
            try:
                io = open(f"/proc/{pid}/io").read()
                read = int(re.search(r"read_bytes:\s*(\d+)", io).group(1))
                written = int(re.search(r"write_bytes:\s*(\d+)", io).group(1))
                elapsed = int(shell(f"ps -o etimes= -p {pid}").strip() or 0)
                wchan = open(f"/proc/{pid}/wchan").read().strip()
            except (OSError, AttributeError, ValueError):
                continue
            rows.append({"name": name, "pid": int(pid), "minutes": elapsed // 60,
                         "read_mb": round(read / 1048576, 1),
                         "written_mb": round(written / 1048576, 1), "wchan": wchan})
    return rows


def minutes_since_changed(key: str, value, history: list):
    """How long this counter has held its current value, in minutes."""
    last = None
    for check in reversed(history):
        if key not in check["counters"]:
            continue
        if check["counters"][key] != value:
            break
        last = check["when"]
    if last is None:
        return 0
    try:
        then = dt.datetime.fromisoformat(last)
    except ValueError:
        return 0
    return int((dt.datetime.now() - then).total_seconds() // 60)


def load_state() -> dict:
    try:
        return json.load(open(STATE))
    except (OSError, ValueError):
        return {"checks": []}


def verdicts(now: dict, history: list) -> list:
    """One line per counter: moving, quiet once, or stalled with the evidence."""
    out = []
    for key, value in sorted(now["counters"].items()):
        previous = [c["counters"].get(key) for c in history[-STALL_CHECKS:]
                    if key in c["counters"]]
        if not previous:
            out.append((key, "FIRST", f"{value} (nothing to compare with yet)"))
            continue
        if value != previous[-1]:
            out.append((key, "MOVING", f"{previous[-1]} -> {value}"))
            continue
        age = minutes_since_changed(key, value, history)
        needed = QUIET_AFTER_MINUTES.get(key, 150)
        unchanged = (len(previous) >= STALL_CHECKS and all(p == value for p in previous)
                     and age >= needed)
        if not unchanged:
            if age is not None and age < needed:
                out.append((key, "QUIET", f"{value}, unchanged for {age} min; this counter is "
                                          f"only written every {needed} min"))
                continue
            out.append((key, "QUIET", f"{value}, unchanged since the last check"))
            continue
        # Unchanged across several checks. Is anything actually working?
        silent = [w for w in now["workers"] if w["written_mb"] == 0 and w["minutes"] >= 10]
        if silent:
            evidence = "; ".join(f"{w['name']} pid {w['pid']} {w['minutes']}min "
                                 f"read {w['read_mb']} MB written 0 MB ({w['wchan']})"
                                 for w in silent[:3])
            out.append((key, "STALLED", f"{value} for {len(previous)} checks; {evidence}"))
        elif not now["workers"]:
            out.append((key, "STALLED", f"{value} for {len(previous)} checks; no worker running"))
        else:
            out.append((key, "SLOW", f"{value} for {len(previous)} checks, but workers are "
                                     f"writing -- load {now['load1']}"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--show", action="store_true", help="print recent checks without taking one")
    ap.add_argument("--force", action="store_true", help="take a snapshot even if one is recent")
    args = ap.parse_args()
    state = load_state()

    if not args.show and state["checks"] and not args.force:
        try:
            last = dt.datetime.fromisoformat(state["checks"][-1]["when"])
            waited = (dt.datetime.now() - last).total_seconds() / 60
            if waited < MIN_SNAPSHOT_MINUTES:
                print(f"last check was {waited:.0f} min ago; a snapshot sooner than "
                      f"{MIN_SNAPSHOT_MINUTES} min says nothing. Use --force to override.")
                return
        except (ValueError, KeyError):
            pass

    if args.show:
        for check in state["checks"][-6:]:
            print(f"{check['when']}  " + ", ".join(f"{k}={v}" for k, v in
                                                   sorted(check["counters"].items())))
        return

    load1 = (shell("cat /proc/loadavg").split() or ["?"])[0]
    now = {"when": dt.datetime.now().isoformat(timespec="seconds"),
           "counters": counters(), "workers": live_workers(), "load1": load1}
    lines = verdicts(now, state["checks"])

    print(f"{now['when']}  load1 {load1}  {len(now['workers'])} worker process(es)")
    for key, verdict, detail in lines:
        print(f"  {verdict:8s} {key:28s} {detail}")
    stalled = [k for k, v, _d in lines if v == "STALLED"]
    if stalled:
        print(f"\nSTALLED: {', '.join(stalled)}. Check the process before assuming it is slow: a "
              f"futex spin\nburns a core exactly as hard as real work, and the number that tells "
              f"them apart is write_bytes.")

    state["checks"] = (state["checks"] + [now])[-200:]
    os.makedirs(os.path.dirname(STATE), exist_ok=True)
    json.dump(state, open(STATE, "w"), indent=1)
    with open(LOG, "a") as fh:
        for key, verdict, detail in lines:
            fh.write(f"{now['when']}\t{verdict}\t{key}\t{detail}\n")


if __name__ == "__main__":
    main()
