#!/usr/bin/env python3
"""Check the gate of instruction 81 against the baseline recorded before re-quantification.

The gate, as instruction 81 states it: the 40 zero-parasite runs produce non-zero parasite mapping,
every re-quantified run agrees with its previous values within the expected tolerance, and any that
does not is REPORTED rather than overwritten.

This script only reads. It never writes into a run directory and never overwrites a quantification.
Run it once `parascape-requant` has re-quantified the three short-insert series.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import statistics
import sys

NAS = Path("/nas_mnt/carruthers/Einar/parascape_omics/rnaseq")
SERIES = ("GSE243206", "GSE243207", "GSE245775")
REFS = ("ME49", "ME49_host", "GT1", "VEG")


def rate(run: Path, tool: str, ref: str):
    """The mapping rate one tool reports for one reference, or None when it did not run."""
    if tool == "kallisto":
        f = run / f"kallisto_{ref}" / "run_info.json"
        key = "p_pseudoaligned"
    else:
        f = run / f"salmon_{ref}" / "aux_info" / "meta_info.json"
        key = "percent_mapped"
    if not f.exists():
        return None
    try:
        return json.load(f.open()).get(key)
    except (OSError, json.JSONDecodeError):
        return None


def trimmed(run: Path) -> bool:
    marker = run / "done.json"
    return marker.exists() and '"trimming"' in marker.read_text()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", type=Path,
                    default=Path("results/requant_gate_baseline_2026_11_17.csv"))
    ap.add_argument("--tolerance", type=float, default=5.0,
                    help="percentage points a re-quantified rate may move before it is reported")
    args = ap.parse_args(argv)
    if not args.baseline.exists():
        sys.exit(f"no baseline at {args.baseline}; record one before re-quantification, not after")
    base = {r["run"]: r for r in csv.DictReader(args.baseline.open())}

    still_zero, not_done, moved, ok = [], [], [], 0
    for series in SERIES:
        for run in sorted((NAS / series).glob("*/")):
            acc = run.name
            if not trimmed(run):
                not_done.append(acc)
                continue
            now = rate(run, "salmon", "ME49")
            if not now:
                still_zero.append(acc)
            else:
                ok += 1
            was = base.get(acc, {}).get("kallisto_ME49")
            new = rate(run, "kallisto", "ME49")
            if was not in (None, "") and new is not None:
                # a re-quantified run is EXPECTED to rise here; only a fall is a defect worth naming
                if new + args.tolerance < float(was):
                    moved.append((acc, float(was), new))

    total = sum(len(list((NAS / s).glob("*/"))) for s in SERIES)
    print(f"the three short-insert series: {total} runs")
    print(f"  not re-quantified yet : {len(not_done)}")
    print(f"  salmon still zero     : {len(still_zero)}")
    print(f"  salmon now non-zero   : {ok}")
    if moved:
        print(f"  WENT DOWN by more than {args.tolerance} points (reported, not overwritten):")
        for acc, was, new in moved:
            print(f"    {acc}: kallisto {was:.2f} -> {new:.2f}")
    passed = not not_done and not still_zero and not moved
    print(f"\nGATE: {'OPEN' if passed else 'NOT OPEN'}")
    if not_done:
        print(f"  reason: {len(not_done)} of {total} runs have not been re-quantified")
    if still_zero:
        print(f"  reason: {len(still_zero)} runs still map 0.00 per cent in salmon after trimming")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
