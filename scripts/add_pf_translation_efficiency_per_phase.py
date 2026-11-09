"""Translation efficiency per blood-stage phase, computed from columns already shipped (instruction 41).

The literature index proposed GSE58402 for the Plasmodium translation-efficiency slot. GSE58402 is **already in
the registry**, as `pf_riboseq`: its ribosome-footprint and mRNA densities per phase (ring, early and late
trophozoite, schizont, merozoite) have been in `pf_nodes.parquet` for months. The slot read empty only because no
column divided one by the other, which is what translational efficiency is.

So nothing is downloaded. Five derived columns are computed from the shipped pairs, plus the span across the
cycle: a gene whose efficiency changes by phase is the interesting case, and the span is how to find it.

Efficiency is log2(footprints / mRNA) in the SAME phase, so a phase a gene was not measured in stays missing
rather than borrowing another phase's mRNA. The source columns are densities (RPKM) and already length-normalised,
so the ratio needs no further correction.

    python scripts/add_pf_translation_efficiency_per_phase.py [--dry-run]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

PHASES = ("ring", "early_trophozoite", "late_trophozoite", "schizont", "merozoite")
TABLE = os.path.join(ROOT, "parascape", "data", "pf_nodes.parquet")


def efficiency(nodes: pd.DataFrame) -> pd.DataFrame:
    """One efficiency column per phase, and the span across the phases a gene was measured in."""
    out = {}
    for phase in PHASES:
        footprints, transcripts = f"riboseq_rpf_{phase}", f"riboseq_mrna_{phase}"
        if footprints not in nodes.columns or transcripts not in nodes.columns:
            continue
        rpf = pd.to_numeric(nodes[footprints], errors="coerce")
        mrna = pd.to_numeric(nodes[transcripts], errors="coerce")
        usable = (rpf > 0) & (mrna > 0)
        column = pd.Series(np.nan, index=nodes.index)
        column[usable] = np.log2(rpf[usable] / mrna[usable])
        out[f"te_{phase}_derived"] = column.round(4)
    frame = pd.DataFrame(out, index=nodes.index)
    if frame.empty:
        return frame
    # The span needs at least two phases to mean anything: one phase is a value, not a change.
    comparable = frame.notna().sum(axis=1) >= 2
    span = frame.max(axis=1) - frame.min(axis=1)
    frame["te_phase_span_derived"] = span.where(comparable).round(4)
    return frame


def main(argv=None) -> int:
    """Compute the columns and write them into the shipped Plasmodium table."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dry-run", action="store_true", help="report the coverage and write nothing")
    args = parser.parse_args(argv)
    nodes = pd.read_parquet(TABLE)
    frame = efficiency(nodes)
    if frame.empty:
        print("no riboseq footprint/mRNA pairs in the table: nothing to derive")
        return 1
    for column in frame.columns:
        print(f"  {column:34s} {int(frame[column].notna().sum()):5d} genes")
    if args.dry_run:
        print("dry run: nothing written")
        return 0
    for column in frame.columns:
        nodes[column] = frame[column].to_numpy()
    nodes.to_parquet(TABLE, index=False)
    print(f"wrote {len(frame.columns)} columns into {TABLE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
