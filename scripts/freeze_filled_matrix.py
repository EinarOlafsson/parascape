"""Freeze and date the filled cells of the information map (instruction 73, phase 6).

Phase 6 asks for the filled matrix to be frozen and dated so papers published afterwards can test it. The
filled cells are the only part of the map that is a PREDICTION -- every other cell is a measurement, and a
measurement needs no prospective test. So the freeze is of those cells, written with the question that produced
each one, the calibrated probability, the partitions its question cleared, and a hash.

Grading later means: verify the hash, find the genes among them that a later paper measures for that trait, and
report how many filled cells were right. The fills claim 90 % precision at their lower bound, so the test is
whether at least that share survives.

    python scripts/freeze_filled_matrix.py
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
MAP = os.path.join(ROOT, "parascape", "data", "information_map")
STABILITY = os.path.join(ROOT, "results", "information_map_fills_2026_11_08", "partition_stability.tsv")
OUT = os.path.join(ROOT, "results", f"filled_matrix_freeze_{_dt.date.today().isoformat().replace('-', '_')}")


def main() -> int:
    """Write the frozen filled cells, their provenance and the manifest."""
    import pandas as pd

    from parascape import __version__, information_map as IM

    os.makedirs(OUT, exist_ok=True)
    published = IM.load_published(MAP)
    cells = published.predicted_cells().copy()
    if cells.empty:
        print("the map holds no filled cells: nothing to freeze")
        return 1
    # Each fill carries the question that produced it in its source tier; the partition record says which
    # partitions that question cleared, which is what a later reader needs to judge the claim.
    stability = pd.read_csv(STABILITY, sep="\t") if os.path.exists(STABILITY) else pd.DataFrame()
    if not stability.empty:
        passed = stability[stability.may_fill].groupby("question").seed.apply(
            lambda s: ", ".join(str(x) for x in sorted(s)))
        bounds = stability.groupby("question").lower90.apply(lambda s: ", ".join(f"{x:.3f}" for x in s))
        question = cells.source_tier.astype(str).str.split().str[-1]
        cells["question"] = question
        cells["partitions_cleared"] = question.map(passed).fillna("not recorded")
        cells["lower_bounds"] = question.map(bounds).fillna("not recorded")
    path = os.path.join(OUT, "filled_cells.tsv")
    cells.sort_values(["column", "gene_id"]).to_csv(path, sep="\t", index=False)
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    manifest = dict(
        frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S"), parascape_version=__version__,
        filled_cells=int(len(cells)), columns=sorted(cells.column.unique()),
        genes=int(cells.gene_id.nunique()), filled_cells_sha256=digest,
        bar="the fills claim 90 % precision at the frozen-test lower bound, on at least 2 of 3 partitions "
            "(Einar's rule, 2026-11-08)",
        grade_with="verify the hash, then for the genes a later paper measures for that trait, report the share "
                   "of filled cells that were right; the claim fails if it falls below 90 %",
        note="Only the PREDICTED cells are frozen. Every other cell of the map is a measurement, and a "
             "measurement needs no prospective test.")
    json.dump(manifest, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    print(f"froze {len(cells)} filled cells over {len(manifest['columns'])} column(s) -> {OUT}")
    print(f"sha256 {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
