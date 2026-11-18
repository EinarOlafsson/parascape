"""The dual perturb-seq data at single-cell resolution: what the parasite does when an effector goes.

Figure 10 of the ParaScape paper (instruction 101), rebuilt on the data rather than on its summary.
The supplementary table the first version used is pseudobulk differential expression per effector:
it says which host genes each effector moves and nothing about the parasite. GSE229505 itself holds,
FOR EACH CELL, the host transcriptome, the parasite transcriptome and the CRISPR guide, so it can
answer a question the table cannot:

    when an effector is knocked out, does the PARASITE change too, in the same cell?

That matters because every host effect attributed to an effector assumes the parasite is otherwise
unchanged. If knocking out an effector also moves the parasite transcriptome, part of the host
response is downstream of a sick parasite rather than of the missing effector.

HOW A CELL GETS A GENOTYPE. A cell is assigned to a guide only when one guide carries at least
`MIN_GUIDE_UMI` reads and at least `GUIDE_PURITY` of that cell's guide reads. A cell with two
guides at similar depth is two parasites or a doublet, and is dropped and counted.

THE COMPARISON IS WITHIN ONE SAMPLE. Cells from different 10x runs are never pooled: the lane is a
batch, and this is the same argument as figure 4, where the laboratory structured the corpus three
times more than the perturbation did.

Writes `results/dual_perturb_single_cell_<date>/`.
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

DATA = "/media/carruthers/mnt1/parascape_downloads/GSE229505"

#: A guide needs this many reads in a cell before the cell is called for it.
MIN_GUIDE_UMI = 3

#: And that guide must hold this share of the guide reads in the cell.
GUIDE_PURITY = 0.75

#: A guide with fewer cells than this is not tested.
MIN_CELLS = 25

#: A gene detected in fewer than this share of the sample is not tested.
MIN_DETECTED = 0.05


def read_sample(sample: str) -> tuple:
    """The matrix, the features and the barcodes of one 10x run, as CSR.

    Read by hand rather than through `scipy.io.mmread`, which materialises the whole coordinate
    list as float64 before building the matrix: these counts are small integers and the largest
    run here has 22,553 cells, so int32 triplets keep it to a few hundred megabytes.
    """
    base = os.path.join(DATA, f"GSE229505_{sample}")
    features = pd.read_csv(base + "_features.tsv.gz", sep="\t", header=None,
                           names=["id", "name", "kind"])
    barcodes = pd.read_csv(base + "_barcodes.tsv.gz", sep="\t", header=None, names=["barcode"])
    with gzip.open(base + "_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                n_features, n_cells, _nnz = (int(x) for x in line.split())
                break
    # Read by pandas, not by a Python loop: the genome-scale run has about a hundred million
    # non-zero entries, and three Python lists of that length is more memory than the machine has.
    triplets = pd.read_csv(base + "_matrix.mtx.gz", sep=r"\s+", comment="%", skiprows=1,
                           header=None, names=["i", "j", "v"],
                           dtype={"i": np.int32, "j": np.int32, "v": np.int32})
    matrix = sp.csr_matrix((triplets["v"].to_numpy(),
                            (triplets["i"].to_numpy() - 1, triplets["j"].to_numpy() - 1)),
                           shape=(n_features, n_cells))
    del triplets
    return matrix, features, barcodes


def call_guides(matrix, features) -> tuple:
    """One guide per cell where the evidence allows, and a count of the cells that got none."""
    is_guide = (features["kind"] == "CRISPR Guide Capture").to_numpy()
    guides = matrix[is_guide]
    total = np.asarray(guides.sum(axis=0)).ravel()
    dense = guides.toarray() if guides.shape[0] < 4000 else None
    if dense is None:
        dense = guides.toarray()
    best = dense.argmax(axis=0)
    top = dense[best, np.arange(dense.shape[1])]
    with np.errstate(invalid="ignore", divide="ignore"):
        purity = np.where(total > 0, top / np.maximum(total, 1), 0.0)
    called = (top >= MIN_GUIDE_UMI) & (purity >= GUIDE_PURITY)
    names = features.loc[is_guide, "name"].to_numpy()
    assignment = np.where(called, names[best], None)
    return assignment, int((~called).sum())


def target_of(guide: str) -> str:
    """The gene a guide targets, from its name: TGME49_068760_12_2152 -> TGME49_068760."""
    parts = str(guide).split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else str(guide)


def main() -> None:
    sample = os.environ.get("GSE229505_SAMPLE", "BUT290A1_BUT290A2")
    print(f"reading {sample} ...", flush=True)
    matrix, features, barcodes = read_sample(sample)
    print(f"  {matrix.shape[0]:,} features x {matrix.shape[1]:,} cells")

    assignment, unassigned = call_guides(matrix, features)
    called = pd.Series(assignment)
    print(f"  {int(called.notna().sum()):,} cells carry one guide, {unassigned:,} do not "
          f"({called.notna().mean():.1%} assigned)")
    if called.notna().sum() < MIN_CELLS * 2:
        raise SystemExit("too few assigned cells in this sample to compare anything")

    host = features["id"].str.startswith("GRCh38_").to_numpy()
    parasite = (~host) & (features["kind"] == "Gene Expression").to_numpy()
    host_total = np.asarray(matrix[host].sum(axis=0)).ravel()
    parasite_total = np.asarray(matrix[parasite].sum(axis=0)).ravel()
    print(f"  median reads per cell: host {np.median(host_total):,.0f}, "
          f"parasite {np.median(parasite_total):,.0f}")

    targets = called.map(lambda g: target_of(g) if isinstance(g, str) else None)
    counts = targets.value_counts()
    testable = counts[counts >= MIN_CELLS]
    print(f"  {len(counts)} targets; {len(testable)} with at least {MIN_CELLS} cells")
    if testable.empty:
        raise SystemExit("no target reaches the minimum number of cells in this sample")

    # THE QUESTION. For each target with enough cells, is the parasite's own transcriptional
    # output different in those cells? Compared with all other assigned cells of the same lane,
    # so the batch is held fixed.
    assigned = targets.notna().to_numpy()
    # Where the experiment carries a neutral control, compare against IT rather than against
    # "everything else". UPRT is the standard neutral target in Toxoplasma, and in a run with two
    # guides "each against the rest" is one comparison printed twice with the signs swapped.
    control = next((t for t in testable.index if t.upper().startswith("UPRT")), None)
    if control:
        print(f"  neutral control present: {control}; every target is compared against it")
    rows = []
    for target in testable.index:
        if target == control:
            continue
        mine = (targets == target).to_numpy()
        others = ((targets == control).to_numpy() if control else assigned & ~mine)
        if mine.sum() < MIN_CELLS or others.sum() < MIN_CELLS:
            continue
        p = mannwhitneyu(parasite_total[mine], parasite_total[others], alternative="two-sided")
        h = mannwhitneyu(host_total[mine], host_total[others], alternative="two-sided")
        rows.append({"target": target, "cells": int(mine.sum()),
                     "parasite_reads_median": float(np.median(parasite_total[mine])),
                     "parasite_reads_median_others": float(np.median(parasite_total[others])),
                     "parasite_p": float(p.pvalue),
                     "host_reads_median": float(np.median(host_total[mine])),
                     "host_reads_median_others": float(np.median(host_total[others])),
                     "host_p": float(h.pvalue),
                     "compared_with": control or "every other assigned cell"})
    table = pd.DataFrame(rows)
    if table.empty:
        raise SystemExit("no target had enough cells on both sides to compare")
    for column in ("parasite_p", "host_p"):
        table[column.replace("_p", "_q")] = (table[column] * len(table)).clip(upper=1.0)
    moved = table[table["parasite_q"] < 0.05]
    print(f"\n{len(moved)} of {len(table)} targets change the PARASITE transcriptional output of "
          f"the cells they are in (q < 0.05).")
    if len(moved):
        for row in moved.sort_values("parasite_q").head(8).itertuples():
            direction = "lower" if row.parasite_reads_median < row.parasite_reads_median_others \
                else "higher"
            print(f"   {row.target:18s} {row.cells:5,} cells  parasite reads {direction}: "
                  f"{row.parasite_reads_median:,.0f} against "
                  f"{row.parasite_reads_median_others:,.0f}  q = {row.parasite_q:.2g}")
    print(f"\nWHAT THIS TESTS, AND WHAT IT DOES NOT. The quantity is each cell's TOTAL "
          f"transcriptional\noutput, host and parasite. It answers 'is the parasite sick when this "
          f"effector is gone',\nwhich is the assumption every pseudobulk host analysis rests on. "
          f"It does NOT test the\nCOMPOSITION of either transcriptome: an effector could change "
          f"which genes are expressed\nwithout changing how many reads there are, and this would "
          f"not see it.")
    if len(table) >= 40:
        print(f"\nAt this scale the assumption largely holds: {len(moved)} of {len(table)} targets "
              f"move the parasite\nand "
              f"{int((table['host_q'] < 0.05).sum())} move the host, which is about what "
              f"correction leaves when there is\nlittle or nothing there.")

    out = os.path.join(ROOT, "results", f"dual_perturb_single_cell_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    table.to_csv(os.path.join(out, f"{sample}_targets.tsv"), sep="\t", index=False)
    json.dump({"sample": sample, "cells": int(matrix.shape[1]),
               "cells_assigned": int(called.notna().sum()), "cells_unassigned": unassigned,
               "min_guide_umi": MIN_GUIDE_UMI, "guide_purity": GUIDE_PURITY,
               "targets": int(len(counts)), "targets_tested": int(len(table)),
               "targets_moving_the_parasite": int(len(moved)),
               "neutral_control": control,
               "median_host_reads": float(np.median(host_total)),
               "median_parasite_reads": float(np.median(parasite_total))},
              open(os.path.join(out, f"{sample}_summary.json"), "w"), indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
