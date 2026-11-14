"""The Anopheles midgut host transcriptome of GSE284537, pseudo-bulked per sample.

Instruction 41. `Pf_host transcriptome · Anopheles midgut` had no source. Every bulk candidate whose
title named the midgut turned out to hold whole bodies or salivary glands, and the one genuinely
midgut series is single-cell: 32 of 32 samples record `tissue: midgut` and `source name: midgut`,
dual-organism *Anopheles gambiae* and *Plasmodium falciparum*, 36 h / 2 d / 4 d / 7 d after an
infectious blood meal, four biological replicates, dsGFP-injected controls beside dsEcR knockdowns.
Einar admitted it as a per-sample pseudobulk on 2026-11-13: the exclusion recorded in
`rnaseq_excluded_single_cell` is about pipeline fit, not about the evidence.

**What the archive holds.** `GSE284537_RAW.tar` is 2,800,435,200 bytes and contains 192 files: a
Cell Ranger 4 triplet (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`) per sample per
organism, named `_AnG_` and `_Pf_`. There is no processed gene table and no combined matrix. The
matrices are the RAW, unfiltered barcode matrices -- 0.9 to 3.7 million barcodes each -- so a
barcode here is a droplet and not yet a cell.

**How host and parasite are separated.** By the reference the submitters quantified against, read
off the features file rather than off the file name: every `_AnG_` features file is the same 13,845
AGAP genes (one checksum across all 32), and every `_Pf_` features file is 3,000-odd PF3D7 genes.
Only the `_AnG_` matrices are read here.

**Cells, not barcodes.** Counts are summed over barcodes carrying at least `MIN_UMI` (500) unique
molecules. Summing every barcode would sum the ambient RNA of two to three million empty droplets
with the cells, and the share of molecules inside called barcodes runs from 2.1 to 95.2 per cent
across these samples, so the choice changes the answer rather than rounding it.

**The slot column, and what is left beside it.** Only the **dsGFP control** arm becomes the slot:
the dsEcR samples are an ecdysone-receptor knockdown, which is a perturbation and a different
question, and they are written to the deposit as their own columns for a perturbation slot to claim.
Counts per million within each sample, then the mean of the four replicates -- the same policy
`host_bmdm_baseline` uses, so two host transcriptomes in this map are the same quantity.

**Identifiers.** The deposit is keyed on the AGAP ids of the Cell Ranger reference the submitters
built (13,845 features); the Anopheles gene space is the Ensembl Metazoa AgamP4 protein-coding build
(12,614 genes). The join is on the identifier itself -- no mapping table, no sequence work -- and it
reaches 12,574 of 12,614 space genes. The 40 it misses, and the 1,271 deposit features that are not
protein-coding genes of this build, stay missing rather than zero.

Writes `host_transcriptome_anopheles_midgut_GSE284537.csv` (the keyed table: the four control
columns and the 32 per-sample counts) and `host_transcriptome_anopheles_midgut_GSE284537_cells.csv`
(barcodes, called cells and molecules per sample) to the dataset root.
`parascape.deposits.anopheles_midgut_transcriptome` reads the first.

    python scripts/build_anopheles_midgut_host_deposit.py
"""
from __future__ import annotations

import glob
import gzip
import hashlib
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

SRC = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"
FOLDER = os.path.join("transcription", "scRNAseq", "GSE284537")
OUT = "host_transcriptome_anopheles_midgut_GSE284537.csv"
CELLS = "host_transcriptome_anopheles_midgut_GSE284537_cells.csv"
#: The Anopheles gene space, built by `build_host_gene_spaces_notebook.py` (instruction 53, WP5).
SPACE = os.path.join(ROOT, "results", "host_gene_spaces_2026_11_09", "ag_gene_nodes.parquet")

#: The control arm, per timepoint, as the deposit names them. `48h` is the series' "2 days".
TIMEPOINTS = ("36h", "48h", "4d", "7d")
CONTROL = "GFP"

#: A barcode with fewer molecules than this is treated as an empty droplet rather than a cell. The
#: deposit ships the raw barcode matrix, so some threshold has to be chosen; 500 is the usual floor
#: below which a 10x barcode is ambient. It is recorded here, not buried, because a different
#: threshold is a different column.
MIN_UMI = 500

#: `GSM8686106_AnG_r1-36h-GFP-NS_matrix.mtx.gz` -> replicate, timepoint, arm.
NAME = re.compile(r"^(GSM\d+)_AnG_(r\d)-(\w+?)-(GFP|EcR)-NS_matrix\.mtx\.gz$")


def pseudobulk(path: str, n_genes: int):
    """Sum one sample over its called barcodes. Returns the gene vector and the quality counts."""
    with gzip.open(path, "rt") as handle:
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        genes, barcodes, _nnz = (int(v) for v in line.split())
        if genes != n_genes:
            raise ValueError(f"{os.path.basename(path)}: {genes} features, expected {n_genes}")
        frame = pd.read_csv(handle, sep=r"\s+", header=None, names=["gene", "cell", "count"],
                            dtype=np.int64)
    per_barcode = np.bincount(frame["cell"].to_numpy(), weights=frame["count"].to_numpy(),
                              minlength=barcodes + 1)
    called = per_barcode >= MIN_UMI
    keep = called[frame["cell"].to_numpy()]
    summed = np.bincount(frame.loc[keep, "gene"].to_numpy(),
                         weights=frame.loc[keep, "count"].to_numpy(), minlength=n_genes + 1)
    return summed[1:].astype(np.int64), dict(
        barcodes=int(barcodes), cells=int(called.sum()), molecules=int(per_barcode.sum()),
        molecules_in_cells=int(per_barcode[called].sum()),
        genes_detected=int((summed[1:] > 0).sum()))


def main() -> int:
    folder = os.path.join(SRC, FOLDER)
    matrices = sorted(glob.glob(os.path.join(folder, "*_AnG_*_matrix.mtx.gz")))
    if not matrices:
        print(f"no Anopheles matrices under {folder}")
        return 1
    print(f"{len(matrices)} Anopheles samples under {folder}")

    features = matrices[0].replace("matrix.mtx.gz", "features.tsv.gz")
    with gzip.open(features, "rt") as handle:
        ids = [line.split("\t")[0] for line in handle.read().splitlines()]
    print(f"{len(ids):,} AGAP features, first {ids[0]}, last {ids[-1]}")

    columns, quality = {"gene_id": ids}, []
    for matrix in matrices:
        name = NAME.match(os.path.basename(matrix))
        if name is None:
            raise ValueError(f"unexpected matrix name: {os.path.basename(matrix)}")
        gsm, replicate, timepoint, arm = name.groups()
        own = matrix.replace("matrix.mtx.gz", "features.tsv.gz")
        with gzip.open(own, "rt") as handle:
            if [line.split("\t")[0] for line in handle.read().splitlines()] != ids:
                raise ValueError(f"{gsm}: features differ from the first sample")
        values, counts = pseudobulk(matrix, len(ids))
        label = f"{replicate}_{timepoint}_{arm}"
        columns[label] = values
        digest = hashlib.sha256()
        with open(matrix, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        quality.append(dict(sample=label, gsm=gsm, replicate=replicate, timepoint=timepoint,
                            arm=arm, matrix_sha256=digest.hexdigest(), **counts))
        share = 100 * counts["molecules_in_cells"] / counts["molecules"]
        print(f"  {label:16s} {gsm}  {counts['barcodes']:>9,} barcodes, "
              f"{counts['cells']:>7,} cells, {counts['molecules_in_cells']:>12,} molecules in "
              f"cells ({share:.1f} per cent of the library), {counts['genes_detected']:,} genes")

    counted = pd.DataFrame(columns).set_index("gene_id")
    cells = pd.DataFrame(quality)
    shared = cells.groupby("matrix_sha256")["sample"].apply(list)
    for samples in shared[shared.str.len() > 1]:
        print(f"  NOTE: {' and '.join(samples)} carry byte-identical matrices in the deposit")

    # Counts per million inside each sample, then the mean of the four replicates of one timepoint.
    # A mean of raw counts would be a mean of sequencing effort: these libraries run from 11 million
    # molecules to 324 million.
    per_million = counted / counted.sum(axis=0) * 1e6
    control = {}
    for timepoint in TIMEPOINTS:
        arms = [c for c in counted.columns if c.endswith(f"_{timepoint}_{CONTROL}")]
        control[f"anopheles_midgut_{timepoint}_cpm"] = per_million[arms].mean(axis=1).round(4)

    space = pd.read_parquet(SPACE, columns=["gene_id", "gene_name"])
    out = space.copy()
    for name, values in control.items():
        out[name] = out.gene_id.map(values)
    for sample in counted.columns:
        out["sample_" + sample] = out.gene_id.map(counted[sample]).astype("Int64")
    reached = int(out[list(control)].notna().all(axis=1).sum())
    out.to_csv(os.path.join(SRC, OUT), index=False)
    cells.to_csv(os.path.join(SRC, CELLS), index=False)
    print(f"wrote {os.path.join(SRC, OUT)}: {len(out):,} space genes x {out.shape[1] - 2} columns")
    print(f"  {reached:,} of {len(out):,} space genes carry a value "
          f"({100 * reached / len(out):.1f} per cent); "
          f"{len(set(counted.index) - set(out.gene_id)):,} deposit features are not "
          f"protein-coding genes of this build")
    print(f"wrote {os.path.join(SRC, CELLS)}")
    print(f"  cells per sample: min {cells.cells.min():,}, median "
          f"{int(cells.cells.median()):,}, max {cells.cells.max():,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
