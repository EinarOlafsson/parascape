"""The oocyst transcriptome of GSE271211, keyed on CryptoDB `cgd` identifiers.

Instruction 41. `Cp_transcription · oocyst` had no source: of the 63 Cryptosporidium slots only nine
graded, all from the gene table itself, and the audit of instruction 65 reported zero Cp sources.

**Why this deposit and not GSE232438**, which an earlier pass accepted for this slot. GSE232438 is a
genuine lifecycle series, but no SAMPLE is an oocyst: every mention of the word sits in the growth
protocol -- "infected with 250,000 oocysts" -- where the oocyst is what was put on the cells. Here the
oocyst is the thing measured. Three `Inoculum` samples at time zero are undigested oocysts in
triplicate, and twelve more follow them through simulated gastric and ileal digestion, child and
adult.

**Only the three inoculum samples become the slot.** The digestion arms are oocysts under a
perturbation -- the excystation trigger -- so averaging them into an unperturbed slot would be the
same mistake instruction 65 found in the Toxoplasma table, where nine AP2XII-1-depleted columns sat
in a slot declaring itself unperturbed. They are written to the deposit as their own columns and left
for a perturbation slot to claim.

**The identifier problem, and what it costs.** The counts are keyed on `CPATCC_` locus tags from
GCA_015245375.1 (IOWA-ATCC, University of Georgia, 2022), while the gene table is keyed on `cgd`
from GCA_000165345.1 (Iowa II). Two assemblies, two independent annotations, no shared identifier
and no NCBI Gene record to join on. `build_cpatcc_cgd_map.py` maps 2,571 genes one-to-one by exact
protein identity -- 66 per cent of the ATCC annotation. **The remaining third is left unmapped rather
than filled by proximity**, so this column covers about two thirds of the genome and says so.

Writes `stage_measured_cryptosporidium_oocyst_GSE271211.csv` to the dataset root.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

SRC = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"
MAP = os.path.join(ROOT, "results", "cryptosporidium_space_2026_11_08", "cpatcc_to_cgd.tsv")
COUNTS = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271211/suppl/"
          "GSE271211_ReadCounts.txt.gz")
OUT = os.path.join(SRC, "stage_measured_cryptosporidium_oocyst_GSE271211.csv")

# The submitter names the undigested inoculum this way; everything else in the series is an oocyst
# that has been through a digestion compartment.
INOCULUM = "Inoculum"


def main() -> int:
    with urllib.request.urlopen(COUNTS, timeout=300) as fh:
        raw = fh.read()
    print(f"counts: {len(raw):,} bytes, sha256 {hashlib.sha256(raw).hexdigest()[:16]}")
    counts = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep="\t")
    counts = counts.rename(columns={counts.columns[0]: "cpatcc_locus_tag"})
    print(f"{len(counts):,} genes x {counts.shape[1] - 1} samples")

    mapping = pd.read_csv(MAP, sep="\t")
    merged = counts.merge(mapping[["cpatcc_locus_tag", "cgd_gene_id"]], on="cpatcc_locus_tag")
    print(f"{len(merged):,} rows carry a one-to-one cgd identifier "
          f"({100 * len(merged) / len(counts):.1f} per cent of the table)")

    samples = [c for c in counts.columns if c != "cpatcc_locus_tag"]
    inoculum = [c for c in samples if c.startswith(INOCULUM)]
    digested = [c for c in samples if c not in inoculum]
    assert len(inoculum) == 3, inoculum
    print(f"inoculum (undigested oocysts): {inoculum}")

    # Counts per million within each sample, then log2, then the mean over the three replicates.
    # Depth differs between samples and a mean of raw counts would be a mean of sequencing effort.
    values = merged[samples].apply(pd.to_numeric, errors="coerce")
    cpm = values.div(values.sum(axis=0), axis=1) * 1e6
    logged = np.log2(cpm + 1)

    out = pd.DataFrame({"gene_id": merged.cgd_gene_id})
    out["oocyst_log2cpm"] = logged[inoculum].mean(axis=1)
    out["oocyst_detected"] = (values[inoculum] > 0).any(axis=1)
    for column in digested:
        out["digested_" + column.lower()] = logged[column]
    out = out.groupby("gene_id", as_index=False).mean(numeric_only=False)

    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}")
    print(f"  {len(out):,} genes, {int(out.oocyst_detected.sum()):,} detected in at least one "
          f"undigested oocyst replicate")
    print(f"  oocyst_log2cpm: median {out.oocyst_log2cpm.median():.2f}, "
          f"max {out.oocyst_log2cpm.max():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
