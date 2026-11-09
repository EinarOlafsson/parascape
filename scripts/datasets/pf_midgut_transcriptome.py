#!/usr/bin/env python3
"""Mosquito-midgut transcription, ookinete to oocyst

How much each gene is transcribed in the mosquito midgut, at four points after the blood meal

    level / kind : transcription / scRNAseq
    provides     : expr_midgut_36h, expr_midgut_48h, expr_midgut_4d, expr_midgut_7d
    coverage     : 5,720 genes at four timepoints
    citation     : GSE284537, single-cell RNA of infected Anopheles midguts
    accession    : GSE284537
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284537
    local path   : datasets/transcription/scRNAseq/GSE284537/

Quirks that cost time once:
    Single-cell matrices summed over every cell, replicates added, log2(CPM + 1) per timepoint.
    Only the dsGFP CONTROL arms: the series pairs each timepoint with an ecdysone-receptor
    knockdown, which is a perturbation and a different question, and its Anopheles arm is not a
    Plasmodium table. 36 h is the late ookinete and 48 h to 7 d the developing oocyst -- there
    is NO sporozoite timepoint, so the mosquito-stage slot is filled for two stages of three.
    Depth differs sharply by timepoint (371,390 summed counts at 48 h against 3,866,156 at 7 d),
    and a zero is a non-detection at that depth rather than a measured absence: 4,281 of 5,720
    genes are seen at all at 36 h against 5,591 at 7 d. Einar took the call to use it on
    2026-11-09 knowing it covers part of the slot.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_midgut_transcriptome.py
"""
from _common import run

KEY = "pf_midgut_transcriptome"

if __name__ == "__main__":
    run(KEY)
