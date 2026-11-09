#!/usr/bin/env python3
"""Liver-stage transcription, late (day 6)

How much each gene is transcribed in the late liver stage, in a humanised mouse

    level / kind : transcription / RNAseq
    provides     : expr_pf_liver_day6
    coverage     : 5,720 genes
    citation     : GSE220039, Plasmodium falciparum and P. vivax liver stages in humanised mice
    PMID         : 41428747
    accession    : GSE220039
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220039
    local path   : datasets/transcription/RNAseq/GSE220039/GSE220039_PfHsMmu_GRCh38_GRCm39_Pf3D7v58_raw_counts.csv.gz

Quirks that cost time once:
    Reads mapped jointly to parasite, human and mouse; only the PF3D7 rows are kept and
    renormalised to the PARASITE total, because the parasite is a small minority of a liver
    biopsy and a joint CPM would describe hepatocytes. ONLY DAY 6 IS SHIPPED, and the reason is
    the measurement: the series holds 408 parasite reads on day 0, 27,260 on day 2 and 71,738 on
    day 4, against 9,165,928 on day 6 -- across 5,720 genes that is twelve reads a gene on day 4
    and under five on day 2, so an early-day column would be noise wearing a number. Day 6 is
    the late liver stage before merozoite release, which is what the slot asks for. The P. vivax
    arm of the same series is not read here: there is no shipped vivax table to put it in.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_liver_transcriptome.py
"""
from _common import run

KEY = "pf_liver_transcriptome"

if __name__ == "__main__":
    run(KEY)
