#!/usr/bin/env python3
"""Anopheles midgut transcription through a Plasmodium infection

How much each mosquito gene is transcribed in the midgut at four points after an infectious blood meal

    level / kind : transcription / scRNAseq
    provides     : anopheles_midgut_36h_cpm, anopheles_midgut_48h_cpm, anopheles_midgut_4d_cpm, anopheles_midgut_7d_cpm
    coverage     : 12,574 of 12,614 AGAP genes
    citation     : GSE284537 (PMID 41125888), single-cell RNA of infected Anopheles gambiae midguts, dual-organism with Plasmodium falciparum
    PMID         : 41125888
    accession    : GSE284537
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284537
    local path   : datasets/transcription/scRNAseq/GSE284537/

Quirks that cost time once:
    A PSEUDOBULK OF SINGLE-CELL DATA, which is not the same thing as a bulk library: a value is
    the average of a cell mixture whose composition changes between timepoints, so a rise can be
    a gene transcribed harder or a cell type becoming commoner, and the two cannot be told apart
    here. The only genuinely midgut series in the archive -- 32 of 32 samples record tissue and
    source name as midgut, where every bulk candidate whose title named the midgut held whole
    bodies or salivary glands. Einar admitted it on 2026-11-13: the single-cell exclusion is
    about pipeline fit, not about the evidence. The deposit ships RAW unfiltered Cell Ranger
    matrices, so counts are summed over barcodes carrying at least 500 molecules rather than
    over every droplet; cells per sample run 355 to 71,918 and the share of molecules inside
    called barcodes runs 2.1 to 97.3 per cent, so the threshold changes the answer rather than
    rounding it. The dsGFP CONTROL arm only, four replicates per timepoint, counts per million
    within a sample then the mean of the replicates; the dsEcR arm is an ecdysone-receptor
    knockdown and answers a different question. The 4 d control arm is the weak one: three of
    its four replicates call 1,153, 1,243 and 1,592 cells and hold under a fifth of their
    molecules in them. Keyed on VectorBase AGAP gene ids, joined on the identifier itself, and
    40 space genes plus 1,271 deposit features do not meet -- those stay missing rather than
    zero. Redistributed with its accession and paper cited, per Einar on 2026-11-13.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_anopheles_midgut_transcriptome.py
"""
from _common import run

KEY = "host_anopheles_midgut_transcriptome"

if __name__ == "__main__":
    run(KEY)
