#!/usr/bin/env python3
"""Human immune-cell transcription (blood atlas)

How much each gene is transcribed in a human monocyte and in a myeloid dendritic cell

    level / kind : transcription / RNAseq
    provides     : monocyte_tpm, myeloid_dc_tpm
    coverage     : 19,135 reviewed proteins
    citation     : Human Protein Atlas blood atlas, consensus immune-cell RNA (rna_immune_cell.tsv)
    url          : https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip
    local path   : datasets/host/immune_cell_rna/rna_immune_cell.tsv.zip

Quirks that cost time once:
    The monocyte is the cell Toxoplasma is carried in and the one the host-transcriptome slot
    asks about, and NO GTEx tissue is a monocyte -- whole blood is a mixture dominated by
    neutrophil and erythroid genes, so that question could not be answered from the tissue atlas
    at all. nTPM rather than TPM, because the Atlas normalises across cell types so two columns
    can be compared, which is what a repertoire question does. Two of the eighteen cell types
    are carried: the rest are lymphocyte subsets no slot of this map names, and eighteen columns
    answering two questions is a denominator problem, not coverage.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_immune_cell_rna.py
"""
from _common import run

KEY = "host_immune_cell_rna"

if __name__ == "__main__":
    run(KEY)
