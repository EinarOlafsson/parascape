#!/usr/bin/env python3
"""Human bone marrow protein abundance (PaxDb, PXD020192)

How abundant each protein is in human bone marrow, in parts per million

    level / kind : translation / proteomics
    provides     : bone_marrow_protein_ppm
    coverage     : 2,059 proteins (93.8% of the 2,195 in the source file) — 9.8% of the 21,020-row host space
    citation     : Di Meo A et al., Proteomic Profiling of the Human Tissue and Biological Fluid Proteome, J Proteome Res 2021 (PMID 33107741), PRIDE PXD020192, via PaxDb v6.1 (https://pax-db.org)
    PMID         : 33107741
    accession    : PXD020192
    url          : https://pax-db.org/downloads/latest/datasets/9606/9606-PXD020192_bone_marrow.txt
    local path   : datasets/host/paxdb/

Quirks that cost time once:
    ONE autopsy study, label-free, NSAF: the shallowest of the host proteomes (2,059 proteins
    against 15,108 for liver) because the sample is a marrow aspirate dominated by a few
    abundant proteins, one of which holds 16.7% of the whole. Whole marrow is a mixture of
    erythroid, myeloid and stromal cells, so a value is not a cell type. Parts per million;
    joined by gene symbol, and this join keeps only 76.7% of the abundance mass, so the dropped
    symbols include abundant proteins and the shipped values do not sum to a million.
    Redistributed with its accession and paper cited.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_bone_marrow_proteome.py
"""
from _common import run

KEY = "host_bone_marrow_proteome"

if __name__ == "__main__":
    run(KEY)
