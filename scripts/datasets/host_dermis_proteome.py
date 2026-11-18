#!/usr/bin/env python3
"""Human dermis protein abundance (PaxDb, PXD027093)

How abundant each protein is in the superficial and the deep dermis, in parts per million

    level / kind : translation / proteomics
    provides     : dermis_superficial_protein_ppm, dermis_deep_protein_ppm
    coverage     : 4,437 superficial and 4,090 deep proteins (4,723 in all) — 22.5% of the 21,020-row host space
    citation     : Li J et al., Spatially resolved proteomic map shows that extracellular matrix regulates epidermal growth, Nat Commun 2022 (PMID 35817779), PXD027093, via PaxDb v6.1 (https://pax-db.org)
    PMID         : 35817779
    accession    : PXD027093
    url          : https://pax-db.org/downloads/latest/datasets/9606/9606-PXD027093_deep_dermis.txt
    local path   : datasets/host/paxdb/

Quirks that cost time once:
    Two layers kept apart, because the sporozoite is deposited in the dermis and the layers
    differ: laser-microdissected human skin, data-independent acquisition, NSAF. The accession
    is not served by PRIDE (the PRIDE API returns 404) and was confirmed from the paper's own
    text through Europe PMC full-text search, where PXD027093 occurs once. The whole-skin
    integrated set is NOT used: it is skin, not dermis. Parts per million; joined by gene
    symbol, so a symbol no reviewed protein claims is dropped (the join keeps 82.3% and 83.8% of
    the abundance mass). Redistributed with its accession and paper cited.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_dermis_proteome.py
"""
from _common import run

KEY = "host_dermis_proteome"

if __name__ == "__main__":
    run(KEY)
