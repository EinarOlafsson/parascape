#!/usr/bin/env python3
"""Human skin fibroblast protein abundance (PaxDb, PXD032893)

How abundant each protein is in untreated human skin fibroblasts, in parts per million

    level / kind : translation / proteomics
    provides     : fibroblast_protein_ppm
    coverage     : 6,390 proteins (95.3% of the 6,703 in the source file) — 30.4% of the 21,020-row host space
    citation     : Martinez-Martinez E et al., Int J Mol Sci 2022 (PMID 35743299), PRIDE PXD032893, via PaxDb v6.1 (https://pax-db.org)
    PMID         : 35743299
    accession    : PXD032893
    url          : https://pax-db.org/downloads/latest/datasets/9606/9606-SKIN_FIBROBLAST-integrated.txt
    local path   : datasets/host/paxdb/

Quirks that cost time once:
    PaxDb calls this set integrated, but it integrates the two VEHICLE arms (DMSO, high and low
    glucose) of ONE study, so it is one deposit and not a consensus. The cells are normal human
    PRIMARY SKIN fibroblasts, serum-starved (0.25% FBS) for 24 h with vehicle, not the foreskin
    line (HFF) Toxoplasma is grown in: the slot asks about a human fibroblast and this is the
    closest unperturbed proteome that exists, which is a judgement and not an identity. NSAF-
    based, in parts per million; joined by gene symbol, so a symbol no reviewed protein claims
    is dropped (the join keeps 91.3% of the abundance mass). Redistributed with its accession
    and paper cited.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_fibroblast_proteome.py
"""
from _common import run

KEY = "host_fibroblast_proteome"

if __name__ == "__main__":
    run(KEY)
