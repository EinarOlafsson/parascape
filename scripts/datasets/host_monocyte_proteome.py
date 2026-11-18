#!/usr/bin/env python3
"""Human monocyte protein abundance (PaxDb, Kim 2014 draft map)

How abundant each protein is in adult human monocytes, in parts per million

    level / kind : translation / proteomics
    provides     : monocyte_protein_ppm
    coverage     : 4,989 proteins (94.8% of the 5,265 in the source file)
    citation     : Kim MS et al., A draft map of the human proteome, Nature 2014 (PMID 24870542, doi 10.1038/nature13302), PRIDE PXD000561, via PaxDb v6.1 (https://pax-db.org)
    PMID         : 24870542
    accession    : PXD000561
    url          : https://pax-db.org/downloads/latest/datasets/9606/9606-Adult_Monocytes_Kim_2014_SEQUEST.txt
    local path   : datasets/host/paxdb/

Quirks that cost time once:
    ONE study, spectral counting, not an integrated set: a value is a count-based estimate and
    the deep tail rests on few spectra. Monocytes from healthy adult donors, not infected and
    not sorted into subsets. Parts per million; joined by gene symbol, so a symbol no reviewed
    protein claims is dropped (the join keeps 95.2% of the abundance mass). Redistributed with
    its accession and paper cited.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_monocyte_proteome.py
"""
from _common import run

KEY = "host_monocyte_proteome"

if __name__ == "__main__":
    run(KEY)
