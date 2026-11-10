#!/usr/bin/env python3
"""Human tissue protein abundance (PaxDb)

How abundant each protein is in a human tissue, in parts per million of the proteome

    level / kind : translation / proteomics
    provides     : astrocyte_protein_ppm, brain_protein_ppm, liver_protein_ppm
    coverage     : 9,930 to 15,108 proteins per tissue
    citation     : PaxDb v6.1 integrated datasets (https://pax-db.org)
    url          : https://pax-db.org/downloads/latest/datasets/9606/
    local path   : datasets/host/paxdb/

Quirks that cost time once:
    Parts per million of the proteome, which is what makes two tissues comparable: a raw
    intensity is a property of the run, a ppm is a property of the sample. The INTEGRATED sets
    are used, a weighted average over the studies PaxDb holds for that tissue, rather than any
    single deposit. Joined by gene SYMBOL, which is the weaker join and is why coverage is
    reported per column: PaxDb identifies proteins by Ensembl PROTEIN id and no ENSP-to-gene map
    is on disk, so a symbol no reviewed protein claims is dropped rather than guessed.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_tissue_proteome.py
"""
from _common import run

KEY = "host_tissue_proteome"

if __name__ == "__main__":
    run(KEY)
