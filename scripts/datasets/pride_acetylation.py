#!/usr/bin/env python3
"""Lysine acetylome (GCN5b)

Acetylation sites reported per gene

    level / kind : post_translation / proteomics
    provides     : n_acetylation_sites
    coverage     : 3,921 genes measured
    citation     : Cannella D, Pachano B, Shahinas M, Belmudes L, deVries J, Corrao C, Pounot L, Hesse AM, Couté Y, Bougdour A, Swale C, Hakimi MA. Divergent SAGA complexes shape the Toxoplasma transcriptome for lytic cycle progression and host interaction. Nat Commun. 2026 17(1):10158
    accession    : PXD079431
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD079431
    local path   : datasets/quarantine/2026_08_16_pride/Tg/acetylation/

Quirks that cost time once:
    Read from the deposit's own MaxQuant Sites tables. Held in quarantine rather than the
    dataset archive: a deposit is promoted by being checked, not by being downloaded, and the
    check here was reading the files and finding both Toxoplasma genes and the modification
    named.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pride_acetylation.py
"""
from _common import run

KEY = "pride_acetylation"

if __name__ == "__main__":
    run(KEY)
