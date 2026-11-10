#!/usr/bin/env python3
"""Mouse tissue protein abundance (PaxDb)

How abundant each protein is in mouse brain, macrophage or muscle, in parts per million

    level / kind : translation / proteomics
    provides     : mouse_brain_protein_ppm, bmdm_protein_ppm, skeletal_muscle_protein_ppm
    coverage     : 4,603 to 10,347 proteins per tissue
    citation     : PaxDb v6.1 integrated datasets (https://pax-db.org)
    url          : https://pax-db.org/downloads/latest/datasets/10090/
    local path   : datasets/host/paxdb/

Quirks that cost time once:
    The three tissues a Toxoplasma question asks about in the mouse: brain, where the cyst sits;
    macrophage, the cell that carries the parasite; muscle, the other tissue cysts are found in.
    Same measure and same symbol join as the human sets, against the mouse index -- a mouse
    protein is not a human one and the two host tables stay separate. Coverage is reported per
    column because the symbol join drops what it cannot resolve. The names differ from the human
    set on purpose: a host slot resolves to ONE organism, and the generator refuses a column
    name two species claim.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_mouse_tissue_proteome.py
"""
from _common import run

KEY = "host_mouse_tissue_proteome"

if __name__ == "__main__":
    run(KEY)
