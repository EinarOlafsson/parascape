#!/usr/bin/env python3
"""STRING v12.0, T. gondii VEG, organism's own evidence

Modules (Louvain communities of 10 or more genes) of STRING pairs with Toxoplasma's own experimental or curated-database evidence

    level / kind : reference / interaction
    provides     : string_module_measured
    coverage     : 199 genes, 8 modules
    citation     : Szklarczyk D et al. Nucleic Acids Res 2023;51:D638-D646
    PMID         : 36370105
    url          : https://stringdb-downloads.org/download/protein.links.full.v12.0/432359.protein.links.full.v12.0.txt.gz
    local path   : parascape/data/interaction_labels/interaction_labels.parquet

Quirks that cost time once:
    Observed evidence only (channels `experiments` and `database`, nothing transferred). VEG
    locus numbers equal ME49 ones. Built by scripts/labels/build_interaction_labels.py; source
    files and SHA-256 in parascape/data/interaction_labels/sources.json (instruction 103).

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/string_tgon_measured.py
"""
from _common import run

KEY = "string_tgon_measured"

if __name__ == "__main__":
    run(KEY)
