#!/usr/bin/env python3
"""STRING v12.0, T. gondii VEG, evidence transferred from orthologs

Modules of STRING pairs whose transferred experimental or database score is at least 700

    level / kind : reference / orthology
    provides     : string_module_orthology
    coverage     : 2,469 genes, 14 modules
    citation     : Szklarczyk D et al. Nucleic Acids Res 2023;51:D638-D646
    PMID         : 36370105
    url          : https://stringdb-downloads.org/download/protein.links.full.v12.0/432359.protein.links.full.v12.0.txt.gz
    local path   : parascape/data/interaction_labels/interaction_labels.parquet

Quirks that cost time once:
    COPIED from other organisms: a label, never a measurement of Toxoplasma. STRING co-
    expression and text mining are not used (co-expression is circular with the transcript
    features).

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/string_tgon_transferred.py
"""
from _common import run

KEY = "string_tgon_transferred"

if __name__ == "__main__":
    run(KEY)
