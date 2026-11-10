#!/usr/bin/env python3
"""Translation efficiency, asexual blood stage

How much each transcript is translated, as ribosome footprints over its own mRNA

    level / kind : translation / riboseq
    provides     : pf_ribosome_footprint_cpm
    coverage     : 5,720 genes
    citation     : GSE286083, ribosome profiling with matched RNA in P. falciparum trophozoites
    PMID         : 41099697
    accession    : GSE286083
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286083
    local path   : datasets/translation/riboseq/GSE286083/GSE286083_raw_counts.txt.gz

Quirks that cost time once:
    Only the DMSO arms are used. The series also holds an i14G1-12 inhibitor arm, and a drug-
    treated sample answers a different question from how well a transcript is translated.
    Cycloheximide footprints, because the harringtonine arms have no matched RNA and so cannot
    give a ratio. A gene with no footprints AND no mRNA is missing rather than zero: it was not
    measured, not measured as untranslated. One stage only (trophozoite), so this is not a
    cycle-wide efficiency.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_translation_efficiency.py
"""
from _common import run

KEY = "pf_translation_efficiency"

if __name__ == "__main__":
    run(KEY)
