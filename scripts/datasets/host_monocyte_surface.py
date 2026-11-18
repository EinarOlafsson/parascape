#!/usr/bin/env python3
"""Human monocyte cell-surface proteome (PXD013832)

Which proteins a primary human monocyte carries at the surface, and how strongly (TMT reporter intensity)

    level / kind : translation / proteomics
    provides     : monocyte_surface_classical_intensity
    coverage     : 432 proteins (436 of 437 quantified rows joined, one more outside the shipped host table)
    citation     : Ravenhill BJ et al., Comprehensive cell surface proteomics defines markers of classical, intermediate and non-classical monocytes, Sci Rep 2020 (PMID 32165698), PRIDE PXD013832
    PMID         : 32165698
    accession    : PXD013832
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD013832
    local path   : datasets/translation/proteomics/PXD013832/

Quirks that cost time once:
    Plasma membrane profiling (surface sialic acid oxidised and biotinylated, then streptavidin
    capture), TMT, FACS-sorted CLASSICAL monocytes from three blood donors; the value is the
    mean reporter intensity over the three donors. The spreadsheet is padded to 9,999 rows and
    437 of them carry a classical measurement. A protein not quantified is missing, which is
    UNKNOWN and not absent: the study reports no negatives. The list includes proteins the
    authors flag as nuclear or intracellular carry-over, so a value means recovered by the
    surface capture, not proven exposed. The intermediate and non-classical subsets belong to a
    second experiment and are not shipped. Keyed by reviewed UniProt accession with the isoform
    suffix removed (429 of 437), and by symbol for 7 more. Redistributed with its accession and
    paper cited.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_monocyte_surface.py
"""
from _common import run

KEY = "host_monocyte_surface"

if __name__ == "__main__":
    run(KEY)
