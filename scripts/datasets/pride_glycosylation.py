#!/usr/bin/env python3
"""O-fucosylated glycoproteins (AAL pulldown)

Peptide identifications in the AAL lectin pulldown, per gene

    level / kind : post_translation / proteomics
    provides     : n_o_fucosyl_peptides
    coverage     : 394 genes
    citation     : Bandini G, Haserick JR, Motari E, Ouologuem DT, Lourido S, Roos DS, Costello CE, Robbins PW, Samuelson J. O-fucosylated glycoproteins form assemblies in close proximity to the nuclear pore complexes of Toxoplasma gondii. Proc Natl Acad Sci U S A. 2016 Oct 11;113(41):11567-11572
    accession    : PXD004426
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD004426
    local path   : datasets/quarantine/2026_08_16_pride/Tg/glycosylation/

Quirks that cost time once:
    Aleuria aurantia lectin affinity purification, so the count is peptide IDENTIFICATIONS and
    not sites -- the same standing as the proximity-labelling column, which is also a claim
    about what came down rather than about a residue. Verified against compartment: the pulldown
    is enriched for nucleus-chromatin at odds 3.71 (p = 1.5e-22) and cytosol at 2.61, and not
    enriched for mitochondrion. O-fucosylation through SPY is a nucleocytoplasmic modification
    and the paper describes punctiform signal beside the nuclei, so that is the right answer.
    Keyed on TGGT1_ accessions throughout.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pride_glycosylation.py
"""
from _common import run

KEY = "pride_glycosylation"

if __name__ == "__main__":
    run(KEY)
