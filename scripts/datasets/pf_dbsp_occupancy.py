#!/usr/bin/env python3
"""PfDBSP occupancy per gene

How much PfDBSP sits on each gene and its promoter, as ChIP coverage over input

    level / kind : DNA / ChIPseq
    provides     : dbsp_occupancy
    coverage     : 5,318 genes
    citation     : GSE337518, PfDBSP ChIP-seq in P. falciparum 3D7
    accession    : GSE337518
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337518
    local path   : datasets/DNA/ChIPseq/GSE337518/GSM9856928_PfDBSP_inputnormalized.bedgraph.gz

Quirks that cost time once:
    An interval track, not a peak list, so the per-gene value is input-normalised coverage
    averaged over the gene and 2 kb upstream of its start on the strand it is read from; the
    gene coordinates come from Ensembl Protists r63, since neither shipped table carries them
    and PlasmoDB is behind a login. ONE REPLICATE: a single ChIP and a single input, so
    replicate agreement cannot be checked -- unlike the AP2-HS series refused the same day,
    whose replicates disagreed 25 peaks against 1,727, here there is nothing to disagree with. A
    continuous occupancy and never a bound/not-bound flag, because a threshold would invent a
    confidence one replicate does not support. Its target is described as a double-stranded-DNA
    binding protein rather than a classical ApiAP2; Einar took the call to use it on 2026-11-09.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_dbsp_occupancy.py
"""
from _common import run

KEY = "pf_dbsp_occupancy"

if __name__ == "__main__":
    run(KEY)
