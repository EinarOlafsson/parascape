#!/usr/bin/env python3
"""Human brain-cell transcription by cell type

How much each gene is transcribed in a neuron, an astrocyte or a microglial cell

    level / kind : transcription / scRNAseq
    provides     : neuron_excitatory_ncpm, neuron_inhibitory_ncpm, astrocyte_ncpm, microglia_ncpm
    coverage     : 18,939 reviewed proteins
    citation     : Human Protein Atlas single-cell RNA (rna_single_cell_type.tsv)
    url          : https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip
    local path   : datasets/host/single_cell_rna/rna_single_cell_type.tsv.zip

Quirks that cost time once:
    The bradyzoite cyst sits inside a NEURON, and the cells around it are astrocytes and
    microglia. Bulk brain tissue cannot answer that: a tissue average mixes the cell the
    parasite lives in with the ones reacting to it, and the question is the difference.
    Excitatory and inhibitory neurons are kept apart because the Atlas measures them separately;
    collapsing them would be our averaging rather than theirs. nCPM is the Atlas's cross-cell-
    type normalisation, so two columns can be compared.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_brain_cell_rna.py
"""
from _common import run

KEY = "host_brain_cell_rna"

if __name__ == "__main__":
    run(KEY)
