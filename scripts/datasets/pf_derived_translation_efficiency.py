#!/usr/bin/env python3
"""Translation efficiency per blood-stage phase (DERIVED)

Ribosome footprints over mRNA in the same phase, and how much that changes across the cycle

    level / kind : translation / riboseq
    provides     : te_ring_derived, te_early_trophozoite_derived, te_late_trophozoite_derived, te_schizont_derived, te_merozoite_derived, te_phase_span_derived
    coverage     : 1,174 to 2,983 genes per phase; 3,020 with a span
    local path   : scripts/add_pf_translation_efficiency_per_phase.py
    derived from : riboseq_rpf_ring, riboseq_rpf_early_trophozoite, riboseq_rpf_late_trophozoite, riboseq_rpf_schizont, riboseq_rpf_merozoite, riboseq_mrna_ring, riboseq_mrna_early_trophozoite, riboseq_mrna_late_trophozoite, riboseq_mrna_schizont, riboseq_mrna_merozoite

Quirks that cost time once:
    Computed, not acquired, and the reason is worth recording: the literature index proposed
    GSE58402 for the translation-efficiency slot, and GSE58402 was ALREADY shipped as
    `pf_riboseq` -- same accession, same PMID 25493618. The slot read empty only because no
    column divided footprints by mRNA, which is what efficiency is. log2(RPF/mRNA) within the
    SAME phase, so a phase a gene was not measured in stays missing rather than borrowing
    another phase's mRNA; the sources are already length-normalised densities, so the ratio
    needs no further correction. The span needs two phases, because one phase is a value and not
    a change.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_derived_translation_efficiency.py
"""
from _common import run

KEY = "pf_derived_translation_efficiency"

if __name__ == "__main__":
    run(KEY)
