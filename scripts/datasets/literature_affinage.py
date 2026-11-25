#!/usr/bin/env python3
"""Experimental findings read from single-gene papers

Per gene: the compartment shown by experiment, whether it is shown to act on the host, and how many papers report location, phenotype, activity, partners, complexes and structure

    level / kind : reference / literature
    provides     : lit_compartment, lit_host_facing, lit_n_papers, lit_n_localisation, lit_n_phenotype, lit_n_activity, lit_n_partner, lit_n_complex, lit_n_host_facing, lit_n_structure
    coverage     : 747 genes, 977 papers
    accession    : PubMed abstracts via NCBI E-utilities
    url          : https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
    local path   : parascape/data/literature_affinage/literature_gene_labels.parquet

Quirks that cost time once:
    READ, not deposited: 2,963 abstracts linked to genes by the release literature scan
    (catalogue papers naming more than 50 genes dropped, at most 20 per gene) were each read
    once by Claude Code subagents under fixed rules
    (scripts/literature_affinage/READING_RULES.md): a finding is kept only when THIS paper shows
    it by a direct experiment for a protein the abstract itself names. Fixed checks A1-A7 with
    no model screened the output; seven ATP synthase subunits credited from their product names
    were excluded by hand, and 46 findings whose protein the abstract names only as a family are
    kept in the findings table but not used here. Where hyperLOPIT also places the gene, 67.5%
    of 514 agree (a secreted destination counts as its organelle of origin). These labels come
    from the same papers as literature co-mention and attention depth, so those must never be
    used to score them. Method after Di Bernardo and Cheeseman, arXiv 2607.02217 (2026).

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/literature_affinage.py
"""
from _common import run

KEY = "literature_affinage"

if __name__ == "__main__":
    run(KEY)
