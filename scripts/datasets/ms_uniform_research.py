#!/usr/bin/env python3
"""Uniform re-search of Toxoplasma PRIDE projects (protein abundance)

Per gene and admitted PRIDE project: mean log2 Top3 protein abundance (median-centred per raw file) from one uniform Sage search of every project's raw files

    level / kind : post_translation / proteomics
    provides     : msu_PXD001042, msu_PXD004895, msu_PXD010622_ftms, msu_PXD024210, msu_PXD044588, msu_PXD068413
    coverage     : PRIDE projects curated as abundance-meaningful (see ms_curated_projects.tsv)
    accession    : PRIDE (ProteomeXchange); per-project list in parascape/data/ms_uniform/
    url          : https://ftp.pride.ebi.ac.uk/pride/data/archive/
    local path   : parascape/data/ms_uniform/ms_project_abundance.parquet

Quirks that cost time once:
    RE-SEARCHED, not deposited: raw files streamed from PRIDE, converted (ThermoRawFileParser),
    searched with Sage against ME49 + human + mouse + contaminants with reversed decoys
    (scripts/ms_uniform_search.py); peptides at q <= 0.01 unique to one Toxoplasma gene, Top3
    (Silva et al. 2006), log2, median-centred per file (parascape/ms_genes.py). Admitted:
    curated as abundance-meaningful (whole proteome, subcellular fraction or secretome; IP,
    proximity labelling, crosslinking and PTM enrichment measure enrichment, not abundance), >=
    100 genes and Spearman >= 0.15 with the existing proteome; the reason for every folder left
    out is in ms_project_abundance_projects.tsv. Isobaric projects are per raw file (channels
    summed). Per-condition values follow when sample designs are curated (instruction 104).

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/ms_uniform_research.py
"""
from _common import run

KEY = "ms_uniform_research"

if __name__ == "__main__":
    run(KEY)
