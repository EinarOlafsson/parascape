#!/usr/bin/env python3
"""Uniform re-search of Toxoplasma PRIDE projects (protein abundance)

Per gene and admitted PRIDE project: mean log2 Top3 protein abundance (median-centred per raw file) from one uniform Sage search of every project's raw files

    level / kind : post_translation / proteomics
    provides     : msu_PXD001042, msu_PXD004083, msu_PXD004895, msu_PXD006235, msu_PXD010622_ftms, msu_PXD011692, msu_PXD012103, msu_PXD014133, msu_PXD015809, msu_PXD016383, msu_PXD016846_ftms, msu_PXD016846_itms, msu_PXD018056, msu_PXD018781, msu_PXD021516, msu_PXD022700, msu_PXD022713, msu_PXD022725, msu_PXD022785, msu_PXD023485, msu_PXD024210, msu_PXD024326_ftms, msu_PXD024326_itms, msu_PXD024656, msu_PXD028781, msu_PXD029455, msu_PXD030932, msu_PXD031116, msu_PXD031161, msu_PXD031164, msu_PXD031649, msu_PXD032102, msu_PXD032360, msu_PXD032362, msu_PXD033109, msu_PXD035106, msu_PXD035654, msu_PXD038158, msu_PXD039431, msu_PXD039432, msu_PXD039434, msu_PXD039896, msu_PXD039919, msu_PXD039922, msu_PXD039979, msu_PXD039983, msu_PXD039985, msu_PXD039986, msu_PXD041352, msu_PXD041887, msu_PXD042658, msu_PXD043808, msu_PXD043924, msu_PXD043994, msu_PXD044080, msu_PXD044081, msu_PXD044091, msu_PXD044588, msu_PXD044592, msu_PXD046083, msu_PXD046583, msu_PXD047027, msu_PXD047907, msu_PXD049694, msu_PXD049975, msu_PXD050011, msu_PXD051495, msu_PXD051551, msu_PXD053452, msu_PXD053703, msu_PXD053932, msu_PXD053944, msu_PXD054719, msu_PXD056911_ftms, msu_PXD057662, msu_PXD058095, msu_PXD059579, msu_PXD059833, msu_PXD060424, msu_PXD062403, msu_PXD063615, msu_PXD063619, msu_PXD063659, msu_PXD065585, msu_PXD066820, msu_PXD068228, msu_PXD068413, msu_PXD071065, msu_PXD072342, msu_PXD079431, msu_PXD079432, msu_PXD079900, msu_PXD083501
    coverage     : 93 PRIDE project folders, 5,422 genes
    accession    : PRIDE (ProteomeXchange); per-project list in parascape/data/ms_uniform/
    url          : https://ftp.pride.ebi.ac.uk/pride/data/archive/
    local path   : parascape/data/ms_uniform/ms_project_abundance.parquet

Quirks that cost time once:
    RE-SEARCHED, not deposited: raw files streamed from PRIDE, converted (ThermoRawFileParser),
    searched with Sage against ME49 + human + mouse + contaminants with reversed decoys
    (scripts/ms_uniform_search.py); peptides at q <= 0.01 unique to one Toxoplasma gene, Top3
    (Silva et al. 2006), log2, median-centred per file (parascape/ms_genes.py). Admitted: >= 100
    genes and Spearman >= 0.15 with the existing proteome; the reason for every folder left out
    is in ms_project_abundance_projects.tsv. Isobaric projects are per raw file (channels
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
