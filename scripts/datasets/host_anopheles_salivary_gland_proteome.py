#!/usr/bin/env python3
"""Anopheles salivary gland proteome, with and without a Plasmodium berghei infection

How abundant each mosquito protein is in the salivary gland (parts per million, uninfected), and how it changes 21 days after a Plasmodium berghei infection

    level / kind : translation / proteomics
    provides     : anopheles_salivary_gland_protein_ppm, anopheles_salivary_gland_infection_log2fc, anopheles_salivary_gland_infection_padj
    coverage     : 2,936 mosquito genes; 2,458 with an abundance, 2,366 with a contrast — 19.5% of the 12,574-row host space
    citation     : Alves e Silva TL et al., High-resolution proteomics unveils salivary gland disruption and saliva-hemolymph protein exchange in Plasmodium-infected mosquitoes, Nat Commun 2025 (PMID 41266344, doi 10.1038/s41467-025-64837-6), PRIDE PXD057628 (CC0)
    PMID         : 41266344
    accession    : PXD057628
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD057628
    local path   : datasets/translation/proteomics/PXD057628/

Quirks that cost time once:
    Glands dissected 21 days after the blood meal, when sporozoites sit in them: the stage the
    slot names. Mosquito: Anopheles gambiae, searched against the PEST gene set. The parasite is
    Plasmodium BERGHEI ANKA (GFPcon 259cl2), not falciparum, which the Plasmodium arm accepts by
    design; the same paper's saliva proteome also covers falciparum but saliva is not the gland.
    Three independent infections, each a pool of 25 gland pairs, uninfected (fed on uninfected
    blood) against infected; three replicates a side, so the test is weak: the smallest adjusted
    p is 0.093 and NO protein passes 0.05. Abundance is the uninfected arm only: each
    replicate's LFQ share of the summed LFQ of the kept host groups, in parts per million, then
    the mean. The contrast is Welch's t on log2 LFQ with Benjamini-Hochberg, needing two
    quantified values a side; the authors imputed missing values and this does not, so a protein
    seen in one arm only has no fold change. Groups that also hold a P. berghei protein (22) or
    span several genes (135) are dropped, and 15 AGAP ids the shipped gene build does not hold.
    Keyed on VectorBase AGAP gene ids. Redistributed with its accession and paper cited.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_anopheles_salivary_gland_proteome.py
"""
from _common import run

KEY = "host_anopheles_salivary_gland_proteome"

if __name__ == "__main__":
    run(KEY)
