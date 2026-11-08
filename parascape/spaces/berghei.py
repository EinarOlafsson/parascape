"""The Plasmodium berghei ANKA gene space (instruction 53, WP6): PlasmoDB PBANKA identifiers.

PlasmoGEM blood-stage knockout phenotypes (Bushell et al. 2017) get their native home here instead of existing
only as transfers into P. falciparum (`pb_transferred_*`). Built locally by
`scripts/build_berghei_space_notebook.py`; the PlasmoDB annotation terms are verified before a pack is distributed.
"""
from __future__ import annotations

from .. import organisms as O

BERGHEI_SPACE = O.Space(O.BERGHEI, 'Plasmodium berghei', 'ANKA', O.PARASITE,
    r'PBANKA_\d{7}', ('PBANKA_',), 'PlasmoDB', 'https://plasmodb.org/plasmo/app/record/gene/{id}',
    'pb_gene_nodes.parquet', graph='pb_gene_graph.npz', orthomcl='pber',
    targets=('plasmogem_phenotype',), distribution='pack')
