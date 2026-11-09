"""The cat gene space: Ensembl cat genes (instruction 53, WP4).

The cat is the definitive host: the only species where Toxoplasma completes its sexual cycle, so the
enteroepithelial stages exist nowhere else. Several slots of the catalogue name feline tissue for that reason
and have never had a species to hang on.

**It holds the gene build and nothing else yet.** No cat expression, localisation or infection data is on disk,
so the space is *declared and installable* rather than calibrated. It can name genes and hold measurements; it
can answer no question until measurements arrive.
"""
from __future__ import annotations

from .. import organisms as O

CAT_SPACE = O.Space(O.CAT, 'Felis catus', 'Ensembl Felis_catus_9.0', O.HOST,
    r'ENSFCAG\d{11}', ('ENSFCAG',), 'Ensembl',
    'https://www.ensembl.org/Felis_catus/Gene/Summary?g={id}',
    'fc_gene_nodes.parquet', graph='fc_gene_graph.npz', distribution='pack')
