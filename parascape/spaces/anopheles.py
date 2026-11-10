"""The Anopheles gambiae gene space (instruction 53, WP5): the mosquito that carries Plasmodium.

A vector, not a host: the parasite completes its sexual cycle here, so the ookinete, oocyst and sporozoite
stages happen inside this animal. `organisms.VECTOR` is its kind, which is why the space exists separately
from the host spaces rather than beside them.

Declared from the Ensembl Metazoa build. It holds the gene build only; no mosquito measurement is on disk, so
the space can name genes and answer nothing yet.
"""
from __future__ import annotations

from .. import organisms as O

ANOPHELES_SPACE = O.Space(O.ANOPHELES, 'Anopheles gambiae', 'Ensembl Metazoa AgamP4', O.VECTOR,
    r'AGAP\d{6}', ('AGAP',), 'VectorBase',
    'https://vectorbase.org/vectorbase/app/record/gene/{id}',
    'ag_gene_nodes.parquet', graph='ag_gene_graph.npz', distribution='pack')
