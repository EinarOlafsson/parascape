"""The Cryptosporidium parvum gene space (instruction 53, WP6): Iowa II genes, CryptoDB identifiers.

Declared like the host spaces and built locally (`scripts/build_cryptosporidium_space_notebook.py`); a distributable
pack waits on verifying the redistribution terms of its sources (CryptoDB tables, the Guerin 2023 LOPIT supplement).
"""
from __future__ import annotations

from .. import organisms as O

CRYPTO_SPACE = O.Space(O.CRYPTOSPORIDIUM, 'Cryptosporidium parvum', 'Iowa II', O.PARASITE,
    r'cgd\d+_\d+', ('cgd',), 'CryptoDB', 'https://cryptodb.org/cryptodb/app/record/gene/{id}',
    'cp_gene_nodes.parquet', graph='cp_gene_graph.npz', orthomcl='cpar',
    targets=('lopit_location',), distribution='pack')
