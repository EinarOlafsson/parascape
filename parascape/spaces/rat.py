"""The rat gene space: Ensembl rat genes (instruction 53, WP4).

Declared like the mouse space (`spaces.mouse.MOUSE_SPACE`). The rat is a Toxoplasma host and the species the
behavioural work is done in, so the questions it will be asked are host questions.

**It holds the gene build and nothing else yet.** No rat expression, localisation or screen data is on disk, so
the space is *declared and installable* rather than calibrated: it can hold measurements and name genes, and it
can answer no question until measurements arrive. `space_catalogue` reports that state rather than hiding it.
"""
from __future__ import annotations

from .. import organisms as O

RAT_SPACE = O.Space(O.RAT, 'Rattus norvegicus', 'Ensembl mRatBN7.2', O.HOST,
    r'ENSRNOG\d{11}', ('ENSRNOG',), 'Ensembl',
    'https://www.ensembl.org/Rattus_norvegicus/Gene/Summary?g={id}',
    'rn_gene_nodes.parquet', graph='rn_gene_graph.npz', distribution='pack')
