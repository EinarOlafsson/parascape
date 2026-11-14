"""The mouse gene space: Ensembl mouse genes, kept apart from the mouse protein table.

Declared like the human space (`spaces.human.HUMAN_SPACE`) and distributed as a pack. Mouse subcellular
locations in the source are transferred from human orthologs (Human Protein Atlas); they are inferred, never
observed, and the column names say so.
"""
from __future__ import annotations

from .. import organisms as O

# The DEFAULT layout is the tissue lens (instruction 53, route 1). The 2026.11.07 pack
# shipped `mm_gene_graph.npz`, whose recipe was not recorded beside it; each lens layout
# records its features, its recipe and the SHA-256 of the table it was built from.
MOUSE_SPACE = O.Space(O.MOUSE, 'Mus musculus', 'Ensembl GRCm39 / Bgee', O.HOST,
    r'ENSMUSG\d{11}', ('ENSMUSG',), 'Ensembl', 'https://www.ensembl.org/Mus_musculus/Gene/Summary?g={id}',
    'mm_gene_nodes.parquet', graph='mm_gene_graph_tissue.npz',
    numbers=('bgee_top_score', 'bgee_tissues_present'), distribution='pack')
