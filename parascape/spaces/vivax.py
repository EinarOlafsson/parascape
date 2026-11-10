"""The Plasmodium vivax gene space (instruction 53, WP6): the other human malaria parasite.

Kept apart from falciparum rather than folded into it. P. vivax forms hypnozoites, a dormant liver stage that
falciparum does not have, and it is the species behind most relapsing malaria -- a slot about dormancy means
something different in each, so one table answering for both would be a false equivalence.

Declared from the Ensembl Protists build of PvP01. It holds the gene build only; the vivax measurements on
disk are a single liver-stage series that was read for falciparum, so the space can name genes and answer
nothing yet.
"""
from __future__ import annotations

from .. import organisms as O

VIVAX_SPACE = O.Space(O.VIVAX, 'Plasmodium vivax', 'Ensembl Protists PvP01 v2', O.PARASITE,
    # The organelle genes are named apart, as in falciparum: PVP01_API05000, PVP01_MIT02100. A pattern of
    # digits alone refuses 33 real genes, 30 apicoplast and 3 mitochondrial.
    r'PVP01_(?:\d{7}|API\d{5}|MIT\d{5})', ('PVP01_',), 'PlasmoDB',
    'https://plasmodb.org/plasmo/app/record/gene/{id}',
    'pv_gene_nodes.parquet', graph='pv_gene_graph.npz', orthomcl='pviv', distribution='pack')
