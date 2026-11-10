"""The Anopheles stephensi gene space (instruction 53, WP5): the second malaria vector.

A separate space from `anopheles.ANOPHELES_SPACE` rather than a strain of it. A. stephensi is a different
species with its own genome: it carries malaria across South Asia and has spread into Africa, where it breeds
in urban containers rather than in rural pools as A. gambiae does. A question about where a vector lives, or
which insecticide reaches it, has a different answer in each, so one table for both would be a false
equivalence -- the same reason P. vivax is kept apart from P. falciparum.

Declared from the Ensembl Metazoa build. Gene build only: no stephensi measurement is on disk.
"""
from __future__ import annotations

from .. import organisms as O

STEPHENSI_SPACE = O.Space(O.STEPHENSI, 'Anopheles stephensi', 'Ensembl Metazoa AsteS1', O.VECTOR,
    r'ASTE\d{6}', ('ASTE',), 'VectorBase',
    'https://vectorbase.org/vectorbase/app/record/gene/{id}',
    'as_gene_nodes.parquet', distribution='pack')
# No graph: no edges have been computed for this space. The attribute is left unset rather than naming a
# file that does not exist -- `build_pack` validates whatever a space declares, so a claimed graph is a
# promise the pack cannot keep. It is set when the first edge layer is built. (no graph yet)
