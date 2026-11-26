#!/usr/bin/env python3
"""CORUM human complexes via OrthoMCL 6.21 orthogroups

The most specific human CORUM complex that keeps 10 Toxoplasma genes, mapped through orthogroups

    level / kind : reference / orthology
    provides     : corum_complex_orthology
    coverage     : 414 genes, 26 complexes
    citation     : Tsitsiridis G et al. Nucleic Acids Res 2023;51:D539-D545; OrthoMCL release 6.21
    PMID         : 36382402
    url          : https://mips.helmholtz-muenchen.de/corum/
    local path   : parascape/data/interaction_labels/interaction_labels.parquet

Quirks that cost time once:
    COPIED by orthology: subunit UniProt -> Ensembl gene (HGNC) -> OrthoMCL 6.21 orthogroup ->
    Toxoplasma genes; orthogroups with more than 3 Toxoplasma genes are not used and a complex
    needs 2 mapped orthogroups.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`parascape.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/corum_human_by_orthology.py
"""
from _common import run

KEY = "corum_human_by_orthology"

if __name__ == "__main__":
    run(KEY)
