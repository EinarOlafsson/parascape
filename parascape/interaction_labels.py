"""STRING and CORUM labels for Toxoplasma genes (instruction 103), shipped as package data.

Built by `scripts/labels/build_interaction_labels.py` and joined to the node table by
`build_graph.load_nodes` from release 0.0.14.0 (`datasets.JOINED_AT_NEXT_RELEASE`). Measured and copied
evidence are kept apart and the column names say which:

- `string_module_measured`: modules of STRING v12.0 T. gondii pairs with Toxoplasma's own experimental or
  curated-database evidence (observed);
- `string_module_orthology`: modules of STRING pairs whose evidence was transferred from orthologs (copied);
- `corum_complex_orthology`: CORUM human complexes mapped through OrthoMCL 6.21 orthogroups (copied).

How each was made, the classes and what is left out of the features are in `label_catalogue`.
"""
from __future__ import annotations

import os

import pandas as pd

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "interaction_labels")
FILE = "interaction_labels.parquet"
#: node column -> the catalogue label it carries
COLUMNS = {"string_module_measured": "interaction module (STRING, measured in Toxoplasma)",
           "string_module_orthology": "interaction module (STRING, copied by orthology)",
           "corum_complex_orthology": "complex (CORUM human, copied by orthology)"}


def table(root: str = HERE) -> pd.DataFrame:
    """The labels as built: gene_id plus one column per catalogue label (named as in the catalogue)."""
    path = os.path.join(root, FILE)
    return pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()


def labels(base: str = "", log=print, resolve=None, root: str = HERE) -> pd.DataFrame:
    """Per-gene node columns, indexed by gene id. `base` is unused, kept for symmetry with the other joins."""
    t = table(root)
    if t.empty:
        log(f"interaction labels: no table at {root}; columns left out")
        return pd.DataFrame()
    t = t.set_index("gene_id")
    if resolve is not None:
        t.index = pd.Index([resolve(g) or g for g in t.index], name="gene_id")
        t = t[~t.index.duplicated()]
    out = pd.DataFrame({col: t[label].astype(object).where(t[label].notna()) for col, label in COLUMNS.items()},
                       index=t.index)
    log("interaction labels: " + ", ".join(f"{c} {int(out[c].notna().sum()):,} genes" for c in out))
    return out
