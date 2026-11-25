"""What single-gene papers established about each gene, read from their abstracts.

The free Toxoplasma Affinage run (Einar, 2026-11-24; after Di Bernardo and Cheeseman 2026,
arXiv 2607.02217). Each of 2,963 abstracts linked to a gene was read once for every gene it
names. Only findings that the paper shows by a direct experiment, for a protein that the
abstract itself names, were kept. Fixed rule checks with no model then screened the result.
The scripts, rules and checks are in `scripts/literature_affinage/`.

The four tables in `data/literature_affinage/`:
- `literature_findings.parquet`: one row per finding, with its paper, method and method tier;
- `literature_gene_labels.parquet`: observed labels per gene;
- `literature_partner_pairs.parquet`: gene pairs from partner and complex findings;
- `literature_doubts.parquet`: gene-matching errors noticed while reading.

These are OBSERVED labels from the source group `literature_affinage`. Literature co-mention
and attention depth come from the same papers, so they must not be used to score these labels
(the circularity rule of the Infer benchmark).
"""
from __future__ import annotations

import os

import pandas as pd

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "literature_affinage")
SOURCE_GROUP = "literature_affinage"

#: Gene-table columns this provides.
COLUMNS = ("lit_compartment", "lit_host_facing", "lit_n_papers", "lit_n_localisation", "lit_n_phenotype",
           "lit_n_activity", "lit_n_partner", "lit_n_complex", "lit_n_host_facing", "lit_n_structure")


def _read(name: str, root: str = HERE) -> pd.DataFrame:
    path = os.path.join(root, name)
    return pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()


def findings(root: str = HERE) -> pd.DataFrame:
    """Every admitted finding, one row each."""
    return _read("literature_findings.parquet", root)


def pairs(root: str = HERE) -> pd.DataFrame:
    """Gene pairs from established partner and complex findings, with their papers."""
    return _read("literature_partner_pairs.parquet", root)


def doubts(root: str = HERE) -> pd.DataFrame:
    """What the readers flagged: genes matched to the wrong paper, proteins not listed."""
    return _read("literature_doubts.parquet", root)


def labels(base: str = "", log=print, resolve=None, root: str = HERE) -> pd.DataFrame:
    """Per-gene columns for the gene table, indexed by gene id. `base` is unused, kept for symmetry."""
    t = _read("literature_gene_labels.parquet", root)
    if t.empty:
        log(f"literature findings: no table at {root}; columns left out")
        return pd.DataFrame()
    t = t.set_index("gene_id")
    if resolve is not None:
        t.index = pd.Index([resolve(g) or g for g in t.index], name="gene_id")
        t = t[~t.index.duplicated()]
    out = pd.DataFrame(index=t.index)
    out["lit_compartment"] = t["lit_compartment"].fillna("").astype(str)
    out["lit_host_facing"] = t["lit_host_facing"].astype(bool)
    for c in COLUMNS[2:]:
        out[c] = pd.to_numeric(t.get(c, 0), errors="coerce").fillna(0).astype(int)
    log(f"literature findings: {len(out):,} genes, {int((out.lit_compartment != '').sum()):,} with a "
        f"compartment shown by experiment, {int(out.lit_host_facing.sum())} host-facing")
    return out
