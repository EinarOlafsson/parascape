"""Gene-level protein abundance from the uniform Sage re-search of PRIDE projects (Einar 2026-11-26).

`scripts/ms_uniform_search.py` searches every Toxoplasma PRIDE project the same way (ME49 + human + mouse +
contaminant proteins, reversed decoys) and writes Sage's label-free quantification, `lfq.tsv`: one row per
peptide and charge, its proteins, a q-value, and an MS1 intensity per raw file. This module turns that into a
gene x raw-file table for Toxoplasma:

1. keep peptides at q <= 0.01;
2. keep only peptides UNIQUE to one Toxoplasma gene: every protein they match is a transcript of that gene (a
   peptide shared with a host or contaminant protein, or with a second parasite gene, says nothing about which);
3. gene intensity per file = Top3 (Silva et al. 2006, Mol Cell Proteomics 5:144): the mean of the three most
   intense unique peptides of the gene in that file (fewer when it has fewer); absent when none was measured;
4. log2, then each file is median-centred on the genes quantified in every file of the project.

Isobaric (iTRAQ / TMT) projects were searched label-free too, so their values are per raw file (the channels
summed), not per labelled sample; `labelled` in the project summary says so.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

Q_MAX = 0.01
TOP = 3
GENE = re.compile(r"^(TGME49_\d+)")
META = ("peptide", "charge", "proteins", "q_value", "score", "spectral_angle")


def gene_of(proteins: str) -> str | None:
    """The one Toxoplasma gene every protein of a peptide belongs to, else None."""
    genes = set()
    for p in str(proteins).split(";"):
        m = GENE.match(p.strip())
        if not m:
            return None                       # a host, contaminant or decoy protein: not unique to the parasite
        genes.add(m.group(1))
    return genes.pop() if len(genes) == 1 else None


def gene_table(lfq: pd.DataFrame, q_max: float = Q_MAX, top: int = TOP) -> tuple:
    """(log2 gene x file table, median-centred; per-gene unique peptide counts) from one Sage lfq table."""
    files = [c for c in lfq.columns if c not in META]
    t = lfq[lfq["q_value"] <= q_max].copy()
    t["gene"] = t["proteins"].map(gene_of)
    t = t[t["gene"].notna()]
    if t.empty or not files:
        return pd.DataFrame(columns=files), pd.Series(dtype=int)
    X = t[files].astype(float).where(lambda v: v > 0)
    X["gene"] = t["gene"].to_numpy()

    def top_mean(col):
        return X.groupby("gene")[col].apply(lambda v: v.dropna().nlargest(top).mean() if v.notna().any() else np.nan)

    G = pd.DataFrame({f: top_mean(f) for f in files})
    G = np.log2(G)
    shared = G.dropna()
    if len(shared) >= 10:
        G = G - shared.median(axis=0)
    else:
        G = G - G.median(axis=0)
    peptides = t.groupby("gene")["peptide"].nunique()
    return G, peptides


def summary(project: str, lfq: pd.DataFrame, G: pd.DataFrame, labelled: bool) -> dict:
    """What one project contributes: files, genes quantified, and how much of its signal is parasite."""
    passed = lfq[lfq["q_value"] <= Q_MAX]
    parasite = passed["proteins"].astype(str).str.contains("TGME49_").mean() if len(passed) else np.nan
    return {"project": project, "files": int(G.shape[1]), "peptides_q01": int(len(passed)),
            "parasite_share_of_peptides": float(parasite) if parasite == parasite else np.nan,
            "genes_quantified": int(G.notna().any(axis=1).sum()),
            "genes_in_every_file": int(G.notna().all(axis=1).sum()) if G.shape[1] else 0,
            "labelled": bool(labelled)}


#: Admission into the gene table (Einar 2026-11-26: per project now, designs later).
MIN_GENES = 100
MIN_RHO_PROTEOME = 0.15     # below this a project tracks enrichment (pull-down), not abundance
TABLE_DIR = "ms_uniform"
TABLE = "ms_project_abundance.parquet"
PROJECTS_TSV = TABLE.replace(".parquet", "_projects.tsv")


def project_abundance(G: pd.DataFrame) -> pd.Series:
    """One value per gene for a project: the mean of its median-centred log2 Top3 over the files that measured it."""
    return G.mean(axis=1, skipna=True)


def column_name(project: str, instrument: str, n_instruments: int) -> str:
    """Node column for one project (and instrument, when a project was searched in more than one)."""
    suffix = "" if n_instruments == 1 else "_" + instrument.replace("sage_", "").lower()
    return f"msu_{project}{suffix}"


def labels(base: str = "", log=print, resolve=None, root: str | None = None) -> pd.DataFrame:
    """Per-gene node columns (one per admitted project), indexed by gene id; joined by build_graph.load_nodes."""
    import os
    root = root or os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", TABLE_DIR)
    path = os.path.join(root, TABLE)
    if not os.path.exists(path):
        log(f"uniform MS re-search: no table at {root}; columns left out")
        return pd.DataFrame()
    t = pd.read_parquet(path).set_index("gene_id")
    if resolve is not None:
        t.index = pd.Index([resolve(g) or g for g in t.index], name="gene_id")
        t = t[~t.index.duplicated()]
    log(f"uniform MS re-search: {t.shape[1]} projects, {int(t.notna().any(axis=1).sum()):,} genes quantified")
    return t
