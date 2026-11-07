"""Cross-space evidence links, organism-qualified (instruction 64.24): only bridges already verified.

One table for every link that leaves an organism's own space: who is linked to whom, in which organism, by what
kind of relation, from which source, with what support. The relation types are kept apart and never merged into
one score -- a measured interaction, an orthology and an inferred link are different claims:

* ``measured interaction (curated literature)`` -- `host_interactions.parquet`, parasite protein -> host target,
  with interaction mode, confidence and reference;
* ``measured interaction (crosslinking MS)`` -- `pf_host_bridges.parquet`, Plasmodium protein -> host protein;
* ``orthology (shared OrthoMCL group)`` -- Toxoplasma <-> Plasmodium genes in one orthogroup of their node tables.

The host species is written where the source records it and "host (species not recorded)" where it does not;
nothing is guessed. Exploring a link claims no inference accuracy; inference through a bridge needs its own
registered benchmark.
"""
from __future__ import annotations

import os

import pandas as pd

from . import organisms

COLUMNS = ["source_organism", "source_id", "target_organism", "target_id", "relation", "context", "support",
           "reference"]
HOST_UNRECORDED = "host (species not recorded)"


def _data(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), "data", name)


def host_links() -> pd.DataFrame:
    """Parasite -> host links from the two measured sources."""
    rows = []
    path = _data("host_interactions.parquet")
    if os.path.exists(path):
        h = pd.read_parquet(path)
        for r in h.itertuples(index=False):
            rows.append(dict(source_organism="Tg" if "Toxoplasma" in str(r.organism) else str(r.organism),
                             source_id=str(r.gene_id or r.parasite_protein), target_organism=HOST_UNRECORDED,
                             target_id=str(r.host_target), relation="measured interaction (curated literature)",
                             context=str(r.interaction_mode), support=float(r.confidence)
                             if r.confidence == r.confidence else float("nan"), reference=str(r.key_reference)))
    path = _data("pf_host_bridges.parquet")
    if os.path.exists(path):
        p = pd.read_parquet(path)
        for r in p.itertuples(index=False):
            rows.append(dict(source_organism="Pf", source_id=str(r.gene_id), target_organism="human",
                             target_id=f"{r.host_name} ({r.host_id})", relation="measured interaction (crosslinking MS)",
                             context=f"{int(r.crosslinks)} crosslink(s)", support=float(r.crosslinks),
                             reference=str(r.evidence)))
    return pd.DataFrame(rows, columns=COLUMNS)


def orthology_links() -> pd.DataFrame:
    """Toxoplasma <-> Plasmodium genes sharing an OrthoMCL group in their own node tables."""
    try:
        tg = pd.read_parquet(organisms.nodes_path("Tg"), columns=["gene_id", "orthogroup"]).dropna()
        pf = pd.read_parquet(organisms.nodes_path("Pf"), columns=["gene_id", "orthogroup"]).dropna()
    except Exception:
        return pd.DataFrame(columns=COLUMNS)
    m = tg.merge(pf, on="orthogroup", suffixes=("_tg", "_pf"))
    return pd.DataFrame({"source_organism": "Tg", "source_id": m["gene_id_tg"].astype(str),
                         "target_organism": "Pf", "target_id": m["gene_id_pf"].astype(str),
                         "relation": "orthology (shared OrthoMCL group)", "context": m["orthogroup"].astype(str),
                         "support": float("nan"), "reference": "node tables (orthogroup column)"}, columns=COLUMNS)


def all_links() -> pd.DataFrame:
    """Every verified cross-space link, relation types kept apart."""
    return pd.concat([host_links(), orthology_links()], ignore_index=True)


def coverage(links: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per relation: links, distinct source genes and distinct targets -- the mapping/coverage evidence."""
    links = all_links() if links is None else links
    return links.groupby("relation").agg(links=("source_id", "size"), sources=("source_id", "nunique"),
                                         targets=("target_id", "nunique")).reset_index()


def links_of(gene: str, links: pd.DataFrame | None = None) -> pd.DataFrame:
    """The cross-space links of one gene, from either side."""
    links = all_links() if links is None else links
    return links[(links["source_id"] == str(gene)) | (links["target_id"] == str(gene))]
