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
            rows.append(dict(source_organism=organisms.TOXOPLASMA if "Toxoplasma" in str(r.organism) else str(r.organism),
                             source_id=str(r.gene_id or r.parasite_protein), target_organism=HOST_UNRECORDED,
                             target_id=str(r.host_target), relation="measured interaction (curated literature)",
                             context=str(r.interaction_mode), support=float(r.confidence)
                             if r.confidence == r.confidence else float("nan"), reference=str(r.key_reference)))
    path = _data("pf_host_bridges.parquet")
    if os.path.exists(path):
        p = pd.read_parquet(path)
        for r in p.itertuples(index=False):
            rows.append(dict(source_organism=organisms.FALCIPARUM, source_id=str(r.gene_id), target_organism="human",
                             target_id=f"{r.host_name} ({r.host_id})", relation="measured interaction (crosslinking MS)",
                             context=f"{int(r.crosslinks)} crosslink(s)", support=float(r.crosslinks),
                             reference=str(r.evidence)))
    return pd.DataFrame(rows, columns=COLUMNS)


def orthology_links() -> pd.DataFrame:
    """Toxoplasma <-> Plasmodium genes sharing an OrthoMCL group in their own node tables."""
    try:
        tg = pd.read_parquet(organisms.nodes_path(organisms.TOXOPLASMA), columns=["gene_id", "orthogroup"]).dropna()
        pf = pd.read_parquet(organisms.nodes_path(organisms.FALCIPARUM), columns=["gene_id", "orthogroup"]).dropna()
    except Exception:
        return pd.DataFrame(columns=COLUMNS)
    m = tg.merge(pf, on="orthogroup", suffixes=("_tg", "_pf"))
    return pd.DataFrame({"source_organism": organisms.TOXOPLASMA, "source_id": m["gene_id_tg"].astype(str),
                         "target_organism": organisms.FALCIPARUM, "target_id": m["gene_id_pf"].astype(str),
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


# --------------------------------------------------------------------------- the other side, in its own space
HOST_SPACES = {"human": "Hs", "mouse": "Mm"}
_HOST_FILES = {"Hs": ("hs_gene_nodes.parquet", "hs_protein_gene_links.parquet"),
               "Mm": ("mm_gene_nodes.parquet", "mm_protein_gene_links.parquet")}


def _accession(target_id: str) -> str:
    """The protein accession of a host endpoint written "NAME (ACCESSION)", or the id itself."""
    text = str(target_id)
    return text[text.rfind("(") + 1:-1] if text.endswith(")") and "(" in text else text


def counterpart(link, from_gene: str = "") -> dict:
    """Where the other end of a link lives: its own space, its genes there and that space's own rows.

    Returns ``{"organism", "id", "space", "genes", "rows", "space_genes", "status"}``. A parasite end opens in its
    own map; a host protein end is resolved to host GENES through the installed host pack's reviewed protein links
    (instructions 64.21/22) -- never by name -- and its rows come from that host table, so every count stays
    species-specific (``space_genes`` is the host space's own denominator). A host whose species the source does
    not record gets no space. Nothing here claims an inference: it is a lookup.
    """
    org, other = str(link["target_organism"]), str(link["target_id"])
    if from_gene and other == str(from_gene):
        org, other = str(link["source_organism"]), str(link["source_id"])
    out = dict(organism=org, id=other, space=None, genes=[], rows=pd.DataFrame(), space_genes=0, status="")
    if org == HOST_UNRECORDED:
        out["status"] = "the source does not say which host species; no space is chosen"
        return out
    if org in organisms.SPACES and organisms.get(org).distribution == "wheel":
        try:
            nodes = pd.read_parquet(organisms.nodes_path(org))
        except Exception:
            out["status"] = f"the {organisms.get(org).species} map is not built here"
            return out
        out.update(space=org, genes=[other], rows=nodes[nodes["gene_id"].astype(str) == other],
                   space_genes=len(nodes), status=f"opens in the {organisms.get(org).species} map")
        return out
    code = HOST_SPACES.get(org)
    if code is None:
        out["status"] = f"no gene space is declared for {org}"
        return out
    from . import packs
    try:
        where = packs.installed(code)
    except Exception as error:
        where, out["status"] = None, f"the {org} gene space could not be read: {error}"
    if where is None:
        out["status"] = out["status"] or f"the {org} gene space is not installed (data pack {code})"
        return out
    nodes_file, links_file = _HOST_FILES[code]
    nodes = pd.read_parquet(where / nodes_file)
    links = pd.read_parquet(where / links_file)
    genes = sorted(links.loc[links["host_id"].astype(str) == _accession(other), "gene_id"].dropna().astype(str).unique())
    out.update(space=code, genes=genes, rows=nodes[nodes["gene_id"].isin(genes)], space_genes=len(nodes))
    out["status"] = (f"{len(genes)} gene{'s' if len(genes) != 1 else ''} in the {org} gene space of {len(nodes):,}"
                     if genes else f"this protein has no reviewed gene in the {org} gene space")
    return out
