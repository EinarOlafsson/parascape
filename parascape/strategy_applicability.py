"""Which strategies can run in which space, decided from that space own data (instruction 53, WP8).

A strategy declares what it `needs` in words -- "a categorical column", "two or more edge layers",
"the other organism's table". Those words are the vocabulary this module checks against a space, so
the answer for a new space is read off its data rather than assumed from Toxoplasma. A space with no
graph cannot run a network strategy; a space with no partner table cannot transfer from one; a space
with one categorical column cannot run a strategy that compares two.

Nothing here runs a strategy. The question is only whether it COULD, and when it could not, which
requirement is missing -- so a new space arrives with an honest list rather than a wall of failures.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

#: Column names that are identity or bookkeeping rather than measurements, so a strategy asking for
#: "a categorical column" is not satisfied by the gene id.
NOT_EVIDENCE = ("gene_id", "product", "description", "name")

#: Layer-name fragments that mark a measured interaction layer. Spelled out because the graphs name
#: the same evidence differently -- `xlms` and `ip_ms` in the shipped tables, `crosslink` elsewhere --
#: and a check that missed one would report Plasmodium as having no interaction evidence when it has
#: both crosslinking and IP-MS.
INTERACTION_LAYERS = ("interaction", "crosslink", "xlms", "ip_ms", "ipms")


@dataclass(frozen=True)
class SpaceEvidence:
    """What one space has, in the terms strategies ask about."""
    code: str
    genes: int
    categorical: int
    numeric: int
    edge_layers: tuple
    has_partner_table: bool
    has_orthogroups: bool
    has_literature_attention: bool

    def satisfies(self, need: str) -> tuple:
        """(can it, why not) for one declared requirement, in the strategy own words."""
        layers = len(self.edge_layers)
        text = need.lower()
        if "two categorical" in text:
            return (self.categorical >= 2, f"{self.categorical} categorical columns")
        if "categorical" in text:
            return (self.categorical >= 1, "no categorical column")
        if "two numeric" in text:
            return (self.numeric >= 2, f"{self.numeric} numeric columns")
        if "numeric" in text:
            return (self.numeric >= 1, "no numeric column")
        if "two or more edge layers" in text:
            return (layers >= 2, f"{layers} edge layers")
        if "measured edge layers" in text or "edge layer" in text:
            return (layers >= 1, "no edge layer")
        if "crosslink or ip-ms" in text:
            return (any(k in l for l in self.edge_layers for k in INTERACTION_LAYERS),
                    "no crosslink or IP-MS layer")
        if "structural-similarity" in text:
            return (any("struct" in l for l in self.edge_layers), "no structural-similarity layer")
        if "interaction layer" in text:
            return (any(k in l for l in self.edge_layers for k in INTERACTION_LAYERS),
                    "no interaction layer")
        if "co-mention" in text:
            return (any("mention" in l or "pubmed" in l for l in self.edge_layers), "no co-mention layer")
        if "measurement and literature" in text:
            return (self.numeric >= 1 and self.has_literature_attention,
                    "no literature-attention column" if self.numeric else "no measurements")
        if "literature-attention" in text:
            return (self.has_literature_attention, "no literature-attention column")
        if "orthogroup" in text:
            return (self.has_orthogroups, "no orthogroup column")
        if "other organism" in text:
            return (self.has_partner_table, "no partner table declared for this space")
        if "gene list" in text:
            return (self.genes > 0, "no genes")        # any space with genes can be given a list
        return (True, "")                              # a requirement this module does not know of


def evidence_of(code: str, nodes: pd.DataFrame, edge_layers=(), partner: bool | None = None) -> SpaceEvidence:
    """Read one space table in the terms strategies ask about.

    `edge_layers` are the layer names of its graph, empty when it has none. `partner` overrides the
    declared partner of the space, for a space whose counterpart table is not built yet.
    """
    from . import organisms

    columns = [c for c in nodes.columns if c not in NOT_EVIDENCE]
    numeric = [c for c in columns if pd.api.types.is_numeric_dtype(nodes[c])
               and not pd.api.types.is_bool_dtype(nodes[c])]
    categorical = [c for c in columns if c not in numeric and nodes[c].nunique(dropna=True) > 1]
    if partner is None:
        try:
            space = organisms.get(code)
            partner = bool(space.partner) and space.partner in organisms.SPACES
        except KeyError:
            partner = False
    attention = any("pubmed" in c or "attention" in c or "citation" in c for c in nodes.columns)
    return SpaceEvidence(code=code, genes=len(nodes), categorical=len(categorical), numeric=len(numeric),
                         edge_layers=tuple(edge_layers), has_partner_table=bool(partner),
                         has_orthogroups="orthogroup" in nodes.columns, has_literature_attention=attention)


def applicability(evidence: SpaceEvidence, registry=None) -> pd.DataFrame:
    """One row per strategy: can it run in this space, and if not, what is missing.

    `registry` defaults to every strategy registered when this is called, so a caller that has not
    imported the strategy modules gets an empty table rather than a wrong one.
    """
    from . import strategies

    registry = strategies.REGISTRY if registry is None else registry
    rows = []
    for strategy in registry.values():
        missing = []
        for need in strategy.needs:
            ok, why = evidence.satisfies(need)
            if not ok:
                missing.append(why or need)
        rows.append(dict(space=evidence.code, number=strategy.number, strategy=strategy.key,
                         family=strategy.family, task=getattr(strategy, "task", ""),
                         needs="; ".join(strategy.needs) or "(nothing in particular)",
                         applicable=not missing, missing="; ".join(missing)))
    frame = pd.DataFrame(rows)
    return frame.sort_values(["applicable", "number"], ascending=[False, True]).reset_index(drop=True)
