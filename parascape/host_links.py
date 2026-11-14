"""Carrying a selection between a host protein space and the host gene space of the same species.

Instruction 53, WP9 route 1, part 3. `hs_protein_gene_links.parquet` and
`mm_protein_gene_links.parquet` ship inside the host gene packs and hold one row per
(UniProt accession, Ensembl gene) pair, with a `link` column saying what the pair is. Measured on
the 2026-11-07 packs:

| | human | mouse |
|---|---:|---:|
| proteins | 20,989 | 15,590 |
| proteins reaching several genes (up to) | 1,572 (24) | 0 |
| genes reached by several proteins (up to) | 145 (5) | 0 |
| proteins reaching no reviewed gene | 1,703 | 1,151 |

So the link is many-to-many in human and one-to-one where it exists in mouse. A rule written against
the mouse table would look correct and be wrong in human, which is why the rule is stated here once
and both directions are tested.

## The rule

1. **A selection carries as the union of everything it names.** Every gene any selected protein
   names is carried, and every protein any selected gene names. Nothing is resolved to a single
   best partner: choosing one of 24 genes for a protein would put an invented fact on screen, and
   the reader would have no way to see it happened.
2. **An endpoint that names nothing carries nothing, and is counted.** It is never dropped in
   silence. `Carried.unresolved` holds those inputs, and `Carried.summary()` says how many.
3. **Nothing numeric travels.** Carrying moves WHICH rows are selected and no values at all. A
   measurement spread over a one-to-many link would be a measurement of the link.
4. **Carrying is not invertible, and says so.** Carry 1 protein out and back and you may get 3
   proteins, because the gene in the middle is named by 3. `Carried.grew` is true when the carried
   set is larger than what was put in, so a caller can tell the reader instead of implying that a
   round trip returns where it started.
5. **The two spaces stay separate windows.** Carrying opens or updates the other space's window; it
   never merges the two tables, because a protein is not a gene.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

#: code -> the link table inside that species' gene pack.
LINK_TABLES = {"Hs": "hs_protein_gene_links.parquet", "Mm": "mm_protein_gene_links.parquet"}

#: Host protein space code -> host gene space code, and back. The protein spaces are `HsP`/`MmP`
#: (see `spaces.host_proteins`), the gene spaces `Hs`/`Mm`.
GENE_OF_PROTEIN = {"HsP": "Hs", "MmP": "Mm"}
PROTEIN_OF_GENE = {gene: protein for protein, gene in GENE_OF_PROTEIN.items()}


@dataclass(frozen=True)
class Carried:
    """What a carried selection is, and what carrying it cost.

    `ids` is what to select on the other side. The rest is what a reader needs in order to trust it:
    which inputs named exactly one partner, which named several, which named none, and whether the
    selection grew on the way across.
    """

    direction: str                      # "proteins to genes" or "genes to proteins"
    code: str                           # the species code of the link table, Hs or Mm
    ids: list = field(default_factory=list)          # what is selected on the other side
    asked: list = field(default_factory=list)        # what was selected on this side
    unresolved: list = field(default_factory=list)   # inputs that name nothing over there
    one_to_one: int = 0                 # inputs naming exactly one partner
    fanned_out: int = 0                 # inputs naming several
    shared: int = 0                     # carried ids named by more than one input

    @property
    def grew(self) -> bool:
        """Whether the carried selection holds more rows than were put in."""
        return len(self.ids) > len(self.asked)

    def summary(self) -> str:
        """One plain sentence about what carried and what did not; "unknown" is never implied."""
        parts = [f"{len(self.asked):,} selected, {len(self.ids):,} carried {self.direction}"]
        if self.fanned_out:
            parts.append(f"{self.fanned_out:,} named several")
        if self.shared:
            parts.append(f"{self.shared:,} of the carried rows are named by more than one")
        if self.unresolved:
            parts.append(f"{len(self.unresolved):,} named none and carried nothing")
        return "; ".join(parts) + "."


def read_links(code: str) -> pd.DataFrame:
    """The (protein, gene) pairs of one species, out of its installed gene pack.

    Rows with no partner are kept, because "this protein names no reviewed gene" is a fact the
    carrying rule reports and not an absence of data.
    """
    from . import paths
    from .spaces import host_genes
    if code not in LINK_TABLES:
        raise KeyError(f"no host link table for {code!r}; there are {', '.join(LINK_TABLES)}")
    table = pd.read_parquet(paths.space_file(host_genes.space(code), LINK_TABLES[code]))
    for column in ("host_id", "gene_id"):
        if column not in table.columns:
            raise ValueError(f"{LINK_TABLES[code]}: expected host_id and gene_id columns")
    return table


def _carry(table: pd.DataFrame, source: str, target: str, asked, direction: str,
           code: str) -> Carried:
    """The rule, applied in one direction. Both directions share it, so neither can drift."""
    asked = [str(a) for a in dict.fromkeys(asked)]           # keep order, drop repeats
    if not asked:
        return Carried(direction, code)
    pairs = table[[source, target]].dropna()
    pairs = pairs[pairs[source].astype(str).isin(set(asked))]
    reached = pairs.groupby(pairs[source].astype(str))[target].apply(
        lambda values: list(dict.fromkeys(str(v) for v in values)))
    ids: list = []
    for one in asked:
        for partner in reached.get(one, ()):
            if partner not in ids:
                ids.append(partner)
    counts = pairs.groupby(pairs[target].astype(str))[source].nunique()
    return Carried(direction=direction, code=code, ids=ids, asked=asked,
                   unresolved=[a for a in asked if not len(reached.get(a, ()))],
                   one_to_one=sum(1 for a in asked if len(reached.get(a, ())) == 1),
                   fanned_out=sum(1 for a in asked if len(reached.get(a, ())) > 1),
                   shared=int((counts > 1).sum()))


def genes_of(code: str, proteins) -> Carried:
    """The genes a set of host proteins names, under the rule in the module docstring.

    `code` is the species code of the link table (Hs or Mm), not the protein space code; pass a
    protein space code through `GENE_OF_PROTEIN` first.
    """
    return _carry(read_links(code), "host_id", "gene_id", proteins, "proteins to genes", code)


def proteins_of(code: str, genes) -> Carried:
    """The host proteins a set of host genes names, under the same rule."""
    return _carry(read_links(code), "gene_id", "host_id", genes, "genes to proteins", code)
