"""Which known perturbation does my mutant look like? (instruction 71, product 2)

The signature store is the within-study contrast columns of the expression traits release: a shrunken
log2 fold change per gene, with its standard error, each computed inside one study against that study's
own control -- never raw values pooled across studies. A query signature is ranked against every stored
signature by rank correlation over the genes both measure; a study can be held out so a perturbation is
only ever matched to OTHER labs' experiments.

    store = load()
    hits = query(my_signature, store, exclude_study="GSE...")

No window here; every rule is testable.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

#: A comparison over fewer shared genes than this is not reported.
MIN_SHARED_GENES = 1000


def default_release() -> str:
    """The newest expression traits release in the repository results."""
    root = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    found = sorted(d for d in os.listdir(root) if d.startswith("expression_traits_20")
                   and os.path.exists(os.path.join(root, d, "expression_traits.parquet")))
    return os.path.join(root, found[-1])


def load(release: str | None = None) -> dict:
    """The signature store: fold changes (genes x signatures), their registry rows, the release."""
    release = release or default_release()
    values = pd.read_parquet(os.path.join(release, "expression_traits.parquet"))
    registry = pd.read_csv(os.path.join(release, "column_registry.tsv"), sep="\t")
    sig = registry[(registry["quantity"] == "log2_fold_change") & registry["column"].isin(values.columns)]
    lfc = values[sig["column"].tolist()]
    keep = lfc.notna().sum() >= MIN_SHARED_GENES
    return {"lfc": lfc.loc[:, keep], "registry": sig[keep.values].set_index("column"),
            "release": os.path.basename(release)}


def _ranks(x: np.ndarray) -> np.ndarray:
    out = np.full(x.shape, np.nan)
    ok = ~np.isnan(x)
    out[ok] = rankdata(x[ok])
    return out


def query(signature: pd.Series, store: dict, exclude_study: str | None = None,
          exclude_columns=()) -> pd.DataFrame:
    """Every stored signature ranked by Spearman correlation with `signature` (indexed by gene id)."""
    lfc, reg = store["lfc"], store["registry"]
    q = pd.Series(signature, dtype=float).reindex(lfc.index)
    rows = []
    for column in lfc.columns:
        study = str(reg.loc[column, "study"])
        if column in exclude_columns or (exclude_study and study == exclude_study):
            continue
        both = q.notna() & lfc[column].notna()
        n = int(both.sum())
        if n < MIN_SHARED_GENES:
            continue
        a, b = _ranks(q[both].to_numpy()), _ranks(lfc.loc[both, column].to_numpy())
        rho = float(np.corrcoef(a, b)[0, 1])
        rows.append((column, study, str(reg.loc[column, "contrast_label"]), rho, n))
    out = pd.DataFrame(rows, columns=["signature", "study", "contrast", "spearman", "genes"])
    out = out.sort_values("spearman", ascending=False, ignore_index=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def recovery(store: dict, groups: dict) -> pd.DataFrame:
    """For each member of each replicate group: the rank of its best partner from ANOTHER study.

    `groups` maps a perturbation name to signature columns that should find each other across
    independent studies. The query's own study is held out entirely, so a perturbation can only be
    matched to another lab's experiment; chance is the expected best rank of the partners among the
    candidates, (candidates + 1) / (partners + 1).
    """
    reg = store["registry"]
    rows = []
    for name, members in groups.items():
        members = [m for m in members if m in reg.index]
        for column in members:
            study = str(reg.loc[column, "study"])
            partners = [m for m in members if str(reg.loc[m, "study"]) != study]
            if not partners:
                continue
            hits = query(store["lfc"][column], store, exclude_study=study)
            ranks = hits.loc[hits["signature"].isin(partners), "rank"]
            if not len(ranks):
                continue
            best = int(ranks.min())
            rows.append(dict(perturbation=name, query=column, study=study, partners=len(partners),
                             candidates=len(hits), best_partner_rank=best,
                             chance_rank=round((len(hits) + 1) / (len(partners) + 1), 1),
                             top1=best == 1, top5=best <= 5,
                             best_partner=hits.loc[hits["rank"] == best, "signature"].iloc[0]))
    return pd.DataFrame(rows)


def host_store(atlas: str | None = None) -> dict:
    """The host response atlas (instruction 71, product 3) as a signature store for `query` and `recovery`.

    Same shape as `load`: host fold changes (host genes x signatures) and their registry, with a readable contrast
    label. Human and mouse signatures share no gene identifiers, so a query only ever meets its own species.
    """
    root = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    if atlas is None:
        found = sorted(d for d in os.listdir(root) if d.startswith("host_response_atlas_20")
                       and os.path.exists(os.path.join(root, d, "host_signatures.parquet")))
        atlas = os.path.join(root, found[-1])
    lfc = pd.read_parquet(os.path.join(atlas, "host_signatures.parquet"))
    reg = pd.read_csv(os.path.join(atlas, "registry.tsv"), sep="\t").set_index("column")
    reg["contrast_label"] = (reg["axis"].astype(str) + ": " + reg["treated"].astype(str) + " vs "
                             + reg["baseline"].astype(str) + " (" + reg["background"].astype(str) + ")")
    keep = [c for c in lfc.columns if c in reg.index and lfc[c].notna().sum() >= MIN_SHARED_GENES]
    return {"lfc": lfc[keep], "registry": reg.loc[keep], "release": os.path.basename(atlas)}


def replicate_groups(store: dict, min_studies: int = 2) -> dict:
    """Signatures naming the same treated condition, in the same host species, in at least `min_studies` studies.

    A data-driven replication set for `recovery`: the treated label is normalised (case, spaces, punctuation);
    a group names a perturbation only when independent studies state it.
    """
    reg = store["registry"]
    key = (reg["species"].astype(str) + "|" + reg["axis"].astype(str) + "|"
           + reg["treated"].astype(str).str.lower().str.replace(r"[^a-z0-9]+", "", regex=True))
    groups = {}
    for k, cols in reg.groupby(key).groups.items():
        if k.rsplit("|", 1)[-1] in ("unstated", "none", "nan", "", "control", "wt", "wildtype"):
            continue                               # a missing or baseline label names no perturbation
        if reg.loc[list(cols), "study"].nunique() >= min_studies:
            groups[k] = list(cols)
    return groups


def host_symbol_map(root: str | None = None) -> pd.DataFrame:
    """Human and mouse host genes keyed by the symbol they share, for crossing the species (instruction 71).

    Returns `species, gene_id, symbol`. Human symbols come from the human gene space (`hs_genes.parquet`), mouse
    ones from the mouse protein-gene links joined to the shipped mouse host table. The symbol is upper-cased, so
    mouse `Bfd1` meets human `BFD1`.

    This is a symbol mapping, not a tested orthology: it cannot resolve a gene family where one species has
    several members, and it says nothing about which member is the functional counterpart. A symbol that names
    more than one gene in a species is dropped rather than guessed, and the count of those is the honest measure
    of what the mapping misses.
    """
    here = os.path.dirname(os.path.dirname(__file__))
    root = root or os.path.join(here, "results")
    human = pd.read_parquet(os.path.join(root, "human_gene_space_2026_11_07", "hs_genes.parquet"),
                            columns=["gene_id", "gene_name"])
    links = pd.read_parquet(os.path.join(root, "mouse_gene_space_2026_11_07", "mm_protein_gene_links.parquet"))
    names = pd.read_parquet(os.path.join(here, "parascape", "data", "mm_host_proteins.parquet"),
                            columns=["host_id", "host_name"])
    mouse = links.merge(names, on="host_id")[["gene_id", "host_name"]].rename(columns={"host_name": "gene_name"})
    frames = []
    for species, frame in (("human", human), ("mouse", mouse)):
        frame = frame.dropna()
        frame = frame[frame.gene_name.astype(str).str.len() > 0]
        frame = frame.assign(symbol=frame.gene_name.astype(str).str.upper(), species=species)
        frame = frame.drop_duplicates(["gene_id", "symbol"])
        counts = frame.symbol.value_counts()
        frames.append(frame[frame.symbol.isin(counts[counts == 1].index)][["species", "gene_id", "symbol"]])
    return pd.concat(frames, ignore_index=True)


def cross_species(store: dict, root: str | None = None) -> dict:
    """The host atlas re-indexed on shared gene symbols, so human and mouse signatures can meet.

    The store's rows are versioned Ensembl ids prefixed `HOST:`; both are reduced to the plain id, mapped to a
    symbol with `host_symbol_map`, and signatures of both species are stacked on that symbol index. A gene with
    no symbol, or a symbol that is not one-to-one in its species, is dropped. The returned store has the same
    shape as `host_store`, with `dropped` recording what the crossing cost.
    """
    mapping = host_symbol_map(root)
    by_species = {s: dict(zip(g.gene_id, g.symbol)) for s, g in mapping.groupby("species")}
    lfc, reg = store["lfc"], store["registry"]
    plain = pd.Index(lfc.index.astype(str)).str.replace("^HOST:", "", regex=True).str.replace(r"\.\d+$", "", regex=True)
    columns, kept, dropped = {}, 0, 0
    for column in lfc.columns:
        species = str(reg.loc[column, "species"])
        symbols = pd.Series(plain.map(by_species.get(species, {})), index=lfc.index)
        values = lfc[column][symbols.notna()]
        index = symbols[symbols.notna()]
        series = values.groupby(index.values).mean()        # a symbol reached twice: the mean of its rows
        columns[column] = series
        kept += len(series)
        dropped += int(symbols.isna().sum())
    crossed = pd.DataFrame(columns)
    return {"lfc": crossed, "registry": reg, "release": store["release"] + " (symbols)",
            "dropped": dict(rows_without_symbol=dropped, symbols=len(crossed))}
