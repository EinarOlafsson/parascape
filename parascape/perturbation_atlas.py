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
