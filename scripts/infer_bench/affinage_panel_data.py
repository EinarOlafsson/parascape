"""Data for the figure panels taken from Affinage (Di Bernardo and Cheeseman 2026; Einar, 2026-11-24).

Writes tables only; the graph forms are chosen with Einar before anything is drawn. Three tables:

1. `named_vs_evidence.tsv`: genes ToxoDB names against genes it calls "hypothetical protein",
   crossed with the kind of OBSERVED evidence ParaScape holds for them: a paper showing a finding by
   experiment (the literature findings), a hyperLOPIT place, an EC number, a fibroblast fitness
   score, an imaging phenotype. Shows whether a product name stands for knowledge, or only for a
   resemblance.
2. `uncharacterised_census.tsv`: genes with no observed function from any source (no finding in a
   paper, no EC number, no imaging phenotype, no hyperLOPIT place), split by named / hypothetical and
   by whether losing them costs growth in fibroblasts. To be recounted at every release.
3. `findings_by_year.tsv`: the literature findings per publication year and kind; when the field
   established what it knows.

    python scripts/infer_bench/affinage_panel_data.py OUT_DIR [NODES_PARQUET]
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from parascape import literature_findings as LF  # noqa: E402

NODES = "/media/carruthers/mnt1/parascape_data/data-0.0.13.1-hostslots/nodes.parquet"
#: The fibroblast fitness score below which losing a gene is taken to cost growth (Sidik 2016 used -2
#: as a conservative cut for "fitness-conferring" genes in its own discussion).
FITNESS_CUT = -2.0


def evidence(nodes: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    """One row per gene: named or hypothetical, and each kind of observed evidence."""
    g = nodes[["gene_id", "product"]].copy()
    g["named"] = ~g["product"].fillna("").str.contains("hypothetical", case=False)
    lit = labels.set_index("gene_id")
    g["paper_finding"] = g["gene_id"].map(lit["lit_n_papers"]).fillna(0).astype(int) > 0
    g["paper_function"] = (g["gene_id"].map(lit[["lit_n_activity", "lit_n_phenotype", "lit_n_host_facing",
                                                 "lit_n_complex"]].sum(axis=1)).fillna(0) > 0)
    comp = nodes["compartment"].astype(str)
    g["hyperlopit_place"] = comp.ne("unassigned") & comp.ne("nan") & comp.ne("None")
    g["ec_number"] = nodes["ec_number"].notna() & nodes["ec_number"].astype(str).ne("")
    g["fitness_measured"] = nodes["fit_invitro_hff"].notna()
    g["needed_in_fibroblasts"] = nodes["fit_invitro_hff"] < FITNESS_CUT
    g["imaging_phenotype"] = nodes["screen_any_phenotype"].notna()
    g["observed_function"] = g["paper_function"] | g["ec_number"] | g["imaging_phenotype"]
    g["any_observed"] = g["observed_function"] | g["paper_finding"] | g["hyperlopit_place"]
    return g


def named_vs_evidence(g: pd.DataFrame) -> pd.DataFrame:
    """Per named / hypothetical: how many genes, and how many carry each kind of observed evidence."""
    kinds = ["paper_finding", "paper_function", "hyperlopit_place", "ec_number", "fitness_measured",
             "imaging_phenotype", "observed_function", "any_observed"]
    rows = []
    for named, part in g.groupby("named"):
        for k in kinds:
            rows.append({"product": "named" if named else "hypothetical protein", "evidence": k,
                         "genes": int(part[k].sum()), "of": len(part),
                         "percent": round(100 * part[k].mean(), 1)})
    return pd.DataFrame(rows)


def census(g: pd.DataFrame) -> pd.DataFrame:
    """Genes with no observed function, by product name and by need in fibroblasts."""
    none = g[~g["observed_function"]]
    rows = []
    for named, part in none.groupby("named"):
        for need, sub in part.groupby(part["needed_in_fibroblasts"]):
            rows.append({"product": "named" if named else "hypothetical protein",
                         "needed_in_fibroblasts": bool(need), "genes": len(sub),
                         "with_a_place": int(sub["hyperlopit_place"].sum())})
    out = pd.DataFrame(rows)
    out.loc[len(out)] = {"product": "all", "needed_in_fibroblasts": None, "genes": len(none),
                         "with_a_place": int(none["hyperlopit_place"].sum())}
    return out


def by_year(findings: pd.DataFrame) -> pd.DataFrame:
    """Admitted findings per publication year and kind."""
    f = findings[findings["named_in_abstract"].astype(bool)].dropna(subset=["year"])
    f = f.assign(year=f["year"].astype(int))
    t = f.groupby(["year", "kind"]).size().unstack(fill_value=0)
    t["all"] = t.sum(axis=1)
    t["genes_first_shown"] = f.groupby("gene_id")["year"].min().value_counts().reindex(t.index).fillna(0).astype(int)
    return t.reset_index()


def main(argv=None) -> int:
    """Write the three tables and print a summary."""
    args = sys.argv[1:] if argv is None else argv
    if not 1 <= len(args) <= 2:
        raise SystemExit(__doc__)
    out = args[0]
    os.makedirs(out, exist_ok=True)
    nodes = pd.read_parquet(args[1] if len(args) > 1 else NODES)
    g = evidence(nodes, LF.labels(log=lambda *a: None).reset_index())
    a, c, y = named_vs_evidence(g), census(g), by_year(LF.findings())
    a.to_csv(os.path.join(out, "named_vs_evidence.tsv"), sep="\t", index=False)
    c.to_csv(os.path.join(out, "uncharacterised_census.tsv"), sep="\t", index=False)
    y.to_csv(os.path.join(out, "findings_by_year.tsv"), sep="\t", index=False)
    print(a.pivot(index="evidence", columns="product", values="percent").to_string())
    print(c.to_string(index=False))
    print(y.tail(8).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
