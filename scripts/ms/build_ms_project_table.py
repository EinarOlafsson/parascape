"""Gene x project protein abundance from the uniform Sage re-search, admitted to the gene table (instruction 104).

Reads prepared_ms_genes/<project>/<instrument>/genes.parquet (scripts/ms/build_ms_gene_tables.py), keeps folders
with >= ms_genes.MIN_GENES Toxoplasma genes quantified and a Spearman correlation with the existing proteome
(`proteome_log2_normalized*` columns) >= ms_genes.MIN_RHO_PROTEOME, and writes
parascape/data/ms_uniform/ms_project_abundance.parquet (gene_id + one `msu_<project>` column per admitted folder:
mean of the median-centred log2 Top3 over the project's raw files) and ms_project_abundance_projects.tsv (every
folder, admitted or not, with the reason).

    python scripts/ms/build_ms_project_table.py --nodes NODES.parquet
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from parascape import ms_genes as M  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--omics", default="/nas_mnt/carruthers/Einar/parascape_omics")
    ap.add_argument("--nodes", required=True)
    args = ap.parse_args(argv)
    nodes = pd.read_parquet(args.nodes).set_index("gene_id")
    ref = nodes[[c for c in nodes.columns if c.startswith("proteome_log2_normalized")]].mean(axis=1)
    summary = pd.read_csv(os.path.join(args.omics, "prepared_ms_genes", "summary.csv"))
    # Curated experiment types (instruction 104, Sonnet curation 2026-11-26): only abundance-meaningful projects
    # (whole proteome, subcellular fraction, secretome) may enter an abundance column.
    cur = pd.read_csv(os.path.join(ROOT, "parascape", "data", M.TABLE_DIR, "ms_curated_projects.tsv"), sep="\t")
    curated = {r.project: {"experiment_type": r.experiment_type, "abundance_meaningful": bool(r.abundance_meaningful)}
               for r in cur.itertuples()}
    per_project = summary.groupby("project")["instrument"].nunique()
    cols, meta = {}, []
    for f in sorted(glob.glob(os.path.join(args.omics, "prepared_ms_genes", "*", "sage_*", "genes.parquet"))):
        instrument, project = os.path.basename(os.path.dirname(f)), os.path.basename(os.path.dirname(os.path.dirname(f)))
        G = pd.read_parquet(f)
        a = M.project_abundance(G)
        both = pd.concat([a, ref], axis=1).dropna()
        rho = float(spearmanr(both.iloc[:, 0], both.iloc[:, 1])[0]) if len(both) > 50 else np.nan
        n = int(a.notna().sum())
        kind = curated.get(project, {})
        why = (f"curated as {kind.get('experiment_type')}: measures enrichment or identification, not abundance"
               if kind and not kind.get("abundance_meaningful") else
               "" if n >= M.MIN_GENES and (np.isfinite(rho) and rho >= M.MIN_RHO_PROTEOME) else
               f"{n} genes (< {M.MIN_GENES})" if n < M.MIN_GENES else
               f"Spearman {rho:.2f} with the existing proteome (< {M.MIN_RHO_PROTEOME}): enrichment-like")
        name = M.column_name(project, instrument, int(per_project.get(project, 1)))
        meta.append({"column": name, "project": project, "instrument": instrument, "files": G.shape[1],
                     "experiment_type": kind.get("experiment_type", "not curated"),
                     "genes": n, "rho_existing_proteome": rho, "admitted": not why, "why_not": why})
        if not why:
            cols[name] = a
    T = pd.DataFrame(cols).astype("float32")
    T.index.name = "gene_id"
    out = os.path.join(ROOT, "parascape", "data", M.TABLE_DIR)
    os.makedirs(out, exist_ok=True)
    T.reset_index().to_parquet(os.path.join(out, M.TABLE), index=False)
    pd.DataFrame(meta).to_csv(os.path.join(out, M.PROJECTS_TSV), sep="\t", index=False)
    m = pd.DataFrame(meta)
    print(f"{int(m.admitted.sum())} of {len(m)} folders admitted; {T.notna().any(axis=1).sum()} genes; "
          f"not admitted: {m.loc[~m.admitted, 'why_not'].str.split(' ').str[1:3].str.join(' ').value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
