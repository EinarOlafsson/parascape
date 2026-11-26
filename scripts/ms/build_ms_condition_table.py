"""Per-condition protein abundance for admitted uniform-MS projects that have an SDRF sample design (instruction 104).

A raw file's condition is the joined `factor value[...]` cells of its SDRF row (PRIDE submitted or community-
annotated design, processing_audit/pride_sample_designs; not independently verified). Writes
parascape/data/ms_uniform/ms_condition_abundance.parquet, long form: gene_id, project, condition, files, mean_log2
(mean of the median-centred log2 Top3 over the condition's raw files), and ms_condition_designs.tsv (one row per
project x condition with its files and the design source).

    python scripts/ms/build_ms_condition_table.py
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from parascape import ms_genes as M  # noqa: E402

OMICS = "/nas_mnt/carruthers/Einar/parascape_omics"
DESIGNS = os.path.join(OMICS, "processing_audit", "pride_sample_designs")


def stem(name) -> str:
    return re.sub(r"\.(mzML\.gz|mzML|raw|RAW|d|wiff|mgf)$", "", os.path.basename(str(name)))


def main() -> int:
    current = json.load(open(os.path.join(DESIGNS, "collection.json")))["directory"]
    adm = pd.read_csv(os.path.join(ROOT, "parascape", "data", M.TABLE_DIR, M.PROJECTS_TSV), sep="\t")
    admitted = set(adm.loc[adm["admitted"], "project"])
    long, designs = [], []
    for f in sorted(glob.glob(os.path.join(DESIGNS, current, "*.json"))):
        j = json.load(open(f))
        if "headers" not in j or "comment[data file]" not in j["headers"]:
            continue
        h = j["headers"]
        di = h.index("comment[data file]")
        fv = [i for i, x in enumerate(h) if x.startswith("factor value")] or \
             [i for i, x in enumerate(h) if x.startswith("characteristics[") and "replicate" not in x]
        cond_of = {stem(r["cells"][di]): " | ".join(str(r["cells"][i]) for i in fv) for r in j["rows"]}
        for p in j.get("primary_file_metadata_sha256", {}):
            if p not in admitted:
                continue
            for gp in sorted(glob.glob(os.path.join(OMICS, "prepared_ms_genes", p, "sage_*", "genes.parquet"))):
                G = pd.read_parquet(gp)
                groups: dict = {}
                for c in G.columns:
                    if stem(c) in cond_of:
                        groups.setdefault(cond_of[stem(c)], []).append(c)
                for cond, cols in groups.items():
                    v = G[cols].mean(axis=1, skipna=True).dropna()
                    long.append(pd.DataFrame({"gene_id": v.index, "project": p, "condition": cond,
                                              "files": len(cols), "mean_log2": v.to_numpy("float32")}))
                    designs.append({"project": p, "instrument": os.path.basename(os.path.dirname(gp)),
                                    "condition": cond, "files": len(cols), "factor_columns": "; ".join(h[i] for i in fv),
                                    "design_source": j.get("source_url", ""), "origin": j.get("metadata_origin", "")})
    out = os.path.join(ROOT, "parascape", "data", M.TABLE_DIR)
    L = pd.concat(long, ignore_index=True)
    L.to_parquet(os.path.join(out, "ms_condition_abundance.parquet"), index=False)
    D = pd.DataFrame(designs)
    D.to_csv(os.path.join(out, "ms_condition_designs.tsv"), sep="\t", index=False)
    print(f"{D['project'].nunique()} projects, {len(D)} conditions, {L['gene_id'].nunique():,} genes, {len(L):,} values")
    return 0


if __name__ == "__main__":
    sys.exit(main())
