"""Bait x prey evidence from the enrichment PRIDE projects (instruction 104).

For every admitted-by-search project curated as enrichment_ip, proximity_labelling or crosslinking, with >= 2
curated conditions: per condition, each gene's mean log2 Top3 over the condition's raw files minus its mean over
the project's other files (missing = not detected there; genes detected in only one side get the side's 5th
percentile as floor), z-scored within the condition. The condition's top-enriched gene is its putative bait; genes
with z >= 2 are its prey. Pairs are graded against STRING pairs measured in Toxoplasma and StarPath crosslinks.

Writes mnt2/parascape/ms_enrichment/enrichment.parquet (project, condition, gene, log2 enrichment, z),
pairs.parquet (bait, prey, project, condition, z) and grading.csv.

    python scripts/ms/build_ms_enrichment_pairs.py
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OMICS = "/nas_mnt/carruthers/Einar/parascape_omics"
OUT = "/media/carruthers/mnt2/parascape/ms_enrichment"
LABELS = "/media/carruthers/mnt2/parascape/labels/interaction_2026_11_25"
KINDS = {"enrichment_ip", "proximity_labelling", "crosslinking"}
Z_MIN = 2.0


def stem(c):
    return re.sub(r"\.(mzML\.gz|mzML)$", "", os.path.basename(c))


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    D = os.path.join(ROOT, "parascape", "data", "ms_uniform")
    proj = pd.read_csv(os.path.join(D, "ms_curated_projects.tsv"), sep="\t")
    files = pd.read_csv(os.path.join(D, "ms_curated_files.tsv"), sep="\t").fillna("")
    enr, pairs = [], []
    for p in proj.loc[proj["experiment_type"].isin(KINDS), "project"]:
        cond = dict(zip(files.loc[files.project == p, "raw_file"], files.loc[files.project == p, "condition"]))
        for gp in glob.glob(os.path.join(OMICS, "prepared_ms_genes", p, "sage_*", "genes.parquet")):
            G = pd.read_parquet(gp)
            G.columns = [stem(c) for c in G.columns]
            groups = {}
            for c in G.columns:
                if cond.get(c):
                    groups.setdefault(cond[c], []).append(c)
            if len(groups) < 2:
                continue
            floor = np.nanpercentile(G.to_numpy(float), 5)
            for k, cols in groups.items():
                rest = [c for c in G.columns if c not in cols and c in cond]
                a = G[cols].mean(axis=1)
                b = G[rest].mean(axis=1)
                seen = a.notna() | b.notna()
                d = (a.fillna(floor) - b.fillna(floor))[seen]
                if len(d) < 50 or d.std() == 0:
                    continue
                z = (d - d.mean()) / d.std()
                enr.append(pd.DataFrame({"project": p, "condition": k, "gene": d.index, "log2_enrichment": d.values,
                                         "z": z.values}))
                bait = z.idxmax()
                prey = z[(z >= Z_MIN) & (z.index != bait)]
                pairs.append(pd.DataFrame({"bait": bait, "prey": prey.index, "project": p, "condition": k,
                                           "z": prey.values}))
    E = pd.concat(enr, ignore_index=True)
    P = pd.concat(pairs, ignore_index=True)
    E.to_parquet(os.path.join(OUT, "enrichment.parquet"), index=False)
    P.to_parquet(os.path.join(OUT, "pairs.parquet"), index=False)
    key = lambda a, b: tuple(sorted((a, b)))
    ref = {}
    s = pd.read_parquet(os.path.join(LABELS, "pairs_string_measured_in_toxoplasma.parquet"))
    ref["STRING measured"] = {key(a, b) for a, b in zip(s.gene_a, s.gene_b)}
    x = pd.read_parquet(os.path.join(ROOT, "parascape", "data", "crosslink_models.parquet"))
    ref["crosslink"] = {key(a, b) for a, b in zip(x.gene_a, x.gene_b)}
    genes = sorted(set(E.gene))
    rng = np.random.default_rng(0)
    rows = []
    for name, R in ref.items():
        hit = np.mean([key(a, b) in R for a, b in zip(P.bait, P.prey)])
        rand = np.mean([key(a, b) in R for a, b in zip(rng.choice(genes, len(P)), rng.choice(genes, len(P)))])
        # same baits, random prey among detected genes: the fair null
        same = np.mean([key(a, b) in R for a, b in zip(P.bait, rng.choice(genes, len(P)))])
        rows.append({"reference": name, "pairs": len(P), "share_in_reference": hit, "random_pairs": rand,
                     "same_bait_random_prey": same, "fold_over_same_bait": hit / same if same else np.nan})
    Gd = pd.DataFrame(rows)
    Gd.to_csv(os.path.join(OUT, "grading.csv"), index=False)
    print(f"{E.project.nunique()} projects, {E.groupby(['project', 'condition']).ngroups} conditions, "
          f"{len(P):,} bait-prey pairs, {P.bait.nunique()} baits")
    print(Gd.round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
