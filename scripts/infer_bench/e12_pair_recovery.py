"""Infer E12: do genes that interact have similar profiles? Pair recovery for every kind of known pair.

Instruction 103 (use the STRING / CORUM pair tables for binding partners). For each pair set:
  - StarPath XL-MS crosslinks (direct physical evidence; `parascape/data/crosslink_models.parquet`);
  - STRING pairs measured in Toxoplasma (experiments or curated database, nothing transferred);
  - STRING pairs copied by orthology (transferred score >= 700);
  - CORUM human complexes copied by orthology (all pairs within a mapped complex),
the cosine similarity of the two genes' profiles is compared with DEGREE-MATCHED decoys: each positive pair
(a, b) gets decoys (a, b') where b' has the same number of known partners as b (same degree bin) and no
recorded interaction with a in any set. Reported: AUROC and the share of positives in the top 1 % of all
scored pairs (enrichment over the decoys), for two embeddings - the E8 reference (per view PCA 30) and the
E10 default (no per view PCA, whitening) - and a negative control with gene rows shuffled.
The interaction view (counts of studies and partners) is LEFT OUT of the profiles: it measures how often a gene
was pulled down, which tracks degree, not who its partner is.

Writes pair_scores.parquet (every positive and decoy with its similarity), summary.csv.

    python scripts/infer_bench/e12_pair_recovery.py OUT_DIR
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import e0_harness as H  # noqa: E402
import e10_ops_ledger as E10  # noqa: E402

LABELS = "/media/carruthers/mnt2/parascape/labels/interaction_2026_11_25"
DECOYS = 5
SEED = 12


def pair_sets(idx: dict) -> dict:
    out = {}
    xl = pd.read_parquet(os.path.join(ROOT, "parascape", "data", "crosslink_models.parquet"))
    out["crosslink (XL-MS, measured)"] = xl[["gene_a", "gene_b"]]
    s = pd.read_parquet(os.path.join(LABELS, "pairs_string_measured_in_toxoplasma.parquet"))
    out["STRING, measured in Toxoplasma"] = s[["gene_a", "gene_b"]]
    s = pd.read_parquet(os.path.join(LABELS, "pairs_string_copied_by_orthology.parquet"))
    out["STRING, copied by orthology"] = s.loc[s["combined_score"] >= 700, ["gene_a", "gene_b"]]
    c = pd.read_parquet(os.path.join(LABELS, "pairs_corum_copied_by_orthology.parquet"))
    out["CORUM, copied by orthology"] = c[["gene_a", "gene_b"]]
    clean = {}
    for k, t in out.items():
        a = t["gene_a"].astype(str).map(idx); b = t["gene_b"].astype(str).map(idx)
        ok = a.notna() & b.notna() & (a != b)
        P = np.unique(np.sort(np.c_[a[ok].astype(int), b[ok].astype(int)], axis=1), axis=0)
        clean[k] = P
    return clean


def decoys(P: np.ndarray, known: set, degree: np.ndarray, eligible: np.ndarray, rng) -> np.ndarray:
    """For each positive (a, b): DECOYS pairs (a, b') with b' in b's degree bin and no known pair with a."""
    bins = np.digitize(degree, [1, 2, 4, 8, 16, 32, 64])
    by_bin = {k: eligible[bins[eligible] == k] for k in np.unique(bins[eligible])}
    out = []
    for a, b in P:
        pool = by_bin.get(bins[b], eligible)
        for _ in range(DECOYS):
            for _try in range(20):
                c = int(pool[rng.integers(len(pool))])
                if c != a and (min(a, c), max(a, c)) not in known:
                    out.append((a, c)); break
    return np.array(out)


def cosine_pairs(E: np.ndarray, P: np.ndarray) -> np.ndarray:
    En = E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)
    return np.einsum("ij,ij->i", En[P[:, 0]], En[P[:, 1]])


def main(argv=None) -> int:
    out = (argv if argv is not None else sys.argv[1:])[0]
    os.makedirs(out, exist_ok=True)
    nodes, views, extra, _ = H.load()
    views = {v: c for v, c in views.items() if v != "interaction"}       # degree, not partner identity
    idx = {g: i for i, g in enumerate(nodes["gene_id"].astype(str))}
    sets = pair_sets(idx)
    n = len(nodes)
    known = set()
    for P in sets.values():
        known |= set(map(tuple, P))
    degree = np.zeros(n, int)
    for a, b in known:
        degree[a] += 1; degree[b] += 1
    embeds = {}
    for name, arm in (("E8 reference", "baseline"), ("E10 default", "no per-view PCA + whitened")):
        raw = E10.views_for(nodes, views, extra, arm)
        E, rows = E10.join(E10.blocks_for(raw, arm), "whitened" if "whitened" in arm else arm)
        full = np.full((n, E.shape[1]), np.nan); full[rows] = E
        embeds[name] = full
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(n)
    embeds["E10 default, genes shuffled (control)"] = embeds["E10 default"][perm]
    rows_out, scores = [], []
    for sname, P in sets.items():
        ok_rows = np.flatnonzero(~np.isnan(embeds["E10 default"]).any(axis=1)
                                 & ~np.isnan(embeds["E8 reference"]).any(axis=1))
        okset = set(ok_rows.tolist())
        Pk = np.array([p for p in P if p[0] in okset and p[1] in okset])
        D = decoys(Pk, known, degree, ok_rows, np.random.default_rng(SEED + 1))
        for ename, E in embeds.items():
            Eu = np.nan_to_num(E)
            sp, sd = cosine_pairs(Eu, Pk), cosine_pairs(Eu, D)
            y = np.r_[np.ones(len(sp)), np.zeros(len(sd))]
            s = np.r_[sp, sd]
            top = s >= np.quantile(s, 0.99)
            rows_out.append({"pairs": sname, "embedding": ename, "positives": len(sp), "decoys": len(sd),
                             "auroc": roc_auc_score(y, s), "median_cos_pos": float(np.median(sp)),
                             "median_cos_decoy": float(np.median(sd)),
                             "top1pct_positive_share": float(y[top].mean()),
                             "top1pct_enrichment": float(y[top].mean() / y.mean())})
            scores.append(pd.DataFrame({"pairs": sname, "embedding": ename, "gene_a": np.r_[Pk[:, 0], D[:, 0]],
                                        "gene_b": np.r_[Pk[:, 1], D[:, 1]], "positive": y, "cosine": s}))
            print(rows_out[-1], flush=True)
    pd.concat(scores).to_parquet(os.path.join(out, "pair_scores.parquet"), index=False)
    S = pd.DataFrame(rows_out)
    S.to_csv(os.path.join(out, "summary.csv"), index=False)
    print(S.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
