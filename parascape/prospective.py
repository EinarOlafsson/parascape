"""Grade frozen prospective claims against truth that arrived later (instructions 63 and 64.38).

A frozen set is a folder holding `claims.tsv` (gene_id, question, a probability, rank) and
`manifest.json` (the claims SHA-256, strategy and evidence hashes, version, time). `grade` refuses a
set whose file no longer matches its hash -- a ranking edited after the fact grades nothing -- and
never refits anything: it reads the frozen file and the new truth and counts.

Every gene is put in exactly one bin: tested (the new truth measures it; then correct or wrong at the
frozen cut), unmeasured (the new truth is silent on it). Precision is reported at the top-k cuts that
were declared before the truth existed, and AUC over the tested genes.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

#: The top-k cuts a grade reports, fixed here so they cannot be chosen after seeing the truth.
CUTS = (50, 100, 200)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(folder) -> tuple:
    """(claims, manifest) of a frozen set, after checking the claims file against its recorded hash."""
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    path = folder / "claims.tsv"
    if _sha256(path) != manifest.get("claims_sha256"):
        raise ValueError(f"{path} does not match the hash frozen in its manifest; it was changed "
                         "after the freeze and cannot be graded")
    return pd.read_csv(path, sep="\t", dtype={"gene_id": str}), manifest


def grade(folder, truth: pd.Series, score: str = "probability_needed", cuts=CUTS) -> dict:
    """Score a frozen set against later truth (a boolean Series indexed by gene id). Fits nothing."""
    claims, manifest = load(folder)
    truth = pd.Series(truth).dropna().astype(bool)
    truth.index = truth.index.astype(str)
    ranked = claims.sort_values("rank") if "rank" in claims else claims.sort_values(score, ascending=False)
    tested = ranked[ranked["gene_id"].isin(truth.index)]
    y = truth.reindex(tested["gene_id"]).to_numpy(dtype=bool)
    out = {"question": manifest.get("question"), "frozen_at": manifest.get("frozen_at"),
           "claims": int(len(claims)), "tested": int(len(tested)),
           "unmeasured": int(len(claims) - len(tested)), "positives": int(y.sum())}
    for k in cuts:
        top = ranked.head(k)
        hit = top["gene_id"].isin(truth.index)
        n = int(hit.sum())
        right = int(truth.reindex(top.loc[hit, "gene_id"]).sum()) if n else 0
        out[f"top{k}_tested"] = n
        out[f"top{k}_correct"] = right
        out[f"top{k}_precision"] = right / n if n else float("nan")
    if 0 < y.sum() < len(y):
        from sklearn.metrics import roc_auc_score
        out["auc"] = float(roc_auc_score(y, tested[score].to_numpy(dtype=float)))
        out["prevalence"] = float(y.mean())
    else:
        out["auc"] = float("nan")
        out["prevalence"] = float(y.mean()) if len(y) else float("nan")
    return out
