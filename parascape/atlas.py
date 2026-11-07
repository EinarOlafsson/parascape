"""The inference atlas of one organism: what every method has precomputed, and what runs on demand (64.26-64.29).

`manifest(organism)` has one row per catalogue method with its role (gene label, value, ranking, module,
replication, relationship) and, for each kind of output, either where it is stored or a specific status:

* **known genes** -- held-out per-gene rows in the shipped track record (each gene hidden with its orthogroup);
* **unknown genes** -- deployment claims made for genes without the label (`claims.parquet`, by generator);
* **card** -- the calibration scorecard measured on this organism.

Known genes are only ever read from held-out rows and unknown genes only from deployment output, so a gene with
a label is never shown a prediction that saw its own answer. `cached_or_run` serves a query from the atlas cache
by its identity (organism, method, label, genes, settings) and runs only what is missing, storing the result.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd

ROLES = {"label calls": "gene label", "values": "numeric value", "ranking": "gene ranking",
         "cluster recovery": "module", "replication": "replication", "set retrieval": "gene set"}


def _data(name: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name)


def manifest(organism: str, ledger: pd.DataFrame | None = None, claims: pd.DataFrame | None = None,
             calibration=None) -> pd.DataFrame:
    """One row per catalogue method: role, held-out record, deployment claims, card, and status of each."""
    from . import calibration as CAL, strategies as S, track_record as T
    ledger = T.shipped(organism) if ledger is None else ledger
    if claims is None:
        path = _data("claims.parquet")
        claims = pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame(columns=["organism", "generator"])
    claims = claims[claims["organism"].astype(str) == organism]
    entry = calibration or (lambda key: CAL.entry(key, organism))
    together = ledger[ledger["mode"].astype(str) == "together"] if len(ledger) else ledger
    rows = []
    for s in S.catalog():
        held = together[together["strategy"].astype(str) == s.key] if len(together) else together
        labels = sorted(held["target"].astype(str).unique()) if len(held) else []
        made = claims[claims["generator"].astype(str) == s.key]
        card = entry(s.key) or {}
        if len(held):
            known = f"held-out rows for {held['gene_id'].nunique():,} genes on {len(labels)} labels"
        elif s.task != "label calls":
            known = (f"no per-gene label rows: a {ROLES.get(s.task, s.task)} method is graded by its own "
                     "held-out scorecard (card)")
        else:
            known = "not in the per-gene record for this organism (the method cannot speak on its labels)"
        if len(made):
            unknown = f"{len(made):,} deployment claims on {made['target'].astype(str).nunique()} labels"
        elif s.task == "label calls" and len(held):
            unknown = ("on demand: deployment claims are built only for the best calibrated recipe of a label "
                       "(claim_recipes.parquet); run it from Infer")
        else:
            unknown = "on demand: run from its module (results are cached by query identity)"
        rows.append(dict(organism=organism, strategy=s.key, name=s.name, role=ROLES.get(s.task, s.task),
                         task=s.task, known_genes=known, unknown_genes=unknown,
                         heldout_genes=int(held["gene_id"].nunique()) if len(held) else 0,
                         heldout_labels=", ".join(labels), deployment_claims=int(len(made)),
                         card="measured" if card else "not measured on this organism",
                         card_grade=str(card.get("grade", "")) if card else ""))
    return pd.DataFrame(rows)


def query_identity(organism: str, strategy: str, target: str = "", genes=(), settings: dict | None = None) -> str:
    """The cache key of a query: everything that changes its answer, and nothing else."""
    body = {"organism": organism, "strategy": strategy, "target": target or "",
            "genes": sorted(map(str, genes)), "settings": settings or {}}
    return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()


def cache_dir() -> Path:
    """Where on-demand atlas results are kept (per user, beside the data packs)."""
    from . import paths
    return Path(paths.user_cache_dir()) / "atlas"


def cached_or_run(organism: str, strategy: str, run, *, target: str = "", genes=(), settings=None,
                  directory=None) -> tuple:
    """(frame, "cached" | "ran"): the stored answer to this exact query, or `run()` now and store it."""
    directory = Path(directory) if directory is not None else cache_dir()
    key = query_identity(organism, strategy, target, genes, settings)
    path = directory / organism / f"{strategy}__{key[:24]}.parquet"
    if path.is_file():
        return pd.read_parquet(path), "cached"
    frame = run()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".partial")
    frame.to_parquet(temporary, index=False)
    os.replace(temporary, path)
    meta = {"organism": organism, "strategy": strategy, "target": target, "genes": len(list(genes)),
            "settings": settings or {}, "identity": key}
    path.with_suffix(".json").write_text(json.dumps(meta, sort_keys=True, default=str))
    return frame, "ran"


def coverage(m: pd.DataFrame) -> dict:
    """How much of the atlas is precomputed: methods with held-out rows, with claims, with a card."""
    return dict(methods=len(m), with_heldout_rows=int((m["heldout_genes"] > 0).sum()),
                with_deployment_claims=int((m["deployment_claims"] > 0).sum()),
                with_card=int((m["card"] == "measured").sum()),
                every_method_has_a_status=bool(m[["known_genes", "unknown_genes"]].notna().all().all()))
