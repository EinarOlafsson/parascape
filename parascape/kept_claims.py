"""Your own kept inferences: a layer of claims you accepted, never mixed with the evidence.

Instruction 78. Einar: "if inference is done the user can choose to add it to their own ParaScape
dataset and then it should turn up in `Predicted`."

The rules this module enforces, each because the alternative turns a guess into evidence:

* **A separate store** in the user's own folder (`kept_claims.jsonl`); nothing is written to the
  installed release, so the information map, the gene table and every training set read exactly what
  they read before. A kept claim is never training data because nothing that builds training data
  reads this file.
* **Reproducible or not kept.** Every claim records the value, its probability, the method and its
  version, the question, the score it was accepted at, whether it cleared the frozen-test bar, the
  data release it was computed against, and when.
* **Shown as kept, never as fact**: the Gene module lists kept claims under "What can be inferred".
* **Reverted individually, or every claim of one method version at once.**
* **Superseded loudly.** When a later release MEASURES a kept cell, `check_against` says so and whether
  the measurement agreed -- a kept guess that happens to be right must not become indistinguishable
  from evidence.
* **Exported with its provenance**, so it can be shared or cited.

No window here; every rule is testable.
"""
from __future__ import annotations

import json
import os
import time
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd

#: A kept claim cleared the bar when its frozen-test precision is at least this (instruction 72).
FROZEN_BAR = 0.90


def store_path() -> Path:
    """Where the user's kept claims live: their own cache folder, never the release."""
    override = os.environ.get("PARASCAPE_KEPT_CLAIMS", "")
    if override:
        return Path(override)
    from . import paths
    return Path(paths.user_cache_dir()) / "kept_claims.jsonl"


def release_id() -> str:
    """The data release the application is reading, from its `data_release.json`."""
    try:
        from . import paths
        with open(Path(paths.data_dir()) / "data_release.json") as fh:
            return str(json.load(fh).get("version", "unknown"))
    except Exception:
        return "unknown"


@dataclass
class KeptClaim:
    """One inference the user chose to keep, with everything needed to reproduce it."""

    gene: str
    trait: str                       # the label or column the claim is about
    value: str
    probability: float = float("nan")
    uncertainty: str = ""            # e.g. "precision 0.83 on 412 held-out calls"
    method: str = ""                 # strategy key or module, e.g. "infer:host_facing" or "map neighbours"
    method_version: str = ""
    question: str = ""
    accepted_score: float = float("nan")
    cleared_bar: bool = False
    release: str = field(default_factory=lambda: release_id())
    kept_at: str = field(default_factory=lambda: time.strftime("%Y-%m-%dT%H:%M:%S"))
    claim_id: str = field(default_factory=lambda: uuid.uuid4().hex[:12])
    reverted: bool = False


def _read() -> list:
    path = store_path()
    if not path.exists():
        return []
    out = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    continue
    return out


def _write(records: list) -> None:
    path = store_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with tmp.open("w") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")
    tmp.replace(path)


def keep(claim: KeptClaim) -> KeptClaim:
    """Add a claim to the store. Refused when it cannot be reproduced (no method or no value)."""
    if not claim.gene or not claim.trait or claim.value in ("", None) or not claim.method:
        raise ValueError("a kept claim needs a gene, a trait, a value and the method that made it")
    claim.cleared_bar = bool(claim.accepted_score == claim.accepted_score
                             and claim.accepted_score >= FROZEN_BAR) or claim.cleared_bar
    records = _read()
    records.append(asdict(claim))
    _write(records)
    return claim


def claims(gene: str | None = None, include_reverted: bool = False) -> pd.DataFrame:
    """The kept claims, for one gene or all."""
    frame = pd.DataFrame(_read())
    if not len(frame):
        return pd.DataFrame(columns=list(KeptClaim.__dataclass_fields__))
    if not include_reverted:
        frame = frame[~frame["reverted"].astype(bool)]
    if gene is not None:
        frame = frame[frame["gene"].astype(str) == str(gene)]
    return frame.reset_index(drop=True)


def revert(claim_id: str) -> bool:
    """Withdraw one kept claim (kept in the file, marked reverted, so the history survives)."""
    records = _read()
    hit = False
    for r in records:
        if r.get("claim_id") == claim_id and not r.get("reverted"):
            r["reverted"] = True
            hit = True
    if hit:
        _write(records)
    return hit


def revert_method(method: str, version: str | None = None) -> int:
    """Withdraw every claim of one method (and version, when given). Returns how many."""
    records = _read()
    n = 0
    for r in records:
        if r.get("method") == method and (version is None or r.get("method_version") == version) \
                and not r.get("reverted"):
            r["reverted"] = True
            n += 1
    if n:
        _write(records)
    return n


def check_against(published) -> pd.DataFrame:
    """Kept claims whose cell a release now MEASURES, and whether the measurement agrees.

    `published` is anything with `gene_facts(gene)` (a `PublishedMap`). The result is meant to be
    shown, not acted on silently: a superseded claim is reported as such, agreeing or not.
    """
    kept = claims()
    rows = []
    for c in kept.itertuples(index=False):
        try:
            facts = published.gene_facts(c.gene)
        except Exception:
            continue
        if facts is None or not len(facts):
            continue
        hit = facts[(facts["column"].astype(str) == c.trait) & facts["observed"].astype(bool)]
        tiers = hit["source_tier"].astype(str).str.lower() if len(hit) else hit
        hit = hit[~tiers.str.startswith(("filled", "exploratory", "predicted"))] if len(hit) else hit
        if not len(hit):
            continue
        from .gene_reading import value_text
        measured = value_text(c.trait, hit.iloc[0]["value"])
        kept_text = value_text(c.trait, c.value)
        rows.append(dict(claim_id=c.claim_id, gene=c.gene, trait=c.trait, kept_value=c.value,
                         measured_value=measured, release=c.release,
                         # Read the same way (True and "yes" are one answer), or a right guess
                         # would be reported as contradicted.
                         agrees=measured.strip().lower() == kept_text.strip().lower()))
    return pd.DataFrame(rows, columns=["claim_id", "gene", "trait", "kept_value",
                                       "measured_value", "release", "agrees"])


def export(path: str) -> int:
    """Write every kept claim, with its provenance, to a CSV. Returns how many rows."""
    frame = claims(include_reverted=True)
    frame.to_csv(path, index=False)
    return len(frame)
