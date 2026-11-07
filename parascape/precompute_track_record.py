"""The production track record as resumable, fingerprinted partitions (instruction 64.25).

`scripts/build_track_record.py` builds the hold-out record in one pass: an interruption loses everything and any
change rebuilds everything. Here the same work is cut into partitions -- one per (organism, label, method) -- each
written once with a fingerprint of exactly what produced it:

* the organism node table (its bytes), so a changed dataset invalidates that organism only;
* the method and the settings it chose for that label, so a changed setting invalidates that method only;
* the code that computes held-out calls (`track_record`, `strategies`, `strategy_catalog`) and the fold/seed/
  curve constants.

`run` builds what is missing or stale, reuses what is valid (file present, bytes match the journal), records a
failure and carries on with the other partitions, and stops cleanly on cancellation or budget (jobs, seconds,
memory) with every earlier partition kept. The journal is written atomically after every partition, with
runtime, peak memory and bytes per partition. `assemble` concatenates the valid partitions into the shipped
table format. Nothing here changes what a strategy computes.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

from .precompute_jobs import _atomic_json, _peak_rss, _rss

FOLDS_SEED_CURVE = {"folds": 5, "seed": 0, "curve_sizes": (1, 5, 20, 100, 500), "curve_repeats": 4,
                    "curve_seed": 1}
_CODE = ("track_record.py", "strategies.py", "strategy_catalog.py", "precompute_track_record.py")


def _sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def code_identity() -> str:
    """Hash of the source that computes held-out calls, plus the fold/seed/curve constants."""
    here = os.path.dirname(os.path.abspath(__file__))
    parts = [_sha_file(os.path.join(here, name)) for name in _CODE]
    return _sha_bytes(json.dumps([parts, FOLDS_SEED_CURVE], sort_keys=True, default=list).encode())


@dataclass(frozen=True)
class Partition:
    """One (organism, label, method) piece of the record and the identities that produced it."""

    organism: str
    target: str
    strategy: str
    inputs: str          # node-table bytes hash
    settings: str        # the method settings chosen for this label, as canonical JSON
    code: str            # code_identity()

    @property
    def fingerprint(self) -> str:
        """What must be equal for an existing partition to be reused."""
        return _sha_bytes(json.dumps(asdict(self), sort_keys=True).encode())

    @property
    def name(self) -> str:
        """File stem: readable, and unique per fingerprint."""
        safe = lambda s: "".join(c if c.isalnum() or c in "-_" else "_" for c in s)[:60]
        return f"{self.organism}__{safe(self.target)}__{safe(self.strategy)}__{self.fingerprint[:16]}"


@dataclass(frozen=True)
class Budget:
    """Stop cleanly before a partition when any limit is reached."""

    max_jobs: int = 10_000
    max_seconds: float = 24 * 3600.0
    memory_limit_bytes: int = 64 << 30


def plan(organism: str, targets=None, strategies=None) -> list:
    """Partitions for one organism: every recorded label (or `targets`) x every supported method."""
    from . import organisms, strategies as S, track_record as T
    ctx = S.Context.shipped(organism)
    targets = list(targets or T.labels(ctx))
    keys = list(strategies or T.supported(organism))
    inputs, code = _sha_file(organisms.nodes_path(organism)), code_identity()
    out = []
    for target in targets:
        for key in keys:
            try:
                chosen = T._permitted(ctx, dict(S.get(key).settings(ctx, target=target)), target)
            except Exception as error:                 # the method cannot be set up for this label
                chosen = {"unavailable": f"{type(error).__name__}: {error}"}
            out.append(Partition(organism, target, key, inputs,
                                 json.dumps(chosen, sort_keys=True, default=str), code))
    return out


_CONTEXTS: dict = {}


def build_partition(p: Partition) -> pd.DataFrame:
    """The production computation for one partition -- exactly `build_track_record.build` for one method."""
    from . import strategies as S, track_record as T
    ctx = _CONTEXTS.get(p.organism)
    if ctx is None:
        ctx = _CONTEXTS[p.organism] = S.Context.shipped(p.organism)
    c = FOLDS_SEED_CURVE
    folds = T.evaluate(ctx, p.strategy, p.target, folds=c["folds"], seed=c["seed"])
    if not len(folds):
        return pd.DataFrame(columns=list(T.COLUMNS))
    return pd.concat([folds,
                      T.evaluate_sets(ctx, p.strategy, p.target, T.class_sets(ctx, p.target)),
                      T.evaluate_sets(ctx, p.strategy, p.target,
                                      T.random_sets(ctx, p.target, c["curve_sizes"], c["curve_repeats"],
                                                    seed=c["curve_seed"]))], ignore_index=True)


def _journal(directory: Path) -> dict:
    path = directory / "journal.json"
    if not path.exists():
        return {"partitions": {}, "runs": []}
    envelope = json.loads(path.read_text())
    body = envelope.get("journal")
    if not isinstance(body, dict) or _sha_bytes(json.dumps(body, sort_keys=True).encode()) != envelope.get("sha256"):
        raise ValueError("journal is corrupt; move it aside to rebuild")
    return body


def _save(directory: Path, journal: dict) -> None:
    body = json.loads(json.dumps(journal, sort_keys=True))
    _atomic_json(directory / "journal.json",
                 {"journal": body, "sha256": _sha_bytes(json.dumps(body, sort_keys=True).encode())})


def valid(directory, p: Partition, journal: dict | None = None) -> bool:
    """A partition is reusable when its journal record is ok for this fingerprint and its bytes match."""
    directory = Path(directory)
    journal = _journal(directory) if journal is None else journal
    record = journal["partitions"].get(p.name)
    path = directory / "partitions" / f"{p.name}.parquet"
    return bool(record and record.get("status") == "ok" and record.get("fingerprint") == p.fingerprint
                and path.is_file() and _sha_file(path) == record.get("sha256"))


def run(partitions, directory, budget: Budget = Budget(), *, builder=build_partition,
        cancelled=lambda: False, log=lambda _m: None) -> dict:
    """Build missing or stale partitions, reuse valid ones; returns {name: status} for this run.

    Statuses: ``reused``, ``built``, ``failed`` (error recorded, other partitions continue), ``pending``
    (stopped by cancellation or budget before it started).
    """
    directory = Path(directory)
    (directory / "partitions").mkdir(parents=True, exist_ok=True)
    journal = _journal(directory)
    started, built = time.monotonic(), 0
    run_record = {"started_unix": time.time(), "budget": asdict(budget), "partitions": len(partitions),
                  "stop_reason": "", "statuses": {}}
    journal["runs"].append(run_record)
    states = {}
    for p in partitions:
        if valid(directory, p, journal):
            states[p.name] = "reused"
            continue
        reason = ("cancelled" if cancelled() else
                  "job budget" if built >= budget.max_jobs else
                  "time budget" if time.monotonic() - started >= budget.max_seconds else
                  "memory budget" if _rss() > budget.memory_limit_bytes else "")
        if reason:
            run_record["stop_reason"] = reason
            states.update({q.name: "pending" for q in partitions if q.name not in states})
            break
        t0 = time.monotonic()
        record = {"fingerprint": p.fingerprint, "partition": asdict(p)}
        try:
            frame = builder(p)
            path = directory / "partitions" / f"{p.name}.parquet"
            temporary = path.with_suffix(".partial")
            frame.to_parquet(temporary, index=False)
            os.replace(temporary, path)
            record.update(status="ok", rows=int(len(frame)), sha256=_sha_file(path), bytes=path.stat().st_size)
            states[p.name] = "built"
        except Exception as error:
            record.update(status="failed", error=f"{type(error).__name__}: {error}")
            states[p.name] = "failed"
        record.update(seconds=round(time.monotonic() - t0, 3), peak_rss_bytes=_peak_rss())
        journal["partitions"][p.name] = record
        built += 1
        run_record["statuses"] = dict(states)
        _save(directory, journal)
        log(f"{p.name}: {states[p.name]} in {record['seconds']}s")
    run_record.update(statuses=dict(states), seconds=round(time.monotonic() - started, 3))
    _save(directory, journal)
    return states


def assemble(directory, partitions) -> pd.DataFrame:
    """The valid partitions concatenated in plan order, in the shipped table format."""
    directory = Path(directory)
    journal = _journal(directory)
    frames = [pd.read_parquet(directory / "partitions" / f"{p.name}.parquet")
              for p in partitions if valid(directory, p, journal)]
    frames = [f for f in frames if len(f)]
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    for column in ("organism", "strategy", "target", "setting_key", "mode", "set_name", "truth", "prediction"):
        if column in out:
            out[column] = out[column].astype("category")
    return out


def report(directory) -> pd.DataFrame:
    """Per partition: status, rows, seconds, peak memory and bytes, from the journal."""
    journal = _journal(Path(directory))
    rows = [{"partition": k, **{f: v.get(f) for f in ("status", "rows", "seconds", "peak_rss_bytes", "bytes", "error")}}
            for k, v in journal["partitions"].items()]
    return pd.DataFrame(rows)
