"""Packed record folders: many small files kept as one Parquet table, with their exact bytes.

Instruction 102 (Einar, 2026-11-25: "Parquet first, no DuckDB yet"). A data release held 10,934 files,
7,432 of them in `processing_audit/`, and about 3,300 of those in folders that no module of the
application reads -- provenance kept for audit, not for use. Each such folder now ships as ONE
`<folder>.pack.parquet` with a row per file:

    path    the file path inside the folder (posix)
    data    the file's exact bytes
    sha256  the hex SHA-256 of those bytes

The bytes are kept exactly, not parsed, because the audit records are content-addressed: their
names, receipts and manifests carry SHA-256 digests of the bytes, and a re-serialised JSON would
fail every one of those checks. `unpack` restores a folder byte for byte.
"""
from __future__ import annotations

import hashlib
import io
import re
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

SUFFIX = ".pack.parquet"


def pack_bytes(folder: Path, keep=lambda p: True) -> bytes:
    """The pack of one folder as Parquet bytes: one row per file, sorted by path (deterministic)."""
    folder = Path(folder)
    files = sorted(p for p in folder.rglob("*") if p.is_file() and keep(p))
    paths = [p.relative_to(folder).as_posix() for p in files]
    data = [p.read_bytes() for p in files]
    table = pa.table({"path": pa.array(paths, pa.string()),
                      "data": pa.array(data, pa.binary()),
                      "sha256": pa.array([hashlib.sha256(d).hexdigest() for d in data], pa.string())})
    sink = io.BytesIO()
    pq.write_table(table, sink, compression="zstd")
    return sink.getvalue()


def read(pack: Path) -> pd.DataFrame:
    """The rows of a pack, each checked against its digest."""
    t = pq.read_table(pack).to_pandas()
    bad = [p for p, d, s in zip(t["path"], t["data"], t["sha256"]) if hashlib.sha256(d).hexdigest() != s]
    if bad:
        raise ValueError(f"{pack}: {len(bad)} files do not match their digest, e.g. {bad[:3]}")
    return t


def file_bytes(pack: Path, path: str) -> bytes:
    """The exact bytes of one file of a pack."""
    t = read(pack)
    hit = t[t["path"] == path]
    if not len(hit):
        raise FileNotFoundError(f"{path} is not in {pack}")
    return bytes(hit["data"].iloc[0])


def unpack(pack: Path, dest: Path) -> int:
    """Write every file of a pack under `dest`, byte for byte; returns how many."""
    dest = Path(dest)
    t = read(pack)
    for path, data in zip(t["path"], t["data"]):
        if path.startswith("/") or ".." in Path(path).parts:
            raise ValueError(f"{pack}: unsafe path {path!r}")
        out = dest / path
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_bytes(bytes(data))
    return len(t)


#: Audit folders whose every reader goes through `read_bytes`, so the release may ship them packed.
#: A folder is added here only together with moving all its readers (tests/test_store.py checks it).
READ_THROUGH_STORE = ("all_geo_samples", "pride_projects", "protein_column_design_links")

_PACKS: dict = {}


def _pack_for(path: Path):
    """(pack file, path inside the pack) for a file under a packed folder, or None."""
    for parent in path.parents:
        if not parent.name:                      # the filesystem root has no sibling to look for
            break
        pack = parent.with_name(parent.name + SUFFIX)
        if pack.is_file():
            return pack, path.relative_to(parent).as_posix()
    return None


def read_bytes(path) -> bytes:
    """A file's exact bytes: from disk when it is there, else from the pack of its folder.

    The pack is read once per process and checked against its digests then; a file that is
    neither on disk nor in a pack raises FileNotFoundError, as `Path.read_bytes` would.
    """
    path = Path(path)
    if path.is_file():
        return path.read_bytes()
    found = _pack_for(path)
    if found is None:
        raise FileNotFoundError(path)
    pack, inner = found
    key = str(pack.resolve())
    if key not in _PACKS:
        t = read(pack)
        _PACKS[key] = dict(zip(t["path"], (bytes(d) for d in t["data"])))
    try:
        return _PACKS[key][inner]
    except KeyError:
        raise FileNotFoundError(path) from None


def exists(path) -> bool:
    """Whether a file is on disk or in the pack of its folder."""
    try:
        read_bytes(path)
        return True
    except FileNotFoundError:
        return False


def read_text(path, encoding: str = "utf-8") -> str:
    """A file's text, from disk or from the pack of its folder."""
    return read_bytes(path).decode(encoding)


def unread_folders(audit: Path, code: Path) -> list:
    """Sub-folders of `audit` whose name appears in no module under `code`: nothing reads them."""
    text = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in Path(code).rglob("*.py"))
    out = []
    for folder in sorted(p for p in Path(audit).iterdir() if p.is_dir() and not p.name.startswith(".")):
        if not re.search(r"(?<![\w-])" + re.escape(folder.name) + r"(?![\w-])", text):
            out.append(folder)
    return out
