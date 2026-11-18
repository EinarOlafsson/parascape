#!/usr/bin/env python3
"""Publish instruction-53 organism and host packs to Hugging Face and record their checksums.

Instruction 99 step 1 (Einar, 2026-11-13: all ParaScape data is hosted on Hugging Face). A pack used
to be fetched by typing its address and its checksum by hand. This uploads each pack to
``packs/<organism>/<version>/`` in the dataset repository and writes the checksum into
``parascape/pack_index.json``, which ships inside the wheel, so `packs.download_pack()` never has to
be told one.

The index is written from the files themselves, never from a value passed in: a checksum that was
typed is a checksum that can be typed wrong.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent.parent
INDEX = ROOT / "parascape" / "pack_index.json"
#: Both spellings a pack has ever had. `build_space_packs.py` writes
#: `parascape_<organism>_<version>.zip`; an earlier hand-made pack was
#: `<organism>-<version>.zip`. The two scripts had never been run against each other, so
#: publishing a freshly built pack failed on its own name (found 2026-11-18).
NAME = re.compile(r"^(?:parascape[_-])?(?P<organism>[A-Za-z]{2})[-_](?P<version>[0-9][0-9.]*)\.zip$")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_index() -> dict:
    if not INDEX.exists():
        return {"format": 1, "packs": {}}
    raw = json.loads(INDEX.read_text())
    if raw.get("format") != 1:
        raise SystemExit(f"{INDEX} is not a format 1 index")
    return raw


def record(archives: list[Path]) -> dict:
    """Add each pack to the index, keyed by organism and version, with its own measured checksum."""
    index = read_index()
    for archive in archives:
        m = NAME.match(archive.name)
        if not m:
            raise SystemExit(f"{archive.name}: expected <organism>-<version>.zip")
        organism, version = m["organism"], m["version"]
        entry = index["packs"].setdefault(organism, {})
        entry[version] = {
            "sha256": digest(archive),
            "size": archive.stat().st_size,
            "path": f"{organism}/{version}/{archive.name}",
        }
        print(f"recorded {organism} {version}: {entry[version]['sha256'][:12]}… "
              f"{entry[version]['size']:,} bytes")
    INDEX.write_text(json.dumps(index, indent=2, sort_keys=True) + "\n")
    return index


def upload(archives: list[Path], repo: str, token: str | None) -> None:
    """Upload each pack to the dataset repository, under the path the index records."""
    from huggingface_hub import HfApi
    api = HfApi(token=token)
    for archive in archives:
        m = NAME.match(archive.name)
        path = f"packs/{m['organism']}/{m['version']}/{archive.name}"
        api.upload_file(path_or_fileobj=str(archive), path_in_repo=path, repo_id=repo,
                        repo_type="dataset")
        print(f"uploaded {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("archives", nargs="+", type=Path)
    ap.add_argument("--repo", default="einarolafsson/parascape-data")
    ap.add_argument("--token", default=None)
    ap.add_argument("--record-only", action="store_true",
                    help="write the index without uploading, for a dry run or a re-record")
    args = ap.parse_args(argv)
    missing = [str(a) for a in args.archives if not a.is_file()]
    if missing:
        raise SystemExit("no such pack: " + ", ".join(missing))
    if not args.record_only:
        upload(args.archives, args.repo, args.token)
    record(args.archives)
    print(f"index written: {INDEX.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
