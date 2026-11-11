"""Fetch the supplementary files recorded for a source but missing from this machine (65.04, batch 4).

`verify_supplementary_tables.py` could not check 21 of the 63 supplementary-table sources because their file
is not here. Every one of them records a URL, so "not on this machine" was an assumption about the work, not a
fact about the data: the files can be fetched from what the registry already says.

**Most of the URLs are a Europe PMC `supplementaryFiles` endpoint, which serves a ZIP of everything a paper
deposited**, not the single file. So the archive is downloaded, opened, and the member matching the recorded
basename is taken out. Where the basename is not in the archive the members are listed instead of a guess
being made -- picking "the spreadsheet that looks right" is how the wrong table gets installed.

**Downloads are untrusted.** Each goes into its own new empty directory under `--into`, nothing is executed
from it, and extraction refuses any member whose path escapes its directory. The files are placed where the
registry says they live only after that.

    python scripts/fetch_missing_supplementary.py --into /path/for/downloads [--dry-run]
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
import urllib.request
import tarfile
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")
#: Where a recorded `datasets/...` path is rooted, matching verify_supplementary_tables.py.
DESTINATION = (os.environ.get("STARPLAST_DATA_PARENT")
               or "/media/carruthers/mnt3/claude/toxoplasma_projects")
STALL_SECONDS = 240
AGENT = {"User-Agent": "parascape/0.0.11.0 (dataset provenance check; einar.olafsson@gmail.com)"}


def download(url: str, into: str) -> str:
    """Fetch one URL into its own directory and return the local path."""
    os.makedirs(into, exist_ok=True)
    name = os.path.basename(url.split("?")[0]) or "download"
    if not os.path.splitext(name)[1]:
        name += ".zip"                                   # the Europe PMC endpoint serves an archive
    target = os.path.join(into, name)
    request = urllib.request.Request(url, headers=AGENT)
    # A deadline as well as a socket timeout. One Europe PMC endpoint trickled bytes for ten minutes without
    # ever tripping the timeout, which stalls the whole pass on a single source; a slow server is a reason to
    # move on, not to wait indefinitely.
    deadline = time.monotonic() + STALL_SECONDS
    with urllib.request.urlopen(request, timeout=120) as response, open(target, "wb") as handle:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            handle.write(chunk)
            if time.monotonic() > deadline:
                raise TimeoutError(f"gave up after {STALL_SECONDS}s, {handle.tell():,} bytes")
    return target


def is_office_file(path: str) -> bool:
    """True when a ZIP is actually an .xlsx/.docx, which is a ZIP with an OOXML part inside.

    This caught me out: five sources were reported as "the archive does not contain the recorded file" and I
    wrote them up as a provenance defect. Two were nothing of the sort -- the download WAS the spreadsheet,
    and the code was looking for the file inside itself. An Office document is a ZIP, so `is_zipfile` says
    yes to one, and only the content tells them apart.
    """
    try:
        with zipfile.ZipFile(path) as bundle:
            return "[Content_Types].xml" in bundle.namelist()
    except Exception:
        return False


def member_named(archive: str, basename: str, into: str):
    """Extract the member whose basename matches, or return None and leave the archive alone.

    A member path is checked before it is written: a name containing `..` or an absolute path would
    otherwise be extracted outside the directory it was meant to stay in.
    """
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        wanted = [n for n in names if os.path.basename(n) == basename]
        if not wanted:
            return None, names
        name = wanted[0]
        if os.path.isabs(name) or ".." in name.replace("\\", "/").split("/"):
            return None, names
        bundle.extract(name, into)
        return os.path.join(into, name), names


def dangling_link(path: str) -> str:
    """The first component of `path` that is a symlink going nowhere, or "" when the way is clear."""
    part = os.path.dirname(path)
    while part and part != os.path.dirname(part):
        if os.path.islink(part) and not os.path.exists(part):
            return os.readlink(part)
        part = os.path.dirname(part)
    return ""


def main(argv=None) -> int:
    """Fetch every missing supplementary file and report what landed where."""
    import pandas as pd
    from parascape import datasets as ds

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--into", required=True, help="a directory for the downloads; one per source")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--only", action="append", default=[])
    args = parser.parse_args(argv)

    table = pd.read_csv(os.path.join(OUT, "supplementary_table_verification.tsv"), sep="\t")
    missing = sorted(table.loc[table.verdict == "file not on this machine", "dataset_id"])
    if args.only:
        missing = [k for k in missing if k in args.only]
    registry = {d.key: d for d in ds.REGISTRY}
    rows, placed = [], {}
    for key in missing:
        entry = registry[key]
        url, recorded = (entry.url or "").strip(), str(entry.path or "")
        destination = os.path.join(DESTINATION, recorded)
        basename = os.path.basename(recorded.rstrip("/"))
        if args.dry_run:
            print(f"{key}: would fetch {url[:70]} -> {recorded}")
            continue
        folder = os.path.join(args.into, key)
        try:
            fetched = download(url, folder)
        except Exception as exc:
            rows.append(dict(dataset_id=key, outcome="download failed", detail=f"{type(exc).__name__}: {exc}"))
            print(f"{key:<30}download failed: {type(exc).__name__}", flush=True)
            continue
        if zipfile.is_zipfile(fetched) and not is_office_file(fetched):
            inner, names = member_named(fetched, basename, folder)
            # One level of nesting: several publishers bundle the supplementary ZIP inside the
            # article ZIP (Malaria Journal's MOESM1_ESM.zip, Wellcome Open Research's .tgz).
            if inner is None:
                for nested in [n for n in names if n.lower().endswith((".zip", ".tgz", ".tar.gz"))]:
                    try:
                        with zipfile.ZipFile(fetched) as outer:
                            outer.extract(nested, folder)
                        candidate = os.path.join(folder, nested)
                        if zipfile.is_zipfile(candidate):
                            inner, _ = member_named(candidate, basename, folder)
                        elif tarfile.is_tarfile(candidate):
                            # Wellcome Open Research bundles its supplementary data as .tgz inside the
                            # article archive. A tar member is checked for escaping its directory the same
                            # way a zip member is -- tar has carried absolute paths since long before zip.
                            with tarfile.open(candidate) as ball:
                                hit = next((m for m in ball.getmembers()
                                            if m.isfile() and os.path.basename(m.name) == basename), None)
                                if hit and not os.path.isabs(hit.name) and ".." not in hit.name.split("/"):
                                    ball.extract(hit, folder)
                                    inner = os.path.join(folder, hit.name)
                        if inner is not None:
                            break
                    except Exception:
                        continue
            if inner is None:
                rows.append(dict(dataset_id=key, outcome="recorded file not in the archive",
                                 detail=f"looked for {basename}; archive holds {', '.join(names[:6])}"))
                print(f"{key:<30}{basename} not in the archive", flush=True)
                continue
            fetched = inner
        # Part of this dataset tree is symlinked to another machine, and those links dangle here.
        # `makedirs` on such a path raises, and if that drive were ever mounted it would write into
        # the OTHER machine's checkout. The file stays in its download directory instead, and the
        # mapping below tells the verifier where it actually is. The tree is never repaired blindly.
        link = dangling_link(destination)
        size = os.path.getsize(fetched)
        if link:
            rows.append(dict(dataset_id=key, outcome="fetched, left beside the download",
                             detail=f"{recorded} is behind a symlink to another machine ({link}); "
                                    f"kept at {fetched} ({size:,} bytes)"))
            placed[key] = fetched
            print(f"{key:<30}fetched {size:>11,} bytes (tree symlinked away)", flush=True)
            continue
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        shutil.copyfile(fetched, destination)
        rows.append(dict(dataset_id=key, outcome="placed", detail=f"{recorded} ({size:,} bytes)"))
        placed[key] = destination
        print(f"{key:<30}placed {size:>12,} bytes", flush=True)
    # MERGE, never replace. This file is the only record of where a fetched file ended up when the tree is
    # symlinked away, and rewriting it per run silently forgot everything an earlier run had fetched: two
    # sources that verified an hour ago went back to "file not on this machine" because a later pass happened
    # to time out on them, while their files sat on disk the whole time. An entry is dropped only when the
    # file it points at is gone -- the same rule as `add_deposits.py`, which refuses a merge that loses values.
    import json
    mapping_path = os.path.join(OUT, "supplementary_fetched_paths.json")
    try:
        with open(mapping_path) as handle:
            known = json.load(handle)
    except (OSError, ValueError):
        known = {}
    kept = {k: v for k, v in known.items() if os.path.exists(v)}
    dropped = sorted(set(known) - set(kept))
    merged = {**kept, **placed}
    if merged != known:
        with open(mapping_path, "w") as handle:
            json.dump(merged, handle, indent=1, sort_keys=True)
    if dropped:
        print(f"dropped from the map (file gone): {', '.join(dropped)}")
    print(f"fetched-file map: {len(merged)} entries ({len(placed)} from this run)")
    if rows:
        frame = pd.DataFrame(rows)
        frame.to_csv(os.path.join(OUT, "supplementary_fetch_log.tsv"), sep="\t", index=False)
        print()
        print(frame.outcome.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
