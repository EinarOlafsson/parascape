"""Versioned, hash-checked data packs stored outside the installed application.

Packs contain data only. Installation validates every member before atomically changing the
active version. Downloading requires an independently supplied SHA256; the embedded manifest
checks integrity, not the publisher's identity. No network access occurs at import or lookup.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import tempfile
from urllib.parse import urlparse
from urllib.request import urlopen
import zipfile

from . import paths

#: Where published packs live, in the same Hugging Face dataset repository as the release archives,
#: under ``packs/<organism>/<version>/<organism>-<version>.zip``. Instruction 99 step 1, 2026-11-17:
#: a pack had to be fetched by typing its address and its checksum by hand, which meant the checksum
#: a user pasted was the only thing standing between them and someone else's file.
PACKS_HOST = ('https://huggingface.co/datasets/einarolafsson/parascape-data/resolve/main/packs')

#: The two settings the release archives already have, for packs: an alternative base address, and a
#: single complete address that overrides it.
ENV_PACKS_HOST = 'PARASCAPE_PACKS_HOST'
ENV_PACKS_URL = 'PARASCAPE_PACKS_URL'

#: Hosts a pack may be fetched from. The configured host is added to this set, so setting
#: ``$PARASCAPE_PACKS_HOST`` is enough; an address on any other host is refused before any request.
PACK_HOSTS = ('huggingface.co', 'github.com')

#: The checksums of the published packs, shipped inside the wheel so that `download_pack` never has
#: to be told one. Written by ``scripts/publish_packs.py`` when a pack is published.
PACK_INDEX = "pack_index.json"

MAX_ARCHIVE_BYTES = 512 * 1024 * 1024
MAX_UNPACKED_BYTES = 1024 * 1024 * 1024
MANIFEST = "manifest.json"

#: Licence codes a distributed pack may carry for a column.
#:
#: ``O``  an open licence, verified.
#: ``SA`` share-alike: redistribution is allowed and the licence travels with it.
#: ``G``  a public deposit, redistributed with its accession and the paper named. Einar decided this
#:        on 2026-11-13 ("yes, with citation"), which closed a question open since 2026-11-09: a GEO
#:        submission states no licence of its own, and ``O`` would have been a claim nobody made.
#:
#: ``V`` (verify the terms first) and ``X`` (local build only) are not distributable and are refused.
PACK_LICENSES = ("O", "SA", "G")

#: The code whose columns must each name where they came from. It is accepted ONLY with both an
#: accession and a citation recorded, because the permission Einar gave is to redistribute a deposit
#: WITH its attribution: a G column with no accession is an unlicensed column with a nicer letter.
CITED_DEPOSIT = "G"

#: An accession as the public archives write them: GSE119855, PXD001647, E-MTAB-3579, phs000424.
ACCESSION = re.compile(r"[A-Za-z][A-Za-z0-9-]{1,19}[0-9]{3,}")


def _component(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}", value):
        raise ValueError(f"unsafe pack component: {value!r}")
    return value


def _filename(value):
    _component(value)
    if value == MANIFEST or Path(value).suffix not in (".parquet", ".npz", ".json", ".tsv", ".gz", ".txt"):
        raise ValueError(f"unsupported pack data file: {value!r}")
    return value


def digest(path) -> str:
    """SHA256 of a file, read in bounded chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def root(code: str) -> Path:
    """Writable directory containing this organism's versions and active pointer."""
    return Path(paths.user_cache_dir()) / "spaces" / _component(code)


def _manifest(raw):
    m = json.loads(raw)
    if not isinstance(m, dict) or m.get("format") != 1:
        raise ValueError("unsupported pack manifest")
    _component(m["organism"])
    _component(m["version"])
    if not isinstance(m.get("files"), dict) or not m["files"]:
        raise ValueError("pack has no files")
    total = 0
    for name, record in m["files"].items():
        _filename(name)
        if (not isinstance(record, dict) or not isinstance(record.get("size"), int)
                or record["size"] < 0 or not re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256", "")))):
            raise ValueError(f"invalid file record: {name}")
        total += record["size"]
    if total > MAX_UNPACKED_BYTES:
        raise ValueError("pack exceeds unpacked size limit")
    for field in ("nodes", "graph"):
        if m.get(field) and m[field] not in m["files"]:
            raise ValueError(f"pack is missing its {field} file")
    if not m.get("nodes"):
        raise ValueError("pack needs a node table")
    licences = m.get("column_licenses")
    if (not isinstance(licences, dict) or not licences
            or any(v not in PACK_LICENSES for v in licences.values())):
        raise ValueError("distributable columns need verified O, SA or G licenses")
    _check_citations(licences, m.get("column_sources") or {})
    return m


def _check_citations(licences: dict, sources) -> None:
    """A cited-deposit column is accepted only when both its accession and its paper are recorded.

    Checked in the manifest rather than only in the builder, because the manifest is what an
    installed pack is read back from: a column that lost its attribution between building and
    installing would otherwise be redistributed with nothing saying whose data it is.
    """
    if not isinstance(sources, dict):
        raise ValueError("column_sources must map a column to its accession and citation")
    for column in sorted(c for c, code in licences.items() if code == CITED_DEPOSIT):
        record = sources.get(column)
        if not isinstance(record, dict):
            raise ValueError(f"cited public deposit column has no recorded source: {column}")
        accession = str(record.get("accession") or "").strip()
        citation = str(record.get("citation") or "").strip()
        if not ACCESSION.fullmatch(accession):
            raise ValueError(f"cited public deposit column needs an accession: {column}")
        if len(citation) < 10:
            raise ValueError(f"cited public deposit column needs a citation: {column}")


def build(space, version: str, files: dict, column_licenses: dict, destination,
          column_sources: dict | None = None) -> str:
    """Create a data archive and return its SHA256; validate tables before calling this function.

    `column_sources` maps a column to `{"accession": ..., "citation": ...}`. Every column licensed
    `G` (a public deposit) needs one, or the build is refused.
    """
    destination = Path(destination)
    manifest = {"format": 1, "organism": space.code, "version": version,
                "nodes": space.nodes, "graph": space.graph, "column_licenses": column_licenses,
                "column_sources": dict(column_sources or {}),
                "files": {name: {"size": Path(p).stat().st_size, "sha256": digest(p)}
                          for name, p in files.items()}}
    raw = json.dumps(manifest, sort_keys=True, indent=2)
    _manifest(raw)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as fh:
        temporary = Path(fh.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(MANIFEST, raw)
            for name, path in sorted(files.items()):
                archive.write(path, name)
        verify(temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return digest(destination)


def verify(archive_path, expected_sha256: str | None = None) -> dict:
    """Validate the entire archive, including paths, member types, sizes and checksums."""
    if Path(archive_path).stat().st_size > MAX_ARCHIVE_BYTES:
        raise ValueError("pack exceeds archive size limit")
    if expected_sha256 is not None:
        if not re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256):
            raise ValueError("expected SHA256 must have 64 hexadecimal characters")
        if digest(archive_path) != expected_sha256.lower():
            raise ValueError("pack SHA256 mismatch")
    with zipfile.ZipFile(archive_path) as archive:
        infos = archive.infolist()
        names = [i.filename for i in infos]
        if len(names) > 10000 or len(set(names)) != len(names) or MANIFEST not in names:
            raise ValueError("duplicate, excessive or missing pack members")
        if archive.getinfo(MANIFEST).file_size > 1024 * 1024:
            raise ValueError("pack manifest exceeds size limit")
        manifest = _manifest(archive.read(MANIFEST))
        if set(names) != set(manifest["files"]) | {MANIFEST}:
            raise ValueError("pack members do not match manifest")
        for info in infos:
            kind = stat.S_IFMT(info.external_attr >> 16)
            if info.is_dir() or kind not in (0, stat.S_IFREG):
                raise ValueError("pack members must be regular files")
            if info.filename == MANIFEST:
                continue
            record = manifest["files"][info.filename]
            if info.file_size != record["size"]:
                raise ValueError(f"pack size mismatch: {info.filename}")
            h = hashlib.sha256()
            with archive.open(info) as fh:
                for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                    h.update(chunk)
            if h.hexdigest() != record["sha256"]:
                raise ValueError(f"pack member SHA256 mismatch: {info.filename}")
    return manifest


def installed(code: str, version: str | None = None) -> Path | None:
    """Resolve an installed version, or None; a malformed active pointer fails explicitly."""
    base = root(code)
    if version is None:
        pointer = base / "active.json"
        if not pointer.exists():
            return None
        version = json.loads(pointer.read_text())["version"]
    directory = base / _component(version)
    if not directory.is_dir():
        return None
    m = _manifest((directory / MANIFEST).read_text())
    if m["organism"] != code or m["version"] != version:
        raise ValueError("installed pack identity mismatch")
    if any(not (directory / name).is_file() for name in m["files"]):
        raise ValueError("installed pack is incomplete")
    return directory


def install(archive_path, expected_sha256: str | None = None, organism: str | None = None) -> Path:
    """Verify and install a version, then atomically activate it; never overwrite a version."""
    # Copy before checking to ensure extraction reads the exact bytes that were verified.
    state = Path(paths.user_cache_dir())
    state.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pack-", dir=state) as work:
        snapshot = Path(work) / "download.zip"
        if Path(archive_path).stat().st_size > MAX_ARCHIVE_BYTES:
            raise ValueError("pack exceeds archive size limit")
        shutil.copyfile(archive_path, snapshot)
        m = verify(snapshot, expected_sha256)
        if organism is not None and m["organism"] != organism:
            raise ValueError("pack organism mismatch")
        base = root(m["organism"])
        base.mkdir(parents=True, exist_ok=True)
        destination = base / m["version"]
        with tempfile.TemporaryDirectory(prefix="stage-", dir=base) as stage:
            staged = Path(stage) / "data"
            staged.mkdir()
            with zipfile.ZipFile(snapshot) as archive:
                for name in archive.namelist():
                    with archive.open(name) as src, open(staged / name, "wb") as dst:
                        shutil.copyfileobj(src, dst, 1024 * 1024)
            if destination.exists():
                existing = _manifest((destination / MANIFEST).read_text())
                if existing != m or any(digest(destination / n) != r["sha256"]
                                        for n, r in m["files"].items()):
                    raise ValueError("installed version differs; choose a new version")
            else:
                os.rename(staged, destination)
            pointer = Path(stage) / "active.json"
            pointer.write_text(json.dumps({"version": m["version"]}) + "\n")
            os.replace(pointer, base / "active.json")
    return destination


def packs_host() -> str:
    """The configured base address for packs, without its trailing slash."""
    return (os.environ.get(ENV_PACKS_HOST) or PACKS_HOST).rstrip('/')


def _checked_pack_url(url: str) -> str:
    """Refuse an address that is not https on a known host before anything is requested."""
    parsed = urlparse(url)
    allowed = set(PACK_HOSTS) | {urlparse(packs_host()).hostname or ''}
    if parsed.scheme != 'https' or (parsed.hostname or '') not in allowed:
        raise ValueError(f"refusing to fetch a pack from {url!r}: https on a known host only")
    return url


def pack_index() -> dict:
    """The published packs this build knows about: organism -> version -> {sha256, size, path}.

    Shipped with the application. An empty index is normal in a source checkout that has published
    nothing yet, and makes `download_pack` refuse rather than guess."""
    here = Path(__file__).resolve().parent / PACK_INDEX
    if not here.exists():
        return {}
    try:
        raw = json.loads(here.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"the pack index shipped with this build is unreadable: {exc}") from exc
    if not isinstance(raw, dict) or raw.get("format") != 1:
        raise ValueError("unsupported pack index format")
    return raw.get("packs", {})


def pack_release(organism: str, version: str | None = None) -> tuple[str, str, str]:
    """The address, checksum and version of a published pack, from the index in this build."""
    packs = pack_index()
    have = packs.get(organism)
    if not have:
        known = ", ".join(sorted(packs)) or "none"
        raise LookupError(f"no published pack for {organism!r} in this build; published: {known}")
    if version is None:
        version = max(have)      # versions are zero-padded dotted numbers, so this is the newest
    record = have.get(version)
    if not record:
        raise LookupError(f"no published pack {organism!r} {version!r}; have: {', '.join(sorted(have))}")
    path = record.get("path") or f"{organism}/{version}/{organism}-{version}.zip"
    return _checked_pack_url(f"{packs_host()}/{path}"), record["sha256"], version


def download_pack(organism: str, version: str | None = None) -> Path:
    """Download and install a published pack, with the checksum taken from this build.

    This is the route a user should take: the checksum comes from the application rather than from
    something they were told to paste, so a substituted file fails the check instead of installing.
    `$PARASCAPE_PACKS_URL` overrides the address for a mirror or a local test, and the checksum from
    the index still has to match."""
    url, sha, _version = pack_release(organism, version)
    override = os.environ.get(ENV_PACKS_URL)
    if override:
        url = _checked_pack_url(override)
    return download(url, sha, organism)


def download(url: str, expected_sha256: str, organism: str | None = None) -> Path:
    """Download an HTTPS pack with a trusted checksum and a bounded size, then install it."""
    if urlparse(url).scheme != "https" or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256):
        raise ValueError("downloads require HTTPS and an independently supplied SHA256")
    # no host allowlist here on purpose: this entry point takes an address AND a checksum the caller
    # supplies, which is the mirror case, and the checksum is what guards it. The allowlist applies
    # where ParaScape builds the address itself -- `pack_release` and the override in `download_pack`.
    with tempfile.TemporaryDirectory(prefix="parascape-download-") as work:
        archive = Path(work) / "pack.zip"
        with urlopen(url, timeout=30) as response, archive.open("wb") as fh:
            if urlparse(response.geturl()).scheme != "https":
                raise ValueError("pack download redirected away from HTTPS")
            total = 0
            while chunk := response.read(1024 * 1024):
                total += len(chunk)
                if total > MAX_ARCHIVE_BYTES:
                    raise ValueError("pack exceeds archive size limit")
                fh.write(chunk)
        return install(archive, expected_sha256, organism)


def main(argv=None):
    """Inspect, install or download a data pack without starting the desktop application."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("verify", "install"):
        command = commands.add_parser(name)
        command.add_argument("archive", type=Path)
        command.add_argument("--sha256")
        if name == "install":
            command.add_argument("--organism")
    command = commands.add_parser("download")
    command.add_argument("url", nargs="?",
                         help="a full https address; omit it and give --organism instead")
    command.add_argument("--sha256",
                         help="required with a URL; never needed with --organism, where the checksum "
                              "comes from the index shipped with this build")
    command.add_argument("--organism")
    command.add_argument("--version")
    commands.add_parser("published")
    command = commands.add_parser("installed")
    command.add_argument("organism")
    command.add_argument("--version")
    args = parser.parse_args(argv)
    if args.command == "verify":
        print(json.dumps(verify(args.archive, args.sha256), indent=2))
    elif args.command == "install":
        print(install(args.archive, args.sha256, args.organism))
    elif args.command == "download":
        if args.url:
            if not args.sha256:
                parser.error("--sha256 is required when a URL is given")
            print(download(args.url, args.sha256, args.organism))
        elif args.organism:
            print(download_pack(args.organism, args.version))
        else:
            parser.error("give a URL with --sha256, or --organism to use the published index")
    elif args.command == "published":
        print(json.dumps(pack_index(), indent=2, sort_keys=True))
    else:
        directory = installed(args.organism, args.version)
        print(directory or "Not installed")
        return 0 if directory else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
