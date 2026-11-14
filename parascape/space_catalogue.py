"""Which organism spaces exist here, and what state each is in (instruction 64.23).

One row per declared space -- the two parasite maps shipped with the program and the human and mouse gene
packs -- with an explicit state, never implied:

* ``ready``        the space opens here (a parasite map whose cache is present);
* ``installed``    a host gene pack is installed (its version and hash are shown);
* ``installable``  a pack archive with its published SHA256 is available locally, not yet installed;
* ``unavailable``  neither; the row says what is missing.

``calibrated`` is True only for a space with a held-out track record; an installed host pack is never shown as
calibrated or as ready for inference until it has one. Nothing here touches the network.
"""
from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass

from . import organisms

#: Host spaces distributed as packs. Rn and Fc joined on 2026-11-09 with their gene builds only: they are
#: listed as unavailable until a pack is built, which is the honest state -- a declared space with no
#: measurement can name genes and answer nothing.
HOST_PACK_SPACES = ("Hs", "Mm", "Rn", "Fc")

#: Spaces declared with a gene build but no pack yet: the mosquito vector and the second malaria parasite
#: (instruction 53, WP5 and WP6). Listed so a reader knows they exist; `unavailable` until a pack is built.
DECLARED_ONLY_SPACES = ("Ag", "Pv", "As")

#: Every space distributed as a pack, whichever list declares it. The state of each is read from disk the same
#: way -- installed, a verified archive here, or nothing -- so a newly built pack shows up without a code
#: change. Before 2026-11-09 the second group was hardcoded `unavailable` and three real packs stayed hidden.
PACK_SPACES = HOST_PACK_SPACES + DECLARED_ONLY_SPACES
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass(frozen=True)
class Entry:
    """One space and its state."""

    code: str
    species: str
    kind: str
    distribution: str
    state: str
    version: str = ""
    archive: str = ""
    sha256: str = ""
    calibrated: bool = False
    note: str = ""


def _host_space(code: str):
    if code == "Hs":
        from .spaces.human import HUMAN_SPACE
        return HUMAN_SPACE
    if code == "Rn":
        from .spaces.rat import RAT_SPACE
        return RAT_SPACE
    if code == "Fc":
        from .spaces.cat import CAT_SPACE
        return CAT_SPACE
    if code == "Ag":
        from .spaces.anopheles import ANOPHELES_SPACE
        return ANOPHELES_SPACE
    if code == "Pv":
        from .spaces.vivax import VIVAX_SPACE
        return VIVAX_SPACE
    if code == "As":
        from .spaces.stephensi import STEPHENSI_SPACE
        return STEPHENSI_SPACE
    from .spaces.mouse import MOUSE_SPACE
    return MOUSE_SPACE


def _version_key(path: str) -> tuple:
    """The pack version in a file name, as numbers, so 2026.11.13 sorts after 2026.11.07."""
    match = re.search(r"parascape_[A-Za-z]+_([0-9.]+)\.zip$", os.path.basename(path))
    return tuple(int(p) for p in match.group(1).split(".")) if match else ()


def local_archives(code: str, roots=None) -> list:
    """(archive, sha256) for every pack of `code` with a published SHA256 file beside it, newest first.

    Newest by the VERSION in the file name, not by the path. Sorting on the path put the
    `mouse_gene_space_2026_11_07` pack ahead of `host_gene_packs_2026_11_13` because an `m` sorts
    after an `h`, so a newly built pack would have been installed only by accident of its folder name.
    """
    roots = roots or [os.path.join(_REPO, "results")]
    found = []
    for root in roots:
        for path in glob.glob(os.path.join(root, "*", f"parascape_{code}_*.zip")):
            sums = os.path.join(os.path.dirname(path), "SHA256")
            if not os.path.exists(sums):
                continue
            for line in open(sums):
                parts = line.split()
                if len(parts) == 2 and parts[1] == os.path.basename(path) and re.fullmatch(r"[0-9a-f]{64}", parts[0]):
                    found.append((path, parts[0]))
    return [p for _k, p in sorted(((_version_key(path), (path, sha)) for path, sha in found),
                                  key=lambda row: row[0], reverse=True)]


def _calibrated(code: str) -> bool:
    try:
        from . import track_record
        return bool(len(track_record.shipped(code)))
    except Exception:
        return False


def catalogue(roots=None) -> list:
    """Every declared space with its state (see the module docstring)."""
    from . import packs
    out = []
    for code, space in organisms.SPACES.items():
        try:
            present = os.path.exists(organisms.nodes_path(code)) and os.path.exists(organisms.graph_path(code))
        except FileNotFoundError:
            present = False
        out.append(Entry(code, space.species, space.kind, space.distribution,
                         "ready" if present else "unavailable", calibrated=present and _calibrated(code),
                         note="" if present else "its map cache is not built here"))
    for code in PACK_SPACES:
        space = _host_space(code)
        archives = local_archives(code, roots)
        try:
            where = packs.installed(code)
        except ValueError as error:
            where, broken = None, str(error)
        else:
            broken = ""
        if where is not None:
            manifest = packs._manifest((where / packs.MANIFEST).read_text())
            out.append(Entry(code, space.species, space.kind, "pack", "installed", version=manifest["version"],
                             note="gene space for links from parasite genes; no held-out track record yet, so "
                                  "nothing here is calibrated for inference"))
        elif archives:
            path, sha = archives[0]
            out.append(Entry(code, space.species, space.kind, "pack", "installable", archive=path, sha256=sha,
                             note="a verified pack is available in this checkout"))
        else:
            out.append(Entry(code, space.species, space.kind, "pack", "unavailable",
                             note=broken or "no pack installed and none available here"))
    return out


def install(code: str, archive: str = "", sha256: str = "", roots=None):
    """Install a host gene pack, checked against its published SHA256; returns the installed folder.

    Without arguments the newest local archive is used. A pack without a published checksum is refused.
    """
    from . import packs
    if code not in PACK_SPACES:
        raise ValueError(f"{code} is not a pack space")
    if not archive:
        archives = local_archives(code, roots)
        if not archives:
            raise FileNotFoundError(f"no {code} pack is available here")
        archive, sha256 = archives[0]
    if not sha256:
        sums = os.path.join(os.path.dirname(archive), "SHA256")
        sha256 = next((l.split()[0] for l in open(sums) if l.split()[1:] == [os.path.basename(archive)]), "") \
            if os.path.exists(sums) else ""
    if not re.fullmatch(r"[0-9a-f]{64}", sha256 or ""):
        raise ValueError("a pack needs its published SHA256 to be installed")
    return packs.install(archive, sha256, code)


def describe(entry: Entry) -> str:
    """One plain line for a menu or a table."""
    state = {"ready": "ready", "installed": f"installed ({entry.version})", "installable": "can be installed",
             "unavailable": "not available"}[entry.state]
    calibrated = "; calibrated on held-out genes" if entry.calibrated else "; not calibrated"
    return f"{entry.species} -- {state}{calibrated}" + (f". {entry.note}" if entry.note else "")
