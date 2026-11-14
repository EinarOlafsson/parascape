"""The lenses of the host GENE spaces (instruction 53, WP9 route 1).

Route 2 gave the host PROTEIN tables their own spaces and their own two layouts. This module does
the same for the host GENE spaces -- `spaces.human.HUMAN_SPACE` (58,988 ENSG) and
`spaces.mouse.MOUSE_SPACE` (18,428 ENSMUSG) -- so the two kinds of host space are read through the
same two claims:

* ``tissue``          expressed in the same tissues and cell types;
* ``parasite_facing`` alike in how they respond to infection and in what the parasite does to them.

Two differences from `host_proteins`, both of which follow from the gene spaces being **packs**
rather than wheel data:

* a lens layout is resolved through `paths.space_file`, out of the installed pack, never from the
  data directory. A gene space with no installed pack has no lenses, which is the honest answer;
* the 2026-11-07 packs shipped one layout, `hs_gene_graph.npz`, whose recipe was not recorded beside
  it. It is not carried forward. Each space's `graph` now names its tissue layout, so what a pack
  calls its default layout is a lens with its features written into it.

**The mouse parasite-facing lens is built here, and was refused for the protein table.** The mouse
protein table carries one parasite-facing measurement (147 of 16,473 proteins); the mouse gene table
carries ten, from six Toxoplasma, Cryptosporidium and Plasmodium series. Same floor, different data.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd

from . import validate_graph, validate_nodes

#: code -> lens -> the layout file inside that space's pack. Order is menu order, and the first is
#: the default lens: the one a space opens with when nobody picks.
LENSES = {
    "Hs": {"tissue": "hs_gene_graph_tissue.npz",
           "parasite_facing": "hs_gene_graph_parasite_facing.npz"},
    "Mm": {"tissue": "mm_gene_graph_tissue.npz",
           "parasite_facing": "mm_gene_graph_parasite_facing.npz"},
}

#: Lenses the builder refused, with the reason. Empty for the gene spaces: both lenses measure well
#: over the tenth of genes the builder requires. Kept rather than dropped, because the next host gene
#: space may well have a thin lens and a menu has to be able to say why one is missing.
REFUSED: dict = {}

#: Plain names for the lenses, as menus show them. Same wording as the protein spaces, deliberately:
#: the two kinds of space make the same two claims and should not read as different questions.
LENS_NAMES = {"tissue": "by tissue", "parasite_facing": "by what the parasite does to them"}


def space(code: str):
    """One host gene space by code, or a KeyError naming the ones that have lenses."""
    if code not in LENSES:
        raise KeyError(f"no host gene space {code!r}; there are {', '.join(LENSES)}")
    if code == "Hs":
        from .human import HUMAN_SPACE
        return HUMAN_SPACE
    from .mouse import MOUSE_SPACE
    return MOUSE_SPACE


def default_lens(code: str) -> str:
    """The lens a space opens with when nobody picks one."""
    return next(iter(LENSES[code]))


def layout_path(code: str, lens: str) -> str:
    """Where one lens layout lives inside the installed pack of `code`.

    Raises `FileNotFoundError` when the pack is not installed or does not carry that lens, which is
    what makes `lenses` able to answer honestly rather than offering a menu entry that fails.
    """
    from .. import paths
    if lens not in LENSES[code]:
        raise KeyError(f"{code}: no {lens} lens; there are {', '.join(LENSES[code])}")
    return paths.space_file(space(code), LENSES[code][lens])


def lenses(code: str) -> list:
    """The lenses of a space whose layout is in the installed pack, in declared order."""
    out = []
    for lens in LENSES.get(code, {}):
        try:
            layout_path(code, lens)
        except (FileNotFoundError, ValueError):
            continue
        out.append(lens)
    return out


def available() -> list:
    """Codes of the host gene spaces that can open here (an installed pack with a lens in it)."""
    return [code for code in LENSES if lenses(code)]


def layout_metadata(path: str) -> dict:
    """What a layout file records about itself: space, lens, features, and the table behind it."""
    with np.load(path, allow_pickle=False) as z:
        return json.loads(str(z["layout_metadata"])) if "layout_metadata" in z.files else {}


def read_nodes(code: str) -> pd.DataFrame:
    """The gene table of an installed pack, with its ids validated against the space pattern."""
    from .. import paths
    target = space(code)
    nodes = pd.read_parquet(paths.space_file(target, target.nodes))
    validate_nodes(target, nodes)
    return nodes


def load(code: str, lens: str = ""):
    """(nodes, xyz, metadata) for one lens of one host gene space.

    The same three checks the protein spaces make, for the same reason -- each has failed somewhere
    in this project before: the layout names the table it was built from by SHA-256, its genes are
    the table rows in the same order, and every coordinate is finite. The order and finiteness come
    from `validate_graph`; the SHA-256 tie is what stops a layout built from an older table being
    drawn over a newer one, which no order check would notice when the gene list did not change.
    """
    lens = lens or default_lens(code)
    if (code, lens) in REFUSED:
        raise ValueError(f"{space(code).species}: the {lens} lens was refused: "
                         f"{REFUSED[(code, lens)]}")
    from .. import paths
    target = space(code)
    nodes = read_nodes(code)
    node_file = paths.space_file(target, target.nodes)
    table_sha = hashlib.sha256(open(node_file, "rb").read()).hexdigest()
    path = layout_path(code, lens)
    meta = layout_metadata(path)
    recorded = meta.get("node_file_sha256")
    if recorded != table_sha:
        raise ValueError(f"{LENSES[code][lens]} was built from a different {target.nodes} "
                         f"(recorded {str(recorded)[:12]}, here {table_sha[:12]}); rebuild it with "
                         "scripts/build_host_gene_lenses.py")
    validate_graph(nodes, path)
    with np.load(path, allow_pickle=False) as z:
        xyz = z["xyz"].astype(np.float32)
    return nodes, xyz, meta
