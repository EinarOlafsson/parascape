"""The human and mouse host PROTEIN spaces (instruction 53, WP9; Einar decided "both" on 2026-11-13).

ParaScape has two kinds of host space, kept apart on purpose:

* the host GENE spaces (`spaces.human.HUMAN_SPACE`, `spaces.mouse.MOUSE_SPACE`): Ensembl gene ids, distributed
  as packs, reached from a parasite gene through its links;
* the host PROTEIN spaces declared here: the tables that back the host bridges, `hs_host_proteins.parquet`
  (21,020 proteins) and `mm_host_proteins.parquet` (16,473), keyed on UniProt accessions and shipped in the data
  release.

A protein is not a gene. Mapping one onto the other inherits every ambiguity of a protein-to-gene mapping (one
gene, several reviewed proteins; one protein, several genes), so the protein tables are served as their own
spaces with their own identifier pattern rather than relabelled as the gene spaces.

Each protein space opens in a map window through one of its **lenses** -- the layouts built by
`scripts/build_host_layouts.py`, each a separate claim about what makes two host proteins similar:

* ``tissue``          expressed in the same tissues and cell types;
* ``parasite_facing`` alike in how they respond to infection and in what the parasite does to them.

A layout is served only beside the exact table it was built from: each file records the SHA-256 of its node
file, and a layout whose recorded table differs from the one on disk is refused rather than drawn over the
wrong proteins. The layouts are not bit-reproducible (umap-learn here does not honour its seed as a guarantee),
so the shipped file is the artefact.

These spaces are DECLARED, not registered in `organisms.SPACES`: registering a further space changes what
every slot, strategy and calibration caller sees, and a UniProt accession does not say which species it
belongs to, so `organisms.detect` could not tell the two apart.
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import pandas as pd

from .. import organisms as O
from . import validate_nodes

#: A UniProtKB accession (the documented pattern), with no isoform suffix: the tables hold one row per entry.
UNIPROT = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
RECORD = "https://www.uniprot.org/uniprotkb/{id}/entry"

HUMAN_PROTEIN_SPACE = O.Space(O.HUMAN_PROTEINS, "Homo sapiens host proteins", "UniProtKB human entries",
    O.HOST, UNIPROT, (), "UniProt", RECORD, O.HOST_TABLES[O.HUMAN], distribution="wheel")
MOUSE_PROTEIN_SPACE = O.Space(O.MOUSE_PROTEINS, "Mus musculus host proteins", "UniProtKB mouse entries",
    O.HOST, UNIPROT, (), "UniProt", RECORD, O.HOST_TABLES[O.MOUSE], distribution="wheel")

#: code -> space. Order is menu order.
SPACES = {s.code: s for s in (HUMAN_PROTEIN_SPACE, MOUSE_PROTEIN_SPACE)}

#: code -> lens -> layout file. Only lenses that were built; a refused lens is in `REFUSED` instead.
LENSES = {
    O.HUMAN_PROTEINS: {"tissue": "hs_host_graph_tissue.npz",
                       "parasite_facing": "hs_host_graph_parasite_facing.npz"},
    O.MOUSE_PROTEINS: {"tissue": "mm_host_graph_tissue.npz"},
}

#: Lenses refused by the layout builder, with the reason, so a menu can say why one is missing rather than
#: silently offering fewer. The mouse table carries one parasite-facing measurement (the macrophage surface),
#: and 147 of 16,473 proteins (0.9 per cent) have a value: below the builder's floor of a tenth.
REFUSED = {
    (O.MOUSE_PROTEINS, "parasite_facing"):
        "only 147 of 16,473 mouse proteins (0.9 per cent) carry a parasite-facing value; a map of that would "
        "be mostly filled-in values wearing the shape of a map",
}

#: Plain names for the lenses, as menus show them.
LENS_NAMES = {"tissue": "by tissue", "parasite_facing": "by what the parasite does to them"}

_PACKAGED = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def get(code: str) -> O.Space:
    """One protein space by code, or a KeyError naming the ones that exist."""
    if code not in SPACES:
        raise KeyError(f"no host protein space {code!r}; there are {', '.join(SPACES)}")
    return SPACES[code]


def nodes_path(code: str) -> str:
    """The protein table of a space, in the data directory (the downloaded release)."""
    from .. import paths
    return os.path.join(paths.data_dir(), get(code).nodes)


def layout_path(code: str, lens: str) -> str:
    """Where a lens layout lives: the data directory first, then the copy in the source tree.

    The second place exists because the layouts were added after the current data release was cut: until the
    next release ships them, the downloaded directory has the protein tables and not their layouts. Falling
    back is safe here, unlike for a table, because `load` refuses any layout whose recorded table is not the
    exact table being served -- a stale or foreign layout cannot be drawn over the wrong proteins.
    """
    from .. import paths
    name = LENSES[code][lens]
    for folder in (paths.data_dir(), _PACKAGED):
        path = os.path.join(folder, name)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(f"{get(code).species}: no {lens} layout ({name})")


def lenses(code: str) -> list:
    """The lenses of a space whose table and layout are both present here, in declared order."""
    if not os.path.exists(nodes_path(code)):
        return []
    out = []
    for lens in LENSES.get(code, {}):
        try:
            layout_path(code, lens)
        except FileNotFoundError:
            continue
        out.append(lens)
    return out


def available() -> list:
    """Codes of the protein spaces that can open here (a table and at least one lens)."""
    return [code for code in SPACES if lenses(code)]


def read_nodes(code: str) -> pd.DataFrame:
    """The protein table as a node table: `host_id` becomes `gene_id`, every value kept, ids validated.

    The rename is the only change. `gene_id` is the column every map, gate and table reads; the values stay
    UniProt accessions, and `validate_nodes` refuses the table if any is not one.
    """
    nodes = pd.read_parquet(nodes_path(code))
    if "host_id" not in nodes.columns or "gene_id" in nodes.columns:
        raise ValueError(f"{get(code).nodes}: expected a host_id column and no gene_id column")
    nodes = nodes.rename(columns={"host_id": "gene_id"})
    validate_nodes(get(code), nodes)
    return nodes


def layout_metadata(path: str) -> dict:
    """What a layout file records about itself: space, lens, features, and the table it was built from."""
    with np.load(path, allow_pickle=False) as z:
        return json.loads(str(z["layout_metadata"])) if "layout_metadata" in z.files else {}


def load(code: str, lens: str = ""):
    """(nodes, xyz, metadata) for one lens of one protein space; refuses a layout built from another table.

    Three checks, each of which has failed somewhere in this project before: the layout names the table it
    was built from by SHA-256, its proteins are the table rows in the same order, and every coordinate is
    finite.
    """
    lens = lens or next(iter(LENSES[code]))
    if (code, lens) in REFUSED:
        raise ValueError(f"{get(code).species}: the {lens} lens was refused: {REFUSED[(code, lens)]}")
    nodes = read_nodes(code)
    table_sha = hashlib.sha256(open(nodes_path(code), "rb").read()).hexdigest()
    path = layout_path(code, lens)
    meta = layout_metadata(path)
    if meta.get("node_file_sha256") != table_sha:
        raise ValueError(f"{os.path.basename(path)} was built from a different {get(code).nodes} "
                         f"(recorded {str(meta.get('node_file_sha256'))[:12]}, here {table_sha[:12]}); "
                         "rebuild it with scripts/build_host_layouts.py")
    with np.load(path, allow_pickle=False) as z:
        ids, xyz = z["gene_ids"].astype(str), z["xyz"].astype(np.float32)
    if not np.array_equal(ids, nodes.gene_id.astype(str).to_numpy()):
        raise ValueError(f"{os.path.basename(path)}: protein order differs from {get(code).nodes}")
    if xyz.shape != (len(nodes), 3) or not np.isfinite(xyz).all():
        raise ValueError(f"{os.path.basename(path)}: needs three finite coordinates per protein")
    return nodes, xyz, meta
