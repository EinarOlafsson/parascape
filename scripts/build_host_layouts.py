"""Instruction 53, WP9: display coordinates for the host spaces, as two selectable lenses.

A map is a layout, and no host layout existed: `hs_host_proteins.parquet` and
`mm_host_proteins.parquet` carry no coordinate column, and the only layouts shipped are the
Toxoplasma and Plasmodium ones.

**Choosing what makes two host genes similar IS the claim the map makes**, so it is not a default to
be inherited from the parasite side. Einar chose **both, as selectable lenses** (2026-11-13), which
is the honest answer to a question with two defensible answers: build both and let the reader say
which claim they want the map to make.

* **`tissue`** -- two genes are similar when they are expressed in the same tissues and cell types.
  Measured, broad, and says nothing in itself about the parasite.
* **`parasite_facing`** -- two genes are similar in how they respond to infection and in what the
  parasite does to them: enrichment at the vacuole, the transcriptional response of a fibroblast or
  a hepatocyte to infection, rhoptry discharge, the erythrocyte surface. Closer to what ParaScape is
  for, and much sparser.

A gene that no lens measures is NOT dropped: every host gene keeps a position in both layouts, so a
map cannot quietly lose the genes with no data. Dropping them would make a map of what has been
studied look like a map of the proteome.

Each layout is tied to the table it was built from by the SHA-256 of the node file, the way
`pf_graph.npz` is, so the two cannot drift apart.

    python scripts/build_host_layouts.py            # both lenses for Hs and Mm
    python scripts/build_host_layouts.py --space Hs
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import paths  # noqa: E402
from parascape.embedding import EmbeddingSpec, embed as project  # noqa: E402

#: What each lens is allowed to look at. Named explicitly rather than discovered, because the point
#: of the choice is that it is a claim: a reader must be able to see exactly which measurements
#: produced the arrangement they are reading.
LENSES = {
    "Hs": {
        "tissue": ("fibroblast_tpm", "hepatocyte_tpm", "monocyte_tpm", "myeloid_dc_tpm",
                   "neuron_excitatory_ncpm", "neuron_inhibitory_ncpm", "astrocyte_ncpm",
                   "microglia_ncpm", "astrocyte_protein_ppm", "brain_protein_ppm",
                   "liver_protein_ppm"),
        "parasite_facing": ("pv_enrichment_log2", "hff_tg_infection_log2fc",
                            "hepatocyte_pf_infection_log2fc", "rhoptry_discharge_score",
                            "rhoptry_discharge_beta", "rbc_membrane_psms", "rbc_cytoplasm_psms",
                            "rbc_surface_copies_uk", "rbc_surface_copies_senegal"),
    },
    "Mm": {
        "tissue": ("brain_tpm", "skeletal_muscle_tpm", "bmdm_tpm", "mouse_brain_protein_ppm",
                   "bmdm_protein_ppm", "skeletal_muscle_protein_ppm"),
        # Thin on purpose, and the thinness is the finding: the mouse table carries one
        # parasite-facing measurement, the macrophage surface. A lens built on one column is
        # reported rather than padded with expression columns that would make it the other lens.
        "parasite_facing": ("bmdm_surface_intensity",),
    },
}

TABLE = {"Hs": "hs_host_proteins.parquet", "Mm": "mm_host_proteins.parquet"}

#: A lens must actually measure a reasonable share of the genes, or it is refused rather than drawn.
#:
#: This is not caution, it is a correction. The mouse parasite-facing lens has ONE column and 147 of
#: 16,473 genes carry a value -- 0.9 per cent. Built anyway, it produced 16,473 DISTINCT positions,
#: because UMAP spreads points apart even when their inputs are identical: 99 per cent of that map
#: would have been the repulsion term acting on median-imputed copies of the same row. It looks like
#: a map, reads like a map, and is noise. A map that is mostly imputation is worse than no map,
#: because nobody can see which it is.
MINIMUM_MEASURED = 0.10


def build(code: str, lens: str, columns: tuple, out_dir: str) -> dict:
    node_path = os.path.join(paths.data_dir(), TABLE[code])
    nodes = pd.read_parquet(node_path)
    present = [c for c in columns if c in nodes.columns]
    missing = [c for c in columns if c not in nodes.columns]
    frame = nodes[["host_id"] + present].rename(columns={"host_id": "gene_id"})

    measured = int(frame[present].notna().any(axis=1).sum())
    share = measured / len(frame) if len(frame) else 0.0
    if share < MINIMUM_MEASURED:
        print(f"  {code}/{lens}: REFUSED -- {measured:,} of {len(frame):,} genes measured "
              f"({100 * share:.1f} per cent), below the {100 * MINIMUM_MEASURED:.0f} per cent floor")
        print(f"    columns: {present or 'none'}. A layout here would be imputation wearing the "
              f"shape of a map.")
        return dict(space=code, lens=lens, refused=True, genes=len(frame),
                    genes_with_any_value=measured, features=list(present),
                    reason=f"only {100 * share:.1f} per cent of genes carry any value")
    spec = EmbeddingSpec(name=f"{code}_{lens}", blocks=(), extra_columns=tuple(present),
                         na_policy="median", max_missing=1.0, method="umap", n_components=3,
                         random_state=42)
    coords, features, rows, metadata = project(frame, spec, log=lambda m: print("   ", m),
                                               return_metadata=True)
    if not np.all(rows):
        raise ValueError(f"{code}/{lens}: a host layout must keep every gene")
    xyz = np.asarray(coords, dtype=np.float32)

    record = dict(metadata)
    record.update(space=code, lens=lens, features=list(features), requested_columns=list(columns),
                  missing_columns=missing, genes=len(frame), genes_with_any_value=measured,
                  node_file=TABLE[code],
                  node_file_sha256=hashlib.sha256(open(node_path, "rb").read()).hexdigest(),
                  coordinates_sha256=hashlib.sha256(xyz.tobytes()).hexdigest())

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{code.lower()}_host_graph_{lens}.npz")
    np.savez_compressed(path, xyz=xyz, gene_ids=frame.gene_id.to_numpy(dtype=str),
                        layout_metadata=np.array(json.dumps(record)))
    print(f"  {code}/{lens}: {len(frame):,} genes, {len(features)} features, "
          f"{measured:,} with at least one value ({100 * measured / len(frame):.1f} per cent)")
    if missing:
        print(f"    columns not in the table: {missing}")
    return record


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--space", choices=sorted(LENSES), help="only this space")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "host_layouts_2026_11_13"))
    args = ap.parse_args()

    summary = {}
    for code in ([args.space] if args.space else sorted(LENSES)):
        for lens, columns in LENSES[code].items():
            summary[f"{code}_{lens}"] = build(code, lens, columns, args.out)
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump({k: {m: v[m] for m in ("space", "lens", "features", "genes",
                                         "genes_with_any_value", "executed_method",
                                         "node_file_sha256", "coordinates_sha256")
                       if m in v} | ({"refused": True, "reason": v["reason"]}
                                     if v.get("refused") else {})
                   for k, v in summary.items()}, fh, indent=1)
    built = sum(1 for v in summary.values() if not v.get("refused"))
    print(f"\nwrote {built} layouts of {len(summary)} lenses to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
