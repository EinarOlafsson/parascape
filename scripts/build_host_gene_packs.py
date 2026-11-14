"""Rebuild the human and mouse gene packs, now carrying the lens layouts and the GEO columns.

Instruction 53, route 1. Two things changed on 2026-11-13 and both land here.

**The lens layouts.** Each pack carried one layout, `hs_gene_graph.npz`, whose recipe was not
recorded beside it. It now carries two, `<code>_gene_graph_tissue.npz` and
`<code>_gene_graph_parasite_facing.npz`, each recording its own features, its recipe and the SHA-256
of the node table it was built from. The space declares the tissue layout as its default.

**GEO-derived columns ship, with citation.** Einar decided on 2026-11-13: "yes, with citation". That
was not the one-line change the earlier note promised, because `packs.build` accepted only `O` and
`SA`. It now accepts a third code, `G` -- a public deposit -- and ONLY when both an accession and a
citation are recorded for that column, which is what the permission actually says. Every `G` column
is named in `NOTICE.txt` with its accession and the paper behind it.

**Publishing the packs is not this script.** A separate agent is moving all ParaScape data hosting to
Hugging Face on branch `data-on-huggingface` (Einar, 2026-11-13). This script writes the archives and
their SHA256 into `results/`, where `space_catalogue` finds them for a local install; pointing the
downloader at a published URL belongs to that work and is deliberately left to it.

    python scripts/build_host_gene_packs.py
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

SOURCE = os.path.join(ROOT, "results", "host_gene_lenses_2026_11_13")
PREVIOUS = {"Hs": os.path.join(ROOT, "results", "human_pack_2026_11_07",
                               "parascape_Hs_2026.11.07.zip"),
            "Mm": os.path.join(ROOT, "results", "mouse_gene_space_2026_11_07",
                               "parascape_Mm_2026.11.07.zip")}
OUT = os.path.join(ROOT, "results", "host_gene_packs_2026_11_13")
VERSION = "2026.11.13"
LINKS = {"Hs": "hs_protein_gene_links.parquet", "Mm": "mm_protein_gene_links.parquet"}

NOTICE_HEAD = {
    "Hs": """ParaScape human gene pack {version}

Genes and expression: GTEx Analysis Release V10 (dbGaP phs000424), median TPM, GENCODE 39 / GRCh38.p13.
  The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of the Office of the
  Director of the NIH. Data retrieved 2026-10-08 from https://gtexportal.org. Use does not imply
  endorsement by GTEx or the NIH. This pack pins GTEx V10.
Gene-protein links: UniProt Knowledgebase release 2026_03, reviewed human entries, CC BY 4.0.
Annotations and single-cell expression (HPA_ and hpa_sc_ columns): Human Protein Atlas
  (https://www.proteinatlas.org), CC BY-SA 4.0; those columns, and works adapted from them, must be
  shared under the same licence.
""",
    "Mm": """ParaScape mouse gene pack {version}

Genes: Ensembl GRCm39 release 116. Expression calls and per-tissue expression scores: Bgee 15.2 (CC0).
Human orthologs and the locations transferred along them: Ensembl Compara. A transferred location is
  inferred, never measured in the mouse, and its column name says so.
Gene-protein links: UniProt Knowledgebase release 2026_03, reviewed mouse entries, CC BY 4.0.
""",
}

GEO_NOTICE = """
Public deposits redistributed with citation (licence code G)
------------------------------------------------------------
Each column below is derived from a public sequence-archive deposit. A deposit states no licence of
its own, so it is redistributed here with its accession and the publication behind it named, which is
the condition under which it may travel in this pack. The values are our own summaries (a log2 ratio
between the arms of that deposit), not the depositor's files.
"""


def _licence_lines(licences: dict, sources: dict) -> str:
    out = []
    for column in sorted(licences):
        code = licences[column]
        if code != "G":
            continue
        record = sources[column]
        out.append(f"  {column}: {record['accession']} -- {record['citation']}")
    return "\n".join(out)


def main() -> int:
    """Build both packs, verify each round-trips, and write SHA256 beside them."""
    import pandas as pd

    from parascape import packs
    from parascape.space_catalogue import _host_space
    from parascape.spaces import build_pack
    from parascape.spaces import host_genes

    with open(os.path.join(SOURCE, "summary.json")) as fh:
        summary = json.load(fh)

    os.makedirs(OUT, exist_ok=True)
    sums, report = [], {}
    for code in ("Hs", "Mm"):
        space = _host_space(code)
        built = summary[code]
        licences = dict(built["column_licenses"])
        sources = dict(built["column_sources"])
        table = pd.read_parquet(os.path.join(SOURCE, space.nodes))
        undeclared = [c for c in table.columns if c not in licences]
        if undeclared:
            raise SystemExit(f"{code}: no licence declared for {undeclared}")

        staging = os.path.join(OUT, f"_stage_{code}")
        shutil.rmtree(staging, ignore_errors=True)
        os.makedirs(staging)
        table.to_parquet(os.path.join(staging, space.nodes), index=False)
        extra = {}
        for lens, name in host_genes.LENSES[code].items():
            shutil.copyfile(os.path.join(SOURCE, f"{code.lower()}_gene_graph_{lens}.npz"),
                            os.path.join(staging, name))
            if name != space.graph:
                extra[name] = os.path.join(staging, name)
        # The link table and the column descriptions travel with the pack, as they did before: the
        # link table is what carries a selection between the two kinds of host space.
        with zipfile.ZipFile(PREVIOUS[code]) as z:
            z.extract(LINKS[code], staging)
        extra[LINKS[code]] = os.path.join(staging, LINKS[code])
        metadata = os.path.join(SOURCE, f"{code.lower()}_column_metadata.json")
        shutil.copyfile(metadata, os.path.join(staging, "column_metadata.json"))
        extra["column_metadata.json"] = os.path.join(staging, "column_metadata.json")

        notice = NOTICE_HEAD[code].format(version=VERSION)
        cited = _licence_lines(licences, sources)
        if cited:
            notice += GEO_NOTICE + cited + "\n"
        notice += ("\nLicence per column: "
                   + ", ".join(f"{c}={licences[c]}" for c in sorted(licences)) + "\n")
        with open(os.path.join(staging, "NOTICE.txt"), "w") as fh:
            fh.write(notice)
        extra["NOTICE.txt"] = os.path.join(staging, "NOTICE.txt")

        destination = os.path.join(OUT, f"parascape_{code}_{VERSION}.zip")
        # The 2026.11.07 table is the baseline: `validate_nodes` refuses the build if any gene,
        # column or measured value of it is missing from the new one. The columns are added, never
        # swapped, and this is what says so rather than a claim in a commit message.
        baseline = os.path.join(SOURCE, "previous", space.nodes)
        sha = build_pack(space, staging, VERSION, licences, destination,
                         previous=baseline if os.path.exists(baseline) else None,
                         extra=extra, column_sources=sources)
        manifest = packs.verify(destination, sha)
        size = os.path.getsize(destination)
        sums.append(f"{sha}  {os.path.basename(destination)}")
        report[code] = dict(version=VERSION, sha256=sha, bytes=size, genes=len(table),
                            columns=len(table.columns),
                            cited_columns=sorted(c for c, v in licences.items() if v == "G"),
                            files=sorted(manifest["files"]))
        print(f"{code}: {len(table):,} genes, {len(table.columns)} columns, "
              f"{size / 1e6:.1f} MB, {len(report[code]['cited_columns'])} cited deposit columns")
        shutil.rmtree(staging, ignore_errors=True)

    with open(os.path.join(OUT, "SHA256"), "w") as fh:
        fh.write("\n".join(sums) + "\n")
    with open(os.path.join(OUT, "packs.json"), "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"\n{OUT}")
    print("Publishing these archives is the Hugging Face hosting work on branch "
          "data-on-huggingface, not this script.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
