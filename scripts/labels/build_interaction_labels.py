"""Toxoplasma labels and partner tables from STRING and CORUM, with measured and copied evidence apart.

Einar, 2026-11-25: add STRING and CORUM to the Toxoplasma labels; keep evidence measured in
Toxoplasma apart from evidence copied from other organisms, and say which in every label name.

Sources (all public, downloaded into `--raw`):

- **STRING v12.0, Toxoplasma gondii VEG (taxon 432359)**, `protein.links.full` (every evidence
  channel split into "measured in this organism" and "transferred" from orthologs) and
  `protein.aliases` (UniProt -> `TGVEG_` locus tag). VEG and ME49 share locus numbers, so
  `TGVEG_123456` is joined to `TGME49_123456`; a number with no ME49 gene is dropped.
- **CORUM (current release, human complexes)**: subunit UniProt ids -> Ensembl genes (HGNC
  complete set) -> OrthoMCL 6.21 orthogroups (`hsap|ENSG...` members) -> the ParaScape genes in
  the same orthogroup (`nodes.orthogroup`). Orthogroups with more than `MAX_TOXO_PER_OG`
  Toxoplasma genes are not used (a kinase family would put dozens of genes in one complex).

What is built (each name says the evidence):

| Label | Evidence | Made from |
|---|---|---|
| interaction module (STRING, measured in Toxoplasma) | observed | STRING `experiments` or `database` (organism's own, not transferred) > 0; Louvain communities |
| interaction module (STRING, copied by orthology) | copied | `experiments_transferred` or `database_transferred` >= `COPIED_MIN`; Louvain communities |
| complex (CORUM human, copied by orthology) | copied | each gene's CORUM complex with the most mapped Toxoplasma genes |

STRING co-expression and text mining are NOT used: co-expression is circular with the transcript
view, and text mining is neither a measurement nor curated. Classes need `MIN_CLASS` genes.
Pair tables (`pairs_*.parquet`, one row per unordered gene pair) are written for the binding-partner
benchmark. `labels.parquet` holds one column per label; `sources.json` holds file checksums and counts.

    python scripts/labels/build_interaction_labels.py --raw DIR --out DIR [--nodes nodes.parquet]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms.community import louvain_communities

MIN_CLASS = 10
COPIED_MIN = 700
MAX_TOXO_PER_OG = 3
RESOLUTION = 1.0
SEED = 0
STRING_TAXON = "432359"
FILES = {"string_links": f"{STRING_TAXON}.protein.links.full.v12.0.txt.gz",
         "string_aliases": f"{STRING_TAXON}.protein.aliases.v12.0.txt.gz",
         "corum": "corum_human.txt",
         "hgnc": "hgnc_complete_set.txt",
         "orthomcl_groups": "groups_OrthoMCL-6.21.txt.gz"}
URLS = {"string_links": f"https://stringdb-downloads.org/download/protein.links.full.v12.0/{FILES['string_links']}",
        "string_aliases": f"https://stringdb-downloads.org/download/protein.aliases.v12.0/{FILES['string_aliases']}",
        "corum": "https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=human&file_format=txt",
        "hgnc": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        "orthomcl_groups": "https://orthomcl.org/common/downloads/release-6.21/groups_OrthoMCL-6.21.txt.gz"}
LABELS = {"string_measured": "interaction module (STRING, measured in Toxoplasma)",
          "string_copied": "interaction module (STRING, copied by orthology)",
          "corum_copied": "complex (CORUM human, copied by orthology)"}


def sha256(path: str) -> str:
    """Checksum of a source file, so a rebuild can say whether its inputs changed."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def string_pairs(raw: str, genes: set) -> pd.DataFrame:
    """STRING VEG links joined to ME49 genes: one row per unordered pair, every channel kept."""
    alias = pd.read_csv(os.path.join(raw, FILES["string_aliases"]), sep="\t", comment=None,
                        names=["protein", "alias", "source"], skiprows=1, dtype=str)
    alias = alias[alias["alias"].str.fullmatch(r"TGVEG_\d+")]
    alias["gene"] = alias["alias"].str.replace("TGVEG_", "TGME49_", regex=False)
    alias = alias[alias["gene"].isin(genes)].drop_duplicates("protein")
    to_gene = dict(zip(alias["protein"], alias["gene"]))
    t = pd.read_csv(os.path.join(raw, FILES["string_links"]), sep=" ")
    t["a"], t["b"] = t["protein1"].map(to_gene), t["protein2"].map(to_gene)
    t = t.dropna(subset=["a", "b"])
    t = t[t["a"] < t["b"]]                                     # each pair is listed both ways
    return t.drop(columns=["protein1", "protein2"]).rename(columns={"a": "gene_a", "b": "gene_b"})


def modules(pairs: pd.DataFrame, name: str, product: dict) -> pd.Series:
    """Louvain communities of a pair network; classes under MIN_CLASS genes are dropped."""
    g = nx.Graph()
    g.add_edges_from(zip(pairs["gene_a"], pairs["gene_b"]))
    comms = sorted(louvain_communities(g, resolution=RESOLUTION, seed=SEED), key=len, reverse=True)
    out = {}
    k = 0
    for members in comms:
        if len(members) < MIN_CLASS:
            continue
        k += 1
        hub = max(members, key=g.degree)
        tag = str(product.get(hub, hub))[:40]
        for m in members:
            out[m] = f"{name} {k:02d} (hub {hub}: {tag})"
    return pd.Series(out, dtype=object)


def corum(raw: str, nodes: pd.DataFrame) -> tuple:
    """CORUM human complexes mapped to Toxoplasma genes through OrthoMCL orthogroups."""
    hg = pd.read_csv(os.path.join(raw, FILES["hgnc"]), sep="\t", dtype=str,
                     usecols=["ensembl_gene_id", "uniprot_ids"]).dropna()
    uni2ens = {u: e for e, us in zip(hg["ensembl_gene_id"], hg["uniprot_ids"]) for u in us.split("|")}
    ens2og = {}
    with gzip.open(os.path.join(raw, FILES["orthomcl_groups"]), "rt") as f:
        for line in f:
            og, members = line.split(":", 1)
            for m in members.split():
                if m.startswith("hsap|"):
                    ens2og[m[5:].split(".")[0]] = og
    og2tg = nodes.groupby(nodes["orthogroup"].astype(str))["gene_id"].apply(list).to_dict()
    c = pd.read_csv(os.path.join(raw, FILES["corum"]), sep="\t", dtype=str, encoding="utf-8")
    rows, pair_rows = [], []
    for cid, cname, subs in zip(c["complex_id"], c["complex_name"], c["subunits_uniprot_id"].fillna("")):
        tg, used_og = set(), set()
        for u in subs.split(";"):
            og = ens2og.get(uni2ens.get(u.strip(), ""), "")
            members = og2tg.get(og, [])
            if og and 0 < len(members) <= MAX_TOXO_PER_OG:
                tg.update(members); used_og.add(og)
        if len(used_og) >= 2:
            rows.append({"complex_id": cid, "complex_name": cname, "toxo_genes": sorted(tg),
                         "human_subunits": len([u for u in subs.split(";") if u.strip()]),
                         "orthogroups_mapped": len(used_og)})
            s = sorted(tg)
            pair_rows += [{"gene_a": a, "gene_b": b, "corum_complex": cname, "corum_id": cid}
                          for i, a in enumerate(s) for b in s[i + 1:]]
    cx = pd.DataFrame(rows)
    # Each gene takes the MOST SPECIFIC complex (fewest mapped genes) that still has MIN_CLASS genes.
    # "Largest wins" (first build, 2026-11-25) put 104 ribosomal and nucleolar proteins under
    # "Nop56p-associated pre-rRNA complex" and left "Ribosome, cytoplasmic" with 19.
    pool = [r for r in rows if len(r["toxo_genes"]) >= MIN_CLASS]
    for _ in range(5):                                         # dropping a class can free genes
        best = {}
        for r in sorted(pool, key=lambda r: len(r["toxo_genes"])):
            for g_ in r["toxo_genes"]:
                best.setdefault(g_, r["complex_name"])
        lab = pd.Series(best, dtype=object)
        counts = lab.value_counts()
        small = set(counts[counts < MIN_CLASS].index)
        if not small:
            break
        pool = [r for r in pool if r["complex_name"] not in small]
    lab = lab[lab.map(lab.value_counts()) >= MIN_CLASS]
    pairs = pd.DataFrame(pair_rows).drop_duplicates(["gene_a", "gene_b"]) if pair_rows else pd.DataFrame()
    return lab, cx, pairs


def main(argv=None) -> int:
    """Build the three labels, the pair tables and the source record."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", required=True, help="folder with the downloaded source files")
    ap.add_argument("--out", required=True)
    ap.add_argument("--nodes", default="", help="nodes.parquet (default: the current release cache)")
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    missing = [f for f in FILES.values() if not os.path.exists(os.path.join(args.raw, f))]
    if missing:
        raise SystemExit(f"missing source files in {args.raw}: {missing}; URLs: {URLS}")
    if not args.nodes:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "infer_bench"))
        import e0_harness as H
        args.nodes = os.path.join(H.IM.CACHE, "nodes.parquet")
    nodes = pd.read_parquet(args.nodes, columns=["gene_id", "orthogroup", "product"])
    genes = set(nodes["gene_id"].astype(str))
    product = dict(zip(nodes["gene_id"], nodes["product"]))
    P = string_pairs(args.raw, genes)
    measured = P[(P["experiments"] > 0) | (P["database"] > 0)]
    copied = P[(P["experiments_transferred"] >= COPIED_MIN) | (P["database_transferred"] >= COPIED_MIN)]
    keep = ["gene_a", "gene_b", "experiments", "database", "experiments_transferred",
            "database_transferred", "combined_score"]
    measured[keep].to_parquet(os.path.join(args.out, "pairs_string_measured_in_toxoplasma.parquet"), index=False)
    copied[keep].to_parquet(os.path.join(args.out, "pairs_string_copied_by_orthology.parquet"), index=False)
    lab_c, cx, cpairs = corum(args.raw, nodes)
    cx.to_parquet(os.path.join(args.out, "corum_complexes_mapped.parquet"), index=False)
    cpairs.to_parquet(os.path.join(args.out, "pairs_corum_copied_by_orthology.parquet"), index=False)
    L = pd.DataFrame(index=nodes["gene_id"].astype(str))
    L[LABELS["string_measured"]] = modules(measured, "module", product).reindex(L.index)
    L[LABELS["string_copied"]] = modules(copied, "module", product).reindex(L.index)
    L[LABELS["corum_copied"]] = lab_c.reindex(L.index)
    L.index.name = "gene_id"
    L.reset_index().to_parquet(os.path.join(args.out, "labels.parquet"), index=False)
    record = {"files": {k: {"file": f, "url": URLS[k], "sha256": sha256(os.path.join(args.raw, f))}
                        for k, f in FILES.items()},
              "settings": {"MIN_CLASS": MIN_CLASS, "COPIED_MIN": COPIED_MIN, "MAX_TOXO_PER_OG": MAX_TOXO_PER_OG,
                           "louvain_resolution": RESOLUTION, "seed": SEED},
              "pairs": {"string_measured": int(len(measured)), "string_copied": int(len(copied)),
                        "corum_copied": int(len(cpairs))},
              "corum_complexes_mapped": int(len(cx)),
              "labels": {name: {"genes": int(L[name].notna().sum()), "classes": int(L[name].nunique())}
                         for name in L.columns}}
    json.dump(record, open(os.path.join(args.out, "sources.json"), "w"), indent=1)
    print(json.dumps({k: record[k] for k in ("pairs", "corum_complexes_mapped", "labels")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
