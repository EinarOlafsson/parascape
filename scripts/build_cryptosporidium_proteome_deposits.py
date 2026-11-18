"""Two Cryptosporidium proteome deposits, keyed on CryptoDB `cgd` identifiers (instruction 41, 2026-11-18).

`Cp_protein abundance · other life stages` and `Cp_secretome / excreted` had no source. Both candidates
were already on disk in the Cryptosporidium acquisition store (MaxQuant `proteinGroups.txt`, extracted from
the depositors' own archives, checksummed), so nothing is fetched here.

* **PXD023935** (John et al. 2023, PMID 36775128): the total proteome of *C. parvum* oocysts, eight FAIMS
  fractions, MaxQuant iBAQ. Searched against the Iowa-ATCC and the Iowa II proteomes TOGETHER, so a protein
  group names the same protein in both annotations, `CPATCC_nnnnnnn-RA-p1;A3FPL8`.
* **PXD028969** (Li et al. 2025, PMID 41137792): proteins of *C. parvum* sporozoites that were released into
  the medium while the sporozoites glided (3 h) and adhered or invaded (6 h) on HCT-8 cells, TMT, searched
  against UniProt entries of *Cryptosporidium* (taxon 5807, every strain), so a group lists UniProt accessions.

**Identifier families, and how each reaches `cgd`.** `cgd` is the CryptoDB name of the Iowa II genome and the
key of the gene table. A UniProt accession of the Iowa II reference proteome (UP000006726) carries its `cgd`
name in the field "Gene Names (ORF)": an exact record, no map. A `CPATCC_` tag belongs to the Iowa-ATCC
annotation and reaches `cgd` only through `cpatcc_to_cgd.tsv`, the 2,571-gene exact-protein map. A group is
assigned to a gene only when every token that CAN be mapped agrees on ONE `cgd`; a group whose tokens name
two genes is dropped, a group none of whose tokens maps stays unmapped, and both are counted in the log.
No identifier is mapped by proximity or by name.

Writes `proteome_measured_cryptosporidium_oocyst_PXD023935.csv` and
`secretome_measured_cryptosporidium_sporozoite_PXD028969.csv` to the dataset root.
"""
from __future__ import annotations

import hashlib
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.environ.get("PARASCAPE_DATA") or "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"
STORE = os.path.join(SRC, "spaces", "Cp")
UNIPROT = os.path.join(STORE, "reference", "UniProt_UP000006726", "uniprot_UP000006726.tsv.gz")
CPATCC_MAP = os.path.join(ROOT, "results", "cryptosporidium_space_2026_11_08", "cpatcc_to_cgd.tsv")
OOCYST_IN = os.path.join(STORE, "abundance", "PXD023935", "proteinGroups.txt.gz")
SECRETOME_IN = os.path.join(STORE, "abundance", "PXD028969", "proteinGroups.txt.gz")
OOCYST_OUT = os.path.join(SRC, "proteome_measured_cryptosporidium_oocyst_PXD023935.csv")
SECRETOME_OUT = os.path.join(SRC, "secretome_measured_cryptosporidium_sporozoite_PXD028969.csv")


def _sha(path: str) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


def _clean(groups: pd.DataFrame) -> pd.DataFrame:
    """Drop reverse hits, contaminants and groups identified only by a modified site."""
    keep = (groups["Reverse"].isna() & groups["Potential contaminant"].isna()
            & groups["Only identified by site"].isna())
    return groups[keep].copy()


def _assign(tokens: str, uniprot: dict, cpatcc: dict) -> tuple:
    """(cgd or None, route) for one protein group; None with a reason when it cannot be one gene."""
    via_u, via_c = set(), set()
    for token in str(tokens).split(";"):
        if token in uniprot:
            via_u.add(uniprot[token])
        locus = token.split("-RA-p")[0]
        if locus.startswith("CPATCC_") and locus in cpatcc:
            via_c.add(cpatcc[locus])
    genes = via_u | via_c
    if not genes:
        return None, "unmapped"
    if len(genes) > 1:
        return None, "ambiguous"
    return next(iter(genes)), ("both agree" if via_u and via_c else "UniProt" if via_u else "CPATCC map")


def main() -> int:
    u = pd.read_csv(UNIPROT, sep="\t")
    uniprot = dict(zip(u["Entry"], u["Gene Names (ORF)"].fillna("")))
    uniprot = {k: v for k, v in uniprot.items() if v.startswith("cgd")}
    cpatcc_table = pd.read_csv(CPATCC_MAP, sep="\t")
    cpatcc = dict(zip(cpatcc_table["cpatcc_locus_tag"], cpatcc_table["cgd_gene_id"]))
    print(f"UniProt Iowa II entries with a cgd name: {len(uniprot):,}; CPATCC map: {len(cpatcc):,} genes")

    # ---- the oocyst proteome ----
    g = _clean(pd.read_csv(OOCYST_IN, sep="\t", low_memory=False))
    print(f"\nPXD023935 (sha {_sha(OOCYST_IN)}): {len(g):,} clean protein groups")
    routes = g["Protein IDs"].map(lambda t: _assign(t, uniprot, cpatcc))
    g["gene_id"] = [r[0] for r in routes]
    g["route"] = [r[1] for r in routes]
    print(g.route.value_counts().to_string())
    g = g[g["iBAQ"] > 0]
    total = g["iBAQ"].sum()          # the denominator keeps the groups that do not map: a share of what was seen
    ok = g[g.gene_id.notna()].sort_values("iBAQ", ascending=False)
    dup = int(ok.gene_id.duplicated().sum())
    ok = ok.drop_duplicates("gene_id")
    out = pd.DataFrame({"gene_id": ok.gene_id.to_numpy(),
                        "oocyst_protein_ppm": (ok["iBAQ"] / total * 1e6).round(4).to_numpy(),
                        "oocyst_unique_peptides": ok["Unique peptides"].to_numpy(),
                        "oocyst_route": ok.route.to_numpy()})
    out.to_csv(OOCYST_OUT, index=False)
    print(f"  {len(g):,} groups with iBAQ; {len(out):,} genes written ({dup} genes had a second group, the "
          f"larger kept); mapped share of iBAQ mass {100 * ok['iBAQ'].sum() / total:.1f}%")

    # ---- the sporozoite secretome ----
    s = _clean(pd.read_csv(SECRETOME_IN, sep="\t", low_memory=False))
    print(f"\nPXD028969 (sha {_sha(SECRETOME_IN)}): {len(s):,} clean protein groups")
    routes = s["Protein IDs"].map(lambda t: _assign(t, uniprot, {}))
    s["gene_id"] = [r[0] for r in routes]
    s["route"] = [r[1] for r in routes]
    print(s.route.value_counts().to_string())
    counts = [c for c in s.columns if c.startswith("Reporter intensity count")]
    s["channels"] = (s[counts] > 0).sum(axis=1)
    quantified = s[s.channels > 0]
    ok = quantified[quantified.gene_id.notna()].sort_values("channels", ascending=False)
    dup = int(ok.gene_id.duplicated().sum())
    ok = ok.drop_duplicates("gene_id")
    out = pd.DataFrame({"gene_id": ok.gene_id.to_numpy(), "secreted_sporozoite_detected": 1.0,
                        "secreted_channels_quantified": ok.channels.to_numpy()})
    out.to_csv(SECRETOME_OUT, index=False)
    print(f"  {len(quantified):,} groups quantified in at least one of {len(counts)} TMT channels; "
          f"{len(out):,} genes written ({dup} had a second group)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
