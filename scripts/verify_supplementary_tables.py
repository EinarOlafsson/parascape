"""Check that a supplementary table really holds the data we ship from it (instruction 65.04, batch 4).

Most sources do not record a deposit. They record a LOCATION -- "eLife 85654 supplementary file 6" -- and for
those `verify_origin_papers.py` can only say that the location is the paper itself, with nothing to cross-check.
Verifying them means opening the table and asking whether the genes we carry a value for are the genes that
table is about. That is the question this script answers, for every such source whose file is on this machine.

**What it checks, and what it deliberately does not.** It does NOT re-derive the column: each dataset has its
own derivation, and reimplementing 39 of them here would test this script rather than the data. It checks the
identity link, which is where this kind of source goes wrong: the gene identifiers in the supplementary table
are read out of every sheet and compared with the genes carrying a non-null value in the installed column. A
source whose table names genes we do not carry, or which carries genes the table never mentions, is reported
with the counts so a person can look. An overlap at or above `AGREES` is the table we say it is.

A low overlap is not automatically an error -- a derivation may filter hard, and a table keyed on protein
accessions rather than genes resolves through the identity layer -- so the verdict says "look at this", never
"this is wrong".

    STARPLAST_DATA=... python scripts/verify_supplementary_tables.py
"""
from __future__ import annotations

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")
VERIFICATION = os.path.join(OUT, "origin_paper_verification.tsv")

#: Where a recorded `datasets/...` path is rooted. `paths.dataset_root()` answers for the packaged caches, not
#: for the raw acquisition tree, which is why it is named here rather than derived.
BASES = (os.environ.get("STARPLAST_DATA_PARENT")
         or "/media/carruthers/mnt3/claude/toxoplasma_projects", ROOT)

#: Gene identifiers these sources are keyed on, plus UniProt accessions for the proteomics tables.
IDENTIFIER = re.compile(
    r"\b(TGME49_\d{6}[A-Z]?|TGGT1_\d{6}[A-Z]?|TGVEG_\d{6}[A-Z]?"
    r"|PF3D7_\d{7}|PVP01_(?:\d{7}|API\d{5}|MIT\d{5})|PBANKA_\d{6}"
    r"|ENS(?:MUS|RNO|FCA)?G\d{11}"
    r"|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")

#: At or above this share of the installed genes found in the table, the table is the one we say it is.
AGREES = 0.5
#: Sheets are capped so a 200 MB workbook cannot stall the pass; the identity link shows up in the first rows.
MAX_ROWS = 60000


#: Sources where a low overlap was read and explained, so it is not re-flagged every run. The identity check
#: compares gene accessions; it cannot apply when the table is keyed on something else, and saying that is more
#: useful than reporting 0 % forever. Each entry names the evidence inside the file.
KEYED_DIFFERENTLY = {
    "pf_hsp90_chemoproteome":
        "the table is keyed on a 'UniProt ID' column, not on PF3D7 gene accessions, so the link to our "
        "gene-keyed column runs through the identity layer and cannot be seen by matching accessions",
    "host_k562_rhoptry_screen":
        "the data sheet is MAGeCK output keyed on gene SYMBOLS while the host table is keyed on UniProt "
        "accessions; the first sheet is a README, which is why so few identifiers appear at all",
    "host_gtex_transcriptome":
        "the GTEx median-TPM matrix is keyed on ENSG gene identifiers while the human host table is keyed "
        "on UniProt accessions, so the two never match as strings; the link is the identity layer",
    "pf_resistome":
        "the table is one row per SNV in compound-selected clones, and two of the six installed columns "
        "(pf6_field_dnds, pf6_field_nonsyn_snvs) come from the Pf6 field data instead, so this table is "
        "expected to cover only part of the column set",
}
#: Sources whose unit is not a gene at all, so a gene-identity check is the wrong question.
NOT_GENE_KEYED = {
    "lipidome_vesicles": "the rows are lipid species in the metabolite table, not genes",
    "metabolome_iron": "the rows are metabolites, not genes",
    "pf_host_bridge_xlms": "it contributes bridge rows between two organisms, not a node column",
}


def resolve(path: str):
    """The recorded path on this machine, or None when the file was never kept here."""
    for base in BASES:
        candidate = os.path.join(base, path)
        if os.path.exists(candidate):
            return candidate
    return None


def identifiers_in(path: str) -> tuple:
    """Every gene or protein identifier the file mentions, and how it was read.

    Workbooks are read sheet by sheet with no header assumption, because a supplementary table puts its real
    header anywhere -- after a title row, after a legend, sometimes twice. Reading every cell as text and
    matching identifiers sidesteps that entirely: it does not matter which column they are in.
    """
    import pandas as pd
    lower = path.lower()
    try:
        if lower.endswith((".xlsx", ".xls")):
            book = pd.read_excel(path, sheet_name=None, header=None, nrows=MAX_ROWS, dtype=str)
            text = "\n".join(frame.astype(str).to_csv(index=False) for frame in book.values())
            return set(IDENTIFIER.findall(text)), f"{len(book)} sheet(s)"
        if lower.endswith((".csv", ".tsv", ".txt", ".gct", ".gz")):
            separator = "," if ".csv" in lower else "\t"
            frame = pd.read_csv(path, sep=separator, header=None, nrows=MAX_ROWS,
                                dtype=str, engine="python", on_bad_lines="skip")
            return set(IDENTIFIER.findall(frame.astype(str).to_csv(index=False))), "text table"
        if lower.endswith(".parquet"):
            frame = pd.read_parquet(path)
            return set(IDENTIFIER.findall(frame.head(MAX_ROWS).astype(str).to_csv(index=False))), "parquet"
        if lower.endswith(".pdf"):
            return set(), "PDF, not a table a program can read"
    except Exception as exc:                                  # a corrupt or exotic workbook is a finding
        return set(), f"could not be read: {type(exc).__name__}"
    return set(), "unsupported file type"


#: A Toxoplasma gene is the same gene in every strain, and a supplementary table may key on any of them:
#: TGGT1_224540 and TGME49_224540 are one gene under two accessions, which the identity layer resolves at
#: ingest. Comparing the strings literally reported 0 % overlap for four CRISPR screens whose tables are
#: entirely GT1 -- the identifiers matched perfectly once the strain was taken off.
def normalise(identifier: str) -> str:
    """One comparable key per gene, with the Toxoplasma strain removed."""
    match = re.fullmatch(r"TG(?:ME49|GT1|VEG)_(\d{6}[A-Z]?)", identifier)
    return f"TG_{match.group(1)}" if match else identifier


#: Which shipped table holds a given organism's genes. Picking the first table that merely CONTAINS the column
#: name compared a Plasmodium source against Toxoplasma genes, because both tables carry `melting_temperature`.
TABLE_FOR = {"Tg": "nodes.parquet", "Pf": "pf_nodes.parquet",
             "Hs": "hs_host_proteins.parquet", "Mm": "mm_host_proteins.parquet"}


def installed_genes(columns, organism: str = "") -> tuple:
    """The genes carrying a value in those columns, and which shipped table they came from."""
    import pandas as pd
    from parascape import paths
    data = paths.data_dir()
    preferred = TABLE_FOR.get(str(organism), "")
    order = ([preferred] if preferred else []) + [
        n for n in ("nodes.parquet", "pf_nodes.parquet", "hs_host_proteins.parquet",
                    "mm_host_proteins.parquet", "metabolites.parquet") if n != preferred]
    for name in order:
        path = os.path.join(data, name)
        if not os.path.exists(path):
            continue
        head = pd.read_parquet(path, columns=None).head(0)
        present = [c for c in columns if c in head.columns]
        if not present:
            continue
        key = "gene_id" if "gene_id" in head.columns else head.columns[0]
        frame = pd.read_parquet(path, columns=[key] + present)
        keep = frame[present].notna().any(axis=1)
        return set(frame.loc[keep, key].astype(str)), name
    return set(), ""


def main() -> int:
    """Verify every supplementary-table source whose file is here, and write the verdicts."""
    import pandas as pd
    from parascape import datasets as ds

    verification = pd.read_csv(VERIFICATION, sep="\t")
    wanted = sorted(set(verification.loc[
        verification.verdict == "supplementary table of a named paper", "dataset_id"]))
    registry = {d.key: d for d in ds.REGISTRY}
    rows = []
    for key in wanted:
        entry = registry[key]
        if key in NOT_GENE_KEYED:
            rows.append(dict(dataset_id=key, verdict="not keyed on genes, so this check does not apply",
                             detail=NOT_GENE_KEYED[key]))
            continue
        if key in KEYED_DIFFERENTLY:
            rows.append(dict(dataset_id=key, verdict="read: keyed on something other than gene accessions",
                             detail=KEYED_DIFFERENTLY[key]))
            continue
        recorded = str(entry.path or "")
        if not recorded:
            rows.append(dict(dataset_id=key, verdict="no file recorded",
                             detail="the source records a paper but no file, so there is nothing to open"))
            continue
        path = resolve(recorded)
        if path is None:
            rows.append(dict(dataset_id=key, verdict="file not on this machine",
                             detail=f"{recorded} is recorded but absent here"))
            continue
        if os.path.isdir(path):
            rows.append(dict(dataset_id=key, verdict="directory, not one table",
                             detail=f"{recorded} is a folder; the source needs a file named before it can be read"))
            continue
        mentioned, how = identifiers_in(path)
        ours, table = installed_genes(list(entry.columns or ()), getattr(entry, 'organism', ''))
        if not ours:
            rows.append(dict(dataset_id=key, verdict="no installed column to compare",
                             detail=f"{', '.join(entry.columns or ()) or 'no columns'} is not in a shipped table"))
            continue
        if not mentioned:
            rows.append(dict(dataset_id=key, verdict="table could not be read", detail=how))
            continue
        keys = {normalise(i) for i in mentioned}
        shared = {g for g in ours if normalise(str(g)) in keys}
        share = len(shared) / len(ours)
        verdict = ("the table names the genes we ship" if share >= AGREES
                   else "LOOK: the table and the column barely overlap")
        rows.append(dict(dataset_id=key, verdict=verdict,
                         detail=(f"{len(shared):,} of {len(ours):,} installed genes appear in the table "
                                 f"({share:.0%}); the table names {len(mentioned):,} identifiers [{how}, "
                                 f"column read from {table}]")))
        print(f"{key:<34}{share:>5.0%}  {verdict}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "supplementary_table_verification.tsv"), sep="\t", index=False)
    print()
    print(frame.verdict.value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
