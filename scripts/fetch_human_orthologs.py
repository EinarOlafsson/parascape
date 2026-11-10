"""Human orthologs of rat and cat genes, from the Ensembl Compara dump (instruction 53, WP4).

The rat space transfers annotation from human the way the mouse space does, which needs a gene-to-gene ortholog
map. BioMart refuses both GET and POST from here (405), and the REST endpoint answers one gene per call, so the
source is the Compara homology dump: one file per species listing its homologies with EVERY other species,
about 620 MB, streamed and split as it arrives.

The file read is the HUMAN one, because the per-species files are not symmetric: the rat file holds 4,058,043
rows over 180 species and no human rows at all, while the human file holds rat and cat. That was found by
filtering the rat file to nothing and then checking the file rather than assuming the filter was at fault.

Written in Python rather than a shell pipeline on purpose: the first two attempts were `curl | awk`, and a
pipeline reports the exit status of its LAST command, so a truncated download looked like success and left a
file holding nothing but its header. Here the stream is counted and the count is checked.

    python scripts/fetch_human_orthologs.py
"""
from __future__ import annotations

import csv
import io
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "host_gene_spaces_2026_11_09")
#: Read from the HUMAN side, not from each species'. The per-species files are not symmetric: the rat file
#: holds 4,058,043 rows over 180 species and **no human rows at all**, while the human file holds both rat and
#: cat. Found by filtering the rat file to nothing and checking instead of assuming the filter was wrong.
SOURCE = ("https://ftp.ensembl.org/pub/release-113/tsv/ensembl-compara/homologies/"
          "homo_sapiens/Compara.113.protein_default.homologies.tsv")
SPECIES = {"rattus_norvegicus": "rn_human_orthologs.tsv", "felis_catus": "fc_human_orthologs.tsv"}


def fetch() -> dict:
    """Stream the human homologies once and split out the rat and cat orthologs. Returns rows per species."""
    kept = {species: 0 for species in SPECIES}
    handles, writers = {}, {}
    read = 0
    with urllib.request.urlopen(SOURCE, timeout=120) as response:
        stream = io.TextIOWrapper(response, encoding="utf-8", errors="replace")
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        column = {name: i for i, name in enumerate(header)}
        needed = ("gene_stable_id", "homology_type", "homology_gene_stable_id", "homology_species",
                  "is_high_confidence")
        missing = [n for n in needed if n not in column]
        if missing:
            raise SystemExit(f"the dump has no column {missing}; header is {header}")
        for species, name in SPECIES.items():
            handles[species] = open(os.path.join(OUT, name), "w", newline="")
            writers[species] = csv.writer(handles[species], delimiter="\t")
            # The gene of THIS space first, then its human partner: a reader of the rat table wants rat ids.
            writers[species].writerow(["gene_id", "homology_type", "human_gene_id", "is_high_confidence"])
        try:
            for row in reader:
                read += 1
                if len(row) <= column["is_high_confidence"]:
                    continue
                species = row[column["homology_species"]]
                if species not in SPECIES or "ortholog" not in row[column["homology_type"]]:
                    continue
                writers[species].writerow([row[column["homology_gene_stable_id"]],
                                           row[column["homology_type"]],
                                           row[column["gene_stable_id"]],
                                           row[column["is_high_confidence"]]])
                kept[species] += 1
        finally:
            for handle in handles.values():
                handle.close()
    print(f"{read:,} human homology rows read; kept " +
          ", ".join(f"{n:,} for {s}" for s, n in kept.items()))
    if read < 1_000_000:
        raise SystemExit(f"only {read} rows arrived -- the stream was cut, not filtered empty")
    return kept


def main() -> int:
    """Fetch both species, refusing to leave a truncated file behind."""
    os.makedirs(OUT, exist_ok=True)
    kept = fetch()
    empty = [species for species, n in kept.items() if not n]
    if empty:
        raise SystemExit(f"no human orthologs found for {empty}, which cannot be right")
    return 0


if __name__ == "__main__":
    sys.exit(main())
