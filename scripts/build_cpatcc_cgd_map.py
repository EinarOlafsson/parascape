"""Instruction 41: correspondence between the two Cryptosporidium parvum annotations.

GSE271211 -- the only candidate found for `Cp_transcription · oocyst` whose oocysts are measured
rather than named in a growth protocol -- reports counts on `CPATCC_` locus tags. The ParaScape Cp
gene table is keyed on CryptoDB `cgd` identifiers. They are two independent annotations of two
assemblies of the same organism:

  * `CPATCC_`: GCA_015245375.1 (ASM1524537v1), *C. parvum* IOWA-ATCC, annotated by the University
    of Georgia, released 2022-10-21. 4,035 genes, which matches the GSE271211 counts row for row.
  * `cgd`:     GCA_000165345.1 (ASM16534v1), *C. parvum* Iowa II, via Ensembl Protists. 3,923 genes
    in the shipped table.

There is no shared identifier and NCBI Gene holds no `CPATCC_` record, so no lookup exists to join
on. What does exist is the protein sequence, and these are two annotations of ONE organism, so most
proteins should be byte-identical.

**This maps only what is identical, and says plainly what it does not reach.** A pair is accepted
only when one ATCC gene and one cgd gene share a protein sequence exactly and neither sequence
appears anywhere else -- so a gene family with repeated sequences is excluded rather than resolved
by a guess. Everything else is left unmapped and counted, because two thirds of a mapping reported
honestly is worth more than a whole one reported hopefully.

Writes ``results/cryptosporidium_space_2026_11_08/cpatcc_to_cgd.tsv``.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import os
import re
import sys
import urllib.request
import zipfile
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

OUT = os.path.join(ROOT, "results", "cryptosporidium_space_2026_11_08")

ATCC = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        "GCA_015245375.1/download?include_annotation_type={t}")
ENSEMBL = ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/current/fasta/"
           "protists_alveolata1_collection/cryptosporidium_parvum_iowa_ii_gca_000165345/pep/"
           "Cryptosporidium_parvum_iowa_ii_gca_000165345.ASM16534v1.pep.all.fa.gz")


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=300) as fh:
        return fh.read()


def fasta(text, identify):
    """(identifier, sequence) per record, skipping records the identifier cannot name."""
    out, cur, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if cur:
                out.append((cur, "".join(seq)))
            cur, seq = identify(line), []
        elif cur:
            seq.append(line)
    if cur:
        out.append((cur, "".join(seq)))
    return [(g, s) for g, s in out if g and s]


def main() -> int:
    # The protein FASTA names the locus tag only for uncharacterized proteins -- 2,719 of 3,896
    # carry a product name instead -- so the tag comes from the GFF, where every CDS has one.
    gff = zipfile.ZipFile(io.BytesIO(fetch(ATCC.format(t="GENOME_GFF"))))
    gff_name = next(n for n in gff.namelist() if n.endswith(".gff"))
    protein_of_tag = {}
    for line in gff.read(gff_name).decode().splitlines():
        if "\tCDS\t" not in line:
            continue
        tag = re.search(r"locus_tag=([^;\n]+)", line)
        pid = re.search(r"protein_id=([^;\n]+)", line) or re.search(r"Name=([^;\n]+)", line)
        if tag and pid:
            protein_of_tag[pid.group(1)] = tag.group(1)

    prot = zipfile.ZipFile(io.BytesIO(fetch(ATCC.format(t="PROT_FASTA"))))
    faa = next(n for n in prot.namelist() if n.endswith(".faa"))
    atcc = fasta(prot.read(faa).decode(),
                 lambda h: protein_of_tag.get(h[1:].split()[0]))
    cgd = fasta(gzip.decompress(fetch(ENSEMBL)).decode(),
                lambda h: (re.search(r"gene:(\S+)", h).group(1) if "gene:" in h else None))

    def by_sequence(rows):
        index = defaultdict(set)
        for gene, seq in rows:
            index[hashlib.sha256(seq.replace("*", "").encode()).hexdigest()].add(gene)
        return index

    left, right = by_sequence(atcc), by_sequence(cgd)
    shared = set(left) & set(right)
    # One gene on each side, or it is not a correspondence. A hash holding several genes is a
    # repeated sequence -- a gene family -- and picking one of them would be invention.
    candidate = sorted((a, c) for h in shared
                       if len(left[h]) == 1 and len(right[h]) == 1
                       for a in left[h] for c in right[h])
    ambiguous = sum(1 for h in shared if len(left[h]) > 1 or len(right[h]) > 1)

    # Uniqueness per SEQUENCE is not uniqueness per GENE. A gene with two protein products appears
    # under two hashes and can pick up a partner from each: `CPATCC_0027530` matched `cgd2_3260`
    # and `cgd2_3270`, two adjacent genes, because the ATCC annotation merges a locus that CryptoDB
    # splits. Keeping either one would be choosing, so the gene is dropped and counted. Found by
    # the test, not by reading the code.
    seen_left = defaultdict(int)
    seen_right = defaultdict(int)
    for a, c in candidate:
        seen_left[a] += 1
        seen_right[c] += 1
    pairs = [(a, c) for a, c in candidate if seen_left[a] == 1 and seen_right[c] == 1]
    split = len(candidate) - len(pairs)

    atcc_genes = {g for g, _ in atcc}
    cgd_genes = {g for g, _ in cgd}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "cpatcc_to_cgd.tsv"), "w") as fh:
        fh.write("cpatcc_locus_tag\tcgd_gene_id\tevidence\n")
        for a, c in pairs:
            fh.write(f"{a}\t{c}\tidentical protein sequence, unique on both sides\n")

    print(f"ATCC genes (GCA_015245375.1): {len(atcc_genes)}")
    print(f"cgd  genes (GCA_000165345.1): {len(cgd_genes)}")
    print(f"one-to-one pairs            : {len(pairs)}")
    print(f"  ATCC mapped: {len({a for a, _ in pairs})} ({100*len({a for a,_ in pairs})/len(atcc_genes):.0f} per cent)")
    print(f"  cgd  mapped: {len({c for _, c in pairs})} ({100*len({c for _,c in pairs})/len(cgd_genes):.0f} per cent)")
    print(f"  sequences shared but repeated on a side (left unmapped): {ambiguous}")
    print(f"  dropped because one gene matched more than one partner: {split}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
