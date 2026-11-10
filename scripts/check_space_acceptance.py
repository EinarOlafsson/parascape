"""WP3 acceptance for the five further spaces (instruction 53).

The acceptance for the human and mouse spaces was: of the reviewed UniProt proteins of that species carrying an
Ensembl cross-reference, at least 95 % must resolve to a gene in the space's universe. It asks the only
question that matters about a gene universe -- is anything real missing from it -- using a source that was not
used to build it.

**It does not transfer to every species, and that is a finding rather than a failure.** The measure needs a
reviewed proteome, and the less-studied species barely have one: rat has 8,239 reviewed proteins, cat 234 and
Anopheles gambiae 303, but Anopheles stephensi has **9** and Plasmodium vivax **11**. A 95 % bar over 9
proteins is a statement about nine proteins, not about a universe of 13,336 genes.

So two checks are run, and each space gets whichever is meaningful:

* **reviewed-proteome resolution**, where enough reviewed proteins exist (at least 100);
* **independent gene-list coverage**, where a list of that species' genes exists that was produced by someone
  else's pipeline -- the quantified genes in a published series -- which is what the parasite and vector
  spaces have instead.

    python scripts/check_space_acceptance.py
"""
from __future__ import annotations

import io
import json
import os
import sys
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SPACES = os.path.join(ROOT, "results", "host_gene_spaces_2026_11_09")
OUT = os.path.join(ROOT, "results", "space_acceptance_2026_11_09")
BAR = 0.95
MIN_REVIEWED = 100

#: space code -> (NCBI taxon id, the Ensembl id prefix its GENE ids use, the build that maps transcript to gene)
#: UniProt's Ensembl cross-reference gives TRANSCRIPT ids (ENSRNOT..., not ENSRNOG...), so the build is read
#: for the transcript-to-gene parentage. Reading the gene ids straight out of the cross-reference finds none,
#: which is what made rat look uncheckable at first.
TAXA = {"Rn": (10116, "ENSRNOG", "reference/ensembl_r113/Rattus_norvegicus.mRatBN7.2.113.gff3.gz"),
        "Fc": (9685, "ENSFCAG", "reference/ensembl_r113/Felis_catus.Felis_catus_9.0.113.gff3.gz"),
        "Ag": (7165, "AGAP", "reference/ensembl_metazoa_r60/Anopheles_gambiae.AgamP4.60.chr.gff3.gz"),
        "As": (30069, "ASTE", "reference/ensembl_metazoa_r60/Anopheles_stephensi.AsteS1.60.gff3.gz"),
        "Pv": (5855, "PVP01_",
               "reference/ensembl_protists_r63/Plasmodium_vivax_gca900093555.GCA900093555v2.63.gff3.gz")}
#: space code -> an independent gene list, produced by someone else's pipeline.
INDEPENDENT = {
    "Ag": ("transcription/scRNAseq/GSE284537/GSM8686106_AnG_r1-36h-GFP-NS_features.tsv.gz",
           "the genes quantified in GSE284537 by its authors' own Cell Ranger run"),
    "Pv": ("transcription/RNAseq/GSE220039/GSE220039_PvHsMmu_GRCh38_GRCm39_P01v58_raw_counts.csv.gz",
           "the genes counted in GSE220039 by its authors' own pipeline"),
}


def transcript_to_gene(relative: str) -> dict:
    """Transcript id -> gene id, from the build, so a UniProt cross-reference can be resolved to a gene."""
    import gzip
    import re

    from parascape import paths

    path = os.path.join(paths.dataset_root(), relative)
    if not os.path.exists(path):
        return {}
    mapping = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            field = line.rstrip("\n").split("\t")
            if len(field) < 9 or field[2] not in ("mRNA", "transcript", "ncRNA", "rRNA", "tRNA"):
                continue
            child = re.search(r"ID=transcript:([^;]+)", field[8])
            parent = re.search(r"Parent=gene:([^;]+)", field[8])
            if child and parent:
                mapping[child.group(1)] = parent.group(1)
    return mapping


def reviewed_proteins(taxon: int, prefix: str, mapping: dict) -> set:
    """Gene ids reached by that species' reviewed UniProt proteins, through the build's parentage."""
    query = urllib.parse.quote(f"reviewed:true AND organism_id:{taxon}")
    url = (f"https://rest.uniprot.org/uniprotkb/stream?query={query}"
           f"&fields=accession,xref_ensembl&format=tsv")
    with urllib.request.urlopen(url, timeout=180) as response:
        text = io.TextIOWrapper(response, encoding="utf-8").read()
    genes, rows = set(), 0
    for line in text.splitlines()[1:]:
        rows += 1
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        for token in parts[1].replace(";", " ").split():
            token = token.split(" ")[0].split(".")[0].strip("[]")
            if token.startswith(prefix):
                genes.add(token)                      # already a gene id (AGAP, ASTE, PVP01_)
            elif token in mapping:
                genes.add(mapping[token])             # a transcript id, resolved through the build
    return genes, rows


def independent_genes(code: str) -> tuple:
    """Gene ids from a list someone else's pipeline produced, or (None, reason)."""
    import gzip

    import pandas as pd
    from parascape import paths

    if code not in INDEPENDENT:
        return None, "no independent gene list on disk"
    relative, what = INDEPENDENT[code]
    path = os.path.join(paths.dataset_root(), relative)
    if not os.path.exists(path):
        return None, f"{relative} is not on disk"
    if path.endswith("features.tsv.gz"):
        with gzip.open(path, "rt") as handle:
            return {line.split("\t")[0] for line in handle.read().splitlines()}, what
    frame = pd.read_csv(path, usecols=[0])
    return set(frame.iloc[:, 0].astype(str).str.strip()), what


def main() -> int:
    """Check every further space and write the verdicts."""
    import pandas as pd

    from parascape import space_catalogue

    os.makedirs(OUT, exist_ok=True)
    rows = []
    for code, (taxon, prefix, build) in TAXA.items():
        space = space_catalogue._host_space(code)
        path = os.path.join(SPACES, space.nodes)
        if not os.path.exists(path):
            continue
        universe = set(pd.read_parquet(path, columns=["gene_id"]).gene_id.astype(str))
        cross, reviewed = reviewed_proteins(taxon, prefix, transcript_to_gene(build))
        row = dict(space=code, species=space.species, genes=len(universe), reviewed_proteins=reviewed,
                   reviewed_with_ensembl=len(cross))
        if reviewed >= MIN_REVIEWED and cross:
            resolved = len(cross & universe)
            row["resolve_rate"] = round(resolved / len(cross), 4)
            row["verdict"] = "passes" if resolved / len(cross) >= BAR else "FAILS the 95 % bar"
            row["check"] = "reviewed proteome"
            row["detail"] = f"{resolved} of {len(cross)} cross-referenced genes are in the universe"
        else:
            genes, what = independent_genes(code)
            if genes:
                shared = len({g for g in genes if g.startswith(prefix)} & universe)
                theirs = len({g for g in genes if g.startswith(prefix)})
                row["resolve_rate"] = round(shared / theirs, 4) if theirs else None
                row["verdict"] = ("passes" if theirs and shared / theirs >= BAR
                                  else "FAILS the 95 % bar")
                row["check"] = "independent gene list"
                row["detail"] = f"{shared} of {theirs} genes in {what} are in the universe"
            else:
                row["verdict"] = "cannot be checked"
                row["check"] = "none available"
                row["detail"] = (f"only {reviewed} reviewed proteins exist, and {what}"
                                 if reviewed < MIN_REVIEWED else what)
        rows.append(row)
        print(f"{code}: {row['verdict']} -- {row['detail']}")
    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "acceptance.tsv"), sep="\t", index=False)
    json.dump(dict(bar=BAR, min_reviewed=MIN_REVIEWED,
                   passing=sorted(frame.loc[frame.verdict == "passes", "space"]),
                   failing=sorted(frame.loc[frame.verdict.str.startswith("FAILS"), "space"]),
                   unchecked=sorted(frame.loc[frame.verdict == "cannot be checked", "space"])),
              open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
