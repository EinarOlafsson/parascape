"""The human gene space, second partition (instruction 64.21), as an executed notebook.

Built on the validated foundation (`results/human_gene_space_foundation_2026_10_08/`: 58,988 GTEx v10 genes and
the reviewed UniProt cross-references). This partition adds the Human Protein Atlas gene annotations already on disk
(by Ensembl gene id) and links the protein table ParaScape ships (`hs_host_proteins.parquet`, keyed by UniProt) to
genes through the reviewed cross-references only. Protein measurements are NOT copied onto genes: the links are a
separate table, every one-to-many and lost join is counted. A distributable pack remains a separate step.

    python scripts/build_human_gene_space_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"human_gene_space_{DATE}")
NB = os.path.join(REPO, "notebooks", f"human_gene_space_{DATE}.ipynb")
FOUNDATION = os.path.join(REPO, "results", "human_gene_space_foundation_2026_10_08")
SRC = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The human gene space: annotations and protein links",
                          namespace={"REPO": REPO, "OUT": OUT, "SRC": SRC, "FOUNDATION": FOUNDATION})
    nb.md("**What this adds (instruction 64.21, second partition).** The foundation holds one row per human gene "
          "(58,988 GTEx v10 genes, Ensembl ids) and the reviewed UniProt cross-references. This partition adds the "
          "Human Protein Atlas annotations already on disk (RNA tissue specificity and distribution, subcellular "
          "location with its reliability, secretome location, protein class) and links the host-protein table "
          "ParaScape ships to genes.",
          "**Rule kept:** a protein is not a gene. Protein measurements stay on proteins; the gene-protein links are "
          "their own table, built only from reviewed cross-references, and every protein with no gene or several "
          "genes is counted.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import hashlib, os, json, numpy as np, pandas as pd",
            "def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()",
            "genes = pd.read_parquet(os.path.join(FOUNDATION, 'hs_gene_nodes.parquet'))",
            "xref = pd.read_parquet(os.path.join(FOUNDATION, 'gene_protein_crossrefs.parquet'))",
            "files = {'transcriptome': os.path.join(SRC, 'hostloc_transcriptome_human.csv'),",
            "         'proteome': os.path.join(SRC, 'hostloc_proteome_human.csv'),",
            "         'shipped_proteins': os.path.join(REPO, 'parascape/data/hs_host_proteins.parquet'),",
            "         'foundation_genes': os.path.join(FOUNDATION, 'hs_gene_nodes.parquet'),",
            "         'foundation_crossrefs': os.path.join(FOUNDATION, 'gene_protein_crossrefs.parquet')}",
            "hashes = {k: sha(p) for k, p in files.items()}",
            "tx, pr = pd.read_csv(files['transcriptome']), pd.read_csv(files['proteome'])",
            "print(len(genes), 'foundation genes;', len(xref), 'reviewed gene-protein pairs;', len(tx), 'HPA rows')")
    nb.md("## 1 · Human Protein Atlas annotations, by gene",
          "HPA rows without an Ensembl id cannot be placed and are listed (`hpa_rows_without_gene_id.tsv`); HPA genes "
          "outside the foundation are counted and listed, never added as new genes.")
    nb.code("no_id = tx[tx.ensembl.isna()]",
            "no_id.to_csv(os.path.join(OUT, 'hpa_rows_without_gene_id.tsv'), sep='\\t', index=False)",
            "keep_tx = ['ensembl', 'HPA_RNA_tissue_specificity', 'HPA_RNA_tissue_distribution', 'HPA_RNA_singlecell_specificity']",
            "keep_pr = ['ensembl', 'cell_localization', 'HPA_subcellular_reliability', 'HPA_secretome_location', 'HPA_protein_class']",
            "hpa = tx[keep_tx].dropna(subset=['ensembl']).merge(pr[keep_pr].dropna(subset=['ensembl']), on='ensembl', how='outer')",
            "dup_hpa = int(hpa.ensembl.duplicated().sum())",
            "hpa = hpa.drop_duplicates('ensembl').rename(columns={'ensembl': 'gene_id', 'cell_localization': 'HPA_subcellular_location'})",
            "outside = hpa[~hpa.gene_id.isin(genes.gene_id)]",
            "outside.to_csv(os.path.join(OUT, 'hpa_genes_outside_foundation.tsv'), sep='\\t', index=False)",
            "space = genes.merge(hpa, on='gene_id', how='left')",
            "annotated = int(space['HPA_RNA_tissue_specificity'].notna().sum())",
            "print(len(no_id), 'HPA rows without an id;', dup_hpa, 'duplicated;', len(outside), 'HPA genes outside the foundation;',",
            "      annotated, 'foundation genes now carry an HPA annotation')")
    nb.md("## 2 · The shipped host proteins, linked to genes",
          "Each shipped protein (UniProt accession) is looked up in the reviewed cross-references. A protein with one "
          "gene, several genes or none is counted; the links table keeps all of them, and the measurements stay on "
          "the protein.")
    nb.code("prot = pd.read_parquet(files['shipped_proteins'])",
            "links = prot[['host_id']].merge(xref[['gene_id', 'protein_id']], left_on='host_id', right_on='protein_id', how='left').drop(columns='protein_id')",
            "n_genes = links.groupby('host_id').gene_id.nunique()",
            "links['link'] = links.host_id.map(n_genes).map(lambda n: 'no reviewed gene' if n == 0 else ('one gene' if n == 1 else 'several genes'))",
            "genes_per_protein = links.drop_duplicates('host_id').link.value_counts().to_dict()",
            "proteins_per_gene = links.dropna(subset=['gene_id']).groupby('gene_id').host_id.nunique()",
            "links.to_parquet(os.path.join(OUT, 'hs_protein_gene_links.parquet'), index=False)",
            "links[links.link == 'no reviewed gene'][['host_id']].to_csv(os.path.join(OUT, 'proteins_without_gene.tsv'), sep='\\t', index=False)",
            "space['shipped_protein_count'] = space.gene_id.map(proteins_per_gene).fillna(0).astype(int)",
            "print(genes_per_protein, '; genes with more than one shipped protein:', int((proteins_per_gene > 1).sum()))")
    nb.md("## 3 · What each column is",
          "Species, unit, source, licence and what the column may be used for. HPA data are distributed under "
          "Creative Commons Attribution-ShareAlike 4.0 (attribution and share-alike apply to a redistributed pack).")
    nb.code("meta = json.load(open(os.path.join(FOUNDATION, 'column_metadata.json')))",
            "hpa_lic = 'Human Protein Atlas, CC BY-SA 4.0 (attribution; share-alike on redistribution)'",
            "for c in [c for c in hpa.columns if c != 'gene_id']:",
            "    src = 'transcriptome' if c in keep_tx else 'proteome'",
            "    meta[c] = dict(organism='Hs', storage_unit='gene', quantity_unit='category (HPA vocabulary)', evidence_grade='curation',",
            "                   source_id='hpa_' + src, source_sha256=hashes[src], context='Human Protein Atlas consensus',",
            "                   license=hpa_lic, redistribution='allowed with attribution and share-alike',",
            "                   missingness='gene not in HPA, or no value', targets='label or input; location columns are labels for host localisation questions')",
            "meta['shipped_protein_count'] = dict(organism='Hs', storage_unit='gene', quantity_unit='count', evidence_grade='derived',",
            "    source_id='parascape hs_host_proteins + reviewed UniProt crossrefs', source_sha256=hashes['shipped_proteins'],",
            "    license='derived', redistribution='as its sources', missingness='0 = no reviewed shipped protein', targets='identity only')",
            "json.dump(meta, open(os.path.join(OUT, 'column_metadata.json'), 'w'), indent=1)",
            "pd.DataFrame({c: {k: meta[c].get(k) for k in ('quantity_unit', 'evidence_grade', 'license')} for c in space.columns}).T")
    nb.md("## 4 · Checks", "Unique ids, every foundation gene kept in its order, no protein row turned into a gene.")
    nb.code("checks = dict(unique_ids=bool(space.gene_id.is_unique),",
            "              every_foundation_gene_kept=bool(len(space) == len(genes) and (space.gene_id.values == genes.gene_id.values).all()),",
            "              foundation_values_unchanged=bool(space[genes.columns].equals(genes)),",
            "              no_protein_ids_as_genes=bool(~space.gene_id.isin(prot.host_id).any()),",
            "              every_shipped_protein_in_links=bool(set(prot.host_id) == set(links.host_id)),",
            "              every_column_described=all(c in meta for c in space.columns))",
            "assert all(checks.values()), checks",
            "space.to_parquet(os.path.join(OUT, 'hs_genes.parquet'), index=False)",
            "report = dict(foundation_genes=len(genes), hpa_rows_without_gene_id=len(no_id), hpa_duplicated_ids=dup_hpa,",
            "              hpa_genes_outside_foundation=len(outside), genes_with_hpa_annotation=annotated,",
            "              shipped_proteins=int(prot.host_id.nunique()), proteins_by_link=genes_per_protein,",
            "              genes_with_several_shipped_proteins=int((proteins_per_gene > 1).sum()), checks=checks, sha256=hashes)",
            "json.dump(report, open(os.path.join(OUT, 'mapping_report.json'), 'w'), indent=1, default=str)",
            "report")
    nb.write(NB)


if __name__ == "__main__":
    main()
