"""The mouse gene space and its distributable pack (instruction 64.22), as an executed notebook.

One row per mouse gene (Ensembl id) from the mouse files already on disk: Bgee expression (per-tissue scores ->
best score, best tissue, tissues present), the human ortholog and its type, and subcellular locations TRANSFERRED
from the human ortholog (Human Protein Atlas; inferred, never observed). The shipped mouse protein table
(`mm_host_proteins.parquet`, keyed by UniProt) is linked through the source gene->protein column in a separate
table, with every loss counted. Then a 3-D layout, the attribution notice, a hash-checked pack, and an
independent opening in an empty cache -- the same steps as the human pack (64.21).

    python scripts/build_mouse_gene_space_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
VERSION = "2026.11.07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"mouse_gene_space_{DATE}")
NB = os.path.join(REPO, "notebooks", f"mouse_gene_space_{DATE}.ipynb")
SRC = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The mouse gene space and pack", namespace={"REPO": REPO, "OUT": OUT, "SRC": SRC,
                                                                     "VERSION": VERSION})
    nb.md("**What this builds (instruction 64.22).** One row per mouse gene (Ensembl id) with what the mouse "
          "sources on disk say: where it is expressed (Bgee), its human ortholog, and the subcellular location of "
          "that ortholog. Then the same pack as for human: layout, notice, licences, opened on its own.",
          "**Two rules kept.** A protein is not a gene: the shipped mouse proteins are linked in their own table. "
          "And a location copied from a human ortholog is not a mouse observation: those columns are named "
          "`*_from_human_ortholog` and graded inferred, so they are never used as observed labels.",
          "**Provenance stated.** The source files were assembled on 2026-07-17 (Bgee, Ensembl Compara, Human "
          "Protein Atlas); the script that assembled them is not in the tree, so their hashes are recorded and "
          "their release versions are a stated gap.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import hashlib, json, os, shutil, numpy as np, pandas as pd",
            "from parascape import organisms as O, packs, paths",
            "from parascape.spaces import build_pack, validate_graph, validate_nodes",
            "from parascape.spaces.mouse import MOUSE_SPACE",
            "def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()",
            "files = {'transcriptome': os.path.join(SRC, 'hostloc_transcriptome_mouse.csv'),",
            "         'proteome': os.path.join(SRC, 'hostloc_proteome_mouse.csv'),",
            "         'shipped_proteins': os.path.join(REPO, 'parascape/data/mm_host_proteins.parquet')}",
            "hashes = {k: sha(p) for k, p in files.items()}",
            "tx, pr = pd.read_csv(files['transcriptome']), pd.read_csv(files['proteome'])",
            "assert len(tx) == len(pr) and (tx.ensembl.fillna('').values == pr.ensembl.fillna('').values).all()",
            "print(len(tx), 'source rows;', int(tx.ensembl.isna().sum()), 'without an Ensembl id')")
    nb.md("## 1 · One row per gene",
          "Bgee lists tissues with a score each (0-100, highest first); kept as the best score, the best tissue and "
          "the number of tissues where the gene is present. Rows without an Ensembl id are listed, not guessed.")
    nb.code("no_id = tx[tx.ensembl.isna()]",
            "no_id[['gene_name', 'uniprot']].to_csv(os.path.join(OUT, 'rows_without_gene_id.tsv'), sep='\\t', index=False)",
            "keep = tx.ensembl.notna()",
            "t, p = tx[keep], pr[keep]",
            "scores = t.Bgee_expression_score.fillna('').astype(str).str.split(';')",
            "genes = pd.DataFrame({'gene_id': t.ensembl.values, 'gene_name': t.gene_name.values, 'uniprot': t.uniprot.values,",
            "    'bgee_top_score': pd.to_numeric(scores.str[0], errors='coerce').values,",
            "    'bgee_top_tissue': t.tissue.fillna('').astype(str).str.split(';').str[0].replace('', np.nan).values,",
            "    'bgee_tissues_present': t.Bgee_n_tissues_present.values,",
            "    'location_from_human_ortholog': p.HPAortho_cell_localization.values,",
            "    'location_from_human_ortholog_probability': p.HPAortho_cell_localization_probability.values,",
            "    'human_ortholog': p.ortholog_human_gene.values, 'human_ortholog_type': p.ortholog_type.values,",
            "    'human_ortholog_confidence': p.ortholog_confidence.values})",
            "genes = genes.sort_values('gene_id').reset_index(drop=True)",
            "validate_nodes(MOUSE_SPACE, genes)",
            "print(len(genes), 'genes;', int(genes.bgee_top_score.notna().sum()), 'with Bgee expression;',",
            "      int(genes.location_from_human_ortholog.notna().sum()), 'with a location transferred from a human ortholog')")
    nb.md("## 2 · The shipped mouse proteins, linked to genes",
          "Through the source gene->UniProt column (one accession per gene in this source, not a reviewed "
          "many-to-many cross-reference; stated). Proteins with no gene are counted and listed.")
    nb.code("prot = pd.read_parquet(files['shipped_proteins'])",
            "acc = genes[['gene_id', 'uniprot']].dropna()",
            "links = prot[['host_id']].merge(acc, left_on='host_id', right_on='uniprot', how='left').drop(columns='uniprot')",
            "n = links.groupby('host_id').gene_id.nunique()",
            "links['link'] = links.host_id.map(n).map(lambda k: 'no gene' if k == 0 else ('one gene' if k == 1 else 'several genes'))",
            "by_link = links.drop_duplicates('host_id').link.value_counts().to_dict()",
            "links.to_parquet(os.path.join(OUT, 'mm_protein_gene_links.parquet'), index=False)",
            "links[links.link == 'no gene'][['host_id']].to_csv(os.path.join(OUT, 'proteins_without_gene.tsv'), sep='\\t', index=False)",
            "nodes = genes.drop(columns='uniprot').assign(shipped_protein_count=genes.gene_id.map(links.dropna().groupby('gene_id').host_id.nunique()).fillna(0).astype(int))",
            "by_link")
    nb.md("## 3 · What each column is, and its licence",
          "Bgee is CC0; Ensembl data carry no restriction on use (open); the Human Protein Atlas is CC BY-SA 4.0, "
          "so the transferred location columns are share-alike (SA). Everything else is open with attribution (O).")
    nb.code("base = dict(organism='Mm', storage_unit='gene')",
            "meta = {",
            " 'gene_id': dict(base, quantity_unit='identifier', evidence_grade='curation', source='Ensembl (via source file)', license='Ensembl, no restrictions', targets='identity only'),",
            " 'gene_name': dict(base, quantity_unit='symbol', evidence_grade='curation', source='Ensembl (via source file)', license='Ensembl, no restrictions', targets='identity only'),",
            " 'bgee_top_score': dict(base, quantity_unit='Bgee expression score 0-100 (best tissue)', evidence_grade='measured (summarised)', source='Bgee', license='CC0', targets='input or numeric label'),",
            " 'bgee_top_tissue': dict(base, quantity_unit='Uberon tissue name', evidence_grade='measured (summarised)', source='Bgee', license='CC0', targets='input or label'),",
            " 'bgee_tissues_present': dict(base, quantity_unit='count of tissues', evidence_grade='measured (summarised)', source='Bgee', license='CC0', targets='input or numeric label'),",
            " 'location_from_human_ortholog': dict(base, quantity_unit='HPA location terms', evidence_grade='inferred (orthology transfer)', source='Human Protein Atlas via Ensembl Compara', license='CC BY-SA 4.0', targets='input only; never an observed label'),",
            " 'location_from_human_ortholog_probability': dict(base, quantity_unit='probability 0-1', evidence_grade='inferred (orthology transfer)', source='Human Protein Atlas via Ensembl Compara', license='CC BY-SA 4.0', targets='input only'),",
            " 'human_ortholog': dict(base, quantity_unit='human gene symbol', evidence_grade='computed', source='Ensembl Compara', license='Ensembl, no restrictions', targets='bridge to the human space'),",
            " 'human_ortholog_type': dict(base, quantity_unit='Compara homology type', evidence_grade='computed', source='Ensembl Compara', license='Ensembl, no restrictions', targets='bridge'),",
            " 'human_ortholog_confidence': dict(base, quantity_unit='0/1', evidence_grade='computed', source='Ensembl Compara', license='Ensembl, no restrictions', targets='bridge'),",
            " 'shipped_protein_count': dict(base, quantity_unit='count', evidence_grade='derived', source='parascape mm_host_proteins + source UniProt column', license='derived', targets='identity only')}",
            "for c, m in meta.items():",
            "    m['source_sha256'] = hashes['proteome' if 'ortholog' in c else 'transcriptome']",
            "assert set(meta) == set(nodes.columns)",
            "json.dump(meta, open(os.path.join(OUT, 'column_metadata.json'), 'w'), indent=1)",
            "lic = {c: ('SA' if 'CC BY-SA' in meta[c]['license'] else 'O') for c in nodes.columns}",
            "pd.DataFrame(meta).T[['evidence_grade', 'license', 'targets']]")
    nb.md("## 4 · Layout, notice and pack",
          "Three coordinates per gene (UMAP, cosine, seed 0) from Bgee score, tissues present, best-tissue and "
          "transferred-location indicators; display only.")
    nb.code("import umap",
            "parts = [nodes[['bgee_top_score', 'bgee_tissues_present']].fillna(0).to_numpy(float) / [100.0, 600.0]]",
            "for c in ('bgee_top_tissue', 'location_from_human_ortholog'):",
            "    tok = nodes[c].fillna('').astype(str).str.split(';').explode().str.strip()",
            "    parts.append(pd.get_dummies(tok).groupby(level=0).max().drop(columns='', errors='ignore').to_numpy(float))",
            "X = np.hstack(parts) + 1e-6 * np.random.default_rng(0).standard_normal((len(nodes), sum(q.shape[1] for q in parts)))",
            "xyz = umap.UMAP(n_components=3, n_neighbors=15, metric='cosine', random_state=0, low_memory=True).fit_transform(X).astype(np.float32)",
            "build = os.path.join(OUT, 'build'); os.makedirs(build, exist_ok=True)",
            "nodes.to_parquet(os.path.join(build, MOUSE_SPACE.nodes), index=False)",
            "np.savez(os.path.join(build, MOUSE_SPACE.graph), xyz=xyz, gene_ids=nodes.gene_id.to_numpy(dtype=str))",
            "validate_graph(nodes, os.path.join(build, MOUSE_SPACE.graph))",
            "notice = '\\n'.join([f'ParaScape mouse gene pack {VERSION}', '',",
            "    'Expression: Bgee (https://www.bgee.org), CC0. Gene ids, symbols and human orthologs: Ensembl / Ensembl Compara (https://www.ensembl.org), no restrictions on use.',",
            "    'Locations (location_from_human_ortholog*): Human Protein Atlas (https://www.proteinatlas.org), CC BY-SA 4.0, transferred from the human ortholog;',",
            "    '  inferred, not observed in mouse. Those columns and adaptations of them must be shared under the same licence.',",
            "    'Source files assembled 2026-07-17; release versions of Bgee and Ensembl not recorded in them.', '',",
            "    'Licence per column: ' + ', '.join(f'{c}={v}' for c, v in lic.items())]) + '\\n'",
            "open(os.path.join(build, 'NOTICE.txt'), 'w').write(notice)",
            "archive = os.path.join(OUT, f'parascape_Mm_{VERSION}.zip')",
            "sha_pack = build_pack(MOUSE_SPACE, build, VERSION, lic, archive, extra={'NOTICE.txt': os.path.join(build, 'NOTICE.txt'),",
            "    'column_metadata.json': os.path.join(OUT, 'column_metadata.json'), 'mm_protein_gene_links.parquet': os.path.join(OUT, 'mm_protein_gene_links.parquet')})",
            "open(os.path.join(OUT, 'SHA256'), 'w').write(f'{sha_pack}  {os.path.basename(archive)}\\n')",
            "print(os.path.getsize(archive) // 1024, 'KiB;', sha_pack)")
    nb.md("## 5 · Opening it on its own, and checks")
    nb.code("state = os.path.join(OUT, '_open_check'); shutil.rmtree(state, ignore_errors=True)",
            "os.environ[paths.ENV_STATE] = state",
            "added = 'Mm' not in O.SPACES",
            "if added: O.SPACES['Mm'] = MOUSE_SPACE",
            "try:",
            "    where = packs.install(archive, sha_pack, 'Mm')",
            "    back = pd.read_parquet(O.nodes_path('Mm'))",
            "    validate_graph(back, O.graph_path('Mm'))",
            "    opened = dict(genes=len(back), identical=bool(back.equals(nodes)), notice=(where / 'NOTICE.txt').exists())",
            "finally:",
            "    if added: O.SPACES.pop('Mm')",
            "    os.environ.pop(paths.ENV_STATE)",
            "    shutil.rmtree(state, ignore_errors=True)",
            "checks = dict(unique_ids=bool(nodes.gene_id.is_unique), every_source_gene_kept=len(nodes) == int(tx.ensembl.notna().sum()),",
            "              sorted_order=bool(nodes.gene_id.is_monotonic_increasing), no_protein_ids_as_genes=bool(~nodes.gene_id.isin(prot.host_id).any()),",
            "              every_shipped_protein_in_links=bool(set(prot.host_id) == set(links.host_id)), opened_identical=opened['identical'],",
            "              notice_in_pack=bool(opened['notice']))",
            "assert all(checks.values()), checks",
            "report = dict(source_rows=len(tx), rows_without_gene_id=len(no_id), genes=len(nodes), shipped_proteins=int(prot.host_id.nunique()),",
            "              proteins_by_link=by_link, licences=lic, version=VERSION, pack_sha256=sha_pack, checks=checks, source_sha256=hashes,",
            "              gaps=['source release versions not recorded', 'gene->protein from one source column, not reviewed cross-references',",
            "                    'no observed mouse labels yet: locations are transferred from human; benchmarks belong to 64.29',",
            "                    'Mm not registered in organisms.SPACES until host spaces are wired in (64.24)'])",
            "json.dump(report, open(os.path.join(OUT, 'report.json'), 'w'), indent=1)",
            "shutil.rmtree(build)",
            "report")
    nb.write(NB)


if __name__ == "__main__":
    main()
