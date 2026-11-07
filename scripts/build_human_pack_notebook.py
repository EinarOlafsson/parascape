"""The distributable human gene pack (instruction 64.21, last partition), as an executed notebook.

Takes the gene table from `build_human_gene_space_notebook.py`, lays the genes out in three dimensions from their
GTEx expression and Human Protein Atlas annotations, writes the attribution notice the source terms require,
builds a hash-checked pack with a licence per column, then installs it into a throwaway cache and opens it there
the way ParaScape opens any pack space.

    python scripts/build_human_pack_notebook.py
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
SPACE_DIR = os.path.join(REPO, "results", f"human_gene_space_{DATE}")
OUT = os.path.join(REPO, "results", f"human_pack_{DATE}")
NB = os.path.join(REPO, "notebooks", f"human_pack_{DATE}.ipynb")
STATE = os.environ.get("HUMAN_PACK_STATE", os.path.join(OUT, "_open_check"))


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The human gene pack", namespace={"REPO": REPO, "OUT": OUT, "SPACE_DIR": SPACE_DIR,
                                                           "VERSION": VERSION, "STATE": STATE})
    nb.md("**What this builds (instruction 64.21).** A data pack for the human gene space: the gene table (58,988 "
          "GTEx v10 genes with HPA annotations), a 3-D layout of the genes, the reviewed gene-protein links, the "
          "description of every column and the attribution notice. Each column carries its licence: GTEx and "
          "UniProt-derived columns are open with attribution (O); HPA columns are share-alike (SA).",
          "**Then it is opened on its own:** installed into an empty cache and read back through the same calls "
          "ParaScape uses for any pack space.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, shutil, numpy as np, pandas as pd",
            "from parascape import organisms as O, packs, paths",
            "from parascape.spaces import build_pack, validate_graph, validate_nodes",
            "from parascape.spaces.human import HUMAN_SPACE",
            "nodes = pd.read_parquet(os.path.join(SPACE_DIR, 'hs_genes.parquet'))",
            "meta = json.load(open(os.path.join(SPACE_DIR, 'column_metadata.json')))",
            "validate_nodes(HUMAN_SPACE, nodes)",
            "print(len(nodes), 'genes,', nodes.shape[1], 'columns')")
    nb.md("## 1 · Where each gene sits",
          "Three coordinates per gene from log expression in the two GTEx contexts and the HPA categories "
          "(tissue specificity and distribution, single-cell specificity, subcellular location, secretome location, "
          "protein class), one indicator per category value. UMAP, cosine distance, 15 neighbours, fixed seed. "
          "Genes without any HPA annotation are placed by expression alone. The layout is for display; it is not an "
          "input to any answer.")
    nb.code("import umap",
            "cats = [c for c in nodes.columns if c.startswith('HPA_')]",
            "parts = [np.log1p(nodes[['gtex_fibroblast_median_tpm', 'gtex_lcm_hepatocyte_median_tpm']].fillna(0)).to_numpy()]",
            "for c in cats:",
            "    tokens = nodes[c].fillna('').astype(str).str.replace(';', ',').str.split(',')",
            "    parts.append(pd.get_dummies(tokens.explode().str.strip()).groupby(level=0).max().drop(columns='', errors='ignore').to_numpy(float))",
            "X = np.hstack(parts)",
            "X = X + 1e-6 * np.random.default_rng(0).standard_normal(X.shape)  # no all-zero rows under cosine distance",
            "xyz = umap.UMAP(n_components=3, n_neighbors=15, metric='cosine', random_state=0, low_memory=True).fit_transform(X).astype(np.float32)",
            "build = os.path.join(OUT, 'build'); os.makedirs(build, exist_ok=True)",
            "nodes.to_parquet(os.path.join(build, HUMAN_SPACE.nodes), index=False)",
            "np.savez(os.path.join(build, HUMAN_SPACE.graph), xyz=xyz, gene_ids=nodes.gene_id.to_numpy(dtype=str))",
            "validate_graph(nodes, os.path.join(build, HUMAN_SPACE.graph))",
            "print(X.shape[1], 'features; layout', xyz.shape)")
    nb.md("## 2 · Licence per column and the attribution notice",
          "From the documentary partition (`results/human_reference_terms_2026_10_08_v5/`): GTEx public data may be "
          "redistributed with source and date acknowledged, no implied endorsement, and a notice when the version is "
          "not current; UniProt is CC BY 4.0; the Human Protein Atlas is CC BY-SA 4.0.")
    nb.code("lic = {c: ('SA' if c.startswith('HPA_') else 'O') for c in nodes.columns}",
            "notice = '\\n'.join([",
            "    f'ParaScape human gene pack {VERSION}',",
            "    '',",
            "    'Genes and expression: GTEx Analysis Release V10 (dbGaP phs000424), median TPM, GENCODE 39 / GRCh38.p13.',",
            "    '  The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund of the Office of the Director of the NIH.',",
            "    '  Data retrieved 2026-10-08 from https://gtexportal.org. Use does not imply endorsement by GTEx or the NIH.',",
            "    '  This pack pins GTEx V10; newer GTEx releases may exist (check https://gtexportal.org/home/releaseInfoPage).',",
            "    'Gene-protein links: UniProt Knowledgebase release 2026_03, reviewed human entries, CC BY 4.0 (https://www.uniprot.org).',",
            "    'Annotations (HPA_ columns): Human Protein Atlas (https://www.proteinatlas.org), CC BY-SA 4.0;',",
            "    '  those columns, and works adapted from them, must be shared under the same licence.',",
            "    '',",
            "    'Licence per column: ' + ', '.join(f'{c}={v}' for c, v in lic.items()),",
            "]) + '\\n'",
            "open(os.path.join(build, 'NOTICE.txt'), 'w').write(notice)",
            "print(notice)")
    nb.md("## 3 · The pack")
    nb.code("extra = {'NOTICE.txt': os.path.join(build, 'NOTICE.txt'),",
            "         'column_metadata.json': os.path.join(SPACE_DIR, 'column_metadata.json'),",
            "         'hs_protein_gene_links.parquet': os.path.join(SPACE_DIR, 'hs_protein_gene_links.parquet')}",
            "archive = os.path.join(OUT, f'parascape_Hs_{VERSION}.zip')",
            "sha = build_pack(HUMAN_SPACE, build, VERSION, lic, archive, extra=extra)",
            "manifest = packs.verify(archive, sha)",
            "open(os.path.join(OUT, 'SHA256'), 'w').write(f'{sha}  {os.path.basename(archive)}\\n')",
            "print(os.path.getsize(archive) // 1024, 'KiB;', sha)",
            "sorted(manifest['files'])")
    nb.md("## 4 · Opening it on its own",
          "An empty cache, the pack installed with its checksum, and the gene table and layout read back through "
          "`organisms.nodes_path` / `graph_path` -- the calls the application uses.")
    nb.code("shutil.rmtree(STATE, ignore_errors=True)",
            "os.environ[paths.ENV_STATE] = STATE",
            "added = 'Hs' not in O.SPACES",
            "if added: O.SPACES['Hs'] = HUMAN_SPACE",
            "try:",
            "    where = packs.install(archive, sha, 'Hs')",
            "    back = pd.read_parquet(O.nodes_path('Hs'))",
            "    validate_graph(back, O.graph_path('Hs'))",
            "    opened = dict(installed_at=os.path.relpath(str(where), STATE), genes=len(back),",
            "                  identical=bool(back.equals(nodes)), notice=os.path.exists(where / 'NOTICE.txt'),",
            "                  protein_rows_as_genes=int(back.gene_id.isin(pd.read_parquet(where / 'hs_protein_gene_links.parquet').host_id).sum()))",
            "finally:",
            "    if added: O.SPACES.pop('Hs')",
            "    os.environ.pop(paths.ENV_STATE)",
            "    shutil.rmtree(STATE, ignore_errors=True)",
            "assert opened['identical'] and opened['notice'] and opened['protein_rows_as_genes'] == 0",
            "json.dump(dict(version=VERSION, sha256=sha, archive=os.path.basename(archive), licences=lic, opened=opened),",
            "          open(os.path.join(OUT, 'pack_report.json'), 'w'), indent=1)",
            "shutil.rmtree(build)",
            "opened")
    nb.write(NB)


if __name__ == "__main__":
    main()
