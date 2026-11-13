"""The Cryptosporidium parvum gene space and a first benchmark on its own truth (instruction 53, WP6).

Genes: the 3,941 Iowa II genes of the OrthoMCL table (CryptoDB `cgd` identifiers), with product, OrthoMCL group and
protein length; measured hyperLOPIT location (Guerin et al. 2023, 1,742 genes); measured stage (Tandel et al.
2019, 1,477 genes); InterPro/Pfam domains; DeepTMHMM signal peptides and transmembrane segments. Every column
records its source and its licence status; the pack is not distributed until the CryptoDB and supplement terms are
verified (licence code V).

Benchmark: does the measured localisation of a gene's Toxoplasma orthologs (hyperLOPIT, through the shared OrthoMCL
group) predict the measured Cryptosporidium location? The mapping between the two species' location vocabularies is
learned on training genes only; folds keep whole OrthoMCL groups together; compared with the commonest class.

    python scripts/build_cryptosporidium_space_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"cryptosporidium_space_{DATE}")
NB = os.path.join(REPO, "notebooks", f"cryptosporidium_space_{DATE}.ipynb")
SRC = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The Cryptosporidium parvum gene space", namespace={"REPO": REPO, "OUT": OUT, "SRC": SRC})
    nb.md("**What this builds (instruction 53, WP6).** One row per Cryptosporidium parvum Iowa II gene from the "
          "sources already on disk, each column with its source and licence status, and a first benchmark on "
          "Cryptosporidium truth: whether the localisation measured on a gene's Toxoplasma orthologs predicts the "
          "localisation measured in Cryptosporidium.",
          "**Not yet distributed:** the CryptoDB tables and the Guerin 2023 supplement need their redistribution "
          "terms verified before a pack (licence code V).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import hashlib, json, os, re, numpy as np, pandas as pd",
            "from parascape.spaces import validate_nodes",
            "from parascape.spaces.crypto import CRYPTO_SPACE",
            "def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()",
            "F = {k: os.path.join(SRC, v) for k, v in dict(genes='orthomcl_cryptosporidium_parvum_IowaII.csv',",
            "     lopit='lopit_cryptosporidium_parvum_MEASURED_Guerin2023.csv', stage='stage_measured_cryptosporidium_Tandel2019.csv',",
            "     interpro='interpro_cpar.csv', tm='deeptmhmm_cp.csv',",
            "     oocyst='stage_measured_cryptosporidium_oocyst_GSE271211.csv').items()}",
            "hashes = {k: sha(p) for k, p in F.items()}",
            "g = pd.read_csv(F['genes'])",
            "bad = g[~g.gene_source_id.str.fullmatch(r'cgd\\d+_\\d+')]",
            "g = g[g.gene_source_id.str.fullmatch(r'cgd\\d+_\\d+')]",
            "print(len(g), 'genes;', len(bad), 'rows with a non-cgd identifier left out:', bad.gene_source_id.tolist()[:5])")
    nb.md("## 1 · One row per gene")
    nb.code("nodes = pd.DataFrame({'gene_id': g.gene_source_id.values, 'product': g.gene_product.values,",
            "                      'orthogroup': g.gene_orthomcl_name.values, 'protein_length': g.protein_length.values})",
            "lopit = pd.read_csv(F['lopit']).drop_duplicates('gene_source_id').set_index('gene_source_id').MAP_location",
            "nodes['lopit_location'] = nodes.gene_id.map(lopit).replace('unassigned', np.nan)",
            "nodes['lopit_assigned'] = nodes.gene_id.map(lopit).notna()",
            "st = pd.read_csv(F['stage']).drop_duplicates('gene_id').set_index('gene_id')",
            "for c in ('measured_stage', 'asexual_pct', 'female_pct', 'logFC_female_vs_asexual', 'sporozoite_specific'):",
            "    nodes['stage_' + c if not c.startswith('measured') else c] = nodes.gene_id.map(st[c])",
            "ip = pd.read_csv(F['interpro']); ip['gene'] = ip.gene_source_id.str.lower()",
            "ip = ip.groupby('gene').agg(interpro_ids=('interpro_id', lambda s: ';'.join(sorted({x for v in s.dropna() for x in str(v).split(';')}))))",
            "nodes['interpro_ids'] = nodes.gene_id.map(ip.interpro_ids)",
            "nodes['n_interpro'] = nodes.interpro_ids.fillna('').map(lambda v: len([x for x in v.split(';') if x]))",
            "tm = pd.read_csv(F['tm'], low_memory=False).drop_duplicates('identifier').set_index('identifier')",
            "for c in ('dtm_type', 'n_signal', 'n_tm'):",
            "    nodes['tm_' + c if c != 'dtm_type' else c] = nodes.gene_id.map(tm[c])",
            "# GSE271211: the undigested oocyst inoculum, three replicates, mean log2 CPM. Only the",
            "# three time-zero samples -- the twelve digestion arms are oocysts under the excystation",
            "# trigger and belong to a perturbation slot, not to this one. Keyed on CPATCC_ locus tags",
            "# from a different assembly, so it arrives through the 2,571-gene exact-protein map; the",
            "# 1,352 genes that map does not reach stay missing rather than being filled with a zero,",
            "# because 'not measured here' and 'measured and absent' are different answers.",
            "oo = pd.read_csv(F['oocyst']).drop_duplicates('gene_id').set_index('gene_id')",
            "nodes['expr_oocyst'] = nodes.gene_id.map(oo.oocyst_log2cpm)",
            "nodes['oocyst_detected'] = nodes.gene_id.map(oo.oocyst_detected)",
            "print(f\"oocyst column: {nodes.expr_oocyst.notna().sum():,} of {len(nodes):,} genes measured, \"",
            "      f\"{int(nodes.oocyst_detected.fillna(False).sum()):,} detected\")",
            "nodes = nodes.sort_values('gene_id').reset_index(drop=True)",
            "validate_nodes(CRYPTO_SPACE, nodes)",
            "nodes.to_parquet(os.path.join(OUT, 'cp_gene_nodes.parquet'), index=False)",
            "cover = nodes.notna().sum().to_frame('genes')",
            "cover")
    nb.md("## 2 · What each column is, and its licence status")
    nb.code("SOURCE = {'gene_id': ('genes', 'identity'), 'product': ('genes', 'annotation'), 'orthogroup': ('genes', 'orthology'),",
            "          'protein_length': ('genes', 'sequence'), 'lopit_location': ('lopit', 'place (measured)'),",
            "          'lopit_assigned': ('lopit', 'place (measured)'), 'measured_stage': ('stage', 'transcript (measured)'),",
            "          'stage_asexual_pct': ('stage', 'transcript (measured)'), 'stage_female_pct': ('stage', 'transcript (measured)'),",
            "          'stage_logFC_female_vs_asexual': ('stage', 'transcript (measured)'), 'stage_sporozoite_specific': ('stage', 'transcript (measured)'),",
            "          'interpro_ids': ('interpro', 'domain (predicted)'), 'n_interpro': ('interpro', 'domain (predicted)'),",
            "          'dtm_type': ('tm', 'structure (predicted)'), 'tm_n_signal': ('tm', 'structure (predicted)'), 'tm_n_tm': ('tm', 'structure (predicted)'),",
            "          'expr_oocyst': ('oocyst', 'transcript (measured)'), 'oocyst_detected': ('oocyst', 'transcript (measured)')}",
            "REF = {'genes': 'CryptoDB / OrthoMCL (VEuPathDB)', 'lopit': 'Guerin et al. 2023 hyperLOPIT supplement',",
            "       'stage': 'Tandel et al. 2019 single-cell RNA-seq', 'interpro': 'InterPro via CryptoDB', 'tm': 'DeepTMHMM run in this project',",
            "       'oocyst': 'GSE271211 (PMID 39354600), undigested oocyst inoculum, mapped CPATCC to cgd by exact protein identity'}",
            "meta = {c: dict(organism='Cp', source=REF[k], kind=kind, source_sha256=hashes[k],",
            "                licence='O (computed here)' if k == 'tm' else 'P (public GEO deposit)' if k == 'oocyst'",
            "                        else 'V (verify redistribution terms)') for c, (k, kind) in SOURCE.items()}",
            "assert set(meta) == set(nodes.columns)",
            "json.dump(meta, open(os.path.join(OUT, 'column_metadata.json'), 'w'), indent=1)",
            "pd.DataFrame(meta).T[['source', 'kind', 'licence']]")
    nb.md("## 3 · Benchmark: Toxoplasma orthologs predict the Cryptosporidium location?",
          "Toxoplasma hyperLOPIT (`lopit_unified`) per OrthoMCL group (commonest label among the group Toxoplasma "
          "genes). Cryptosporidium locations with at least 10 genes. Five folds over whole OrthoMCL groups; in each, "
          "the Toxoplasma-to-Cryptosporidium location mapping is the commonest Cryptosporidium location per "
          "Toxoplasma location among the TRAINING genes; a held-out gene gets its group Toxoplasma location mapped.")
    nb.code("tg = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet'), columns=['gene_id', 'orthogroup', 'lopit_unified'])",
            "tg_loc = tg.dropna().groupby('orthogroup').lopit_unified.agg(lambda s: s.value_counts().index[0])",
            "truth = nodes.set_index('gene_id').lopit_location.dropna()",
            "keep = truth.value_counts()[lambda c: c >= 10].index",
            "truth = truth[truth.isin(keep)]",
            "og = nodes.set_index('gene_id').orthogroup.reindex(truth.index)",
            "src = og.map(tg_loc)",
            "from sklearn.model_selection import GroupKFold",
            "rows = []",
            "groups = og.fillna(pd.Series(truth.index, index=truth.index))",
            "for fold, (tr, te) in enumerate(GroupKFold(5).split(truth, truth, groups)):",
            "    t_tr, s_tr = truth.iloc[tr], src.iloc[tr]",
            "    both = s_tr.notna()",
            "    mapping = t_tr[both].groupby(s_tr[both]).agg(lambda s: s.value_counts().index[0])",
            "    common = t_tr.value_counts().index[0]",
            "    for gene in truth.index[te]:",
            "        p = mapping.get(src[gene]) if pd.notna(src[gene]) else None",
            "        rows.append(dict(gene=gene, fold=fold, truth=truth[gene], tg_location=src[gene], predicted=p,",
            "                         commonest=common))",
            "B = pd.DataFrame(rows)",
            "ans = B.dropna(subset=['predicted'])",
            "bench = dict(labelled=len(B), classes=int(truth.nunique()), with_tg_ortholog_location=int(B.tg_location.notna().sum()),",
            "             answered=len(ans), right=round(float((ans.predicted == ans.truth).mean()), 3),",
            "             chance_on_answered=round(float((ans.commonest == ans.truth).mean()), 3))",
            "B.to_csv(os.path.join(OUT, 'benchmark_rows.tsv'), sep='\\t', index=False)",
            "bench")
    nb.code("per = ans.assign(right=ans.predicted == ans.truth).groupby('truth').agg(genes=('gene', 'size'), right=('right', 'mean')).round(2)",
            "json.dump(dict(benchmark=bench, genes=len(nodes), sha256=hashes), open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "per.sort_values('genes', ascending=False)")
    nb.write(NB)


if __name__ == "__main__":
    main()
