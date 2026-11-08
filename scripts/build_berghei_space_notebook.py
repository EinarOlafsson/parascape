"""The Plasmodium berghei ANKA gene space, PlasmoGEM in its native home (instruction 53, WP6), as a notebook.

Genes: every PBANKA gene of the PlasmoDB-68 annotated proteome on disk. Measurements: the PlasmoGEM blood-stage
knockout screen (Bushell et al. 2017, Cell, PMC5509546 Table S1): phenotype class and relative growth rate with its
interval, keyed by current PBANKA identifiers, with the P. falciparum ortholog the authors name.

Checks: (1) WP6 acceptance -- the P. falciparum `pb_transferred_phenotype` column equals the projection of this
space through the named orthologs; (2) a benchmark on P. berghei truth -- does the P. falciparum saturation
mutagenesis index of the named ortholog (piggyBac MIS, Zhang 2018) separate essential from dispensable P. berghei
knockouts (AUC)?

    python scripts/build_berghei_space_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"berghei_space_{DATE}")
NB = os.path.join(REPO, "notebooks", f"berghei_space_{DATE}.ipynb")
SRC = "/media/carruthers/mnt3/claude/toxoplasma_projects/datasets"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The Plasmodium berghei gene space", namespace={"REPO": REPO, "OUT": OUT, "SRC": SRC})
    nb.md("**What this builds (instruction 53, WP6).** One row per P. berghei ANKA gene, with the PlasmoGEM "
          "blood-stage knockout screen in its native home. Until now PlasmoGEM existed in ParaScape only as a "
          "transfer into P. falciparum (`pb_transferred_*`).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import hashlib, json, os, re, numpy as np, pandas as pd",
            "from sklearn.metrics import roc_auc_score",
            "from parascape.spaces import validate_nodes",
            "from parascape.spaces.berghei import BERGHEI_SPACE",
            "def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()",
            "FASTA = os.path.join(SRC, 'pan_apicomplexan/PbergheiANKA/PbergheiANKA_annotation/PlasmoDB-68_PbergheiANKA_AnnotatedProteins.fasta')",
            "GEM = os.path.join(SRC, 'recovered_processed_sources/pf_berghei_transfer/mmc1.xlsx')",
            "hashes = {'proteome': sha(FASTA), 'plasmogem': sha(GEM)}",
            "genes = {}",
            "for line in open(FASTA):",
            "    if line.startswith('>'):",
            "        f = dict(kv.split('=', 1) for kv in (p.strip() for p in line[1:].split('|')[1:]) if '=' in kv)",
            "        g = f.get('gene')",
            "        if g and g not in genes:",
            "            genes[g] = dict(product=f.get('gene_product'), protein_length=int(f.get('protein_length', 0) or 0))",
            "nodes = pd.DataFrame([dict(gene_id=g, **v) for g, v in genes.items()]).sort_values('gene_id').reset_index(drop=True)",
            "odd = nodes[~nodes.gene_id.str.fullmatch(r'PBANKA_\\d{7}')]",
            "nodes = nodes[nodes.gene_id.str.fullmatch(r'PBANKA_\\d{7}')].reset_index(drop=True)",
            "print(len(nodes), 'PBANKA genes;', len(odd), 'non-PBANKA identifiers left out:', odd.gene_id.tolist()[:5])")
    nb.md("## 1 · PlasmoGEM in its own space",
          "Rows keyed by the current PBANKA identifier; 'Insufficient data' is kept as stated, not turned into a class.")
    nb.code("gem = pd.read_excel(GEM)",
            "gem = gem.rename(columns={'P. berghei current ID': 'gene_id', 'P. falciparum ID': 'pf_ortholog', 'Phenotype': 'plasmogem_phenotype',",
            "                          'Relative growth rate': 'plasmogem_growth_rate', 'RGR CI low': 'plasmogem_rgr_ci_low',",
            "                          'RGR CI high': 'plasmogem_rgr_ci_high', 'Confidence': 'plasmogem_confidence'})",
            "no_id = int(gem.gene_id.isna().sum())",
            "gem = gem.dropna(subset=['gene_id']).assign(gene_id=lambda d: d.gene_id.astype(str).str.strip())",
            "dup = int(gem.gene_id.duplicated().sum())",
            "gem = gem.drop_duplicates('gene_id').set_index('gene_id')",
            "missing = sorted(set(gem.index) - set(nodes.gene_id))",
            "for c in ('pf_ortholog', 'plasmogem_phenotype', 'plasmogem_growth_rate', 'plasmogem_rgr_ci_low', 'plasmogem_rgr_ci_high', 'plasmogem_confidence'):",
            "    nodes[c] = nodes.gene_id.map(gem[c])",
            "validate_nodes(BERGHEI_SPACE, nodes)",
            "nodes.to_parquet(os.path.join(OUT, 'pb_gene_nodes.parquet'), index=False)",
            "print(int(nodes.plasmogem_phenotype.notna().sum()), 'genes with a PlasmoGEM phenotype;', no_id, 'rows without a PBANKA id;', dup, 'duplicate rows;',",
            "      len(missing), 'PlasmoGEM genes not in the PlasmoDB-68 annotation:', missing[:5])",
            "nodes.plasmogem_phenotype.value_counts(dropna=False).to_frame('genes')")
    nb.md("## 2 · WP6 acceptance: the P. falciparum transfer is the projection of this space",
          "`pf_nodes.pb_transferred_phenotype` was built from the same table through the named orthologs. For every "
          "P. falciparum gene carrying it, the P. berghei gene naming that ortholog must have the same phenotype.")
    nb.code("pf = pd.read_parquet(os.path.join(REPO, 'parascape/data/pf_nodes.parquet'), columns=['gene_id', 'pb_transferred_phenotype', 'piggybac_mis'])",
            "proj = nodes.dropna(subset=['pf_ortholog']).assign(pf=lambda d: d.pf_ortholog.astype(str).str.split(r'[;,\\s]+')).explode('pf')",
            "proj = proj[proj.pf.str.match(r'PF3D7_')]",
            "per_pf = proj.groupby('pf').plasmogem_phenotype.agg(lambda s: s.iloc[0] if s.nunique() == 1 else 'AMBIGUOUS')",
            "carried = pf.dropna(subset=['pb_transferred_phenotype']).set_index('gene_id').pb_transferred_phenotype.astype(str)",
            "both = carried.index.intersection(per_pf.index)",
            "agree = (carried[both] == per_pf[both].astype(str))",
            "acceptance = dict(pf_genes_carrying_transfer=len(carried), reached_from_pb_space=len(both),",
            "                  identical=int(agree.sum()), differing=int((~agree).sum()),",
            "                  differing_examples=[(g, carried[g], per_pf[g]) for g in both[~agree.values][:5]])",
            "acceptance")
    nb.md("## 3 · Benchmark on P. berghei truth: does the P. falciparum mutagenesis index predict essentiality?",
          "Essential (1) against dispensable (0) P. berghei knockouts, slow and insufficient data left out; the "
          "piggyBac mutagenesis index (low = the gene resists disruption) of the named P. falciparum ortholog as "
          "the score, negated so higher means more likely essential.")
    nb.code("ess = nodes[nodes.plasmogem_phenotype.isin(['Essential', 'Dispensable'])].dropna(subset=['pf_ortholog']).copy()",
            "ess['pf'] = ess.pf_ortholog.astype(str).str.split(r'[;,\\s]+').str[0]",
            "mis = pf.set_index('gene_id').piggybac_mis",
            "ess['mis'] = ess.pf.map(mis)",
            "ok = ess.dropna(subset=['mis'])",
            "y = (ok.plasmogem_phenotype == 'Essential').astype(int)",
            "bench = dict(genes=len(ok), essential=int(y.sum()), auc=round(float(roc_auc_score(y, -ok.mis)), 3))",
            "json.dump(dict(acceptance={k: v for k, v in acceptance.items() if k != 'differing_examples'}, benchmark=bench,",
            "               genes=len(nodes), sha256=hashes), open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "bench")
    nb.write(NB)


if __name__ == "__main__":
    main()
