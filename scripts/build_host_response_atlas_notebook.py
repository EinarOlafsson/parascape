"""Host response atlas (instruction 71, product 3), as an executed notebook.

The same within-study contrasts as the parasite-side signatures (the contrast ladder), read on the HOST side of
the competitive salmon reference: host transcripts summed to host genes, counts per million over the host
library of that run, log2(CPM + 1), treated minus baseline group means, with a standard error where both arms
are replicated. Per host species (human or mouse, whichever the run's reads are), never pooled.

    python scripts/build_host_response_atlas_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"host_response_atlas_{DATE}")
NB = os.path.join(REPO, "notebooks", f"host_response_atlas_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Host response atlas: what each parasite perturbation does to the host cell",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this is.** For every within-study contrast of the contrast ladder (a mutant against its "
          "parent, a treatment against its control, in the same host cells), how the HOST transcriptome "
          "changes -- read from the host side of the same competitive quantification that gave the "
          "parasite values (instruction 71, product 3).",
          "**The check.** GSE125120 compares infected with bystander cells for wild-type and MYR1-knockout "
          "parasites. MYR1 exports the dense-granule effectors that rewire host transcription, so the "
          "infected-cell response should shrink without it -- measured below, not assumed.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code(
        "import gzip, os, re",
        "from pathlib import Path",
        "import numpy as np, pandas as pd",
        "BASE = Path('/nas_mnt/carruthers/Einar/parascape_omics')",
        "snap = sorted(p for p in (BASE / 'app_snapshots').iterdir()",
        "              if re.fullmatch(r'\\d{8}T\\d{6}Z-\\d+', p.name) and (p / 'rna/sample_design.tsv').exists())[-1]",
        "design = pd.read_csv(snap / 'rna/sample_design.tsv', sep='\\t', dtype=str)",
        "design = design[design.processing_status == 'verified'].drop_duplicates('run_accession')",
        "ladder = pd.read_csv(os.path.join(REPO, 'results/contrast_ladder_2026_11_06/contrasts.tsv'), sep='\\t')",
        "tx2gene = pd.read_csv(BASE / 'references/ME49_host/tx2gene.tsv', sep='\\t', header=None, names=['tx', 'gene'])",
        "host_tx = tx2gene[tx2gene.gene.str.startswith('HOST:')].set_index('tx').gene",
        "def members(text):",
        "    return [x.strip() for x in str(text).split('||') if x.strip()]",
        "groups = {g for col in ('baseline_members', 'treated_members') for v in ladder[col].dropna() for g in members(v)}",
        "runs = design[design.condition_group.isin(groups)]",
        "print(snap.name, len(ladder), 'contrasts;', len(runs), 'runs to read;', host_tx.nunique(), 'host genes in the reference')")

    nb.md("## 1 · Host counts per run", "Salmon NumReads of host transcripts summed to host genes.")
    nb.code(
        "def host_counts(study, run):",
        "    path = BASE / 'rnaseq' / study / run / 'salmon_ME49_host' / 'quant.sf.gz'",
        "    if not path.exists():",
        "        return None",
        "    q = pd.read_csv(path, sep='\\t', usecols=['Name', 'NumReads'])",
        "    q = q[q.Name.isin(host_tx.index)]",
        "    return q.groupby(q.Name.map(host_tx)).NumReads.sum()",
        "counts = {}",
        "for r in runs.itertuples():",
        "    c = host_counts(r.study, r.run_accession)",
        "    if c is not None:",
        "        counts[r.run_accession] = c",
        "H = pd.DataFrame(counts).fillna(0.0)",
        "human = H.index.str.startswith('HOST:ENSG'); mouse = H.index.str.startswith('HOST:ENSMUSG')",
        "species = pd.Series(np.where(H[human].sum() >= H[mouse].sum(), 'human', 'mouse'), index=H.columns)",
        "host_reads = pd.Series(np.where(species == 'human', H[human].sum(), H[mouse].sum()), index=H.columns)",
        "print(H.shape[1], 'runs read;', species.value_counts().to_dict(), '; median host reads', int(host_reads.median()))")

    nb.md("## 2 · Host signatures, one per contrast",
          "Within the run's own species; CPM over that species' genes; genes kept with mean CPM >= 1 across the "
          "contrast's runs. A contrast with fewer than 100,000 host reads in any run is skipped.")
    nb.code(
        "group_runs = design.groupby('condition_group').run_accession.apply(list)",
        "sigs, ses, rows = {}, {}, []",
        "for i, c in ladder.iterrows():",
        "    base = [r for g in members(c.baseline_members) for r in group_runs.get(g, []) if r in H.columns]",
        "    treat = [r for g in members(c.treated_members) for r in group_runs.get(g, []) if r in H.columns]",
        "    if not base or not treat:",
        "        continue",
        "    sp = species[base + treat]",
        "    if sp.nunique() != 1 or host_reads[base + treat].min() < 1e5:",
        "        continue",
        "    mask = human if sp.iloc[0] == 'human' else mouse",
        "    X = H.loc[mask, base + treat]",
        "    cpm = X / X.sum() * 1e6",
        "    cpm = cpm[cpm.mean(axis=1) >= 1]",
        "    L = np.log2(cpm + 1)",
        "    lfc = L[treat].mean(axis=1) - L[base].mean(axis=1)",
        "    se = np.sqrt(L[treat].var(axis=1, ddof=1) / len(treat) + L[base].var(axis=1, ddof=1) / len(base)) \\",
        "        if len(treat) > 1 and len(base) > 1 else pd.Series(np.nan, index=L.index)",
        "    name = f'host_lfc_{c.study}_{i}'",
        "    sigs[name], ses[name] = lfc, se",
        "    rows.append(dict(column=name, study=c.study, species=sp.iloc[0], axis=c.axis, baseline=c.baseline_level,",
        "                     treated=c.treated_level, background=c.background, n_baseline=len(base), n_treated=len(treat),",
        "                     host_genes=len(lfc), changed_2fold=int((lfc.abs() >= 1).sum())))",
        "registry = pd.DataFrame(rows)",
        "signatures = pd.DataFrame(sigs)",
        "signatures.to_parquet(os.path.join(OUT, 'host_signatures.parquet'))",
        "pd.DataFrame(ses).to_parquet(os.path.join(OUT, 'host_signature_se.parquet'))",
        "registry.to_csv(os.path.join(OUT, 'registry.tsv'), sep='\\t', index=False)",
        "print(len(registry), 'host signatures in', registry.study.nunique(), 'studies;', registry.species.value_counts().to_dict())",
        "registry.sort_values('changed_2fold', ascending=False).head(12)")

    nb.md("## 3 · The check: the infected-cell response without MYR1",
          "GSE125120, infected against bystander cells. If MYR1-exported effectors drive the host response, the "
          "wild-type contrast changes more host genes than the MYR1-knockout one, and the two agree only "
          "partly. MYC, whose induction by infection is reported to need MYR1, is shown for both.")
    nb.code(
        "m = registry[(registry.study == 'GSE125120') & (registry.axis == 'condition')]",
        "check = m[['column', 'background', 'n_baseline', 'n_treated', 'host_genes', 'changed_2fold']].copy()",
        "myc = [g for g in signatures.index if g.startswith('HOST:ENSG00000136997')]",
        "check['MYC_lfc'] = [round(float(signatures.loc[myc[0], c]), 2) if myc and pd.notna(signatures.loc[myc[0], c]) else np.nan for c in check.column]",
        "if len(check) == 2:",
        "    a, b = check.column.tolist()",
        "    both = signatures[[a, b]].dropna()",
        "    check['spearman_between'] = round(float(both[a].corr(both[b], method='spearman')), 3)",
        "check.to_csv(os.path.join(OUT, 'myr1_check.tsv'), sep='\\t', index=False)",
        "check")
    nb.write(NB)


if __name__ == "__main__":
    main()
