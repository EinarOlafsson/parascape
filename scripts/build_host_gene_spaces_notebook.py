"""The rat and cat gene spaces (instruction 53, WP4), as an executed notebook.

Both are Toxoplasma hosts, and both were recorded for weeks as blocked behind the VEuPathDB login. They are
not: a rat and a cat are vertebrates, their gene builds are on Ensembl, and nothing about them needs a parasite
database. That assumption was checked on 2026-11-09 and was wrong.

**What this builds is a gene space, not an answer.** One row per gene, with its id, symbol, biotype and
coordinates from the Ensembl build. No rat or cat expression, localisation, screen or infection data is on
disk, so both spaces are *declared and installable* and neither is calibrated: a space can name genes and hold
measurements long before it can answer anything, and `space_catalogue` reports that state rather than hiding
it. The cat matters out of proportion to its data, because it is the only species where Toxoplasma completes
its sexual cycle -- several slots naming feline tissue have never had a species to hang on.

    python scripts/build_host_gene_spaces_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_09"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"host_gene_spaces_{DATE}")
NB = os.path.join(REPO, "notebooks", f"host_gene_spaces_{DATE}.ipynb")
BUILDS = os.path.join("reference", "ensembl_r113")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The rat and cat gene spaces",
                          namespace={"REPO": REPO, "OUT": OUT, "BUILDS": BUILDS})
    nb.md("**What this builds (instruction 53, WP4-6).** One row per gene for four species, from their Ensembl "
          "builds: id, symbol, biotype, chromosome, start, end, strand. Rat and cat are Toxoplasma hosts; "
          "Anopheles gambiae is the VECTOR that carries Plasmodium through its sexual cycle; P. vivax is the "
          "other human malaria parasite, kept apart from falciparum because it forms hypnozoites and "
          "falciparum does not.",
          "**Both were recorded as blocked behind the VEuPathDB login, and that was wrong.** A rat and a cat "
          "are vertebrates; their builds are on Ensembl and nothing about them needs a parasite database. The "
          "assumption went unchecked for weeks.",
          "**A gene space is not an answer.** No rat or cat measurement is on disk, so both spaces are "
          "declared and installable, and neither is calibrated. The cat matters out of proportion to its data: "
          "it is the only species where Toxoplasma completes its sexual cycle, so the enteroepithelial stages "
          "exist nowhere else.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import gzip, hashlib, json, os, re, pandas as pd",
            "from parascape import paths",
            "from parascape.spaces.anopheles import ANOPHELES_SPACE",
            "from parascape.spaces.cat import CAT_SPACE",
            "from parascape.spaces.rat import RAT_SPACE",
            "from parascape.spaces.stephensi import STEPHENSI_SPACE",
            "from parascape.spaces.vivax import VIVAX_SPACE",
            "root = paths.dataset_root()",
            "SOURCES = {RAT_SPACE: ('reference/ensembl_r113', 'Rattus_norvegicus.mRatBN7.2.113.gff3.gz'),",
            "           CAT_SPACE: ('reference/ensembl_r113', 'Felis_catus.Felis_catus_9.0.113.gff3.gz'),",
            "           ANOPHELES_SPACE: ('reference/ensembl_metazoa_r60', 'Anopheles_gambiae.AgamP4.60.chr.gff3.gz'),",
            "           STEPHENSI_SPACE: ('reference/ensembl_metazoa_r60', 'Anopheles_stephensi.AsteS1.60.gff3.gz'),",
            "           VIVAX_SPACE: ('reference/ensembl_protists_r63',",
            "                         'Plasmodium_vivax_gca900093555.GCA900093555v2.63.gff3.gz')}",
            "def sha(path):",
            "    return hashlib.sha256(open(path, 'rb').read()).hexdigest()",
            "paths_by_space = {s: os.path.join(root, where, f) for s, (where, f) in SOURCES.items()}",
            "for space, path in paths_by_space.items():",
            "    print(space.code, os.path.basename(path), 'present:', os.path.exists(path))",
            "hashes = {s.code: sha(p) for s, p in paths_by_space.items() if os.path.exists(p)}",
            "pd.Series(hashes).to_frame('sha256 of the build')")
    nb.md("## 1 · One row per gene",
          "Genes only, not transcripts: a gene space is keyed by gene. The id is the Ensembl stable id with no "
          "version, which is what the space's pattern matches.")
    nb.code("def genes_of(path):",
            "    rows = []",
            "    with gzip.open(path, 'rt') as handle:",
            "        for line in handle:",
            "            if line.startswith('#'):",
            "                continue",
            "            field = line.rstrip('\\n').split('\\t')",
            "            if len(field) < 9 or field[2] != 'gene':",
            "                continue",
            "            attributes = field[8]",
            "            found = re.search(r'ID=gene:([^;]+)', attributes)",
            "            name = re.search(r'Name=([^;]+)', attributes)",
            "            biotype = re.search(r'biotype=([^;]+)', attributes)",
            "            rows.append(dict(gene_id=(found.group(1) if found else ''),",
            "                             gene_name=(name.group(1) if name else ''),",
            "                             biotype=(biotype.group(1) if biotype else ''),",
            "                             seqid=field[0], start=int(field[3]), end=int(field[4]),",
            "                             strand=field[6]))",
            "    frame = pd.DataFrame(rows)",
            "    return frame[frame.gene_id != ''].drop_duplicates('gene_id').reset_index(drop=True)",
            "tables = {space: genes_of(path) for space, path in paths_by_space.items() if os.path.exists(path)}",
            "summary = pd.DataFrame([dict(space=s.code, species=s.species, genes=len(t),",
            "                             protein_coding=int((t.biotype == 'protein_coding').sum()),",
            "                             named=int((t.gene_name != '').sum()),",
            "                             ids_canonical=bool(t.gene_id.map(s.matches).all()))",
            "                        for s, t in tables.items()])",
            "summary")
    nb.md("## 2 · What Bgee says each gene is expressed in",
          "The same two columns the mouse space carries, so the three host spaces answer the same question the "
          "same way: the best expression score over the tissues a gene is called present in, and how many "
          "tissues those are. A gene called absent everywhere keeps a score and a count of zero -- Bgee tested "
          "it, so that is a measurement, not a gap. A gene Bgee never tested stays missing.")
    nb.code("BGEE = os.path.join('host', 'bgee_expression')",
            "bgee_files = {RAT_SPACE: 'Rattus_norvegicus_expr_simple.tsv.gz',",
            "              CAT_SPACE: 'Felis_catus_expr_simple.tsv.gz'}",
            "for space, table in list(tables.items()):",
            "    if space not in bgee_files:",
            "        continue                       # Bgee covers vertebrates; the mosquito and the parasite are not in it",
            "    path = os.path.join(root, BGEE, bgee_files[space])",
            "    if not os.path.exists(path):",
            "        print(space.code, 'no Bgee file'); continue",
            "    calls = pd.read_csv(path, sep='\t')",
            "    calls.columns = [c.strip('\"') for c in calls.columns]",
            "    calls['score'] = pd.to_numeric(calls['Expression score'], errors='coerce')",
            "    present = calls[calls.Expression == 'present']",
            "    top = present.groupby('Gene ID').score.max()",
            "    count = present.groupby('Gene ID').size()",
            "    tested = calls['Gene ID'].drop_duplicates()",
            "    table['bgee_top_score'] = table.gene_id.map(top)",
            "    table['bgee_tissues_present'] = table.gene_id.map(count)",
            "    # Tested and never present is a measured zero; never tested stays missing.",
            "    was_tested = table.gene_id.isin(set(tested))",
            "    for column in ('bgee_top_score', 'bgee_tissues_present'):",
            "        table.loc[was_tested & table[column].isna(), column] = 0",
            "    hashes[space.code + '_bgee'] = sha(path)",
            "    print(f'{space.code}: {int(was_tested.sum())} of {len(table)} genes tested by Bgee, '",
            "          f'{int(table.bgee_tissues_present.gt(0).sum())} expressed somewhere, '",
            "          f'{calls[\'Anatomical entity ID\'].nunique()} tissues')",
            "pd.DataFrame({s.code: t[['bgee_top_score', 'bgee_tissues_present']].notna().sum()",
            "              for s, t in tables.items() if 'bgee_top_score' in t}).T")
    nb.md("## 3 · The human ortholog of each gene",
          "What lets annotation transfer from human the way the mouse space does. One-to-one only: a gene with "
          "several human counterparts has no single one to copy from, and picking one would invent a certainty "
          "the homology does not have. The count of the other kinds is kept so the loss is visible.")
    nb.code("ORTHO = {RAT_SPACE: 'rn_human_orthologs.tsv', CAT_SPACE: 'fc_human_orthologs.tsv'}",
            "for space, table in tables.items():",
            "    if space not in ORTHO:",
            "        continue                       # a human ortholog means nothing for a mosquito or a parasite",
            "    path = os.path.join(OUT, ORTHO[space])",
            "    if not os.path.exists(path):",
            "        print(space.code, 'no ortholog file'); continue",
            "    links = pd.read_csv(path, sep='\t')",
            "    one = links[links.homology_type == 'ortholog_one2one'].drop_duplicates('gene_id')",
            "    table['human_ortholog'] = table.gene_id.map(one.set_index('gene_id').human_gene_id)",
            "    table['human_ortholog_high_confidence'] = table.gene_id.map(",
            "        one.set_index('gene_id').is_high_confidence).astype('Float64')",
            "    other = links[links.homology_type != 'ortholog_one2one'].gene_id.nunique()",
            "    print(f'{space.code}: {int(table.human_ortholog.notna().sum())} genes with a one-to-one human '",
            "          f'ortholog; {other} more have only one-to-many or many-to-many and are left empty')",
            "pd.DataFrame({s.code: t[['human_ortholog', 'human_ortholog_high_confidence']].notna().sum()",
            "              for s, t in tables.items() if 'human_ortholog' in t}).T")
    nb.md("## 4 · What the parasite space already has on disk",
          "The P. vivax arm of GSE220039 -- the same humanised-liver series read for falciparum, whose vivax "
          "half had never been opened. **Controls only**: the series pairs them with a PI4K-inhibitor arm, "
          "which is a drug experiment and a different question. Counts are renormalised to the PARASITE total, "
          "because the parasite is a small minority of a liver biopsy.")
    nb.code("import numpy as np",
            "PV_LIVER = os.path.join('transcription', 'RNAseq', 'GSE220039',",
            "                        'GSE220039_PvHsMmu_GRCh38_GRCm39_P01v58_raw_counts.csv.gz')",
            "pv_path = os.path.join(root, PV_LIVER)",
            "if VIVAX_SPACE in tables and os.path.exists(pv_path):",
            "    counts = pd.read_csv(pv_path)",
            "    counts = counts.rename(columns={counts.columns[0]: 'gene_id'})",
            "    counts['gene_id'] = counts.gene_id.astype(str).str.strip()",
            "    parasite = counts[counts.gene_id.str.startswith('PVP01_')].set_index('gene_id')",
            "    controls = [c for c in parasite.columns if 'Cont' in c]",
            "    total = parasite[controls].to_numpy().sum()",
            "    cpm = 1e6 * parasite[controls].sum(axis=1) / total",
            "    table = tables[VIVAX_SPACE]",
            "    table['expr_pv_liver'] = table.gene_id.map(np.log2(cpm + 1).round(4))",
            "    hashes['Pv_liver'] = sha(pv_path)",
            "    print(f'Pv: {len(controls)} control samples, {int(total):,} parasite reads, '",
            "          f'{int(table.expr_pv_liver.notna().sum())} of {len(table)} genes with a value')",
            "    print('   drug arm excluded:', [c for c in parasite.columns if 'PI4K' in c])",
            "else:",
            "    print('no vivax liver counts on disk')",
            "tables[VIVAX_SPACE].head(3) if VIVAX_SPACE in tables else 'no vivax table'")
    nb.md("## 5 · What the vector space already has on disk",
          "The Anopheles arm of GSE284537 -- the same infected-midgut series read for the Plasmodium "
          "mosquito-stage slot, whose mosquito half had never been opened. The **control (dsGFP)** arms only, "
          "pseudo-bulked per timepoint: the series pairs each with an ecdysone-receptor knockdown, which is a "
          "perturbation and a different question. These are the mosquito's own genes responding to infection.",
          "**Named `anopheles_midgut_*`, not `expr_midgut_*`.** The Plasmodium table already carries "
          "`expr_midgut_*` from the PARASITE arm of the same files, and one name meaning two species' "
          "expression is how a join goes quietly wrong -- the same collision the mouse and human protein "
          "columns hit an hour earlier.")
    nb.code("import glob, gzip, re",
            "AG_DIR = os.path.join('transcription', 'scRNAseq', 'GSE284537')",
            "ag_dir = os.path.join(root, AG_DIR)",
            "matrices = sorted(glob.glob(os.path.join(ag_dir, '*_AnG_*matrix.mtx.gz')))",
            "if ANOPHELES_SPACE in tables and matrices:",
            "    per_timepoint = {}",
            "    for matrix in matrices:",
            "        tag = re.search(r'_AnG_r\\d-(\\w+)-GFP', os.path.basename(matrix))",
            "        features = matrix.replace('matrix.mtx.gz', 'features.tsv.gz')",
            "        if tag is None or not os.path.exists(features):",
            "            continue",
            "        with gzip.open(features, 'rt') as handle:",
            "            genes = [line.split('\\t')[0] for line in handle.read().splitlines()]",
            "        totals = np.zeros(len(genes))",
            "        with gzip.open(matrix, 'rt') as handle:",
            "            line = handle.readline()",
            "            while line.startswith('%'):",
            "                line = handle.readline()",
            "            for line in handle:",
            "                gene, _cell, value = line.split()",
            "                totals[int(gene) - 1] += int(value)",
            "        per_timepoint.setdefault(tag.group(1), []).append(pd.Series(totals, index=genes))",
            "    table = tables[ANOPHELES_SPACE]",
            "    for timepoint in ('36h', '48h', '4d', '7d'):",
            "        replicates = per_timepoint.get(timepoint)",
            "        if not replicates:",
            "            continue",
            "        summed = sum(replicates)",
            "        if summed.sum() <= 0:",
            "            continue",
            "        table[f'anopheles_midgut_{timepoint}'] = table.gene_id.map(",
            "            np.log2(1e6 * summed / summed.sum() + 1).round(4))",
            "    made = [c for c in table.columns if c.startswith('anopheles_midgut_')]",
            "    print(f'Ag: {len(made)} timepoints, {int(table[made].notna().all(axis=1).sum())} of {len(table)}'",
            "          f' genes with a value in each')",
            "else:",
            "    print('no Anopheles matrices on disk')",
            "tables[ANOPHELES_SPACE].head(3) if ANOPHELES_SPACE in tables else 'no vector table'")
    nb.md("## 6 · Validate and write",
          "`spaces.validate_nodes` is the same check every space passes: unique, non-missing, canonical ids. "
          "A space whose ids do not match its own pattern is refused here rather than at first use.")
    nb.code("from parascape.spaces import validate_nodes",
            "written = {}",
            "for space, table in tables.items():",
            "    validate_nodes(space, table)",
            "    where = os.path.join(OUT, space.nodes)",
            "    table.to_parquet(where, index=False)",
            "    written[space.code] = dict(rows=len(table), columns=list(table.columns), file=space.nodes,",
            "                               build=SOURCES[space], build_sha256=hashes.get(space.code, ''))",
            "json.dump(written, open(os.path.join(OUT, 'manifest.json'), 'w'), indent=1)",
            "print('written:', {k: v['rows'] for k, v in written.items()})",
            "pd.DataFrame(written).T[['rows', 'file', 'build']]")
    nb.md("## 7 · What these spaces cannot do yet",
          "Neither holds a measurement, so neither can answer a question or be calibrated. What each needs "
          "first is stated here so the gap is a list rather than a feeling.")
    nb.code("needs = pd.DataFrame([",
            "    dict(space='Rn', first_need='an infection series: no rat Toxoplasma transcriptome is on disk',",
            "         then='a pack, so the space can be installed and leave results/'),",
            "    dict(space='Fc', first_need='the enteroepithelial series -- the stages that exist in no other host',",
            "         then='a pack, so the space can be installed and leave results/'),",
            "    dict(space='Ag', first_need='the salivary gland, where the sporozoite waits; the midgut is now covered',",
            "         then='a pack, so the space can be installed and leave results/'),",
            "    dict(space='As', first_need='any stephensi measurement at all: none is on disk',",
            "         then='the urban-breeding and insecticide-resistance data this species is studied for'),",
            "    dict(space='Pv', first_need='hypnozoite data -- the dormant stage falciparum does not have',",
            "         then='a pack, so the space can be installed and leave results/'),",
            "])",
            "needs.to_csv(os.path.join(OUT, 'what_each_space_needs.tsv'), sep='\\t', index=False)",
            "needs")
    nb.write(NB)


if __name__ == "__main__":
    main()
