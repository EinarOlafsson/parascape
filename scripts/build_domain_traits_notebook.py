"""Instruction 72, protein-domain (InterPro / Pfam) tier-2 traits of the ParaScape information map.

Writes notebooks/domain_traits_2026_11_04.ipynb (executed) and results/domain_traits_2026_11_04/.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('Protein-domain traits of the information map')
nb.md('# Protein-domain traits (InterPro and Pfam)',
      'Instruction 72, tier 2 (Einar, 2026-11-04: "a very wide matrix, with each measurement of the gene as a tier 2 column, '
      'and then columns like Kinase domain, host facing"). One row per Toxoplasma gene (the 8,140 genes of `nodes.parquet`).',
      'A **domain** is a recognisable protein module (for example the kinase module, or the EF-hand that binds calcium). '
      'InterPro and Pfam are databases of such modules; software scans each protein sequence and says which modules it contains. '
      'This notebook turns that into columns in three groups:',
      '1. **One yes/no column per domain family** that occurs in enough genes (InterPro entries and Pfam families separately).\n'
      '2. **Named classes** a biologist asks for by name: kinase, AP2 transcription factor, rhoptry kinase, SRS surface antigen, '
      'EF-hand, WD40, ankyrin, RNA-binding, proteasome/ubiquitin, transmembrane-rich.\n'
      '3. **Counts**, plus `domain_novel`: the protein has no InterPro and no Pfam domain at all. This matters for the host-facing '
      'question because secreted effectors are often unlike anything in the databases.',
      '**Important caveat.** A domain call is itself a **prediction from the amino-acid sequence**, not an experiment. A gene with no '
      'domain is "not recognised", which is not the same as "has no domain".')

nb.md('## 1 · Inputs and what the delimited columns look like')
nb.code('import re',
        'import numpy as np, pandas as pd',
        'from pathlib import Path',
        'from scipy.stats import fisher_exact',
        "OUT = Path('results/domain_traits_2026_11_04'); OUT.mkdir(parents=True, exist_ok=True)",
        "nodes = pd.read_parquet('parascape/data/nodes.parquet').set_index('gene_id')",
        "print(len(nodes), 'genes')",
        "for c in ['interpro_id', 'interpro_desc', 'pfam_id', 'n_interpro', 'has_domain']:",
        "    print(c, '| dtype', nodes[c].dtype, '| examples:', [v for v in nodes[c].head(40) if v != ''][:3])",
        "print('Delimiter: a semicolon (;) between entries. IDs and descriptions are separate columns in the same order. Empty string = none.')",
        "print('interpro_id empty:', int((nodes.interpro_id == '').sum()), '| pfam_id empty:', int((nodes.pfam_id == '').sum()),",
        "      '| has_domain==0:', int((nodes.has_domain == 0).sum()))",
        "split = lambda s: [x for x in s.split(';') if x]",
        "ip = nodes.interpro_id.map(split); pf = nodes.pfam_id.map(split)",
        "print('n_interpro equals the number of listed InterPro IDs for', int((ip.map(len) == nodes.n_interpro).sum()), 'of', len(nodes), 'genes')",
        "print('max InterPro IDs per gene:', int(ip.map(len).max()), '| max Pfam per gene:', int(pf.map(len).max()))")

nb.md('## 2 · Where the data come from (provenance) and whether the lists are complete',
      'Both columns come from the one dataset registered as `interpro`. The registry holds no PMID or accession for it: it is a '
      'ToxoDB gene-attribute download (InterPro scan results, which include Pfam), file `datasets/interpro_tgon.csv`. '
      'That raw file also carries Pfam descriptions, which `nodes.parquet` does not, so it is used (read-only) for the '
      'descriptions and to check that `nodes.parquet` is complete.')
nb.code('from parascape.datasets import provenance',
        "for col in ('interpro_id', 'pfam_id'):",
        "    p = provenance(col, 'Tg')",
        "    print(col, '->', p.key, '|', p.name, '| level', p.level, '| pmid', p.pmid, '| accession', p.accession, '| path', p.path)",
        "    print('   note:', p.note)",
        "print('url:', provenance('interpro_id', 'Tg').url[:110], '...')")
nb.code("RAW = '/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/interpro_tgon.csv'",
        "raw = pd.read_csv(RAW, keep_default_na=False)",
        "print(len(raw), 'rows,', raw.gene_source_id.nunique(), 'distinct genes (a few genes have several transcript rows)')",
        "clean = lambda s: {x for x in s.split(';') if x and x != 'N/A'}",
        "union = lambda col: raw.groupby('gene_source_id')[col].agg(lambda s: set().union(*map(clean, s)))",
        "rip, rpf = union('interpro_id').reindex(nodes.index), union('pfam_id').reindex(nodes.index)",
        "print('genes of nodes.parquet missing from the raw file:', int(rip.isna().sum()))",
        "d_ip = [g for g, a, b in zip(nodes.index, ip, rip) if set(a) != b]",
        "d_pf = [g for g, a, b in zip(nodes.index, pf, rpf) if set(a) != b]",
        "print('genes where nodes.parquet InterPro list differs from the raw file:', len(d_ip), d_ip[:5])",
        "print('genes where nodes.parquet Pfam list differs from the raw file:', len(d_pf), d_pf[:5])",
        "at_cap = ip[ip.map(len) >= 20]",
        "print('genes whose InterPro ID list reaches 20 (a possible cap):', list(at_cap.index))",
        "for g in at_cap.index:",
        "    print('  ', g, nodes.loc[g, 'product'], '| IDs', len(ip[g]), '| descriptions', len(nodes.loc[g, 'interpro_desc'].split(';')))",
        "COVERAGE_NOTE = f'{len(d_ip)} InterPro and {len(d_pf)} Pfam lists differ from the raw file; {len(at_cap)} gene(s) at 20 IDs'")
nb.md('**Reading of the check.** If the differing counts above are zero the per-gene lists in `nodes.parquet` are the complete raw lists. '
      'A gene whose ID list stops at exactly 20 while its description list is longer has a **truncated** ID list in the source '
      'itself; that single case is reported and left as is (nothing is invented). Everything below uses `nodes.parquet`.')

nb.md('## 3 · Accession to description tables',
      'Descriptions are paired with IDs by position, using only rows where the two lists have the same length (the rare unequal '
      'row is skipped, never guessed).')
nb.code("ip_desc, pf_desc = {}, {}",
        "for col_id, col_desc, store in (('interpro_id', 'interpro_desc', ip_desc), ('pfam_id', 'pfam_desc', pf_desc)):",
        "    for a, b in zip(raw[col_id], raw[col_desc]):",
        "        A, B = [x for x in a.split(';') if x != 'N/A'], [x for x in b.split(';') if x != 'N/A']",
        "        if A and len(A) == len(B):",
        "            for x, y in zip(A, B): store.setdefault(x, y)",
        "for g in nodes.index:   # same pairing from nodes.parquet fills any InterPro accession the raw pass missed",
        "    A, B = ip[g], split(nodes.loc[g, 'interpro_desc'])",
        "    if len(A) == len(B):",
        "        for x, y in zip(A, B): ip_desc.setdefault(x, y)",
        "all_ip = {a for l in ip for a in l}; all_pf = {a for l in pf for a in l}",
        "print('InterPro accessions in nodes:', len(all_ip), '| with a description:', len(all_ip & set(ip_desc)))",
        "print('Pfam accessions in nodes:', len(all_pf), '| with a description:', len(all_pf & set(pf_desc)))")

nb.md('## 4 · Frequency cutoff',
      'A family seen in only one or two genes cannot be learned from and just widens the matrix. The cutoff is the minimum number of '
      'genes that must carry the family. The number of columns at four cutoffs is shown so the choice is visible; **10 is used**.')
nb.code("ip_n = pd.Series([a for l in ip for a in set(l)]).value_counts()",
        "pf_n = pd.Series([a for l in pf for a in set(l)]).value_counts()",
        "cut = pd.DataFrame({c: {'InterPro families': int((ip_n >= c).sum()), 'Pfam families': int((pf_n >= c).sum()),",
        "                        'total columns': int((ip_n >= c).sum() + (pf_n >= c).sum())} for c in (5, 10, 20, 50)}).T",
        "cut.index.name = 'minimum genes'",
        "CUTOFF = 10",
        "cut")

nb.md('## 5 · Group 1: one yes/no column per family',
      'Column `dom_<accession>` is True when the gene has the family. False means "not found in the sequence scan", not "proven absent".')
nb.code("keep_ip = sorted(ip_n[ip_n >= CUTOFF].index); keep_pf = sorted(pf_n[pf_n >= CUTOFF].index)",
        "sets_ip, sets_pf = ip.map(set), pf.map(set)",
        "fam = {}",
        "for a in keep_ip: fam['dom_' + a] = sets_ip.map(lambda s, a=a: a in s)",
        "for a in keep_pf: fam['dom_' + a] = sets_pf.map(lambda s, a=a: a in s)",
        "fam = pd.DataFrame(fam, index=nodes.index)",
        "lookup = pd.DataFrame(",
        "    [('dom_' + a, a, ip_desc.get(a, ''), int(ip_n[a]), 'InterPro') for a in keep_ip] +",
        "    [('dom_' + a, a, pf_desc.get(a, ''), int(pf_n[a]), 'Pfam') for a in keep_pf],",
        "    columns=['column', 'accession', 'description', 'gene_count', 'source'])",
        "print(len(keep_ip), 'InterPro +', len(keep_pf), 'Pfam =', fam.shape[1], 'family columns')",
        "print('families kept without a description:', int((lookup.description == '').sum()))",
        "lookup.head(8)")

nb.md('## 6 · Group 2: named domain classes, with the exact matching rule',
      'Each class is defined by an explicit list of InterPro and Pfam accessions (some also by words in the description). '
      'The rules were written from the lists of what exists in this dataset and then FIXED; known genes are checked afterwards and '
      'a miss is reported, not repaired. Superfamily entries that are very broad (for example "Protein kinase-like domain '
      'superfamily", which also covers non-kinase folds) are deliberately left out.',
      '* **protein_kinase**: InterPro IPR000719 (Protein kinase domain), IPR001245 and IPR020635 (Ser/Thr/Tyr kinase catalytic), '
      'IPR027916 (Protein kinase-like, Apicomplexa); Pfam PF00069, PF07714, PF14531.\n'
      '* **ap2_apiap2**: InterPro IPR001471 (AP2/ERF domain); Pfam PF00847.\n'
      '* **rhoptry_kinase_ropk**: the Apicomplexa kinase-like entries (IPR027916, PF14531), OR any protein_kinase gene whose product '
      'name contains "rhoptry" or a ROP name (regex `rhoptry|\\bROP\\d+`). The second arm relies on the product name, which is an '
      'annotation, not a domain.\n'
      '* **srs_sag**: InterPro IPR007226, IPR036755 (SRS domain), IPR028352 (SAG1 family); Pfam PF04092.\n'
      '* **ef_hand**: InterPro IPR002048, IPR011992; Pfam PF00036, PF13499, PF13202, PF13405, PF13833.\n'
      '* **wd40**: InterPro IPR001680, IPR036322, IPR019775, IPR020472; Pfam PF00400.\n'
      '* **ankyrin**: InterPro IPR002110, IPR036770; Pfam PF00023, PF12796, PF13606, PF13637, PF13857.\n'
      '* **rna_binding**: RRM (IPR000504, IPR003954, PF00076), Pumilio/PUF (IPR001313, IPR033133, PF00806) and KH (IPR004087, '
      'IPR004088, IPR036612, PF00013) when present in the data.\n'
      '* **proteasome_ubiquitin**: any InterPro or Pfam entry whose description contains (case-insensitive) proteasome, ubiquitin, '
      'HECT, F-box, cullin, SUMO, NEDD8 or neddylation.\n'
      '* **tm_rich**: DeepTMHMM predicts at least 3 transmembrane helices (`n_tm >= 3`).')
nb.code("RULES = {",
        " 'protein_kinase': dict(ip=['IPR000719', 'IPR001245', 'IPR020635', 'IPR027916'], pf=['PF00069', 'PF07714', 'PF14531']),",
        " 'ap2_apiap2': dict(ip=['IPR001471'], pf=['PF00847']),",
        " 'srs_sag': dict(ip=['IPR007226', 'IPR036755', 'IPR028352'], pf=['PF04092']),",
        " 'ef_hand': dict(ip=['IPR002048', 'IPR011992'], pf=['PF00036', 'PF13499', 'PF13202', 'PF13405', 'PF13833']),",
        " 'wd40': dict(ip=['IPR001680', 'IPR036322', 'IPR019775', 'IPR020472'], pf=['PF00400']),",
        " 'ankyrin': dict(ip=['IPR002110', 'IPR036770'], pf=['PF00023', 'PF12796', 'PF13606', 'PF13637', 'PF13857']),",
        " 'rna_binding': dict(ip=['IPR000504', 'IPR003954', 'IPR001313', 'IPR033133', 'IPR004087', 'IPR004088', 'IPR036612'],",
        "                     pf=['PF00076', 'PF00806', 'PF00013']),",
        "}",
        "UB_RX = re.compile(r'proteasome|ubiquitin|HECT|F-box|cullin|SUMO|NEDD8|neddylation', re.I)",
        "ub_ip = sorted(a for a in all_ip if UB_RX.search(ip_desc.get(a, '')))",
        "ub_pf = sorted(a for a in all_pf if UB_RX.search(pf_desc.get(a, '')))",
        "RULES['proteasome_ubiquitin'] = dict(ip=ub_ip, pf=ub_pf)",
        "hit = lambda sets, accs: sets.map(lambda s: bool(s & set(accs)))",
        "cls = {}",
        "for name, r in RULES.items():",
        "    cls['domclass_' + name] = hit(sets_ip, r['ip']) | hit(sets_pf, r['pf'])",
        "    missing = [a for a in r['ip'] if a not in all_ip] + [a for a in r['pf'] if a not in all_pf]",
        "    if missing: print(name, ': accessions in the rule that no gene carries:', missing)",
        "ROP_RX = re.compile(r'rhoptry|\\bROP\\d+', re.I)",
        "apicomplexa_kinase = hit(sets_ip, ['IPR027916']) | hit(sets_pf, ['PF14531'])",
        "cls['domclass_rhoptry_kinase_ropk'] = apicomplexa_kinase | (cls['domclass_protein_kinase'] & nodes['product'].fillna('').str.contains(ROP_RX))",
        "cls['domclass_tm_rich'] = nodes.n_tm.fillna(0) >= 3",
        "classes = pd.DataFrame(cls, index=nodes.index)",
        "print('proteasome_ubiquitin matches', len(ub_ip), 'InterPro and', len(ub_pf), 'Pfam entries')",
        "classes.sum().rename('genes').to_frame()")

nb.md('### Spot checks on known genes',
      'Four genes whose biology is established: CDPK1 (a calcium-dependent kinase with EF-hands), ROP18 (a secreted rhoptry kinase), '
      'AP2XI-5 (an AP2 transcription factor) and SAG1 (the major surface antigen, ToxoDB name SRS29B). Expected True is what the '
      'biology says; a mismatch is **reported, the rule is not changed**.')
nb.code("known = {'CDPK1': ('TGME49_301440', ['protein_kinase', 'ef_hand']),",
        "         'ROP18': ('TGME49_205250', ['protein_kinase', 'rhoptry_kinase_ropk']),",
        "         'AP2XI-5': ('TGME49_216220', ['ap2_apiap2']),",
        "         'SAG1 (SRS29B)': ('TGME49_233460', ['srs_sag'])}",
        "rows, MISSED = [], []",
        "for label, (g, expect) in known.items():",
        "    for c in expect:",
        "        ok = bool(classes.loc[g, 'domclass_' + c])",
        "        rows.append(dict(gene=label, gene_id=g, product=nodes.loc[g, 'product'], expected_class=c, caught=ok))",
        "        if not ok: MISSED.append(f'{label} ({g}) not caught by {c}')",
        "print('Missed known genes:', MISSED if MISSED else 'none')",
        "for label, (g, _) in known.items(): print(label, '| InterPro', ip[g], '| Pfam', pf[g])",
        "pd.DataFrame(rows)")

nb.md('## 7 · Group 3: counts and the "no domain" flag')
nb.code("counts = pd.DataFrame({'n_interpro_domains': ip.map(len).astype(int), 'n_pfam_domains': pf.map(len).astype(int)}, index=nodes.index)",
        "counts['has_any_domain'] = (counts.n_interpro_domains > 0) | (counts.n_pfam_domains > 0)",
        "counts['domain_novel'] = ~counts.has_any_domain",
        "print('genes with no InterPro and no Pfam domain at all (domain_novel):', int(counts.domain_novel.sum()), 'of', len(counts),",
        "      f'({counts.domain_novel.mean():.1%})')",
        "print('has an InterPro entry but no Pfam:', int(((counts.n_interpro_domains > 0) & (counts.n_pfam_domains == 0)).sum()))",
        "print('has a Pfam entry but no InterPro:', int(((counts.n_interpro_domains == 0) & (counts.n_pfam_domains > 0)).sum()))",
        "print('consistency with the has_domain column of nodes.parquet:', bool((counts.has_any_domain.astype(int) == nodes.has_domain).all()))",
        "counts.describe().T")

nb.md('## 8 · The 30 most common families')
nb.code("top = pd.concat([lookup.sort_values('gene_count', ascending=False).head(30)]).reset_index(drop=True)",
        "top[['accession', 'source', 'description', 'gene_count']]")

nb.md('## 9 · Cross-check with the map: are host-facing genes more often domain-novel?',
      'The host-facing labels (65 True genes in the v0 map) are few, and this is a question that matters to Einar, so the result '
      'is read as a **hint with wide uncertainty**, not a finding. The Fisher test compares the 65 host-facing genes with all other '
      'genes.')
nb.code("v = pd.read_parquet('results/information_map_v0_2026_11_04/values.parquet')",
        "hf = v.host_facing.eq(True).reindex(nodes.index).fillna(False).astype(bool)",
        "neg = v.host_facing.eq(False).reindex(nodes.index).fillna(False).astype(bool)",
        "nov = counts.domain_novel",
        "a, b = int(nov[hf].sum()), int(hf.sum())",
        "c, d = int(nov[~hf].sum()), int((~hf).sum())",
        "odds, p = fisher_exact([[a, b - a], [c, d - c]])",
        "from scipy.stats import beta",
        "ci = (beta.ppf(0.025, a, b - a + 1), beta.ppf(0.975, a + 1, b - a))",
        "print(f'host-facing genes: {b}; domain_novel {a} ({a / b:.1%}; 95% CI {ci[0]:.1%}-{ci[1]:.1%})')",
        "print(f'all genes: {len(nov)}; domain_novel {int(nov.sum())} ({nov.mean():.1%})')",
        "print(f'all non-host-facing genes: {d}; domain_novel {c} ({c / d:.1%})')",
        "print(f'Fisher exact, host-facing vs all other genes: odds ratio {odds:.2f}, two-sided p = {p:.3g}')",
        "n2 = int(neg.sum()); a2 = int(nov[neg].sum())",
        "o2, p2 = fisher_exact([[a, b - a], [a2, n2 - a2]])",
        "print(f'curated NOT host-facing genes: {n2}; domain_novel {a2} ({a2 / n2:.1%}); Fisher vs host-facing: odds {o2:.2f}, p = {p2:.3g}')",
        "HF_RESULT = dict(a=a, b=b, c=c, d=d, odds=odds, p=p, a2=a2, n2=n2, o2=o2, p2=p2)")

nb.md('## 10 · Register the columns in the information map',
      'Every column is tier 2, `normalization="not_applicable"`, `role="measurement"`, kind `measure`. The family is `interpro` or '
      '`pfam` for the family columns; the counts and classes are `interpro`. The registration is done here, not by editing '
      '`parascape/information_map.py`.')
nb.code("from parascape.information_map import Trait, register, TRAITS",
        "traits = pd.concat([fam, classes, counts], axis=1)",
        "desc = {}",
        "for _, r in lookup.iterrows(): desc[r.column] = f\"Has {r.source} domain {r.accession}: {r.description}\"",
        "CLASS_DESC = {'protein_kinase': 'Has a protein kinase domain (InterPro IPR000719/IPR001245/IPR020635/IPR027916 or Pfam PF00069/PF07714/PF14531)',",
        "  'ap2_apiap2': 'Has an AP2/ApiAP2 DNA-binding domain (IPR001471 or PF00847)',",
        "  'rhoptry_kinase_ropk': 'Rhoptry kinase (ROPK): Apicomplexa kinase-like domain, or a kinase whose product name says rhoptry/ROP',",
        "  'srs_sag': 'Has an SRS (SAG-related sequence) surface antigen domain',",
        "  'ef_hand': 'Has an EF-hand calcium-binding domain', 'wd40': 'Has WD40 repeats', 'ankyrin': 'Has ankyrin repeats',",
        "  'rna_binding': 'Has an RNA-binding domain (RRM, Pumilio or KH)',",
        "  'proteasome_ubiquitin': 'Has a proteasome or ubiquitin-system domain (by description words)',",
        "  'tm_rich': 'Predicted to cross the membrane at least 3 times (DeepTMHMM n_tm >= 3)'}",
        "for k, t in CLASS_DESC.items(): desc['domclass_' + k] = t",
        "desc.update(n_interpro_domains='Number of InterPro entries found in the protein', n_pfam_domains='Number of Pfam families found in the protein',",
        "  has_any_domain='Has at least one InterPro or Pfam domain', domain_novel='No InterPro and no Pfam domain at all (structurally unrecognised by sequence)')",
        "lines = []",
        "for c in traits.columns:",
        "    fam_name = 'pfam' if c.startswith('dom_PF') else 'interpro'",
        "    register(Trait(c, 'measure', desc[c], tier=2, family=fam_name, normalization='not_applicable', role='measurement'))",
        "    lines.append(f\"register(Trait({c!r}, 'measure', {desc[c]!r}, tier=2, family={fam_name!r}, normalization='not_applicable', role='measurement'))\")",
        "print(len(traits.columns), 'columns registered;', sum(t.tier == 2 and t.family in ('interpro', 'pfam') for t in TRAITS.values()), 'interpro/pfam tier-2 traits in the registry')",
        "print('\\n'.join(lines[:3] + ['...'] + [l for l in lines if 'domclass_' in l or 'domain_novel' in l or \"'n_\" in l or 'has_any' in l]))",
        "REGISTER_LINES = lines")


README_CODE = '''
rule_lines = []
for k, v in RULES.items():
    if k == 'proteasome_ubiquitin':
        rule_lines.append(f"* **{k}**: {len(v['ip'])} InterPro and {len(v['pf'])} Pfam entries whose description matches the regex below.")
    else:
        rule_lines.append(f"* **{k}**: InterPro {v['ip']}; Pfam {v['pf']}")
cut_lines = [f"* {k}: {int(x['InterPro families'])} InterPro / {int(x['Pfam families'])} Pfam / {int(x['total columns'])} total" for k, x in cut.iterrows()]
class_lines = [f"* {c}: {int(classes[c].sum())} genes" for c in classes.columns]
readme = "\\n".join([
 "# Domain traits (InterPro and Pfam), tier 2 of the information map", "",
 "Built by scripts/build_domain_traits_notebook.py (notebook notebooks/domain_traits_2026_11_04.ipynb). One row per gene, the 8,140 genes of parascape/data/nodes.parquet.", "",
 "## Files",
 f"* domain_traits.parquet: {traits.shape[1]} columns (index gene_id): {fam.shape[1]} family flags dom_<accession>, {classes.shape[1]} named classes domclass_*, 4 counts/flags.",
 "* domain_lookup.tsv: column, accession, description, gene count, source for every family column.",
 "* coverage.tsv: how many genes have what.", "",
 "## Frequency cutoff", f"A family must occur in at least {CUTOFF} genes. Columns at each cutoff:"] + cut_lines + [
 f"Chosen: {CUTOFF}, giving {len(keep_ip)} InterPro + {len(keep_pf)} Pfam = {fam.shape[1]} family columns.", "",
 "## Matching rules for the named classes (verbatim)"] + rule_lines + [
 f"* proteasome_ubiquitin regex (case-insensitive, on the InterPro and Pfam descriptions): {UB_RX.pattern}",
 f"* rhoptry_kinase_ropk: IPR027916 or PF14531 (Apicomplexa protein kinase-like), OR protein_kinase AND product name matches regex {ROP_RX.pattern} (case-insensitive). The second arm uses the product annotation, not a domain.",
 "* tm_rich: DeepTMHMM n_tm >= 3.",
 "* Broad superfamily entries (for example IPR011009 Protein kinase-like superfamily) are deliberately excluded.", "",
 "Gene counts per class:"] + class_lines + ["",
 "Known-gene spot checks (CDPK1, ROP18, AP2XI-5, SAG1), missed: " + (str(MISSED) if MISSED else "none"), "",
 "## Coverage",
 f"{int(counts.has_any_domain.sum())} of {len(nodes)} genes have at least one domain; {int(counts.domain_novel.sum())} ({counts.domain_novel.mean():.1%}) have none (domain_novel). Check against the raw file: {COVERAGE_NOTE}.",
 "Source: ToxoDB gene attributes (datasets/interpro_tgon.csv, registry key interpro; the registry holds no PMID or accession).", "",
 "## Host-facing cross-check",
 f"domain_novel in {r['a']} of {r['b']} host-facing genes ({r['a']/r['b']:.1%}) versus {int(nov.sum())} of {len(nov)} genes overall ({nov.mean():.1%}); Fisher exact, host-facing vs all other genes: odds ratio {r['odds']:.2f}, p = {r['p']:.3g}. Only {r['b']} labelled genes: a hint, not a finding.", "",
 "## Caveats",
 "* An InterPro or Pfam assignment is itself a prediction from the amino-acid sequence, not an experiment. False means not recognised, not proven absent.",
 "* Domain annotation partly records study effort, not only conserved architecture (registry note).",
 "* Missing is not zero: domain_novel genes may carry real domains that no database knows.", ""])
(OUT / 'README.md').write_text(readme)
print(readme)
cov
'''

nb.md('## 11 · Outputs')
nb.code("traits.index.name = 'gene_id'",
        "traits.to_parquet(OUT / 'domain_traits.parquet')",
        "lookup.to_csv(OUT / 'domain_lookup.tsv', sep='\\t', index=False)",
        "cov = [dict(item='genes in universe', genes=len(nodes)),",
        "       dict(item='genes with at least one InterPro entry', genes=int((counts.n_interpro_domains > 0).sum())),",
        "       dict(item='genes with at least one Pfam family', genes=int((counts.n_pfam_domains > 0).sum())),",
        "       dict(item='genes with any domain', genes=int(counts.has_any_domain.sum())),",
        "       dict(item='genes with no domain (domain_novel)', genes=int(counts.domain_novel.sum())),",
        "       dict(item=f'InterPro families at >= {CUTOFF} genes', genes=len(keep_ip)),",
        "       dict(item=f'Pfam families at >= {CUTOFF} genes', genes=len(keep_pf)),",
        "       dict(item='genes with at least one kept family column True', genes=int(fam.any(axis=1).sum())),",
        "       dict(item='InterPro lists differing from raw file', genes=len(d_ip)),",
        "       dict(item='InterPro ID lists at the 20-ID cap', genes=len(at_cap))]",
        "cov += [dict(item='class ' + c, genes=int(classes[c].sum())) for c in classes.columns]",
        "cov = pd.DataFrame(cov); cov.to_csv(OUT / 'coverage.tsv', sep='\\t', index=False)",
        "r = HF_RESULT",
        README_CODE)

nb.write(str(ROOT / 'notebooks/domain_traits_2026_11_04.ipynb'))
print('wrote notebook')
