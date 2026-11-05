"""Instruction 72, v0 build: the information map from OBSERVED facts only.

Writes notebooks/information_map_v0_2026_11_04.ipynb (executed) and results/information_map_v0_2026_11_04/.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('Information map v0: observed facts per gene, with a probability for each')
nb.md('Instruction 72 (Einar, 2026-11-04), v0 build. One row per Toxoplasma ME49 gene (the 8,140 genes of the ParaScape '
      'node table), one column per trait, and two matching tables of the same shape: how probable the fact is, and where '
      'it came from.',
      '**Observed facts only.** Nothing predicted is written here: the TAGM / `lopit_*` columns of the node table are '
      'classifier outputs and are left out. **Missing is empty, never zero, never "no".** A gene nobody tested for a '
      'trait is simply blank; unlabeled is not negative.',
      '**The probabilities are provisional (v0).** They are tier defaults, not measured: direct experiment 0.95, '
      'curated marker 0.90, single high-throughput screen call 0.80, proximity labelling 0.70; StarPath crosslinks use '
      '1 minus their own local false-discovery rate; members of interaction studies get 0.5 until a background frequency-style '
      'background exists. In v1 these defaults are replaced by measured agreement between independent sources.')

nb.md('## 1 · Setup, gene universe and tier defaults',
      'The gene universe is the 8,140 ME49 genes of `parascape/data/nodes.parquet`. The tier defaults are written down '
      'once here so every later step uses the same numbers.')
nb.code(r'''
import json, re, time, hashlib, datetime as dt
from pathlib import Path
import numpy as np, pandas as pd, requests
OUT = Path('results/information_map_v0_2026_11_04'); SRC = OUT / 'sources'
SRC.mkdir(parents=True, exist_ok=True)
nodes = pd.read_parquet('parascape/data/nodes.parquet')
GENES = pd.Index(nodes['gene_id'].astype(str), name='gene_id')
assert GENES.is_unique
GENESET = set(GENES)
TIER_P = {'observed-direct': 0.95, 'observed-curated': 0.90, 'observed-screen': 0.80, 'observed-proximity': 0.70}
STUDY_MEMBER_P = 0.5
EV = []   # long-form evidence rows
def add_ev(gene, trait, value, p, tier, ref, dataset):
    EV.append((gene, trait, value, float(p), tier, ref, dataset))
print(len(GENES), 'genes;', 'tier defaults', TIER_P)
''')

nb.md('## 2 · QuickGO experimental annotations for Toxoplasma (taxon 5811 and its strains)',
      'GO annotations with an experimental evidence code (ECO:0000269 and relatives: direct assay, immunofluorescence, '
      'physical interaction, mutant phenotype and so on) are the closest thing to a published, citable fact per protein. '
      'We ask QuickGO for four sets: host-cell cellular-component terms (for the host-facing trait, including their '
      'child terms), and all cellular-component, biological-process and molecular-function annotations. Responses are '
      'fetched page by page, politely (half a second between calls), and the raw JSON of every page is saved with the '
      'retrieval time so the step can be audited. Annotations with a `NOT` qualifier assert absence of a term and are '
      'set aside, not counted as positive facts.')
nb.code(r'''
BASE = 'https://www.ebi.ac.uk/QuickGO/services/annotation/search'
ECO = 'ECO:0000269,ECO:0000314,ECO:0000353,ECO:0000315,ECO:0000316,ECO:0000270'
HOST_GO = 'GO:0033643,GO:0020005,GO:0030430,GO:0042025,GO:0020002,GO:0043657'
sess = requests.Session(); sess.headers['User-Agent'] = 'starplast-information-map-v0 (research; einarolafssonbirnir@gmail.com)'
def get(url, params, tries=4):
    for k in range(tries):
        r = sess.get(url, params=params, timeout=60)
        if r.status_code == 200:
            return r
        time.sleep(2 * (k + 1))
    r.raise_for_status()
def quickgo(name, extra):
    params = dict(taxonId=5811, taxonUsage='descendants', evidenceCode=ECO, evidenceCodeUsage='descendants', limit=100, **extra)
    rows, page, stamps = [], 1, []
    while True:
        params['page'] = page
        r = get(BASE, params)
        d = r.json()
        stamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')
        (SRC / f'quickgo_{name}_page{page:03d}.json').write_text(json.dumps({'retrieved_utc': stamp, 'url': r.url, 'response': d}))
        rows += d['results']; stamps.append(stamp)
        if page >= d['pageInfo']['total']:
            break
        page += 1; time.sleep(0.5)
    print(f'{name}: {len(rows)} annotations in {page} pages (hits reported {d["numberOfHits"]}), retrieved {stamps[0]}')
    return rows
qg = {
    'hostcell_cc': quickgo('hostcell_cc', dict(goId=HOST_GO, goUsage='descendants')),
    'cc': quickgo('cc', dict(aspect='cellular_component')),
    'bp': quickgo('bp', dict(aspect='biological_process')),
    'mf': quickgo('mf', dict(aspect='molecular_function')),
}
''')

nb.md('## 3 · Map UniProt accessions to ME49 gene IDs',
      'QuickGO identifies proteins by UniProt accession, often from a strain other than ME49 (RH, GT1, VEG). UniProt '
      'lists cross-references to every Toxoplasma genome build at ToxoDB, so we read the `TGME49_` identifiers from '
      'those and keep only those that exist in our 8,140-gene universe. An accession that maps to more than one '
      'universe gene is attached to each (and listed as ambiguous); an accession with none is written to '
      '`sources/unmapped_uniprot.tsv` with the traits it carried, never silently dropped.')
nb.code(r'''
accs = sorted({r['geneProductId'].split(':', 1)[1] for rows in qg.values() for r in rows if r['geneProductId'].startswith('UniProtKB:')})
other = sorted({r['geneProductId'] for rows in qg.values() for r in rows if not r['geneProductId'].startswith('UniProtKB:')})
print(len(accs), 'UniProt accessions;', len(other), 'non-UniProt product ids', other[:5])
acc2genes = {}
for i in range(0, len(accs), 40):
    chunk = accs[i:i + 40]
    q = ' OR '.join(f'accession:{a}' for a in chunk)
    r = get('https://rest.uniprot.org/uniprotkb/search', dict(query=q, fields='accession,xref_veupathdb', format='tsv', size=500))
    for line in r.text.splitlines()[1:]:
        a, _, x = line.partition('\t')
        acc2genes[a] = sorted({t.split(':', 1)[1] for t in x.split(';') if t.startswith('ToxoDB:TGME49_')} & GENESET)
    time.sleep(0.5)
mapped = {a: g for a, g in acc2genes.items() if g}
ambiguous = {a: g for a, g in mapped.items() if len(g) > 1}
unmapped = sorted(set(accs) - set(mapped))
print(len(mapped), 'mapped;', len(ambiguous), 'ambiguous (several universe genes);', len(unmapped), 'unmapped')
''')
nb.code(r'''
rows = []
for name, lst in qg.items():
    for r in lst:
        a = r['geneProductId'].split(':', 1)[1]
        if a in unmapped or not r['geneProductId'].startswith('UniProtKB:'):
            rows.append((a, name, r['goId'], r['reference'], r['symbol']))
um = pd.DataFrame(rows, columns=['uniprot', 'query_set', 'go_id', 'reference', 'symbol'])
um.to_csv(SRC / 'unmapped_uniprot.tsv', sep='\t', index=False)
pd.Series({'unmapped accessions': um.uniprot.nunique(), 'annotation rows lost to unmapped': len(um)})
''')

nb.md('## 4 · GO term names',
      'QuickGO returns term identifiers without names, so the names are fetched from its ontology service in batches '
      'and used to make the matrix readable (a top term shows as `GO:0020002 host cell plasma membrane`).')
nb.code(r'''
all_go = sorted({r['goId'] for lst in qg.values() for r in lst})
go_name = {}
for i in range(0, len(all_go), 100):
    ids = ','.join(all_go[i:i + 100])
    d = get(f'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{ids}', {}).json()
    for t in d['results']:
        go_name[t['id']] = t['name']
    time.sleep(0.5)
(SRC / 'go_term_names.json').write_text(json.dumps(go_name, indent=0))
print(len(go_name), 'of', len(all_go), 'GO terms named')
''')

nb.md('## 5 · Turn GO annotations into evidence rows',
      'Each (gene, GO term, reference) becomes one evidence row, tier `observed-direct`, probability 0.95 (provisional). '
      'Host-cell terms feed the **host_facing** trait (a protein with experimental evidence of being in a host-cell '
      'compartment is exported there); cellular-component terms feed **compartment_go_cc**; biological-process and '
      'molecular-function terms feed **function_go_bp** and **function_go_mf**. `NOT` annotations are counted and skipped.')
nb.code(r'''
TRAIT_OF = {'hostcell_cc': 'host_facing', 'cc': 'compartment_go_cc', 'bp': 'function_go_bp', 'mf': 'function_go_mf'}
not_skipped = {}
for name, lst in qg.items():
    n_not = 0
    for r in lst:
        if str(r.get('qualifier', '')).upper().startswith('NOT'):
            n_not += 1; continue
        a = r['geneProductId'].split(':', 1)[1]
        for g in mapped.get(a, []):
            ds = f'QuickGO/GOA experimental ({r["assignedBy"]}; {r["evidenceCode"]}; query {name})'
            if name == 'hostcell_cc':
                add_ev(g, 'host_facing', 'yes', TIER_P['observed-direct'], 'observed-direct', r['reference'], ds + f'; {r["goId"]}')
            else:
                add_ev(g, TRAIT_OF[name], f'{r["goId"]} {go_name.get(r["goId"], "")}'.strip(), TIER_P['observed-direct'], 'observed-direct', r['reference'], ds)
    not_skipped[name] = n_not
print('NOT-qualified annotations skipped:', not_skipped)
pd.DataFrame(EV, columns=['gene_id','trait','value','probability','source_tier','reference','dataset']).groupby('trait').gene_id.nunique()
''')

nb.md('## 6 · hyperLOPIT markers (compartment, curated)',
      'The 718 marker proteins of Barylyuk et al. 2020 (Cell Host and Microbe) were chosen from prior knowledge and '
      'microscopy, then used to train the spatial classifier. They are the curation-grade compartment truth (tier '
      '`observed-curated`, 0.90) -- not independent experiments, which is why they sit below direct GO evidence. Only '
      'the `markers` column is used; the `tagm.map.*` columns are classifier predictions and are ignored. The file '
      'is checked against the pinned SHA-256 before it is read.')
nb.code(r'''
MMC5 = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/sources/barylyuk2020_mmc5/mmc5.xls')
sha = hashlib.sha256(MMC5.read_bytes()).hexdigest()
assert sha == 'a58a8778714c42e344099de9429daa3be4e9c70021b0aeab3de65578c165b8f2', sha
mk = pd.read_excel(MMC5, sheet_name='S4a - TAGM-MAP classification')[['Accession', 'markers']]
mk = mk[mk.markers.notna() & (mk.markers != 'unknown')]
mk_in = mk[mk.Accession.isin(GENESET)]
mk_out = mk[~mk.Accession.isin(GENESET)]
for g, m in zip(mk_in.Accession, mk_in.markers):
    add_ev(g, 'compartment_hyperlopit', m, TIER_P['observed-curated'], 'observed-curated', 'PMID:32619502', 'Barylyuk 2020 mmc5.xls S4a markers')
mk_out.to_csv(SRC / 'hyperlopit_markers_not_in_universe.tsv', sep='\t', index=False)
print(len(mk), 'markers;', len(mk_in), 'in universe;', len(mk_out), 'not in the 8,140-gene universe:', list(mk_out.Accession))
mk_in.markers.value_counts().head(8)
''')

nb.md('## 7 · Curated host interactions and PVM proximity (host-facing)',
      'Two more host-facing sources. (a) The 18 curated parasite-to-host interactions (ROP16-STAT3 and so on): a '
      'protein that acts on a host protein is host-facing; tier `observed-curated` 0.90, reference = the paper in the '
      'table. (b) PVM proximity labelling (PMID 34749525): the 73 proteins the study called positive are host-facing '
      'candidates at the vacuole membrane, tier `observed-proximity` 0.70. The same study also lists 1,201 "likely '
      'negatives", but not being near the vacuole membrane is **not** evidence of not being exported, so those stay '
      'empty. Provenance of the proximity column is read from `parascape/datasets.py`.')
nb.code(r'''
from parascape import datasets as pds
prov = pds.provenance('pvm_proximity_positive', 'Tg')
print('pvm_proximity_positive ->', prov.key, '| PMID', prov.pmid, '|', prov.accession)
hi = pd.read_parquet('parascape/data/host_interactions.parquet')
for _, r in hi.iterrows():
    if r.gene_id in GENESET:
        add_ev(r.gene_id, 'host_facing', 'yes', TIER_P['observed-curated'], 'observed-curated', r.key_reference,
               f'host_interactions.parquet ({r.parasite_protein} -> {r.host_target})')
pv = nodes.loc[nodes.pvm_proximity_positive == 1, 'gene_id']
for g in pv:
    add_ev(g, 'host_facing', 'yes', TIER_P['observed-proximity'], 'observed-proximity', f'PMID:{prov.pmid}', f'{prov.key} ({prov.accession}); pvm_proximity_positive')
print(len(hi), 'curated interactions;', hi.gene_id.nunique(), 'distinct parasite genes;', len(pv), 'PVM proximity positives')
''')

nb.md('## 7b · Proteins that are definitely NOT host-facing',
      'Einar, 2026-11-04: "you also know which proteins are definitely not host facing, so that information should be '
      'there as well in this column." Without real negatives, precision cannot be measured honestly, so `host_facing` '
      'is three-valued: True, **False**, or empty for untested.',
      '* **Curated gold negatives** (16 proteins, `hostfacing_goldsets.json`): each was localised and found to stay '
      'inside the parasite or the vacuole -- GRA9 by immunogold electron microscopy in the dense granules and then the '
      'intravacuolar network, GRA2 shaping the network from the vacuole side, and so on. Tier `observed-curated`, 0.90.',
      '* **hyperLOPIT markers of compartments physically inside the parasite**: apicoplast, mitochondrion, nucleus, '
      'nucleolus, ribosome, proteasome, Golgi, ER, cytosol, tubulin cytoskeleton and IMC. Probability **0.80**, below '
      'the curated tier, because this is an inference from a resident marker rather than a direct export test.',
      '* **Deliberately not negatives:** dense granules, rhoptries, micronemes, apical, plasma membrane and '
      'endomembrane vesicles. These are the secretory and boundary compartments whose proteins can be exported, so a '
      'marker there says nothing either way. The PVM-proximity "likely negatives" also stay empty.')
nb.code(r'''
GOLD = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/hostfacing_goldsets.json')
gold = json.loads(GOLD.read_text())
n_gold = 0
for r in gold['gold_negative']:
    if r['gene_id'] in GENESET:
        add_ev(r['gene_id'], 'host_facing', 'no', TIER_P['observed-curated'], 'observed-curated',
               r['evidence'][:200], f"hostfacing_goldsets.json gold_negative ({r['name']}; confidence {r['confidence']})")
        n_gold += 1
INSIDE = ('apicoplast', 'mitochondrion', 'nucleus', 'nucleolus', 'ribosome', 'proteasome', 'Golgi', 'ER',
          'cytosol', 'tubulin cytoskeleton', 'IMC')
NOT_NEGATIVE = ('dense granules', 'rhoptries', 'micronemes', 'apical', 'PM -', 'endomembrane vesicles')
def inside_parasite(niche):
    if any(k in niche for k in NOT_NEGATIVE):
        return False
    return any(niche.startswith(k) or k in niche for k in INSIDE)
inside = {g: m for g, m in zip(mk_in.Accession, mk_in.markers) if inside_parasite(m)}
for g, m in inside.items():
    add_ev(g, 'host_facing', 'no', 0.80, 'observed-curated', 'PMID:32619502',
           f'Barylyuk 2020 hyperLOPIT marker of {m} (a compartment inside the parasite; inference from residence)')
print(n_gold, 'curated gold negatives;', len(inside), 'internal-compartment markers')
pd.Series(sorted({m for m in inside.values()})).to_frame('niches counted as inside the parasite').T
''')

nb.md('## 7c · Einar\'s review of the 140 literature claims (2026-11-04)',
      'Einar reviewed every published claim behind the host-facing candidates: 107 claims host-facing, 33 not, '
      'over 98 genes. Collapsed per gene by the rule he set -- strongest evidence wins, with three explicit '
      'overrides -- this gives **66 host-facing and 32 not**. A gene decided by a direct experimental claim is '
      '`observed-direct` 0.95; one decided by a proximity or inferred claim is `observed-curated` 0.90, because '
      'his judgement is worth more than the raw screen call and less than a published localisation.',
      'This supersedes the raw proximity calls wherever he ruled: he **rejected 22 dense granule proteins** that '
      'the PVM proximity screen had called positive (GRA1, GRA4, GRA22, GRA30, GRA36, GRA37, GRA40, GRA43, GRA46, '
      'GRA48, GRA49, GRA51-GRA59, GRA12B, ROP34). Those proximity rows stay in the evidence table as the record '
      'but lose on probability, so the matrix shows "not host-facing". This is the distinction the first model '
      'could not learn: in the vacuole versus exported into the host.')
nb.code(r'''
REVIEW = Path('results/host_facing_truth_2026_11_04/evidence.parquet')
review = pd.read_parquet(REVIEW)
for r in review.itertuples():
    add_ev(r.gene_id, 'host_facing', r.value, r.probability, r.source_tier, r.reference, r.dataset)
print(len(review), 'reviewed genes;', (review.value == 'yes').sum(), 'host-facing,', (review.value == 'no').sum(), 'not')
review.value.value_counts().to_frame('genes').T
''')

nb.md('## 7d · Function classes from the reviewed truth set',
      'The map held experimental GO function for 24-27 genes. The function truth set (2026-11-04) adds a coarse '
      'class for **242 genes**: 32 from experimental GO, 113 from reviewed Swiss-Prot curator comments and 97 '
      'from the imaging screen, with 57 % of genes carrying more than one class. It enters as a `terms` trait, '
      'so the matrix shows the best-supported class and `function_class_n` counts how many a gene has, while '
      'every class stays in the evidence table.',
      'The probability is the tier of the claim that decided the gene: 0.95 for experimental GO, 0.90 for a '
      'curator comment, 0.80 for a screen phenotype. A screen phenotype supports a *process* claim only, never '
      'a molecular function, which is why those genes sit at the lowest tier.')
nb.code(r'''
FUNC = Path('results/function_truth_2026_11_04/evidence.parquet')
func = pd.read_parquet(FUNC)
for r in func.itertuples():
    for cls in str(r.value).split('|'):          # multi-label genes carry several classes
        if cls.strip():
            add_ev(r.gene_id, 'function_class', cls.strip(), r.probability, r.source_tier, r.reference, r.dataset)
print(len(func), 'genes;', sum(len(str(r.value).split('|')) for r in func.itertuples()), 'class assignments')
func.source_tier.value_counts().to_frame('genes').T
''')

nb.md('## 8 · Screen phenotypes (function, observed-screen)',
      'The arrayed splitCas9 imaging screen (PMID 35538310) scored 319 genes for egress, actin, apicoplast and '
      'replication phenotypes. Here a measured "no phenotype" is a real observation of that screen (the gene was '
      'tested), so both True and False are recorded, tier `observed-screen` 0.80. Genes that were not screened stay '
      'empty. Provenance of every `screen_` column is read from `datasets.py`.')
nb.code(r'''
screen_cols = [c for c in nodes.columns if c.startswith('screen_')]
for c in screen_cols:
    p = pds.provenance(c, 'Tg')
    print(c, '->', p.key, 'PMID', p.pmid)
    sub = nodes.loc[nodes[c].notna(), ['gene_id', c]]
    for g, v in zip(sub.gene_id, sub[c]):
        add_ev(g, c, bool(v), TIER_P['observed-screen'], 'observed-screen', f'PMID:{p.pmid}', f'{p.key} ({p.accession}); {c}')
''')

nb.md('## 9 · Binding partners (pairs table)',
      'Pair-valued facts do not fit one cell per gene, so they live in `pairs.parquet`: gene_a, gene_b, which side the '
      'partner is on (parasite or host), probability, tier and reference.',
      '* **StarPath crosslinks** (cross-linking mass spectrometry): probability = 1 minus the pair local FDR, tier '
      '`observed-direct`. Reference is the StarPath preprint (PMC12505969).',
      '* **Interaction studies** (BioID / pull-down tables parsed from papers): pairs are bait to every listed member, '
      'for the 31 studies where the bait gene is known. Tier `observed-screen`, probability **0.5, provisional** until a '
      'background frequency-style background says which members are sticky contaminants.',
      '* **Curated host interactions**: parasite gene to the named host protein, side `host`, tier `observed-curated` 0.90.')
nb.code(r'''
cl = pd.read_parquet('parascape/data/crosslink_models.parquet')
cl = cl[cl.gene_a.isin(GENESET) & cl.gene_b.isin(GENESET) & (cl.gene_a != cl.gene_b)]
p1 = pd.DataFrame({'gene_a': cl.gene_a, 'gene_b': cl.gene_b, 'partner_side': 'parasite',
                   'probability': (1 - cl.local_fdr).clip(0, 1), 'source_tier': 'observed-direct',
                   'reference': 'PMC12505969 (StarPath XL-MS); starpath_id ' + cl.starpath_id.astype(str)})
st = pd.read_parquet('parascape/data/interaction_studies.parquet')
mem = pd.read_parquet('parascape/data/interaction_study_members.parquet')
baits = st.loc[st.bait_gene.notna(), ['pmid', 'bait_gene', 'method']]
j = mem.merge(baits[['pmid', 'bait_gene']], on='pmid')
j = j[j.bait_gene.isin(GENESET) & j.gene_id.isin(GENESET) & (j.bait_gene != j.gene_id)]
p2 = pd.DataFrame({'gene_a': j.bait_gene, 'gene_b': j.gene_id, 'partner_side': 'parasite',
                   'probability': STUDY_MEMBER_P, 'source_tier': 'observed-screen',
                   'reference': 'PMID:' + j.pmid + ' (' + j.method + '; provisional 0.5 pending background table)'})
p3 = pd.DataFrame({'gene_a': hi.gene_id, 'gene_b': hi.host_target, 'partner_side': 'host',
                   'probability': TIER_P['observed-curated'], 'source_tier': 'observed-curated',
                   'reference': hi.key_reference})
pairs = pd.concat([p1, p2, p3], ignore_index=True)
pairs = pairs.sort_values(['source_tier', 'probability'], ascending=[True, False]).reset_index(drop=True)
print(len(st), 'studies,', st.bait_gene.notna().sum(), 'with a known bait')
pairs.groupby(['source_tier', 'partner_side']).agg(pairs=('gene_a', 'size'), min_p=('probability', 'min'), max_p=('probability', 'max'))
''')

nb.md('## 10 · Assemble the three aligned tables',
      'From the evidence rows we build the matrix. Rules, so nothing is hidden:',
      '* **host_facing**: `yes` where any source says so; its probability is the highest among the sources and its '
      'source is that source tier. All sources stay in the evidence table (conflicts and agreement are visible there).',
      '* **compartment_hyperlopit**: the marker niche. **compartment_go_cc**: the most frequently annotated '
      'experimental CC term (ties by GO id) plus `compartment_go_cc_n`, the number of distinct terms. These two are '
      'kept as separate columns on purpose, so disagreement between sources is never averaged away.',
      '* **function_go_bp / function_go_mf**: top term the same way, plus counts.',
      '* **partner_***: best observed partner (highest probability, ties by id), its probability, the number of '
      'partners with probability at least 0.95 and at least 0.5. Host partners are shown as `host:<name>`.',
      '* **screen_***: the screen result (True/False), blank if not screened.',
      'Every cell of `probability` and `source` is empty exactly where `values` is empty.',
      'The rules live in `parascape/information_map.py`: each column is one registered trait (name, kind, description) '
      'and the matrix is assembled from the evidence rows, so a new column is one `register` line plus its evidence.')
nb.code(r'''
from parascape import information_map as im
im.register(im.Trait('function_class', 'terms', 'Coarse function class from the reviewed truth set',
                     family='function_truth', depends_on=('splitcas9_imaging_screen',)))
for c in screen_cols:                               # one registry line per new column; nothing else changes
    im.register(im.Trait(c, 'measure', f'Arrayed screen result ({c}); a measured no-phenotype is False'))
ev = im.evidence_frame(pd.DataFrame(EV, columns=list(im.EVIDENCE_COLUMNS)))
order = ['host_facing', 'compartment_hyperlopit', 'compartment_go_cc', 'partners', 'function_class',
         'function_go_bp', 'function_go_mf'] + screen_cols
values, prob, srcs = im.build(ev, GENES, pairs=pairs, order=order)
cols = list(values.columns)
values.notna().sum().rename('genes with a value').to_frame()
''')

nb.md('## 11 · Coverage per trait and tier',
      'How many genes have at least one observed fact, per trait and source tier. A gene can appear under several '
      'tiers of the same trait, so the "any tier" row is the union, not the sum. Partner counts come from the pairs '
      'table (a gene counts if it is on either side of a parasite pair, or is the parasite side of a host pair).')
nb.code(r'''
rows = []
def cov(trait, genes_by_tier):
    for t, g in genes_by_tier.items():
        rows.append((trait, t, len(g)))
    rows.append((trait, 'any tier', len(set().union(*genes_by_tier.values())) if genes_by_tier else 0))
for trait in ['host_facing', 'compartment_hyperlopit', 'compartment_go_cc', 'function_go_bp', 'function_go_mf'] + screen_cols:
    e = ev[ev.trait == trait]
    cov(trait, {t: set(g.gene_id) for t, g in e.groupby('source_tier')})
cov('partners', {t: set(g.gene_a) | set(g.gene_b[g.partner_side == 'parasite']) for t, g in pairs.groupby('source_tier')})
coverage = pd.DataFrame(rows, columns=['trait', 'source_tier', 'genes_with_observed_fact'])
coverage.to_csv(OUT / 'coverage.tsv', sep='\t', index=False)
coverage
''')
nb.code('pairs.groupby("source_tier").size().rename("pairs").to_frame()')

nb.md('## 11b · Where sources disagree',
      'A gene with both a positive and a negative host-facing source is shown with the stronger source, and both are '
      'kept. These need a human decision, so they are written to `conflicts.tsv`. Every one of them is a PVM-proximity '
      'positive (0.70) contradicted by a curated negative (0.80-0.90): the proximity screen labelled the host-cytosolic '
      'face of the vacuole membrane, while the curated call places the protein on the vacuole side or in the '
      'intravacuolar network. Which of the two counts as "host-facing" is a definition question for Einar, not '
      'something to average away.')
nb.code(r'''
conf = im.conflicts(ev)
conf.to_csv(OUT / 'conflicts.tsv', sep='\t', index=False)
print(conf.gene_id.nunique(), 'genes with sources that disagree')
conf[conf.trait == 'host_facing'][['gene_id', 'value', 'probability', 'dataset']].sort_values(['gene_id', 'probability'], ascending=[True, False])
''')

nb.md('## 12 · Write the outputs',
      'Parquet files for the three aligned tables, the long evidence table and the pairs table, plus the README that '
      'states what v0 is and is not. Every file has the same gene index; empty means no observed fact.')
nb.code(r'''
values.to_parquet(OUT / 'values.parquet'); prob.to_parquet(OUT / 'probability.parquet'); srcs.to_parquet(OUT / 'source.parquet')
ev.to_parquet(OUT / 'evidence.parquet', index=False); pairs.to_parquet(OUT / 'pairs.parquet', index=False)
cc = coverage.pivot(index='trait', columns='source_tier', values='genes_with_observed_fact').fillna(0).astype(int)
readme = f"""# Information map v0 (2026-11-04)

Instruction 72, v0 build. Observed facts only; no predictions, no gap filling. Built by
`scripts/build_information_map_v0_notebook.py`, executed notebook `notebooks/information_map_v0_2026_11_04.ipynb`.

## Files
* `values.parquet`, `probability.parquet`, `source.parquet` -- {len(GENES)} ME49 genes (index `gene_id`, the rows of
  `parascape/data/nodes.parquet`) x {len(cols)} trait columns, identical shape. Empty (NaN) = no observed fact. Never zero, never "no".
  `host_facing` is only ever `True`; screen columns hold both True and False because a screened gene with no phenotype is a measurement.
* `evidence.parquet` -- long form, one row per fact: gene_id, trait, value, probability, source_tier, reference, dataset ({len(ev)} rows).
* `pairs.parquet` -- binding partners: gene_a (parasite), gene_b (parasite gene, or host protein name when partner_side = host), partner_side, probability, source_tier, reference ({len(pairs)} rows).
* `coverage.tsv` -- genes with an observed fact per trait and tier.
* `sources/` -- raw QuickGO pages with retrieval time, GO term names, unmapped UniProt accessions, hyperLOPIT markers outside the gene universe.

## Provisional probabilities (v0)
Tier defaults, NOT measured: observed-direct 0.95, observed-curated 0.90, observed-screen 0.80, observed-proximity 0.70.
StarPath crosslinks: 1 - local_fdr. Interaction-study members: 0.5, provisional until a background table exists.
v1 replaces these by measured agreement with independent sources.

## Rules followed
* TAGM / `lopit_*` node columns are classifier predictions and are excluded; only the 718 hyperLOPIT markers (curated) are used.
* PVM proximity "likely negatives" are NOT written as `host_facing = no` (not near the PVM is not proof of not exported).
* GO annotations with a NOT qualifier are skipped. Experimental evidence codes: ECO:0000269, 314, 353, 315, 316, 270 (with descendants).
* When several sources speak about one gene, the matrix cell takes the highest-probability one; all remain in `evidence.parquet`.
  Hyperlopit and GO CC are separate columns so disagreement stays visible.
* UniProt to ME49 mapping via UniProt VEuPathDB cross-references, restricted to the 8,140-gene universe; ambiguous (>1 gene): {len(ambiguous)}; unmapped accessions: {len(unmapped)} (listed in `sources/unmapped_uniprot.tsv`).
* The literature host-facing candidate list is not included; it enters after Einar reviews it.

## Coverage (genes with an observed fact)
```
{cc.to_string()}
```
"""
(OUT / 'README.md').write_text(readme)
print(sorted(p.name for p in OUT.iterdir()))
''')

nb.write('notebooks/information_map_v0_2026_11_04.ipynb')
print('wrote notebook')
