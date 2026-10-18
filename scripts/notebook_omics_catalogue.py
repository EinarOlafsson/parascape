"""Write notebooks/omics_catalogue_2026_10_17.ipynb: every Toxoplasma omics study in GEO, PRIDE and the
Metabolomics Workbench, classified by what was perturbed and what was read out.

Instruction 69. Every accession comes from the repositories' own APIs; every query is in the notebook.

    python scripts/notebook_omics_catalogue.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from notebook_runner import ExecutedNotebook  # noqa: E402

OUT = os.path.join(ROOT, "results", "omics_catalogue_2026_10_17")

nb = ExecutedNotebook("Toxoplasma omics catalogue")
nb.md("# Toxoplasma omics studies: what was perturbed, what was measured",
      "Instruction 69. The request (2026-10-17): compare expression between Toxoplasma strains and under "
      "perturbations such as a thioredoxin knockout; as many gene-knockout and drug studies with RNA-seq "
      "as possible; and other omics readouts. This notebook lists every candidate study from the "
      "repositories' own APIs -- no accession typed by hand -- and classifies each, so loading can be "
      "prioritised by what is actually there.",
      "Sources: NCBI GEO via E-utilities (sequencing and arrays), PRIDE (proteomics), Metabolomics "
      "Workbench (metabolomics). Classification is by keyword on title and summary, and is a first pass "
      "for a person to check: the table keeps the text it matched on.")
nb.code("import json, re, time, os, urllib.parse, urllib.request",
        "import pandas as pd",
        f"OUT = {OUT!r}",
        "os.makedirs(OUT, exist_ok=True)",
        "def get(url, tries=4):",
        "    for i in range(tries):",
        "        try:",
        "            with urllib.request.urlopen(url, timeout=60) as r:",
        "                return r.read().decode()",
        "        except Exception as e:",
        "            if i == tries - 1: raise",
        "            time.sleep(2 * (i + 1))",
        "QUERIES = []  # every query sent, kept for the record",
        "def logged(url):",
        "    QUERIES.append(url); return get(url)")
nb.md("## 1. GEO series",
      "Two searches, unioned: the organism field (series whose samples are Toxoplasma), and all fields "
      "(also series indexed under the host only, where the parasite is the perturbation).")
nb.code("E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'",
        "terms = ['\"Toxoplasma gondii\"[Organism] AND gse[ETYP]', 'Toxoplasma[All Fields] AND gse[ETYP]']",
        "ids = {}",
        "for t in terms:",
        "    q = f'{E}/esearch.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'term': t, 'retmax': 5000, 'retmode': 'json'})",
        "    r = json.loads(logged(q))['esearchresult']",
        "    ids[t] = set(r['idlist']); time.sleep(0.4)",
        "    print(f\"{r['count']:>5} {t}\")",
        "uids = sorted(set().union(*ids.values()))",
        "len(uids)")
nb.code("docs = []",
        "for i in range(0, len(uids), 200):",
        "    q = f'{E}/esummary.fcgi?' + urllib.parse.urlencode({'db': 'gds', 'id': ','.join(uids[i:i + 200]), 'retmode': 'json'})",
        "    r = json.loads(logged(q))['result']",
        "    docs += [r[u] for u in r['uids']]; time.sleep(0.4)",
        "geo = pd.DataFrame([{",
        "    'source': 'GEO', 'accession': d.get('accession'), 'title': d.get('title', ''),",
        "    'summary': d.get('summary', ''), 'type': d.get('gdstype', ''), 'samples': int(d.get('n_samples') or 0),",
        "    'taxa': d.get('taxon', ''), 'pubmed': ';'.join(map(str, d.get('pubmedids', []))),",
        "    'supplementary': d.get('suppfile', ''), 'date': d.get('pdat', ''),",
        "    'parasite_measured': 'Toxoplasma' in d.get('taxon', ''),",
        "} for d in docs])",
        "geo['parasite_measured'].value_counts()")
nb.md("## 1b. SRA-only RNA-seq (BioProject)",
      "Some studies deposit reads in SRA without a GEO series -- the Trx1 thioredoxin mutant is one "
      "(PRJNA1008424). Transcriptome projects on the organism, minus those already linked to a GEO series.")
nb.code("Q = chr(34)",
        "# Two searches, unioned: the organism field, and all fields -- the Trx1 thioredoxin mutant project",
        "# (PRJNA1008424) has an EMPTY organism field and is found only by the second.",
        "terms_bp = [f'{Q}Toxoplasma gondii{Q}[Organism] AND {Q}transcriptome or gene expression{Q}[Project Data Type]',",
        "            f'Toxoplasma[All Fields] AND ({Q}transcriptome or gene expression{Q}[Project Data Type] OR {Q}raw sequence reads{Q}[Project Data Type]) AND (RNA-seq[All Fields] OR transcriptome[All Fields] OR RNA sequencing[All Fields])']",
        "bp_ids = set()",
        "for t in terms_bp:",
        "    q = f'{E}/esearch.fcgi?' + urllib.parse.urlencode({'db': 'bioproject', 'term': t, 'retmax': 5000, 'retmode': 'json'})",
        "    r = json.loads(logged(q))['esearchresult']; bp_ids |= set(r['idlist']); time.sleep(0.4)",
        "    print(f\"{r['count']:>5} {t}\")",
        "bp_ids = sorted(bp_ids)",
        "bp = []",
        "for i in range(0, len(bp_ids), 200):",
        "    q = f'{E}/esummary.fcgi?' + urllib.parse.urlencode({'db': 'bioproject', 'id': ','.join(bp_ids[i:i + 200]), 'retmode': 'json'})",
        "    r = json.loads(logged(q))['result']; bp += [r[u] for u in r['uids']]; time.sleep(0.4)",
        "# A GEO series names its BioProject in the summary record's text only sometimes; match on the",
        "# GSE accession through elink instead, which is exact.",
        "linked = set()",
        "for i in range(0, len(bp_ids), 100):",
        "    # One id= per project, so each comes back in its own linkset; a comma-joined id list returns",
        "    # ONE linkset for all of them, which once marked all 125 projects as linked.",
        "    q = f'{E}/elink.fcgi?' + urllib.parse.urlencode([('dbfrom', 'bioproject'), ('db', 'gds'), ('retmode', 'json')] + [('id', x) for x in bp_ids[i:i + 100]])",
        "    for ls in json.loads(logged(q)).get('linksets', []):",
        "        if any(db.get('links') for db in ls.get('linksetdbs', [])): linked.update(map(str, ls.get('ids', [])))",
        "    time.sleep(0.4)",
        "sra_only = pd.DataFrame([{",
        "    'source': 'SRA', 'accession': d.get('project_acc'), 'title': d.get('project_title', ''),",
        "    'summary': d.get('project_description', ''), 'type': 'expression profiling by high throughput sequencing (SRA only)',",
        "    'samples': None, 'taxa': d.get('organism_name', '') or 'Toxoplasma gondii', 'pubmed': '',",
        "    'supplementary': 'raw reads only', 'date': d.get('registration_date', ''), 'parasite_measured': True,",
        "} for d in bp if str(d.get('uid')) not in linked])",
        "len(bp), len(linked), len(sra_only)")
nb.md("## 2. PRIDE (proteomics)")
nb.code("rows, page = [], 0",
        "while True:",
        "    q = 'https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?' + urllib.parse.urlencode({'keyword': 'Toxoplasma', 'pageSize': 100, 'page': page})",
        "    batch = json.loads(logged(q))",
        "    if not batch: break",
        "    rows += batch; page += 1; time.sleep(0.3)",
        "    if len(batch) < 100: break",
        "pride = pd.DataFrame([{",
        "    'source': 'PRIDE', 'accession': p.get('accession'), 'title': p.get('title', ''),",
        "    'summary': p.get('projectDescription') or '',  # not the protocol: 'treated with trypsin' is not a drug",
        "    'protocol': p.get('sampleProcessingProtocol') or '',",
        "    'type': 'proteomics: ' + ', '.join(p.get('experimentTypes') or []),",
        "    'samples': None, 'taxa': ', '.join(p.get('organisms') or []),",
        "    'pubmed': '', 'supplementary': '', 'date': p.get('submissionDate', ''),",
        "    'parasite_measured': any('toxoplasma' in str(o).lower() for o in (p.get('organisms') or [])),",
        "} for p in rows])",
        "len(pride), int(pride['parasite_measured'].sum())")
nb.md("## 3. Metabolomics Workbench")
nb.code("mw = []",
        "for field in ('study_title', 'species'):",
        "    q = f'https://www.metabolomicsworkbench.org/rest/study/{field}/Toxoplasma/summary'",
        "    try:",
        "        r = json.loads(logged(q))",
        "    except Exception as e:",
        "        print(field, 'failed:', e); continue",
        "    r = [r] if isinstance(r, dict) and 'study_id' in r else list(r.values()) if isinstance(r, dict) else r",
        "    mw += r; time.sleep(0.3)",
        "meta = pd.DataFrame([{",
        "    'source': 'Metabolomics Workbench', 'accession': s.get('study_id'), 'title': s.get('study_title', ''),",
        "    'summary': s.get('study_summary', ''), 'type': 'metabolomics: ' + str(s.get('analysis_type', '')),",
        "    'samples': int(s['number_of_samples']) if str(s.get('number_of_samples', '')).isdigit() else None,",
        "    'taxa': s.get('species', ''), 'pubmed': '', 'supplementary': '', 'date': s.get('release_date', ''),",
        "    'parasite_measured': 'toxoplasma' in str(s.get('species', '')).lower(),",
        "} for s in mw]).drop_duplicates('accession')",
        "len(meta)")
nb.md("## 4. Classification",
      "What was perturbed, by keyword on title and summary. A study can be several things at once (a "
      "knockout in two strains). The matched text is kept so each call can be checked.")
nb.code("cat = pd.concat([geo, sra_only, pride, meta], ignore_index=True).drop_duplicates('accession')",
        "TEXT = (cat['title'] + ' ' + cat['summary']).str.replace(r'\\s+', ' ', regex=True)",
        "RULES = {",
        "    'gene perturbation': r'knock-?out|knock-?down|\\bdelet(?:ion|ed)\\b|\\bnull mutant|\\bmutant\\b|\\u0394\\w|\\bKO\\b|auxin|\\bm?AID\\b|degron|conditional (?:knock|deplet)|deplet(?:ion|ed) of|CRISPR|tet-?off|ATc|overexpress',",
        "    'drug or compound': r'\\bdrugs?\\b|\\bcompounds?\\b|inhibitors?\\b|treated with (?!trypsin)|treatment with (?!trypsin)|pyrimethamine|sulfadiazine|atovaquone|artemisinin|clindamycin|endochin|buparvaquone|\\bIC50\\b|antiparasitic|anti-?toxoplasma',",
        "    'strain comparison': r'\\bstrains?\\b.{0,80}?(?:type I\\b|type II\\b|type III|\\bRH\\b|ME49|PRU\\b|VEG\\b|GT1|CTG|CEP\\b)|(?:type I\\b|type II\\b|type III).{0,60}?(?:type I\\b|type II\\b|type III)|\\b(?:RH|GT1)\\b.{0,60}?\\b(?:ME49|PRU|VEG|CTG)\\b',",
        "    'stage or differentiation': r'bradyzoite|tachyzoite.{0,40}bradyzoite|differentiation|cyst|sporozoite|oocyst|merozoite|sexual|gametocyte|life cycle',",
        "    'host response': r'host cell|infected (?:cells|mice|fibroblast|macrophage)|human foreskin|HFF|mouse brain|macrophage|dendritic',",
        "}",
        "for name, rx in RULES.items():",
        "    m = TEXT.str.extract('(' + rx + ')', flags=re.I)[0]",
        "    cat[name] = m.notna(); cat[name + ' (matched)'] = m.fillna('')",
        "def readout(t):",
        "    t = str(t).lower()",
        "    for k, v in (('proteomics', 'proteome'), ('metabolomics', 'metabolome'), ('methylation', 'methylation'),",
        "                 ('genome binding', 'chromatin / binding'), ('genome variation', 'genome'), ('non-coding', 'non-coding RNA'),",
        "                 ('high throughput sequencing', 'RNA-seq / sequencing'), ('array', 'microarray')):",
        "        if k in t: return v",
        "    return 'other'",
        "cat['readout'] = cat['type'].map(readout)",
        "cat.groupby('readout')[list(RULES)].sum().assign(studies=cat.groupby('readout').size())")
nb.md("## 5. The thioredoxin knockout and the strain comparisons asked for",
      "Exact-word matches, listed whole so a person can read them.")
nb.code("trx = cat[TEXT.str.contains(r'thioredoxin|\\bTRX\\d?\\b|\\bTrxL', case=False, regex=True)]",
        "trx[['source', 'accession', 'readout', 'samples', 'title']]")
nb.code("strains = cat[cat['strain comparison'] & cat['parasite_measured']]",
        "strains[['source', 'accession', 'readout', 'samples', 'strain comparison (matched)', 'title']].sort_values('samples', ascending=False)")
nb.md("## 6. Gene perturbations and drugs measured on the parasite",
      "The pool for loading: the parasite is measured (not only the host), and something was knocked "
      "out, knocked down or dosed. Samples is the series total; replicate structure is checked per "
      "study when it is loaded.")
nb.code("pool = cat[cat['parasite_measured'] & (cat['gene perturbation'] | cat['drug or compound'])]",
        "print(len(pool), 'studies;', pool.groupby('readout').size().to_dict())",
        "pool[['source', 'accession', 'readout', 'samples', 'gene perturbation (matched)', 'drug or compound (matched)', 'title']].sort_values(['readout', 'samples'], ascending=[True, False]).head(60)")
nb.code("cat.to_csv(os.path.join(OUT, 'catalogue.csv'), index=False)",
        "pd.Series(QUERIES).to_csv(os.path.join(OUT, 'queries.txt'), index=False, header=False)",
        "len(cat), len(QUERIES)")
nb.write(os.path.join(ROOT, "notebooks", "omics_catalogue_2026_10_17.ipynb"))
