"""Instruction 72: the FUNCTION truth set for the information map.

Writes `notebooks/function_truth_2026_11_04.ipynb` (executed) and the tables under
`results/function_truth_2026_11_04/`. Raw API responses come from `scripts/fetch_function_sources.py`,
which pins them with their retrieval times before this script runs.

The third target trait of instruction 72 is "what is the function of all Toxoplasma proteins". The map
held experimental GO function annotations for ~25 genes, far too few to measure anything at the 90 %
precision bar. This builds the labelled set the function column has to be scored on, from three
sources kept apart and graded, the way the host-facing truth set was built.
"""
import sys
import textwrap  # noqa: F401  (kept: cell() documents why dedent is not used)
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from notebook_runner import ExecutedNotebook  # noqa: E402


def cell(nb, source: str):
    """Run one cell. Strips the 4-space block indent this file uses, line by line, rather than with
    `textwrap.dedent`: one cell embeds a Markdown literal whose own lines start at column 0, which
    would leave dedent with no common prefix to remove."""
    body = '\n'.join(line[4:] if line.startswith('    ') else line for line in source.split('\n'))
    return nb.code(body.strip('\n'))


nb = ExecutedNotebook('The function truth set for the information map')

nb.md('# The function truth set for the information map',
      'Instruction 72, target trait 3 (Einar, 2026-11-04). The information map must predict **what a '
      'Toxoplasma protein does**, and a prediction can only be scored against labels that were not '
      'themselves produced by a predictor. Before this notebook the map held experimental GO '
      'biological-process / molecular-function annotations for about 25 genes: too few to measure '
      'precision at the 0.90 fill bar, let alone its lower confidence bound.',
      '**Three sources, in descending evidence strength, kept separate and recorded per row:**',
      '1. **Experimental GO** -- QuickGO, taxon 5811 with descendants, aspects `biological_process` and '
      '`molecular_function`, evidence restricted to the experimental codes (ECO:0000269 EXP, 0000314 IDA, '
      '0000353 IPI, 0000315 IMP, 0000316 IGI, 0000270 IEP). Tier `observed-direct`, probability 0.95.',
      '2. **Reviewed UniProt (Swiss-Prot)** entries for *Toxoplasma gondii*: the curator-written `FUNCTION` '
      'comment, assembled from published work and citing its PubMed ids. Tier `observed-curated`, 0.90.',
      '3. **Phenotype-backed process claims** -- the arrayed splitCas9 imaging screen in `nodes.parquet` '
      '(318 genes scored for egress, F-actin, apicoplast and replication; PMID 35538310). These support a '
      '*process* claim and never a molecular function. Tier `observed-screen`, 0.80.',
      '**Excluded on purpose, and this is why the set is small:** InterPro/Pfam domain assignments, ToxoDB '
      '`product` names, and every GO annotation whose evidence code is IEA, ISS, ISO or another '
      'electronic/similarity code. All of those are *predictions from sequence*. They are the thing the '
      'information map is supposed to predict; training or scoring on them would measure agreement with a '
      'predictor rather than with biology (instruction 70) and would repeat the failure diagnosed in both '
      'host-facing scoring runs, where the strongest features turned out to be downstream of the conclusion.',
      '**The probabilities below are provisional tier defaults** (instruction 72 v0), not measured error '
      'rates. v1 replaces them by how often each source is confirmed by an independent one.')

nb.md('## 1 · Inputs and the pinned raw responses',
      'Every response was written to `results/function_truth_2026_11_04/sources/` with its URL, HTTP status, '
      'retrieval time in UTC and sha256 (`sources/manifest.json`). The 32 MB GO ontology is pinned beside the '
      'other dataset snapshots instead of inside the repository. Requests were made one at a time with a '
      '0.5 s pause, well inside the 3 requests/s budget.')
cell(nb, """
    import json, re, csv, collections, math, hashlib, urllib.request, time
    from datetime import datetime, timezone
    from pathlib import Path
    import numpy as np, pandas as pd

    OUT = Path('results/function_truth_2026_11_04'); SRC = OUT / 'sources'
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((SRC / 'manifest.json').read_text())
    pd.DataFrame(manifest)[['source', 'file', 'status', 'retrieved_utc', 'bytes', 'sha256']]
""")

nb.md('### The gene universe',
      'The map has 8,140 ME49 rows (`nodes.parquet`). A label that cannot be placed on one of those rows is '
      'useless to the map, so unmapped and ambiguous accessions are listed rather than forced onto a gene.')
cell(nb, """
    nodes = pd.read_parquet('parascape/data/nodes.parquet')
    universe = set(nodes.gene_id)
    print(f'{len(universe):,} ME49 genes in the map')
""")

nb.md('## 2 · Accession to TGME49 gene id',
      'Both GO and UniProt key on UniProt accessions, and Toxoplasma has reviewed entries from several strains '
      '(ME49, GT1, VEG, RH88, ARI...), so a mapping step is unavoidable. Two independent routes are used and '
      'compared:',
      '* **route A -- UniProt own cross-references**: the `VEuPathDB` `ToxoDB:TGME49_*` ids of each reviewed entry '
      '(fetched in a second pass, because the field list of the main query suppresses every cross-reference but GO).',
      '* **route B -- the installed DeepTMHMM table**: `datasets/deeptmhmm_toxo.csv`, whose `accession` column holds '
      'the comma-separated UniProt accessions of each of the 8,140 genes.',
      'Where both routes fire they are checked against each other. An accession reaching more than one ME49 gene is '
      'recorded as **ambiguous** and its labels are held back, not spread over the candidates.')
cell(nb, """
    xref = json.loads((SRC / 'uniprot_reviewed_5811_xref_page1.json').read_text())['results']
    route_a = collections.defaultdict(set)
    for entry in xref:
        for ref in entry.get('uniProtKBCrossReferences', []):
            if ref['database'] == 'VEuPathDB' and ref['id'].startswith('ToxoDB:TGME49_'):
                route_a[entry['primaryAccession']].add(ref['id'].split(':', 1)[1])

    dtm = pd.read_csv('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/deeptmhmm_toxo.csv',
                      usecols=['identifier', 'accession'])
    route_b = collections.defaultdict(set)
    for gene, accs in zip(dtm.identifier, dtm.accession):
        if isinstance(accs, str):
            for acc in accs.split(','):
                route_b[acc.strip()].add(gene)

    # Route D, the fallback: many reviewed entries are the GT1, VEG or RH88 copy of a gene and carry no
    # TGME49 cross-reference at all. ToxoDB keeps the numeric suffix of a gene id across strains for the
    # syntenic ortholog (O00834 is cross-referenced as TGARI_291890, TGGT1_291890, TGME49_291890 ... all
    # 291890), so an entry whose only ToxoDB ids share one suffix is placed on TGME49_<suffix> when that
    # gene exists. This is an inference from an id convention, not a measurement, so it is labelled and
    # every gene reached this way is penalised in the review queue.
    route_d = {}
    for entry in xref:
        ids = [ref['id'].split(':', 1)[1] for ref in entry.get('uniProtKBCrossReferences', [])
               if ref['database'] == 'VEuPathDB' and ref['id'].startswith('ToxoDB:')]
        suffixes = {i.rsplit('_', 1)[-1] for i in ids}
        if len(suffixes) == 1:
            candidate = f'TGME49_{suffixes.pop()}'
            if candidate in universe:
                route_d[entry['primaryAccession']] = {candidate}

    def resolve(acc):
        a, b = route_a.get(acc, set()) & universe, route_b.get(acc, set()) & universe
        if a and b and a != b:
            return sorted(a | b), 'A+B disagree'
        if a:
            return sorted(a), 'A (UniProt ToxoDB xref)'
        if b:
            return sorted(b), 'B (DeepTMHMM accession)'
        if acc in route_d:
            return sorted(route_d[acc]), 'D (ToxoDB strain-suffix, inferred)'
        return [], 'unmapped'

    print(f'route A covers {len(route_a):,} reviewed accessions, route B covers {len(route_b):,} '
          f'accessions overall, route D adds {len(set(route_d) - set(route_a)):,} reviewed accessions')
""")

nb.md('## 3 · Source 1 -- experimental GO (tier `observed-direct`, 0.95)',
      'QuickGO returns one row per annotation with its ECO id, GO evidence code and reference. Only the six '
      'experimental codes asked for are kept; the high-throughput experimental codes (HDA, HMP, HGI, HEP) and '
      'every electronic or similarity code are dropped, and the counts of what was dropped are printed so the '
      'exclusion is visible rather than silent.')
cell(nb, """
    EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'}
    go_rows = []
    dropped = collections.Counter()
    for aspect in ('biological_process', 'molecular_function'):
        with open(SRC / f'quickgo_{aspect}_experimental.tsv') as fh:
            for row in csv.DictReader(fh, delimiter='\\t'):
                code = row['GO EVIDENCE CODE']
                if code not in EXPERIMENTAL:
                    dropped[code] += 1
                    continue
                go_rows.append(dict(accession=row['GENE PRODUCT ID'], symbol=row['SYMBOL'],
                                    qualifier=row['QUALIFIER'], go_id=row['GO TERM'],
                                    aspect=row['GO ASPECT'], eco=row['ECO ID'], code=code,
                                    reference=row['REFERENCE'], taxon=row['TAXON ID'],
                                    assigned_by=row['ASSIGNED BY'], date=row['DATE']))
    go = pd.DataFrame(go_rows)
    print(f'{len(go)} experimental GO annotations, {go.accession.nunique()} gene products, '
          f'{go.go_id.nunique()} distinct terms; dropped by evidence code: {dict(dropped) or "none"}')
    go.code.value_counts().to_frame('annotations').T
""")

nb.md('For GO rows the `SYMBOL` field gives two more routes, used only when the accession routes find '
      'nothing: the symbol is sometimes already a TGME49 id (route C), and sometimes another strain id whose '
      'suffix identifies the ME49 gene (route D again, and labelled as the inference it is).')
cell(nb, """
    def go_gene(row):
        genes, route = resolve(row.accession)
        if genes:
            return genes, route
        symbol = str(row.symbol)
        if re.fullmatch(r'TGME49_\\w+', symbol) and symbol in universe:
            return [symbol], 'C (QuickGO symbol)'
        if re.fullmatch(r'TG[A-Z0-9]+_\\d+', symbol):
            candidate = f'TGME49_{symbol.rsplit("_", 1)[-1]}'
            if candidate in universe:
                return [candidate], 'D (ToxoDB strain-suffix, inferred)'
        return [], 'unmapped'

    go['genes'], go['map_route'] = zip(*[go_gene(r) for r in go.itertuples()])
    go['gene_id'] = [g[0] if len(g) == 1 else '' for g in go.genes]
    go['ambiguous'] = [len(g) > 1 for g in go.genes]
    print(f'{(go.gene_id != "").sum()} annotations placed on a single ME49 gene '
          f'({go.loc[go.gene_id != "", "gene_id"].nunique()} genes); '
          f'{int(go.ambiguous.sum())} ambiguous, {(go.map_route == "unmapped").sum()} unmapped')
    go.map_route.value_counts().to_frame('annotations')
""")

nb.md('## 4 · Source 2 -- reviewed UniProt curator function (tier `observed-curated`, 0.90)',
      'A reviewed (Swiss-Prot) entry is read and written by a curator from the published experiments, and its '
      '`FUNCTION` comment carries the PubMed ids and ECO codes of each statement. That is the largest machine-readable '
      'body of real Toxoplasma function knowledge there is. The `GO` cross-references of the same entries are also '
      'read, again keeping only experimental evidence types -- they turn out to be a **subset** of the QuickGO rows, '
      'which is the expected relationship (UniProt-GOA feeds QuickGO) and a useful consistency check.',
      'The recommended protein name is kept as a secondary text field, because in a reviewed entry it is a '
      'curator statement from the literature (ECO:0000303) rather than a sequence prediction. **UniProt keywords '
      'are recorded but never used to assign a class**: some keywords (Transmembrane, Signal) come from automatic '
      'annotation, which is exactly the excluded kind of evidence.')
cell(nb, """
    entries = json.loads((SRC / 'uniprot_reviewed_5811_page1.json').read_text())['results']
    up_rows, go_xref = [], []
    for e in entries:
        acc = e['primaryAccession']
        genes, route = resolve(acc)
        desc = e.get('proteinDescription', {})
        name = desc.get('recommendedName', {}).get('fullName', {}).get('value', '')
        gene_names = [g.get('geneName', {}).get('value') for g in e.get('genes', []) if g.get('geneName')]
        orf = [o['value'] for g in e.get('genes', []) for o in g.get('orfNames', [])]
        texts, pmids, ecos = [], set(), set()
        for comment in e.get('comments', []):
            if comment.get('commentType') != 'FUNCTION':
                continue
            for t in comment.get('texts', []):
                texts.append(t['value'])
                for ev in t.get('evidences', []):
                    ecos.add(ev.get('evidenceCode', ''))
                    if ev.get('source') == 'PubMed':
                        pmids.add(ev['id'])
        domains = [f.get('description', '') for f in e.get('features', []) if f.get('type') == 'Domain']
        up_rows.append(dict(accession=acc, uniprot_id=e['uniProtkbId'], protein_name=name,
                            gene_names='|'.join(filter(None, gene_names)), orf_names='|'.join(orf),
                            function_text=' '.join(texts), function_pmids='|'.join(sorted(pmids)),
                            function_eco='|'.join(sorted(ecos)),
                            keywords='|'.join(k['name'] for k in e.get('keywords', [])),
                            domains='|'.join(filter(None, domains)),
                            lit_pmids='|'.join(sorted({r['citation']['id'] for r in e.get('references', [])
                                                       if r.get('citation', {}).get('citationCrossReferences')
                                                       and r['citation'].get('id', '').isdigit()})),
                            genes=genes, map_route=route))
        for ref in e.get('uniProtKBCrossReferences', []):
            if ref['database'] != 'GO':
                continue
            props = {p['key']: p['value'] for p in ref.get('properties', [])}
            evidence = props.get('GoEvidenceType', '').split(':')[0]
            aspect = props.get('GoTerm', '')[:1]
            if evidence in EXPERIMENTAL and aspect in 'PF':
                go_xref.append(dict(accession=acc, go_id=ref['id'], aspect=aspect, code=evidence,
                                    term=props.get('GoTerm', '')[2:]))
    up = pd.DataFrame(up_rows)
    up['gene_id'] = [g[0] if len(g) == 1 else '' for g in up.genes]
    xref_go = pd.DataFrame(go_xref)
    print(f'{len(up)} reviewed entries; {(up.function_text != "").sum()} carry a FUNCTION comment; '
          f'{(up.gene_id != "").sum()} map to one ME49 gene, '
          f'{(up.map_route == "unmapped").sum()} unmapped, {sum(len(g) > 1 for g in up.genes)} ambiguous')
    print(f'{len(xref_go)} experimental GO cross-references on the reviewed entries, '
          f'{len(set(zip(xref_go.go_id, xref_go.aspect)) - set(zip(go.go_id, go.aspect)))} of their '
          f'(term, aspect) pairs absent from the QuickGO download')
    up.map_route.value_counts().to_frame('entries')
""")

nb.md('### What the unmapped accessions are',
      'Anything left without an ME49 gene is printed in full. These are real function labels that the map cannot '
      'yet use, and guessing a gene for them would put curated biology on the wrong row.')
cell(nb, """
    unmapped = up[(up.gene_id == '') | (up.map_route == 'unmapped')][
        ['accession', 'uniprot_id', 'protein_name', 'gene_names', 'orf_names', 'map_route']]
    unmapped.to_csv(OUT / 'unmapped_accessions.tsv', sep='\\t', index=False)
    print(f'{len(unmapped)} reviewed entries with no single ME49 gene')
    unmapped.head(30)
""")

nb.md('## 5 · The GO ontology: term names, and how specific each term is',
      "QuickGO's TSV gives only the GO id, so names come from the `go-basic.obo` snapshot. The same parse gives the "
      '**information content** instruction 72 asks for, so that a vague answer cannot score like a specific one.',
      'Two measures are reported, because they answer different questions:',
      '* `ic_structural` = −log2((1 + descendants) / terms in the aspect) -- how specific the term is *in the '
      'ontology*, independent of our small set. "protein binding" sits near the root and scores near 0; '
      '"actin filament severing activity" is a leaf and scores high.',
      '* `ic_corpus` = −log2(genes carrying the label / labelled genes) -- how surprising the label is *in this '
      'truth set*. With a few hundred genes this is noisy, and is reported as such.',
      'Ancestors are followed over `is_a` and `part_of` only, which is what `go-basic` is filtered for.')
cell(nb, """
    OBO = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/sources/go_basic_2026_11_04/go-basic.obo')
    names, namespace, parents, obsolete = {}, {}, collections.defaultdict(set), set()
    by_name = {}
    current = None
    for line in OBO.read_text().splitlines():
        if line in ('[Term]', '[Typedef]'):
            current = {} if line == '[Term]' else None
        elif current is None:
            continue
        elif line.startswith('id: GO:'):
            current['id'] = line[4:]
        elif line.startswith('name: '):
            names[current['id']] = line[6:]
        elif line.startswith('namespace: '):
            namespace[current['id']] = line[11:]
        elif line.startswith('is_a: GO:'):
            parents[current['id']].add(line[6:16])
        elif line.startswith('relationship: part_of GO:'):
            parents[current['id']].add(line[22:32])
        elif line.startswith('is_obsolete: true'):
            obsolete.add(current['id'])
    for go_id, name in names.items():
        if go_id not in obsolete:
            by_name.setdefault(name, go_id)

    children = collections.defaultdict(set)
    for child, ps in parents.items():
        for p in ps:
            children[p].add(child)

    def closure(seed, graph):
        out, stack = set(), [seed]
        while stack:
            node = stack.pop()
            for nxt in graph.get(node, ()):
                if nxt not in out:
                    out.add(nxt)
                    stack.append(nxt)
        return out

    aspect_size = collections.Counter(ns for go_id, ns in namespace.items() if go_id not in obsolete)
    ancestors = {}
    def anc(go_id):
        if go_id not in ancestors:
            ancestors[go_id] = closure(go_id, parents)
        return ancestors[go_id]

    def ic_structural(go_id):
        ns = namespace.get(go_id)
        if not ns:
            return float('nan')
        return -math.log2((1 + len(closure(go_id, children))) / aspect_size[ns])

    print(f'{len(names) - len(obsolete):,} current GO terms '
          f'({aspect_size["biological_process"]:,} BP, {aspect_size["molecular_function"]:,} MF)')
    print('examples:',
          {names[t]: round(ic_structural(t), 2) for t in ('GO:0005515', 'GO:0003789', 'GO:0035891')})
""")

nb.md('## 6 · The label space, and why free text had to become classes',
      'Instruction 72 needs a function column that can be *scored*. A curator paragraph cannot be scored, and the '
      'raw GO terms are too sparse to score: 85 distinct experimental terms over about 50 genes means most terms '
      'have one or two genes, so held-out precision could never be measured per term. So the gene gets one or more '
      '**coarse classes** as the scoring convenience, while the fine labels -- the GO terms and the curator text -- '
      'are kept in the long-form evidence table as the actual truth.',
      '15 classes, chosen to match the way Toxoplasma biology is divided in the literature and the way this '
      'particular evidence actually falls (every class below is populated; nothing was invented for symmetry):',
      '`host_cell_modulation`, `invasion_motility`, `egress`, `cell_cycle_division`, '
      '`protein_trafficking_secretion`, `signalling`, `transport_homeostasis`, `metabolism`, '
      '`lipid_membrane_biology`, `apicoplast_mitochondrion`, `translation`, `transcription_chromatin`, '
      '`stress_differentiation`, `proteostasis`, `cytoskeleton`.',
      'Rules are of two kinds. A **GO rule** names a seed term by its ontology name; a gene matches the class if '
      'one of its experimental GO terms is that seed or a descendant of it (`is_a`/`part_of` closure). A **text '
      'rule** is a regular expression over the curator `FUNCTION` comment and the recommended protein name. '
      'GO rules are the stronger kind and are applied first; a class reached by a GO rule keeps the GO tier.',
      'Two deliberate choices inside the rules: `host_cell_modulation` uses only the *perturbation/modulation of '
      'host process* branch, not the broad "biological process involved in interaction with host" (which would '
      'swallow invasion and egress); and `(?<!host )mitochondri` keeps GRA/MAF proteins that tether the vacuole to '
      '**host** mitochondria out of `apicoplast_mitochondrion`, where they would be plain wrong.')
cell(nb, """
    CLASS_RULES = {
      'host_cell_modulation': dict(
        go=['modulation of process of another organism',
            'symbiont-mediated perturbation of host process',
            'symbiont-mediated response to host defenses',
            'Toll-like receptor binding'],
        text=[r'suppress\\w*[^.]{0,60}host', r'modulat\\w*[^.]{0,60}host',
              r'host (innate )?immun', r'host (interferon|cytokine|transcription|gene expression)',
              r'immunity-related GTPase', r'\\bIRG[AB]?\\d*\\b', r'\\bIFN-gamma\\b',
              r'export\\w*[^.]{0,40}(into|to) the host', r'tether\\w*[^.]{0,60}host',
              r'host (mitochondri|nucleus|nuclear|cytosol|endocytic|lysosom)',
              r'\\bvirulence (factor|in mice)\\b', r'NF-kappa-?B', r'\\bSTAT[136]\\b',
              r'effector[^.]{0,40}(export|translocat)', r'host cell cycle']),
      'invasion_motility': dict(
        go=['symbiont entry into host', 'adhesion of symbiont to host', 'host cell surface binding',
            'regulation of cell motility', 'cell adhesion molecule production', 'cell-cell adhesion'],
        text=[r'invasion', r'invad\\w+', r'adhesin', r'attachment (of|to)', r'moving junction',
              r'motility', r'glid\\w+', r'conoid', r'apical complex', r'binding to host cells?',
              r'microneme secretion', r'rhoptry secretion', r'\\bMJ complex\\b']),
      'egress': dict(
        go=['exit from host cell'],
        text=[r'egress', r'exit from (the )?host', r'lysis of the host cell',
              r'permeabiliz\\w+[^.]{0,60}vacuol', r'parasite exit']),
      'cell_cycle_division': dict(
        go=['cell cycle', 'cell division', 'chromosome segregation',
            'regulation of cell population proliferation', 'DNA replication'],
        text=[r'cell (division|cycle)', r'daughter (cell|parasite)', r'centrosome', r'cytokinesis',
              r'\\bbudding\\b', r'mitosis', r'mitotic', r'G[12]/[MS]', r'intracellular (multiplication|replication)',
              r'parasite (division|replication|multiplication)', r'chromosome', r'\\bcell polarity\\b']),
      'protein_trafficking_secretion': dict(
        go=['intracellular protein localization', 'regulated exocytosis', 'vesicle-mediated transport',
            'regulation of regulated secretory pathway', 'vesicle fusion to plasma membrane',
            'protein secretion', 'Golgi organization'],
        text=[r'traffick\\w+', r'secretion', r'secretory', r'exocytosis', r'export of',
              r'import of[^.]{0,60}into the', r'vesicul\\w+', r'localization to[^.]{0,40}membrane',
              r'escorter', r'sort\\w+ (of|into)', r'\\bGolgi\\b', r'endoplasmic reticulum',
              r'dense granule (transport|exocytosis)']),
      'signalling': dict(
        go=['protein kinase activity', 'protein phosphorylation', 'signal transduction',
            'cyclic nucleotide biosynthetic process', 'cyclic nucleotide binding',
            'phosphoprotein phosphatase activity', 'cellular response to chemical stimulus'],
        text=[r'\\bkinase\\b', r'\\bcyclase\\b', r'phosphorylat\\w+', r'phosphatase',
              r'second messenger', r'\\bcGMP\\b', r'\\bcAMP\\b', r'signal transduc\\w+',
              r'calcium[- ]dependent', r'\\bpseudokinase\\b', r'senses? ', r'detects and transduces']),
      'transport_homeostasis': dict(
        go=['transmembrane transporter activity', 'transmembrane transport', 'channel activity',
            'monoatomic ion homeostasis', 'cellular homeostasis', 'export across plasma membrane'],
        text=[r'transporter', r'\\bchannel\\b', r'permeable to', r'homeostasis',
              r'import across', r'uptake of', r'\\befflux\\b', r'extrusion of',
              r'\\btransports?\\b[^.]{0,40}?\\b(across|into|out of)\\b',
              r'\\btransports?\\b[^.]{0,30}?\\b(ions?|water|glycerol|urea|pyruvate|glucose|hexose|'
              r'amino acids?|lysine|arginine|methionine|nucleotides?|sugars?|lipids?|cholesterol)\\b']),
      'metabolism': dict(
        go=['small molecule metabolic process', 'carbohydrate metabolic process',
            'nucleobase-containing small molecule metabolic process', 'amino acid metabolic process',
            'generation of precursor metabolites and energy',
            'catalytic activity, acting on a nucleic acid'],
        text=[r'catalyz\\w+', r'biosynthe\\w+', r'metaboli\\w+', r'glycolysis', r'glycolytic',
              r'salvage', r'\\bnucleotide', r'\\bpurine\\b', r'\\bpyrimidine\\b', r'tricarboxylic acid',
              r'pentose phosphate', r'\\bfolate\\b', r'\\bdehydrogenase\\b', r'\\bsynthase\\b',
              r'\\btransferase\\b(?! )', r'amino acid (biosynthesis|metabolism)']),
      'lipid_membrane_biology': dict(
        go=['lipid metabolic process', 'lipid binding', 'membrane organization',
            'lipid transport', 'phospholipid biosynthetic process'],
        text=[r'\\blipid', r'phospholipid', r'phosphatidyl', r'fatty acid', r'\\bsterol\\b',
              r'cholesterol', r'membrane (fusion|organization|remodel)', r'binds lipid',
              r'\\bflippase\\b', r'\\bpore[- ]forming\\b', r'(membranous|intravacuolar) network']),
      'apicoplast_mitochondrion': dict(
        go=[],   # cellular-component terms are out of scope here; this class is text-derived only
        text=[r'apicoplast', r'(?<!host )mitochondri\\w+', r'\\bplastid\\b', r'isoprenoid',
              r'methylerythritol', r'\\bMEP pathway\\b', r'\\bFASII\\b', r'respiratory chain']),
      'translation': dict(
        go=['translation', 'ribosome biogenesis', 'tRNA metabolic process',
            'structural constituent of ribosome', 'translation regulator activity'],
        text=[r'\\btranslation (of|initiation|elongation|termination)\\b',
              r'translational (control|regulation|repression)', r'ribosom\\w+', r'\\btRNA\\b',
              r'aminoacyl', r'\\belongation factor\\b', r'\\binitiation factor\\b']),
      'transcription_chromatin': dict(
        go=['DNA-binding transcription factor activity', 'regulation of DNA-templated transcription',
            'chromatin organization', 'mRNA processing', 'regulation of gene expression',
            'DNA binding', 'RNA binding'],
        text=[r'(?<!host )transcription', r'chromatin', r'histone', r'promoter', r'\\bmRNA\\b',
              r'\\bsplic\\w+', r'(?<!host )gene expression', r'DNA[- ]binding',
              r'\\bRNA polymerase\\b', r'\\bcapping\\b']),
      'stress_differentiation': dict(
        go=['response to stress', 'cellular developmental process', 'dormancy process',
            'cell fate commitment', 'regulation of cell differentiation'],
        text=[r'bradyzoite', r'differentiat\\w+', r'\\bcyst\\b', r'cyst wall', r'dormanc\\w+',
              r'\\bstress\\b', r'heat shock', r'oxidative', r'starvation', r'chronic[- ]stage',
              r'\\bpersistence\\b', r'reactivation']),
      'proteostasis': dict(
        go=['proteolysis', 'protein folding', 'autophagy', 'peptidase activity',
            'protein maturation', 'protein modification by small protein conjugation or removal'],
        text=[r'\\bprotease\\b', r'peptidase', r'proteoly\\w+', r'(?<!host )autophag\\w+', r'chaperon\\w+',
              r'\\bfolding\\b', r'ubiquitin', r'proteasom\\w+', r'degradation of',
              r'processing[- ]mediated maturation', r'\\bcleav\\w+\\b', r'\\bisomerase\\b']),
      'cytoskeleton': dict(
        go=['cytoskeleton organization', 'actin filament-based process', 'actin binding',
            'microtubule-based process', 'structural constituent of cytoskeleton', 'cytoskeletal motor activity',
            'actin filament severing activity', 'actin monomer binding'],
        text=[r'\\bactin\\b', r'tubulin', r'microtubule', r'\\bmyosin\\b', r'cytoskelet\\w+',
              r'inner membrane complex', r'\\bIMC\\b', r'\\bspectrin\\b', r'\\bseptin\\b',
              r'peripheral annuli', r'\\bF-actin\\b']),
    }
    SCREEN_RULES = {   # phenotype column -> the process class it supports, and nothing more
      'screen_egress_phenotype': 'egress',
      'screen_actin_phenotype': 'cytoskeleton',
      'screen_apicoplast_phenotype': 'apicoplast_mitochondrion',
      'screen_replication_phenotype': 'cell_cycle_division',
    }
    CLASSES = list(CLASS_RULES)
    print(len(CLASSES), 'classes;',
          sum(len(v['go']) for v in CLASS_RULES.values()), 'GO seed terms,',
          sum(len(v['text']) for v in CLASS_RULES.values()), 'text patterns')
""")

nb.md('Every GO seed must resolve to a current ontology term, or the rule is silently dead. That is checked here '
      'and the resolved id, its own structural information content and the size of its descendant set are printed '
      'in full -- this table plus the regular expressions above is what `class_rules.md` records verbatim.')
cell(nb, """
    seeds, missing = {}, []
    for cls, rule in CLASS_RULES.items():
        for name in rule['go']:
            go_id = by_name.get(name)
            if go_id is None:
                missing.append((cls, name))
            else:
                seeds.setdefault(cls, []).append((name, go_id))
    assert not missing, f'GO seed names not found in the ontology: {missing}'
    seed_table = pd.DataFrame([dict(cls=cls, seed=name, go_id=go_id, namespace=namespace[go_id],
                                    descendants=len(closure(go_id, children)),
                                    ic_structural=round(ic_structural(go_id), 2))
                               for cls, rows in seeds.items() for name, go_id in rows])
    print(f'{len(seed_table)} GO seeds, all resolved')
    seed_table
""")

nb.md('## 7 · Applying the rules',
      'For each gene the three sources are applied independently so that the deciding source of every class is '
      'recorded. A class reached from an experimental GO term is `observed-direct` (0.95); reached only from the '
      'curator text, `observed-curated` (0.90); reached only from a screen phenotype, `observed-screen` (0.80).')
cell(nb, """
    SEED_IDS = {cls: {go_id for _, go_id in rows} for cls, rows in seeds.items()}

    def classes_from_go(go_ids):
        hits = collections.defaultdict(set)
        for go_id in go_ids:
            reach = {go_id} | anc(go_id)
            for cls, seed_ids in SEED_IDS.items():
                if reach & seed_ids:
                    hits[cls].add(go_id)
        return hits

    PATTERNS = {cls: [(p, re.compile(p, re.I)) for p in rule['text']] for cls, rule in CLASS_RULES.items()}
    #: Curators write what a protein does NOT do as often as what it does: ASP1 is "dispensable for
    #: ... tachyzoite and bradyzoite viability", MAF1a1 "does not play a role in tethering the
    #: parasitophorous vacuole to host mitochondria", FER2 is "Not required for microneme secretion".
    #: Taking those at face value would hand out exactly the wrong class. A match is therefore
    #: discarded when a negation cue appears BEFORE it inside the same sentence.
    NEGATION = re.compile(r'\\b(not|non|never|cannot|dispensab\\w*|without|lacks?|absence of|'
                          r'fails? to|neither|unable)\\b', re.I)
    #: Second guard, for the three classes where the HOST has the same process. ROP17 "May modulate
    #: gene expression in human cells" and "Promotes autophagy in human cells" are statements about
    #: the host, not about parasite transcription or parasite proteostasis -- that gene is already
    #: host_cell_modulation from its own patterns. A match of these three classes inside a sentence
    #: that names the host, human cells, mice or macrophages is therefore discarded.
    HOST_CONTEXT = re.compile(r'\\b(host|human cells?|mouse|mice|murine|macrophages?|monocytes?|'
                              r'infected cells?)\\b', re.I)
    HOST_AMBIGUOUS = {'transcription_chromatin', 'translation', 'proteostasis'}
    suppressed = []

    def sentences(text):
        offset = 0
        for part in re.split(r'(?<=[.;])\\s+', text or ''):
            yield offset, part
            offset += len(part) + 1

    def classes_from_text(text, gene=''):
        hits = collections.defaultdict(list)
        for cls, pats in PATTERNS.items():
            for raw, pat in pats:
                for _, sentence in sentences(text):
                    found = pat.search(sentence)
                    if not found:
                        continue
                    cue = NEGATION.search(sentence[:found.start()])
                    host = HOST_CONTEXT.search(sentence) if cls in HOST_AMBIGUOUS else None
                    if cue or host:
                        suppressed.append(dict(gene_id=gene, cls=cls, pattern=raw,
                                               cue=cue.group(0) if cue else f'host context: {host.group(0)}',
                                               guard='negation' if cue else 'host context',
                                               sentence=sentence[:160]))
                        continue
                    hits[cls].append(raw)
                    break
        return hits

    assignments = []     # one row per gene x class x source
    for gene, rows in go[go.gene_id != ''].groupby('gene_id'):
        for cls, terms in classes_from_go(sorted(set(rows.go_id))).items():
            refs = sorted({r for r in rows.loc[rows.go_id.isin(terms), 'reference'] if r})
            assignments.append(dict(gene_id=gene, cls=cls, source='experimental_go',
                                    tier='observed-direct', probability=0.95,
                                    detail='; '.join(f'{t} {names.get(t, "?")}' for t in sorted(terms)),
                                    reference='; '.join(refs), n_hits=len(terms)))
    for row in up[up.gene_id != ''].itertuples():
        text = f'{row.protein_name}. {row.function_text}'
        for cls, pats in classes_from_text(text, row.gene_id).items():
            assignments.append(dict(gene_id=row.gene_id, cls=cls, source='uniprot_reviewed',
                                    tier='observed-curated', probability=0.90,
                                    detail=f'{row.accession} ({row.uniprot_id}) matched: ' + ', '.join(pats),
                                    reference='; '.join(f'PMID:{p}' for p in row.function_pmids.split('|') if p),
                                    n_hits=len(pats)))
    screen = nodes[nodes.screen_any_phenotype.notna()]
    for column, cls in SCREEN_RULES.items():
        for row in screen[screen[column] == True].itertuples():   # noqa: E712
            assignments.append(dict(gene_id=row.gene_id, cls=cls, source='phenotype_screen',
                                    tier='observed-screen', probability=0.80,
                                    detail=f'{column} = True (arrayed splitCas9 imaging screen)',
                                    reference='PMID:35538310', n_hits=1))
    assign = pd.DataFrame(assignments)
    print(f'{len(assign)} gene x class x source assignments over {assign.gene_id.nunique()} genes')
    assign.groupby(['source', 'tier']).agg(assignments=('cls', 'size'), genes=('gene_id', 'nunique'))
""")

nb.md('### What the two guards threw away',
      'Every match discarded -- because a negation cue preceded it in the same sentence, or because the '
      'sentence was about the host and the class was one the host also has -- is listed, so the guards are '
      'auditable both ways: a wrong suppression is as visible as a wrong assignment.')
cell(nb, """
    negated = pd.DataFrame(suppressed)
    negated.to_csv(OUT / 'negated_matches.tsv', sep='\\t', index=False)
    print(f'{len(negated)} text matches discarded over '
          f'{negated.gene_id.nunique() if len(negated) else 0} genes: '
          f'{dict(negated.guard.value_counts()) if len(negated) else {}}')
    negated.head(25)
""")

nb.md('### Essentiality is not a function',
      'The CRISPR fitness columns were considered as a fourth source and **rejected as a class source**. A fitness '
      'score says a gene is needed for growth; it does not say what the protein does, and none of the 15 classes '
      'is "required for growth". Mapping a low fitness score onto a function class would manufacture labels. The '
      'fitness value is carried into the evidence table as context for the review queue and nothing else.')
cell(nb, """
    FITNESS = 'fit_invitro_hff'
    fitness = nodes[['gene_id', FITNESS]].dropna()
    print(f'{len(fitness):,} genes have {FITNESS}; recorded as context, used for no class. '
          f'{(fitness[FITNESS] < -2).sum():,} are below -2 (needed for growth in fibroblasts).')
""")

nb.md('## 8 · Per-gene truth, multi-label by design',
      'A gene keeps every class it genuinely fits, because proteins do several things: ROP18 is a kinase that '
      'suppresses host immunity, and both halves are true. The gene-level `tier` is the strongest tier among its '
      'classes, with the deciding source and reference taken from that tier.')
cell(nb, """
    TIER_RANK = {'observed-direct': 0, 'observed-curated': 1, 'observed-screen': 2}
    truth_rows = []
    for gene, rows in assign.groupby('gene_id'):
        rows = rows.assign(rank=rows.tier.map(TIER_RANK)).sort_values(['rank', 'n_hits'],
                                                                      ascending=[True, False])
        best = rows.iloc[0]
        classes = sorted(set(rows.cls))
        truth_rows.append(dict(gene_id=gene, classes='|'.join(classes), n_classes=len(classes),
                               multi_label=len(classes) > 1, tier=best.tier, probability=best.probability,
                               deciding_source=best.source, reference=best.reference,
                               deciding_class=best.cls, deciding_detail=best.detail,
                               sources='|'.join(sorted(set(rows.source)))))
    truth = pd.DataFrame(truth_rows).sort_values(['tier', 'gene_id'])
    # Several reviewed entries (different strains of the same gene) can reach one ME49 id, so the
    # names are collapsed per gene before the merge; otherwise the truth table would gain rows.
    named = (up[up.gene_id != ''].groupby('gene_id')
             .agg(protein_name=('protein_name', lambda s: ' / '.join(sorted({x for x in s if x}))),
                  gene_names=('gene_names', lambda s: ' / '.join(sorted({x for x in s if x}))),
                  accessions=('accession', lambda s: '|'.join(sorted(set(s))))).reset_index())
    truth = truth.merge(named, on='gene_id', how='left')
    truth[['protein_name', 'gene_names', 'accessions']] = truth[
        ['protein_name', 'gene_names', 'accessions']].fillna('')
    print(f'{len(truth)} genes carry at least one function class; '
          f'{int(truth.multi_label.sum())} are multi-label ({truth.multi_label.mean():.0%})')
    truth.tier.value_counts().to_frame('genes')
""")

nb.md('## 9 · Information content of each label',
      'A class that half the labelled genes carry is not much of an answer. `ic_corpus` makes that explicit, and '
      'for the fine GO labels `ic_structural` says how specific the term is in the ontology regardless of our set. '
      'Any scoring of the function column should weight a correct call by the information content of the label it '
      'got right, exactly so that "protein binding" does not count as much as "actin filament severing activity".')
cell(nb, """
    long = assign[['gene_id', 'cls', 'source', 'tier', 'probability', 'detail', 'reference']].copy()
    n_labelled = truth.gene_id.nunique()
    class_ic = []
    for cls in CLASSES:
        genes = set(assign.loc[assign.cls == cls, 'gene_id'])
        p = len(genes) / n_labelled if n_labelled else float('nan')
        class_ic.append(dict(label=cls, kind='coarse_class', genes=len(genes), prevalence=round(p, 4),
                             ic_corpus=round(-math.log2(p), 2) if genes else float('nan'),
                             ic_structural=float('nan')))
    term_genes = collections.defaultdict(set)
    for row in go[go.gene_id != ''].itertuples():
        term_genes[(row.go_id, row.aspect)].add(row.gene_id)
    for (go_id, aspect), genes in term_genes.items():
        p = len(genes) / n_labelled
        class_ic.append(dict(label=f'{go_id} {names.get(go_id, "?")}', kind=f'go_{aspect}', genes=len(genes),
                             prevalence=round(p, 4), ic_corpus=round(-math.log2(p), 2),
                             ic_structural=round(ic_structural(go_id), 2)))
    ic = pd.DataFrame(class_ic).sort_values(['kind', 'ic_corpus'])
    ic.to_csv(OUT / 'label_information_content.tsv', sep='\\t', index=False)
    print(f'{len(ic)} labels scored for information content '
          f'({len(CLASSES)} coarse classes, {len(term_genes)} fine GO terms)')
    ic[ic.kind == 'coarse_class'].sort_values('ic_corpus')
""")
cell(nb, """
    vague = ic[(ic.kind != 'coarse_class')].sort_values('ic_structural').head(8)
    print('the least specific fine labels in the set (structural IC near the root):')
    vague[['label', 'genes', 'ic_structural', 'ic_corpus']]
""")

nb.md('## 10 · Verifying the PubMed ids rather than trusting the field',
      'Both QuickGO and UniProt hand us PubMed ids in a field; a reference that does not exist, or that is a '
      'different paper, would quietly make a row unverifiable. A deterministic sample is checked against NCBI '
      'E-utilities `esummary` and the raw response is pinned beside the other sources.')
cell(nb, """
    pmids = set()
    for ref in go.loc[go.gene_id != '', 'reference']:
        if str(ref).startswith('PMID:'):
            pmids.add(str(ref).split(':', 1)[1])
    for field in up.loc[up.gene_id != '', 'function_pmids']:
        pmids.update(p for p in str(field).split('|') if p)
    pmids.add('35538310')   # the arrayed screen behind every observed-screen row
    pmids = sorted(pmids)
    rng = np.random.default_rng(72)
    sample = sorted(rng.choice(pmids, size=min(30, len(pmids)), replace=False).tolist())
    url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id='
           + ','.join(sample))
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    with urllib.request.urlopen(urllib.request.Request(
            url, headers={'User-Agent': 'ParaScape/function-truth (starplast)'}), timeout=120) as response:
        body = response.read()
    (SRC / 'pubmed_verification.json').write_bytes(body)
    manifest.append(dict(url=url, status=200, retrieved_utc=stamp, bytes=len(body),
                         sha256=hashlib.sha256(body).hexdigest(), file='pubmed_verification.json',
                         source='ncbi_eutils_esummary', sampled=len(sample), pmids_total=len(pmids)))
    (SRC / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\\n')
    result = json.loads(body)['result']
    checked = pd.DataFrame([dict(pmid=p, found='error' not in result.get(p, {'error': 1}),
                                 journal=result.get(p, {}).get('source', ''),
                                 year=str(result.get(p, {}).get('pubdate', ''))[:4],
                                 title=result.get(p, {}).get('title', '')[:90]) for p in sample])
    checked.to_csv(OUT / 'pubmed_verification.tsv', sep='\\t', index=False)
    print(f'{len(pmids)} distinct PubMed ids cited by the set; {len(sample)} sampled (seed 72), '
          f'{int(checked.found.sum())} resolve at NCBI, {int((~checked.found).sum())} do not')
    checked
""")

nb.md('## 11 · The review queue, one row per gene',
      'Einar reviewed the host-facing candidates per claim and said the duplication was tedious, so this queue is '
      '**one row per gene** with its whole evidence folded into the row. The score is the sum of explicit '
      'penalties, worst first, so the top of the list is where a wrong class is most likely:',
      '* **+3** the only evidence is a screen phenotype (a process claim with no molecular statement);',
      '* **+3** the gene reaches 4 or more classes -- usually a broad text pattern (`catalyz`, `transport`) firing '
      'on a sentence about something else;',
      '* **+2** every class came from a single text pattern match, with no GO support;',
      '* **+2** the gene has experimental GO but all of it is uninformative (structural IC below 4, e.g. only '
      '`protein binding`), so nothing specific is actually known here;',
      '* **+2** the accession route disagreed or was ambiguous;',
      '* **+1** no PubMed id travels with the deciding evidence;',
      '* **+1** screen-only and the two scorers of the screen disagreed on this gene.')
cell(nb, """
    term_ic = {go_id: ic_structural(go_id) for (go_id, _) in term_genes}
    go_by_gene = go[go.gene_id != ''].groupby('gene_id')
    scorers = nodes.set_index('gene_id').screen_scorers_agree.to_dict()
    route_by_gene = {}
    for frame in (up, go):
        for row in frame.itertuples():
            if getattr(row, 'gene_id', ''):
                route_by_gene.setdefault(row.gene_id, row.map_route)

    queue = []
    for row in truth.itertuples():
        rows = assign[assign.gene_id == row.gene_id]
        score, why = 0, []
        if set(rows.source) == {'phenotype_screen'}:
            score += 3; why.append('screen phenotype only: a process claim, no molecular function')
        if row.n_classes >= 4:
            score += 3; why.append(f'{row.n_classes} classes -- check for a broad pattern firing wrongly')
        text_only = rows[rows.source == 'uniprot_reviewed']
        if not rows.source.eq('experimental_go').any() and len(text_only) and text_only.n_hits.max() == 1:
            score += 2; why.append('every class rests on a single text pattern')
        mine = go_by_gene.get_group(row.gene_id) if row.gene_id in go_by_gene.groups else None
        if mine is not None and len(mine):
            best_ic = max(term_ic.get(t, 0) for t in mine.go_id)
            if best_ic < 4:
                score += 2; why.append(f'experimental GO is all vague (best structural IC {best_ic:.1f})')
        route = route_by_gene.get(row.gene_id, '')
        if 'disagree' in route:
            score += 2; why.append('accession mapping routes disagree')
        if route.startswith('D '):
            score += 2; why.append('gene id inferred from the ToxoDB strain-suffix convention, '
                                   'not from a cross-reference')
        if not str(row.reference).strip():
            score += 1; why.append('no PubMed id on the deciding evidence')
        if set(rows.source) == {'phenotype_screen'} and scorers.get(row.gene_id) is False:
            score += 1; why.append('screen scorers disagreed on this gene')
        queue.append(dict(gene_id=row.gene_id, protein=row.protein_name or row.gene_names or '',
                          classes=row.classes, n_classes=row.n_classes, tier=row.tier,
                          uncertainty=score, reasons='; '.join(why) or 'none',
                          evidence=' || '.join(f'[{r.source}] {r.detail}' for r in rows.itertuples())[:1200],
                          reference=row.reference,
                          question='Are these the right function classes for this gene?'))
    # Genes with real experimental evidence that no rule could turn into a class -- usually because the
    # only annotation is something like "protein binding". They have fine-grained truth and no coarse
    # label, which is a rule gap rather than an absence of knowledge, so they head the queue.
    no_class = sorted((set(go.loc[go.gene_id != '', 'gene_id']) |
                       set(up.loc[up.gene_id != '', 'gene_id'])) - set(truth.gene_id))
    for gene in no_class:
        bits = [f'[experimental_go] {r.go_id} {names.get(r.go_id, "?")} ({r.code})'
                for r in go[go.gene_id == gene].itertuples()]
        bits += [f'[uniprot_reviewed] {r.protein_name}. {r.function_text}'
                 for r in up[up.gene_id == gene].itertuples()]
        queue.append(dict(gene_id=gene, protein='', classes='', n_classes=0, tier='no class assigned',
                          uncertainty=5,
                          reasons='has experimental or curated function evidence that no class rule '
                                  'matched -- either the evidence is too vague to classify or a rule is missing',
                          evidence=' || '.join(bits)[:1200], reference='',
                          question='Which class does this gene belong to, or is the evidence too vague?'))
    print(f'{len(no_class)} genes have function evidence that matched no class rule')
    review = pd.DataFrame(queue).sort_values(['uncertainty', 'n_classes'], ascending=[False, False])
    review.to_csv(OUT / 'review_queue.tsv', sep='\\t', index=False)
    print(f'{len(review)} genes in the queue, {int((review.uncertainty > 0).sum())} with at least one '
          f'penalty; the top 100 are what Einar would rule on first')
    review.uncertainty.value_counts().sort_index(ascending=False).to_frame('genes')
""")

nb.md('## 12 · Tables written',
      '`function_truth.tsv` (per gene), `function_evidence.tsv` (long form: every gene x source x term/text with '
      'its evidence code and PMID), `coverage.tsv` (genes per class per tier), `class_rules.md` (the rules '
      'verbatim), `label_information_content.tsv`, `review_queue.tsv`, `unmapped_accessions.tsv`, '
      '`pubmed_verification.tsv`, `evidence.parquet` (information-map rows) and `sources/`.')
cell(nb, """
    evidence_rows = []
    for row in go[go.gene_id != ''].itertuples():
        evidence_rows.append(dict(gene_id=row.gene_id, source='experimental_go', aspect=row.aspect,
                                  term_or_text=f'{row.go_id} {names.get(row.go_id, "?")}',
                                  qualifier=row.qualifier, evidence_code=f'{row.code} ({row.eco})',
                                  pmid=row.reference, tier='observed-direct', probability=0.95,
                                  ic_structural=round(ic_structural(row.go_id), 2),
                                  accession=row.accession, map_route=row.map_route,
                                  assigned_by=row.assigned_by))
    for row in up[up.gene_id != ''].itertuples():
        if row.function_text:
            evidence_rows.append(dict(gene_id=row.gene_id, source='uniprot_reviewed', aspect='FUNCTION',
                                      term_or_text=row.function_text, qualifier='',
                                      evidence_code=row.function_eco, pmid=row.function_pmids,
                                      tier='observed-curated', probability=0.90,
                                      ic_structural=float('nan'), accession=row.accession,
                                      map_route=row.map_route, assigned_by='UniProt curator'))
    for column, cls in SCREEN_RULES.items():
        for row in screen[screen[column] == True].itertuples():   # noqa: E712
            evidence_rows.append(dict(gene_id=row.gene_id, source='phenotype_screen', aspect='P',
                                      term_or_text=f'{column} -> required for {cls}', qualifier='involved_in',
                                      evidence_code='arrayed splitCas9 imaging screen, scored by eye',
                                      pmid='PMID:35538310', tier='observed-screen', probability=0.80,
                                      ic_structural=float('nan'), accession='', map_route='n/a',
                                      assigned_by='Smith et al. 2022'))
    long_evidence = pd.DataFrame(evidence_rows)
    long_evidence.to_csv(OUT / 'function_evidence.tsv', sep='\\t', index=False)
    # Context that is recorded but assigns no class: UniProt names/keywords/domains (keywords and
    # domains can come from automatic annotation) and the CRISPR fitness scores (essentiality).
    context = [dict(gene_id=r.gene_id, source='uniprot_reviewed_context',
                    term_or_text=f'name={r.protein_name}; keywords={r.keywords}; domains={r.domains}',
                    note='keywords and domains may be automatic annotation; not used for class assignment',
                    pmid=r.lit_pmids) for r in up[up.gene_id != ''].itertuples()]
    context += [dict(gene_id=r.gene_id, source='crispr_fitness_context',
                     term_or_text=f'{FITNESS}={getattr(r, FITNESS):.3f}',
                     note='essentiality is not a function class', pmid='')
                for r in fitness.itertuples()]
    pd.DataFrame(context).to_csv(OUT / 'context_not_used_for_classes.tsv', sep='\\t', index=False)
    truth[['gene_id', 'protein_name', 'gene_names', 'classes', 'n_classes', 'multi_label', 'tier',
           'probability', 'deciding_source', 'deciding_class', 'reference', 'deciding_detail',
           'sources']].to_csv(OUT / 'function_truth.tsv', sep='\\t', index=False)
    assign.to_csv(OUT / 'function_classes_long.tsv', sep='\\t', index=False)
    coverage = (assign.groupby(['cls', 'tier']).gene_id.nunique().unstack(fill_value=0)
                .reindex(CLASSES).fillna(0).astype(int))
    coverage['genes_any_tier'] = [assign.loc[assign.cls == c, 'gene_id'].nunique() for c in coverage.index]
    coverage['ic_corpus'] = [ic.loc[ic.label == c, 'ic_corpus'].iloc[0] for c in coverage.index]
    coverage.to_csv(OUT / 'coverage.tsv', sep='\\t')
    print(f'{len(long_evidence):,} evidence rows written')
    coverage
""")

nb.md('The information-map rows use the same shape as the host-facing truth set, so the function column can be '
      'added with one `register(Trait(...))` line and these evidence rows.')
cell(nb, """
    map_rows = [dict(gene_id=r.gene_id, trait='function', value=r.classes, probability=r.probability,
                     source_tier=r.tier, reference=r.reference,
                     dataset=f'function truth 2026-11-04 ({r.deciding_source}; {r.deciding_class})')
                for r in truth.itertuples()]
    pd.DataFrame(map_rows).to_parquet(OUT / 'evidence.parquet', index=False)
    print(f'{len(map_rows)} information-map rows for trait "function"')
""")

nb.md('## 13 · Is this enough to score function prediction at 90 % precision?',
      'The honest arithmetic, printed rather than asserted: per class, how many positives exist at the two '
      'non-screen tiers, and how many held-out calls a Clopper-Pearson lower 90 % bound would need to reach 0.90 '
      'if the true precision were 0.95.')
cell(nb, """
    from scipy.stats import beta
    def lower_bound(correct, total, alpha=0.10):
        return 0.0 if correct == 0 else beta.ppf(alpha, correct, total - correct + 1)
    need = {}
    for n in range(5, 400):
        k = round(0.95 * n)
        if lower_bound(k, n) >= 0.90:
            need['at true precision 0.95'] = n
            break
    for n in range(5, 2000):
        k = round(0.92 * n)
        if lower_bound(k, n) >= 0.90:
            need['at true precision 0.92'] = n
            break
    strong = assign[assign.tier != 'observed-screen']
    per_class = (strong.groupby('cls').gene_id.nunique().reindex(CLASSES).fillna(0).astype(int)
                 .sort_values(ascending=False).to_frame('non_screen_genes'))
    needed = need['at true precision 0.95']
    per_class['enough_alone'] = per_class.non_screen_genes >= needed
    print('held-out calls a lower 90% bound needs to reach 0.90:', need)
    print(f'{strong.gene_id.nunique()} genes at observed-direct/observed-curated; '
          f'{int(per_class["enough_alone"].sum())} of {len(CLASSES)} classes have, on their own, the '
          f'{needed} held-out calls a 0.90 lower bound needs even at a true precision of 0.95')
    per_class
""")

nb.md('With five folds the held-out calls per fold are a fifth of those counts, so **no single class can '
      'demonstrate a 0.90 lower bound today**. The one question the set can support is the pooled, '
      'multi-label one: over all classes at once, of the confident class calls on held-out genes, what '
      'fraction is right. That has a few hundred calls, which is enough to *measure* a bound -- and the '
      'bound will only clear 0.90 if the true precision is well above it.')

nb.md('## 14 · The rules and the README, written verbatim',
      '`class_rules.md` contains the GO seeds with their resolved ids and the regular expressions exactly as '
      'applied, so a reader can reproduce any class on any gene by hand.')
cell(nb, """
    lines = ['# Function class rules (2026-11-04), verbatim', '',
             'Applied by `scripts/build_function_truth_notebook.py`. Two rule kinds per class.', '',
             '* **GO rule**: a gene matches the class if one of its *experimental* GO terms is the seed term or '
             'a descendant of it over `is_a`/`part_of` in the `go-basic.obo` snapshot of 2026-11-04. '
             'A class reached this way is tier `observed-direct`, probability 0.95.',
             '* **Text rule**: a Python regular expression (`re.I`) searched in '
             '`"<recommended protein name>. <FUNCTION comment>"` of a reviewed UniProt entry. A class reached '
             'only this way is tier `observed-curated`, probability 0.90.',
             '* **Screen rule**: a True phenotype column of the arrayed imaging screen supports one process '
             'class and nothing else. Tier `observed-screen`, probability 0.80.', '',
             'Never used to assign a class: InterPro/Pfam domains, ToxoDB product names, IEA/ISS/ISO GO '
             'annotations, UniProt keywords and domain features, and CRISPR fitness scores.', '']
    for cls in CLASSES:
        rule = CLASS_RULES[cls]
        lines += [f'## `{cls}`', '', 'GO seeds:', '']
        if rule['go']:
            lines += [f'* `{go_id}` {names[go_id]} ({namespace[go_id]}, '
                      f'{len(closure(go_id, children))} descendants, structural IC '
                      f'{ic_structural(go_id):.2f})'
                      for _, go_id in seeds.get(cls, [])]
        else:
            lines += ['* none -- this class is text-derived only (the matching GO terms are '
                      'cellular-component terms, which are a different trait of the map)']
        lines += ['', 'Text patterns:', '']
        lines += [f'* `{p}`' for p in rule['text']]
        screen_cols = [c for c, k in SCREEN_RULES.items() if k == cls]
        lines += ['', 'Screen columns: ' + (', '.join(f'`{c}`' for c in screen_cols) or 'none'), '']
    (OUT / 'class_rules.md').write_text('\\n'.join(lines) + '\\n')
    print(f'class_rules.md: {len(lines)} lines, {len(CLASSES)} classes')
""")

cell(nb, """
    counts = truth.tier.value_counts().to_dict()
    top_classes = coverage.sort_values('genes_any_tier', ascending=False)
    readme = f'''# The function truth set (2026-11-04)

Built by `scripts/build_function_truth_notebook.py`; raw responses fetched by
`scripts/fetch_function_sources.py` and pinned in `sources/` with their retrieval times and sha256.
The executed notebook is `notebooks/function_truth_2026_11_04.ipynb`.

Instruction 72, target trait 3. The map previously held experimental GO function annotations for about
25 genes. This set labels **{len(truth)} genes** with at least one coarse function class.

| tier | probability | genes | what it is |
|---|---|---|---|
| `observed-direct` | 0.95 | {counts.get('observed-direct', 0)} | experimental GO (EXP/IDA/IPI/IMP/IGI/IEP), BP + MF, QuickGO taxon 5811 |
| `observed-curated` | 0.90 | {counts.get('observed-curated', 0)} | reviewed UniProt (Swiss-Prot) curator `FUNCTION` comment |
| `observed-screen` | 0.80 | {counts.get('observed-screen', 0)} | arrayed splitCas9 imaging screen phenotype (PMID 35538310); a process claim only |

**These three probabilities are provisional tier defaults** (instruction 72 v0), recorded as such. They
are not measured error rates. v1 replaces each one by how often that source is confirmed by an
independent direct source on the genes both cover, with a confidence interval, and no source keeps a
probability above what its replication supports.

{int(truth.multi_label.sum())} of {len(truth)} genes carry more than one class ({truth.multi_label.mean():.0%}),
because proteins do several things; the classes are a set, never one winner.

## Excluded, and this is why the set is small
InterPro/Pfam domain assignments, ToxoDB `product` names, and every GO annotation with an IEA, ISS, ISO
or other electronic/similarity evidence code are **not used**, for either labels or classes. They are
predictions from sequence -- the thing the information map is supposed to predict. Scoring against them
would measure agreement with a predictor rather than with biology (instruction 70), which is the
failure both host-facing scoring runs diagnosed: run 1 cleared the bar only because its strongest
features were downstream of the conclusion. UniProt keywords and domain features and the CRISPR fitness
scores are recorded in `context_not_used_for_classes.tsv` and assign no class either; essentiality is
not a function.

## The label space
{len(CLASSES)} coarse classes, derived from the curated text and GO terms by the rules printed verbatim in
`class_rules.md`. The coarse classes are a **scoring convenience**; the truth is the fine-grained
evidence in `function_evidence.tsv` (GO terms with their codes and PMIDs, and the curator function text).
Information content is reported per label in `label_information_content.tsv`: `ic_corpus` for how common
a label is in this set and `ic_structural` for how specific a GO term is in the ontology, so that a vague
answer ("protein binding", structural IC {ic.loc[ic.label.str.startswith('GO:0005515'), 'ic_structural'].iloc[0] if (ic.label.str.startswith('GO:0005515')).any() else float('nan')})
cannot score like a specific one.

Largest classes: {', '.join(f'{i} ({int(r.genes_any_tier)})' for i, r in top_classes.head(6).iterrows())}.

## Two things a reader should know before trusting a row
* **Gene ids.** Many reviewed entries are the GT1, VEG or RH88 copy of a gene and carry no TGME49
  cross-reference. Those are placed on ME49 through the ToxoDB convention that the numeric suffix of a
  gene id is shared across strains (`TGGT1_254370` -> `TGME49_254370`). That is an inference from an id
  convention, labelled `D (ToxoDB strain-suffix, inferred)` in every table and penalised in the review
  queue. 1 reviewed entry and 1 accession remain unplaced or ambiguous (`unmapped_accessions.tsv`).
* **Two guards on the text rules**, both auditable in `negated_matches.tsv`. A match is discarded when a
  negation cue precedes it in the same sentence (ASP1 is "dispensable for ... bradyzoite viability", so it
  is not a differentiation gene), and a match of `transcription_chromatin`, `translation` or
  `proteostasis` is discarded when its sentence names the host (ROP17 "Promotes autophagy in human cells"
  is host modulation, which it already has, not parasite proteostasis).

Genes whose only function evidence is uninformative -- an experimental `protein binding` annotation, or a
reviewed entry with no FUNCTION comment at all (RON4, MIC21, the enolases, GRA1) -- get no class and head
the review queue, because that is a gap in the rules or in the evidence, not an absence of knowledge.

## Is this enough to score function prediction at 90 % precision?
Not per class. A Clopper-Pearson lower 90 % bound needs {needed} held-out calls to reach 0.90 even if the
true precision is 0.95, and {int(need['at true precision 0.92'])} if it is 0.92. At the two non-screen tiers
{strong.gene_id.nunique()} genes are labelled in total and {int(per_class['enough_alone'].sum())} of
{len(CLASSES)} classes reach that count on their own -- before splitting into folds. What the set *can*
support is the pooled multi-label question (all classes at once, a few hundred held-out class calls), plus
the CAFA no-knowledge rule: genes with no prior annotation scored separately.

## Files
* `function_truth.tsv` -- one row per gene: classes, tier, probability, deciding source, reference.
* `function_evidence.tsv` -- long form: gene x source x GO term or curator text, evidence code, PMID.
* `function_classes_long.tsv` -- gene x class x source, the per-class tiers behind the per-gene row.
* `class_rules.md` -- the rules verbatim (GO seeds with resolved ids, regular expressions, screen columns).
* `coverage.tsv` -- genes per class per tier, with each class information content.
* `label_information_content.tsv` -- per-label specificity, coarse classes and fine GO terms.
* `review_queue.tsv` -- **one row per gene**, worst-first, for Einar's ruling on the top 100.
* `unmapped_accessions.tsv` -- reviewed entries with no single ME49 gene; curated function the map cannot place.
* `pubmed_verification.tsv` -- a 30-id sample of the cited PMIDs checked against NCBI E-utilities.
* `negated_matches.tsv` -- every text match the two guards discarded, with the cue and the sentence.
* `context_not_used_for_classes.tsv` -- UniProt keywords/domains and CRISPR fitness, recorded only.
* `evidence.parquet` -- information-map rows for the `function` trait.
* `sources/` -- every raw API response with URL, status, retrieval time and sha256 (`manifest.json`).
'''
    (OUT / 'README.md').write_text(readme)
    print(readme[:1500])
""")

nb.write('notebooks/function_truth_2026_11_04.ipynb')
print('notebook written')
