"""Instruction 72, sequence and structure (3Di) tier-2 traits of the ParaScape information map.

Writes notebooks/sequence_traits_2026_11_04.ipynb (executed) and results/sequence_traits_2026_11_04/.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
from notebook_runner import ExecutedNotebook  # noqa: E402

nb = ExecutedNotebook('Sequence and 3Di traits of the information map')
nb.md('# Sequence and structure (3Di) traits',
      'Instruction 72, tier 2 (Einar, 2026-11-04). One row per Toxoplasma gene, giving the protein sequence and its '
      'Foldseek 3Di structural-alphabet string cut into four regions: the whole protein, the mature protein (signal '
      'peptide removed), the stretch before the first transmembrane (TM) helix, and the stretch after the last one. '
      'The point is that only a region facing outward can reach the host, so the flanks are kept separately and each '
      'records which side of the membrane DeepTMHMM puts it on.',
      '**The no-TM rule.** A protein with no TM domain has no membrane to split it, so BOTH flank columns are the whole '
      'mature sequence and are IDENTICAL (`flanks_identical` is True). This is intended, and asserted below.',
      '**Known gap.** Propeptide / TEXEL cleavage (e.g. of rhoptry and microneme proteins) is NOT modelled: there is no '
      'reliable per-gene table here. `seq_mature` removes only the DeepTMHMM signal peptide.')

nb.md('## 1 · Inputs')
nb.code('import json, os',
        'import numpy as np, pandas as pd',
        'from pathlib import Path',
        "OUT = Path('results/sequence_traits_2026_11_04'); OUT.mkdir(parents=True, exist_ok=True)",
        "nodes = pd.read_parquet('parascape/data/nodes.parquet', columns=['gene_id', 'sequence']).set_index('gene_id')",
        "dt = pd.read_csv('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/deeptmhmm_toxo.csv', low_memory=False)",
        "print(len(nodes), 'genes in nodes.parquet;', len(dt), 'DeepTMHMM rows,', dt.identifier.nunique(), 'unique identifiers')",
        "dup = dt[dt.identifier.duplicated(keep=False)].identifier.unique()",
        "print('duplicated identifiers (first row kept):', list(dup))",
        "N_DT_ROWS = len(dt); dt = dt.drop_duplicates('identifier').set_index('identifier')",
        "print('nodes genes without a DeepTMHMM row:', len(nodes.index.difference(dt.index)))",
        "print('topology string alphabet:', sorted(set(''.join(dt.topology_string.dropna().iloc[:500]))))")

nb.md('## 2 · Check 1: the DeepTMHMM coordinate convention',
      'The spec assumes 1-based inclusive coordinates. Rather than trust that, compare `sequence[start-1:end]` with the '
      'recorded `tm_<i>_sequence` and `sp_sequence` over every domain of every gene, and also compare the recorded '
      'DeepTMHMM sequence with the nodes sequence.')
nb.code('same_seq = (nodes.sequence == dt.sequence.reindex(nodes.index))',
        "print('nodes sequence == DeepTMHMM sequence:', int(same_seq.sum()), 'of', len(nodes))",
        'tot = ok = bad = 0; bad_examples = []',
        'for g, r in dt.iterrows():',
        '    s = r.sequence',
        "    pairs = [('sp', r.sp_start, r.sp_end, r.sp_sequence)] if pd.notna(r.sp_start) else []",
        "    pairs += [(f'tm_{i}', r[f'tm_{i}_start'], r[f'tm_{i}_end'], r[f'tm_{i}_sequence']) for i in range(1, 25) if pd.notna(r[f'tm_{i}_start'])]",
        '    for name, a, b, rec in pairs:',
        '        tot += 1',
        '        if s[int(a) - 1:int(b)] == rec: ok += 1',
        '        else:',
        '            bad += 1; bad_examples.append((g, name, int(a), int(b)))',
        "print(f'1-based inclusive sequence[start-1:end] == recorded: {ok} of {tot} domains; mismatches {bad}', bad_examples[:5])",
        "assert bad == 0, 'convention disagrees: work out the real one'",
        "CONVENTION = f'1-based inclusive confirmed on {ok} of {tot} domains (signal peptides and TM helices, all genes)'",
        'CONVENTION')

nb.md('## 3 · Reading the 3Di databases',
      'Reused from `af3_folding/analysis_04_3di.py` (`read_db`): the `.lookup` file maps record index to name and `_ss` '
      'holds the NUL-separated 3Di strings. Only TGME49 entries are kept (the databases hold other species too). '
      '`naive_db` names are bare ids; `cleaned_db` names carry suffixes such as `__sol` or `__ecto_24-288`, which are '
      'stripped. Some cleaned entries are domain fragments of a protein, so several entries can map to one gene; the '
      'entry whose length equals the annotated sequence is chosen. Naive is preferred, then cleaned. An AlphaFold model '
      'can cover a different length from the annotated protein; those genes are reported and excluded (empty 3Di).')
nb.code('def read_db(prefix):',
        '    """name -> list of 3Di strings (a gene can have several records); logic of analysis_04_3di.read_db."""',
        '    names = {}',
        '    for ln in open(f"{prefix}.lookup"):',
        '        f = ln.rstrip("\\n").split("\\t")',
        '        if len(f) >= 2: names[int(f[0])] = f[1]',
        '    with open(f"{prefix}_ss", "rb") as fh:',
        '        recs = fh.read().split(b"\\x00")',
        '    out = {}',
        '    for i, r in enumerate(recs):',
        '        s = r.decode("ascii", "ignore").strip()',
        '        if not s: continue',
        '        n = names.get(i)',
        '        if n is None: continue',
        '        n = n[4:] if n.startswith("AF3-") else n',
        '        if n.endswith(".cif"): n = n[:-4]',
        '        n = n.split("__")[0]',
        '        if n.startswith("TGME49_"): out.setdefault(n, []).append(s)',
        '    return out',
        "db = {'naive_db': read_db('/home/carruthers/foldseek/naive_db'), 'cleaned_db': read_db('/home/carruthers/foldseek/cleaned_db')}",
        "for k, v in db.items(): print(k, len(v), 'TGME49 genes,', sum(len(x) for x in v.values()), 'records')")
nb.code('tdi, tdi_source, mismatch = {}, {}, {}',
        'for g, seq in nodes.sequence.items():',
        "    for src in ('naive_db', 'cleaned_db'):",
        '        cands = db[src].get(g, [])',
        '        hit = [c for c in cands if len(c) == len(seq)]',
        '        if hit:',
        '            tdi[g], tdi_source[g] = hit[0], src; break',
        '        if cands: mismatch.setdefault(g, []).append((src, [len(c) for c in cands]))',
        "no_db = [g for g in nodes.index if g not in tdi and g not in mismatch]",
        "mm = [g for g in nodes.index if g not in tdi and g in mismatch]",
        "print('3Di length matches protein length:', len(tdi))",
        "print('in a database but every record has a different length (excluded):', len(mm))",
        "print('in neither database:', len(no_db))",
        "print('examples of length mismatch (gene, [(db, record lengths)], protein length):')",
        "for g in mm[:5]: print(' ', g, mismatch[g], len(nodes.sequence[g]))",
        "print('used cleaned_db although naive_db had a record of the wrong length:', sum(1 for g, s in tdi_source.items() if s == 'cleaned_db' and g in mismatch))")


nb.md('## 3b · Filling the gaps from AlphaFold3 models already on disk',
      '**Source.** `/media/carruthers/mnt1/af3/af_output/` holds AlphaFold3 outputs. Folders named `af3_<UniProt>` are '
      'whole-protein models; `af3_<UniProt>_f<i>of<n>_<start>-<end>` are FRAGMENTS of one long protein, the name giving the '
      '1-based inclusive residue range. (The same tree also holds interaction jobs named `bait_`, `pep_`, `xl_` ... that carry '
      'TGME49 ids but are domains or complexes, not whole proteins; they are NOT used.) UniProt accessions are mapped to '
      'genes with the `accession` column of the DeepTMHMM table; an accession that maps to more than one gene is attached to '
      'none and listed below.',
      '**Only genes lacking a database 3Di are converted**: the Foldseek databases above stay the first choice wherever their '
      'length already matches. `foldseek createdb` is run on copies of the `model.cif` files in batches of at most 2,000 '
      'in a scratch directory on mnt1 (outside the repo and the NAS), the `_ss` file is read with the same NUL-separated '
      'logic, and each batch is deleted afterwards. The AF3 tree itself is never written to.',
      '**Fragment stitching.** Each fragment 3Di is placed at its residue range inside a string of the full protein length. '
      'Residues covered by no fragment get the gap character `X` (not a 3Di letter). Where fragments overlap, the fragment '
      'with the higher mean pLDDT (mean of `atom_plddts` in `af3_<ID>_confidences.json`) keeps the residue; if no score can be '
      'read, the earlier fragment (lower start) wins. `tdi_coverage` is the fraction of residues with a real letter and '
      '`tdi_n_fragments` the number of models used.',
      '**Acceptance rule.** A 3Di string is accepted only when its length equals the protein length in `nodes.parquet` (for a '
      'whole model) or when every fragment 3Di length equals its stated range and the ranges fit inside the protein (for a '
      'stitched string, which is then exactly the protein length by construction). Otherwise the gene goes to the mismatch '
      'list with both lengths and its 3Di stays empty: nothing is padded or trimmed to fit.')
nb.code('import re, json, shutil, subprocess, collections, datetime',
        "AF3 = Path('/media/carruthers/mnt1/af3/af_output')",
        "folders = sorted(x for x in os.listdir(AF3) if x.startswith('af3_'))",
        "pat = re.compile(r'^af3_(.+?)(?:_f(\\d+)of(\\d+)_(\\d+)-(\\d+))?$')",
        'a2g = collections.defaultdict(set)',
        "for g, a in zip(dt.index, dt.accession.fillna('')):",
        "    for x in a.split(','):",
        '        if x.strip(): a2g[x.strip()].add(g)',
        'entries = collections.defaultdict(list)       # accession -> list of dict(folder, start, end, frag)',
        'n_nonmatch = 0',
        'for f in folders:',
        '    m = pat.match(f)',
        '    if not m: n_nonmatch += 1; continue',
        '    acc, i, n, s, e = m.groups()',
        '    entries[acc].append(dict(folder=f, frag=i is not None, i=int(i) if i else 0, n=int(n) if n else 1,',
        '                             start=int(s) if s else None, end=int(e) if e else None))',
        'ambiguous = sorted(a for a in entries if len(a2g.get(a, ())) > 1)',
        'unmapped = sorted(a for a in entries if a not in a2g)',
        'gene_acc = collections.defaultdict(list)',
        'for a in entries:',
        '    if len(a2g.get(a, ())) == 1: gene_acc[next(iter(a2g[a]))].append(a)',
        'gene_acc = {g: v for g, v in gene_acc.items() if g in nodes.index}',
        "print(len(folders), 'af3_* folders;', len(entries), 'distinct ids; TGME49-named whole-protein ids:', sum(a.startswith('TGME49_') for a in entries))",
        "print('accessions mapping to MORE than one gene (attached to none):', len(ambiguous), ambiguous)",
        "print('accessions not found in the DeepTMHMM accession column (not attached):', len(unmapped))",
        "print('genes with at least one AF3 folder:', len(gene_acc), '; genes with >1 accession:', sum(len(v) > 1 for v in gene_acc.values()))",
        'need = [g for g in gene_acc if g not in tdi]',
        "print('of these, genes still lacking a database 3Di (the ones AF3 must fill):', len(need))",
        'amb_report = {a: sorted(a2g[a]) for a in ambiguous}')
nb.code('# ---- cell B: run foldseek on the needed models, batches of <= 2000',
        "FS_BIN = '/home/carruthers/foldseek/foldseek/bin/foldseek'",
        "FS_VERSION = subprocess.run([FS_BIN, 'version'], capture_output=True, text=True).stdout.strip()",
        "print('foldseek version:', FS_VERSION)",
        "SCR = Path('/media/carruthers/mnt1/parascape_omics_scratch') / 'af3_3di_20261104'",
        "todo = sorted({e['folder'] for g in need for a in gene_acc[g] for e in entries[a]})",
        "print('model.cif files to convert:', len(todo))",
        'af3_3di = {}',
        'def read_ss(prefix):',
        '    names = {}',
        "    for ln in open(f'{prefix}.lookup'):",
        "        f = ln.rstrip('\\n').split('\\t')",
        '        if len(f) >= 2: names[int(f[0])] = f[1]',
        "    recs = open(f'{prefix}_ss', 'rb').read().split(b'\\x00')",
        '    out = {}',
        '    for i, r in enumerate(recs):',
        "        s = r.decode('ascii', 'ignore').strip()",
        '        if s and i in names: out.setdefault(names[i], []).append(s)',
        '    return out',
        'fs_log = []',
        'for b in range(0, len(todo), 2000):',
        '    batch = todo[b:b + 2000]',
        "    d = SCR / f'batch_{b // 2000:03d}'; (d / 'in').mkdir(parents=True, exist_ok=True)",
        '    for f in batch:',
        "        src = AF3 / f / f'{f}_model.cif'",
        "        if src.exists(): shutil.copy(src, d / 'in' / f'{f}.cif')",
        "    r = subprocess.run([FS_BIN, 'createdb', str(d / 'in'), str(d / 'db'), '--threads', '8'], capture_output=True, text=True)",
        '    fs_log.append(r.stdout[-300:] + r.stderr[-300:])',
        "    for name, ss in read_ss(str(d / 'db')).items():",
        "        name = re.sub(r'\\.cif$', '', name)",
        "        key = name if name in batch else re.sub(r'_[A-Za-z0-9]$', '', name)    # strip a chain suffix if foldseek added one",
        '        if key in batch: af3_3di.setdefault(key, []).extend(ss)',
        '    shutil.rmtree(d)                                                         # scratch removed per batch',
        "print('foldseek produced 3Di for', len(af3_3di), 'of', len(todo), 'models; scratch batches deleted')",
        'multi_chain = {k: [len(s) for s in v] for k, v in af3_3di.items() if len(v) > 1}',
        "print('models with more than one chain in the 3Di (first chain not guessed; excluded):', len(multi_chain))")
nb.code('# ---- cell C: stitching and acceptance',
        'def mean_plddt(folder):',
        '    try:',
        "        j = json.load(open(AF3 / folder / f'{folder}_confidences.json'))",
        "        return float(np.mean(j['atom_plddts']))",
        '    except Exception:',
        '        return None',
        'L_of = nodes.sequence.str.len()',
        'af3_acc, af3_src, af3_cov, af3_nfrag, af3_mismatch, af3_gaps, af3_overlap_n = {}, {}, {}, {}, {}, {}, 0',
        'for g in need:',
        '    L = int(L_of[g]); tried = []',
        '    for a in sorted(gene_acc[g]):',
        '        es = entries[a]',
        "        whole = [e for e in es if not e['frag']]",
        "        frags = sorted((e for e in es if e['frag']), key=lambda e: e['start'])",
        '        cands = []',
        '        for e in whole:',
        "            ss = af3_3di.get(e['folder'], [])",
        "            if len(ss) == 1: cands.append(('af3', [(1, len(ss[0]), ss[0], e['folder'])], ss[0], None))",
        "            elif ss: tried.append((e['folder'], 'multi-chain', [len(s) for s in ss]))",
        '        if cands:',
        '            kind, parts, ss, _ = cands[0]',
        '            if len(ss) == L:',
        "                af3_acc[g], af3_src[g], af3_cov[g], af3_nfrag[g] = ss, 'af3', 1.0, 1; break",
        "            tried.append((whole[0]['folder'], 'whole model', len(ss))); continue",
        '        if frags:',
        '            parts = []; bad = False',
        '            for e in frags:',
        "                ss = af3_3di.get(e['folder'], [])",
        "                if len(ss) != 1: bad = True; tried.append((e['folder'], 'no single-chain 3Di', len(ss))); continue",
        "                if len(ss[0]) != e['end'] - e['start'] + 1: bad = True; tried.append((e['folder'], 'fragment length != range', len(ss[0]), e['start'], e['end'])); continue",
        "                parts.append((e['start'], e['end'], ss[0], e['folder']))",
        '            if bad or not parts: continue',
        '            top = max(p[1] for p in parts)',
        "            if top > L: tried.append((a, 'fragment range beyond protein length', top, L)); continue",
        "            arr = ['X'] * L; owner = [None] * L; score = {}",
        '            for st, en, ss, f in parts:',
        '                p = mean_plddt(f) if any(st <= q[1] and en >= q[0] and q[3] != f for q in parts) else None',
        '                score[f] = p',
        '            for st, en, ss, f in parts:            # earlier-listed fragment wins ties; higher mean pLDDT wins an overlap',
        '                for k, ch in zip(range(st - 1, en), ss):',
        '                    o = owner[k]',
        '                    if o is None or (score.get(f) is not None and score.get(o) is not None and score[f] > score[o]):',
        '                        if o is not None: af3_overlap_n += 1',
        '                        arr[k], owner[k] = ch, f',
        "            s = ''.join(arr)",
        "            cov = 1 - s.count('X') / L",
        "            af3_acc[g], af3_src[g], af3_cov[g], af3_nfrag[g] = s, 'af3_fragments', cov, len(parts)",
        "            af3_gaps[g] = s.count('X'); break",
        '    else:',
        "        if g not in af3_acc: af3_mismatch[g] = tried or [('no usable model', )]",
        '# a gene whose break fired is accepted; others fall through to mismatch',
        'af3_mismatch = {g: v for g, v in af3_mismatch.items() if g not in af3_acc}',
        "print('AF3 accepted:', len(af3_acc), '| whole:', sum(v == 'af3' for v in af3_src.values()), '| stitched from fragments:', sum(v == 'af3_fragments' for v in af3_src.values()))",
        "print('AF3 model(s) present but length differs from the protein (recorded, 3Di left empty):', len(af3_mismatch))")

nb.code("tdi, tdi_source, tdi_cov, tdi_nfrag = dict(tdi), dict(tdi_source), {g: 1.0 for g in tdi}, {g: 1 for g in tdi}",
        "for g, v in af3_acc.items(): tdi[g], tdi_source[g], tdi_cov[g], tdi_nfrag[g] = v, af3_src[g], af3_cov[g], af3_nfrag[g]",
        "n_before = len(tdi_source) - len(af3_acc)",
        "mm_db = list(mm); no_db_db = list(no_db)                       # state before AF3",
        "still = [g for g in nodes.index if g not in tdi]",
        "mm = [g for g in still if g in mismatch or g in af3_mismatch]   # a model exists but its length differs",
        "no_db = [g for g in still if g not in mm]                       # no usable model anywhere",
        "stitched = [g for g, v in af3_src.items() if v == 'af3_fragments']",
        "print('3Di before AF3:', n_before, '| after:', len(tdi), '| still missing:', len(still), '(length mismatch', len(mm), ', no model', len(no_db), ')')",
        "print('stitched from fragments:', len(stitched)); print(pd.Series([af3_cov[g] for g in stitched]).describe().round(3).to_string())",
        "print('fully covered stitched genes:', sum(af3_cov[g] == 1.0 for g in stitched), '; stitched with gaps:', sum(af3_cov[g] < 1.0 for g in stitched))",
        "print('overlapping residues re-assigned by the pLDDT rule:', af3_overlap_n)",
        "print('AF3 length-mismatch examples (gene: details, protein length):')",
        "for g in list(af3_mismatch)[:8]: print(' ', g, af3_mismatch[g], int(L_of[g]))")
nb.code("mismatch_rows = [dict(gene_id=g, protein_len=int(L_of[g]), detail=str(af3_mismatch.get(g, mismatch.get(g)))) for g in mm]",
        "pd.DataFrame(mismatch_rows).to_csv(OUT / 'tdi_length_mismatch.tsv', sep='\\t', index=False)",
        "pd.DataFrame([dict(accession=a, genes=','.join(v)) for a, v in amb_report.items()]).to_csv(OUT / 'af3_ambiguous_accessions.tsv', sep='\\t', index=False)",
        "print(len(mismatch_rows), 'length-mismatch genes written; ambiguous accessions:', len(amb_report))")

nb.md('## 4 · Cutting the regions',
      'All coordinates are 1-based inclusive in the DeepTMHMM file; below they become Python slices on the FULL sequence '
      'and the same slices are applied to the 3Di string.',
      '* **mature** = `[sp_end : L]` (signal peptide removed), else the whole sequence.',
      '* **N-terminal flank** = from the start of the mature sequence to the residue before the FIRST TM helix: '
      '`[sp_end : first_tm_start-1]`. **C-terminal flank** = from the residue after the LAST TM helix: '
      '`[last_tm_end : L]`. With no TM domain both are the whole mature sequence.',
      '* **Side** (`nterm_side`, `cterm_side`) is read from the per-residue topology string (I inside, O outside, S signal '
      'peptide, M membrane). A flank is `inside` or `outside` when at least 90% of its residues carry that one call, '
      'otherwise (mixed, empty, or no I/O call) `unknown`. Note that for a TM-less protein the side describes the whole '
      'protein as called by DeepTMHMM (globular proteins are mostly called inside or outside as a whole).')
nb.code('rows = []',
        'for g, seq in nodes.sequence.items():',
        '    r = dt.loc[g] if g in dt.index else None',
        '    L = len(seq)',
        '    sp_end = int(r.sp_end) if r is not None and pd.notna(r.sp_end) else 0',
        '    tms = sorted((int(r[f"tm_{i}_start"]), int(r[f"tm_{i}_end"])) for i in range(1, 25) if r is not None and pd.notna(r[f"tm_{i}_start"])) if r is not None else []',
        '    n_tm = len(tms)',
        '    if n_tm:',
        '        n0, n1 = sp_end, max(sp_end, tms[0][0] - 1)',
        '        c0, c1 = tms[-1][1], L',
        '    else:',
        '        n0, n1, c0, c1 = sp_end, L, sp_end, L',
        '    topo = r.topology_string if r is not None and isinstance(r.topology_string, str) and len(r.topology_string) == L else None',
        '    def side(a, b):',
        '        if topo is None or b <= a: return "unknown"',
        '        seg = topo[a:b]; i, o = seg.count("I"), seg.count("O")',
        '        if i + o == 0: return "unknown"',
        '        if i / len(seg) >= 0.9: return "inside"',
        '        if o / len(seg) >= 0.9: return "outside"',
        '        return "unknown"',
        '    t = tdi.get(g, "")',
        '    cut = lambda x, a, b: x[a:b] if x else ""',
        '    rows.append(dict(gene_id=g,',
        '        seq_full=seq, seq_mature=seq[sp_end:], seq_nterm_flank=seq[n0:n1], seq_cterm_flank=seq[c0:c1],',
        '        tdi_full=t, tdi_mature=cut(t, sp_end, L), tdi_nterm_flank=cut(t, n0, n1), tdi_cterm_flank=cut(t, c0, c1),',
        '        n_tm=n_tm, has_signal_peptide=sp_end > 0, nterm_side=side(n0, n1), cterm_side=side(c0, c1),',
        '        seq_len=L, mature_len=L - sp_end, nterm_len=n1 - n0, cterm_len=c1 - c0,',
        '        tdi_source=tdi_source.get(g, "none"), tdi_coverage=round(tdi_cov.get(g, 0.0), 4), tdi_n_fragments=tdi_nfrag.get(g, 0), flanks_identical=seq[n0:n1] == seq[c0:c1]))',
        "traits = pd.DataFrame(rows).set_index('gene_id')",
        'print(traits.shape)',
        'traits.drop(columns=[c for c in traits if c.startswith(("seq_full","seq_mat","seq_n","seq_c","tdi_full","tdi_mat","tdi_n","tdi_c"))]).describe(include="all").T')

nb.md('## 5 · Check 2: five named proteins with a signal peptide and a TM domain',
      'For eyeballing: length and first/last 20 residues of every region, sequence and 3Di. The first five proteins that '
      'carry both features and have 3Di are shown (with a short name where the product description is well known).')
nb.code("pick = traits[(traits.has_signal_peptide) & (traits.n_tm >= 1) & (traits.tdi_source != 'none') & (traits.n_tm <= 3)].index[:5]",
        'for g in pick:',
        '    t = traits.loc[g]',
        '    print(f"== {g}: n_tm={t.n_tm} seq_len={t.seq_len} mature_len={t.mature_len} nterm_len={t.nterm_len} cterm_len={t.cterm_len} nterm_side={t.nterm_side} cterm_side={t.cterm_side} 3Di={t.tdi_source}")',
        '    for col in ("seq_full", "seq_mature", "seq_nterm_flank", "seq_cterm_flank", "tdi_full", "tdi_mature", "tdi_nterm_flank", "tdi_cterm_flank"):',
        '        v = t[col]; print(f"  {col:16s} len={len(v):5d}  {v[:20]} ... {v[-20:]}")')

nb.md('## 6 · Checks 3 and 4',
      'Genes with no TM domain must have identical flanks equal to the mature sequence; every 3Di region must have the '
      'same length as its sequence region.')
nb.code('no_tm = traits[traits.n_tm == 0]',
        'assert (no_tm.seq_nterm_flank == no_tm.seq_cterm_flank).all() and (no_tm.seq_nterm_flank == no_tm.seq_mature).all()',
        'assert (no_tm.tdi_nterm_flank == no_tm.tdi_cterm_flank).all() and no_tm.flanks_identical.all()',
        "print('check 3 passed:', len(no_tm), 'genes with n_tm == 0 have identical flanks equal to the mature sequence')",
        'has = traits[traits.tdi_source != "none"]',
        'for a, b in (("seq_full", "tdi_full"), ("seq_mature", "tdi_mature"), ("seq_nterm_flank", "tdi_nterm_flank"), ("seq_cterm_flank", "tdi_cterm_flank")):',
        '    assert (has[a].str.len() == has[b].str.len()).all(), a',
        '    assert (traits.loc[traits.tdi_source == "none", b] == "").all(), b',
        "print('check 4 passed: every 3Di region has the length of its sequence region for', len(has), 'genes; empty for the other', len(traits) - len(has))",
        'assert (traits.nterm_len >= 0).all() and (traits.cterm_len >= 0).all()',
        'print("TM proteins whose flank is empty (TM at the very end/start):", int(((traits.n_tm > 0) & ((traits.nterm_len == 0) | (traits.cterm_len == 0))).sum()))')

nb.md('## 7 · Register the columns in the information map',
      'Every column is a tier-2 trait with `normalization="not_applicable"` and `role="measurement"`. The string columns are '
      'registered as kind `category` and the numeric/boolean ones as `measure`. `category` does not describe a long '
      'sequence string well (its rule is "one label per gene"), but it is the only kind that holds a single value per gene '
      'without implying a number or a yes-only flag, so it is used here. The families are held out together: '
      '`proteome_sequence`, `foldseek_3di` and `deeptmhmm`.')
nb.code('from parascape.information_map import Trait, register, TRAITS',
        "desc = {",
        " 'seq_full': 'Annotated protein sequence', 'seq_mature': 'Protein sequence with the DeepTMHMM signal peptide removed (no propeptide cleavage)',",
        " 'seq_nterm_flank': 'Mature sequence from its start to just before the first TM helix (whole mature sequence if no TM)',",
        " 'seq_cterm_flank': 'Sequence from just after the last TM helix to the C-terminus (whole mature sequence if no TM)',",
        " 'tdi_full': 'Foldseek 3Di string of the whole protein', 'tdi_mature': '3Di string of the mature region',",
        " 'tdi_nterm_flank': '3Di string of the N-terminal flank', 'tdi_cterm_flank': '3Di string of the C-terminal flank',",
        " 'n_tm': 'Number of DeepTMHMM transmembrane helices', 'has_signal_peptide': 'DeepTMHMM predicts a signal peptide',",
        " 'nterm_side': 'Side of the membrane of the N-terminal flank: inside, outside or unknown', 'cterm_side': 'Side of the membrane of the C-terminal flank: inside, outside or unknown',",
        " 'seq_len': 'Protein length', 'mature_len': 'Mature length', 'nterm_len': 'N-terminal flank length', 'cterm_len': 'C-terminal flank length',",
        " 'tdi_source': 'Where the 3Di string came from (naive_db, cleaned_db, af3, af3_fragments, none)', 'tdi_coverage': 'Fraction of residues with a real 3Di letter (0-1); gaps of stitched fragments are X', 'tdi_n_fragments': 'Number of AF3 fragment models stitched (1 for a single model, 0 when no 3Di)', 'flanks_identical': 'N- and C-terminal flanks are identical (always true without a TM domain)'}",
        "cats = {'seq_full','seq_mature','seq_nterm_flank','seq_cterm_flank','tdi_full','tdi_mature','tdi_nterm_flank','tdi_cterm_flank','nterm_side','cterm_side','tdi_source'}",
        'def fam(c): return "proteome_sequence" if c.startswith("seq_") and c not in ("seq_len",) else "foldseek_3di" if c.startswith("tdi_") and c != "tdi_source" else "deeptmhmm"',
        'for c in traits.columns:',
        '    register(Trait(c, "category" if c in cats else "measure", desc[c], tier=2, family=fam(c), normalization="not_applicable", role="measurement"))',
        "pd.DataFrame([(c, TRAITS[c].kind, TRAITS[c].family, TRAITS[c].tier) for c in traits.columns], columns=['trait', 'kind', 'family', 'tier'])")

nb.md('## 8 · Coverage and outputs')
nb.code("traits.to_parquet(OUT / 'sequence_traits.parquet')",
        'cov = []',
        'for c in traits.columns:',
        '    s = traits[c]',
        '    n = int((s.astype(str) != "").sum()) if not (pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s)) else int(s.notna().sum())',
        '    if c in ("nterm_side", "cterm_side"): n = int((s != "unknown").sum())',
        '    if c == "tdi_source": n = int((s != "none").sum())',
        '    cov.append(dict(item=f"column {c} non-empty", genes=n))',
        'cov += [dict(item="no TM domain (n_tm == 0)", genes=int((traits.n_tm == 0).sum())),',
        '        dict(item="with TM domain", genes=int((traits.n_tm > 0).sum())),',
        '        dict(item="signal peptide", genes=int(traits.has_signal_peptide.sum())),',
        '        dict(item="signal peptide and TM", genes=int((traits.has_signal_peptide & (traits.n_tm > 0)).sum())),',
        '        dict(item="3Di from naive_db", genes=int((traits.tdi_source == "naive_db").sum())),',
        '        dict(item="3Di from cleaned_db", genes=int((traits.tdi_source == "cleaned_db").sum())),',
        '        dict(item="3Di from af3 (whole model)", genes=int((traits.tdi_source == "af3").sum())),',
        '        dict(item="3Di from af3_fragments (stitched)", genes=int((traits.tdi_source == "af3_fragments").sum())),',
        '        dict(item="stitched with gaps (coverage < 1)", genes=int(((traits.tdi_source == "af3_fragments") & (traits.tdi_coverage < 1)).sum())),',
        '        dict(item="no 3Di: model(s) of different length (excluded)", genes=len(mm)),',
        '        dict(item="no 3Di: no usable model anywhere", genes=len(no_db)),',
        '        dict(item="flanks identical", genes=int(traits.flanks_identical.sum())),',
        '        dict(item="TM genes with nterm_side known", genes=int(((traits.n_tm > 0) & (traits.nterm_side != "unknown")).sum())),',
        '        dict(item="TM genes with cterm_side known", genes=int(((traits.n_tm > 0) & (traits.cterm_side != "unknown")).sum()))]',
        'cov = pd.DataFrame(cov); cov["fraction"] = (cov.genes / len(traits)).round(3)',
        "cov.to_csv(OUT / 'coverage.tsv', sep='\\t', index=False)",
        "pd.set_option('display.width', 200); print(cov.to_string())",
        "print(traits.nterm_side.value_counts().to_dict(), traits.cterm_side.value_counts().to_dict())")
nb.code("readme = f'''# Sequence and 3Di traits (instruction 72, 2026-11-04)\n\n"
        "One row per gene of parascape/data/nodes.parquet ({len(traits):,} genes), file sequence_traits.parquet, index gene_id. "
        "Built by scripts/build_sequence_traits_notebook.py (notebook notebooks/sequence_traits_2026_11_04.ipynb).\n\n"
        "## Rules\n"
        "* seq_mature: DeepTMHMM signal peptide removed (sequence[sp_end:]); otherwise equal to seq_full.\n"
        "* seq_nterm_flank: start of the mature sequence to the residue before the FIRST TM helix; seq_cterm_flank: residue after the LAST TM helix to the C-terminus.\n"
        "* NO-TM RULE: a protein with no TM domain has BOTH flanks equal to the whole mature sequence, so the two flank columns are identical (intended; asserted).\n"
        "* tdi_*: the same four regions with the same coordinates, cut from the Foldseek 3Di string, only where its length equals the protein length; empty otherwise.\n"
        "* nterm_side / cterm_side: inside or outside when >= 90% of the flank residues carry that I/O call in the DeepTMHMM per-residue topology string; unknown if mixed, empty or absent.\n"
        "* tdi_source: naive_db preferred, then cleaned_db (suffixes such as __sol stripped; the record matching the protein length is chosen), then af3, then af3_fragments, else none.\n"
        "* af3: Foldseek 3Di computed here (version {FS_VERSION}) from the whole-protein AlphaFold3 model af3_<UniProt>_model.cif in /media/carruthers/mnt1/af3/af_output, UniProt mapped to the gene through the accession column of deeptmhmm_toxo.csv; used only where no database 3Di of the right length exists.\n"
        "* af3_fragments: folders af3_<UniProt>_f<i>of<n>_<start>-<end> stitched at their 1-based inclusive ranges into a string of the full protein length; residues in no fragment are 'X' (not a 3Di letter); overlaps are resolved by the higher mean pLDDT (atom_plddts of the confidences json), else the earlier fragment. tdi_coverage = fraction of residues with a real letter; tdi_n_fragments = models used.\n"
        "* ACCEPTANCE: a 3Di string is kept only when its length equals the nodes.parquet protein length (stitched: every fragment length equals its range and the ranges fit inside the protein). Otherwise empty and listed in tdi_length_mismatch.tsv; nothing padded or trimmed.\n"
        "* Stitched strings contain 'X'; downstream code must treat X as missing.\n\n"
        "## Coordinate convention check\n{CONVENTION}. sequence[start-1:end] equals the recorded tm_<i>_sequence / sp_sequence.\n\n"
        "## Gaps\n"
        "* No propeptide / TEXEL cleavage: no reliable per-gene table; only the signal peptide is removed.\n"
        "* 3Di is missing for {int((traits.tdi_source == 'none').sum()):,} genes: {len(no_db):,} with no usable model anywhere (no database record and no af3_<UniProt> folder that maps unambiguously to the gene) and {len(mm):,} whose model(s) differ in length from the annotated protein (tdi_length_mismatch.tsv; excluded, not trimmed).\n"
        "* The AF3 tree has {len(entries):,} distinct af3_<UniProt> ids, all UniProt (none TGME49-named); the TGME49-named folders (bait_, pep_, xl_, ...) are domain or complex jobs and are not used. {len(unmapped)} accessions are not in the DeepTMHMM accession column and {len(ambiguous)} map to more than one gene (attached to none; af3_ambiguous_accessions.tsv).\n"
        "* AF3 model sequence was not compared with the annotated sequence (only the length), so an isoform of identical length would pass.\n"
        "* DeepTMHMM file has {N_DT_ROWS} rows (the 8,141 quoted in the spec counts the header line), all identifiers unique, all present in nodes.parquet.\n"
        "* Long strings are registered with kind 'category' (no better-fitting kind).\n\n"
        "## Coverage\n\n```\n{cov.to_string(index=False)}\n```\n'''",
        "(OUT / 'README.md').write_text(readme)", "print(readme[:400])")

nb.write(str(ROOT / 'notebooks/sequence_traits_2026_11_04.ipynb'))
