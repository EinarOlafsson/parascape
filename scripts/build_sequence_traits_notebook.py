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
        "lacking = set(nodes.index) - set(tdi)      # genes with no database 3Di of the right length (before AF3)",
        "# folders of accessions that are ambiguous (>1 gene) or absent from the DeepTMHMM table are converted too: section 3c",
        "# decides by sequence alignment which gene (if any) they belong to",
        "amb_unmapped = set(ambiguous) | set(unmapped)",
        "todo = sorted({e['folder'] for g in need for a in gene_acc[g] for e in entries[a]} | {e['folder'] for a in amb_unmapped for e in entries[a]})",
        "print('model.cif files to convert:', len(todo), '(of which for ambiguous / unmapped accessions:', len({e['folder'] for a in amb_unmapped for e in entries[a]}), ')')",
        'af3_3di, af3_seq = {}, {}',
        "def read_ss(prefix, suffix='_ss'):          # suffix '' reads the SEQUENCE database that sits beside <prefix>_ss",
        '    names = {}',
        "    for ln in open(f'{prefix}.lookup'):",
        "        f = ln.rstrip('\\n').split('\\t')",
        '        if len(f) >= 2: names[int(f[0])] = f[1]',
        "    recs = open(f'{prefix}{suffix}', 'rb').read().split(b'\\x00')",
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
        "    seqs = read_ss(str(d / 'db'), '')",
        "    for name, ss in read_ss(str(d / 'db')).items():",
        "        name = re.sub(r'\\.cif$', '', name)",
        "        key = name if name in batch else re.sub(r'_[A-Za-z0-9]$', '', name)    # strip a chain suffix if foldseek added one",
        '        if key in batch:',
        '            af3_3di.setdefault(key, []).extend(ss); af3_seq.setdefault(key, []).extend(seqs[name])',
        "            assert [len(x) for x in ss] == [len(x) for x in seqs[name]], name      # the sequence database is in the same order as the 3Di one",
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
nb.md('## 3c · Filling the rest by transferring 3Di through a sequence alignment',
      '**Why the sequences differ.** The protein in `nodes.parquet` is the ToxoDB gene model of ME49. The structure that '
      'AlphaFold predicted was folded from a UniProt entry, and the UniProt entry is often NOT the same protein: ToxoDB may '
      'have merged two neighbouring predicted genes into one long gene (for example `TGME49_500303`, 6,319 residues, is '
      'covered by two UniProt entries of 2,083 and 1,861 residues), extended or shortened a start or a stop, or chosen a '
      'different exon structure. Downloading the structure again would fetch the same mismatched model, and a 3Di string of '
      'the wrong length cannot be laid on our protein residue by residue.',
      '**What transferring 3Di through an alignment means.** A 3Di string has one letter per residue of the structure, '
      'describing the local 3D neighbourhood of that residue. If most of the structure\'s own sequence is identical to a '
      'stretch of our protein, then residue *k* of the structure is residue *j* of our protein, and its 3Di letter belongs to '
      'position *j*. We align the two sequences (global alignment, but a free overhang at both ends, because the models '
      'usually start or stop at a different place), copy the 3Di letter of every aligned residue onto our residue, and write '
      '`X` (missing) wherever no structure residue lines up. With several structures for one gene (two halves of a merged '
      'gene) all of them are transferred; where they overlap, the one with the higher identity keeps the position.',
      '**Why at least 95 % identity.** 3Di describes a specific fold. A different isoform or a paralogue may fold differently '
      'even when 80 % of its residues are shared, and a 3Di letter copied over from it would then be wrong. At 95 % the '
      'structure and our protein are, for practical purposes, the same polypeptide with a few annotation differences. The '
      'gene must also be at least half covered (a short matching stretch tells us little about the whole protein). The '
      'accepted rule is 95 % identity (matches over alignment columns, internal gaps counted as mismatches) per structure and '
      '50 % of our protein covered by the merged result. The 90 % / 30 % trade-off is printed for comparison but not used.',
      '**An aligned 3Di is weaker evidence than a direct one.** Direct 3Di (`naive_db`, `cleaned_db`, `af3`) is computed '
      'from a structure whose sequence is exactly ours. An aligned one is inferred from a similar sequence and may contain '
      '`X` gaps. Any model that uses the `tdi_*` columns must carry `tdi_source` (and `tdi_align_identity`, '
      '`tdi_align_coverage`) with it, so that aligned genes can be down-weighted, flagged or held out. Aligned strings are '
      'full length, like the others; `X` means missing and must never be read as a 3Di letter.')
nb.code(r'''import numpy as np
from numba import njit

# Semi-global (free end gaps on both sequences) affine-gap alignment, written here because Biopython is not installed.
# Score: identity scoring (match +5, mismatch -4; gap of length k costs 12 + 2(k-1)); no substitution matrix is available offline.
@njit
def _sg(a, b, match, mismatch, gopen, gext):
    n, m = len(a), len(b)
    NEG = -10**8
    tb = np.zeros((n + 1, m + 1), dtype=np.uint8)
    Hp = np.zeros(m + 1, dtype=np.int32); Fp = np.full(m + 1, NEG, dtype=np.int32)
    Hc = np.zeros(m + 1, dtype=np.int32); Fc = np.full(m + 1, NEG, dtype=np.int32)
    best = NEG; bi = 0; bj = 0
    for i in range(1, n + 1):
        Hc[0] = 0; Fc[0] = NEG; E = NEG
        for j in range(1, m + 1):
            e_open = Hc[j - 1] + gopen; e_ext = E + gext
            if e_ext >= e_open:
                E = e_ext; bits = 4
            else:
                E = e_open; bits = 0
            f_open = Hp[j] + gopen; f_ext = Fp[j] + gext
            if f_ext >= f_open:
                Fc[j] = f_ext; bits |= 8
            else:
                Fc[j] = f_open
            s = match if a[i - 1] == b[j - 1] else mismatch
            h = Hp[j - 1] + s; src = 0
            if E > h:
                h = E; src = 1
            if Fc[j] > h:
                h = Fc[j]; src = 2
            Hc[j] = h
            tb[i, j] = bits | src
        if Hc[m] > best:
            best = Hc[m]; bi = i; bj = m
        for j in range(m + 1):
            Hp[j] = Hc[j]; Fp[j] = Fc[j]
    for j in range(1, m + 1):
        if Hp[j] > best:
            best = Hp[j]; bi = n; bj = j
    return tb, best, bi, bj

@njit
def _trace(tb, i, j):
    pi = np.empty(i + j + 2, dtype=np.int32); pj = np.empty(i + j + 2, dtype=np.int32)
    k = 0; state = 0
    while i > 0 and j > 0:
        t = tb[i, j]
        if state == 0:
            state = t & 3
            if state == 0:
                pi[k] = i - 1; pj[k] = j - 1; k += 1; i -= 1; j -= 1
        elif state == 1:
            j -= 1
            if (t & 4) == 0: state = 0
        else:
            i -= 1
            if (t & 8) == 0: state = 0
    return pi[:k], pj[:k]

def semiglobal(cand, ours, match=5, mismatch=-4, gopen=-12, gext=-2):
    """Align cand (the 3Di donor) to ours. Returns 0-based aligned pairs, number of identical pairs, alignment columns
    (first aligned pair to last, internal gaps of either sequence included)."""
    a = np.frombuffer(cand.encode('ascii'), dtype=np.uint8); b = np.frombuffer(ours.encode('ascii'), dtype=np.uint8)
    tb, best, bi, bj = _sg(a, b, match, mismatch, gopen, gext)
    pi, pj = _trace(tb, bi, bj); pi, pj = pi[::-1], pj[::-1]
    if len(pi) == 0: return pi, pj, 0, 0
    ident = int((a[pi] == b[pj]).sum())
    cols = int((pi[-1] - pi[0] + 1) + (pj[-1] - pj[0] + 1) - len(pi))
    return pi, pj, ident, cols

# unit tests on a real protein: a fragment, a fragment with a deletion, and the identity of a protein with itself
s0 = nodes.sequence['TGME49_500303']
for name, c in (('fragment 1000-2500', s0[1000:2500]), ('two pieces with a 20-residue deletion', s0[100:300] + s0[320:900]), ('self 0-2000', s0[:2000])):
    pi, pj, idn, cols = semiglobal(c, s0)
    print(f'{name}: aligned {len(pi)}, identical {idn}, columns {cols}, first/last target position {pj[0]}-{pj[-1]}')
assert semiglobal(s0[1000:2500], s0)[1][0] == 1000 and semiglobal(s0[1000:2500], s0)[2] == 1500''')

nb.code(r'''# ---- candidates: every structure sequence we can pair with a gene (database records, AF3 models), with its own 3Di
def read_db_full(prefix):
    """name -> list of (sequence, 3Di, record name). The sequence database <prefix> sits beside <prefix>_ss; same order."""
    names = {}
    for ln in open(f'{prefix}.lookup'):
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 2: names[int(f[0])] = f[1]
    ss = open(f'{prefix}_ss', 'rb').read().split(b'\x00'); sq = open(prefix, 'rb').read().split(b'\x00')
    out = {}; n_len_ok = n_tot = 0
    for i, nm in names.items():
        a, b = ss[i].decode('ascii', 'ignore').strip(), sq[i].decode('ascii', 'ignore').strip()
        n_tot += 1; n_len_ok += (len(a) == len(b))
        if not a or len(a) != len(b): continue
        raw = nm
        nm = nm[4:] if nm.startswith('AF3-') else nm
        nm = nm[:-4] if nm.endswith('.cif') else nm
        nm = nm.split('__')[0]
        if nm.startswith('TGME49_'): out.setdefault(nm, []).append((b, a, raw))
    print(prefix.split('/')[-1], ': stored sequence has the length of its 3Di string for', n_len_ok, 'of', n_tot, 'records')
    return out
dbfull = {k: read_db_full(f'/home/carruthers/foldseek/{k}') for k in ('naive_db', 'cleaned_db')}
# verify the sequence database really is in the same order as the 3Di one: wherever the stored sequence is as long as ours it should BE ours
chk = [(g, r[0] == nodes.sequence[g]) for g, recs in dbfull['naive_db'].items() if g in nodes.index for r in recs if len(r[0]) == len(nodes.sequence[g])]
print('naive_db records of the same length as the annotated protein:', len(chk), '; with identical sequence:', sum(ok for _, ok in chk))
assert sum(ok for _, ok in chk) / len(chk) > 0.99, 'sequence database order does not match the 3Di order'

# AF3 candidates: every chain of every converted model
af3_cands = {}      # folder -> list of (sequence, 3Di, folder#chain)
for fo in af3_3di:
    for k, (ss_, sq_) in enumerate(zip(af3_3di[fo], af3_seq[fo])):
        af3_cands.setdefault(fo, []).append((sq_, ss_, f'{fo}#{k + 1}'))
print('AF3 models with sequence and 3Di:', len(af3_cands))
assert all(len(a) == len(b) for v in af3_cands.values() for a, b, _ in v)

lacking = [g for g in nodes.index if g not in tdi]       # genes still lacking 3Di after sections 3 and 3b
print('genes still lacking 3Di before the alignment:', len(lacking))

# gene -> candidate list: records of the gene in the databases, AF3 folders of its (unambiguously or ambiguously) mapped accessions
pool = {g: [] for g in lacking}
for g in lacking:
    for src in ('naive_db', 'cleaned_db'):
        for sq_, ss_, raw in dbfull[src].get(g, []): pool[g].append((sq_, ss_, f'{src}:{raw}', src))
    accs = set(gene_acc.get(g, [])) | {a for a in ambiguous if g in a2g[a]}
    for a in sorted(accs):
        for e in entries[a]:
            for sq_, ss_, cid in af3_cands.get(e['folder'], []): pool[g].append((sq_, ss_, cid, 'af3'))

# accessions the DeepTMHMM table does not know: a 6-mer screen finds which lacking genes (if any) they resemble
K = 6
kidx = collections.defaultdict(set)
for g in lacking:
    s_ = nodes.sequence[g]
    for i in range(len(s_) - K + 1): kidx[s_[i:i + K]].add(g)
n_screen = n_screen_hit = 0
for a in sorted(unmapped):
    for e in entries[a]:
        for sq_, ss_, cid in af3_cands.get(e['folder'], []):
            n_screen += 1
            ks = {sq_[i:i + K] for i in range(len(sq_) - K + 1)}
            hits = collections.Counter(g for k_ in ks for g in kidx.get(k_, ()))
            sel = [g for g, c in hits.items() if c >= max(10, 0.5 * len(ks))]
            n_screen_hit += bool(sel)
            for g in sel: pool[g].append((sq_, ss_, cid, 'af3_unmapped'))
print('structures of accessions unknown to DeepTMHMM:', n_screen, '; sharing >= 50% of their 6-mers with a gene that lacks 3Di:', n_screen_hit)
# the same sequence from several sources is aligned once (first kept)
for g in pool:
    seen, keep = set(), []
    for c in pool[g]:
        if c[0] not in seen: seen.add(c[0]); keep.append(c)
    pool[g] = keep
print('genes with at least one candidate:', sum(bool(v) for v in pool.values()), '; candidates in total:', sum(len(v) for v in pool.values()))''')

nb.code(r'''# ---- align every candidate to its gene
import time
t0 = time.time(); arows = []; aln = {}
for g in lacking:
    ours = nodes.sequence[g]; L = len(ours)
    for sq_, ss_, cid, kind in pool[g]:
        pi, pj, idn, cols = semiglobal(sq_, ours)
        aln[(g, cid)] = (pi, pj, ss_)
        arows.append(dict(gene_id=g, candidate=cid, kind=kind, protein_len=L, candidate_len=len(sq_), aligned_residues=len(pi),
                          identical=idn, aln_columns=cols, identity=idn / cols if cols else 0.0, identity_ungapped=idn / len(pi) if len(pi) else 0.0,
                          coverage=len(pi) / L))
A = pd.DataFrame(arows)
print(f'{len(A)} alignments in {time.time() - t0:.0f} s; genes aligned:', A.gene_id.nunique())
A.round(4).to_csv(OUT / 'tdi_alignment_candidates.tsv', sep='\t', index=False)
A.groupby('gene_id').identity.max().describe().round(3)''')

nb.code(r'''# ---- merge the candidates of each gene and apply the acceptance rule
def merge(g, min_id):
    """Fill a full-length 3Di string from every candidate with identity >= min_id; higher identity first, no overwriting."""
    L = len(nodes.sequence[g]); arr = ['X'] * L; used = []
    sub = A[(A.gene_id == g) & (A.identity >= min_id)].sort_values(['identity', 'aligned_residues'], ascending=False)
    ident_sum = cols_sum = 0
    for r in sub.itertuples():
        pi, pj, ss_ = aln[(g, r.candidate)]
        new = 0
        for i, j in zip(pi, pj):
            if arr[j] == 'X': arr[j] = ss_[i]; new += 1
        if new:
            used.append(r.candidate); ident_sum += r.identical; cols_sum += r.aln_columns
    s = ''.join(arr)
    return s, (1 - s.count('X') / L), used, (ident_sum / cols_sum if cols_sum else 0.0)
res = {}
for rule, (mi, mc) in {'95/50': (0.95, 0.50), '90/30': (0.90, 0.30)}.items():
    ok = {}
    for g in A.gene_id.unique():
        s, cov, used, idn = merge(g, mi)
        if used and cov >= mc: ok[g] = (s, cov, used, idn)
    res[rule] = ok
print('genes accepted at identity >= 95% and coverage >= 50% (the rule used):', len(res['95/50']))
print('genes that would be accepted at 90% / 30% (comparison only):', len(res['90/30']), '; extra over 95/50:', len(set(res['90/30']) - set(res['95/50'])))
best = A.groupby('gene_id').agg(best_identity=('identity', 'max'), best_candidate_coverage=('coverage', 'max'))
print('best single-candidate identity per gene (', len(best), 'genes ): quantiles')
print(best.best_identity.quantile([0, .1, .25, .5, .75, .9, 1]).round(3).to_string())
acc = res['95/50']
aligned_ok = set(acc)
al_tab = pd.DataFrame([dict(gene_id=g, identity=round(v[3], 4), coverage=round(v[1], 4), n_candidates=len(v[2]), candidates=','.join(v[2]))
                       for g, v in acc.items()]).set_index('gene_id')
al_tab['alignment_length'] = [int(A[(A.gene_id == g) & A.candidate.isin(acc[g][2])].aln_columns.sum()) for g in al_tab.index]
al_tab.to_csv(OUT / 'tdi_aligned_genes.tsv', sep='\t')
print('accepted genes: identity and coverage'); print(al_tab[['identity', 'coverage', 'n_candidates', 'alignment_length']].describe().round(3).to_string())
print('candidates used per gene:', al_tab.n_candidates.value_counts().sort_index().to_dict())
print('from sources:', pd.Series([A[(A.gene_id == g) & (A.candidate == v[2][0])].kind.iloc[0] for g, v in acc.items()]).value_counts().to_dict())
cmp = pd.DataFrame({'rule': ['identity >= 95%, coverage >= 50% (accepted)', 'identity >= 90%, coverage >= 30% (comparison)'],
                    'genes': [len(res['95/50']), len(res['90/30'])]})
cmp''')

nb.code(r'''# ---- apply: aligned strings join the 3Di dictionaries
tdi_align_id, tdi_align_cov = {}, {}
for g, (s, cov, used, idn) in acc.items():
    assert g not in tdi and len(s) == len(nodes.sequence[g])
    tdi[g] = s; tdi_source[g] = 'aligned' if len(used) == 1 else 'aligned_multi'
    tdi_cov[g], tdi_nfrag[g], tdi_align_id[g], tdi_align_cov[g] = cov, len(used), idn, cov
print('3Di strings now:', len(tdi), 'of', len(nodes), '; by source:', pd.Series(tdi_source).value_counts().to_dict())
# genes still empty and why
mm_before = list(mm)
mm = [g for g in mm if g not in aligned_ok]
no_db = [g for g in nodes.index if g not in tdi and g not in mm]
rej = [g for g in mm if g in A.gene_id.values]
print('still without 3Di:', len(mm) + len(no_db), '=', len(mm), 'with a model whose alignment fails the 95/50 rule +', len(no_db), 'with no usable model')
print('of the', len(mm), ': candidates aligned but best identity < 95%:', sum(1 for g in mm if g in best.index and best.best_identity[g] < 0.95),
      '; identity >= 95% but merged coverage < 50%:', sum(1 for g in mm if g in best.index and best.best_identity[g] >= 0.95))''')


nb.code("mismatch_rows = [dict(gene_id=g, protein_len=int(L_of[g]), detail=str(af3_mismatch.get(g, mismatch.get(g))), best_alignment_identity=round(float(best.best_identity.get(g, np.nan)), 4), best_single_candidate_coverage=round(float(best.best_candidate_coverage.get(g, np.nan)), 4)) for g in mm]",
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
        '        tdi_source=tdi_source.get(g, "none"), tdi_coverage=round(tdi_cov.get(g, 0.0), 4), tdi_n_fragments=tdi_nfrag.get(g, 0), tdi_align_identity=tdi_align_id.get(g, np.nan), tdi_align_coverage=tdi_align_cov.get(g, np.nan), flanks_identical=seq[n0:n1] == seq[c0:c1]))',
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
        " 'tdi_source': 'Where the 3Di string came from (naive_db, cleaned_db, af3, af3_fragments, aligned, aligned_multi, none); aligned = transferred through a sequence alignment, weaker evidence', 'tdi_coverage': 'Fraction of residues with a real 3Di letter (0-1); gaps of stitched fragments are X', 'tdi_n_fragments': 'Number of AF3 fragment models stitched, or of aligned structures merged (1 for a single model, 0 when no 3Di)', 'tdi_align_identity': 'Identity (matches over alignment columns) of the structure sequence to the annotated protein for an aligned 3Di; NaN otherwise', 'tdi_align_coverage': 'Fraction of the protein covered by aligned structures for an aligned 3Di; NaN otherwise', 'flanks_identical': 'N- and C-terminal flanks are identical (always true without a TM domain)'}",
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
        '        dict(item="3Di from aligned (one structure)", genes=int((traits.tdi_source == "aligned").sum())),',
        '        dict(item="3Di from aligned_multi (several structures)", genes=int((traits.tdi_source == "aligned_multi").sum())),',
        '        dict(item="aligned with gaps (coverage < 1)", genes=int((traits.tdi_source.isin(["aligned", "aligned_multi"]) & (traits.tdi_coverage < 1)).sum())),',
        '        dict(item="would also pass at 90% identity / 30% coverage (not used)", genes=len(res["90/30"]) - len(res["95/50"])),',
        '        dict(item="stitched with gaps (coverage < 1)", genes=int(((traits.tdi_source == "af3_fragments") & (traits.tdi_coverage < 1)).sum())),',
        '        dict(item="no 3Di: model(s) of different length, alignment fails 95/50", genes=len(mm)),',
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
        "* tdi_*: the same four regions with the same coordinates, cut from the Foldseek 3Di string (full protein length); empty only for genes with no 3Di.\n"
        "* nterm_side / cterm_side: inside or outside when >= 90% of the flank residues carry that I/O call in the DeepTMHMM per-residue topology string; unknown if mixed, empty or absent.\n"
        "* tdi_source: naive_db preferred, then cleaned_db (suffixes such as __sol stripped; the record matching the protein length is chosen), then af3, then af3_fragments, then aligned / aligned_multi (section below), else none.\n"
        "* aligned / aligned_multi: for genes whose structure(s) have a DIFFERENT sequence from the ToxoDB gene model (ToxoDB model merged or extended relative to the UniProt entry AlphaFold folded). The structure own sequence is aligned to ours (semi-global, free end gaps, identity scoring match +5 / mismatch -4 / gap 12+2(k-1); Biopython not installed, own numba implementation) and the 3Di letter of each aligned residue is copied to the residue it aligns to; positions without an aligned residue are X. Accepted only when a structure has >= 95% identity (matches over alignment columns, internal gaps counted) and the merged result covers >= 50% of the protein; several structures are merged, higher identity first, no overwriting. Candidates: naive_db and cleaned_db records of the gene, AF3 models of its accessions, AF3 models of accessions ambiguous between genes or unknown to DeepTMHMM (screened by shared 6-mers, then aligned). tdi_align_identity / tdi_align_coverage hold the identity (over the structures that contributed) and coverage; NaN for other sources; tdi_n_fragments = number of structures merged. Per-candidate numbers: tdi_alignment_candidates.tsv; accepted genes: tdi_aligned_genes.tsv. AN ALIGNED 3Di IS WEAKER EVIDENCE THAN A DIRECT ONE: any model using tdi_* must carry tdi_source.\n"
        "* af3: Foldseek 3Di computed here (version {FS_VERSION}) from the whole-protein AlphaFold3 model af3_<UniProt>_model.cif in /media/carruthers/mnt1/af3/af_output, UniProt mapped to the gene through the accession column of deeptmhmm_toxo.csv; used only where no database 3Di of the right length exists.\n"
        "* af3_fragments: folders af3_<UniProt>_f<i>of<n>_<start>-<end> stitched at their 1-based inclusive ranges into a string of the full protein length; residues in no fragment are 'X' (not a 3Di letter); overlaps are resolved by the higher mean pLDDT (atom_plddts of the confidences json), else the earlier fragment. tdi_coverage = fraction of residues with a real letter; tdi_n_fragments = models used.\n"
        "* ACCEPTANCE (direct sources): a 3Di string is kept only when its length equals the nodes.parquet protein length (stitched: every fragment length equals its range and the ranges fit inside the protein). Otherwise empty and listed in tdi_length_mismatch.tsv; nothing padded or trimmed.\n"
        "* Stitched strings contain 'X'; downstream code must treat X as missing.\n\n"
        "## Coordinate convention check\n{CONVENTION}. sequence[start-1:end] equals the recorded tm_<i>_sequence / sp_sequence.\n\n"
        "## Gaps\n"
        "* No propeptide / TEXEL cleavage: no reliable per-gene table; only the signal peptide is removed.\n"
        "* 3Di is missing for {int((traits.tdi_source == 'none').sum()):,} genes: {len(no_db):,} with no usable model anywhere (no database record and no af3_<UniProt> folder that maps unambiguously to the gene) and {len(mm):,} whose model(s) differ in length from the annotated protein and could not be aligned at >= 95% identity with >= 50% coverage (tdi_length_mismatch.tsv, with the best alignment identity; excluded, not trimmed). {len(aligned_ok):,} genes of the former mismatch/no-model group were filled by alignment (tdi_source aligned / aligned_multi); at 90%/30% {len(res['90/30']) - len(res['95/50']):,} more would pass.\n"
        "* The AF3 tree has {len(entries):,} distinct af3_<UniProt> ids, all UniProt (none TGME49-named); the TGME49-named folders (bait_, pep_, xl_, ...) are domain or complex jobs and are not used. {len(unmapped)} accessions are not in the DeepTMHMM accession column and {len(ambiguous)} map to more than one gene (attached to none; af3_ambiguous_accessions.tsv).\n"
        "* Whole-length direct sources (naive_db, cleaned_db, af3) are accepted on length only: the sequence was not compared with the annotated one, so an isoform of identical length would pass. Aligned sources are checked by identity.\n"
        "* DeepTMHMM file has {N_DT_ROWS} rows (the 8,141 quoted in the spec counts the header line), all identifiers unique, all present in nodes.parquet.\n"
        "* Long strings are registered with kind 'category' (no better-fitting kind).\n\n"
        "## Coverage\n\n```\n{cov.to_string(index=False)}\n```\n'''",
        "(OUT / 'README.md').write_text(readme)", "print(readme[:400])")

nb.write(str(ROOT / 'notebooks/sequence_traits_2026_11_04.ipynb'))
