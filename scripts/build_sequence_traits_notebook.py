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
        '        tdi_source=tdi_source.get(g, "none"), flanks_identical=seq[n0:n1] == seq[c0:c1]))',
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
        " 'tdi_source': 'Which Foldseek database gave the 3Di string (naive_db, cleaned_db, none)', 'flanks_identical': 'N- and C-terminal flanks are identical (always true without a TM domain)'}",
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
        '        dict(item="no 3Di: record of different length (excluded)", genes=len(mm)),',
        '        dict(item="no 3Di: not in either database", genes=len(no_db)),',
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
        "* tdi_source: naive_db preferred, then cleaned_db (suffixes such as __sol stripped; the record matching the protein length is chosen, since cleaned_db also holds domain fragments).\n\n"
        "## Coordinate convention check\n{CONVENTION}. sequence[start-1:end] equals the recorded tm_<i>_sequence / sp_sequence.\n\n"
        "## Gaps\n"
        "* No propeptide / TEXEL cleavage: no reliable per-gene table; only the signal peptide is removed.\n"
        "* 3Di is missing for {int((traits.tdi_source == 'none').sum()):,} genes: {len(no_db):,} in neither Foldseek database and {len(mm):,} whose AlphaFold-derived 3Di length differs from the annotated protein (excluded, not trimmed).\n"
        "* DeepTMHMM file has {N_DT_ROWS} rows (the 8,141 quoted in the spec counts the header line), all identifiers unique, all present in nodes.parquet.\n"
        "* Long strings are registered with kind 'category' (no better-fitting kind).\n\n"
        "## Coverage\n\n```\n{cov.to_string(index=False)}\n```\n'''",
        "(OUT / 'README.md').write_text(readme)", "print(readme[:400])")

nb.write(str(ROOT / 'notebooks/sequence_traits_2026_11_04.ipynb'))
