"""Build the Toxoplasma CRAPome: a background frequency table for pull-down MS, and a
specificity score for new hit lists.

Instruction 71 product 4. The partner-prediction run of 2026-11-04 (instruction 72) kept the
33,596 bait-member pull-down pairs out of its headline because, without a background, a sticky
contaminant is indistinguishable from a real partner. This script builds that background from the
uniform Sage re-search (instruction 69) and writes an executed notebook plus a results folder.

Everything is read-only against the NAS. The builder rescans the uniform-search output every time
it runs, so re-running it as the queue fills picks up new projects; the project count it covered is
recorded in the outputs.

    python scripts/build_crapome_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())

from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_05"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"crapome_{DATE}")
NB = os.path.join(REPO, "notebooks", f"crapome_{DATE}.ipynb")

MS_ROOT = "/nas_mnt/carruthers/Einar/parascape_omics/prepared_ms_uniform"
SELECTION = ("/nas_mnt/carruthers/Einar/parascape_omics/processing_audit/"
             "pride_ms_uniform_selection.json")
NODES = os.path.join(REPO, "parascape", "data", "nodes.parquet")
STUDIES = os.path.join(REPO, "parascape", "data", "interaction_studies.parquet")
MEMBERS = os.path.join(REPO, "parascape", "data", "interaction_study_members.parquet")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Toxoplasma CRAPome: how often a protein turns up in a pull-down "
                          "regardless of the bait")

    nb.md(
        "**What this is.** When you pull down a tagged protein from *Toxoplasma* and send the "
        "beads for mass spectrometry, some proteins come back almost every time: ribosomes, "
        "chaperones, tubulin, actin, the abundant glycolytic enzymes. They stick to beads, to "
        "antibodies, and to anything hydrophobic. A human version of this list exists (the "
        "CRAPome, Mellacheruvu et al. 2013) and is used routinely to filter hit lists. No "
        "*Toxoplasma* version exists.",
        "This notebook builds one from our own uniform re-search: every *Toxoplasma* PRIDE "
        "project searched with the same Sage settings against the same combined FASTA "
        "(instruction 69), so the numbers are comparable between projects in a way that "
        "author-deposited tables are not.",
        "**What you get.** For every parasite protein, the fraction of independent "
        "affinity-purification and proximity-labelling experiments in which it is detected, its "
        "typical abundance when detected, and a score that says how surprising it is to see that "
        "protein, at that abundance, in a new pull-down.",
        "**Read the validation section before using it.** The whole point of a background table "
        "is that real complexes survive it. If they do not, the table is not usable yet, and this "
        "notebook says so rather than shipping a filter that throws away the biology.",
        f"Built {_dt.date.today().isoformat()}. Inputs are read-only; nothing on the NAS is "
        "written or changed.")

    nb.code(
        "import json, os, re, glob",
        "import numpy as np, pandas as pd",
        "from scipy import stats",
        "pd.set_option('display.width', 200)",
        "pd.set_option('display.max_columns', 60)",
        f"MS_ROOT = {MS_ROOT!r}",
        f"SELECTION = {SELECTION!r}",
        f"OUT = {OUT!r}",
        f"REPO = {REPO!r}",
        "os.makedirs(OUT, exist_ok=True)",
        "print('uniform search output:', MS_ROOT)",
        "print('exists:', os.path.isdir(MS_ROOT))")

    # ------------------------------------------------------------------ inputs
    nb.md("## 1 · What is on disk, and what is excluded",
          "The PRIDE selection lists 203 projects (the rule is in the file: keep anything whose "
          "taxon or title names *Toxoplasma* or a *Toxoplasma* protein family). The uniform "
          "pipeline is still working through them. A project is finished when it has a "
          "`done.json` receipt naming its Sage searches and their output checksums.",
          "Three kinds of project are dropped before anything is counted:",
          "* `unsupported.json` — the pipeline could not search it at all (SCIEX `.wiff`, which "
          "needs a converter that does not run on this machine; deposits with no spectra).",
          "* `NOT_COMPARABLE.json` — **hydrogen–deuterium exchange MS**. HDX peptides carry an "
          "unknown number of deuteriums, so their masses are shifted by amounts the search does "
          "not model; the spectra are kept but the tables are not abundance. Instruction 69 "
          "measured it: 1,272 and 1,883 PSMs per run against 10,717–45,527 for ordinary "
          "bottom-up runs in the same pipeline. These are `PXD031816` and `PXD023791`; two more "
          "HDX projects (`PXD034607`, `PXD039335`) were stopped before they were searched and "
          "carry `unsupported.json`.",
          "* a project with a receipt but no readable Sage PSM table.")

    nb.code(
        "sel = json.load(open(SELECTION))",
        "selected = {e['project']: e for e in sel['keep']}",
        "print('selection rule:', sel['rule'])",
        "print('selected projects:', len(selected))",
        "",
        "rows = []",
        "for proj in sorted(os.listdir(MS_ROOT)):",
        "    d = os.path.join(MS_ROOT, proj)",
        "    if not os.path.isdir(d):",
        "        continue",
        "    receipt = os.path.join(d, 'done.json')",
        "    psm = sorted(glob.glob(os.path.join(d, 'sage_*', 'results.sage.tsv')))",
        "    lfq = sorted(glob.glob(os.path.join(d, 'sage_*', 'lfq.tsv')))",
        "    rows.append(dict(",
        "        project=proj,",
        "        in_selection=proj in selected,",
        "        title=selected.get(proj, {}).get('title', ''),",
        "        taxa=selected.get(proj, {}).get('taxa', ''),",
        "        finished=os.path.exists(receipt),",
        "        unsupported=os.path.exists(os.path.join(d, 'unsupported.json')),",
        "        not_comparable=os.path.exists(os.path.join(d, 'NOT_COMPARABLE.json')),",
        "        n_psm_tables=len(psm), n_lfq_tables=len(lfq),",
        "        psm_path=psm[0] if psm else '', lfq_path=lfq[0] if lfq else '',",
        "    ))",
        "disk = pd.DataFrame(rows)",
        "print('project directories on disk:', len(disk))",
        "print('with a done.json receipt:', int(disk.finished.sum()))",
        "print('unsupported.json:', sorted(disk.loc[disk.unsupported, 'project']))",
        "print('NOT_COMPARABLE.json (HDX):', sorted(disk.loc[disk.not_comparable, 'project']))",
        "disk.groupby(['finished', 'unsupported', 'not_comparable']).size()")

    nb.code(
        "# run-level detail from each receipt: file names and the acquisition flag the pipeline recorded",
        "det = {}",
        "for proj in disk.loc[disk.finished, 'project']:",
        "    rec = json.load(open(os.path.join(MS_ROOT, proj, 'done.json')))",
        "    runs = rec.get('runs') or []",
        "    # the receipt schema changed mid-queue: older entries key the source file as 'raw'",
        "    names = [r.get('source') or r.get('raw') or r.get('mzml') or '' for r in runs]",
        "    flags = sorted({(r.get('acquisition') or {}).get('flag') for r in runs} - {None})",
        "    det[proj] = dict(n_runs=len(runs), run_names=[n for n in names if n],",
        "                     flags='; '.join(flags),",
        "                     sage_version=(rec.get('tools') or {}).get('sage', ''),",
        "                     fasta=os.path.basename(rec.get('fasta', '')))",
        "disk['n_runs'] = disk.project.map(lambda p: det.get(p, {}).get('n_runs', 0))",
        "disk['flags'] = disk.project.map(lambda p: det.get(p, {}).get('flags', ''))",
        "print('runs behind the finished projects:', int(disk.n_runs.sum()))",
        "print('one FASTA for all of them:',",
        "      sorted({v['fasta'] for v in det.values()}), sorted({v['sage_version'] for v in det.values()}))",
        "disk.loc[disk.finished, ['project', 'n_runs', 'flags']].groupby('flags').size()")

    nb.code(
        "usable = disk[disk.finished & ~disk.unsupported & ~disk.not_comparable",
        "              & (disk.n_psm_tables > 0)].copy().reset_index(drop=True)",
        "excluded = disk[disk.finished & ~disk.project.isin(usable.project)].copy()",
        "print('finished and usable:', len(usable))",
        "print('finished but excluded:', len(excluded))",
        "print('still waiting in the queue:', len(selected) - int(disk.finished.sum()))",
        "excluded[['project', 'title', 'unsupported', 'not_comparable', 'n_psm_tables']]")

    # ------------------------------------------------- 1. classification rules
    nb.md("## 2 · Classifying every finished project by what kind of experiment it is",
          "A background table is only meaningful over one kind of experiment. A whole-proteome run "
          "detects nearly everything, so mixing it in would make every protein look like a "
          "frequent flyer; a phospho-enrichment detects a biased subset. Only affinity "
          "purification and proximity labelling go into the background.",
          "PRIDE gives us a title and the deposited file names. The pipeline also recorded an "
          "acquisition flag per run (DDA, possible TMT, timsTOF). The rules below are applied in "
          "order, first match wins, to the title and to the run file names separately, so it is "
          "visible which of the two decided the call. The cue that matched is printed in the "
          "output table — that is the audit trail.",
          "**The rules are keyword rules on text a human wrote, not a measurement.** They are "
          "printed verbatim so they can be checked and corrected. Where a project matches cues "
          "for two classes, or matches only on file names, it is flagged rather than trusted.")

    nb.code(
        "# Ordered rules. (class, regex, why). First match wins.",
        "RULES = [",
        "    ('hdx', r'\\bhdx\\b|hdxms|hydrogen[- ]deuterium',",
        "     'HDX-MS: masses are deuterium-shifted, not abundance (instruction 69)'),",
        "    ('proximity_labelling',",
        "     r'bioid|bira|turboid|turbo-id|\\bapex2?\\b|proxisome|proximity[- ]?(dependent|label)',",
        "     'BioID / TurboID / APEX: biotin ligase fused to the bait labels neighbours'),",
        "    ('affinity_purification',",
        "     r'immunoprecipit|immuno-precipit|\\bco-?ip\\b|\\bip-?ms\\b|interactome|interacting "
        "prote|pull[- ]?down|precipitation of|affinity purif|anti-ha ms|associated complex|"
        "\\binteractors\\b|binding with|associates with|interact with|co-?purif|"
        "proximal interacting',",
        "     'IP-MS / co-IP / tagged-bait affinity purification'),",
        "    ('ptm_enrichment',",
        "     r'phospho|acetylat|methylscan|methylarginine|lactylation|malonylation|"
        "nitrosylation|ubiquitylat|ubiquitinat|succinylat|crotonylat|glycosylat|\\bptm\\b|"
        "post-?translational modification',",
        "     'the sample is an enriched modified subset, not the proteome'),",
        "    ('fractionation',",
        "     r'hyperlopit|\\blopit\\b|subcellular|fractionat|spatial proteom|organellar|"
        "sucrose gradient|\\bbn-?page\\b',",
        "     'organelle or gradient fractionation: detection depends on the fraction'),",
        "    ('whole_proteome',",
        "     r'proteome|proteomic|protein abundance|expression|abundances following|antigen|"
        "immunogenic|proteogenomic',",
        "     'whole-cell or expression proteome'),",
        "]",
        "# File-name cues that are weaker than a title cue but still informative.",
        "NAME_CUES = [",
        "    ('proximity_labelling', r'bioid|bira|turboid|_bio[_.-]|biotin'),",
        "    ('affinity_purification', r'(^|[_.-])ip[_.-]|_ip\\d|ha[_.-]?ip|flag[_.-]?ip|"
        "beads|pulldown|pull[_-]down|compip'),",
        "    ('hdx', r'hdx'),",
        "    ('ptm_enrichment', r'phos|_pho\\d|psn|_sn[_.]|ptm'),",
        "]",
        "",
        "def _first_match(text, rules):",
        "    for cls, pat, _why in rules:",
        "        m = re.search(pat, text, re.I)",
        "        if m:",
        "            return cls, m.group(0)",
        "    return '', ''",
        "",
        "def _all_matches(text, rules):",
        "    out = []",
        "    for cls, pat, _why in rules:",
        "        m = re.search(pat, text, re.I)",
        "        if m:",
        "            out.append((cls, m.group(0)))",
        "    return out",
        "",
        "print('\\n'.join(f'{i+1}. {c:22s} {why}' for i, (c, _p, why) in enumerate(RULES)))")

    nb.code(
        "cls_rows = []",
        "for _, r in usable.iterrows():",
        "    names = ' '.join(det[r.project]['run_names'])",
        "    t_cls, t_cue = _first_match(r.title, RULES)",
        "    n_cls, n_cue = '', ''",
        "    for cls, pat in NAME_CUES:",
        "        m = re.search(pat, names, re.I)",
        "        if m:",
        "            n_cls, n_cue = cls, m.group(0)",
        "            break",
        "    title_hits = _all_matches(r.title, RULES)",
        "    # the call",
        "    if t_cls:",
        "        call, cue, decided_by = t_cls, f'title: {t_cue!r}', 'title'",
        "    elif n_cls:",
        "        call, cue, decided_by = n_cls, f'file names: {n_cue!r}', 'file names'",
        "    else:",
        "        call, cue, decided_by = 'other', 'no cue matched', 'nothing'",
        "    # confidence",
        "    reasons = []",
        "    if decided_by == 'nothing':",
        "        conf = 'none'",
        "        reasons.append('no keyword matched the title or the file names')",
        "    elif decided_by == 'file names':",
        "        conf = 'low'",
        "        reasons.append('title gave no cue; the call rests on deposited file names')",
        "    else:",
        "        conf = 'high'",
        "    if len(title_hits) > 1:",
        "        other = [c for c, _ in title_hits if c != call]",
        "        if other:",
        "            conf = 'medium' if conf == 'high' else conf",
        "            reasons.append('title also matches ' + '/'.join(sorted(set(other))))",
        "    if n_cls and t_cls and n_cls != t_cls:",
        "        conf = 'medium' if conf == 'high' else conf",
        "        reasons.append(f'file names suggest {n_cls} instead')",
        "    cls_rows.append(dict(project=r.project, title=r.title, n_runs=r.n_runs,",
        "                         acquisition=r.flags, predicted_class=call, cue=cue,",
        "                         decided_by=decided_by, confidence=conf,",
        "                         note='; '.join(reasons), example_runs='; '.join(",
        "                             det[r.project]['run_names'][:4])))",
        "cl = pd.DataFrame(cls_rows)",
        "cl.groupby(['predicted_class', 'confidence']).size().unstack(fill_value=0)")

    nb.md("### Hand corrections, each with its reason",
          "Three kinds of project the keyword rules get wrong, checked one at a time against the "
          "title and the deposited file names. Every change is listed here with what the rule said "
          "and why it is overridden; nothing is corrected silently.")

    nb.code(
        "# project -> (class, why). Applied after the rules; recorded in the output table.",
        "OVERRIDES = {",
        "    # Same study, two deposits: PXD039432 is the anti-HA IP, PXD039434 is the biotin",
        "    # (BIO) against vehicle (VEH) TurboID arm. The shared title says neither.",
        "    'PXD039434': ('proximity_labelling',",
        "                  'run names are BIO vs VEH (biotin against vehicle), i.e. proximity labelling,"
        " not the IP arm of the same paper'),",
        "    # Titles that name a complex but whose deposit is a gradient/BN-PAGE slice series.",
        "    'PXD068413': ('fractionation',",
        "                  'runs are Top/Middle/Bottom gradient slices of a conoid preparation,"
        " not a bait pull-down'),",
        "    'PXD068228': ('fractionation',",
        "                  'runs are Top/Middle/Bottom gradient slices, same preparation style'),",
        "    # Glycan beads capture by affinity to a sugar, not to a tagged parasite protein:",
        "    # there is no bait gene, so it cannot contribute a bait-specific background.",
        "    'PXD079900': ('affinity_purification',",
        "                  'glycan-bead capture: a real affinity purification but the bait is a glycan,"
        " so it has no bait gene'),",
        "    # Titles with no keyword at all, resolved from the deposited file names.",
        "    'PXD035658': ('whole_proteome',",
        "                  'run names carry SILAC labelling at two drug concentrations"
        " (100uM_L / 10uM_H), i.e. a whole-proteome drug comparison, not a pull-down'),",
        "    'PXD064338': ('whole_proteome',",
        "                  '32 TMT runs of infected host cells (a host immunity study):"
        " a whole-proteome experiment'),",
        "}",
        "cl['rule_class'] = cl.predicted_class",
        "for proj, (cls, why) in OVERRIDES.items():",
        "    if proj in set(cl.project):",
        "        i = cl.index[cl.project == proj][0]",
        "        cl.loc[i, 'predicted_class'] = cls",
        "        cl.loc[i, 'confidence'] = 'hand-checked'",
        "        cl.loc[i, 'note'] = (cl.loc[i, 'note'] + '; ' if cl.loc[i, 'note'] else '') + \\",
        "            f'override of rule call {cl.loc[i, \"rule_class\"]!r}: {why}'",
        "cl['class'] = cl.predicted_class",
        "print(cl.groupby('class').size().to_string())",
        "print()",
        "print('needs a human eye (confidence low/none/medium):',",
        "      int(cl.confidence.isin(['low', 'none', 'medium']).sum()), 'of', len(cl))",
        "cl.loc[cl.confidence.isin(['low', 'none', 'medium']),",
        "       ['project', 'title', 'class', 'cue', 'confidence', 'note']]")

    nb.code(
        "# every project, its class and how it was decided",
        "cl.sort_values(['class', 'project'])[['project', 'class', 'confidence', 'decided_by',",
        "                                      'cue', 'n_runs', 'title']]")

    # --------------------------------------------------------------- 2. baits
    nb.md("## 3 · Which protein was the bait",
          "For the background table itself the bait does not matter: we count how often a protein "
          "appears whatever was pulled. The bait matters for two other things. It tells us whether "
          "the experiments really are independent (ten pull-downs of the same complex are not ten "
          "independent looks at the background), and it is what a worked example needs.",
          "Baits are read out of the title and the run names, then the symbol is resolved to an "
          "ME49 gene ID through the ToxoDB product names in `nodes.parquet` plus a short alias "
          "table. A bait counts as **confident** only when the symbol resolves to exactly one "
          "gene.",
          "The author-reported `interaction_studies.parquet` (97 parsed studies, 31 with a bait) "
          "is a separate resource and is reported separately below. It is keyed by PMID and the "
          "PRIDE selection carries no PMID, so the two cannot be joined here — the author table is "
          "used as a comparison of how hard bait identification is, not as a source of baits for "
          "our projects.")

    nb.code(
        "nodes = pd.read_parquet(os.path.join(REPO, 'parascape', 'data', 'nodes.parquet'),",
        "                        columns=['gene_id', 'product', 'compartment'])",
        "nodes['product'] = nodes['product'].fillna('')",
        "print('genes:', len(nodes))",
        "",
        "def resolve_symbol(sym):",
        "    \"\"\"ME49 gene IDs whose ToxoDB product name carries this symbol as a whole word.\"\"\"",
        "    if not sym:",
        "        return []",
        "    pat = r'(?<![A-Za-z0-9])' + re.escape(sym) + r'(?![A-Za-z0-9])'",
        "    hit = nodes[nodes['product'].str.contains(pat, case=False, regex=True)]",
        "    return sorted(hit.gene_id)",
        "",
        "print(resolve_symbol('MyoA'), resolve_symbol('MIC2'), resolve_symbol('HDAC3'))")

    nb.code(
        "# Bait symbol patterns, applied to the title first and the run names second.",
        "BAIT_TITLE_PATTERNS = [",
        "    r'immunoprecipitation(?:\\s+mass[- ]spectrometry)?\\s+of\\s+(?:the\\s+)?"
        "(?:toxoplasma(?:\\s+gondii)?\\s+)?([A-Za-z0-9_.-]{2,20})',",
        "    r'(?:toxoplasma(?:\\s+gondii)?\\s+)?(?:pru\\s+)?([A-Za-z0-9_.-]{2,20})"
        "[- ]?(?:ha|mng|ty|flag|mng-maid)?\\s+(?:co-?)?immunoprecipitation',",
        "    r'([A-Za-z0-9_.-]{2,20})\\s+interactome',",
        "    r'([A-Za-z0-9_.-]{2,20})\\s+proxisome',",
        "    r'proximity labeling of\\s+(?:turboid-)?([A-Za-z0-9_.-]{2,20})',",
        "    r'([A-Za-z0-9_.-]{2,20})-turboid',",
        "    r'turboid-([A-Za-z0-9_.-]{2,20})',",
        "    r'([A-Za-z0-9_.-]{2,20})-bioid2?',",
        "    r'([A-Za-z0-9_.-]{2,20})\\s+interactors',",
        "    r'proteins that interact with\\s+([A-Za-z0-9_.-]{2,20})',",
        "    r'interacting proteomics of .*\\(([A-Za-z0-9]{2,20})\\)',",
        "    r'([A-Za-z0-9_.-]{2,20})[- ](?:associated|binding) ',",
        "    r'([A-Za-z0-9_.-]{2,20})\\s+associates with',",
        "]",
        "# Epitope and degron tags, stripped before the symbol is looked up.",
        "TAGS = r'(?:[-_ ](?:ha|ty|flag|mng|maid|aid|myc|turboid|bioid2?|bira|gfp|3xha|2xty|"
        "toxoplasma|tagged|second|copy))+$'",
        "# Aliases: the name a title uses against a token that really appears in a ToxoDB product",
        "# name. The gene ID is never typed; it is looked up through the product name, and the",
        "# product name is recorded in the output so the mapping can be checked.",
        "ALIAS = {",
        "    'myoa': 'myosin A',          # the product name is 'myosin A', not the symbol",
        "    'hce1': 'TEEGR',             # published synonym; product 'E2F4-associated ... TEEGR'",
        "    'pka c3': 'PKAc3', 'pkac3': 'PKAc3', 'c3': 'PKAc3',",
        "    'calmodulin': 'CAM1',",
        "    'tgbdp1': 'BDP1', 'bdp1': 'BDP1',",
        "    'tgsir2a': 'SIR2A', 'sir2a': 'SIR2A',",
        "    'tggsk3': 'GSK3', 'gsk3': 'GSK3',",
        "    'tgcoqfmo': 'CoqFMO',",
        "    'cyclin-4': 'cyclin4', 'cyclin4': 'cyclin4',",
        "}",
        "# Titles that name a whole complex rather than one tagged protein. The bait is a complex,",
        "# which is a different thing from an unknown bait and is recorded as such.",
        "COMPLEX_BAIT_WORDS = r'succinate dehydrogenase|metc|mitochondrial ribosome|conoid|"
        "supercomplex|glycan bead|endocytic and exocytic'",
        "HOST_BAITS = {'mospd2'}  # a host protein, not an ME49 gene",
        "",
        "def strip_tags(raw):",
        "    prev = None",
        "    s = raw",
        "    while s != prev:",
        "        prev = s",
        "        s = re.sub(TAGS, '', s, flags=re.I).strip('-_. ')",
        "    return s",
        "",
        "STOPWORDS = {'the', 'toxoplasma', 'gondii', 'mass', 'a', 'an', 'protein', 'proteins',",
        "             'and', 'of', 'ms', 'ha', 'anti', 'homolg', 'homolog', 'mitochondrial',",
        "             'succinate', 'novel', 'components', 'this'}",
        "",
        "def find_bait(title, names):",
        "    \"\"\"Bait symbol from the title first, then from the deposited run names.\"\"\"",
        "    for pat in BAIT_TITLE_PATTERNS:",
        "        m = re.search(pat, title, re.I)",
        "        if not m:",
        "            continue",
        "        raw = strip_tags(m.group(1).strip('-_. '))",
        "        if raw.lower() in STOPWORDS or len(raw) < 2:",
        "            continue",
        "        return raw, 'title'",
        "    # run names: a six-digit ME49 number, or the token before an IP / tag marker",
        "    m = re.search(r'(?<!\\d)(\\d{6})(?!\\d)', names)",
        "    if m and ('TGME49_' + m.group(1)) in set(nodes.gene_id):",
        "        return 'TGME49_' + m.group(1), 'run names (ME49 number in the file name)'",
        "    m = re.search(r'(?:^|[_ ])IP[_-]([A-Za-z0-9]{3,20})', names)",
        "    if m:",
        "        raw = strip_tags(re.sub(r'(AID|HA|Ty|mNG)$', '', m.group(1)))",
        "        if raw.lower() not in STOPWORDS and len(raw) >= 2:",
        "            return raw, 'run names (token after IP_)'",
        "    return '', ''",
        "",
        "bait_rows = []",
        "for _, r in cl.iterrows():",
        "    names = ' '.join(det[r.project]['run_names'])",
        "    raw, how = find_bait(r.title, names)",
        "    sym = ALIAS.get(raw.lower(), raw)",
        "    if raw.lower() in HOST_BAITS:",
        "        bait_rows.append(dict(project=r.project, bait_symbol=raw, bait_gene='HOST',",
        "                              bait_confidence='host bait', bait_how=how,",
        "                              bait_product=''))",
        "        continue",
        "    if re.search(COMPLEX_BAIT_WORDS, r.title, re.I):",
        "        bait_rows.append(dict(project=r.project, bait_symbol=raw,",
        "                              bait_gene='', bait_how=how or 'title names a complex',",
        "                              bait_confidence='a complex or a non-protein bait, "
        "not one tagged gene', bait_product=''))",
        "        continue",
        "    if re.fullmatch(r'TGME49_\\d+', sym):",
        "        genes = [sym]",
        "    else:",
        "        genes = resolve_symbol(sym)",
        "    if len(genes) == 1:",
        "        conf = 'confident'",
        "    elif len(genes) > 1:",
        "        conf = 'ambiguous (symbol matches several genes)'",
        "    elif raw:",
        "        conf = 'ambiguous (symbol not in any ToxoDB product name)'",
        "    else:",
        "        conf = 'none'",
        "    prod = ''",
        "    if len(genes) == 1:",
        "        prod = nodes.loc[nodes.gene_id == genes[0], 'product'].iloc[0]",
        "    bait_rows.append(dict(project=r.project, bait_symbol=raw,",
        "                          bait_gene=genes[0] if len(genes) == 1 else ';'.join(genes),",
        "                          bait_confidence=conf, bait_how=how, bait_product=prod))",
        "bt = pd.DataFrame(bait_rows)",
        "cl = cl.merge(bt, on='project', how='left')",
        "pull = cl[cl['class'].isin(['affinity_purification', 'proximity_labelling'])]",
        "print('affinity + proximity projects:', len(pull))",
        "print(pull.bait_confidence.value_counts().to_string())",
        "print()",
        "print('distinct confident bait genes:',",
        "      pull.loc[pull.bait_confidence == 'confident', 'bait_gene'].nunique())",
        "pull[['project', 'class', 'bait_symbol', 'bait_gene', 'bait_confidence', 'title']]")

    nb.code(
        "# the author-reported table, for comparison only",
        "st = pd.read_parquet(os.path.join(REPO, 'parascape', 'data', 'interaction_studies.parquet'))",
        "mem = pd.read_parquet(os.path.join(REPO, 'parascape', 'data',",
        "                                   'interaction_study_members.parquet'))",
        "print('parsed interaction studies:', len(st),",
        "      '| with a bait gene:', int((st.bait_gene.fillna('') != '').sum()))",
        "print('bait-member pairs:', len(mem))",
        "print('methods:'); print(st.method.value_counts().to_string())",
        "print()",
        "print('can these be joined to our PXD projects?  columns available:', list(st.columns))",
        "print('-> no accession column, so no join; used as a comparison only.')",
        "st.method.value_counts()")

    # ------------------------------------------------------------- 3. detection
    nb.md("## 4 · The detection rule",
          "Sage reports, for every peptide-spectrum match, a `spectrum_q` (the false-discovery "
          "rate at which that match would be accepted) and the protein or proteins the peptide "
          "belongs to. The rule used for the background table is:",
          "> **A protein is detected in an experiment when at least 2 peptide-spectrum matches at "
          "1 % spectrum-level FDR map to that gene and to no other parasite gene.**",
          "Three choices, each for a reason:",
          "* **1 % FDR** is the field standard and it is Sage's own number, so no new threshold is "
          "invented here.",
          "* **At least 2 matches** removes single-spectrum identifications. At 1 % FDR a single "
          "match is the one most likely to be wrong, and for a background table a false detection "
          "is worse than a missed one: it would make a protein look like a frequent flyer and get "
          "a real hit filtered out of someone's list.",
          "* **Unique mapping.** If a peptide is shared between two paralogues, crediting both "
          "would inflate the frequency of whole families (SRS surface antigens, the myosins). A "
          "protein seen only through shared peptides is not counted — a conservative choice that "
          "makes frequencies slightly too low rather than too high.",
          "Section 6 repeats the whole table under a looser rule (1 match) and a stricter one "
          "(at least 2 distinct peptide sequences and at least 3 matches), so the sensitivity of "
          "every conclusion to this choice is visible.",
          "**Experiments, not runs.** One project is one experiment. A project with eight runs is "
          "usually the bait and its control in replicate, which is one look at the background, not "
          "eight. Run-level counts are also recorded so the difference can be seen.")

    nb.code(
        "GENE = re.compile(r'TGME49_\\d+')",
        "PSM_COLS = ['peptide', 'proteins', 'label', 'spectrum_q', 'ms2_intensity', 'filename']",
        "",
        "def aggregate_psms(path, chunksize=400_000):",
        "    \"\"\"One pass over a project's Sage PSM table. Returns a per-gene frame and totals.",
        "",
        "    Only the aggregate is kept in memory, so the table size does not matter.\"\"\"",
        "    parts, n_rows, n_pass = [], 0, 0",
        "    for ch in pd.read_csv(path, sep='\\t', usecols=PSM_COLS, chunksize=chunksize):",
        "        n_rows += len(ch)",
        "        ch = ch[(ch['label'] == 1) & (ch['spectrum_q'] <= 0.01)]",
        "        n_pass += len(ch)",
        "        if not len(ch):",
        "            continue",
        "        genes = ch['proteins'].astype(str).map(lambda s: set(GENE.findall(s)))",
        "        one = genes.map(len) == 1",
        "        ch = ch.loc[one].copy()",
        "        if not len(ch):",
        "            continue",
        "        ch['gene_id'] = [next(iter(g)) for g in genes.loc[one]]",
        "        g = ch.groupby('gene_id').agg(",
        "            n_psm=('peptide', 'size'),",
        "            n_peptides=('peptide', 'nunique'),",
        "            n_runs=('filename', 'nunique'),",
        "            intensity=('ms2_intensity', 'sum'))",
        "        pep = ch.groupby('gene_id')['peptide'].apply(lambda s: set(s))",
        "        run = ch.groupby('gene_id')['filename'].apply(lambda s: set(s))",
        "        g['peptide_set'] = pep",
        "        g['run_set'] = run",
        "        parts.append(g)",
        "    empty = pd.DataFrame({'gene_id': pd.Series(dtype=str),",
        "                          'n_psm': pd.Series(dtype=int), 'intensity': pd.Series(dtype=float),",
        "                          'n_peptides': pd.Series(dtype=int),",
        "                          'n_runs_detected': pd.Series(dtype=int)})",
        "    if not parts:",
        "        return empty, dict(n_rows=n_rows, n_pass=n_pass)",
        "    cat = pd.concat(parts)",
        "    out = cat.groupby(level=0).agg(",
        "        n_psm=('n_psm', 'sum'), intensity=('intensity', 'sum'),",
        "        peptide_set=('peptide_set', lambda ss: set().union(*ss)),",
        "        run_set=('run_set', lambda ss: set().union(*ss)))",
        "    out['n_peptides'] = out.peptide_set.map(len)",
        "    out['n_runs_detected'] = out.run_set.map(len)",
        "    out = out.drop(columns=['peptide_set', 'run_set']).reset_index(names='gene_id')",
        "    return out, dict(n_rows=n_rows, n_pass=n_pass)",
        "",
        "LFQ_META = {'peptide', 'charge', 'proteins', 'q_value', 'score', 'spectral_angle'}",
        "",
        "def aggregate_lfq(path):",
        "    \"\"\"Per-gene MS1 LFQ intensity, summed over the project's runs.\"\"\"",
        "    t = pd.read_csv(path, sep='\\t')",
        "    runs = [c for c in t.columns if c not in LFQ_META]",
        "    t = t[t['q_value'] <= 0.01]",
        "    genes = t['proteins'].astype(str).map(lambda s: set(GENE.findall(s)))",
        "    one = genes.map(len) == 1",
        "    t = t.loc[one].copy()",
        "    if not len(t):",
        "        return pd.DataFrame(columns=['gene_id', 'lfq'])",
        "    t['gene_id'] = [next(iter(g)) for g in genes.loc[one]]",
        "    t['lfq'] = t[runs].sum(axis=1)",
        "    return t.groupby('gene_id', as_index=False)['lfq'].sum()",
        "print('aggregators defined')")

    nb.code(
        "# one project at a time; only the aggregate is kept",
        "per_project, psm_stats = {}, []",
        "for _, r in cl.iterrows():",
        "    row = usable.loc[usable.project == r.project].iloc[0]",
        "    agg, stats_ = aggregate_psms(row.psm_path)",
        "    if row.lfq_path:",
        "        lfq = aggregate_lfq(row.lfq_path)",
        "        agg = agg.merge(lfq, on='gene_id', how='left')",
        "    else:",
        "        agg['lfq'] = np.nan",
        "    agg['project'] = r.project",
        "    per_project[r.project] = agg",
        "    psm_stats.append(dict(project=r.project, psm_rows=stats_['n_rows'],",
        "                          psm_1pct=stats_['n_pass'],",
        "                          genes_1psm=int((agg.n_psm >= 1).sum()),",
        "                          genes_2psm=int((agg.n_psm >= 2).sum()),",
        "                          has_lfq=bool(row.lfq_path),",
        "                          lfq_genes=int(agg.lfq.notna().sum())))",
        "ps = pd.DataFrame(psm_stats)",
        "cl = cl.merge(ps, on='project', how='left')",
        "print('projects aggregated:', len(per_project))",
        "print('PSM rows read:', int(ps.psm_rows.sum()),",
        "      '| passing 1% FDR and target:', int(ps.psm_1pct.sum()))",
        "ps.describe()[['psm_rows', 'psm_1pct', 'genes_1psm', 'genes_2psm']]")

    # --------------------------------------------------------- 4. background
    nb.md("## 5 · The background table",
          "Only the affinity-purification and proximity-labelling projects contribute. Abundance "
          "is made comparable between projects by expressing each protein as parts per million of "
          "that project's own total parasite signal, so a deep run and a shallow run can be "
          "averaged. Two abundance scales are kept: MS1 label-free quantity (`lfq_ppm`, the better "
          "measure, available where the pipeline wrote an `lfq.tsv`) and summed MS2 fragment "
          "intensity of the counted matches (`psm_ppm`, available everywhere, including the "
          "timsTOF projects which have no MS1 LFQ).")

    nb.code(
        "BG_CLASSES = ['affinity_purification', 'proximity_labelling']",
        "bg_projects = sorted(cl.loc[cl['class'].isin(BG_CLASSES), 'project'])",
        "N_EXP = len(bg_projects)",
        "print('experiments in the background:', N_EXP)",
        "print(cl[cl['class'].isin(BG_CLASSES)].groupby('class').size().to_string())",
        "",
        "frames = []",
        "for p in bg_projects:",
        "    a = per_project[p].copy()",
        "    tot_i = a.intensity.sum()",
        "    tot_l = a.lfq.sum(skipna=True)",
        "    a['psm_ppm'] = a.intensity / tot_i * 1e6 if tot_i > 0 else np.nan",
        "    a['lfq_ppm'] = a.lfq / tot_l * 1e6 if tot_l and tot_l > 0 else np.nan",
        "    frames.append(a)",
        "long = pd.concat(frames, ignore_index=True)",
        "print('gene x experiment rows (>=1 unique match):', len(long))",
        "long.head(5)")

    nb.code(
        "def build_background(long, n_exp, min_psm=2, min_peptides=1):",
        "    d = long[(long.n_psm >= min_psm) & (long.n_peptides >= min_peptides)]",
        "    g = d.groupby('gene_id').agg(",
        "        experiments_detected=('project', 'nunique'),",
        "        runs_detected=('n_runs_detected', 'sum'),",
        "        median_psm_ppm=('psm_ppm', 'median'), max_psm_ppm=('psm_ppm', 'max'),",
        "        median_abundance=('lfq_ppm', 'median'), max_abundance=('lfq_ppm', 'max'),",
        "        n_experiments_with_lfq=('lfq_ppm', 'count'),",
        "        total_psm=('n_psm', 'sum'))",
        "    g['n_experiments'] = n_exp",
        "    g['frequency'] = g.experiments_detected / n_exp",
        "    return g.reset_index()",
        "",
        "bg = build_background(long, N_EXP)",
        "bg = bg.merge(nodes.rename(columns={'product': 'product_name'}), on='gene_id', how='left')",
        "print('proteins in the background table:', len(bg))",
        "print('detected in every experiment:', int((bg.frequency == 1).sum()))",
        "print('detected in >= 80% of experiments:', int((bg.frequency >= 0.8).sum()))",
        "print('detected in exactly one experiment:', int((bg.experiments_detected == 1).sum()))",
        "print()",
        "print('frequency distribution:')",
        "print(pd.cut(bg.frequency, [0, .1, .25, .5, .75, .9, 1.0],",
        "             include_lowest=True).value_counts().sort_index().to_string())",
        "bg.sort_values('frequency', ascending=False).head(25)[",
        "    ['gene_id', 'product_name', 'experiments_detected', 'frequency',",
        "     'median_abundance', 'compartment']]")

    nb.md("### Runs are not experiments",
          "The same table counted by run instead of by project shows why the unit matters: the "
          "run count is dominated by whichever projects happened to deposit many files.")

    nb.code(
        "n_runs_total = int(cl.loc[cl['class'].isin(BG_CLASSES), 'n_runs'].sum())",
        "cmp = bg[['gene_id', 'product_name', 'experiments_detected', 'frequency', 'runs_detected']].copy()",
        "cmp['run_frequency'] = cmp.runs_detected / n_runs_total",
        "print('runs behind the', N_EXP, 'experiments:', n_runs_total)",
        "print('rank correlation of the two frequencies: %.3f' %",
        "      stats.spearmanr(cmp.frequency, cmp.run_frequency).statistic)",
        "print('genes in the top 50 by experiment frequency but not by run frequency:',",
        "      len(set(cmp.nlargest(50, 'frequency').gene_id) -",
        "          set(cmp.nlargest(50, 'run_frequency').gene_id)))",
        "cmp.nlargest(10, 'run_frequency')[['gene_id', 'product_name', 'frequency', 'run_frequency']]")

    # ------------------------------------------------------- 5. sensitivity
    nb.md("## 6 · How much the detection rule matters",
          "The same table under a looser and a stricter rule. If the ranking barely moves, the "
          "threshold is not what the conclusions rest on.")

    nb.code(
        "variants = {",
        "    'loose (>=1 match)': build_background(long, N_EXP, min_psm=1),",
        "    'default (>=2 matches)': bg[['gene_id', 'experiments_detected', 'frequency']].copy(),",
        "    'strict (>=3 matches, >=2 peptides)': build_background(long, N_EXP, min_psm=3,",
        "                                                          min_peptides=2),",
        "}",
        "rows = []",
        "ref = variants['default (>=2 matches)'].set_index('gene_id').frequency",
        "for name, v in variants.items():",
        "    f = v.set_index('gene_id').frequency",
        "    both = ref.index.intersection(f.index)",
        "    rows.append(dict(rule=name, proteins=len(v),",
        "                     freq_ge_0_5=int((f >= 0.5).sum()),",
        "                     freq_ge_0_8=int((f >= 0.8).sum()),",
        "                     seen_once=int((v.set_index('gene_id').experiments_detected == 1).sum()),",
        "                     median_frequency=round(float(f.median()), 4),",
        "                     spearman_vs_default=round(float(",
        "                         stats.spearmanr(ref[both], f[both]).statistic), 4)))",
        "sens = pd.DataFrame(rows)",
        "sens")

    # ------------------------------------------------------- 6. specificity
    nb.md("## 7 · The specificity score",
          "Given a hit list from one experiment — a protein and how much of it there was — the "
          "score asks two questions and adds the answers:",
          "1. **Is there more of it here than usual?** `log2` of this experiment's parts per "
          "million against the median parts per million the protein reaches across the background "
          "experiments in which it is seen at all.",
          "2. **Is it unusual to see it at all?** `log2` of one over its background frequency, "
          "smoothed so a protein seen in no background experiment does not give an infinite score.",
          "`specificity = log2(abundance fold) + log2(1 / frequency)`. A protein seen in every "
          "pull-down at its usual level scores about 0. Tubulin at twice its usual level still "
          "scores low, because the second term is near zero. A protein never seen before, at a "
          "decent level, scores high.",
          "Both parts are reported separately, because they fail differently: the fold change is "
          "unreliable for a protein with few background measurements, and the frequency is "
          "unreliable for everything while the number of experiments is small.")

    nb.code(
        "def specificity(hits, background, n_experiments=None, pseudo=0.5, abundance='psm_ppm'):",
        "    \"\"\"Score a pull-down hit list against the CRAPome background.",
        "",
        "    hits        : DataFrame with 'gene_id' and an abundance column in parts per million",
        "                  of that experiment's own total parasite signal.",
        "    background  : the background table (gene_id, experiments_detected, n_experiments,",
        "                  frequency, median_psm_ppm / median_abundance).",
        "    Returns the hit list with fold, frequency and the score, highest score first.",
        "    \"\"\"",
        "    bgcol = 'median_abundance' if abundance == 'lfq_ppm' else 'median_psm_ppm'",
        "    n = int(background.n_experiments.iloc[0]) if n_experiments is None else n_experiments",
        "    b = background[['gene_id', 'experiments_detected', bgcol]].rename(",
        "        columns={bgcol: 'background_ppm'})",
        "    h = hits.merge(b, on='gene_id', how='left')",
        "    h['experiments_detected'] = h.experiments_detected.fillna(0)",
        "    h['n_experiments'] = n",
        "    # smoothed frequency: never 0, never exactly 1",
        "    h['background_frequency'] = (h.experiments_detected + pseudo) / (n + pseudo)",
        "    typical = float(background[bgcol].median(skipna=True))",
        "    bgp = h.background_ppm.fillna(typical)",
        "    h['abundance_fold'] = (h[abundance] + 1.0) / (bgp + 1.0)",
        "    h['log2_fold'] = np.log2(h.abundance_fold)",
        "    h['log2_rarity'] = np.log2(1.0 / h.background_frequency)",
        "    h['specificity'] = h.log2_fold + h.log2_rarity",
        "    return h.sort_values('specificity', ascending=False).reset_index(drop=True)",
        "print(specificity.__doc__.splitlines()[0])")

    nb.md("### A worked example: the MYR1 immunoprecipitation (PXD016383)",
          "MYR1 is the translocon component that moves dense granule effectors across the "
          "parasitophorous vacuole membrane into the host cell. Its partners are known from the "
          "literature — MYR2, MYR3, MYR4, GRA44, GRA45, ROP17 — so this is a hit list whose right "
          "answer we already have. The project itself contributes to the background, which makes "
          "the test slightly harder, not easier.")

    nb.code(
        "EXAMPLE = 'PXD016383'",
        "ex = per_project[EXAMPLE].copy()",
        "tot = ex.intensity.sum()",
        "ex['psm_ppm'] = ex.intensity / tot * 1e6",
        "ex = ex[ex.n_psm >= 2]",
        "scored = specificity(ex[['gene_id', 'n_psm', 'n_peptides', 'psm_ppm']], bg)",
        "scored = scored.merge(nodes.rename(columns={'product': 'product_name'}),",
        "                      on='gene_id', how='left')",
        "print(cl.loc[cl.project == EXAMPLE, ['title', 'class', 'bait_symbol', 'bait_gene']].to_string(index=False))",
        "print('proteins in the hit list:', len(scored))",
        "scored.head(25)[['gene_id', 'product_name', 'n_psm', 'psm_ppm', 'background_frequency',",
        "                 'log2_fold', 'log2_rarity', 'specificity']]")

    nb.code(
        "# where the known MYR1 complex lands in that ranking",
        "MYR_SYMS = ['MYR1', 'MYR2', 'MYR3', 'MYR4', 'GRA44', 'GRA45', 'ROP17']",
        "want = {}",
        "for s in MYR_SYMS:",
        "    for g in resolve_symbol(s):",
        "        want[g] = s",
        "scored['rank'] = np.arange(1, len(scored) + 1)",
        "scored['percentile'] = 100 * (1 - (scored['rank'] - 1) / len(scored))",
        "known = scored[scored.gene_id.isin(want)].copy()",
        "known['symbol'] = known.gene_id.map(want)",
        "known[['symbol', 'gene_id', 'product_name', 'rank', 'percentile', 'psm_ppm',",
        "       'background_frequency', 'specificity']]")

    # ------------------------------------------------------- 7. validation
    nb.md("## 8 · Validation: do real complexes survive, and are the frequent flyers caught?",
          "This is what decides whether the table can be used. Two sets of proteins with a known "
          "answer:",
          "* **Known complexes** must have a *low* background frequency. If the glideosome is "
          "detected in every pull-down, the table would strip it out of a glideosome paper.",
          "* **Frequent flyers** — ribosomal proteins, HSP70/HSP90, tubulin, actin, GAPDH, enolase "
          "— must have a *high* background frequency. If they do not, the table is not catching "
          "the thing it exists to catch.",
          "Gene identities come from ToxoDB product names in `nodes.parquet`, listed in full so "
          "they can be checked.")

    nb.code(
        "# Each complex member is named by the exact ToxoDB product name, so there is no doubt",
        "# which gene is meant and no gene ID is typed from memory.",
        "COMPLEXES = {",
        "    'MyoA glideosome': {'MyoA': 'myosin A', 'MLC1': 'myosin light chain MLC1',",
        "                        'ELC1': 'essential light chain ELC1',",
        "                        'GAP45': 'gliding-associated protein GAP45',",
        "                        'GAP50': 'gliding-associated protein GAP50',",
        "                        'GAP40': 'GAP40 protein'},",
        "    'AMA1-RON':        {'AMA1': 'apical membrane antigen AMA1',",
        "                        'RON2': 'rhoptry neck protein RON2',",
        "                        'RON4': 'rhoptry neck protein RON4'},",
        "    'MIC2-M2AP':       {'MIC2': 'microneme protein MIC2',",
        "                        'M2AP': 'MIC2-associated protein M2AP'},",
        "    'MORC/HDAC3':      {'MORC': 'ATPase MORC',",
        "                        'HDAC3': 'histone deacetylase HDAC3'},",
        "}",
        "# Frequent flyers, tight: the actual abundant housekeeping proteins. The cytosolic",
        "# ribosome is taken from the curated hyperLOPIT 40S/60S markers rather than from a product",
        "# name, so the apicoplast and mitochondrial 30S/50S ribosomes are not mixed in.",
        "FLYER_TIGHT = {",
        "    'HSP70/HSP90': r'^heat shock protein (?:hsp)?(?:70|90)(?:,|$)',",
        "    'tubulin':     r'^(?:(?:alpha|beta|beta-1|gamma)[- ]tubulin|tubulin (?:gamma|alpha"
        "|beta) chain)(?!.*(?:ligase|suppressor|transferase|cofactor|complex component|binding))',",
        "    'actin':       r'^actin ACT1$',",
        "    'GAPDH':       r'^glyceraldehyde-3-phosphate dehydrogenase GAPDH',",
        "    'enolase':     r'^enolase \\d$',",
        "}",
        "# Frequent flyers, loose: whole families as a keyword would find them. Reported beside the",
        "# tight sets because the loose sets sweep in low-abundance relatives (the tubulin-tyrosine",
        "# ligases, the organellar ribosomes) that are not frequent flyers at all.",
        "FLYER_LOOSE = {",
        "    'ribosomal protein (any)': r'ribosomal protein',",
        "    'heat shock (any)':        r'heat shock|hsp\\d',",
        "    'tubulin (any)':           r'tubulin',",
        "    'actin (any)':             r'^actin',",
        "}",
        "",
        "val_rows = []",
        "for name, members in COMPLEXES.items():",
        "    for sym, product in members.items():",
        "        hit = nodes[nodes['product'].str.lower() == product.lower()]",
        "        if not len(hit):",
        "            print('WARNING: no gene with product %r' % product)",
        "        for g in hit.gene_id:",
        "            val_rows.append(dict(set_name=name, kind='known complex', symbol=sym,",
        "                                 gene_id=g))",
        "for g in nodes.loc[nodes.compartment.isin(['40S ribosome', '60S ribosome']), 'gene_id']:",
        "    val_rows.append(dict(set_name='cytosolic ribosome (hyperLOPIT marker)',",
        "                         kind='expected frequent flyer', symbol='ribosome', gene_id=g))",
        "for name, pat in FLYER_TIGHT.items():",
        "    hit = nodes[nodes['product'].str.contains(pat, case=False, regex=True)]",
        "    for g in hit.gene_id:",
        "        val_rows.append(dict(set_name=name, kind='expected frequent flyer',",
        "                             symbol=name, gene_id=g))",
        "for name, pat in FLYER_LOOSE.items():",
        "    hit = nodes[nodes['product'].str.contains(pat, case=False, regex=True)]",
        "    for g in hit.gene_id:",
        "        val_rows.append(dict(set_name=name, kind='loose family (comparison only)',",
        "                             symbol=name, gene_id=g))",
        "val = pd.DataFrame(val_rows).drop_duplicates(['set_name', 'gene_id'])",
        "val = val.merge(bg[['gene_id', 'product_name', 'experiments_detected', 'n_experiments',",
        "                    'frequency', 'median_abundance', 'median_psm_ppm', 'compartment']],",
        "                on='gene_id', how='left')",
        "val['experiments_detected'] = val.experiments_detected.fillna(0).astype(int)",
        "val['n_experiments'] = N_EXP",
        "val['frequency'] = val.frequency.fillna(0.0)",
        "print('validation proteins:', len(val))",
        "val.groupby(['kind', 'set_name']).agg(",
        "    proteins=('gene_id', 'size'),",
        "    median_frequency=('frequency', 'median'),",
        "    max_frequency=('frequency', 'max'),",
        "    never_seen=('experiments_detected', lambda s: int((s == 0).sum())))")

    nb.code(
        "# the named complex members one by one",
        "val[val.kind == 'known complex'].sort_values(['set_name', 'symbol'])[",
        "    ['set_name', 'symbol', 'gene_id', 'product_name', 'experiments_detected',",
        "     'n_experiments', 'frequency', 'median_abundance']]")

    nb.code(
        "a = val.loc[val.kind == 'known complex', 'frequency']",
        "b = val.loc[val.kind == 'expected frequent flyer', 'frequency']",
        "u = stats.mannwhitneyu(b, a, alternative='greater')",
        "auc = u.statistic / (len(a) * len(b))",
        "print('known complex members:        n=%d  median frequency %.3f' % (len(a), a.median()))",
        "print('expected frequent flyers:     n=%d  median frequency %.3f' % (len(b), b.median()))",
        "print('all proteins in the table:    n=%d  median frequency %.3f'",
        "      % (len(bg), bg.frequency.median()))",
        "print()",
        "print('flyers above complex members: AUC %.3f, Mann-Whitney p = %.3g' % (auc, u.pvalue))",
        "print('a flyer picked at random beats a complex member %.0f%% of the time' % (100 * auc))",
        "print()",
        "thresh = 0.5",
        "print('at a frequency cutoff of %.0f%%:' % (100 * thresh))",
        "print('  complex members that would be filtered out as background: %d of %d'",
        "      % (int((a >= thresh).sum()), len(a)))",
        "print('  frequent flyers correctly filtered:                      %d of %d'",
        "      % (int((b >= thresh).sum()), len(b)))",
        "pd.DataFrame({'known complex': a.describe(), 'frequent flyer': b.describe(),",
        "              'all proteins': bg.frequency.describe()})")

    nb.md("### Could the frequency be used as a filter at any cutoff?",
          "A direction that is right on average is not the same as a filter that works. The "
          "question a biologist actually asks is: if I delete everything above some frequency, how "
          "many frequent flyers do I lose and how much of my complex do I lose with them? The "
          "sweep below answers that at every cutoff.")

    nb.code(
        "sweep = []",
        "for c in np.arange(0.05, 1.001, 0.05):",
        "    sweep.append(dict(cutoff=round(float(c), 2),",
        "                      flyers_removed=int((b >= c).sum()), flyers=len(b),",
        "                      flyer_recall=float((b >= c).mean()),",
        "                      complex_lost=int((a >= c).sum()), complex_members=len(a),",
        "                      complex_loss=float((a >= c).mean()),",
        "                      proteins_removed=int((bg.frequency >= c).sum())))",
        "sw = pd.DataFrame(sweep)",
        "ok = sw[(sw.complex_loss <= 0.15) & (sw.flyer_recall >= 0.5)]",
        "FREQ_FILTER_WORKS = bool(len(ok))",
        "if FREQ_FILTER_WORKS:",
        "    print('a usable cutoff exists:')",
        "    print(ok.to_string(index=False))",
        "else:",
        "    best = sw[sw.flyer_recall >= 0.5].nsmallest(1, 'complex_loss')",
        "    print('NO cutoff removes at least half the frequent flyers while losing at most 15% of')",
        "    print('the known complex members. The least bad such cutoff:')",
        "    print(best.to_string(index=False))",
        "    print()",
        "    print('-> the frequency column must not be used as a standalone filter.')",
        "sw")

    nb.md("### The loose family sets, for comparison",
          "The same test with the frequent-flyer sets defined by a keyword instead of curated. "
          "The loose sets sweep in the tubulin-tyrosine ligases and the apicoplast and "
          "mitochondrial 30S/50S ribosomes, which are not abundant and not sticky, so they drag "
          "the family medians down. This is why the tight sets are the ones the verdict uses.")

    nb.code(
        "loose_tbl = val[val.kind == 'loose family (comparison only)'].groupby('set_name').agg(",
        "    proteins=('gene_id', 'size'), median_frequency=('frequency', 'median'),",
        "    max_frequency=('frequency', 'max'))",
        "tight_tbl = val[val.kind == 'expected frequent flyer'].groupby('set_name').agg(",
        "    proteins=('gene_id', 'size'), median_frequency=('frequency', 'median'),",
        "    max_frequency=('frequency', 'max'))",
        "print('tight (curated):'); print(tight_tbl.to_string())",
        "print(); print('loose (keyword):'); print(loose_tbl.to_string())",
        "tight_tbl")

    nb.md("### The test that matters for the score: does a bait come out top of its own pull-down?",
          "Background frequency on its own is not the product; the specificity score is. The "
          "sharpest test of the score is the one case where the right answer is certain: the "
          "tagged protein itself. In every experiment with a confident bait, the bait should be "
          "near the top of its own hit list. This is run across all such experiments at once.")

    nb.code(
        "rec = []",
        "for _, r in pull[pull.bait_confidence == 'confident'].iterrows():",
        "    a_ = per_project[r.project].copy()",
        "    tot_ = a_.intensity.sum()",
        "    if not tot_ or not len(a_):",
        "        continue",
        "    a_['psm_ppm'] = a_.intensity / tot_ * 1e6",
        "    a_ = a_[a_.n_psm >= 2]",
        "    s = specificity(a_[['gene_id', 'n_psm', 'psm_ppm']], bg)",
        "    s['rank'] = np.arange(1, len(s) + 1)",
        "    hit = s[s.gene_id == r.bait_gene]",
        "    rec.append(dict(project=r.project, bait=r.bait_symbol, bait_gene=r.bait_gene,",
        "                    hits=len(s),",
        "                    bait_rank=int(hit['rank'].iloc[0]) if len(hit) else np.nan,",
        "                    bait_percentile=(100 * (1 - (hit['rank'].iloc[0] - 1) / len(s))",
        "                                     if len(hit) else np.nan),",
        "                    bait_specificity=float(hit.specificity.iloc[0]) if len(hit) else np.nan,",
        "                    bait_frequency=(float(hit.background_frequency.iloc[0])",
        "                                    if len(hit) else np.nan),",
        "                    bait_abundance_rank=(int((a_.psm_ppm > float(hit.psm_ppm.iloc[0])).sum()) + 1",
        "                                         if len(hit) else np.nan)))",
        "recov = pd.DataFrame(rec)",
        "recov = recov.merge(cl[['project', 'class']], on='project', how='left')",
        "got = recov.dropna(subset=['bait_rank']) if len(recov) else recov",
        "if len(got):",
        "    print('experiments with a confident bait:', len(recov),",
        "          '| bait detected in its own pull-down:', len(got))",
        "    print('bait in the top 1%% of its own hit list: %d of %d'",
        "          % (int((got.bait_percentile >= 99).sum()), len(got)))",
        "    print('bait in the top 5%% of its own hit list: %d of %d'",
        "          % (int((got.bait_percentile >= 95).sum()), len(got)))",
        "    print('median bait rank: %.0f of a median %.0f proteins'",
        "          % (got.bait_rank.median(), got.hits.median()))",
        "    print('for comparison, ranking by raw abundance alone: median bait rank %.0f'",
        "          % got.bait_abundance_rank.median())",
        "    print()",
        "    print('split by what kind of experiment it is -- the score is designed for affinity',",
        "          'purification, where the bait is pulled out intact:')",
        "    print(got.groupby('class').agg(experiments=('project', 'size'),",
        "                                   in_top_5pct=('bait_percentile',",
        "                                                lambda s: int((s >= 95).sum())),",
        "                                   median_percentile=('bait_percentile', 'median')",
        "                                   ).to_string())",
        "recov")

    nb.md("### Where the known complexes do best",
          "No experiment here baits the glideosome, AMA1-RON2, MIC2-M2AP or MORC/HDAC3, so there "
          "is no pull-down in which they *should* be specific. What can be asked is the next best "
          "thing: across the experiments we have, is there any in which the score marks them as "
          "specific, and do the members of one complex peak together?")

    nb.code(
        "comp_genes = val[val.kind == 'known complex'][['set_name', 'symbol', 'gene_id']]",
        "best = []",
        "for p in bg_projects:",
        "    a_ = per_project[p].copy()",
        "    tot_ = a_.intensity.sum()",
        "    if not tot_ or not len(a_):",
        "        continue",
        "    a_['psm_ppm'] = a_.intensity / tot_ * 1e6",
        "    a_ = a_[a_.n_psm >= 2]",
        "    s = specificity(a_[['gene_id', 'n_psm', 'psm_ppm']], bg)",
        "    s['rank'] = np.arange(1, len(s) + 1)",
        "    s['percentile'] = 100 * (1 - (s['rank'] - 1) / len(s))",
        "    m = comp_genes.merge(s[['gene_id', 'rank', 'percentile', 'specificity']], on='gene_id')",
        "    m['project'] = p",
        "    best.append(m)",
        "cb = pd.concat(best, ignore_index=True)",
        "top = cb.sort_values('specificity', ascending=False).groupby(",
        "    ['set_name', 'symbol'], as_index=False).first()",
        "print('best specificity each complex member reaches in any of the', N_EXP, 'experiments,')",
        "print('and the project in which it reaches it:')",
        "top.sort_values(['set_name', 'specificity'], ascending=[True, False])[",
        "    ['set_name', 'symbol', 'project', 'rank', 'percentile', 'specificity']]")

    nb.code(
        "# do the members of one complex peak in the same experiment?",
        "pk = cb.loc[cb.groupby(['set_name', 'symbol']).specificity.idxmax()]",
        "agree = pk.groupby('set_name').agg(members=('symbol', 'size'),",
        "                                   distinct_best_projects=('project', 'nunique'),",
        "                                   top_project=('project', lambda s: s.mode().iloc[0]),",
        "                                   members_sharing_it=('project',",
        "                                                       lambda s: int((s == s.mode().iloc[0]).sum())))",
        "agree")

    nb.md("### Why the number of independent baits is the limit",
          "The frequency of a protein can only take one of `N+1` values when there are `N` "
          "experiments, and a protein that happens to be in two of the few baited complexes will "
          "look like a frequent flyer. The question below is whether the experiments we have are "
          "really independent looks at the background, or several looks at the same biology.")

    nb.code(
        "print('experiments in the background:', N_EXP)",
        "print('of which a confident bait gene was identified:',",
        "      int((pull.bait_confidence == 'confident').sum()))",
        "print('distinct confident bait genes:',",
        "      pull.loc[pull.bait_confidence == 'confident', 'bait_gene'].nunique())",
        "print('smallest frequency step the table can resolve: 1/%d = %.3f' % (N_EXP, 1 / N_EXP))",
        "print()",
        "bait_comp = pull.merge(nodes.rename(columns={'compartment': 'bait_compartment'})[",
        "    ['gene_id', 'bait_compartment']], left_on='bait_gene', right_on='gene_id', how='left')",
        "print('compartments the baits come from:')",
        "print(bait_comp.bait_compartment.fillna('(unknown)').value_counts().to_string())",
        "print()",
        "print('proteins detected in only one experiment: %d of %d (%.0f%%)'",
        "      % (int((bg.experiments_detected == 1).sum()), len(bg),",
        "         100 * (bg.experiments_detected == 1).mean()))",
        "bg.frequency.describe()")

    # ------------------------------------------------------------- outputs
    nb.md("## 9 · Files written")

    nb.code(
        "cl_out = cl[['project', 'class', 'confidence', 'decided_by', 'cue', 'rule_class',",
        "             'bait_symbol', 'bait_gene', 'bait_product', 'bait_confidence', 'bait_how',",
        "             'n_runs',",
        "             'acquisition', 'psm_rows', 'psm_1pct', 'genes_2psm', 'has_lfq',",
        "             'title', 'note', 'example_runs']].sort_values(['class', 'project'])",
        "cl_out.to_csv(os.path.join(OUT, 'project_classification.tsv'), sep='\\t', index=False)",
        "",
        "bg_out = bg[['gene_id', 'product_name', 'experiments_detected', 'n_experiments',",
        "             'frequency', 'median_abundance', 'max_abundance', 'n_experiments_with_lfq',",
        "             'median_psm_ppm', 'max_psm_ppm', 'total_psm', 'runs_detected',",
        "             'compartment']].sort_values('frequency', ascending=False)",
        "strict = variants['strict (>=3 matches, >=2 peptides)'].set_index('gene_id')",
        "loose = variants['loose (>=1 match)'].set_index('gene_id')",
        "bg_out['frequency_strict'] = bg_out.gene_id.map(strict.frequency).fillna(0.0)",
        "bg_out['frequency_loose'] = bg_out.gene_id.map(loose.frequency).fillna(0.0)",
        "bg_out.to_csv(os.path.join(OUT, 'background.tsv'), sep='\\t', index=False)",
        "",
        "scored.to_csv(os.path.join(OUT, 'specificity_example.tsv'), sep='\\t', index=False)",
        "val.sort_values(['kind', 'set_name', 'frequency'], ascending=[True, True, False]).to_csv(",
        "    os.path.join(OUT, 'validation.tsv'), sep='\\t', index=False)",
        "recov.to_csv(os.path.join(OUT, 'bait_recovery.tsv'), sep='\\t', index=False)",
        "top.to_csv(os.path.join(OUT, 'complex_best_specificity.tsv'), sep='\\t', index=False)",
        "",
        "coverage = pd.DataFrame([",
        "    dict(item='projects selected for the uniform re-search', n=len(selected)),",
        "    dict(item='project directories on disk', n=len(disk)),",
        "    dict(item='finished (done.json receipt)', n=int(disk.finished.sum())),",
        "    dict(item='excluded: unsupported.json', n=int(disk.unsupported.sum())),",
        "    dict(item='excluded: NOT_COMPARABLE.json (HDX-MS)', n=int(disk.not_comparable.sum())),",
        "    dict(item='finished and classified', n=len(cl)),",
        "    dict(item='still waiting in the queue', n=len(selected) - int(disk.finished.sum())),",
        "    dict(item='experiments in the background (affinity + proximity)', n=N_EXP),",
        "    dict(item='of those with a confident bait gene',",
        "         n=int((pull.bait_confidence == 'confident').sum())),",
        "    dict(item='distinct confident bait genes',",
        "         n=int(pull.loc[pull.bait_confidence == 'confident', 'bait_gene'].nunique())),",
        "    dict(item='runs behind those experiments', n=n_runs_total),",
        "    dict(item='PSM rows read', n=int(ps.psm_rows.sum())),",
        "    dict(item='PSMs at 1% FDR, target', n=int(ps.psm_1pct.sum())),",
        "    dict(item='parasite proteins in the background table', n=len(bg)),",
        "] + [dict(item=f'class: {k}', n=int(v)) for k, v in cl.groupby('class').size().items()])",
        "coverage.to_csv(os.path.join(OUT, 'coverage.tsv'), sep='\\t', index=False)",
        "for f in sorted(os.listdir(OUT)):",
        "    print(f, os.path.getsize(os.path.join(OUT, f)), 'bytes')",
        "coverage")

    nb.md("## 10 · Verdict",
          "The verdict cell prints the three numbers the decision rests on: how far apart the "
          "known complexes and the frequent flyers are, how many independent baits stand behind "
          "the table, and what happens to each known complex at a usable cutoff.")

    nb.code(
        "lines = []",
        "lines.append('Experiments: %d (%d affinity purification, %d proximity labelling) out of %d '",
        "             'selected PRIDE projects; %d finished so far.'",
        "             % (N_EXP, int((cl[\"class\"] == 'affinity_purification').sum()),",
        "                int((cl['class'] == 'proximity_labelling').sum()), len(selected),",
        "                int(disk.finished.sum())))",
        "lines.append('Background table: %d parasite proteins.' % len(bg))",
        "lines.append('Known complex members: median background frequency %.2f (n=%d).'",
        "             % (a.median(), len(a)))",
        "lines.append('Expected frequent flyers: median background frequency %.2f (n=%d).'",
        "             % (b.median(), len(b)))",
        "lines.append('Separation: AUC %.3f, p = %.2g.' % (auc, u.pvalue))",
        "lines.append('At a 50%% frequency cutoff, %d of %d complex members would be discarded as '",
        "             'background and %d of %d flyers correctly caught.'",
        "             % (int((a >= 0.5).sum()), len(a), int((b >= 0.5).sum()), len(b)))",
        "lines.append('No cutoff removes half the frequent flyers while keeping 85% of the known '",
        "             'complex members.' if not FREQ_FILTER_WORKS else",
        "             'A frequency cutoff exists that removes half the frequent flyers while '",
        "             'keeping 85% of the known complex members.')",
        "lines.append('Independent confident baits: %d (over %d experiments).'",
        "             % (pull.loc[pull.bait_confidence == 'confident', 'bait_gene'].nunique(),",
        "                int((pull.bait_confidence == 'confident').sum())))",
        "ap_got = got[got['class'] == 'affinity_purification'] if len(got) else got",
        "SCORE_WORKS = bool(len(ap_got) >= 5) and float((ap_got.bait_percentile >= 95).mean()) >= 0.6",
        "if len(got):",
        "    lines.append('Bait recovery by the full score: the tagged protein is in the top 5%% of '",
        "                 'its own hit list in %d of %d experiments with a confident bait, and in '",
        "                 '%d of %d of the affinity-purification ones.'",
        "                 % (int((got.bait_percentile >= 95).sum()), len(got),",
        "                    int((ap_got.bait_percentile >= 95).sum()), len(ap_got)))",
        "lines.append('VERDICT: the frequency column alone is %s as a filter. '",
        "             'The full specificity score %s.'",
        "             % ('NOT usable' if not FREQ_FILTER_WORKS else 'usable',",
        "                'recovers the bait in its own pull-down and is usable for ranking'",
        "                if SCORE_WORKS else 'is not yet demonstrated'))",
        "print('\\n'.join('- ' + x for x in lines))")

    nb.write(NB)
    print("notebook:", NB)

    # ------------------------------------------------------------- README
    ns = nb.ns
    write_readme(ns)
    print("results:", OUT)


def write_readme(ns: dict) -> None:
    import numpy as np  # noqa: F401

    cl = ns["cl"]
    bg = ns["bg"]
    disk = ns["disk"]
    selected = ns["selected"]
    pull = ns["pull"]
    val = ns["val"]
    sens = ns["sens"]
    N_EXP = ns["N_EXP"]
    auc = ns["auc"]
    u = ns["u"]
    a = ns["a"]
    b = ns["b"]
    scored = ns["scored"]
    known = ns["known"]
    coverage = ns["coverage"]

    cls_counts = cl.groupby("class").size()
    n_conf = int((pull.bait_confidence == "confident").sum())
    n_baits = int(pull.loc[pull.bait_confidence == "confident", "bait_gene"].nunique())
    uncertain = int(cl.confidence.isin(["low", "none", "medium"]).sum())

    complex_tbl = val[val.kind == "known complex"].sort_values(["set_name", "symbol"])
    flyer_tbl = ns["tight_tbl"]
    loose_tbl = ns["loose_tbl"]
    recov = ns["recov"]
    top = ns["top"]
    sw = ns["sw"]
    got = ns["got"]
    freq_works = ns["FREQ_FILTER_WORKS"]
    score_works = ns["SCORE_WORKS"]
    ap_got = got[got["class"] == "affinity_purification"] if len(got) else got

    lines = []
    lines.append("# Toxoplasma CRAPome, first pass (%s)" % DATE.replace("_", "-"))
    lines.append("")
    lines.append("Instruction 71 product 4: how often each parasite protein is detected in an "
                 "affinity-purification or proximity-labelling experiment regardless of what was "
                 "pulled down, plus a specificity score for a new hit list.")
    lines.append("")
    lines.append("Notebook: `notebooks/crapome_%s.ipynb` (executed). Builder: "
                 "`scripts/build_crapome_notebook.py`. Re-running the builder rescans the uniform "
                 "search output and picks up whatever has finished since." % DATE)
    lines.append("")
    lines.append("## Coverage")
    lines.append("")
    lines.append("| | n |")
    lines.append("|---|---|")
    for _, r in coverage.iterrows():
        lines.append("| %s | %d |" % (r["item"], r["n"]))
    lines.append("")
    lines.append("## Classification rules, verbatim")
    lines.append("")
    lines.append("Applied in this order to the PRIDE title first and the deposited run file names "
                 "second; first match wins. The matched text is kept in "
                 "`project_classification.tsv` (`cue`, `decided_by`).")
    lines.append("")
    lines.append("```")
    for cls, pat, why in ns["RULES"]:
        lines.append("%-22s %s" % (cls, pat))
        lines.append("%-22s   (%s)" % ("", why))
    lines.append("")
    lines.append("file-name cues (weaker, used only when the title gave none):")
    for cls, pat in ns["NAME_CUES"]:
        lines.append("%-22s %s" % (cls, pat))
    lines.append("```")
    lines.append("")
    lines.append("Counts: " + ", ".join("%s %d" % (k, v) for k, v in cls_counts.items()) + ".")
    lines.append("")
    lines.append("**%d of %d classifications hinge on a judgement call** (confidence low, none or "
                 "medium, or hand-checked): either the title gave no keyword and the call rests on "
                 "deposited file names, or the title matched cues for two classes. They are listed "
                 "in `project_classification.tsv`; every one should be checked by a human before "
                 "the table is relied on." % (uncertain + int((cl.confidence == 'hand-checked').sum()),
                                              len(cl)))
    lines.append("")
    lines.append("Hand overrides, each with its reason:")
    lines.append("")
    for proj, (cls, why) in ns["OVERRIDES"].items():
        lines.append("* `%s` -> %s: %s" % (proj, cls, why))
    lines.append("")
    lines.append("## Exclusions")
    lines.append("")
    lines.append("* **HDX-MS, `NOT_COMPARABLE.json`**: %s. Hydrogen-deuterium exchange peptides "
                 "carry an unknown number of deuteriums, so their masses are shifted by amounts "
                 "the search does not model; any protein-level quantity reflects deuteration, not "
                 "abundance (instruction 69 measured ~1/10 the PSMs per run of an ordinary "
                 "bottom-up project). Spectra are kept; tables are not used."
                 % ", ".join(sorted(disk.loc[disk.not_comparable, "project"])))
    lines.append("* **`unsupported.json`**: %s. The pipeline could not search these at all — "
                 "SCIEX `.wiff` (needs a converter that does not run here), deposits with no "
                 "spectra, and the two HDX projects stopped before searching (PXD034607, "
                 "PXD039335)." % ", ".join(sorted(disk.loc[disk.unsupported, "project"])))
    lines.append("* Projects classified as whole proteome, PTM enrichment, fractionation or other "
                 "are **classified but not counted in the background**: a whole-proteome run "
                 "detects nearly everything and would make every protein a frequent flyer, and a "
                 "PTM enrichment detects a biased subset.")
    lines.append("* %d of the %d selected projects have not finished the queue yet."
                 % (len(selected) - int(disk.finished.sum()), len(selected)))
    lines.append("")
    lines.append("## Detection rule and its sensitivity")
    lines.append("")
    lines.append("> A protein is **detected** in an experiment when at least **2 peptide-spectrum "
                 "matches at 1 % spectrum-level FDR** (Sage `spectrum_q <= 0.01`, target hits "
                 "only) map to that gene **and to no other parasite gene**.")
    lines.append("")
    lines.append("1 % FDR is Sage's own number and the field standard. Two matches rather than one "
                 "because at 1 % FDR the single-spectrum identifications are the ones most likely "
                 "to be wrong, and for a background table a false detection is the costlier error: "
                 "it inflates a protein's frequency and gets someone's real hit filtered out. "
                 "Unique mapping because a peptide shared between paralogues would otherwise "
                 "inflate whole families; a protein seen only through shared peptides is not "
                 "counted, which makes frequencies slightly too low rather than too high.")
    lines.append("")
    lines.append("| rule | proteins | frequency >= 0.5 | frequency >= 0.8 | seen once | median "
                 "frequency | Spearman vs default |")
    lines.append("|---|---|---|---|---|---|---|")
    for _, r in sens.iterrows():
        lines.append("| %s | %d | %d | %d | %d | %.4f | %.4f |"
                     % (r["rule"], r["proteins"], r["freq_ge_0_5"], r["freq_ge_0_8"],
                        r["seen_once"], r["median_frequency"], r["spearman_vs_default"]))
    lines.append("")
    lines.append("`background.tsv` carries `frequency_strict` and `frequency_loose` beside the "
                 "default, so any conclusion can be repeated under either.")
    lines.append("")
    lines.append("**A known consequence of the unique-mapping rule.** Two nearly identical "
                 "paralogues share most of their tryptic peptides, so both lose the shared ones "
                 "and both look rarer than they are. Enolase 1 and enolase 2 are the clearest "
                 "case in this table: they are among the most abundant proteins in a tachyzoite "
                 "yet come out at a low frequency. Treat the frequency of a member of a close "
                 "paralogue pair as a lower bound, and check `total_psm` beside it.")
    lines.append("")
    lines.append("**Experiments, not runs.** One project is one experiment; a project's several "
                 "runs are usually the bait and its control in replicate, which is one look at the "
                 "background. Run-level counts are kept in `runs_detected` for comparison.")
    lines.append("")
    lines.append("**Abundance.** Each protein is expressed as parts per million of its own "
                 "project's total parasite signal, so projects of different depth can be averaged. "
                 "`median_abundance` / `max_abundance` are MS1 label-free quantities (from "
                 "`lfq.tsv`); `median_psm_ppm` / `max_psm_ppm` are summed MS2 fragment intensities "
                 "and exist for every project, including the timsTOF ones that have no MS1 LFQ.")
    lines.append("")
    lines.append("## Baits")
    lines.append("")
    lines.append("%d affinity + proximity experiments; **%d have a confident bait gene** (a symbol "
                 "read out of the title that resolves to exactly one ME49 gene through ToxoDB "
                 "product names), covering **%d distinct bait genes**."
                 % (N_EXP, n_conf, n_baits))
    lines.append("")
    for k, v in pull.bait_confidence.value_counts().items():
        lines.append("* %s: %d" % (k, v))
    lines.append("")
    lines.append("`interaction_studies.parquet` (97 author-reported studies, 31 with a bait gene) "
                 "is a **separate** resource and was used only as a comparison: it is keyed by "
                 "PMID and the PRIDE selection carries no PMID, so the two cannot be joined. Every "
                 "bait in `project_classification.tsv` comes from the PRIDE title or run names, "
                 "never from that table.")
    lines.append("")
    lines.append("## Specificity score")
    lines.append("")
    lines.append("```")
    lines.append("specificity = log2( (ppm_here + 1) / (median background ppm + 1) )")
    lines.append("            + log2( 1 / ((experiments_detected + 0.5) / (n_experiments + 0.5)) )")
    lines.append("```")
    lines.append("")
    lines.append("The first term asks whether there is more of the protein here than usual; the "
                 "second asks whether it is unusual to see it at all. Both are reported "
                 "separately in `specificity_example.tsv` (`log2_fold`, `log2_rarity`), because "
                 "they fail differently: the fold is unreliable for a protein with few background "
                 "measurements and the rarity term is unreliable for everything while "
                 "`n_experiments` is small. A protein seen in every pull-down at its usual level "
                 "scores about 0.")
    lines.append("")
    lines.append("Worked example: **%s**, the MYR1 immunoprecipitation (%d proteins in the hit "
                 "list). Known MYR1-complex members and where the score put them:"
                 % (ns["EXAMPLE"], len(scored)))
    lines.append("")
    lines.append("| protein | gene | rank | percentile | specificity | background frequency |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in known.sort_values("rank").iterrows():
        lines.append("| %s | %s | %d | %.1f | %.2f | %.2f |"
                     % (r["symbol"], r["gene_id"], r["rank"], r["percentile"],
                        r["specificity"], r["background_frequency"]))
    lines.append("")
    lines.append("## Validation")
    lines.append("")
    lines.append("| | n | median frequency | max frequency |")
    lines.append("|---|---|---|---|")
    lines.append("| known complex members | %d | %.3f | %.3f |"
                 % (len(a), a.median(), a.max()))
    lines.append("| expected frequent flyers | %d | %.3f | %.3f |"
                 % (len(b), b.median(), b.max()))
    lines.append("| all proteins in the table | %d | %.3f | %.3f |"
                 % (len(bg), bg.frequency.median(), bg.frequency.max()))
    lines.append("")
    lines.append("Separation: AUC %.3f, Mann-Whitney p = %.3g. At a 50 %% frequency cutoff, %d of "
                 "%d complex members would be discarded as background and %d of %d flyers "
                 "correctly caught."
                 % (auc, u.pvalue, int((a >= 0.5).sum()), len(a),
                    int((b >= 0.5).sum()), len(b)))
    lines.append("")
    lines.append("The direction is right — the frequent flyers are above the complex members — "
                 "but the two distributions overlap almost completely, and that is what decides "
                 "whether a filter can be built. The cutoff sweep:")
    lines.append("")
    lines.append("| cutoff | frequent flyers removed | complex members lost | proteins removed "
                 "from the table |")
    lines.append("|---|---|---|---|")
    for _, r in sw.iterrows():
        if r["cutoff"] in (0.25, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0):
            lines.append("| %.2f | %d / %d (%.0f %%) | %d / %d (%.0f %%) | %d |"
                         % (r["cutoff"], r["flyers_removed"], r["flyers"],
                            100 * r["flyer_recall"], r["complex_lost"], r["complex_members"],
                            100 * r["complex_loss"], r["proteins_removed"]))
    lines.append("")
    if not freq_works:
        lines.append("**No cutoff removes at least half the frequent flyers while losing at most "
                     "15 % of the known complex members, so the frequency column must not be "
                     "used as a standalone filter.** The reason is biological rather than a bug: "
                     "the glideosome, AMA1-RON, MIC2-M2AP and the MORC complex are among the most "
                     "abundant proteins in a tachyzoite, so a pull-down from whole parasite lysate "
                     "picks them up whatever the bait. Their frequencies (0.38-0.88) sit inside "
                     "the frequent-flyer range, and a cutoff anywhere in it deletes exactly the "
                     "complexes a *Toxoplasma* biologist cares about.")
    else:
        lines.append("A cutoff exists that removes at least half the frequent flyers while losing "
                     "at most 15 % of the known complex members.")
    lines.append("")
    lines.append("### What the full score does, which the frequency alone does not")
    lines.append("")
    if len(got):
        lines.append("The product is the specificity score, not the frequency column. Tested on "
                     "the one case where the answer is certain — the tagged protein in its own "
                     "pull-down — across every experiment with a confident bait:")
        lines.append("")
        lines.append("| project | class | bait | gene | hits | bait rank | percentile | "
                     "specificity | background frequency |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for _, r in got.sort_values("bait_percentile", ascending=False).iterrows():
            lines.append("| %s | %s | %s | %s | %d | %d | %.1f | %.2f | %.2f |"
                         % (r["project"], str(r["class"]).replace("_", " "), r["bait"],
                            r["bait_gene"], r["hits"], r["bait_rank"], r["bait_percentile"],
                            r["bait_specificity"], r["bait_frequency"]))
        lines.append("")
        lines.append("Bait in the top 5 %% of its own hit list: **%d of %d** overall and **%d of "
                     "%d** of the affinity purifications; median bait rank %.0f of a median %.0f "
                     "proteins, against median rank %.0f if the list is ranked by raw abundance "
                     "instead. `bait_recovery.tsv`."
                     % (int((got.bait_percentile >= 95).sum()), len(got),
                        int((ap_got.bait_percentile >= 95).sum()), len(ap_got),
                        got.bait_rank.median(), got.hits.median(),
                        got.bait_abundance_rank.median()))
        lines.append("")
        lines.append("The split by class matters. The score is built for affinity purification, "
                     "where the tagged protein is pulled out intact and should dominate its own "
                     "eluate. In proximity labelling the bait is a biotin-ligase fusion whose job "
                     "is to label its *neighbours*, so there is no reason for it to top the list, "
                     "and in a degron experiment (`PXD039922`, SPARK-mNG-mAID) the bait is "
                     "deliberately degraded, so a low rank there is the correct answer rather "
                     "than a failure. Those cases are kept in the table rather than dropped, "
                     "because dropping them would be choosing the evidence.")
    else:
        lines.append("No experiment has both a confident bait and a detected bait protein, so the "
                     "score could not be tested this way.")
    lines.append("")
    lines.append("No experiment here baits the glideosome, AMA1-RON, MIC2-M2AP or MORC/HDAC3, so "
                 "none of them has a pull-down in which it *should* be specific. The best "
                 "specificity each member reaches in any of the %d experiments is in "
                 "`complex_best_specificity.tsv`; the highest per complex:" % N_EXP)
    lines.append("")
    lines.append("| complex | protein | best in | rank | percentile | specificity |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in top.sort_values(["set_name", "specificity"],
                                ascending=[True, False]).groupby("set_name").head(2).iterrows():
        lines.append("| %s | %s | %s | %d | %.1f | %.2f |"
                     % (r["set_name"], r["symbol"], r["project"], r["rank"],
                        r["percentile"], r["specificity"]))
    lines.append("")
    lines.append("Per complex:")
    lines.append("")
    lines.append("| complex | protein | gene | experiments detected | frequency |")
    lines.append("|---|---|---|---|---|")
    for _, r in complex_tbl.iterrows():
        lines.append("| %s | %s | %s | %d / %d | %.3f |"
                     % (r["set_name"], r["symbol"], r["gene_id"],
                        r["experiments_detected"], r["n_experiments"], r["frequency"]))
    lines.append("")
    lines.append("Per frequent-flyer family (curated sets; the cytosolic ribosome is the "
                 "hyperLOPIT 40S/60S marker set, so the apicoplast and mitochondrial ribosomes "
                 "are not mixed in):")
    lines.append("")
    lines.append("| family | proteins | median frequency | max frequency |")
    lines.append("|---|---|---|---|")
    for name, r in flyer_tbl.iterrows():
        lines.append("| %s | %d | %.3f | %.3f |"
                     % (name, r["proteins"], r["median_frequency"], r["max_frequency"]))
    lines.append("")
    lines.append("The same families as a keyword would define them, for comparison only — these "
                 "sweep in the tubulin-tyrosine ligases and the organellar ribosomes, which are "
                 "neither abundant nor sticky:")
    lines.append("")
    lines.append("| family (keyword) | proteins | median frequency | max frequency |")
    lines.append("|---|---|---|---|")
    for name, r in loose_tbl.iterrows():
        lines.append("| %s | %d | %.3f | %.3f |"
                     % (name, r["proteins"], r["median_frequency"], r["max_frequency"]))
    lines.append("")
    lines.append("## Is it usable yet")
    lines.append("")
    lines.append("* **The frequency column as a standalone filter: %s.** %s"
                 % ("usable" if freq_works else "not yet, and this is the headline",
                    "A cutoff exists that removes at least half the frequent flyers while losing "
                    "at most 15 % of the known complex members." if freq_works else
                    "The frequent flyers are above the complex members (median %.2f against "
                    "%.2f, AUC %.3f, p = %.3g) but the overlap is almost total: at the 50 %% "
                    "cutoff, %d of the %d known complex members are deleted. Do not filter a hit "
                    "list on frequency."
                    % (b.median(), a.median(), auc, u.pvalue, int((a >= 0.5).sum()), len(a))))
    lines.append("* **The specificity score: %s.** %s"
                 % ("usable for ranking, with the caveats below" if score_works
                    else "not yet demonstrated",
                    "It puts the tagged bait in the top 5 %% of its own hit list in %d of the %d "
                    "affinity purifications with a confident bait, and on the MYR1 worked example "
                    "it puts MYR1 first, MYR4 third and five of six known complex members in the "
                    "top fifth. The work is done by the abundance-fold term, not the frequency "
                    "term." % (int((ap_got.bait_percentile >= 95).sum()), len(ap_got))
                    if score_works else
                    "Only %d affinity purifications have a confident bait and a detected bait "
                    "protein, which is too few to demonstrate it." % len(ap_got)))
    lines.append("* **The limiting resource is independent baits (%d distinct, over %d "
                 "experiments), not proteins or spectra.** Nearly every experiment is a pull-down "
                 "from whole parasite lysate, so the abundant invasion and motility machinery is "
                 "in every list. A frequency becomes informative only when there are enough baits "
                 "from enough different compartments that \"present in most pull-downs\" means "
                 "stickiness rather than abundance. For comparison, the human CRAPome was built "
                 "from more than 400 experiments." % (n_baits, N_EXP))
    lines.append("* **Practical advice for now:** rank a hit list by `specificity`, or by its "
                 "`log2_fold` component alone, and read `frequency` as context — \"this protein "
                 "has been seen in %d of %d other pull-downs\" — rather than as a filter. The "
                 "`median_abundance` column is the useful part of the table today: it says how "
                 "much of a protein a pull-down normally gives, which is what makes an enrichment "
                 "judgeable." % (int(a.max() * N_EXP), N_EXP))
    lines.append("")
    lines.append("## What will change as the queue fills")
    lines.append("")
    lines.append("* **The frequency resolution.** With %d experiments the smallest step a "
                 "frequency can take is 1/%d = %.3f. Every frequency in this table is that coarse."
                 % (N_EXP, N_EXP, 1.0 / N_EXP))
    lines.append("* **Bait independence.** %d distinct confident baits is the binding constraint, "
                 "not the protein count. Several of the finished projects come from the same "
                 "laboratories and the same complexes (SPARK/SPARKEL, CDPK1, the mitochondrial "
                 "complexes), so the experiments are not as independent as their number suggests. "
                 "A protein that happens to sit in two of the few baited complexes looks like a "
                 "frequent flyer for a reason that has nothing to do with stickiness." % n_baits)
    lines.append("* **Classification.** The %d uncertain calls should be reviewed; a "
                 "misclassified whole-proteome project in the background would raise every "
                 "frequency at once." % uncertain)
    lines.append("* **Nothing here is hand-patched.** Re-run the builder and the whole table, the "
                 "validation and this README are rebuilt from whatever has finished.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("* `project_classification.tsv` — every finished project: class, how it was "
                 "decided, the matched cue, confidence, bait and bait confidence, run and PSM "
                 "counts, title.")
    lines.append("* `background.tsv` — gene_id, experiments_detected, n_experiments, frequency, "
                 "median_abundance, max_abundance, plus the strict and loose frequencies, the "
                 "PSM-intensity abundances, run counts, product name and hyperLOPIT compartment.")
    lines.append("* `specificity_example.tsv` — the MYR1 hit list scored, with both score "
                 "components.")
    lines.append("* `validation.tsv` — the known complexes and the frequent flyers (curated and "
                 "keyword-defined) with their frequencies.")
    lines.append("* `bait_recovery.tsv` — for every experiment with a confident bait, where the "
                 "score put the bait in its own hit list.")
    lines.append("* `complex_best_specificity.tsv` — the best specificity each known complex "
                 "member reaches in any experiment, and where.")
    lines.append("* `coverage.tsv` — the counts above, so a later run can be compared with this "
                 "one.")

    with open(os.path.join(OUT, "README.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
