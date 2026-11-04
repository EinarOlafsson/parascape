"""Build the review page for the host-facing truth candidates (instruction 72).

Einar reviews 140 literature candidates before they become the truth set that host-facing predictions are
scored against. The page shows one piece of evidence at a time with its PubMed reference and a diagram of the
compartment claimed, and writes each decision to the artifact's shared store so it can be read back here.

    python scripts/build_host_facing_review_page.py            # writes the page into the scratchpad
"""
import json
import os
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / 'results/host_facing_truth_candidates_2026_11_04/candidates.tsv'
OUT = Path(os.environ.get('PARASCAPE_PAGE_OUT',
                          '/tmp/claude-1000/-home-carruthers-Documents-repo/'
                          '12a1088e-6d4a-4c6e-87e9-dc1b5e58c660/scratchpad')) / 'host_facing_review.html'
#: Where each claimed compartment sits, outward from the parasite; drives the diagram and the sort order.
PLACES = ['PV_lumen_host_interacting', 'moving_junction', 'PVM', 'host_cytosol', 'host_organelle',
          'host_plasma_membrane', 'host_nucleus']
LABELS = {'PV_lumen_host_interacting': 'vacuole lumen, binds host factors', 'moving_junction': 'moving junction',
          'PVM': 'vacuole membrane', 'host_cytosol': 'host cytosol', 'host_organelle': 'host organelle',
          'host_plasma_membrane': 'host plasma membrane', 'host_nucleus': 'host nucleus'}


def rows() -> list:
    frame = pd.read_csv(CANDIDATES, sep='\t')
    frame['gene_id_me49'] = frame.gene_id_me49.fillna('')     # the id_unresolved rows carry no gene id
    counts = frame.gene_id_me49.value_counts()
    out = []
    for row in frame.itertuples():
        pmid = '' if pd.isna(row.pmid) else str(int(row.pmid))
        key = re.sub(r'[^A-Za-z0-9]+', '_', f'{row.gene_id_me49}_{pmid}_{row.host_compartment}').strip('_')
        out.append(dict(key=key, gene=str(row.gene_id_me49), name=str(row.protein_name),
                        aliases='' if pd.isna(row.aliases) else str(row.aliases),
                        place=str(row.host_compartment), strength=str(row.strength),
                        method='' if pd.isna(row.method) else str(row.method),
                        pmid=pmid, year='' if pd.isna(row.year) else str(int(row.year)),
                        author='' if pd.isna(row.first_author) else str(row.first_author),
                        evidence='' if pd.isna(row.evidence_text) else str(row.evidence_text),
                        source=str(row.source), go='' if pd.isna(row.go_term) else str(row.go_term),
                        status=str(row.status), siblings=int(counts[row.gene_id_me49])))
    order = {p: i for i, p in enumerate(PLACES)}
    out.sort(key=lambda r: (r['status'] != 'ok', order.get(r['place'], 9), r['strength'] != 'direct', r['name']))
    return out


PAGE = """<title>Host-Facing Evidence Review</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,600;1,6..72,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: one piece of evidence at a time, centred under a sticky progress bar; a compartment diagram places
   the claim between parasite and host nucleus. */
:root {
  --ink: #16202b; --ink-soft: #55657a; --ink-faint: #8798ab;
  --paper: #eef1f4; --card: #ffffff; --line: #d3dae2;
  --yes: #0f7a6c; --no: #a8432f; --maybe: #8a6212; --flag: #7b3f9d;
  --display: "Newsreader", Georgia, serif;
  --body: "IBM Plex Sans", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ink: #e6ecf3; --ink-soft: #a3b2c4; --ink-faint: #74869a;
    --paper: #0f161e; --card: #18222e; --line: #2c3a4a;
    --yes: #4fd1b5; --no: #f0836a; --maybe: #e0b34a; --flag: #c59be0;
    color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --ink: #e6ecf3; --ink-soft: #a3b2c4; --ink-faint: #74869a;
  --paper: #0f161e; --card: #18222e; --line: #2c3a4a;
  --yes: #4fd1b5; --no: #f0836a; --maybe: #e0b34a; --flag: #c59be0;
  color-scheme: dark;
}
* { box-sizing: border-box; }
body { background: var(--paper); color: var(--ink); font-family: var(--body); margin: 0; line-height: 1.5; }
.bar { position: sticky; top: env(safe-area-inset-top, 0px); z-index: 5; background: var(--paper);
       border-bottom: 1px solid var(--line); }
.bar-in { max-width: 760px; margin: 0 auto; padding: 10px 16px; display: flex; gap: 12px; align-items: center;
          flex-wrap: wrap; }
h1 { font-family: var(--display); font-size: 1.05rem; font-weight: 600; margin: 0; letter-spacing: -0.01em; }
.count { font-family: var(--mono); font-size: 0.8rem; color: var(--ink-soft); font-variant-numeric: tabular-nums; }
.track { flex: 1 1 120px; min-width: 0; height: 5px; background: var(--line); border-radius: 3px; overflow: hidden; }
.fill { height: 100%; width: 0; background: var(--yes); transition: width .25s ease; }
main { max-width: 760px; margin: 0 auto; padding-block: 20px 56px; padding-inline: 16px; }
.card { background: var(--card); border: 1px solid var(--line); border-radius: 10px; padding: 20px; }
.head { display: flex; gap: 10px; align-items: baseline; flex-wrap: wrap; }
.name { font-family: var(--display); font-size: 1.7rem; font-weight: 600; margin: 0; letter-spacing: -0.015em; }
.gene { font-family: var(--mono); font-size: 0.82rem; color: var(--ink-soft); }
.tags { display: flex; gap: 6px; flex-wrap: wrap; margin: 12px 0 0; }
.tag { font-size: 0.72rem; letter-spacing: .04em; text-transform: uppercase; padding: 3px 8px; border-radius: 999px;
       border: 1px solid var(--line); color: var(--ink-soft); }
.tag.direct { color: var(--yes); border-color: var(--yes); }
.tag.proximity { color: var(--maybe); border-color: var(--maybe); }
.tag.inferred { color: var(--ink-faint); }
.tag.warn { color: var(--flag); border-color: var(--flag); }
figure { margin: 18px 0 0; }
svg { display: block; width: 100%; height: auto; max-width: 100%; }
.place { font-family: var(--mono); font-size: 0.78rem; color: var(--ink-soft); margin-top: 6px; }
h2 { font-size: 0.72rem; letter-spacing: .07em; text-transform: uppercase; color: var(--ink-faint);
     margin: 20px 0 6px; font-weight: 600; }
.evidence { font-family: var(--display); font-size: 1.08rem; margin: 0; }
.method { color: var(--ink-soft); font-size: 0.9rem; margin: 6px 0 0; }
.ref { margin-top: 6px; font-size: 0.9rem; }
.ref a { color: var(--yes); }
.acts { display: flex; gap: 8px; margin-top: 22px; flex-wrap: wrap; }
button { font: inherit; font-weight: 500; padding: 10px 14px; border-radius: 8px; border: 1px solid var(--line);
         background: var(--card); color: var(--ink); cursor: pointer; }
button:hover { border-color: var(--ink-faint); }
button:focus-visible { outline: 2px solid var(--yes); outline-offset: 2px; }
button.on[data-d="yes"] { background: var(--yes); border-color: var(--yes); color: var(--paper); }
button.on[data-d="no"] { background: var(--no); border-color: var(--no); color: var(--paper); }
button.on[data-d="unsure"] { background: var(--maybe); border-color: var(--maybe); color: var(--paper); }
.key { font-family: var(--mono); font-size: 0.72rem; opacity: .65; }
textarea { width: 100%; margin-top: 12px; min-height: 52px; padding: 9px 11px; border-radius: 8px;
           border: 1px solid var(--line); background: var(--paper); color: var(--ink); font: inherit;
           font-size: 0.9rem; resize: vertical; }
nav { display: flex; gap: 8px; align-items: center; margin-top: 16px; flex-wrap: wrap; }
.spacer { flex: 1 1 0; min-width: 0; }
select { font: inherit; font-size: 0.85rem; padding: 7px 9px; border-radius: 8px; border: 1px solid var(--line);
         background: var(--card); color: var(--ink); }
.note { color: var(--ink-soft); font-size: 0.85rem; margin-top: 18px; }
table { width: 100%; border-collapse: collapse; font-size: 0.85rem; }
th, td { text-align: left; padding: 6px 8px; border-bottom: 1px solid var(--line); }
th { font-size: 0.7rem; letter-spacing: .06em; text-transform: uppercase; color: var(--ink-faint); }
td.g { font-family: var(--mono); font-size: 0.78rem; }
.wrap { overflow-x: auto; }
.tally { display: flex; gap: 16px; flex-wrap: wrap; margin: 0 0 14px; font-family: var(--mono);
         font-size: 0.85rem; font-variant-numeric: tabular-nums; }
.tally b { font-weight: 500; }
.loading { color: var(--ink-faint); font-size: 0.9rem; }
@media (max-width: 460px) { .name { font-size: 1.4rem; } .card { padding: 16px; } }
</style>

<div class="bar"><div class="bar-in">
  <h1>Host-facing evidence</h1>
  <span class="count" id="count">—</span>
  <div class="track"><div class="fill" id="fill"></div></div>
  <select id="filter">
    <option value="undecided">Undecided</option>
    <option value="all">All 140</option>
    <option value="flagged">Needs a check (10)</option>
    <option value="decided">Decided</option>
  </select>
</div></div>

<main>
  <section id="review" class="card">
    <p class="loading" id="boot">Opening the store&hellip;</p>
  </section>

  <nav>
    <button id="prev">&larr; Back <span class="key">k</span></button>
    <button id="next">Skip <span class="key">j</span></button>
    <span class="spacer"></span>
    <button id="toggle">Summary</button>
  </nav>

  <section id="summary" hidden>
    <h2>Decisions so far</h2>
    <div class="tally" id="tally"></div>
    <div class="wrap"><table><thead><tr><th>Protein</th><th>Gene</th><th>Claim</th><th>Evidence</th>
      <th>Decision</th><th>Note</th></tr></thead><tbody id="rows"></tbody></table></div>
  </section>

  <p class="note">Each row is one published claim, so a protein with several papers appears more than once.
  <b>Host-facing</b> means the protein, or a cleaved part of it, reaches beyond the parasite plasma membrane.
  A dense granule or rhoptry protein that stays inside the vacuole is <b>not</b> host-facing. Press
  <span class="key">y</span>, <span class="key">n</span> or <span class="key">u</span>; decisions save as you go.</p>
</main>

<script>
const ROWS = __ROWS__, PLACES = __PLACES__, LABELS = __LABELS__;
const BY_KEY = Object.fromEntries(ROWS.map(r => [r.key, r]));
let decisions = {}, order = [], at = 0, db = null, saving = {};

const el = id => document.getElementById(id);
const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));

function diagram(place) {
  const i = PLACES.indexOf(place), stops = [
    ['parasite', 52], ['vacuole', 132], ['PVM', 190], ['host cytosol', 290], ['host nucleus', 400]];
  const x = {PV_lumen_host_interacting: 132, moving_junction: 190, PVM: 190, host_cytosol: 290,
             host_organelle: 300, host_plasma_membrane: 470, host_nucleus: 400}[place] || 190;
  return `<svg viewBox="0 0 520 96" role="img" aria-label="Claimed location: ${esc(LABELS[place] || place)}">
    <rect x="1" y="18" width="518" height="60" rx="8" fill="none" stroke="var(--line)"/>
    <rect x="455" y="18" width="64" height="60" rx="8" fill="none" stroke="var(--line)"/>
    <ellipse cx="400" cy="48" rx="34" ry="22" fill="none" stroke="var(--line)"/>
    <ellipse cx="160" cy="48" rx="76" ry="26" fill="none" stroke="var(--ink-faint)" stroke-dasharray="4 3"/>
    <ellipse cx="120" cy="48" rx="46" ry="19" fill="none" stroke="var(--ink-soft)"/>
    <circle cx="${x}" cy="48" r="7" fill="var(--yes)"/>
    <circle cx="${x}" cy="48" r="13" fill="none" stroke="var(--yes)" opacity=".5"/>
    ${stops.map(([t, cx]) => `<text x="${cx}" y="90" fill="var(--ink-faint)" font-size="9"
      font-family="var(--mono)" text-anchor="middle">${t}</text>`).join('')}
  </svg>`;
}

function visible() {
  const mode = el('filter').value;
  return ROWS.filter(r => mode === 'all' ? true
    : mode === 'flagged' ? r.status !== 'ok'
    : mode === 'decided' ? decisions[r.key]
    : !decisions[r.key]).map(r => r.key);
}

function render() {
  const total = ROWS.length, done = Object.keys(decisions).filter(k => BY_KEY[k]).length;
  el('count').textContent = `${done} of ${total}`;
  el('fill').style.width = (100 * done / total).toFixed(1) + '%';
  order = visible();
  if (!order.length) {
    el('review').innerHTML = `<p class="evidence">Nothing left under this filter.</p>
      <p class="method">Switch the filter above, or open the summary to check the decisions.</p>`;
    return;
  }
  at = Math.max(0, Math.min(at, order.length - 1));
  const r = BY_KEY[order[at]], d = decisions[r.key] || {};
  el('review').innerHTML = `
    <div class="head"><h2 class="name" style="text-transform:none;letter-spacing:-0.015em;color:var(--ink);font-size:1.7rem;margin:0">${esc(r.name)}</h2>
      <span class="gene">${esc(r.gene)}</span></div>
    <div class="tags">
      <span class="tag ${esc(r.strength)}">${esc(r.strength)} evidence</span>
      <span class="tag">${esc(LABELS[r.place] || r.place)}</span>
      <span class="tag">${esc(r.source)}</span>
      ${r.siblings > 1 ? `<span class="tag">${r.siblings} claims for this gene</span>` : ''}
      ${r.status !== 'ok' ? `<span class="tag warn">${esc(r.status.replace(/_/g, ' '))}</span>` : ''}
      ${r.aliases ? `<span class="tag">${esc(r.aliases)}</span>` : ''}
    </div>
    <figure>${diagram(r.place)}<figcaption class="place">Claimed location: ${esc(LABELS[r.place] || r.place)}</figcaption></figure>
    <h2>What the paper reports</h2>
    <p class="evidence">${esc(r.evidence) || '<span class="method">No quotation recorded.</span>'}</p>
    ${r.method ? `<p class="method">${esc(r.method)}</p>` : ''}
    <p class="ref">${r.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(r.pmid)}/" target="_blank" rel="noopener">${esc(r.author)} ${esc(r.year)} &middot; PMID ${esc(r.pmid)}</a>`
      : '<span class="method">No verified reference.</span>'}${r.go ? ` <span class="gene">${esc(r.go)}</span>` : ''}</p>
    <div class="acts">
      <button data-d="yes" class="${d.decision === 'yes' ? 'on' : ''}">Host-facing <span class="key">y</span></button>
      <button data-d="no" class="${d.decision === 'no' ? 'on' : ''}">Not host-facing <span class="key">n</span></button>
      <button data-d="unsure" class="${d.decision === 'unsure' ? 'on' : ''}">Unsure <span class="key">u</span></button>
    </div>
    <textarea id="note" placeholder="Note (optional) — a correction, a better reference, why it is unclear">${esc(d.note || '')}</textarea>`;
  el('review').querySelectorAll('button[data-d]').forEach(b =>
    b.onclick = () => decide(r.key, b.dataset.d));
  el('note').onchange = e => decide(r.key, (decisions[r.key] || {}).decision || 'unsure', e.target.value);
}

async function decide(key, decision, note) {
  const row = BY_KEY[key];
  const value = {decision, note: note != null ? note : ((decisions[key] || {}).note || ''),
                 gene: row.gene, protein: row.name, place: row.place, pmid: row.pmid,
                 strength: row.strength, at: new Date().toISOString()};
  decisions[key] = value;
  if (note == null) { at += 1; }
  render(); drawSummary();
  if (!db || saving[key]) return;
  saving[key] = true;
  try { await db.collection('decisions').doc(key).set(value); }
  catch (e) { el('count').textContent += ' · not saved'; }
  finally { saving[key] = false; }
}

function drawSummary() {
  const counts = {yes: 0, no: 0, unsure: 0};
  Object.entries(decisions).forEach(([k, v]) => { if (BY_KEY[k] && counts[v.decision] != null) counts[v.decision]++; });
  el('tally').innerHTML = `<span><b>${counts.yes}</b> host-facing</span><span><b>${counts.no}</b> not</span>
    <span><b>${counts.unsure}</b> unsure</span><span><b>${ROWS.length - counts.yes - counts.no - counts.unsure}</b> left</span>`;
  el('rows').innerHTML = ROWS.filter(r => decisions[r.key]).map(r => {
    const d = decisions[r.key];
    return `<tr><td>${esc(r.name)}</td><td class="g">${esc(r.gene)}</td><td>${esc(LABELS[r.place] || r.place)}</td>
      <td>${esc(r.strength)}</td><td>${esc(d.decision)}</td><td>${esc(d.note || '')}</td></tr>`;
  }).join('') || '<tr><td colspan="6" class="method">No decisions yet.</td></tr>';
}

el('filter').onchange = () => { at = 0; render(); };
el('prev').onclick = () => { at -= 1; render(); };
el('next').onclick = () => { at += 1; render(); };
el('toggle').onclick = () => { const s = el('summary'); s.hidden = !s.hidden;
  el('toggle').textContent = s.hidden ? 'Summary' : 'Hide summary'; };
document.addEventListener('keydown', e => {
  if (e.target.tagName === 'TEXTAREA' || e.target.tagName === 'SELECT' || e.metaKey || e.ctrlKey) return;
  const r = BY_KEY[order[at]];
  if (!r) return;
  const map = {y: 'yes', n: 'no', u: 'unsure'};
  if (map[e.key]) { e.preventDefault(); decide(r.key, map[e.key]); }
  else if (e.key === 'j' || e.key === 'ArrowRight') { e.preventDefault(); at += 1; render(); }
  else if (e.key === 'k' || e.key === 'ArrowLeft') { e.preventDefault(); at -= 1; render(); }
});

render(); drawSummary();
let settled = false;
const timer = setTimeout(() => { if (!settled) el('count').textContent += ' · offline'; }, 10000);
claude.use('db').then(store => {
  db = store;
  if (!db) { settled = true; clearTimeout(timer); el('count').textContent += ' · sign in to save'; return; }
  db.collection('decisions').onSnapshot(snap => {
    settled = true; clearTimeout(timer);
    const next = {};
    (snap.docs || []).forEach(doc => { next[doc.id] = doc.data ? doc.data() : doc; });
    Object.entries(next).forEach(([k, v]) => { if (!saving[k]) decisions[k] = v; });
    render(); drawSummary();
  }, () => { settled = true; clearTimeout(timer); });
}).catch(() => { settled = true; clearTimeout(timer); });
</script>
"""


def main():
    data = rows()
    page = (PAGE.replace('__ROWS__', json.dumps(data, separators=(',', ':')))
            .replace('__PLACES__', json.dumps(PLACES))
            .replace('__LABELS__', json.dumps(LABELS)))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(page)
    print(json.dumps(dict(rows=len(data), flagged=sum(r['status'] != 'ok' for r in data),
                          genes=len({r['gene'] for r in data}), bytes=len(page), path=str(OUT))))


if __name__ == '__main__':
    main()
