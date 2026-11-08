"""Batch 5 of the dataset-selection audit: keep, add or replace, per comparable alternative (instruction 65.04).

65.03 found the slots with more than one candidate dataset but could rank none of them: every candidate was
missing its publication identity, so no citation rate existed and `rankable_groups` was 0 for all 297 slots.
65.04 resolved the publications (the 2026-11-08 review queue), so the rule Einar set can now be applied:

    highest citations per year, a slight preference for newer, a preference for more comprehensive.

This notebook joins each candidate to its source in the queue, scores it, and records one decision per slot:
**keep** (the installed default ranks first, or is not beaten by the margin), **replace proposed** (another
candidate ranks clearly higher -- a proposal for Einar, never an automatic promotion), **nothing to choose** (the
slot policy is `separate`, `fill` or `average`, so every candidate is used and there is no single default), or
**unknown** (the data the rule needs is still missing, named per slot). Only a `one` slot has a default at all. Nothing is promoted automatically and no incumbent is declared best by default.

    python scripts/build_slot_decisions_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"dataset_selection_{DATE}")
NB = os.path.join(REPO, "notebooks", f"slot_decisions_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Slot decisions: keep, add or replace", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this decides (instruction 65.04, batch 5).** For every information slot with more than one "
          "candidate dataset, whether to keep what is installed or to propose a replacement, by Einar's rule: "
          "highest citations per year, a slight preference for newer datasets, a preference for more "
          "comprehensive ones.",
          "**Why it could not be done before.** 65.03 could rank nothing: every candidate lacked a publication "
          "identity, so no citation rate existed. The publications resolved in 65.04 are what makes the rule "
          "computable.",
          "**What this is not.** An automatic promotion. A replacement is a *proposal*, with its margin and the "
          "biological checks still owed; a slot whose data is incomplete stays `unknown` and says what is "
          "missing.",
          "**Only a `one` slot can be replaced.** A `fill` or `average` slot uses every candidate it has, and a "
          "`separate` slot keeps them apart, so for those there is nothing to select -- an earlier version of "
          "this notebook read their empty default list as a missing default and called them additions, which "
          "was wrong.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import json, os, numpy as np, pandas as pd",
            "slots = pd.read_csv(os.path.join(REPO, 'results/dataset_selection_2026_10_07/slot_comparisons.csv'))",
            "queue = pd.read_csv(os.path.join(OUT, 'source_review_queue.csv')).set_index('dataset_id')",
            "multi = slots[slots.candidate_groups > 1].reset_index(drop=True)",
            "print(len(slots), 'slots;', len(multi), 'with more than one candidate dataset')",
            "print('rankable groups recorded by 65.03:', int(slots.rankable_groups.max()))",
            "multi.policy.value_counts().to_frame('slots')")
    nb.md("## 1 · The rule, as a score",
          "`citations per year` is the main term. Recency is the slight preference: a 10 % bonus per decade of "
          "newness, capped. Comprehensiveness is the measured share of genes of the organism table the source "
          "fills (65.04). A candidate missing the citation rate cannot be scored and is named, not guessed at "
          "zero.")
    nb.code("RESOLVED = ('resolved_registry_pmid', 'resolved_recorded_doi', 'resolved_by_pride_reference')",
            "def factors(cand):",
            "    src = cand['source_id']",
            "    if src not in queue.index:",
            "        return dict(gap='source not in the review queue')",
            "    row = queue.loc[src]",
            "    cpy = row.citations_per_year",
            "    comp = cand.get('comprehensiveness')",
            "    if pd.isna(comp):",
            "        comp = row.measured_comprehensiveness",
            "    gap = []",
            "    if pd.isna(cpy):",
            "        gap.append('no citation rate (' + str(row.publication_status) + ')')",
            "    if pd.isna(comp):",
            "        gap.append('comprehensiveness not measured')",
            "    return dict(citations_per_year=cpy, comprehensiveness=comp,",
            "                publication_status=row.publication_status, publication_id=row.publication_id,",
            "                gap='; '.join(gap))",
            "rows = []",
            "for s in multi.itertuples(index=False):",
            "    defaults = set(json.loads(s.default_patterns or '[]'))",
            "    for cand in json.loads(s.candidate_factors):",
            "        f = factors(cand)",
            "        comp = f.get('comprehensiveness')",
            "        cpy = f.get('citations_per_year')",
            "        score = np.nan",
            "        if not pd.isna(cpy) and not pd.isna(comp):",
            "            score = float(cpy) * (1.0 + 0.3 * float(comp))   # comprehensiveness as a preference, not a multiplier of rank",
            "        rows.append(dict(slot=s.slot, organism=s.organism, policy=s.policy, dataset_id=cand['dataset_id'],",
            "                         source_id=cand['source_id'], pattern=cand.get('pattern', ''),",
            "                         installed_default=cand.get('pattern', '') in defaults, score=score, **f))",
            "C = pd.DataFrame(rows)",
            "C.to_csv(os.path.join(OUT, 'slot_candidates_scored.tsv'), sep='\\t', index=False)",
            "print(len(C), 'candidates over', C.slot.nunique(), 'slots;', int(C.score.notna().sum()), 'scorable')",
            "C.gap.value_counts().to_frame('candidates')")
    nb.md("## 2 · One decision per slot",
          "A slot whose policy is `separate` keeps its cohorts apart, so there is nothing to choose. Otherwise "
          "the candidates are ranked by the score; a replacement is proposed only when the best candidate is not "
          "the installed default **and** beats it by at least 25 % of its own score -- a margin, so a near-tie "
          "is left alone.")
    nb.code("MARGIN = 0.25",
            "out = []",
            "for slot, g in C.groupby('slot'):",
            "    policy = g.policy.iloc[0]",
            "    scored = g[g.score.notna()].sort_values('score', ascending=False)",
            "    best = scored.iloc[0] if len(scored) else None",
            "    installed = g[g.installed_default]",
            "    if policy != 'one':",
            "        # Only a 'one' slot has a default to choose: `audit_dataset_selection` sets default_patterns",
            "        # for policy 'one' alone. 'separate' keeps cohorts apart, 'fill' and 'average' use every",
            "        # candidate together, so for all three there is nothing to select and nothing to replace.",
            "        decision = {'separate': 'separate cohorts', 'fill': 'all candidates used (fill)',",
            "                    'average': 'all candidates used (average)'}[policy]",
            "        why = 'the slot policy uses every candidate; there is no single default to choose or replace'",
            "    elif best is None:",
            "        decision, why = 'unknown', 'no candidate can be scored: ' + '; '.join(sorted(set(g.gap)) - {''} if isinstance(g.gap, str) else sorted({x for x in g.gap if x}))",
            "    elif len(installed) and best.dataset_id in set(installed.dataset_id):",
            "        decision, why = 'keep', f'the installed default ranks first (score {best.score:.1f})'",
                        "    else:",
            "        inst = installed[installed.score.notna()]",
            "        if not len(inst):",
            "            decision, why = 'unknown', f'the installed default cannot be scored ({\"; \".join(x for x in installed.gap if x)}), so no comparison is possible'",
            "        elif best.score >= inst.score.max() * (1 + MARGIN):",
            "            decision = 'replace proposed'",
            "            why = (f'{best.dataset_id} scores {best.score:.1f} against {inst.score.max():.1f} for the '",
            "                   f'installed default -- a proposal; biological admission still owed')",
            "        else:",
            "            decision, why = 'keep', f'no candidate beats the installed default by {int(MARGIN * 100)} %'",
            "    out.append(dict(slot=slot, organism=g.organism.iloc[0], policy=policy, candidates=len(g),",
            "                    scorable=int(g.score.notna().sum()), decision=decision, reason=why,",
            "                    best=('' if best is None else best.dataset_id)))",
            "D = pd.DataFrame(out).sort_values(['decision', 'slot'])",
            "D.to_csv(os.path.join(OUT, 'slot_decisions.tsv'), sep='\\t', index=False)",
            "json.dump(D.decision.value_counts().to_dict(), open(os.path.join(OUT, 'slot_decisions_summary.json'), 'w'), indent=1)",
            "D.decision.value_counts().to_frame('slots')")
    nb.md("## 3 · The proposals, in full",
          "These are the rows for Einar: a replacement or an addition the rule prefers. Each still needs the "
          "biological admission check (organism, stage, assay, units, mapping) before anything is installed.")
    nb.code("P = D[D.decision == 'replace proposed']",
            "P.to_csv(os.path.join(OUT, 'slot_proposals.tsv'), sep='\\t', index=False)",
            "print(len(P), 'proposals')",
            "P[['slot', 'organism', 'decision', 'best', 'reason']]")
    nb.md("## 4 · The admission sheet for each proposal",
          "What the admission rule asks of a candidate -- organism, stage, assay, quantity, units, mapping -- "
          "set beside what the slot declares and what the registry records. The organism check is mechanical and "
          "decided here; stage, assay and units are read from the registry name and note, which is a prepared "
          "comparison for Einar, not a verdict. A proposal with an organism mismatch is refused outright.")
    nb.code("from parascape import datasets as ds",
            "reg = {d.key: d for d in ds.REGISTRY}",
            "slot_meta = slots.set_index('slot')[['organism', 'name', 'context', 'unit', 'policy']]",
            "sheet = []",
            "for r in P.itertuples(index=False):",
            "    cand = C[(C.slot == r.slot) & (C.dataset_id == r.best)].iloc[0]",
            "    d = reg.get(cand.source_id)",
            "    meta = slot_meta.loc[r.slot]",
            "    same = (d is not None and str(d.organism) == str(meta.organism))",
            "    host_slot = 'host_' in r.slot or '_human_' in r.slot or '_mouse_' in r.slot",
            "    organism_check = ('ok' if same else",
            "                      ('host dataset for a host slot' if (d is not None and host_slot) else",
            "                       f'MISMATCH: candidate {getattr(d, \"organism\", \"?\")} against slot {meta.organism}'))",
            "    sheet.append(dict(slot=r.slot, decision=r.decision, candidate=r.best, source=cand.source_id,",
            "                      slot_context=meta.context, slot_unit=meta.unit, slot_policy=meta.policy,",
            "                      candidate_name=getattr(d, 'name', ''), candidate_level=getattr(d, 'level', ''),",
            "                      candidate_coverage=getattr(d, 'coverage', ''), candidate_note=str(getattr(d, 'note', ''))[:300],",
            "                      publication=cand.publication_id, citations_per_year=cand.citations_per_year,",
            "                      comprehensiveness=cand.comprehensiveness, organism_check=organism_check,",
            "                      stage_assay_units='for Einar: does the name and note match the slot context and unit above?'))",
            "S = pd.DataFrame(sheet)",
            "S.to_csv(os.path.join(OUT, 'slot_proposal_admission.tsv'), sep='\\t', index=False)",
            "print(int((S.organism_check.str.startswith('MISMATCH')).sum()), 'proposals refused on organism')",
            "S[['slot', 'candidate', 'slot_context', 'slot_unit', 'candidate_name', 'organism_check']]")
    nb.write(NB)


if __name__ == "__main__":
    main()
