"""The release evidence matrix (instruction 64.40), as an executed notebook.

One place that says, for the two goals of ParaScape (make what is known about a gene or a list accessible and
synthesised; make inference from it possible and graded), what evidence exists and what does not:

1. every catalogue method on each parasite -- role, held-out per-gene rows, scorecard grade, whether it beats chance
   on its labels, deployment output;
2. every organism space and its state (ready / installed / installable / unavailable; calibrated or not);
3. every scorecard level (gene, class, label, method, meta) and where it is stored;
4. the host pilots (human, mouse) measured on host truth;
5. every unresolved gap: the open progress items and the gaps stated by the finished 64.x items.

    python scripts/build_release_evidence_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"release_evidence_{DATE}")
NB = os.path.join(REPO, "notebooks", f"release_evidence_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Release evidence matrix", namespace={"REPO": REPO, "OUT": OUT})
    nb.md("**What this is (instruction 64.40).** The evidence behind a release, gathered from the stored results "
          "rather than restated: per method, per organism, per scorecard level, and every gap still open. A row "
          "without evidence says so; nothing is filled in by assumption.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import glob, json, os, re, pandas as pd",
            "from parascape import atlas as AT, space_catalogue as SC",
            "R = lambda *p: os.path.join(REPO, *p)")
    nb.md("## 1 · Every method on each parasite")
    nb.code("metrics = pd.read_csv(R('results/label_metrics_2026_11_07/overall.tsv'), sep='\\t')",
            "beats = metrics.assign(beats=metrics.accuracy_all > metrics.chance).groupby(['organism', 'strategy']).agg(",
            "    labels=('target', 'size'), beats_chance_on=('beats', 'sum'), best_macro_f1=('macro_f1', 'max')).reset_index()",
            "frames = []",
            "for org in ('Tg', 'Pf'):",
            "    m = pd.read_csv(R(f'results/inference_atlas_2026_11_07/manifest_{org}.tsv'), sep='\\t')",
            "    frames.append(m.merge(beats[beats.organism == org].drop(columns='organism'), on='strategy', how='left'))",
            "methods = pd.concat(frames, ignore_index=True)",
            "methods['evidence'] = methods.apply(lambda r: 'held-out per-gene rows + card' if r.heldout_genes > 0 else",
            "                                    ('card (own held-out scorecard)' if r.card == 'measured' else 'none'), axis=1)",
            "methods.to_csv(os.path.join(OUT, 'methods.tsv'), sep='\\t', index=False)",
            "methods.groupby(['organism', 'evidence']).size().unstack(fill_value=0)")
    nb.code("methods[methods.heldout_genes > 0][['organism', 'strategy', 'labels', 'beats_chance_on', 'best_macro_f1', 'deployment_claims']].round(3)")
    nb.md("## 2 · Organism spaces")
    nb.code("spaces = pd.DataFrame([dict(code=e.code, species=e.species, kind=e.kind, distribution=e.distribution, state=e.state,",
            "                            calibrated=e.calibrated, note=e.note) for e in SC.catalogue()])",
            "spaces.to_csv(os.path.join(OUT, 'spaces.tsv'), sep='\\t', index=False)",
            "spaces")
    nb.md("## 3 · Scorecard levels and where each is stored")
    nb.code("levels = pd.DataFrame([",
            "  dict(level='gene', what='each gene hidden with its orthogroup: the call, right or wrong', where='parascape/data/track_record.parquet; HarnessResult.rows',",
            "       present=os.path.exists(R('parascape/data/track_record.parquet'))),",
            "  dict(level='class', what='recall / precision / F1 per class, confusion', where='results/label_metrics_2026_11_07/by_class.tsv',",
            "       present=os.path.exists(R('results/label_metrics_2026_11_07/by_class.tsv'))),",
            "  dict(level='label', what='coverage, accuracy against chance, macro F1', where='results/label_metrics_2026_11_07/overall.tsv',",
            "       present=os.path.exists(R('results/label_metrics_2026_11_07/overall.tsv'))),",
            "  dict(level='method', what='calibration scorecard per method and organism', where='parascape/data/strategy_calibration.json',",
            "       present=os.path.exists(R('parascape/data/strategy_calibration.json'))),",
            "  dict(level='meta', what='a model over methods, strictly nested, admitted per label', where='results/meta_inference_nested_2026_11_07*',",
            "       present=bool(glob.glob(R('results/meta_inference_nested_2026_11_07*')))),",
            "  dict(level='question', what='harness frozen tests, touched once and logged', where='results/harness_holdout_log.jsonl',",
            "       present=os.path.exists(R('results/harness_holdout_log.jsonl')))])",
            "levels.to_csv(os.path.join(OUT, 'levels.tsv'), sep='\\t', index=False)",
            "levels")
    nb.md("## 4 · Host pilots, measured on host truth")
    nb.code("hs = pd.read_csv(R('results/human_inference_pilot_2026_11_07/scores.tsv'), sep='\\t').assign(organism='Hs')",
            "mm = pd.read_csv(R('results/mouse_inference_pilot_2026_11_07/scores.tsv'), sep='\\t')",
            "mm = mm[mm.partition.str.startswith('frozen')].assign(organism='Mm', location='BMDM surface', frozen_auc=lambda d: d.roc_auc.round(3),",
            "                                                       frozen_ap=lambda d: d.average_precision.round(3), may_fill=False)",
            "hosts = pd.concat([hs[['organism', 'location', 'frozen_auc', 'frozen_ap', 'may_fill']], mm[['organism', 'location', 'frozen_auc', 'frozen_ap', 'may_fill']]])",
            "hosts.to_csv(os.path.join(OUT, 'host_pilots.tsv'), sep='\\t', index=False)",
            "hosts")
    nb.md("## 5 · Every unresolved gap",
          "Open progress items (from `instructions/INDEX.md` and the open 64.x items), and the gaps the finished 64.x "
          "items stated in their own write-ups.")
    nb.code("index = open(R('instructions/INDEX.md')).read()",
            "open_items = [dict(item=m.group(1), title=re.sub(r'\\s+', ' ', m.group(3))[:120]) for m in",
            "              re.finditer(r'^\\| (\\d+) \\| \\[([^\\]]+)\\]\\(open/[^)]+\\)(.*?)\\|\\s*$', index, flags=re.M)]",
            "sub_open = [dict(item=os.path.basename(p)[:5].replace('_', '.'), title=open(p).readline().strip('# \\n'))",
            "            for p in sorted(glob.glob(R('instructions/open/64_*.md'))) if 'Status: DONE' not in open(p).read()]",
            "stated = []",
            "for p in sorted(glob.glob(R('instructions/done/64_*.md'))) + sorted(glob.glob(R('instructions/open/64_*.md'))):",
            "    text = open(p).read()",
            "    for m in re.finditer(r'Gaps, stated:(.*?)(?:\\n\\n|$)', text, flags=re.S):",
            "        stated.append(dict(item=os.path.basename(p)[:5].replace('_', '.'), gap=re.sub(r'\\s+', ' ', m.group(1)).strip()[:400]))",
            "gaps = pd.DataFrame(open_items + sub_open)",
            "pd.DataFrame(stated).to_csv(os.path.join(OUT, 'stated_gaps.tsv'), sep='\\t', index=False)",
            "gaps.to_csv(os.path.join(OUT, 'open_items.tsv'), sep='\\t', index=False)",
            "print(len(open_items), 'open top-level items;', len(sub_open), 'open 64.x items;', len(stated), 'stated gaps in 64.x write-ups')",
            "pd.DataFrame(stated)")
    nb.write(NB)


if __name__ == "__main__":
    main()
