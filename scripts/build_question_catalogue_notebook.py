"""The question catalogue (instruction 73, phase 3), enumerated by rule, as an executed notebook.

Every column of the Toxoplasma table that can serve as a label becomes one candidate question, worded from its
information category ("Is the gene needed under ...?", "Where is the protein?", "Does it respond to ...?").
Each row carries its labelled-gene count, whether it is already registered (and scored) in the question
registry, its measured redundancy group (|Spearman| >= 0.8 with other columns, `redundancy.groups`), and the
admission tests of instruction 73: answerable (enough labelled genes) and non-redundant (first of its measured
group). The effective number of independent questions is the number of distinct measured groups among the
admitted ones. Nothing is run here: this is the list for Einar to prune (phase 3 gate).

    python scripts/build_question_catalogue_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"question_catalogue_{DATE}")
NB = os.path.join(REPO, "notebooks", f"question_catalogue_{DATE}.ipynb")
MIN_LABELLED = 150


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The question catalogue", namespace={"REPO": REPO, "OUT": OUT, "MIN_LABELLED": MIN_LABELLED})
    nb.md("**What this is (instruction 73, phase 3).** Every column of the Toxoplasma table that can be a label, "
          "turned into one candidate question by its information category, with the admission tests recorded. "
          "A list to prune, not results: no question is scored here.",
          f"**Admission tests recorded per row.** *Answerable*: at least {MIN_LABELLED} labelled genes (the "
          "harness default for a frozen test). *Non-redundant*: the first column of its measured redundancy group "
          "(every pair in a group has |Spearman| >= 0.8); the others are listed as restating it. *Registered*: "
          "already in the question registry (scored through the harness).",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, re, numpy as np, pandas as pd",
            "from parascape import information_categories as IC, question_registry as Q, redundancy as R",
            "nodes = pd.read_parquet(os.path.join(REPO, 'parascape/data/nodes.parquet'))",
            "groups = R.groups('Tg').set_index('column')",
            "registered = {}",
            "for q in Q.QUESTIONS.values():",
            "    col = str(q.label_source).rsplit(':', 1)[-1].split(' ')[0]",
            "    registered.setdefault(col, []).append(q.id)",
            "print(len(nodes.columns), 'columns;', len(Q.QUESTIONS), 'registered questions')")
    nb.md("## 1 · One candidate question per label column",
          "Wording by category. Identity, sequence, study-effort and pure annotation-bookkeeping columns are not "
          "questions about biology and are skipped.")
    nb.code("WORDING = {",
            "  'fitness': 'Is the gene needed for growth under {x}?',",
            "  'phenotype': 'What happens when the gene is perturbed ({x})?',",
            "  'place': 'Where is the protein ({x})?',",
            "  'function': 'What does the gene do ({x})?',",
            "  'transcript': 'Does the gene respond in {x}?',",
            "  'translation': 'Is the gene translated differently in {x}?',",
            "  'protein_abundance': 'How much protein is there in {x}?',",
            "  'modification': 'Is the protein modified ({x})?',",
            "  'interaction': 'What does the protein bind ({x})?',",
            "  'host_response': 'What does the host do in response ({x})?',",
            "}",
            "SKIP = {'identity', 'sequence', 'study_effort', 'conservation', 'structure', 'domain', 'annotation'}",
            "rows = []",
            "for c in nodes.columns:",
            "    try:",
            "        cats = tuple(IC.categories_of(c))",
            "    except Exception:",
            "        cats = ()",
            "    cat = next((k for k in cats if k in WORDING), None)",
            "    if cat is None or set(cats) & SKIP and cat not in WORDING:",
            "        continue",
            "    s = nodes[c]",
            "    labelled = int(s.notna().sum())",
            "    numeric = pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s)",
            "    kind = 'value' if numeric and s.nunique() > 12 else 'class'",
            "    rows.append(dict(column=c, category=cat, kind=kind, labelled=labelled,",
            "                     classes=int(s.nunique()) if kind == 'class' else None,",
            "                     question=WORDING[cat].format(x=c.replace('_', ' ')),",
            "                     registered=', '.join(registered.get(c, [])),",
            "                     group=int(groups.loc[c, 'group']) if c in groups.index else None))",
            "cat_df = pd.DataFrame(rows)",
            "print(len(cat_df), 'candidate questions')",
            "cat_df.category.value_counts().to_frame('questions')")
    nb.md("## 2 · Admission: answerable and non-redundant",
          "Within a measured group, the column with the most labelled genes is kept as the question; the others "
          "are marked as restating it. Columns outside any measured group (categorical, or not measured on enough "
          "genes to correlate) stand on their own.")
    nb.code("cat_df['answerable'] = cat_df.labelled >= MIN_LABELLED",
            "cat_df = cat_df.sort_values(['group', 'labelled'], ascending=[True, False])",
            "first = cat_df.dropna(subset=['group']).groupby('group').column.first()",
            "cat_df['restates'] = [first.get(g) if pd.notna(g) and first.get(g) != c else '' for g, c in zip(cat_df.group, cat_df.column)]",
            "cat_df['non_redundant'] = cat_df.restates == ''",
            "cat_df['admitted'] = cat_df.answerable & cat_df.non_redundant",
            "admitted = cat_df[cat_df.admitted]",
            "independent = admitted.group.nunique() + int(admitted.group.isna().sum())",
            "summary = dict(candidates=len(cat_df), answerable=int(cat_df.answerable.sum()),",
            "               restating_another=int((~cat_df.non_redundant).sum()), admitted=len(admitted),",
            "               effective_independent_questions=int(independent),",
            "               already_registered=int((cat_df.registered != '').sum()))",
            "cat_df.sort_values(['category', 'labelled'], ascending=[True, False]).to_csv(os.path.join(OUT, 'catalogue.tsv'), sep='\\t', index=False)",
            "summary")
    nb.code("admitted.groupby('category').agg(questions=('column', 'size'), median_labelled=('labelled', 'median'),",
            "                                 registered=('registered', lambda s: int((s != '').sum())))")
    nb.md("## 3 · The gaps that must be declared",
          "Questions Einar asked (instruction 73) that no column in the table can label: declared, not fudged.")
    nb.code("gaps = pd.DataFrame([",
            "  dict(question='favourable immune response in humans', why='no human screen exists'),",
            "  dict(question='involvement in autophagy', why='no genome-wide autophagy assay in the table'),",
            "  dict(question='protein scavenging from the host', why='no scavenging screen in the table'),",
            "  dict(question='sporozoite-to-tachyzoite switching', why='no stage-switch screen; sporozoite RNA in one study only')])",
            "gaps.to_csv(os.path.join(OUT, 'declared_gaps.tsv'), sep='\\t', index=False)",
            "gaps")
    nb.code("import json",
            "json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "admitted.sort_values('labelled', ascending=False).head(40)[['question', 'category', 'labelled', 'registered']]")
    nb.write(NB)


if __name__ == "__main__":
    main()
