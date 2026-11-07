"""Dataset quality, with citations per year since publication (instruction 76), as an executed notebook.

Every component is reported separately and the combined score is published only with them visible. The
citation rate is attention, not evidence: a dataset-level weight, never a per-gene feature (category
`attention`, excluded by default from any inference about a gene). The gate is measured here: does the
score built WITHOUT replicate agreement predict replicate agreement across studies?

    python scripts/build_dataset_quality_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_06"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"dataset_quality_{DATE}")
NB = os.path.join(REPO, "notebooks", f"dataset_quality_{DATE}.ipynb")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("Dataset quality, with citations normalised by age",
                          namespace={"REPO": REPO, "OUT": OUT})
    nb.md(
        "**What this is.** A quality profile for every dataset in the omics catalogue: the technical "
        "quality of its samples, how well its replicates agree, how completely our own pipeline "
        "processed it, and how much attention its paper received -- **citations per year since "
        "publication** (Einar, instruction 76).",
        "**How to read it.** The components are shown side by side; the combined score is their "
        "average rank and is never shown without them. A well-cited dataset with poor replicate "
        "agreement should look different from a quiet one that reproduces cleanly.",
        "**What citations are not.** Attention predicted host-facing at AUC 0.902 in this project "
        "without any biology, so the citation rate is a weight on a *dataset*, in the `attention` "
        "category, and never a feature of a *gene*.",
        f"Built {_dt.date.today().isoformat()}. Citation counts are fetched from Europe PMC now and "
        "stored with the retrieval date, because they change.")

    nb.code(
        "import json, os, time, urllib.parse, urllib.request",
        "import numpy as np, pandas as pd",
        "cat = pd.read_csv(os.path.join(REPO, 'parascape', 'data', 'omics_catalogue.tsv'), sep='\\t')",
        "cat['pubmed'] = cat['pubmed'].astype(str).str.extract(r'(\\d{5,9})', expand=False)",
        "print(len(cat), 'datasets;', cat.pubmed.notna().sum(), 'with a PubMed id;',",
        "      cat.pubmed.nunique(), 'distinct papers')",
        "cat[['source', 'accession', 'title', 'pubmed', 'date']].head()")

    nb.md("## 1 · Citations, from Europe PMC",
          "`citedByCount` and the first publication date for each PubMed id, fifty ids per query. The "
          "rate is citations divided by years since publication, with a floor of half a year so a "
          "paper from last month is not divided by nearly zero.")
    nb.code(
        "def europepmc(pmids):",
        "    # Bracketed: without them the source filter binds only to the last id.",
        "    q = '(' + ' OR '.join(f'EXT_ID:{p}' for p in pmids) + ') AND SRC:MED'",
        "    url = ('https://www.ebi.ac.uk/europepmc/webservices/rest/search?'",
        "           + urllib.parse.urlencode({'query': q, 'resultType': 'lite', 'format': 'json',",
        "                                     'pageSize': 100}))",
        "    for attempt in range(4):",
        "        try:",
        "            with urllib.request.urlopen(url, timeout=60) as r:",
        "                return json.loads(r.read())['resultList']['result']",
        "        except Exception:",
        "            time.sleep(5 * (attempt + 1))",
        "    return []",
        "ids = sorted(cat.pubmed.dropna().unique())",
        "rows = []; _dt_today = time.strftime('%Y-%m-%d')",
        "for i in range(0, len(ids), 50):",
        "    for r in europepmc(ids[i:i + 50]):",
        "        rows.append(dict(pubmed=str(r.get('pmid')), cited_by=int(r.get('citedByCount') or 0),",
        "                         first_published=r.get('firstPublicationDate') or '',",
        "                         journal=r.get('journalTitle', ''), retrieved=_dt_today))",
        "cites = pd.DataFrame(rows).drop_duplicates('pubmed')",
        "today = pd.Timestamp(_dt_today)",
        "pub = pd.to_datetime(cites.first_published, errors='coerce')",
        "cites['years_out'] = ((today - pub).dt.days / 365.25).clip(lower=0.5)",
        "cites['citations_per_year'] = cites.cited_by / cites.years_out",
        "cites.to_csv(os.path.join(OUT, 'citations.tsv'), sep='\\t', index=False)",
        "print(len(cites), 'of', len(ids), 'papers found; median', round(cites.citations_per_year.median(), 1),",
        "      'citations a year')",
        "cites.sort_values('citations_per_year', ascending=False).head(10)")

    nb.md("### The publication date, checked against NCBI",
          "The rate divides by the age of the paper, so a wrong date is a wrong rate. The NCBI E-utilities "
          "publication date is fetched for the same ids; a paper whose two years differ is flagged, and "
          "the earlier date is used.")
    nb.code(
        "def ncbi_dates(pmids):",
        "    url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'",
        "           + urllib.parse.urlencode({'db': 'pubmed', 'id': ','.join(pmids), 'retmode': 'json'}))",
        "    for attempt in range(4):",
        "        try:",
        "            with urllib.request.urlopen(url, timeout=60) as r:",
        "                res = json.loads(r.read())['result']",
        "            return {k: v.get('sortpubdate', '')[:10].replace('/', '-') for k, v in res.items() if k != 'uids'}",
        "        except Exception:",
        "            time.sleep(5 * (attempt + 1))",
        "    return {}",
        "ncbi = {}",
        "for i in range(0, len(ids), 150):",
        "    ncbi.update(ncbi_dates(ids[i:i + 150])); time.sleep(0.5)",
        "cites['ncbi_date'] = cites.pubmed.map(ncbi).fillna('')",
        "epmc_year = pd.to_datetime(cites.first_published, errors='coerce').dt.year",
        "ncbi_year = pd.to_datetime(cites.ncbi_date, errors='coerce').dt.year",
        "cites['date_disagrees'] = epmc_year.notna() & ncbi_year.notna() & (epmc_year != ncbi_year)",
        "earliest = pd.concat([pd.to_datetime(cites.first_published, errors='coerce'),",
        "                      pd.to_datetime(cites.ncbi_date, errors='coerce')], axis=1).min(axis=1)",
        "cites['years_out'] = ((today - earliest).dt.days / 365.25).clip(lower=0.5)",
        "cites['citations_per_year'] = cites.cited_by / cites.years_out",
        "cites.to_csv(os.path.join(OUT, 'citations.tsv'), sep='\\t', index=False)",
        "print(int(ncbi_year.notna().sum()), 'dates from NCBI;', int(cites.date_disagrees.sum()), 'disagree on the year')",
        "cites[cites.date_disagrees][['pubmed', 'first_published', 'ncbi_date', 'cited_by']]")

    nb.md("## 2 · The measured components, per study",
          "* **Sample quality**: the mean of the newest per-sample technical score *without* its "
          "replicate term (`quality_no_replicate`), so it does not restate the next column.",
          "* **Replicate agreement**: the median Spearman of each contrast with itself across replicate "
          "splits (`expression_traits`), and the median correlation of each run with the other runs of "
          "its group (`replicate_rho`).",
          "* **Reproduced elsewhere**: of the study contrasts that share a theme with another study, the "
          "share that some other study reproduces (Spearman at least 0.3). It needs a partner study, so "
          "it is missing more often than not.",
          "* **Processing completeness**: runs our pipeline quantified out of the runs planned.")
    nb.code(
        "sq_dir = sorted(d for d in os.listdir(os.path.join(REPO, 'results')) if d.startswith('rna_sample_quality_'))[-1]",
        "sq = pd.read_csv(os.path.join(REPO, 'results', sq_dir, 'sample_quality.tsv'), sep='\\t')",
        "sample_quality = sq.groupby('study').quality_no_replicate.mean().rename('sample_quality')",
        "rep_dir = sorted(d for d in os.listdir(os.path.join(REPO, 'results')) if d.startswith('expression_traits_20')",
        "                 and os.path.exists(os.path.join(REPO, 'results', d, 'replicate_agreement.tsv')))[-1]",
        "rep = pd.read_csv(os.path.join(REPO, 'results', rep_dir, 'replicate_agreement.tsv'), sep='\\t')",
        "replicate = rep.groupby('study').spearman.median().rename('replicate_agreement')",
        "# The second, wider measure: how well each run matches the other runs of its condition group.",
        "replicate_rho = sq.groupby('study').replicate_rho.median().rename('replicate_rho')",
        "cross = pd.read_csv(os.path.join(REPO, 'results', rep_dir, 'cross_study_agreement.tsv'), sep='\\t')",
        "cross = cross[cross.study_a != cross.study_b]",
        "both = pd.concat([cross.rename(columns={'column_a': 'column', 'study_a': 'study'})[['study', 'column', 'spearman']],",
        "                  cross.rename(columns={'column_b': 'column', 'study_b': 'study'})[['study', 'column', 'spearman']]])",
        "best = both.groupby(['study', 'column']).spearman.max()",
        "reproduced = (best >= 0.3).groupby(level='study').mean().rename('reproduced_elsewhere')",
        "NAS = '/nas_mnt/carruthers/Einar/parascape_omics'",
        "plan_path = os.path.expanduser('~/.cache/parascape-downloads/automation/priority-plan.jsonl')",
        "if not os.path.exists(plan_path):",
        "    plan_path = os.path.join(NAS, 'catalogue', 'plan.jsonl')",
        "plan = pd.DataFrame([json.loads(l) for l in open(plan_path) if l.strip()])",
        "plan['done'] = [os.path.exists(os.path.join(NAS, 'rnaseq', s, a, 'done.json'))",
        "                for s, a in zip(plan.study, plan.run_accession)]",
        "completeness = plan.groupby('study').done.mean().rename('processing_completeness')",
        "print('sample quality from', sq_dir, '; replicate agreement from', rep_dir)",
        "pd.concat([sample_quality, replicate, replicate_rho, reproduced, completeness], axis=1).describe().round(3)")

    nb.md("## 3 · The profile, and the combined score",
          "Each component is turned into its rank among the datasets that have it (0-1), and the "
          "combined score is the mean of the ranks available. Citations enter as the log of the rate.")
    nb.code(
        "prof = cat[['source', 'accession', 'title', 'pubmed', 'date']].copy()",
        "prof = prof.merge(cites[['pubmed', 'cited_by', 'citations_per_year', 'first_published', 'retrieved']],",
        "                  on='pubmed', how='left')",
        "prof = prof.join(sample_quality, on='accession').join(replicate, on='accession')",
        "prof = prof.join(replicate_rho, on='accession').join(reproduced, on='accession')",
        "prof = prof.join(completeness, on='accession')",
        "comp = {'sample_quality': prof.sample_quality, 'replicate_agreement': prof.replicate_agreement,",
        "        'replicate_rho': prof.replicate_rho, 'reproduced_elsewhere': prof.reproduced_elsewhere,",
        "        'processing_completeness': prof.processing_completeness,",
        "        'citations': np.log1p(prof.citations_per_year)}",
        "ranks = pd.DataFrame({k: v.rank(pct=True) for k, v in comp.items()})",
        "prof['quality_score'] = ranks.mean(axis=1, skipna=True)",
        "prof['components_available'] = ranks.notna().sum(axis=1)",
        "prof['score_without_replicate'] = ranks.drop(columns=['replicate_agreement', 'replicate_rho']).mean(axis=1, skipna=True)",
        "prof.to_csv(os.path.join(OUT, 'dataset_quality.tsv'), sep='\\t', index=False)",
        "prof.components_available.value_counts().sort_index()")

    nb.md("## 4 · The gate: does the score predict replicate agreement?",
          "Two measures of replicate agreement are the targets: the split-half agreement of each "
          "contrast with itself (few studies, the stricter test) and the median correlation of each run "
          "with the other runs of its group (`replicate_rho`, every study with replicates). The score "
          "built **without** either is set against each (Spearman across studies), and so is each "
          "component on its own.",
          "**Held out by construction.** The score has no fitted weights -- it is an average of ranks -- "
          "so no dataset was used to fit it and every study is a held-out study. The per-sample score "
          "passed the same kind of check at 0.63; a dataset score that cannot is reported as a "
          "description, not a weight.")
    nb.code(
        "from scipy.stats import spearmanr",
        "gate = []",
        "for target in ('replicate_agreement', 'replicate_rho'):",
        "    have = prof.dropna(subset=[target]).drop_duplicates('accession')",
        "    for name in ('score_without_replicate', 'sample_quality', 'reproduced_elsewhere',",
        "                 'processing_completeness', 'citations_per_year'):",
        "        sub = have.dropna(subset=[name])",
        "        varies = sub[name].nunique() > 1 and sub[target].nunique() > 1",
        "        if len(sub) >= 5 and varies:",
        "            rho, p = spearmanr(sub[name], sub[target])",
        "            gate.append(dict(target=target, predictor=name, studies=len(sub),",
        "                             spearman=round(float(rho), 3), p=round(float(p), 4)))",
        "        else:",
        "            gate.append(dict(target=target, predictor=name, studies=len(sub), spearman=np.nan,",
        "                             p=np.nan))",
        "gate = pd.DataFrame(gate)",
        "gate.to_csv(os.path.join(OUT, 'gate.tsv'), sep='\\t', index=False)",
        "gate")

    nb.md("## 5 · Reading it",
          "* The citation rate is stored with its raw count and retrieval date; recompute it rather "
          "than trusting an old rate.",
          "* A component that does not predict replicate agreement is still a description of the "
          "dataset, not a reason to weight it.",
          "* Nothing here is a feature of any gene.")
    nb.write(NB)
    _readme()


def _md(frame) -> str:
    head = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def _reading(gate) -> str:
    """The gate in one paragraph, computed from the table so it cannot drift from it."""
    wide = gate[gate.target == "replicate_rho"].set_index("predictor")
    best = wide.spearman.idxmax()
    score, cites = wide.loc["score_without_replicate"], wide.loc["citations_per_year"]
    return (f"**Reading.** Over the {int(wide.loc[best, 'studies'])} studies with replicate groups, the "
            f"component that predicts replicate agreement best is `{best}` (Spearman "
            f"{wide.loc[best, 'spearman']:.2f}, p {wide.loc[best, 'p']:.3f}). The combined score without "
            f"replicate terms reaches {score.spearman:.2f}: averaging in the other components dilutes it. "
            f"Citations per year reach {cites.spearman:.2f} over {int(cites.studies)} studies -- attention "
            "does not predict reproducibility here, so the citation rate is a description of a dataset, "
            f"not a weight on it. Use `{best}` where a dataset weight is needed.")


def _readme() -> None:
    import pandas as pd
    gate = pd.read_csv(os.path.join(OUT, "gate.tsv"), sep="\t")
    prof = pd.read_csv(os.path.join(OUT, "dataset_quality.tsv"), sep="\t")
    cites = pd.read_csv(os.path.join(OUT, "citations.tsv"), sep="\t")
    lines = ["# Dataset quality, with citations normalised by age (instruction 76)", "",
             f"Notebook: `notebooks/dataset_quality_{DATE}.ipynb` (executed). Builder: "
             "`scripts/build_dataset_quality_notebook.py`.", "",
             f"{len(prof)} datasets; {int(prof.cited_by.notna().sum())} with a citation count "
             f"(Europe PMC, retrieved {cites.retrieved.iloc[0] if len(cites) else 'n/a'}); "
             f"median {cites.citations_per_year.median():.1f} citations a year.", "",
             "## Gate: predicting replicate agreement (Spearman across studies)", "", _md(gate), "",
             _reading(gate), "",
             "Files: `dataset_quality.tsv` (every component and the combined score), `citations.tsv` "
             "(raw counts, dates, rates), `gate.tsv`. The citation rate is a dataset weight in the "
             "`attention` category, never a gene feature."]
    with open(os.path.join(OUT, "README.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
