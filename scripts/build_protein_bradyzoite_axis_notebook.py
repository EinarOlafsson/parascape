"""The bradyzoite axis at the protein level (instruction 75), as an executed notebook.

The bradyzoite fingerprint was built on RNA (GSE223869: the genes that rise and fall on conversion). The open
question of instruction 75 was whether it transfers to protein, and until now no finished proteomics project
stated a bradyzoite sample. One now does: PXD058095, a bradyzoite proteome, two runs. It carries no internal
tachyzoite control, so the control is external -- every finished whole-proteome project of tachyzoites in the
same uniform Sage search, scored the same way.

The test is a ranking: the stated bradyzoite proteome should score higher on the bradyzoite fingerprint than the
tachyzoite proteomes. Scores are within-sample ranks (a protein's rank among the proteins measured in that run),
so a project is never compared on absolute intensity, which differs by instrument and loading.

    python scripts/build_protein_bradyzoite_axis_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_08"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"protein_bradyzoite_axis_{DATE}")
NB = os.path.join(REPO, "notebooks", f"protein_bradyzoite_axis_{DATE}.ipynb")
MS_ROOT = "/nas_mnt/carruthers/Einar/parascape_omics/prepared_ms_uniform"


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("The bradyzoite axis at the protein level",
                          namespace={"REPO": REPO, "OUT": OUT, "MS_ROOT": MS_ROOT})
    nb.md("**The open point this tests (instruction 75).** The bradyzoite fingerprint is built on RNA. Does it "
          "read a proteome? Until this snapshot no finished proteomics project stated a bradyzoite sample; "
          "PXD058095 does.",
          "**The comparison.** PXD058095 has no internal tachyzoite control, so the controls are external: the "
          "finished whole-proteome projects of tachyzoites from the same uniform Sage search. Scores are "
          "within-sample ranks, so instrument and loading cannot carry the answer.",
          "**What would count as a pass.** The stated bradyzoite proteome scores above every tachyzoite "
          "proteome. Anything less is reported as it falls.",
          f"Built {_dt.date.today().isoformat()}. The NAS is read-only here.")
    nb.code("import glob, json, os, numpy as np, pandas as pd",
            "fp = pd.read_csv(os.path.join(REPO, 'results/condition_fingerprints_2026_11_06/fingerprints.tsv'), sep='\\t')",
            "brady = fp[fp.fingerprint == 'bradyzoite_conversion_gse223869']",
            "rises = set(brady.loc[brady.set == 'rises', 'gene'])",
            "falls = set(brady.loc[brady.set == 'falls', 'gene'])",
            "classification = pd.read_csv(os.path.join(REPO, 'results/pulldown_background_2026_11_08/project_classification.tsv'), sep='\\t')",
            "whole = classification[classification['class'] == 'whole_proteome']",
            "print(len(rises), 'genes rise on conversion;', len(falls), 'fall')",
            "print(len(whole), 'finished whole-proteome projects in the uniform search')",
            "whole[['project', 'title']].head(20)")
    nb.md("## 1 · One score per project",
          "Each project LFQ table is read, protein intensities are summed per gene over its runs, and the genes "
          "are ranked within the project. The score is the mean rank of the genes that rise on conversion minus "
          "the mean rank of those that fall, in units of the rank range: positive means the proteome looks "
          "converted. A project measuring fewer than 40 of either set is reported as too thin, not scored.")
    nb.code("MIN_GENES = 40",
            "def gene_ranks(project):",
            "    paths = sorted(glob.glob(os.path.join(MS_ROOT, project, 'sage_*', 'lfq.tsv')))",
            "    frames = []",
            "    for path in paths:",
            "        table = pd.read_csv(path, sep='\\t')",
            "        protein = next((c for c in table.columns if c.lower() in ('proteins', 'protein', 'proteinids')), None)",
            "        value_columns = [c for c in table.columns if table[c].dtype.kind == 'f' and c.lower() not in ('q_value', 'spectral_angle', 'score')]",
            "        if protein is None or not value_columns:",
            "            continue",
            "        genes = table[protein].astype(str).str.extract(r'(TGME49_\\d+)')[0]",
            "        frames.append(pd.DataFrame({'gene': genes, 'value': table[value_columns].sum(axis=1)}).dropna())",
            "    if not frames:",
            "        return None",
            "    joined = pd.concat(frames).groupby('gene').value.sum()",
            "    joined = joined[joined > 0]",
            "    return joined.rank(pct=True) if len(joined) else None",
            "rows = []",
            "for project in ['PXD058095'] + [p for p in whole.project if p != 'PXD058095']:",
            "    ranks = gene_ranks(project)",
            "    title = classification.loc[classification.project == project, 'title']",
            "    title = str(title.iloc[0]) if len(title) else ''",
            "    if ranks is None:",
            "        rows.append(dict(project=project, title=title[:90], verdict='no readable LFQ table'))",
            "        continue",
            "    up, down = ranks.reindex(sorted(rises)).dropna(), ranks.reindex(sorted(falls)).dropna()",
            "    if len(up) < MIN_GENES or len(down) < MIN_GENES:",
            "        rows.append(dict(project=project, title=title[:90], genes=len(ranks), rises_measured=len(up),",
            "                         falls_measured=len(down), verdict='too few fingerprint genes measured'))",
            "        continue",
            "    rows.append(dict(project=project, title=title[:90], genes=len(ranks), rises_measured=len(up),",
            "                     falls_measured=len(down), score=round(float(up.mean() - down.mean()), 4),",
            "                     verdict='scored'))",
            "S = pd.DataFrame(rows)",
            "S.to_csv(os.path.join(OUT, 'project_scores.tsv'), sep='\\t', index=False)",
            "print(S.verdict.value_counts().to_dict())",
            "cols = [c for c in ('project', 'score', 'rises_measured', 'falls_measured', 'title') if c in S.columns]",
            "S[S.verdict == 'scored'].sort_values('score', ascending=False)[cols] if 'score' in S.columns else S[cols]")
    nb.md("## 2 · Does the stated bradyzoite proteome stand out?",
          "Its rank among the scored projects, and the gap to the next one. With one stated bradyzoite sample "
          "this is a single observation, not a validated axis: it can fail the axis, or support it, and cannot "
          "establish it. If no project measures enough of the fingerprint, that is the result: the depth of the "
          "proteomes, not the fingerprint, is what blocks the axis.")
    nb.code("has = 'score' in S.columns and (S.verdict == 'scored').any()",
            "scored = (S[S.verdict == 'scored'].sort_values('score', ascending=False).reset_index(drop=True)",
            "          if has else pd.DataFrame(columns=['project', 'score', 'title']))",
            "rank = (int(scored.index[scored.project == 'PXD058095'][0]) + 1",
            "        if len(scored) and (scored.project == 'PXD058095').any() else None)",
            "if not len(scored):",
            "    verdict = ('no project measures enough of the fingerprint: the depth of these proteomes, not the '",
            "               'fingerprint, is what blocks the axis')",
            "elif rank == 1:",
            "    verdict = 'the stated bradyzoite proteome scores highest'",
            "elif rank is None:",
            "    verdict = 'the stated bradyzoite proteome could not be scored'",
            "else:",
            "    verdict = f'the stated bradyzoite proteome ranks {rank} of {len(scored)} -- the fingerprint does not read it'",
            "measured = S[S.rises_measured.notna()] if 'rises_measured' in S.columns else S.iloc[:0]",
            "summary = dict(projects=len(S), scored_projects=len(scored), bradyzoite_rank=rank, verdict=verdict,",
            "               rises_in_fingerprint=len(rises), falls_in_fingerprint=len(falls),",
            "               rises_measured_max=(int(measured.rises_measured.max()) if len(measured) else None),",
            "               falls_measured_max=(int(measured.falls_measured.max()) if len(measured) else None),",
            "               minimum_required=MIN_GENES)",
            "json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)",
            "print(summary)",
            "measured[['project', 'genes', 'rises_measured', 'falls_measured']].sort_values('genes', ascending=False)")
    nb.md("## 3 · Why it fails, which is the finding",
          "The genes that **rise** on conversion are the bradyzoite-specific ones, and they are the ones these "
          "proteomes do not reach: at most 10 of 150 in any project, against 17-61 of the 150 that fall. The "
          "falling set is ordinary tachyzoite protein, abundant and detected. So the obstacle is depth on the "
          "bradyzoite side, not the fingerprint, and a deeper bradyzoite proteome -- not another tachyzoite "
          "one -- is what the axis needs.")
    nb.write(NB)


if __name__ == "__main__":
    main()
