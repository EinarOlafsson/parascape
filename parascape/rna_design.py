"""Study, sample and run relationships with explicit evidence and uncertainty.

An accession links submitted samples and sequencing runs; it does not establish
independent biological replicates. Title-derived groupings are suggestions only.
"""
from __future__ import annotations

import json
from pathlib import Path
import re

import pandas as pd

from . import expression_contrast as X


def annotate(rows: list[dict], metadata: dict | None = None) -> pd.DataFrame:
    """Annotate source-declared factors and link repeated samples without inventing controls."""
    metadata = metadata or {}
    output = []
    factor_tags = {'strain', 'genotype', 'treatment', 'condition', 'time', 'time point',
                   'timepoint', 'developmental stage', 'stage', 'drug', 'dose', 'tissue', 'host'}
    for row in rows:
        item = dict(row)
        sample = str(row.get('sample_accession') or '')
        source = metadata.get(sample, {})
        attrs = {str(k).lower().replace('_', ' '): str(v) for k, v in source.get('attributes', {}).items()}
        title = str(row.get('sample_title') or row.get('sample_alias') or row['run_accession'])
        factors = [(k, v) for k, v in attrs.items() if k in factor_tags and v.lower() not in ('missing', 'not provided', 'not applicable')]
        group = '; '.join(f'{k}: {v}' for k, v in factors)
        item['condition_group'] = (group + ' | ' if group else '') + X._condition(row)
        item['group_evidence'] = 'source factors with title-suggested grouping; review required' if group else 'suggested from sample title; review required'
        reps = [(k, v) for k, v in attrs.items() if 'replicate' in k]
        text = str(row.get('sample_alias') or '') + ' ' + title
        match = re.search(r'(?i)(?:^|[ _-])(?:biological[ _-]*)?rep(?:licate)?[ _-]*(\d+)(?=$|[\W_])', text)
        item['replicate_label'] = '; '.join(f'{k}: {v}' for k, v in reps) or (match.group(0) if match else 'not specified')
        item['replicate_evidence'] = 'source attribute' if reps else 'source sample label' if match else 'unknown'
        item['replicate_type'] = ('biological (source attribute)' if any('biological' in (k + ' ' + v).lower() for k, v in reps)
                                  else 'technical (source attribute)' if any('technical' in (k + ' ' + v).lower() for k, v in reps)
                                  else 'not specified; do not assume independence')
        controls = ' '.join([title, str(row.get('sample_alias') or '')] + [v for _, v in factors])
        if re.search(r'(?i)(?:^|[\W_])control(?=$|[\W_])', controls):
            item['control_role'] = 'control (source label)'
        elif re.search(r'(?i)\b(wild[ _-]*type|wt|untreated|vehicle|mock)\b', controls):
            item['control_role'] = 'candidate baseline; review required'
        else:
            item['control_role'] = 'not assigned'
        item['control_link'] = 'No verified pairing; choose the matching control within this study'
        item['sample_key'] = sample or f"unknown:{row['study']}:{row['run_accession']}"
        item['metadata_url'] = source.get('url', f'https://www.ebi.ac.uk/ena/browser/view/{sample}' if sample else '')
        item['metadata_status'] = source.get('status', 'run metadata only')
        item['attributes_json'] = json.dumps(source.get('attributes', {}), sort_keys=True)
        output.append(item)
    d = pd.DataFrame(output)
    if d.empty:
        return d
    related = d.groupby('sample_key').run_accession.agg(lambda x: ', '.join(sorted(set(x))))
    d['linked_runs'] = d.sample_key.map(related)
    n = d.groupby('sample_key').run_accession.nunique()
    d['run_relationship'] = d.sample_key.map(n).map(lambda v: 'multiple runs from one submitted sample; independence unconfirmed' if v > 1 else 'one run for this submitted sample')
    reused = d.groupby('run_accession').study.nunique()
    d['shared_across_studies'] = d.run_accession.map(reused) > 1
    d.loc[d.shared_across_studies, 'run_relationship'] += '; same raw run listed in multiple studies'
    return d


def samples(root: str = X.NAS) -> pd.DataFrame:
    """Read the shipped sample design, including explicitly pending priority-study runs."""
    p = Path(root) / 'sample_design.tsv'
    if p.exists():
        return pd.read_csv(p, sep='\t', keep_default_na=False)
    rows = []
    available = X.available(root).drop_duplicates('study')
    for study in available.study:
        first = X.available(root).query('study == @study').iloc[0]
        _, _, _, sheet = X.load(study, first.quantifier, first.reference, root)
        rows.extend(sheet.assign(study=study, processing_status='verified').to_dict('records'))
    return annotate(rows)


def independent_runs(sheet: pd.DataFrame, chosen: list[str]) -> bool:
    """Refuse to count repeated runs of the same submitted sample as independent replicates."""
    if 'sample_accession' not in sheet:
        return True
    ids = sheet[sheet.run_accession.isin(chosen)].sample_accession.fillna('').astype(str)
    return not ids[ids != ''].duplicated().any()
