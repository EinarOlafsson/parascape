"""Propose, for review, how literal GEO sample annotations map onto comparison context fields.

The literal registry (scripts/build_sample_registry.py) keeps author annotations verbatim. To compare
studies fairly each sample needs explicit context: parasite strain, parasite genotype or perturbation,
parasite stage, host species, host strain, host cell or tissue, time after infection and treatment.
Authors use the same tag for different things ("strain" is the host in some studies and the parasite in
others), so a tag name alone is not enough. This script proposes a field for each distinct tag/value
pair only through a named, narrow rule; every other pair stays "review" with no value. Nothing here is
final: the proposal is meant to be reviewed and the reviewed table, not this output, is what later
stages use.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import time

import pandas as pd

REGISTRY = Path('/nas_mnt/carruthers/Einar/parascape_omics/processing_audit/sample_registry')

#: Tags whose name alone fixes the field. Values are kept literally.
TAG_FIELDS = {
    'parasite strain': 'parasite_strain', 'toxoplasma strain': 'parasite_strain',
    'mouse strain': 'host_strain', 'host strain': 'host_strain', 'rat strain': 'host_strain',
    'cell line': 'host_cell_or_tissue', 'tissue': 'host_cell_or_tissue', 'host cell': 'host_cell_or_tissue',
    'organ': 'host_cell_or_tissue', 'host': 'host_species', 'organism': 'organism_declared',
}
#: Recognized T. gondii strain or lineage names (whole-word, case-insensitive) for the ambiguous "strain" tag.
PARASITE_STRAINS = r'RH|RHΔku80|ME-?49|PRU|Prugniaud|VEG|GT-?1|CTG|CEP|CZ1|TgCkUg2|COUG|MAS|FOU|GUY-\w+|Type ?(I|II|III|1|2|3)'
#: Common laboratory host strain names for the same tag.
HOST_STRAINS = (r'C57BL/?6\w*|BALB/?c\w*|CBA/?J|CD-?1|Swiss( Webster)?|DBA/?2\w*|129\w*|NOD\w*|A/J|'
                r'Lewis|Brown Norway|Sprague[- ]Dawley|Wistar|Fischer ?344|F344')
#: Parasite life stages, for stage-like tags and for cell type values that name a stage.
STAGES = r'tachyzoites?|bradyzoites?|sporozoites?|merozoites?|oocysts?|tissue cysts?|cysts?|zoites|schizonts?|gametes?|enteroepithelial'
TIME = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|hpi|d|day|days|dpi|w|wk|week|weeks|min|mins|minutes)\b'
                  r'(?:\s*(?:post[- ]?infection|p\.?i\.?))?\s*$', re.I)
HOURS = dict(h=1, hr=1, hrs=1, hour=1, hours=1, hpi=1, d=24, day=24, days=24, dpi=24, w=168, wk=168, week=168,
             weeks=168, min=1 / 60, mins=1 / 60, minutes=1 / 60)
TIME_TAGS = {'time', 'timepoint', 'time point', 'time post infection', 'time post-infection', 'infection time',
             'time of infection', 'hours post infection', 'days post infection'}
STAGE_TAGS = {'stage', 'developmental stage', 'life stage', 'parasite stage', 'life cycle stage'}


def whole(pattern, value):
    return re.fullmatch(r'\s*(?:T\.? ?gondii\s+)?(?:' + pattern + r')(?:\s+strain)?\s*', value, re.I) is not None


def propose(tag: str, value: str) -> tuple[str, str, str]:
    """(field, proposed value, rule) for one literal pair; field 'review' means no proposal."""
    key = tag.lower()
    if not value or value.lower() in {'n/a', 'na', 'none', 'not applicable', 'unknown', 'missing', '-'}:
        return 'review', '', 'empty or explicit not-available value'
    if key in TAG_FIELDS:
        return TAG_FIELDS[key], value, f'tag name "{key}"'
    if key == 'strain':
        if whole(PARASITE_STRAINS, value):
            return 'parasite_strain', value, 'strain value matches a T. gondii strain/lineage name'
        if whole(HOST_STRAINS, value):
            return 'host_strain', value, 'strain value matches a laboratory host strain name'
        return 'review', '', 'strain value matches neither list'
    if key in STAGE_TAGS or key == 'cell type':
        if re.fullmatch(r'\s*(?:T\.? ?gondii\s+)?(?:' + STAGES + r')\s*', value, re.I):
            return 'parasite_stage', value, f'"{key}" value names a parasite life stage'
        if key == 'cell type':
            return 'host_cell_or_tissue', value, 'cell type value is not a parasite stage'
        return 'review', '', 'stage value is not a recognized life stage'
    if key in TIME_TAGS:
        match = TIME.match(value)
        if match:
            hours = float(match.group(1)) * HOURS[match.group(2).lower()]
            # The clock's start (infection, treatment, plating) is only known when the value says so.
            infection = bool(re.search(r'post[- ]?infection|p\.?i\.?|hpi|dpi', value, re.I))
            return ('time_declared_hours', f'{hours:g}',
                    'explicit number and unit; ' + ('stated as after infection' if infection else 'reference event needs review'))
        return 'review', '', 'time without an explicit number and unit'
    if key in {'genotype', 'genotype/variation'}:
        return 'genotype_declared', value, 'genotype tag; parasite or host owner needs review'
    if key in {'treatment', 'infection', 'condition', 'infection agent', 'agent'}:
        return 'condition_declared', value, f'"{key}" tag; infection, drug and control roles need review'
    return 'review', '', 'tag outside the proposal rules'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--registry', type=Path, default=None, help='registry build folder (default: newest)')
    args = parser.parse_args()
    folder = args.registry or sorted(p for p in REGISTRY.iterdir() if p.is_dir())[-1]
    annotations = pd.read_parquet(folder / 'annotations.parquet')
    literal = annotations[annotations.origin.isin(['characteristics', 'channel organism'])]
    pairs = (literal.groupby(['tag', 'value']).agg(samples=('geo_sample', 'nunique'), studies=('study', 'nunique'))
             .reset_index())
    proposals = pairs.apply(lambda r: propose(r.tag, r.value), axis=1, result_type='expand')
    pairs[['field', 'proposed_value', 'rule']] = proposals
    pairs['review_decision'] = ''                   # to be filled by the reviewer: accept / change / reject
    pairs = pairs.sort_values(['field', 'samples'], ascending=[True, False])
    out = folder / 'context_mapping_proposal.csv'
    pairs.to_csv(out, index=False)
    context = pairs[~pairs.field.isin(['review', 'organism_declared'])][['tag', 'value']]
    covered = literal.merge(context, on=['tag', 'value']).geo_sample.nunique()   # organism alone is not context
    summary = dict(operation='propose_context_mapping', epoch=time.time(), registry=str(folder), pairs=len(pairs),
                   proposed=int((pairs.field != 'review').sum()), review=int((pairs.field == 'review').sum()),
                   samples_with_a_context_proposal=int(covered), samples=int(annotations.geo_sample.nunique()),
                   by_field=pairs.groupby('field').samples.sum().astype(int).to_dict(),
                   proposal_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   scope='Proposal for review only; no proposal is used until reviewed.')
    (folder / 'context_mapping_proposal.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
