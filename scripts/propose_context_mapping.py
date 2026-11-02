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
HOURS = dict(h=1, hr=1, hrs=1, hour=1, hours=1, hpi=1, d=24, day=24, days=24, dpi=24, w=168, wk=168, wks=168,
             week=168, weeks=168, min=1 / 60, mins=1 / 60, minute=1 / 60, minutes=1 / 60)
TIME_TAGS = {'time', 'timepoint', 'time point', 'time post infection', 'time post-infection', 'infection time',
             'time of infection', 'hours post infection', 'days post infection'}
STAGE_TAGS = {'stage', 'developmental stage', 'life stage', 'parasite stage', 'life cycle stage', 'develpomental stage',
              'developmentalstage', 'dev. stage', 'toxoplasma stage', 'parasite life cycle stage'}

#: First-pass tag review (2026-11-02, Claude, for Einar to spot-check). Spelling variants join the
#: field of their canonical tag; tags describing assay mechanics are marked not_context so they leave
#: the review queue without a value. Replicate and individual identifiers are kept because they decide
#: which samples are independent.
REVIEWED_TAGS = {
    **dict.fromkeys(['toxoplasma parasite strain', 't. gondii strain', 'strain of toxoplasma', 'toxplasma gondii strain',
                     'strain, toxoplasma gondii', 'tg type strain', 't. gondii strain background'], 'parasite_strain'),
    **dict.fromkeys(['parasite genotype', 't. gondii genotype', 'toxoplasma gondii genotype', 'toxoplasma genetic modification',
                     't. gondii strain and genotype', 'complementation', 'overexpression', 'plasmid transfection',
                     'sirna transfection', 'transfection', 'transfected with', 'flp recombinase'], 'parasite_or_cell_perturbation_declared'),
    **dict.fromkeys(['strain of mouse', 'mice strain'], 'host_strain'),
    **dict.fromkeys(['mouse genotype/variation', 'genetics', 'host cell genotype'], 'host_genotype_declared'),
    **dict.fromkeys(['host cell type', 'host cell line', 'host cells', 'host cell_line', 'host cell_type', 'cellline', 'cell',
                     'cell source', 'tissue/cell line', 'monocyte subset', 'host growth'], 'host_cell_or_tissue'),
    **dict.fromkeys(['sex', 'gender', 'gender (for cell source)'], 'host_sex'),
    **dict.fromkeys(['age', 'age weeks', 'age (cell source)'], 'host_age_declared'),
    **dict.fromkeys(['replicate', 'biological replicate', 'biol. replicate', 'replicate information', 'batch', 'donor',
                     'donor id', 'mouse id', 'mouse number', 'clone'], 'unit_or_replicate_identifier'),
    **dict.fromkeys(['treat agent', 'treat dose', 'treat mode', 'infection dose', 'infection mode', 'infection status',
                     'infected with', 'infected', 'stimulus', 'agonist', 'stimulation', 'stimulated', 'drug treatment',
                     'treatment 1', 'treatment 2', 'treated', 'ifn-g priming', 'infrection condition', 'pathogen', 'disease',
                     'disease state', 'mouse model', 'multiplicity of infection', 'moi', 'culture condition', 'growth condition',
                     'culture conditions', 'incubation media', 'ph', 'growth temperature', 'parasite strain or treatment',
                     'infection protocol', 'parasite', 'parasites', 'differentiation protocol'], 'condition_declared'),
    **dict.fromkeys(['chip antibody', 'chip antibody vendor', 'chip antibody cat.#', 'chip antibody cat. #',
                     'chip antibody lot/batch#', 'chip antibody source', 'chip antibody name', 'chip antibody manufacturer',
                     'chip antibody catalog #', 'chip antobody catalogue number', 'antibody', 'library type', 'array design',
                     'rna population', 'fraction', 'sorted fraction', 'molecule subtype', 'gating strategy', 'hashtag number',
                     'immgent igt_10xlane_id', 'cite seq', 'lab', 'exp name', 'organism name', 'organ uberon_id', 'plasmid',
                     'protein_sequence', 'uniprot id', 'expressed protein', 'surface marker', 'sort population',
                     'sorting strategy', 'dna fragmentation method', 'prep description', 'target cells', 't specificity',
                     'description for exp id', 'jax', 'pool of_mice', 'material type', 'passage', 'passages', 'body weight',
                     'hematocrit', 'platlets', 'wbc count', 'cd4 count', 'viral load copies/ml', 'siv status', 'habitat',
                     'ecotype', 'developmental stage (boyes et al. plant cell 2001)', 'tumor'], 'not_context'),
}
#: Time tags whose name carries the unit (value is a bare number) and whether the start is infection.
UNIT_TIME_TAGS = {'days post-infection': (24, True), 'days postinfection': (24, True), 'day of infection': (24, True),
                  'months post infection': (730.5, True), 'week-post treatment': (168, False)}
TIME_TAGS |= {'infection time point', 'infection time', 'culture time', 'incubation time', 'growth time',
              'collection timepoint', 'time collected_after_initial_stimulation', 'duration', 'treatment duration',
              'length of treatment', 'time stimulated_with_lps+ifny'}


STRAIN_TAGS = {'strain', 'strain background', 'genetic background', 'strain/background', 'cell line/strain',
               'strain/genotype'}
CONDITION_TAGS = {'treatment', 'infection', 'condition', 'infection agent', 'agent', 'drug treatment', 'agonist',
                  'stimulus', 'disease state', 'group', 'groups', 'sample group', 'phenotype', 'condition detailed_organ'}
COMMENT_TAGS = {'mouse comments', 'sample comments'}
NOT_AVAILABLE = {'n/a', 'na', 'not applicable', 'unknown', 'missing', '-', 'not available', 'not determined'}
WORD_NUMBERS = dict(one=1, two=2, three=3, four=4, five=5, six=6, seven=7, eight=8, nine=9, ten=10)
#: "day 5", "Day10", "at day 5 post-infection", "post inoculation day 25", "27 days after treatment", "12h later".
TIME_ANY = re.compile(r'(?:\bday\s*(?P<d>\d+(?:\.\d+)?)\b)|(?P<n>\d+(?:\.\d+)?|' + '|'.join(WORD_NUMBERS) +
                      r')\s*(?P<u>h|hr|hrs|hours?|hpi|d|days?|dpi|w|wks?|weeks?|min|mins|minutes)\b', re.I)
REFERENCE = re.compile(r'post[- ]?infection|p\.?i\.?|hpi|dpi|post[- ]?inoculation|after infection|post[- ]?infect', re.I)


def stage_in(value):
    return re.search(r'\b(?:' + STAGES + r')\b|enter?o?pithelial', value, re.I) is not None


def parse_time(value):
    """Hours and whether the value names infection as the start; None without an explicit number and unit."""
    match = TIME_ANY.search(value)
    if not match:
        return None
    if match.group('d'):
        hours = float(match.group('d')) * 24
    else:
        number = match.group('n').lower()
        number = WORD_NUMBERS.get(number) or float(number)
        hours = float(number) * HOURS[match.group('u').lower().rstrip('s') if match.group('u').lower() not in HOURS else match.group('u').lower()]
    return hours, bool(REFERENCE.search(value))


def whole(pattern, value):
    return re.fullmatch(r'\s*(?:T\.? ?gondii\s+)?(?:' + pattern + r')(?:\s+strain)?\s*', value, re.I) is not None


def propose(tag: str, value: str) -> tuple[str, str, str]:
    """(field, proposed value, rule) for one literal pair; field 'review' means no proposal."""
    key = tag.lower()
    if key in CONDITION_TAGS and value.strip().lower() in {'none', 'no', 'untreated', 'uninfected', 'naive', 'mock'}:
        return 'condition_declared', value, f'"{key}" declares no treatment or infection'
    if key == 'infection time point' and value.strip().lower() == 'naive':
        return 'condition_declared', value, 'naive (uninfected) declared'
    if not value or value.strip().lower() in NOT_AVAILABLE | {'none'}:
        return 'not_available', '', 'explicit not-available or empty value'
    if key in COMMENT_TAGS:
        return 'free_text_comment', value, 'free-text comment kept, not interpreted'
    if key in TAG_FIELDS:
        return TAG_FIELDS[key], value, f'tag name "{key}"'
    if key in UNIT_TIME_TAGS and re.fullmatch(r'\s*\d+(?:\.\d+)?\s*', value):
        factor, infection = UNIT_TIME_TAGS[key]
        return ('time_declared_hours', f'{float(value) * factor:g}',
                f'unit from tag "{key}"; ' + ('stated as after infection' if infection else 'reference event needs review'))
    if key in REVIEWED_TAGS and key not in STAGE_TAGS:
        field = REVIEWED_TAGS[key]
        return field, ('' if field == 'not_context' else value), f'first-pass tag review: "{key}"'
    if key in STRAIN_TAGS:
        if whole(PARASITE_STRAINS, value):
            return 'parasite_strain', value, 'strain value matches a T. gondii strain/lineage name'
        if whole(HOST_STRAINS + r'|B6|C3H\w*|AJ|NOD/ShiLtJ', value):
            return 'host_strain', value, 'strain value matches a laboratory host strain name'
        if re.match(r'\s*(?:Tg ?)?(?:type\s+(?:I|II|III)\s+)?(?:' + PARASITE_STRAINS + r')(?:[\s_:.\-(]|$)', value, re.I) \
                or re.search(r'\bRH\b|RHKU80|RH::', value):
            return 'parasite_strain', value, 'strain value begins with or names a T. gondii strain (modified line kept literally)'
        if re.match(r'\s*(?:B6|C57)[\s./\w-]*', value):
            return 'host_strain', value, 'strain value is a C57BL/6 derivative or cross'
        return 'review', '', 'strain value matches neither list'
    if key == 'type' and value.strip() in {'I', 'II', 'III'}:
        return 'parasite_strain', 'Type ' + value.strip(), 'type tag with a T. gondii lineage numeral'
    if key in STAGE_TAGS | {'sample type'}:
        if stage_in(value):
            return 'parasite_stage', value, f'"{key}" value names a parasite life stage (literal kept)'
        timed = parse_time(value)
        if timed and timed[1]:
            return 'time_declared_hours', f'{timed[0]:g}', f'"{key}" value is a time stated as after infection'
    if key in STAGE_TAGS or key == 'cell type':
        if re.fullmatch(r'\s*(?:T\.? ?gondii\s+)?(?:' + STAGES + r')\s*', value, re.I):
            return 'parasite_stage', value, f'"{key}" value names a parasite life stage'
        if key == 'cell type':
            return 'host_cell_or_tissue', value, 'cell type value is not a parasite stage'
        return 'review', '', 'stage value is not a recognized life stage'
    if key in TIME_TAGS:
        match = TIME.match(value)
        timed = None if match else parse_time(value)
        if timed:
            return ('time_declared_hours', f'{timed[0]:g}', 'number and unit within the value; ' +
                    ('stated as after infection' if timed[1] else 'reference event needs review'))
        if match:
            hours = float(match.group(1)) * HOURS[match.group(2).lower()]
            # The clock's start (infection, treatment, plating) is only known when the value says so.
            infection = bool(re.search(r'post[- ]?infection|p\.?i\.?|hpi|dpi', value, re.I))
            return ('time_declared_hours', f'{hours:g}',
                    'explicit number and unit; ' + ('stated as after infection' if infection else 'reference event needs review'))
        return 'review', '', 'time without an explicit number and unit'
    if key in {'genotype', 'genotype/variation'}:
        return 'genotype_declared', value, 'genotype tag; parasite or host owner needs review'
    if key in CONDITION_TAGS:
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
