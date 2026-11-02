"""Build the literal sample registry: every validated GEO sample, its author annotations and its runs.

Fair comparison across studies needs to know, per biological sample, its strain, stage, host, time
and treatment. Authors record these as free-text GEO "characteristics" tags, so this first stage
keeps every tag and value exactly as submitted (only surrounding whitespace is removed) and adds no
interpretation: no controlled vocabulary, no inferred replicate or control relationships. It also
links each sample to the queued RNA runs through its BioSample accession and writes an inventory of
the tag names authors actually use, which is what a later, explicit vocabulary mapping is built from.

Inputs (read-only): NAS processing_audit/all_geo_samples/<GSE>/geo-samples.json and validation.json,
and the automation RNA plan. Outputs: processing_audit/sample_registry/<stamp>/ with samples.parquet,
annotations.parquet (one row per sample tag), tag_inventory.csv, run_links.parquet and a receipt.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path
import re
import time

import pandas as pd

NAS = Path('/nas_mnt/carruthers/Einar/parascape_omics')
GEO = NAS / 'processing_audit' / 'all_geo_samples'
PLAN = Path.home() / '.cache/parascape-downloads/automation/plan.jsonl'
OUT = NAS / 'processing_audit' / 'sample_registry'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(value) -> str:
    """Remove the surrounding whitespace MINiML adds; internal text is kept exactly."""
    return re.sub(r'^\s+|\s+$', '', str(value or ''))


def study_records(folder: Path):
    """Yield sample rows and annotation rows for one validated GEO series."""
    validation = json.loads((folder / 'validation.json').read_text())
    if not validation.get('exact_sample_ids_reconciled'):
        return [], [], dict(study=folder.name, status='sample_ids_not_reconciled')
    source = json.loads((folder / 'geo-samples.json').read_text())
    samples, annotations = [], []
    for sample in source['samples']:
        gsm = sample['geo_accession']
        relations = {r.get('type'): r.get('target', '') for r in sample.get('relations', [])}
        biosample = relations.get('BioSample', '').rstrip('/').rsplit('/', 1)[-1] if relations.get('BioSample') else ''
        sra = re.search(r'term=([SED]RX\d+)', relations.get('SRA', '') or '')
        organisms = sorted({o for c in sample.get('channels', []) for o in c.get('organisms', [])})
        samples.append(dict(study=folder.name, geo_sample=gsm, title=clean(sample.get('title')),
                            sample_type=clean(sample.get('sample_type')), platforms=';'.join(sample.get('platform_refs', [])),
                            channels=len(sample.get('channels', [])), organisms=';'.join(organisms),
                            biosample=biosample, sra_experiment=sra.group(1) if sra else '',
                            source_metadata_sha256=sha(folder / 'geo-samples.json')))
        for channel in sample.get('channels', []):
            position = clean(channel.get('position')) or '1'
            annotations.append(dict(study=folder.name, geo_sample=gsm, channel=position, tag='source_name',
                                    value=clean(channel.get('source')), origin='channel source'))
            for organism in channel.get('organisms', []):
                annotations.append(dict(study=folder.name, geo_sample=gsm, channel=position, tag='organism',
                                        value=clean(organism), origin='channel organism'))
            for item in channel.get('characteristics', []):
                annotations.append(dict(study=folder.name, geo_sample=gsm, channel=position,
                                        tag=clean(item.get('tag')) or '(untagged)', value=clean(item.get('value')),
                                        origin='characteristics'))
    return samples, annotations, dict(study=folder.name, status='included', samples=len(samples))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--out', type=Path, default=None)
    args = parser.parse_args()
    out = args.out or OUT / time.strftime('%Y%m%dT%H%M%S')
    out.mkdir(parents=True, exist_ok=False)
    samples, annotations, studies = [], [], []
    for folder in sorted(p for p in GEO.iterdir() if p.is_dir() and p.name.startswith('GSE')):
        if not (folder / 'geo-samples.json').exists():
            studies.append(dict(study=folder.name, status='no_sample_metadata')); continue
        s, a, info = study_records(folder)
        samples += s; annotations += a; studies.append(info)
    samples, annotations = pd.DataFrame(samples), pd.DataFrame(annotations)
    # A SuperSeries re-lists the samples of its SubSeries, so one GEO sample can belong to several
    # series. Keep one sample row listing every series; annotations are kept once when every series
    # carries the same ones, and the sample is flagged for review when they differ.
    fields = ['channel', 'tag', 'value', 'origin']
    signature = (annotations.sort_values(['geo_sample', 'study'] + fields)
                 .groupby(['geo_sample', 'study'])[fields]
                 .apply(lambda f: hashlib.sha256(f.to_csv(index=False).encode()).hexdigest()))
    differing = signature.groupby(level=0).nunique()
    conflicts = set(differing[differing > 1].index)
    series = samples.groupby('geo_sample').study.apply(lambda v: ';'.join(sorted(set(v))))
    first = samples.sort_values('study').drop_duplicates('geo_sample').set_index('geo_sample')
    first['series'] = series
    first['annotation_status'] = ['differs between series; review' if g in conflicts else 'consistent'
                                  for g in first.index]
    samples = first.reset_index()
    keep = samples.set_index('geo_sample').study
    annotations = annotations[annotations.study.eq(annotations.geo_sample.map(keep))].reset_index(drop=True)
    plan = [json.loads(line) for line in PLAN.read_text().splitlines() if line.strip()]
    links = pd.DataFrame([dict(study=r['study'], run_accession=r['run_accession'], biosample=r.get('sample_accession', ''),
                               sample_alias=r.get('sample_alias', ''), sample_title=r.get('sample_title', ''),
                               scientific_name=r.get('scientific_name', ''), library_layout=r.get('library_layout', ''))
                          for r in plan])
    by_biosample = samples[samples.biosample != ''].set_index('biosample').geo_sample
    if by_biosample.index.duplicated().any():
        raise ValueError('A BioSample is linked to more than one GEO sample')
    links['geo_sample'] = links.biosample.map(by_biosample).fillna('')
    links['link'] = links.geo_sample.ne('').map({True: 'exact BioSample match', False: 'no GEO sample (SRA-only study or unlinked)'})
    inventory = (annotations[annotations.origin == 'characteristics']
                 .assign(tag_key=lambda f: f.tag.str.lower())
                 .groupby('tag_key').agg(samples=('geo_sample', 'nunique'), studies=('study', 'nunique'),
                                         spellings=('tag', lambda t: ' | '.join(sorted(set(t))[:5])),
                                         example_values=('value', lambda v: ' | '.join(list(dict.fromkeys(v))[:5])))
                 .sort_values('samples', ascending=False))
    samples.to_parquet(out / 'samples.parquet', index=False)
    annotations.to_parquet(out / 'annotations.parquet', index=False)
    links.to_parquet(out / 'run_links.parquet', index=False)
    inventory.to_csv(out / 'tag_inventory.csv')
    receipt = dict(operation='build_literal_sample_registry', epoch=time.time(), studies=len(studies),
                   studies_included=sum(1 for s in studies if s['status'] == 'included'),
                   study_gaps=[s for s in studies if s['status'] != 'included'], samples=len(samples),
                   annotation_rows=len(annotations), distinct_tags=len(inventory),
                   samples_in_several_series=int(samples.series.str.contains(';').sum()),
                   samples_with_differing_annotations=int(samples.annotation_status.ne('consistent').sum()),
                   plan_rows=len(links), plan_rows_linked=int(links.geo_sample.ne('').sum()),
                   plan_sha256=sha(PLAN), builder_sha256=sha(Path(__file__)),
                   outputs={p.name: sha(p) for p in out.iterdir()},
                   scope='Literal author annotations only; no vocabulary mapping, replicate, control or independence inference.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('outputs', 'study_gaps')}))


if __name__ == '__main__':
    main()
