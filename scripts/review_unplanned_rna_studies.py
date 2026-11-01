"""Dry-run review of catalogue RNA studies that are absent from the RNA plan; writes no plan rows.

For every catalogue study with status labelled_parasite_rna_not_in_current_plan, resolve its
sequencing project exactly as omics_quant.plan() does, list every ENA run, and record which runs
the current plan filter would keep and why the others are dropped (host scientific name, library
strategy). Output is a review table and a hashed receipt; admission is a separate decision.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import omics_quant as Q  # noqa: E402

AUDIT = Path('/nas_mnt/carruthers/Einar/parascape_omics/processing_audit')
OUT = AUDIT / 'unplanned_rna_study_review'
KEEP_STRATEGIES = ('RNA-Seq', 'OTHER', 'Ribo-Seq', '')


def main():
    OUT.mkdir(exist_ok=True)
    rows = list(csv.DictReader(open(AUDIT / 'catalogue_coverage.tsv'), delimiter='\t'))
    studies = [r for r in rows if r['coverage_status'] == 'labelled_parasite_rna_not_in_current_plan']
    planned = set()
    with open(Q.PLAN) as fh:
        planned = {json.loads(line)['run_accession'] for line in fh if line.strip()}
    review, runs_out = [], []
    previous, previous_runs = {}, {}
    if (OUT / 'studies.tsv').exists():
        previous = {r['study']: r for r in csv.DictReader(open(OUT / 'studies.tsv'), delimiter='\t') if not r.get('error')}
        for r in csv.DictReader(open(OUT / 'runs.tsv'), delimiter='\t'):
            previous_runs.setdefault(r['study'], []).append(r)
    for row in studies:
        study = row['accession']
        record = dict(study=study, source=row['source'], catalogue_taxa=row['taxa'],
                      catalogue_samples=row['samples'], title=row['title'][:200])
        if study in previous:  # already resolved in an earlier pass; keep its exact rows
            review.append(previous[study]); runs_out.extend(previous_runs.get(study, [])); continue
        try:
            project = Q.geo_project(study) if study.startswith('GSE') else study
            record['project'] = project or ''
            runs = []
            for attempt in range(4):  # ENA returns transient HTTP 500s under load
                try:
                    runs = Q.ena_runs(project) if project else []
                    break
                except Exception:
                    if attempt == 3:
                        raise
                    time.sleep(20 * (attempt + 1))
        except Exception as error:  # network or service failure stays an explicit retry gap
            record.update(project=record.get('project', ''), error=str(error)[:200])
            review.append(record)
            continue
        time.sleep(0.4)
        counts = dict(total=0, already_planned=0, keep_current_filter=0, host_scientific_name=0,
                      other_strategy=0, no_fastq=0)
        names, strategies, keep_bytes, host_bytes = set(), set(), 0, 0
        for run in runs:
            counts['total'] += 1
            name, strategy = run.get('scientific_name', ''), run.get('library_strategy', '')
            names.add(name); strategies.add(strategy)
            size = sum(int(x) for x in run.get('fastq_bytes', '').split(';') if x.isdigit())
            if run['run_accession'] in planned:
                reason = 'already_planned'
            elif strategy not in KEEP_STRATEGIES:
                reason = 'other_strategy'
            elif 'Toxoplasma' not in name:
                reason = 'host_scientific_name'; host_bytes += size
            elif not run.get('fastq_ftp'):
                reason = 'no_fastq'
            else:
                reason = 'keep_current_filter'; keep_bytes += size
            counts[reason] += 1
            runs_out.append(dict(study=study, project=project, run_accession=run['run_accession'],
                                 scientific_name=name, library_strategy=strategy,
                                 library_layout=run.get('library_layout', ''), fastq_bytes=size,
                                 sample_accession=run.get('sample_accession', ''), decision=reason))
        record.update(counts, scientific_names=';'.join(sorted(names)), strategies=';'.join(sorted(strategies)),
                      keep_bytes=keep_bytes, host_named_bytes=host_bytes, error='')
        review.append(record)
        print(json.dumps({k: record.get(k) for k in ('study', 'project', 'total', 'keep_current_filter',
                                                      'host_scientific_name', 'other_strategy')}), flush=True)
    for name, table in (('studies.tsv', review), ('runs.tsv', runs_out)):
        keys = sorted({k for r in table for k in r}) if table else ['study']
        with open(OUT / name, 'w', newline='') as fh:
            w = csv.DictWriter(fh, keys, delimiter='\t'); w.writeheader(); w.writerows(table)
    receipt = dict(epoch=time.time(), operation='dry_run_review_of_unplanned_RNA_studies', studies=len(review),
                   runs=len(runs_out), plan_modified=False,
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   outputs={n: hashlib.sha256((OUT / n).read_bytes()).hexdigest() for n in ('studies.tsv', 'runs.tsv')})
    (OUT / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    main()
