"""Stage, and with --apply append, reviewed bulk RNA studies to the bounded automation queue.

Rows are rebuilt with omics_quant.ena_runs and the exact omics_quant.plan() filter
(Toxoplasma scientific name, RNA-Seq/OTHER/Ribo-Seq/blank strategy), so new entries
match existing ones. Studies are taken only from the dry-run review in
processing_audit/unplanned_rna_study_review. Single-cell or other per-cell
collections are excluded by name with a reason, because bulk salmon/kallisto
quantification of one run per cell is not the same measurement. The supervisor
lock is held while both plan files and the pinned checksum change together.
"""
import argparse
import csv
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import omics_quant as Q  # noqa: E402

AUDIT = Path('/nas_mnt/carruthers/Einar/parascape_omics/processing_audit')
REVIEW = AUDIT / 'unplanned_rna_study_review'
AUTOMATION = Path.home() / '.cache/parascape-downloads/automation'
EXCLUDED = {'GSE145800': 'Per-cell libraries (11,419 runs for 23 GEO samples); single-cell processing required, not bulk quantification',
            'GSE145080': 'Per-cell libraries (5,499 runs for 8 GEO samples); single-cell processing required, not bulk quantification'}
KEEP_STRATEGIES = ('RNA-Seq', 'OTHER', 'Ribo-Seq', '')
#: Hosts whose transcriptomes are in omics_quant.HOST_TRANSCRIPTS; reads of other hosts would have no
#: host transcripts to compete with, so conserved host genes could be counted as parasite.
REFERENCE_HOSTS = ('Homo sapiens', 'Mus musculus')
HOST_MODE = '--host-labelled' in sys.argv
STAGED = 'staged-host-labelled-rows.jsonl' if HOST_MODE else 'staged-plan-rows.jsonl'
RECEIPT = 'staged-host-labelled-receipt.json' if HOST_MODE else 'staged-receipt.json'


def runs_with_retry(project):
    for attempt in range(4):
        try:
            return Q.ena_runs(project)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(20 * (attempt + 1))


def stage():
    studies = [r for r in csv.DictReader(open(REVIEW / 'studies.tsv'), delimiter='\t') if not r.get('error')]
    catalogue = {r['accession']: r for r in csv.DictReader(open(AUDIT / 'catalogue_coverage.tsv'), delimiter='\t')}
    existing = [json.loads(line) for line in (AUTOMATION / 'plan.jsonl').read_text().splitlines() if line.strip()]
    planned = {(r['study'], r['run_accession']) for r in existing}
    planned_runs = {r['run_accession'] for r in existing}
    rows, decisions = [], []
    for study in studies:
        name = study['study']
        if name in EXCLUDED:
            decisions.append(dict(study=name, decision='excluded', reason=EXCLUDED[name])); continue
        if not study.get('project'):
            decisions.append(dict(study=name, decision='no_project')); continue
        refs = ['ME49', 'ME49_host'] + (['GT1', 'VEG'] if catalogue.get(name, {}).get('strain comparison') == 'True' else [])
        added = 0
        for run in runs_with_retry(study['project']):
            if (name, run['run_accession']) in planned:
                continue
            if HOST_MODE:
                if run.get('scientific_name') not in REFERENCE_HOSTS:
                    continue
                run['run_label'] = 'host_labelled_run: ENA scientific name ' + run['scientific_name'] + '; parasite read fraction may be small'
            elif 'Toxoplasma' not in run.get('scientific_name', 'Toxoplasma'):
                continue
            if run.get('library_strategy', 'RNA-Seq') not in KEEP_STRATEGIES:
                continue
            run.update({'study': name, 'project': study['project'], 'references': refs})
            rows.append(run); planned.add((name, run['run_accession'])); added += 1
        decisions.append(dict(study=name, decision='staged', runs=added, references=refs))
        time.sleep(0.4)
    shared = sorted({r['run_accession'] for r in rows} & planned_runs)
    staged = ''.join(json.dumps(r) + '\n' for r in rows).encode()
    (REVIEW / STAGED).write_bytes(staged)
    receipt = dict(epoch=time.time(), operation='stage_reviewed_bulk_RNA_runs', staged_rows=len(rows),
                   staged_sha256=hashlib.sha256(staged).hexdigest(), studies=decisions,
                   runs_already_planned_under_other_studies=shared,
                   fastq_bytes=sum(sum(int(x) for x in r.get('fastq_bytes', '').split(';') if x.isdigit()) for r in rows),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   scope='Bulk RNA runs with the existing plan filter; host-named dual-RNA runs, non-RNA strategies and single-cell collections stay separate review gaps.')
    (REVIEW / RECEIPT).write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'studies'}), flush=True)
    return receipt


def apply(receipt):
    staged = (REVIEW / STAGED).read_bytes()
    if hashlib.sha256(staged).hexdigest() != receipt['staged_sha256']:
        raise ValueError('Staged rows changed after review')
    stamp = time.strftime('%Y%m%dT%H%M%S')
    with (AUTOMATION / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        config = json.loads((AUTOMATION / 'config.json').read_text())
        plan = (AUTOMATION / 'plan.jsonl').read_bytes()
        if hashlib.sha256(plan).hexdigest() != config['plan_sha256']:
            raise ValueError('Current plan does not match its pinned checksum')
        backup = AUTOMATION / ('before-unplanned-rna-' + stamp); backup.mkdir()
        for name in ('config.json', 'plan.jsonl', 'priority-plan.jsonl'):
            shutil.copy2(AUTOMATION / name, backup / name)
        for name in ('plan.jsonl', 'priority-plan.jsonl'):
            raw = (AUTOMATION / name).read_bytes()
            raw = raw if raw.endswith(b'\n') else raw + b'\n'
            temporary = AUTOMATION / (name + '.tmp'); temporary.write_bytes(raw + staged); temporary.replace(AUTOMATION / name)
        new = (AUTOMATION / 'plan.jsonl').read_bytes()
        rows = [line for line in new.splitlines() if line.strip()]
        config.update(plan_sha256=hashlib.sha256(new).hexdigest(), planned_runs=len(rows),
                      extended=config.get('extended', []) + [dict(epoch=time.time(), added_rows=receipt['staged_rows'],
                          staged_sha256=receipt['staged_sha256'], backup=str(backup),
                          reason='User 2026-11-01: finish all downloads; reviewed unplanned RNA studies appended (' + STAGED + ')')])
        temporary = AUTOMATION / 'config.json.tmp'; temporary.write_text(json.dumps(config, indent=2) + '\n')
        temporary.replace(AUTOMATION / 'config.json')
    with (REVIEW / 'operations.jsonl').open('a') as log:
        log.write(json.dumps(dict(operation='append_reviewed_bulk_RNA_runs', epoch=time.time(), added_rows=receipt['staged_rows'],
                                  staged_sha256=receipt['staged_sha256'], new_plan_sha256=config['plan_sha256'],
                                  planned_runs=config['planned_runs'], backup=str(backup))) + '\n')
    print(json.dumps(dict(applied=True, planned_runs=config['planned_runs'], backup=str(backup))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', help='append the staged rows to the automation queue')
    parser.add_argument('--host-labelled', action='store_true', help='stage host-named runs of reference hosts instead')
    options = parser.parse_args()
    if options.apply:
        apply(json.loads((REVIEW / RECEIPT).read_text()))
    else:
        stage()
