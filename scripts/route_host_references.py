"""Give dual RNA-seq from hosts outside ME49_host its own competitive host reference.

ME49_host holds human and mouse transcripts only. Reads from another host would have no host
transcripts to compete with, so conserved host genes could be counted as parasite. This tool
(1) adds the matching omics_quant.HOST_REFERENCES entry to queued rows of studies whose catalogue
taxa name such a host, and (2) stages the host-labelled runs that were held back because their
host was missing, with references ME49 + that host. Applying requires every named index to exist.
Rows already finished keep their outputs; the supervisor sets their receipt aside as a changed plan
row (renamed, never deleted), and the quantifier reruns them with every reference.
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
EXCLUDED = ('GSE145800', 'GSE145080')
OUT = REVIEW / 'host-routing'


def host_of(taxa):
    """The single non-human/mouse host reference named by catalogue taxa, if any."""
    names = {Q.HOST_SPECIES[x.strip()] for x in taxa.replace(',', ';').split(';')
             if Q.HOST_SPECIES.get(x.strip(), 'ME49_host') != 'ME49_host'}
    if len(names) > 1:
        raise ValueError('More than one non-reference host declared: ' + taxa)
    return next(iter(names), None)


def stage():
    OUT.mkdir(exist_ok=True)
    catalogue = {r['accession']: r for r in csv.DictReader(open(AUDIT / 'catalogue_coverage.tsv'), delimiter='\t')}
    plan = [json.loads(line) for line in (AUTOMATION / 'plan.jsonl').read_text().splitlines() if line.strip()]
    amended = []
    for row in plan:
        ref = host_of(catalogue.get(row['study'], {}).get('taxa', ''))
        if ref and ref not in row['references']:
            amended.append(dict(study=row['study'], run_accession=row['run_accession'], add_reference=ref))
    held = [r for r in csv.DictReader(open(REVIEW / 'runs.tsv'), delimiter='\t')
            if r['decision'] == 'host_scientific_name' and r['study'] not in EXCLUDED
            and Q.HOST_SPECIES.get(r['scientific_name'], 'ME49_host') != 'ME49_host']
    planned = {(r['study'], r['run_accession']) for r in plan}
    projects, new_rows = {}, []
    for r in held:
        projects.setdefault((r['study'], r['project']), set()).add(r['run_accession'])
    for (study, project), wanted in projects.items():
        for run in Q.ena_runs(project):
            if run['run_accession'] in wanted and (study, run['run_accession']) not in planned:
                ref = Q.HOST_SPECIES[run['scientific_name']]
                run.update(study=study, project=project, references=['ME49', ref],
                           run_label='host_labelled_run: ENA scientific name ' + run['scientific_name'] + '; parasite read fraction may be small')
                new_rows.append(run); planned.add((study, run['run_accession']))
        time.sleep(0.4)
    staged = ''.join(json.dumps(r) + '\n' for r in new_rows).encode()
    (OUT / 'new-rows.jsonl').write_bytes(staged)
    receipt = dict(epoch=time.time(), operation='stage_host_reference_routing', amended_rows=amended,
                   new_rows=len(new_rows), new_rows_sha256=hashlib.sha256(staged).hexdigest(),
                   references_needed=sorted({a['add_reference'] for a in amended} | {r['references'][1] for r in new_rows}),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: (len(v) if k == 'amended_rows' else v) for k, v in receipt.items()}), flush=True)


def apply():
    receipt = json.loads((OUT / 'receipt.json').read_text())
    staged = (OUT / 'new-rows.jsonl').read_bytes()
    if hashlib.sha256(staged).hexdigest() != receipt['new_rows_sha256']:
        raise ValueError('Staged rows changed after review')
    missing = [ref for ref in receipt['references_needed']
               if not (Path(Q.NAS) / 'references' / ref / 'reference.json').exists()]
    if missing:
        raise ValueError('Build these references first (omics_quant.py host-references): ' + ', '.join(missing))
    amend = {(a['study'], a['run_accession']): a['add_reference'] for a in receipt['amended_rows']}
    stamp = time.strftime('%Y%m%dT%H%M%S')
    with (AUTOMATION / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        config = json.loads((AUTOMATION / 'config.json').read_text())
        if hashlib.sha256((AUTOMATION / 'plan.jsonl').read_bytes()).hexdigest() != config['plan_sha256']:
            raise ValueError('Current plan does not match its pinned checksum')
        backup = AUTOMATION / ('before-host-routing-' + stamp); backup.mkdir()
        for name in ('config.json', 'plan.jsonl', 'priority-plan.jsonl'):
            shutil.copy2(AUTOMATION / name, backup / name)
        for name in ('plan.jsonl', 'priority-plan.jsonl'):
            rows = [json.loads(line) for line in (AUTOMATION / name).read_text().splitlines() if line.strip()]
            for row in rows:
                ref = amend.get((row['study'], row['run_accession']))
                if ref and ref not in row['references']:
                    row['references'] = row['references'] + [ref]
            raw = ''.join(json.dumps(r) + '\n' for r in rows).encode() + staged
            temporary = AUTOMATION / (name + '.tmp'); temporary.write_bytes(raw); temporary.replace(AUTOMATION / name)
        new = (AUTOMATION / 'plan.jsonl').read_bytes()
        config.update(plan_sha256=hashlib.sha256(new).hexdigest(), planned_runs=len(new.splitlines()),
                      extended=config.get('extended', []) + [dict(epoch=time.time(), amended_rows=len(amend),
                          added_rows=receipt['new_rows'], backup=str(backup), reason='Competitive host references for hosts outside ME49_host')])
        temporary = AUTOMATION / 'config.json.tmp'; temporary.write_text(json.dumps(config, indent=2) + '\n')
        temporary.replace(AUTOMATION / 'config.json')
    with (OUT / 'operations.jsonl').open('a') as log:
        log.write(json.dumps(dict(operation='apply_host_reference_routing', epoch=time.time(), amended_rows=len(amend),
                                  added_rows=receipt['new_rows'], new_plan_sha256=config['plan_sha256'], backup=str(backup))) + '\n')
    print(json.dumps(dict(applied=True, planned_runs=config['planned_runs'], backup=str(backup))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    apply() if parser.parse_args().apply else stage()
