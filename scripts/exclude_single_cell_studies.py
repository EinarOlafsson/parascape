"""Take single-cell collections out of the bulk RNA queue, under the supervisor lock with a backup.

Bulk salmon/kallisto quantification of a droplet library (barcode/UMI reads beside one cDNA read) is not the
same measurement as bulk RNA-seq, and the paired-read split writes nothing for such spots (GSE229505,
dual perturb-seq, failed 2026-11-03 with 'no reads written'). Einar chose deposited single-cell matrices for
these studies; this removes their pending rows from plan.jsonl and priority-plan.jsonl and re-pins the
checksum. Runs already quantified keep their results; a removed run is never deleted from the store.
"""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import time

AUTOMATION = Path.home() / '.cache/parascape-downloads/automation'
AUDIT = Path('/nas_mnt/carruthers/Einar/parascape_omics/processing_audit')
NAS = Path('/nas_mnt/carruthers/Einar/parascape_omics')
#: Study -> why it is single-cell (from the GEO/BioProject title or the archive layout).
SINGLE_CELL = {
    'GSE229505': 'dual perturb-seq (pooled CRISPR + 10x single-cell RNA-seq); 1 biological + 3 technical reads per spot',
    'GSE346252': 'single-cell transcriptomic profiling (title)',
    'PRJNA1047766': 'single-cell RNA-seq (title)',
    'PRJDB12314': 'scRNA-seq (title)',
    'PRJNA1053644': 'single-cell survey (title)',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--apply', action='store_true', help='rewrite the queue (default: report only)')
    args = parser.parse_args()
    stamp = time.strftime('%Y%m%dT%H%M%S')
    with (AUTOMATION / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        config = json.loads((AUTOMATION / 'config.json').read_text())
        if hashlib.sha256((AUTOMATION / 'plan.jsonl').read_bytes()).hexdigest() != config['plan_sha256']:
            raise ValueError('Current plan does not match its pinned checksum')
        backup = AUTOMATION / ('before-single-cell-exclusion-' + stamp)
        if args.apply:
            backup.mkdir()
            for name in ('config.json', 'plan.jsonl', 'priority-plan.jsonl'):
                if (AUTOMATION / name).exists():
                    shutil.copy2(AUTOMATION / name, backup / name)
        removed = {}
        for name in ('plan.jsonl', 'priority-plan.jsonl'):
            path = AUTOMATION / name
            if not path.exists():
                continue
            rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
            done = lambda r: (NAS / 'rnaseq' / r['study'] / r['run_accession'] / 'done.json').exists()
            drop = [r for r in rows if r['study'] in SINGLE_CELL and not done(r)]
            removed[name] = len(drop)
            if args.apply and drop:
                keep = [r for r in rows if r not in drop]
                temporary = path.with_suffix('.tmp')
                temporary.write_text(''.join(json.dumps(r) + '\n' for r in keep)); temporary.replace(path)
        print(json.dumps(dict(apply=args.apply, removed=removed)))
        if not args.apply:
            return
        new = (AUTOMATION / 'plan.jsonl').read_bytes()
        config.update(plan_sha256=hashlib.sha256(new).hexdigest(), planned_runs=len(new.splitlines()),
                      extended=config.get('extended', []) + [dict(epoch=time.time(), removed_rows=removed,
                          backup=str(backup), reason='Single-cell collections use deposited matrices, not bulk quantification',
                          studies=SINGLE_CELL)])
        temporary = AUTOMATION / 'config.json.tmp'; temporary.write_text(json.dumps(config, indent=2) + '\n')
        temporary.replace(AUTOMATION / 'config.json')
    with (AUDIT / 'single_cell_exclusions.jsonl').open('a') as log:
        log.write(json.dumps(dict(epoch=time.time(), removed=removed, studies=SINGLE_CELL,
                                  plan_sha256=config['plan_sha256'], backup=str(backup))) + '\n')


if __name__ == '__main__':
    main()
