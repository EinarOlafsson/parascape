"""Work through every selected PRIDE project with ms_uniform_search.py, smallest first, until all are done.

Selection: processing_audit/pride_ms_uniform_selection.json (Toxoplasma-related projects; see its rule).
A project starts only while the NAS keeps NAS_RESERVE_GB free after its RAW files (+30%) land; a failed
project is recorded and retried on a later pass, never blocking the others.
"""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
NAS = Path('/nas_mnt/carruthers/Einar/parascape_omics')
SELECTION = NAS / 'processing_audit/pride_ms_uniform_selection.json'
OUT = NAS / 'prepared_ms_uniform'
STATE = Path.home() / '.local/state/parascape-ms-uniform'
PYTHON = str(Path.home() / 'anaconda3/envs/parascape-ms/bin/python')
NAS_RESERVE_GB = 200


def status(**fields):
    STATE.mkdir(parents=True, exist_ok=True)
    (STATE / 'status.json').write_text(json.dumps(dict(updated=time.strftime('%Y-%m-%dT%H:%M:%S'), **fields), indent=2) + '\n')


def main():
    threads = sys.argv[1] if len(sys.argv) > 1 else '4'
    while True:
        projects = json.loads(SELECTION.read_text())['keep']
        handled = lambda q: (OUT / q['project'] / 'done.json').exists() or (OUT / q['project'] / 'unsupported.json').exists()
        pending = [p for p in projects if not handled(p)]
        if not pending:
            status(phase='complete', projects=len(projects)); return
        failures = {}
        for p in pending:
            need = p['raw_bytes'] * 1.3 / 1e9 + NAS_RESERVE_GB
            while shutil.disk_usage(NAS).free / 1e9 < need:
                status(phase='waiting_for_nas_space', project=p['project'], need_gb=round(need),
                       free_gb=round(shutil.disk_usage(NAS).free / 1e9))
                time.sleep(600)
            status(phase='processing', project=p['project'], raw_gb=round(p['raw_bytes'] / 1e9, 1),
                   done=len(projects) - len(pending), total=len(projects))
            result = subprocess.run([PYTHON, '-u', str(HERE / 'ms_uniform_search.py'), p['project'], '--threads', threads])
            if result.returncode:
                failures[p['project']] = result.returncode
                with (STATE / 'failures.jsonl').open('a') as log:
                    log.write(json.dumps(dict(project=p['project'], code=result.returncode, epoch=time.time())) + '\n')
            pending_now = [q for q in projects if not handled(q)]
            status(phase='processing', last=p['project'], done=len(projects) - len(pending_now), total=len(projects),
                   failed_this_pass=len(failures))
        time.sleep(3600)                                   # retry failures on the next pass


if __name__ == '__main__':
    main()
