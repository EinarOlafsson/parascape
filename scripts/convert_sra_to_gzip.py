"""Convert one cached SRA archive to gzipped FASTQ for runs too large to extract on local scratch.

fasterq-dump needs several times the archive size in uncompressed scratch and failed on the NFS-mounted
NAS while joining its temporary files; ENA's gzip copy can be unreachable. fastq-dump --gzip writes
compressed reads sequentially, which works on the NAS. The output is used by omics_quant.py only after
this script has validated the archive, checked every gzip stream and written done.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

BIN = Path.home() / 'anaconda3/envs/parascape-quant/bin'
CACHE = Path.home() / '.cache/parascape-downloads'
OUT = Path('/nas_mnt/carruthers/Einar/parascape_omics/preconverted_fastq')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 22), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('run')
    args = parser.parse_args()
    run = args.run
    sra = (CACHE / f'{run}.sra').resolve(strict=True)
    folder = OUT / run
    if (folder / 'done.json').exists():
        print(json.dumps(dict(run=run, status='already converted'))); return
    folder.mkdir(parents=True, exist_ok=True)
    for stale in folder.glob(f'{run}*'):
        stale.unlink()
    started = time.time()
    subprocess.run([str(BIN / 'vdb-validate'), str(sra)], check=True, capture_output=True)
    subprocess.run([str(BIN / 'fastq-dump'), '--split-3', '--gzip', '-O', str(folder), str(sra)],
                   check=True, stdout=subprocess.DEVNULL)
    pair = [folder / f'{run}_{k}.fastq.gz' for k in (1, 2)]
    single = folder / f'{run}.fastq.gz'
    files = pair if all(p.exists() for p in pair) else [single]
    if not all(p.exists() and p.stat().st_size > 0 for p in files):
        raise RuntimeError('fastq-dump produced no reads')
    for path in files:
        subprocess.run([str(BIN / 'pigz'), '-t', str(path)], check=True)
    record = dict(run=run, source_archive=str(sra), source_archive_bytes=sra.stat().st_size,
                  files={p.name: dict(bytes=p.stat().st_size, sha256=sha256(p)) for p in files},
                  layout='paired' if files == pair else 'single', tool='fastq-dump --split-3 --gzip 3.4.1',
                  checks='vdb-validate on the archive; pigz -t on every output', seconds=round(time.time() - started),
                  source_sha256=sha256(Path(__file__)))
    (folder / 'done.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record))


if __name__ == '__main__':
    main()
