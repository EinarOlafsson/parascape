"""Keep a private copy of the derived results ToxoHub v4 serves for the Gubbels/Zarringhalam preprint.

Wu et al. 2026, "A MultiOmics time-resolved systems map of gene regulation in Toxoplasma gondii",
bioRxiv 10.64898/2026.10.29.762658. Its reads (22 scRNA-seq libraries, scATAC-seq, CUT&RUN) are not
deposited; the paper points only to https://umbibio.math.umb.edu/toxohub/v4/. Einar (2026-11-03) chose to
harvest the app's JSON API, one request a second, flagged as unreleased preprint data: kept on the NAS for
our own analysis, never shipped in a release, and replaced by the uniform reprocessing once GEO has the reads.

Every response is stored verbatim with its URL, time and sha256; finished paths are skipped on a rerun.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

API = 'https://umbibio.math.umb.edu/toxohub/v4/api/'
OUT = Path('/nas_mnt/carruthers/Einar/parascape_omics/preprint_results/toxohub_v4')
FLAG = dict(status='unreleased preprint results; not for redistribution or release',
            source='ToxoHub v4 API', preprint='10.64898/2026.10.29.762658',
            decision='Einar 2026-11-03: harvest, flag, keep out of releases until GEO deposit')


def fetch(path: str, pause: float) -> dict | None:
    """One API path; None for a 404 (the app treats it as 'not found')."""
    for attempt in range(4):
        try:
            request = urllib.request.Request(API + path, headers={'User-Agent': 'ParaScape research harvest (einar.olafsson@gmail.com)'})
            with urllib.request.urlopen(request, timeout=120) as response:
                body = response.read()
            time.sleep(pause)
            return dict(body=body)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                time.sleep(pause); return None
            if attempt == 3:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
        time.sleep(30 * (attempt + 1))


def store(path: str, pause: float) -> dict | None:
    target = OUT / 'api' / (urllib.parse.quote(path, safe='') + '.json')
    if target.exists():
        return json.loads(target.read_text())['data']
    got = fetch(path, pause)
    if got is None:
        return None
    data = json.loads(got['body'])
    target.parent.mkdir(parents=True, exist_ok=True)
    record = dict(url=API + path, fetched=time.strftime('%Y-%m-%dT%H:%M:%S'),
                  sha256=hashlib.sha256(got['body']).hexdigest(), flag=FLAG['status'], data=data)
    part = target.with_suffix('.part'); part.write_text(json.dumps(record) + '\n'); part.replace(target)
    return data


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--pause', type=float, default=1.0, help='seconds between requests')
    ap.add_argument('--no-genes', action='store_true', help='skip the per-gene pages (about 7,500 requests)')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'FLAG.json').write_text(json.dumps(FLAG, indent=2) + '\n')
    meta = store('meta', args.pause)
    for path in ('overview', 'regulators', 'communities', 'pairs', 'genes', 'umap/meta'):
        store(path, args.pause)
    windows = int(meta['stats'].get('windows', 24))
    for r in store('regulators', args.pause)['rows']:
        rid = urllib.parse.quote(r['id'], safe='')
        store(f'regulator/{rid}', args.pause)
        store(f'regulator/{rid}/targets', args.pause)
        store(f'regulator/{rid}/de', args.pause)                       # whole cycle
        for t in range(windows):
            store(f'regulator/{rid}/de?t={t}', args.pause)
    for c in store('communities', args.pause)['rows']:
        store(f"community/{urllib.parse.quote(str(c['cid']), safe='')}", args.pause)
    for p in store('pairs', args.pause)['rows']:
        store(f"pair/{urllib.parse.quote(p['id'], safe='')}", args.pause)
    if not args.no_genes:
        genes = store('genes', args.pause)
        column = genes['cols'].index('gene')
        for row in genes['rows']:
            store(f"gene/{urllib.parse.quote(str(row[column]), safe='')}", args.pause)
    files = sorted((OUT / 'api').glob('*.json'))
    (OUT / 'manifest.json').write_text(json.dumps(dict(FLAG, app=meta.get('app'), data=meta.get('data'),
        bundle=meta.get('bundle'), files=len(files), finished=time.strftime('%Y-%m-%dT%H:%M:%S')), indent=2) + '\n')
    print(json.dumps(dict(files=len(files), out=str(OUT))))


if __name__ == '__main__':
    main()
