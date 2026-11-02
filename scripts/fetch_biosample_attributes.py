"""Fetch the submitted attributes of queued RNA runs' BioSamples that have no GEO sample record.

About half of the queued RNA runs come from SRA-only studies, so their strain, host, stage, time and
treatment are recorded only as BioSample attributes. This downloads the ENA sample XML for every such
BioSample (in batches), keeps each raw response with its SHA-256, and writes a flat table of the
attributes exactly as submitted, in the same shape as the GEO annotations of the sample registry.
Nothing is interpreted here; scripts/propose_context_mapping.py proposes context fields from them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd

NAS = Path('/nas_mnt/carruthers/Einar/parascape_omics/processing_audit')
REGISTRY = NAS / 'sample_registry'
OUT = NAS / 'biosample_attributes'
URL = 'https://www.ebi.ac.uk/ena/browser/api/xml/{}'


def fetch(accessions, tries=4):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(URL.format(','.join(accessions)), timeout=120) as response:
                return response.read()
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(15 * (attempt + 1))


def parse(raw: bytes):
    """One row per submitted attribute, plus title and organism; values only stripped of edge whitespace."""
    rows = []
    for sample in ET.fromstring(raw).findall('SAMPLE'):
        accession = sample.get('accession')
        def add(tag, value, origin):
            rows.append(dict(biosample=accession, tag=tag, value=re.sub(r'^\s+|\s+$', '', value or ''), origin=origin))
        add('title', sample.findtext('TITLE'), 'sample title')
        add('organism', sample.findtext('SAMPLE_NAME/SCIENTIFIC_NAME'), 'sample organism')
        for attribute in sample.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE'):
            add(re.sub(r'^\s+|\s+$', '', attribute.findtext('TAG') or '') or '(untagged)',
                attribute.findtext('VALUE'), 'biosample attribute')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--batch', type=int, default=50)
    args = parser.parse_args()
    registry = sorted(p for p in REGISTRY.iterdir() if p.is_dir())[-1]
    links = pd.read_parquet(registry / 'run_links.parquet')
    wanted = sorted({b for b in links.loc[links.geo_sample.eq(''), 'biosample'] if re.fullmatch(r'SAM[NED]\w+', b or '')})
    out = OUT / time.strftime('%Y%m%dT%H%M%S'); (out / 'raw').mkdir(parents=True)
    rows, receipts = [], []
    for start in range(0, len(wanted), args.batch):
        batch = wanted[start:start + args.batch]
        raw = fetch(batch)
        digest = hashlib.sha256(raw).hexdigest()
        (out / 'raw' / f'{digest}.xml').write_bytes(raw)
        parsed = parse(raw)
        found = {r['biosample'] for r in parsed}
        receipts.append(dict(requested=batch, returned=sorted(found), sha256=digest, url=URL.format(','.join(batch))))
        rows += parsed
        time.sleep(0.5)
    table = pd.DataFrame(rows)
    table.to_parquet(out / 'attributes.parquet', index=False)
    missing = sorted(set(wanted) - set(table.biosample))
    receipt = dict(operation='fetch_biosample_attributes', epoch=time.time(), registry=str(registry),
                   requested=len(wanted), returned=int(table.biosample.nunique()), missing=missing,
                   attribute_rows=len(table), batches=receipts,
                   attributes_sha256=hashlib.sha256((out / 'attributes.parquet').read_bytes()).hexdigest(),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if k not in ('batches',)})[:600])


if __name__ == '__main__':
    main()
