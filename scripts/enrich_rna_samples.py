"""Record public ENA sample annotations and a transparent study/sample/run design."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parascape.rna_design import annotate


def fetch(sample, cache):
    """Cache the primary sample XML and retain retrieval failures as explicit metadata gaps."""
    url = f'https://www.ebi.ac.uk/ena/browser/api/xml/{sample}'
    p = cache / f'{sample}.xml'
    for attempt in range(3):
        try:
            if p.exists():
                data = p.read_bytes()
            else:
                with urllib.request.urlopen(url, timeout=45) as response:
                    data = response.read()
                ET.fromstring(data)
                p.write_bytes(data)
            root = ET.fromstring(data)
            attributes = {x.findtext('TAG'): x.findtext('VALUE') for x in root.findall('.//SAMPLE_ATTRIBUTE')}
            return sample, {'url': url, 'sha256': hashlib.sha256(data).hexdigest(),
                            'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                            'attributes': attributes, 'status': 'source attributes retrieved'}
        except Exception as exc:
            if attempt == 2:
                return sample, {'url': url, 'status': 'unavailable', 'error': str(exc), 'attributes': {}}
            time.sleep(attempt + 1)


def main():
    """Enrich completed runs and the user-prioritized Trx1 study from frozen plan identities."""
    root = Path(sys.argv[1])
    plan = Path(sys.argv[2])
    cache = Path(sys.argv[3]); cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((root / 'snapshot.json').read_text())
    verified = {'/'.join(key.split('/')[:2]) for key in manifest['inputs']}
    rows = [json.loads(line) for line in plan.read_text().splitlines() if line]
    rows = [{**r, 'processing_status': 'verified' if r['study'] + '/' + r['run_accession'] in verified else 'pending'} for r in rows]
    ids = sorted({r['sample_accession'] for r in rows if r.get('sample_accession')
                  and (r['processing_status'] == 'verified' or r['study'] == 'PRJNA1008424')})
    with ThreadPoolExecutor(max_workers=4) as pool:
        metadata = dict(pool.map(lambda s: fetch(s, cache), ids))
    (root / 'sample_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    d = annotate(rows, metadata)
    d.to_csv(root / 'sample_design.tsv', sep='\t', index=False)
    for filename in ['sample_metadata.json', 'sample_design.tsv']:
        manifest['outputs'][filename] = hashlib.sha256((root / filename).read_bytes()).hexdigest()
    manifest['sample_design'] = {'rows': len(d), 'samples': len(ids),
                                'source_attributes_retrieved': sum(x['status'] == 'source attributes retrieved' for x in metadata.values()),
                                'pending_priority_study': 'PRJNA1008424', 'control_pairing': 'unverified; requires within-study review'}
    (root / 'snapshot.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest['sample_design']), flush=True)


if __name__ == '__main__':
    main()
