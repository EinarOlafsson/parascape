"""Acquire reviewed UniProt activity evidence without treating missing labels as negatives."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import organisms as O  # noqa: E402


def sha(path):
    """Hash a bounded source snapshot."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def acquire(output, archive):
    """Pin complete reviewed species queries, preserving source headers and accession lists."""
    output, archive = Path(output), Path(archive)
    if output.exists() or archive.exists():
        raise ValueError('Use fresh immutable output and source archive directories')
    output.mkdir(parents=True)
    archive.mkdir(parents=True)
    receipts = {}
    for organism in (O.TOXOPLASMA, O.FALCIPARUM):
        species = O.get(organism).species
        query = '(organism_name:"' + species + '") AND (reviewed:true)'
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urlencode({'query': query, 'format': 'json'})
        request = Request(url, headers={'User-Agent': 'Starplast source review; public bounded organism query'})
        with urlopen(request, timeout=60) as response:
            payload = response.read(30 * 1024 * 1024 + 1)
            headers = dict(response.headers.items())
            final_url = response.geturl()
        if len(payload) > 30 * 1024 * 1024:
            raise ValueError('Refuse unexpectedly large reviewed-organism source response')
        data = json.loads(payload)
        if not isinstance(data.get('results'), list):
            raise ValueError('Expected complete UniProt JSON stream results')
        path = archive / (organism + '_reviewed.json')
        path.write_bytes(payload)
        receipts[organism] = {'url': url, 'final_url': final_url, 'query': query, 'path': str(path),
            'sha256': sha(path), 'bytes': len(payload), 'headers': headers,
            'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'entries': len(data['results']),
            'taxa': sorted({(item['organism']['taxonId'], item['organism']['scientificName']) for item in data['results']})}
    (archive / 'receipts.json').write_text(json.dumps(receipts, indent=2) + '\n')
    (archive / 'URLS.txt').write_text('\n'.join(item['url'] for item in receipts.values()) + '\n')
    (archive / 'SHA256SUMS.txt').write_text('\n'.join(item['sha256'] + '  ' + Path(item['path']).name for item in receipts.values()) + '\n')
    (output / 'source_receipts.json').write_text(json.dumps(receipts, indent=2) + '\n')
    return receipts


def main():
    """Execute and save the actual acquisition notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--archive', required=True, type=Path)
    args = parser.parse_args()
    nb = ExecutedNotebook('Reviewed parasite catalytic activity evidence source review')
    nb.md('Complete reviewed UniProt species queries. Experimental catalytic evidence must be attached to the reaction/EC, not inferred from reviewed status or from a protein-level literature reference. Source absence is unknown; reviewed entries are an ascertainment-selected evidence source, not a genome-wide negative set. No installed annotations or inferred claims are modified.')
    nb.code('from scripts.review_functional_activity_truth import acquire',
        f'receipts=acquire({str(args.out)!r}, {str(args.archive)!r})',
        "{organism: {'entries': receipt['entries'], 'taxa': receipt['taxa']} for organism, receipt in receipts.items()}")
    nb.write(str(args.out / 'review.ipynb'))
    print({organism: receipt['entries'] for organism, receipt in receipts.items()})


if __name__ == '__main__':
    main()
