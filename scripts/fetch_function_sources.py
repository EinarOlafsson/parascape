"""Download the raw evidence for the FUNCTION truth set (instruction 72) and pin it with retrieval times.

Three sources, kept apart on purpose (descending evidence strength):

1. **Experimental GO** -- QuickGO annotation search for taxon 5811 (descendants), aspects
   `biological_process` and `molecular_function`, restricted to experimental evidence codes
   (ECO:0000269 and the specific children we were asked for: 314 IDA, 353 IPI, 315 IMP, 316 IGI,
   270 IEP), each code with `evidenceCodeUsage=descendants`.
2. **Reviewed UniProt (Swiss-Prot)** entries for taxon 5811, with the curator-written `cc_function`
   comment, domains, GO terms, keywords and the PubMed ids UniProt cites. Paged through the
   `Link` header exactly as the REST API asks.
3. Phenotype columns come from `nodes.parquet` and need no download.

Nothing else is fetched. InterPro/Pfam, ToxoDB product names and IEA/ISS/ISO GO annotations are
deliberately **not** downloaded: they are sequence-based predictions, which is what the information
map must predict rather than learn from.

Politeness: at most 3 requests/s (we sleep 0.5 s between calls), one process, no retries in a loop.
Every response body is written to disk verbatim together with the URL, the HTTP status and the
retrieval time in UTC, plus a sha256, so the set can be rebuilt or audited later.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/function_truth_2026_11_04/sources'

QUICKGO = 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch'
#: Experimental evidence codes we accept. ECO:0000269 is the root of "experimental evidence"; the
#: others are named explicitly in the request so the query records our intent even though several
#: are reached by `evidenceCodeUsage=descendants` anyway.
EXPERIMENTAL_ECO = ['ECO:0000269', 'ECO:0000314', 'ECO:0000353', 'ECO:0000315', 'ECO:0000316',
                    'ECO:0000270']
UNIPROT = ('https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:5811+AND+reviewed:true'
           '&fields=accession,id,protein_name,gene_names,cc_function,ft_domain,go_p,go_f,keyword,'
           'lit_pubmed_id&format=json&size=500')
#: Second pass over the same entries, asking only for the identifiers needed to reach a TGME49 gene
#: id: UniProt's own VEuPathDB/ToxoDB cross-references, the ordered-locus names and the strain taxon.
#: The `fields` list of the main query suppresses every cross-reference except GO, so the mapping
#: evidence has to be requested separately rather than guessed.
UNIPROT_XREF = ('https://rest.uniprot.org/uniprotkb/search?query=taxonomy_id:5811+AND+reviewed:true'
                '&fields=accession,id,gene_oln,gene_names,organism_id,organism_name,'
                'xref_veupathdb&format=json&size=500')
#: The GO ontology itself. QuickGO's TSV download carries only the GO id, so term names come from
#: here, and so does the structural specificity (depth, descendant count) used for information
#: content. `go-basic.obo` is the filtered release whose relations are all safe to traverse up.
GO_OBO = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
#: 32 MB of ontology does not belong in the repository, so it is pinned next to the other dataset
#: snapshots (where the 2026-10-31 QuickGO download already lives) and `sources/manifest.json`
#: records the path together with its sha256.
OBO_SNAPSHOT = Path('/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/sources/'
                    'go_basic_2026_11_04/go-basic.obo')
PAUSE = 0.5  # seconds between requests: <= 2/s, inside the <= 3/s budget we were given


def get(url: str, accept: str | None = None) -> tuple[bytes, dict]:
    """One polite GET. Returns the body and a provenance record (url, status, time, sha256)."""
    request = urllib.request.Request(url, headers={'User-Agent': 'ParaScape/function-truth (starplast)'})
    if accept:
        request.add_header('Accept', accept)
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    with urllib.request.urlopen(request, timeout=300) as response:
        body = response.read()
        record = dict(url=url, status=response.status, retrieved_utc=stamp,
                      sha256=hashlib.sha256(body).hexdigest(), bytes=len(body),
                      link=response.headers.get('Link', ''))
    time.sleep(PAUSE)
    return body, record


def quickgo(aspect: str, manifest: list) -> None:
    url = (f'{QUICKGO}?taxonId=5811&taxonUsage=descendants&aspect={aspect}'
           f'&evidenceCode={",".join(EXPERIMENTAL_ECO)}&evidenceCodeUsage=descendants'
           f'&downloadLimit=50000&includeFields=goName,taxonName')
    body, record = get(url, accept='text/tsv')
    name = f'quickgo_{aspect}_experimental.tsv'
    (OUT / name).write_bytes(body)
    record.update(file=name, source='quickgo', aspect=aspect,
                  rows=max(0, body.decode().count('\n') - 1))
    manifest.append(record)
    print(f'{name}: {record["rows"]} annotation rows')


def uniprot(manifest: list, url: str, stem: str, source: str) -> None:
    page = 0
    while url:
        body, record = get(url)
        page += 1
        name = f'{stem}_page{page}.json'
        (OUT / name).write_bytes(body)
        entries = len(json.loads(body).get('results', []))
        record.update(file=name, source=source, page=page, entries=entries)
        manifest.append(record)
        print(f'{name}: {entries} reviewed entries')
        link = record['link']
        url = link.split('<', 1)[1].split('>', 1)[0] if 'rel="next"' in link else ''


def main() -> None:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    manifest: list = []
    for aspect in ('biological_process', 'molecular_function'):
        quickgo(aspect, manifest)
    uniprot(manifest, UNIPROT, 'uniprot_reviewed_5811', 'uniprot_reviewed')
    uniprot(manifest, UNIPROT_XREF, 'uniprot_reviewed_5811_xref', 'uniprot_reviewed_xref')
    body, record = get(GO_OBO)
    OBO_SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
    OBO_SNAPSHOT.write_bytes(body)
    record.update(file=str(OBO_SNAPSHOT), source='gene_ontology')
    manifest.append(record)
    print(f'{OBO_SNAPSHOT}: {len(body)} bytes')
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(files=len(manifest), out=str(OUT))))


if __name__ == '__main__':
    main()
