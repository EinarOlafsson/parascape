"""Uniform re-search of one PRIDE project's RAW spectra, keeping salmon-like outputs without the RAW.

Proteomics is the RNA pipeline's counterpart: every project is processed with the same tools, protein
database and FDR, and what is kept lets it be revisited without the instrument files (Einar 2026-11-03):
centroided MS1+MS2 spectra as gzipped indexed mzML, the full Sage PSM table (every peptide-spectrum
match with scores and q-values), the label-free peptide quantities, the search parameters, the
ThermoRawFileParser metadata and checksums. The RAW files are deleted only after all of that is written.

Fragment tolerance follows the MS2 analyser read from the converted spectra (Orbitrap 20 ppm, ion trap
0.5 Da); runs with different analysers are searched separately. Label-free quantification only: TMT
and DIA projects are refused here and left for their own settings.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import urllib.request

BIN = Path.home() / 'anaconda3/envs/parascape-ms/bin'
NAS = Path('/nas_mnt/carruthers/Einar/parascape_omics')
PLAN = NAS / 'published_pride/acquisition_plan.json'
FASTA = NAS / 'references/ms_fasta/parascape_ms_ME49_human_mouse_contaminants.fasta'
SCRATCH = NAS / 'ms_raw_scratch'
OUT = NAS / 'prepared_ms_uniform'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 22), b''):
            digest.update(block)
    return digest.hexdigest()


def raw_files(project: str) -> list[dict]:
    plan = json.loads(PLAN.read_text())
    return [f for f in plan['all_source_files']
            if "RAW" in str(f['categories']) and str(f['studies']).strip() == project]


def fetch(url: str, dest: Path, size: int) -> None:
    part = dest.with_suffix(dest.suffix + '.part')
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as response, part.open('wb') as out:
                shutil.copyfileobj(response, out, 1 << 22)
            if size and part.stat().st_size != size:
                raise IOError(f'{dest.name}: {part.stat().st_size} bytes, declared {size}')
            part.replace(dest)
            return
        except Exception:
            if attempt == 4:
                raise
            time.sleep(30 * (attempt + 1))


def acquisition(mzml: Path) -> dict:
    """Scan types seen in the filter strings: MS3 scans suggest TMT-SPS, MS2 without data-dependent markers DIA."""
    ms2 = ms3 = dependent = 0
    with gzip.open(mzml, 'rt', errors='replace') as stream:
        for line in stream:
            if 'filter string' not in line:
                continue
            if ' ms3 ' in line:
                ms3 += 1
            elif ' ms2 ' in line:
                ms2 += 1
                dependent += ' d ' in line
            if ms2 + ms3 >= 2000:
                break
    return dict(ms2_scans_checked=ms2, ms3_scans_checked=ms3, data_dependent_ms2=dependent,
                flag='possible TMT (MS3 scans)' if ms3 else 'possible DIA (no data-dependent MS2)' if ms2 and not dependent else 'DDA')


def ms2_analyser(mzml: Path) -> str:
    """'FTMS' or 'ITMS', from the filter strings of the first MS2 spectra."""
    seen, checked = set(), 0
    with gzip.open(mzml, 'rt', errors='replace') as stream:
        for line in stream:
            if 'filter string' in line and ' ms2 ' in line:
                seen.add('ITMS' if 'ITMS' in line else 'FTMS' if 'FTMS' in line else '?')
                checked += 1
                if checked >= 200:                     # enough scans to see a mixed acquisition
                    break
    if seen == {'FTMS'} or seen == {'ITMS'}:
        return seen.pop()
    raise ValueError(f'{mzml.name}: MS2 analyser not determined ({seen})')


def sage_config(mzmls: list[Path], analyser: str, out: Path) -> dict:
    fragment = {'ppm': [-20.0, 20.0]} if analyser == 'FTMS' else {'da': [-0.5, 0.5]}
    return {
        'database': {'fasta': str(FASTA), 'generate_decoys': True, 'decoy_tag': 'rev_',
                     'enzyme': {'missed_cleavages': 2, 'min_len': 7, 'max_len': 50, 'cleave_at': 'KR', 'restrict': 'P'},
                     'static_mods': {'C': 57.021464},
                     'variable_mods': {'M': [15.994915], '[': [42.010565]}, 'max_variable_mods': 2},
        'precursor_tol': {'ppm': [-10.0, 10.0]}, 'fragment_tol': fragment, 'isotope_errors': [0, 2],
        'deisotope': analyser == 'FTMS', 'min_peaks': 15, 'max_peaks': 150,
        'quant': {'lfq': True}, 'predict_rt': True, 'wide_window': False,
        'output_directory': str(out), 'mzml_paths': [str(p) for p in mzmls],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('project')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    project = args.project
    out = OUT / project
    if (out / 'done.json').exists():
        print(json.dumps(dict(project=project, status='already done'))); return
    files = raw_files(project)
    if not files:
        raise ValueError(f'{project}: no RAW files in the PRIDE plan')
    scratch = SCRATCH / project
    shutil.rmtree(scratch, ignore_errors=True); scratch.mkdir(parents=True)
    spectra = out / 'spectra'; spectra.mkdir(parents=True, exist_ok=True)
    started = time.time(); receipts = []
    # PRIDE serves one file at roughly 1 MB/s from here; several transfers at once fill the line.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda f: fetch(f['url'], scratch / f['filename'], int(f.get('declared_bytes') or 0)), files))
    for f in files:
        raw = scratch / f['filename']
        raw_sha = sha256(raw)
        # Plain mzML + gzip (-f 1 -g): with an indexed mzML (-f 2) and -g, ThermoRawFileParser 2.0.0.dev reopens an
        # uncompressed file to write the index and fails; an index has no use inside a gzip stream anyway.
        subprocess.run([str(BIN / 'thermorawfileparser'), '-i', str(raw), '-o', str(spectra), '-f', '1', '-g',
                        '-m', '0', '-l', '1'], check=True, stdout=subprocess.DEVNULL)
        mzml = spectra / (raw.stem + '.mzML.gz')
        if not mzml.exists():
            raise RuntimeError(f'{raw.name}: no mzML written')
        receipts.append(dict(raw=f['filename'], url=f['url'], raw_bytes=raw.stat().st_size, raw_sha256=raw_sha,
                             mzml=mzml.name, mzml_bytes=mzml.stat().st_size, mzml_sha256=sha256(mzml),
                             metadata=f'{raw.stem}-metadata.json',
                             metadata_sha256=sha256(spectra / f'{raw.stem}-metadata.json'),
                             ms2_analyser=ms2_analyser(mzml), acquisition=acquisition(mzml)))
        raw.unlink()                                   # spectra and checksums kept; the instrument file is not
    groups = {}
    for r in receipts:
        groups.setdefault(r['ms2_analyser'], []).append(spectra / r['mzml'])
    searches = []
    for analyser, mzmls in sorted(groups.items()):
        target = out / f'sage_{analyser}'
        target.mkdir(exist_ok=True)
        config = sage_config(mzmls, analyser, target)
        (target / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
        subprocess.run([str(BIN / 'sage'), '--batch-size', '4', str(target / 'config.json')], check=True,
                       env={'RAYON_NUM_THREADS': str(args.threads), 'PATH': '/usr/bin:/bin'})
        outputs = {p.name: dict(bytes=p.stat().st_size, sha256=sha256(p)) for p in target.iterdir() if p.is_file()}
        searches.append(dict(analyser=analyser, runs=[p.name for p in mzmls], outputs=outputs))
    shutil.rmtree(scratch, ignore_errors=True)
    record = dict(project=project, runs=receipts, searches=searches, fasta=str(FASTA), fasta_sha256=sha256(FASTA),
                  tools=dict(thermorawfileparser=subprocess.run([str(BIN / 'thermorawfileparser'), '--version'],
                                                                capture_output=True, text=True).stdout.strip(),
                             sage=subprocess.run([str(BIN / 'sage'), '--version'], capture_output=True, text=True).stdout.strip()),
                  kept='centroided MS1+MS2 mzML.gz per run (all scans); instrument metadata JSON; Sage PSM table, LFQ table, parameters; checksums',
                  raw_kept=False, seconds=round(time.time() - started),
                  source_sha256=sha256(Path(__file__)))
    (out / 'done.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(dict(project=project, runs=len(receipts), searches=[s['analyser'] for s in searches],
                          raw_gb=round(sum(r['raw_bytes'] for r in receipts) / 1e9, 2),
                          mzml_gb=round(sum(r['mzml_bytes'] for r in receipts) / 1e9, 2), seconds=record['seconds'])))


if __name__ == '__main__':
    main()
