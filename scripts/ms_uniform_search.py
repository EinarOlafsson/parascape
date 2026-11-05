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
EXTRA = NAS / 'processing_audit/pride_ms_uniform_extra_files.json'
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
    records = json.loads(PLAN.read_text())['all_source_files']
    if EXTRA.exists():                                  # projects found after the acquisition plan was built
        records = records + json.loads(EXTRA.read_text())['files']
    return [f for f in records
            if "RAW" in str(f['categories']) and str(f['studies']).strip() == project]


#: Waits before each retry of a download. PRIDE refuses a connection outright when it has had enough, and a
#: refusal that is retried in 30 s simply comes back; two projects were lost that way on 2026-11-04. The last
#: attempts use a single plain connection, which is slow but is what the server is most willing to serve.
RETRY_WAITS = (30, 120, 300, 600, 1200)


def fetch(url: str, dest: Path, size: int) -> None:
    part = dest.with_suffix(dest.suffix + '.part')
    for attempt in range(len(RETRY_WAITS) + 1):
        try:
            if (BIN / 'aria2c').exists() and attempt < 3:
                # One PRIDE connection gives ~1 MB/s; a few ranged connections give several times that. More
                # than about a dozen at once, counting every file in flight, earns a refusal.
                subprocess.run([str(BIN / 'aria2c'), '-q', '-x6', '-s6', '-k8M', '--file-allocation=none',
                                '--allow-overwrite=true', '--auto-file-renaming=false', '--connect-timeout=60',
                                '--max-tries=1', '-d', str(part.parent), '-o', part.name, url], check=True)
            else:
                with urllib.request.urlopen(url, timeout=600) as response, part.open('wb') as out:
                    shutil.copyfileobj(response, out, 1 << 22)
            if size and part.stat().st_size != size:
                raise IOError(f'{dest.name}: {part.stat().st_size} bytes, declared {size}')
            part.replace(dest)
            return
        except Exception:
            part.unlink(missing_ok=True)                  # a partial file would fail the size check next time
            if attempt == len(RETRY_WAITS):
                raise
            time.sleep(RETRY_WAITS[attempt])


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


def file_kind(name: str) -> str:
    """'raw' (Thermo), 'mzml', 'bruker' (zipped timsTOF .d), 'archive' (zip/tar/rar/7z) or 'unsupported' (SCIEX)."""
    n = name.lower()
    if n.endswith('.raw'):
        return 'raw'
    if n.endswith(('.mzml', '.mzml.gz')):
        return 'mzml'
    if n.endswith('.d.zip'):
        return 'bruker'
    if n.endswith(('.raw.zip', '.zip', '.tar.gz', '.tgz', '.tar', '.rar', '.7z')):
        return 'archive'
    return 'unsupported'


def unpack(archive: Path, target: Path) -> str:
    """Unpack into target; returns '' or a note when a deposited tar is truncated (complete members are kept)."""
    import tarfile
    import zipfile
    target.mkdir(parents=True, exist_ok=True)
    root = target.resolve()
    if archive.name.lower().endswith(('.rar', '.7z')):
        subprocess.run([str(BIN / '7zz'), 'x', '-y', f'-o{target}', str(archive)], check=True, stdout=subprocess.DEVNULL)
        for path in target.rglob('*'):                 # 7-Zip strips unsafe paths; confirm nothing landed outside
            if not path.resolve().is_relative_to(root):
                raise ValueError(f'{archive.name}: unsafe member {path}')
    elif archive.name.lower().endswith('.zip'):
        with zipfile.ZipFile(archive) as z:
            for member in z.namelist():                # refuse paths that would leave the folder
                if not (root / member).resolve().is_relative_to(root):
                    raise ValueError(f'{archive.name}: unsafe member {member}')
            z.extractall(target)
    else:
        member = None
        try:
            with tarfile.open(archive, 'r|*') as t:    # streaming: members before a truncation stay usable
                for member in t:
                    t.extract(member, target, filter='data')   # refuses links and paths outside
                    member = None
        except (EOFError, tarfile.ReadError, OSError) as error:
            if member is not None and member.isfile():
                (target / member.name).unlink(missing_ok=True)  # the member cut off mid-way
            return f'{archive.name} is truncated as deposited ({error}); members before the cut kept'
    return ''


def bruker(d_folder: Path, spectra: Path, record: dict, keep_zip: Path | None = None) -> dict:
    """A timsTOF .d folder: stored zipped as the run's spectra (4D data has no compact mzML form), searched as is."""
    if keep_zip is None:
        keep_zip = Path(shutil.make_archive(str(spectra / d_folder.name), 'zip', d_folder.parent, d_folder.name))
    return dict(record, conversion='Bruker timsTOF .d kept zipped as deposited; Sage reads it directly (MS2)',
                stored=keep_zip.name, stored_bytes=keep_zip.stat().st_size, stored_sha256=sha256(keep_zip),
                search_path=str(d_folder), ms2_analyser='TIMS', acquisition=dict(flag='timsTOF DDA-PASEF or DIA (not read here)'))


def convert(source: Path, spectra: Path, record: dict) -> dict:
    """Centroided mzML.gz for one Thermo RAW (deleted afterwards), or an author mzML kept as deposited."""
    if file_kind(source.name) == 'raw':
        # Plain mzML + gzip (-f 1 -g): with an indexed mzML (-f 2) and -g, ThermoRawFileParser 2.0.0.dev reopens an
        # uncompressed file to write the index and fails; an index has no use inside a gzip stream anyway.
        subprocess.run([str(BIN / 'thermorawfileparser'), '-i', str(source), '-o', str(spectra), '-f', '1', '-g',
                        '-m', '0', '-l', '1'], check=True, stdout=subprocess.DEVNULL)
        mzml = spectra / (source.stem + '.mzML.gz')
        if not mzml.exists():
            raise RuntimeError(f'{source.name}: no mzML written')
        meta = spectra / f'{source.stem}-metadata.json'
        record = dict(record, conversion='ThermoRawFileParser centroided MS1+MS2', metadata=meta.name,
                      metadata_sha256=sha256(meta) if meta.exists() else None)
    else:
        stem = source.name[:-len('.gz')] if source.name.lower().endswith('.gz') else source.name
        mzml = spectra / (stem[:-len('.mzml')] + '.mzML.gz')
        if source.name.lower().endswith('.gz'):
            shutil.copyfile(source, mzml)
        else:
            with source.open('rb') as src, gzip.open(mzml, 'wb') as dst:
                shutil.copyfileobj(src, dst, 1 << 22)
        record = dict(record, conversion='author-deposited mzML kept as deposited (peak picking as submitted)')
    try:
        analyser = ms2_analyser(mzml)
    except ValueError:
        analyser = None
    source.unlink()                                    # spectra and checksums kept; the instrument file is not
    return dict(record, mzml=mzml.name, mzml_bytes=mzml.stat().st_size, mzml_sha256=sha256(mzml),
                ms2_analyser=analyser, acquisition=acquisition(mzml))


def sage_config(mzmls: list[Path], analyser: str, out: Path) -> dict:
    fragment = {'da': [-0.5, 0.5]} if analyser == 'ITMS' else {'ppm': [-20.0, 20.0]}
    return {
        'database': {'fasta': str(FASTA), 'generate_decoys': True, 'decoy_tag': 'rev_',
                     'enzyme': {'missed_cleavages': 2, 'min_len': 7, 'max_len': 50, 'cleave_at': 'KR', 'restrict': 'P'},
                     'static_mods': {'C': 57.021464},
                     'variable_mods': {'M': [15.994915], '[': [42.010565]}, 'max_variable_mods': 2},
        'precursor_tol': {'ppm': [-15.0, 15.0]} if analyser == 'TIMS' else {'ppm': [-10.0, 10.0]},
        'fragment_tol': fragment, 'isotope_errors': [0, 2],
        'deisotope': analyser == 'FTMS', 'min_peaks': 15, 'max_peaks': 150,
        # timsrust in this Sage reads MS2 only, so timsTOF runs get identifications without MS1 label-free quantities.
        'quant': {'lfq': analyser != 'TIMS'}, 'predict_rt': True, 'wide_window': False,
        'output_directory': str(out), 'mzml_paths': [str(p) for p in mzmls],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('project')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    project = args.project
    out = OUT / project
    if (out / 'done.json').exists() or (out / 'unsupported.json').exists():
        print(json.dumps(dict(project=project, status='already handled'))); return
    files = raw_files(project)
    if not files:
        raise ValueError(f'{project}: no RAW files in the PRIDE plan')
    scratch = SCRATCH / project
    shutil.rmtree(scratch, ignore_errors=True); scratch.mkdir(parents=True)
    spectra = out / 'spectra'; spectra.mkdir(parents=True, exist_ok=True)
    started = time.time(); receipts = []; unsupported = []; working = []
    usable = [f for f in files if file_kind(f['filename']) != 'unsupported']
    unsupported += [dict(file=f['filename'], bytes=int(f.get('declared_bytes') or 0),
                         reason='SCIEX wiff or other format without a Linux converter (ProteoWizard needs Windows, Wine or Docker)')
                    for f in files if file_kind(f['filename']) == 'unsupported']
    # One file at a time over six ranged connections: PRIDE refused two whole projects on 2026-11-04 when two
    # files were in flight at six connections each, and a single file at six measured as fast anyway.
    for f in usable:
        fetch(f['url'], scratch / f['filename'], int(f.get('declared_bytes') or 0))
    for f in usable:
        source = scratch / f['filename']
        base = dict(source=f['filename'], url=f['url'], source_bytes=source.stat().st_size, source_sha256=sha256(source))
        if file_kind(f['filename']) == 'bruker':
            kept = spectra / source.name
            shutil.move(str(source), kept)              # the zipped .d is the stored spectra
            unpacked = scratch / (source.name + '_unpacked'); unpack(kept, unpacked)
            folders = sorted(p for p in unpacked.rglob('*.d') if (p / 'analysis.tdf').exists())
            if not folders:                             # e.g. older Bruker .baf data: kept, not searched
                unsupported.append(dict(file=f['filename'], bytes=base['source_bytes'],
                                        reason='Bruker .d without analysis.tdf (not timsTOF); zip kept, not searched'))
            for d_folder in folders:
                receipts.append(bruker(d_folder, spectra, base, keep_zip=kept))
            working.append(unpacked)
        elif file_kind(f['filename']) == 'archive':
            unpacked = scratch / (source.name + '_unpacked')
            note = unpack(source, unpacked)
            if note:
                unsupported.append(dict(file=f['filename'], bytes=base['source_bytes'], reason=note))
            source.unlink()
            members = sorted(p for p in unpacked.rglob('*') if p.is_file() and file_kind(p.name) in ('raw', 'mzml'))
            folders = sorted(p for p in unpacked.rglob('*.d') if (p / 'analysis.tdf').exists())
            if not members and not folders:
                unsupported.append(dict(file=f['filename'], bytes=base['source_bytes'],
                                        reason='archive holds no Thermo RAW, mzML or Bruker .d'))
            for member in members:
                receipts.append(convert(member, spectra, dict(base, archive_member=str(member.relative_to(unpacked)))))
            for d_folder in folders:
                receipts.append(bruker(d_folder, spectra, dict(base, archive_member=str(d_folder.relative_to(unpacked)))))
            if folders:
                working.append(unpacked)                # searched in place, removed after the search
            else:
                shutil.rmtree(unpacked, ignore_errors=True)
        else:
            try:
                receipts.append(convert(source, spectra, base))
            except subprocess.CalledProcessError as error:   # one unreadable deposit must not sink the project
                unsupported.append(dict(file=f['filename'], bytes=base['source_bytes'],
                                        reason=f'conversion failed (exit {error.returncode}); RAW as deposited not readable'))
                for partial in spectra.glob(source.stem + '*'):
                    partial.unlink()
                source.unlink(missing_ok=True)
    if not receipts:
        (out / 'unsupported.json').write_text(json.dumps(dict(project=project, unsupported=unsupported,
            note='No spectra this pipeline can convert; kept as a recorded gap.'), indent=2) + '\n')
        shutil.rmtree(scratch, ignore_errors=True)
        print(json.dumps(dict(project=project, status='unsupported formats only', files=len(unsupported)))); return
    groups = {}
    for r in receipts:
        if r['ms2_analyser'] is None:
            unsupported.append(dict(file=r['mzml'], reason='MS2 analyser not determined from filter strings; spectra kept, not searched'))
            continue
        groups.setdefault(r['ms2_analyser'], []).append(Path(r['search_path']) if 'search_path' in r else spectra / r['mzml'])
    searches = []
    for analyser, mzmls in sorted(groups.items()):
        target = out / f'sage_{analyser}'
        target.mkdir(exist_ok=True)
        config = sage_config(mzmls, analyser, target)
        (target / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
        # Sage holds a whole batch of runs' spectra in memory; very large runs go one at a time.
        largest = max((p.stat().st_size if p.is_file() else 0) for p in mzmls)
        batch = '1' if largest > 400e6 else '4'
        subprocess.run([str(BIN / 'sage'), '--batch-size', batch, str(target / 'config.json')], check=True,
                       env={'RAYON_NUM_THREADS': str(args.threads), 'PATH': '/usr/bin:/bin'})
        outputs = {p.name: dict(bytes=p.stat().st_size, sha256=sha256(p)) for p in target.iterdir() if p.is_file()}
        searches.append(dict(analyser=analyser, runs=[p.name for p in mzmls], outputs=outputs))
    for folder in working:
        shutil.rmtree(folder, ignore_errors=True)
    shutil.rmtree(scratch, ignore_errors=True)
    record = dict(project=project, runs=receipts, searches=searches, unsupported=unsupported,
                  fasta=str(FASTA), fasta_sha256=sha256(FASTA),
                  tools=dict(thermorawfileparser=subprocess.run([str(BIN / 'thermorawfileparser'), '--version'],
                                                                capture_output=True, text=True).stdout.strip(),
                             sage=subprocess.run([str(BIN / 'sage'), '--version'], capture_output=True, text=True).stdout.strip()),
                  kept='centroided MS1+MS2 mzML.gz per run (all scans), or the zipped timsTOF .d; instrument metadata JSON; Sage PSM table, LFQ table, parameters; checksums',
                  raw_kept=False, seconds=round(time.time() - started),
                  source_sha256=sha256(Path(__file__)))
    (out / 'done.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(dict(project=project, runs=len(receipts), searches=[s['analyser'] for s in searches],
                          raw_gb=round(sum({r['source']: r['source_bytes'] for r in receipts}.values()) / 1e9, 2),
                          mzml_gb=round(sum(r.get('mzml_bytes', r.get('stored_bytes', 0)) for r in receipts) / 1e9, 2), seconds=record['seconds'])))


if __name__ == '__main__':
    main()
