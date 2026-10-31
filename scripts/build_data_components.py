"""Prepare deterministic bounded archive parts for independently checked source packs.

The builder copies existing verified artifacts without changing any measurements.
It does not upload or publish the provisional release index.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
import shutil
import sys
import tarfile
import time
import uuid

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parascape import data_components as C


def signature(path):
    stat = path.stat()
    return [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def save(path, value):
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


class _Parts:
    def __init__(self, folder, basename, limit):
        self.folder, self.basename, self.limit = folder, basename, limit
        self.parts = []
        self.archive_hash = hashlib.sha256()
        self.archive_size = 0
        self.stream = None
        self.part_hash = None
        self.part_size = 0

    def _close(self):
        if self.stream is not None:
            self.stream.close()
            self.parts.append(dict(file=self.stream.name.rsplit('/', 1)[-1], size=self.part_size,
                                   sha256=self.part_hash.hexdigest()))
            self.stream = None

    def write(self, data):
        original_size = len(data)
        self.archive_hash.update(data)
        self.archive_size += len(data)
        while data:
            if self.stream is None:
                self.stream = (self.folder / f'{self.basename}.part-{len(self.parts):04d}').open('wb')
                self.part_hash = hashlib.sha256()
                self.part_size = 0
            count = min(len(data), self.limit - self.part_size)
            self.stream.write(data[:count])
            self.part_hash.update(data[:count])
            self.part_size += count
            data = data[count:]
            if self.part_size == self.limit:
                self._close()
        return original_size


class _Source:
    def __init__(self, stream):
        self.stream, self.digest = stream, hashlib.sha256()

    def read(self, size):
        block = self.stream.read(size)
        self.digest.update(block)
        return block


def deposited_identity(files):
    """Address the complete quantitative pack, including its original receipts."""
    return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()


def deposited_files(source):
    """Require full source-cell validation before archiving any quantity payload."""
    from parascape import deposited_values as D
    source = Path(source)
    info = D.manifest(source)
    names = ['manifest.json', 'validation.json', *info['outputs']]
    for name in ('preparation.json', 'technical_qc.csv'):
        if (source / name).exists():
            names.append(name)
    files = {}
    for name in names:
        if Path(name).name != name or (source / name).is_symlink():
            raise ValueError('Deposited component files must be ordinary pack files')
        files[name] = dict(size=(source / name).stat().st_size,
                          sha256=info['outputs'].get(name) or C._hash(source / name))
    receipt = json.loads((source / 'validation.json').read_text())
    if ('technical_qc.csv' in files
            and receipt.get('technical_qc_sha256') != files['technical_qc.csv']['sha256']):
        raise ValueError('Deposited numerical QC checksum differs')
    if 'preparation.json' in files:
        preparation = json.loads((source / 'preparation.json').read_text())
        if preparation['manifest_sha256'] != files['manifest.json']['sha256']:
            raise ValueError('Deposited preparation manifest checksum differs')
    return files


def write_deposited_metadata(source, component, core, relative_directory):
    """Copy only original identity receipts into the core; payload stays optional."""
    from parascape.optional_deposited_data import MARKER
    source, core = Path(source), Path(core)
    C.validate(component)
    C._file_name(relative_directory)
    pieces = Path(relative_directory).parts
    if (len(pieces) not in (2, 3) or pieces[0] != 'deposited_values'
            or pieces[-1] != source.name or any(p.startswith('.') for p in pieces)):
        raise ValueError('Invalid deposited metadata directory')
    files = deposited_files(source)
    identity = deposited_identity(files)
    if component['directory'] != 'deposited_values/' + identity or component['files'] != files:
        raise ValueError('Deposited metadata and prepared component differ')
    target = core.joinpath(*pieces)
    for parent in [core, *[core.joinpath(*pieces[:n]) for n in range(1, len(pieces) + 1)]]:
        if parent.is_symlink():
            raise ValueError('Deposited metadata cannot use linked directories')
    target.mkdir(parents=True, exist_ok=True)
    if any(path.name not in {'manifest.json', 'validation.json', MARKER} for path in target.iterdir()):
        raise ValueError('Metadata-only destination already contains payload or unrelated files')
    marker = dict(schema_version=1, component_id=identity,
                  manifest_sha256=files['manifest.json']['sha256'],
                  validation_sha256=files['validation.json']['sha256'])
    payloads = {name: (source / name).read_bytes() for name in ('manifest.json', 'validation.json')}
    if any(hashlib.sha256(data).hexdigest() != files[name]['sha256'] for name, data in payloads.items()):
        raise ValueError('Deposited metadata changed after source validation')
    payloads[MARKER] = (json.dumps(marker, sort_keys=True, indent=2) + '\n').encode()
    for name, data in payloads.items():
        output = target / name
        if output.is_symlink() or (output.exists() and output.read_bytes() != data):
            raise ValueError('Existing deposited metadata differs; preserve the old generation')
        if not output.exists():
            with output.open('xb') as stream:
                stream.write(data)
    manifest = json.loads(payloads['manifest.json'])
    return dict(components=[component['directory']], source_url=manifest['source_url'],
                input_sha256=manifest['input_sha256'], manifest_sha256=marker['manifest_sha256'],
                validation_sha256=marker['validation_sha256'], metadata_directory=relative_directory)


def build_component(source, directory, destination, *, tag='v0.0.8.0', part_limit=1024**3):
    """Copy one checked pack into deterministic archive parts, retaining all bytes."""
    if type(part_limit) is not int or not 0 < part_limit <= 1024**3:
        raise ValueError('Archive part limit must be at most one GiB')
    source, destination = Path(source), Path(destination)
    kind, identity = directory.split('/')
    if kind not in C.KINDS or not re.fullmatch('[a-f0-9]{64}', identity) or not re.fullmatch('[A-Za-z0-9._-]+', tag) or tag in ('.', '..'):
        raise ValueError('Invalid component directory or release tag')
    manifest = json.loads((source / 'manifest.json').read_text())
    if kind == 'deposited_values':
        files = deposited_files(source)
        if identity != deposited_identity(files):
            raise ValueError('Deposited component identity differs from its original files')
    elif kind == 'workbook_cells':
        from parascape import workbook_values as W
        W.manifest(source, verify=True)
        receipt = W.validation(source, manifest)
        if receipt['status'] not in ('all_OOXML_source_cells_independently_reconciled',
                                      'BIFF_checksums_verified_independent_reader_review_required',
                                      'all_BIFF_physical_cells_and_decoded_values_reconciled'):
            raise ValueError('Workbook component requires a bound verification record with explicit scope')
        files = {name: dict(size=(source/name).stat().st_size, sha256=C._hash(source/name))
                 for name in ('manifest.json', 'independent_validation.json', *manifest['outputs'])}
    else:
        output = 'axis.parquet' if kind == 'sparse_axes' else 'entries.parquet'
        receipt = json.loads((source / 'validation.json').read_text())
        expected_state = 'all_original_axis_cells_reconciled' if kind == 'sparse_axes' else 'all_original_coordinate_records_reconciled'
        if kind == 'sparse_coordinates' and manifest.get('source_format') in ('HDF5_CSC', 'AnnData_HDF5'):
            from parascape.sparse_values import coordinate_reconciliation
            coordinate_reconciliation(manifest, receipt)
            expected_state = ('all_original_HDF5_CSC_records_reconciled' if manifest['source_format'] == 'HDF5_CSC'
                              else 'all_original_AnnData_layer_records_reconciled')
        if (receipt['state'] != expected_state or receipt['manifest_sha256'] != C._hash(source / 'manifest.json')
                or receipt['input_sha256'] != manifest['input_sha256'] or set(manifest['outputs']) != {output}):
            raise ValueError('Original component lacks matching independent validation')
        files = {name: dict(size=(source/name).stat().st_size, sha256=(manifest['outputs'][output] if name == output else C._hash(source/name)))
                 for name in ('manifest.json', 'validation.json', output)}
    for name in files:
        if (source / name).is_symlink():
            raise ValueError('Source artifacts must be regular pack files')
    signatures = {name: signature(source/name) for name in files}
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / (kind + '-' + identity)
    bindings = dict(directory=directory, files=files, tag=tag, part_limit=part_limit)
    if target.is_symlink():
        raise ValueError('Prepared component cannot be a symbolic link')
    if target.exists():
        cached = json.loads((target/'build.json').read_text())
        component = json.loads((target/'component.json').read_text())
        if (cached['bindings'] != bindings or component['directory'] != directory
                or component['files'] != files or len(component['parts']) != len(cached['parts'])):
            raise ValueError('Existing component source bindings differ')
        for number, part in enumerate(cached['parts']):
            if part['file'] != f'parascape-{kind}-{identity}.tar.part-{number:04d}':
                raise ValueError('Unsafe prepared component part name')
            remote = component['parts'][number]
            if not remote['url'].endswith('/' + part['file']) or remote['size'] != part['size'] or remote['sha256'] != part['sha256']:
                raise ValueError('Prepared component part description differs')
            path = target/part['file']
            if path.is_symlink() or (signature(path) != cached['part_signatures'][part['file']] and C._hash(path) != part['sha256']):
                raise ValueError('Prepared component part changed')
        if signatures != cached['source_signatures']:
            for name, item in files.items():
                if C._hash(source/name) != item['sha256']:
                    raise ValueError('Original source component changed')
        C.validate(component)
        return dict(directory=directory, component=component, path=str(target), reused=True)
    staging = destination / ('.building-' + kind + '-' + identity + '-' + uuid.uuid4().hex[:8])
    staging.mkdir()
    writer = _Parts(staging, 'parascape-' + kind + '-' + identity + '.tar', part_limit)
    builder_sha256 = C._hash(Path(__file__))
    try:
        with tarfile.open(fileobj=writer, mode='w|', format=tarfile.USTAR_FORMAT) as bundle:
            for name, item in sorted(files.items()):
                info = tarfile.TarInfo(name)
                info.size, info.mode, info.mtime = item['size'], 0o644, 0
                with (source/name).open('rb') as stream:
                    checked = _Source(stream)
                    bundle.addfile(info, checked)
                    if checked.digest.hexdigest() != item['sha256']:
                        raise ValueError('Source artifact changed while archiving')
        writer._close()
        if signatures != {name: signature(source/name) for name in files} or C._hash(Path(__file__)) != builder_sha256:
            raise ValueError('Source component or builder code changed during preparation')
        component = dict(directory=directory, files=files, archive_size=writer.archive_size,
                         archive_sha256=writer.archive_hash.hexdigest(),
                         parts=[dict(url=f'https://github.com/EinarOlafsson/parascape/releases/download/{tag}/{part["file"]}',
                                     size=part['size'], sha256=part['sha256']) for part in writer.parts])
        C.validate(component)
        save(staging/'component.json', component)
        save(staging/'build.json', dict(bindings=bindings, epoch=time.time(), builder_sha256=builder_sha256,
             source_signatures=signatures, parts=writer.parts,
             part_signatures={part['file']:signature(staging/part['file']) for part in writer.parts},
             measurement_processing='No transformation; existing independently reconciled files retained byte-for-byte'))
        staging.replace(target)
        return dict(directory=directory, component=component, path=str(target), reused=False)
    finally:
        writer._close()
        if staging.exists():
            shutil.rmtree(staging)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('directory')
    parser.add_argument('destination')
    parser.add_argument('--tag', default='v0.0.8.0')
    args = parser.parse_args()
    print(json.dumps(build_component(args.source, args.directory, args.destination, tag=args.tag)))
