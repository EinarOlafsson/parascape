"""Check built package metadata, cache files, icons, and size before uploading."""
from __future__ import annotations

import argparse
from email.parser import BytesParser
from pathlib import Path
import zipfile


def check_wheel(path: Path) -> None:
    """Reject incomplete wheels, accidental saved analyses, and oversized uploads."""
    if path.stat().st_size > 100_000_000:
        raise ValueError(f"{path.name} exceeds the default PyPI 100 MB file limit")
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        metadata = next(n for n in names if n.endswith(".dist-info/METADATA"))
        info = BytesParser().parsebytes(archive.read(metadata))
        if info["Name"] != "parascape":
            raise ValueError(f"Unexpected distribution: {info['Name']}")
        required = {'parascape/app.py', 'parascape/paths.py', 'parascape/data_manager.py',
                    'parascape/data_download_dialog.py', 'parascape/data_manifest.json',
                    'parascape/rna_profiles.py', 'parascape/rna_panel.py',
                    'parascape/rna_design.py', 'parascape/rna_design_panel.py',
                    'parascape/data/icons/parascape.svg', 'parascape/spaces/__init__.py'}
        missing = required - names
        if missing:
            raise ValueError(f"Incomplete wheel: {sorted(missing)}")
        import json
        descriptor = json.loads(archive.read('parascape/data_manifest.json'))
        if not descriptor.get('sha256') or not descriptor.get('files') or descriptor.get('version') != info['Version']:
            raise ValueError('Missing or mismatched GitHub data release descriptor')
        if not descriptor['url'].startswith('https://github.com/EinarOlafsson/parascape/releases/download/'):
            raise ValueError('Data must be hosted on the official GitHub release')
        if any(n.startswith('parascape/data/') and not n.startswith('parascape/data/icons/') for n in names):
            raise ValueError('Scientific data must be downloaded, not bundled in the wheel')
        for requirement in info.get_all("Requires-Dist", []):
            if requirement.lower().startswith("parascape"):
                raise ValueError(f"Unexpected ParaScape dependency: {requirement}")
            if requirement.lower().startswith(("cuml", "cupy")) and 'extra == "gpu"' not in requirement:
                raise ValueError(f"CUDA must be optional: {requirement}")
    print(f"Checked {path.name}: {path.stat().st_size / 1_000_000:.1f} MB")


def main():
    """Validate each wheel in a distribution directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=Path("dist"))
    directory = parser.parse_args().directory
    wheels = list(directory.glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError(f"Expected one wheel, found {len(wheels)}")
    expected_sdist = wheels[0].name.split("-py", 1)[0] + ".tar.gz"
    if {p.name for p in directory.iterdir()} != {wheels[0].name, expected_sdist}:
        raise ValueError("Upload directory must contain only the matching parascape wheel and sdist")
    for wheel in wheels:
        check_wheel(wheel)


if __name__ == "__main__":
    main()
