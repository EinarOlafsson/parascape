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
        required = {"parascape/app.py", "parascape/data/nodes.parquet", "parascape/data/graph.npz",
                    "parascape/data/toxodb_identity.tsv", "parascape/data/pf_nodes.parquet",
                    "parascape/data/pf_graph.npz", "parascape/data/icons/parascape.svg",
                    "parascape/data/icons/Apicomplexa_cells.svg", "parascape/data/icons/Animal_cells.svg",
                    "parascape/data/toxodb_cds.tsv.gz", "parascape/data/plasmodb_cds.tsv.gz",
                    "parascape/workflows.py", "parascape/prediction.py", "parascape/evidence.py",
                    "parascape/data/af3_features.parquet", "parascape/data/af3_manifest.parquet",
                    "parascape/data/af3_features.json", "parascape/data/esm_features.parquet",
                    "parascape/data/esm_manifest.parquet", "parascape/data/pf_mentions.parquet",
                    "parascape/data/hs_host_proteins.parquet",
                    "parascape/data/mm_host_proteins.parquet",
                    "parascape/data/host_species_migration.json", "parascape/packs.py",
                    "parascape/spaces/__init__.py"}
        missing = required - names
        if missing:
            raise ValueError(f"Incomplete wheel: {sorted(missing)}")
        if "parascape/data/host_proteins.parquet" in names:
            raise ValueError("The wheel contains the legacy mixed-species host table; clean the build")
        if any(n.startswith("parascape/data/embeddings/") for n in names):
            raise ValueError("The wheel contains local saved embeddings")
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
