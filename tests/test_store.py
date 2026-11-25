"""Packed record folders (instruction 102): one Parquet per unread audit folder, exact bytes kept."""
from __future__ import annotations

import hashlib
import io
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parascape import store  # noqa: E402


def _tree(root: Path) -> Path:
    f = root / "audit" / "records"
    (f / "sub").mkdir(parents=True)
    (f / "a.json").write_bytes(b'{"b": 1,  "a": 2}\n')           # odd spacing must survive
    (f / "sub" / "c.tsv").write_bytes(b"x\ty\n1\t2\n")
    return f


def test_a_pack_restores_its_folder_byte_for_byte(tmp_path):
    folder = _tree(tmp_path)
    pack = tmp_path / "records.pack.parquet"
    pack.write_bytes(store.pack_bytes(folder))
    assert store.unpack(pack, tmp_path / "out") == 2
    for p in folder.rglob("*"):
        if p.is_file():
            assert (tmp_path / "out" / p.relative_to(folder)).read_bytes() == p.read_bytes()
    assert store.file_bytes(pack, "a.json") == b'{"b": 1,  "a": 2}\n'


def test_packing_is_deterministic(tmp_path):
    """The release descriptor is verified by rebuilding, so the same folder must give the same bytes."""
    folder = _tree(tmp_path)
    assert store.pack_bytes(folder) == store.pack_bytes(folder)


def test_a_damaged_pack_is_refused(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    t = pa.table({"path": ["a.json"], "data": [b"changed"], "sha256": [hashlib.sha256(b"orig").hexdigest()]})
    pack = tmp_path / "bad.pack.parquet"
    pq.write_table(t, pack)
    with pytest.raises(ValueError, match="do not match their digest"):
        store.read(pack)


def test_an_unsafe_path_is_never_written(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    t = pa.table({"path": ["../escape"], "data": [b"x"], "sha256": [hashlib.sha256(b"x").hexdigest()]})
    pack = tmp_path / "evil.pack.parquet"
    pq.write_table(t, pack)
    with pytest.raises(ValueError, match="unsafe path"):
        store.unpack(pack, tmp_path / "out")


def test_only_folders_no_module_names_are_unread(tmp_path):
    audit = tmp_path / "audit"
    for name in ("used_folder", "unused_folder", "used_folder_v2"):
        (audit / name).mkdir(parents=True)
    code = tmp_path / "code"
    code.mkdir()
    (code / "m.py").write_text("ROOT = 'processing_audit/used_folder'\n")
    assert [p.name for p in store.unread_folders(audit, code)] == ["unused_folder", "used_folder_v2"]


def test_the_application_still_names_every_folder_it_reads():
    """A reader module must name its folder literally, or the builder packs what the app reads."""
    code = Path(__file__).resolve().parents[1] / "parascape"
    text = "\n".join(p.read_text(encoding="utf-8") for p in code.rglob("*.py"))
    for folder in ("pride_projects", "all_geo_samples", "protein_column_design_links",
                   "reference_identifiers"):
        assert folder in text, folder


def test_a_file_reads_the_same_from_disk_and_from_its_pack(tmp_path):
    folder = _tree(tmp_path)
    expected = (folder / "sub" / "c.tsv").read_bytes()
    (folder.parent / ("records" + store.SUFFIX)).write_bytes(store.pack_bytes(folder))
    assert store.read_bytes(folder / "sub" / "c.tsv") == expected
    import shutil
    shutil.rmtree(folder)
    store._PACKS.clear()
    assert store.read_bytes(folder / "sub" / "c.tsv") == expected
    with pytest.raises(FileNotFoundError):
        store.read_bytes(folder / "absent.json")
    assert not store.exists(tmp_path / "nowhere" / "at_all.json")   # walks up to / without failing


def test_the_coverage_panel_reads_a_packed_geo_folder(tmp_path):
    """all_geo_samples ships packed; the panel's checksum-verified reader must not notice."""
    import json
    import shutil
    from parascape import dataset_coverage_panel as P
    root = tmp_path / "processing_audit"
    geo = root / "all_geo_samples"
    report = json.dumps({"accession": "GSE1", "source_samples": 1, "samples": [{"gsm": "GSM1"}]}).encode()
    (geo / "GSE1").mkdir(parents=True)
    (geo / "GSE1" / "geo-samples.json").write_bytes(report)
    (geo / "coverage.json").write_text(json.dumps([{"accession": "GSE1", "status": "validated_sample_metadata",
                                                    "sample_metadata_sha256": hashlib.sha256(report).hexdigest()}]))
    loose = P.source_samples(root, "GSE1")
    (root / ("all_geo_samples" + store.SUFFIX)).write_bytes(store.pack_bytes(geo))
    shutil.rmtree(geo)
    store._PACKS.clear()
    assert P.source_samples(root, "GSE1") == loose


def test_every_reader_of_a_packed_read_folder_goes_through_the_store():
    """A folder in READ_THROUGH_STORE ships packed, so any module that names it and still opens files
    with Path.read_bytes / read_text would find nothing in a release."""
    code = Path(__file__).resolve().parents[1] / "parascape"
    naming = sorted(p.name for p in code.rglob("*.py") if p.name != "store.py"
                    and any(folder in p.read_text(encoding="utf-8") for folder in store.READ_THROUGH_STORE))
    assert naming == ["dataset_coverage_panel.py", "source_design_links.py"], \
        f"move these readers to the store: {naming}"
