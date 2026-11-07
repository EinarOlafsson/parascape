"""The organism space catalogue (instruction 64.23): explicit states, never shown as calibrated unless it is."""
from __future__ import annotations

import os
import shutil

import pytest

from parascape import paths, space_catalogue as SC

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def empty_cache(tmp_path, monkeypatch):
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path / "state"))
    return tmp_path


def _by_code(entries):
    return {e.code: e for e in entries}


def test_every_declared_space_has_one_explicit_state(empty_cache):
    entries = SC.catalogue()
    assert {"Tg", "Pf", "Hs", "Mm"} <= {e.code for e in entries}
    assert all(e.state in ("ready", "installed", "installable", "unavailable") for e in entries)
    assert all(SC.describe(e) for e in entries)


def test_host_packs_install_with_their_checksum_and_are_never_calibrated(empty_cache):
    archives = SC.local_archives("Mm")
    if not archives:
        pytest.skip("no mouse pack in this checkout")
    assert _by_code(SC.catalogue())["Mm"].state == "installable"
    SC.install("Mm")
    mm = _by_code(SC.catalogue())["Mm"]
    assert mm.state == "installed" and mm.version and not mm.calibrated
    assert "not calibrated" in SC.describe(mm)


def test_a_pack_without_its_published_checksum_is_refused(empty_cache):
    archives = SC.local_archives("Mm")
    if not archives:
        pytest.skip("no mouse pack in this checkout")
    copy = empty_cache / "elsewhere" / os.path.basename(archives[0][0])
    copy.parent.mkdir()
    shutil.copy(archives[0][0], copy)
    with pytest.raises(ValueError, match="published SHA256"):
        SC.install("Mm", str(copy))
    with pytest.raises(ValueError, match="mismatch"):
        SC.install("Mm", str(copy), "0" * 64)
    assert _by_code(SC.catalogue(roots=[str(empty_cache / "nothing")]))["Mm"].state == "unavailable"
