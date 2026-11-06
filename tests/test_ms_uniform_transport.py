"""PRIDE downloads (instruction 88): the faster transport per project, files in parallel, stalls restarted."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
ms = pytest.importorskip("ms_uniform_search")


def test_the_faster_transport_is_kept(monkeypatch):
    speeds = {"ftp": 0.14, "https": 0.71}
    monkeypatch.setattr(ms, "_probe", lambda url: speeds["ftp" if url.startswith("ftp") else "https"])
    assert ms.fastest_transport("https://ftp.pride.ebi.ac.uk/pride/data/x.raw") == "https"
    speeds.update(ftp=9.2, https=1.0)
    assert ms.fastest_transport("https://ftp.pride.ebi.ac.uk/pride/data/x.raw") == "ftp"


def test_when_both_probes_fail_ftp_is_kept(monkeypatch):
    monkeypatch.setattr(ms, "_probe", lambda url: 0.0)
    assert ms.fastest_transport("https://ftp.pride.ebi.ac.uk/pride/data/x.raw") == "ftp"


def test_parallel_files_stay_within_what_pride_tolerates():
    assert ms.FILES_AT_ONCE * 6 <= 12
