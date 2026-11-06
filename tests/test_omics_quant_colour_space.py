"""The quantifier and colour-space reads (instruction 89).

GSE60305 was sequenced on an AB SOLiD 4, whose reads are colour space ("T1330222..."). kallisto
failed on every run, 49 times in a day, and salmon "succeeded" with 0 of 47 million reads
mapped -- so the study was absent and the cause invisible. Colour-space runs are now recorded
before any download and never retried, and a quantifier that maps nothing is a failure.
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import omics_quant as Q  # noqa: E402


def test_solid_runs_are_colour_space_and_illumina_runs_are_not():
    assert Q.colour_space({"instrument_model": "AB SOLiD 4 System"})
    assert Q.colour_space({"instrument_model": "AB 5500xl Genetic Analyzer SOLiD"})
    assert not Q.colour_space({"instrument_model": "Illumina HiSeq 2500"})
    assert not Q.colour_space({})


def test_salmon_mapped_reads_are_read_from_its_meta_info(tmp_path):
    (tmp_path / "aux_info").mkdir()
    (tmp_path / "aux_info" / "meta_info.json").write_text(json.dumps({"num_mapped": 1234}))
    assert Q.salmon_mapped(str(tmp_path)) == 1234
    assert Q.salmon_mapped(str(tmp_path / "missing")) is None


def test_a_colour_space_run_is_recorded_without_downloading(tmp_path, monkeypatch):
    monkeypatch.setattr(Q, "NAS", str(tmp_path))
    monkeypatch.setattr(Q, "_reads", lambda r: (_ for _ in ()).throw(AssertionError("downloaded")))
    monkeypatch.setattr(Q, "_still_planned", lambda acc: True)
    monkeypatch.setattr(Q, "_claim", lambda acc: True)
    monkeypatch.setattr(Q, "log", lambda msg: None)
    run = {"run_accession": "SRR1542920", "study": "GSE60305", "references": ["ME49"],
           "instrument_model": "AB SOLiD 4 System"}
    Q._quantify_pass([run], None, 0)
    record = json.loads((tmp_path / "rnaseq" / "GSE60305" / "SRR1542920" /
                         Q.COLOUR_SPACE_FILE).read_text())
    assert "colour-space" in record["reason"]
    assert not (tmp_path / "rnaseq" / "GSE60305" / "SRR1542920" / "done.json").exists()
    Q._quantify_pass([run], None, 0)                  # recorded runs are skipped, not retried
