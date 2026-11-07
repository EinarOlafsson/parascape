"""Re-quantifying the earlier runs (instruction 81): beside the verified result, compared with it,
the runs that mapped nothing first, and only after the new downloads."""
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import omics_quant as Q  # noqa: E402


def _salmon(folder: Path, tpm: dict, mapped: int) -> None:
    folder.mkdir(parents=True)
    rows = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n" + "".join(
        f"{k}\t1000\t900\t{v}\t{v}\n" for k, v in tpm.items())
    with gzip.open(folder / "quant.sf.gz", "wt") as fh:
        fh.write(rows)
    (folder / "aux_info").mkdir()
    (folder / "aux_info" / "meta_info.json").write_text(json.dumps({"num_mapped": mapped}))


def test_a_close_requantification_agrees(tmp_path, monkeypatch):
    monkeypatch.setattr(Q, "_cache_module", lambda: type("C", (), {"save": staticmethod(
        lambda path, data: Path(path).write_text(json.dumps(data)))}))
    tpm = {f"t{i}": float(i + 2) for i in range(40)}
    _salmon(tmp_path / "salmon_ME49", tpm, 1000)
    _salmon(tmp_path / "requant" / "salmon_ME49", {k: v * 1.05 for k, v in tpm.items()}, 1010)
    verdict = Q.compare_requant(str(tmp_path), str(tmp_path / "requant"), ["ME49"])
    assert verdict["verdict"] == "agrees" and verdict["spearman"] > 0.99
    assert (tmp_path / "requant" / "compare.json").exists()


def test_a_run_that_mapped_nothing_is_recovered_not_compared(tmp_path, monkeypatch):
    monkeypatch.setattr(Q, "_cache_module", lambda: type("C", (), {"save": staticmethod(
        lambda path, data: Path(path).write_text(json.dumps(data)))}))
    _salmon(tmp_path / "salmon_ME49", {f"t{i}": 0.0 for i in range(40)}, 0)
    _salmon(tmp_path / "requant" / "salmon_ME49", {f"t{i}": float(i) for i in range(40)}, 5000)
    assert Q.compare_requant(str(tmp_path), str(tmp_path / "requant"),
                             ["ME49"])["verdict"].startswith("recovered")


def test_untrimmed_runs_are_queued_zero_mapped_first(tmp_path, monkeypatch):
    plan = tmp_path / "plan.jsonl"
    rows = [{"run_accession": "A", "study": "S", "references": ["ME49"]},
            {"run_accession": "B", "study": "S", "references": ["ME49"]},
            {"run_accession": "C", "study": "S", "references": ["ME49"]}]
    plan.write_text("".join(json.dumps(r) + "\n" for r in rows))
    monkeypatch.setattr(Q, "PLAN", str(plan))
    monkeypatch.setattr(Q, "NAS", str(tmp_path))
    for acc, mapped, trimmed in (("A", 900, False), ("B", 0, False), ("C", 900, True)):
        d = tmp_path / "rnaseq" / "S" / acc
        _salmon(d / "salmon_ME49", {"t": 5.0}, mapped)
        (d / "done.json").write_text(json.dumps({"trimming": {}} if trimmed else {}))
    assert [r["run_accession"] for r in Q.requant_todo()] == ["B", "A"]
    assert Q.pending_new_runs() == 0
