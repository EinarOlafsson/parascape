"""Host salmon later (Einar 2026-11-13, instruction 69).

Salmon against the combined ME49_host index maps some projects about 4.6 times slower than others, and one
deep run held a worker for more than 7 hours. The cheap stages of every queued run now go first; salmon against
ME49_host follows in a later pass. Settings and outputs are unchanged; only the order is. done.json is written
only when every stage is done, so nothing that reads it sees a half-finished run as finished.
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import omics_quant as Q  # noqa: E402

RUN = {"run_accession": "SRR0000001", "study": "GSE1", "references": ["ME49", "ME49_host"]}


def _fake_tools(tmp_path, monkeypatch, calls, reads_after=90):
    monkeypatch.setattr(Q, "NAS", str(tmp_path / "nas"))
    monkeypatch.setattr(Q, "SCRATCH", str(tmp_path / "scratch"))
    monkeypatch.setattr(Q, "DOWNLOAD_CACHE", "")
    monkeypatch.setattr(Q, "NAS_SCRATCH", "")
    monkeypatch.setattr(Q, "DRAIN_FILE", str(tmp_path / "DRAIN"))
    monkeypatch.setattr(Q, "HOST_SALMON_LATER", {"ME49_host"})
    monkeypatch.setattr(Q, "HOST_SALMON_LATER_STUDIES", set())
    monkeypatch.setattr(Q, "log", lambda msg: None)
    monkeypatch.setattr(Q, "version", lambda tool: f"{tool} 1")
    monkeypatch.setattr(Q, "_still_planned", lambda acc: True)
    monkeypatch.setattr(Q, "_claim", lambda acc: True)
    monkeypatch.setattr(Q, "_release", lambda acc: None)
    os.makedirs(Q.SCRATCH, exist_ok=True)
    monkeypatch.setattr(Q, "_reads", lambda r: [f"{r['run_accession']}_1.fq", f"{r['run_accession']}_2.fq"])
    trim = {"reads_before": 100, "reads_after": reads_after}
    monkeypatch.setattr(Q, "_trim", lambda acc, local: (local, dict(trim)))

    def run(cmd, **kw):
        tool = os.path.basename(cmd[0])
        out = Path(cmd[cmd.index("-o") + 1])
        ref = Path(cmd[cmd.index("-i") + 1]).parent.name
        calls.append((tool, ref, out.parent.name))
        (out / "aux_info").mkdir(parents=True, exist_ok=True)
        if tool == "salmon":
            (out / "quant.sf").write_text("Name\tTPM\n")
            (out / "aux_info" / "meta_info.json").write_text(json.dumps({"num_mapped": 5}))
        else:
            (out / "abundance.tsv").write_text("target_id\ttpm\n")
    monkeypatch.setattr(Q, "run", run)


def test_cheap_pass_leaves_host_salmon_pending_and_writes_no_done_json(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    assert Q._quantify_pass([RUN], None, 0, Q.CHEAP) == []
    assert calls == [("salmon", "ME49", "SRR0000001"), ("kallisto", "ME49", "SRR0000001"),
                     ("kallisto", "ME49_host", "SRR0000001")]
    run_dir = tmp_path / "nas" / "rnaseq" / "GSE1" / "SRR0000001"
    assert not (run_dir / "done.json").exists()
    staged = json.loads((run_dir / Q.STAGES_FILE).read_text())
    assert staged["stages_pending"] == ["salmon_ME49_host"]
    assert staged["stages_done"] == ["archive", "reads", "fastp", "salmon_ME49", "kallisto_ME49", "kallisto_ME49_host"]
    assert staged["trimming"]["reads_after"] == 90
    calls.clear()
    assert Q._quantify_pass([RUN], None, 0, Q.CHEAP) == []     # cheap stages are not repeated
    assert calls == []


def test_host_pass_finishes_the_run_with_the_same_record(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    Q._quantify_pass([RUN], None, 0, Q.CHEAP)
    calls.clear()
    Q._quantify_pass([RUN], None, 0, Q.HOST)
    assert calls == [("salmon", "ME49_host", "SRR0000001")]
    run_dir = tmp_path / "nas" / "rnaseq" / "GSE1" / "SRR0000001"
    done = json.loads((run_dir / "done.json").read_text())
    assert not (run_dir / Q.STAGES_FILE).exists()
    assert done["stages_pending"] == [] and "salmon_ME49_host" in done["stages_done"]
    assert done["finished"] and done["salmon"] and done["kallisto"]          # what queue verification needs
    assert all(done[key] == value for key, value in RUN.items())             # the plan row, unchanged
    assert (run_dir / "salmon_ME49_host" / "quant.sf.gz").exists()
    assert (run_dir / "kallisto_ME49_host" / "abundance.tsv.gz").exists()


def test_host_pass_refuses_reads_that_differ_from_the_first_pass(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    Q._quantify_pass([RUN], None, 0, Q.CHEAP)
    _fake_tools(tmp_path, monkeypatch, calls, reads_after=89)
    calls.clear()
    Q._quantify_pass([RUN], None, 0, Q.HOST)
    run_dir = tmp_path / "nas" / "rnaseq" / "GSE1" / "SRR0000001"
    assert calls == [] and not (run_dir / "done.json").exists() and (run_dir / Q.STAGES_FILE).exists()


def test_every_run_gets_its_cheap_stages_before_any_host_salmon(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    runs = [dict(RUN, run_accession=f"SRR000000{i}") for i in (1, 2)]
    Q._run_phases(runs, None)
    order = [(tool, ref) for tool, ref, _ in calls]
    assert order.index(("salmon", "ME49_host")) > max(i for i, c in enumerate(order) if c != ("salmon", "ME49_host"))
    assert all((tmp_path / "nas" / "rnaseq" / "GSE1" / r["run_accession"] / "done.json").exists() for r in runs)


def test_a_run_without_a_deferred_reference_finishes_in_the_first_pass(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    cow = dict(RUN, references=["ME49", "ME49_cow"])
    Q._quantify_pass([cow], None, 0, Q.CHEAP)
    assert (tmp_path / "nas" / "rnaseq" / "GSE1" / "SRR0000001" / "done.json").exists()
    assert ("salmon", "ME49_cow", "SRR0000001") in calls


def test_the_later_pass_can_be_narrowed_to_named_studies(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    monkeypatch.setattr(Q, "HOST_SALMON_LATER_STUDIES", {"GSE210671"})
    Q._quantify_pass([RUN], None, 0, Q.CHEAP)                 # GSE1 is not named: every stage at once
    assert (tmp_path / "nas" / "rnaseq" / "GSE1" / "SRR0000001" / "done.json").exists()
    assert ("salmon", "ME49_host", "SRR0000001") in calls


def test_one_pass_when_nothing_is_deferred(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    monkeypatch.setattr(Q, "HOST_SALMON_LATER", set())
    Q._run_phases([RUN], None)
    assert [c[:2] for c in calls] == [("salmon", "ME49"), ("kallisto", "ME49"),
                                      ("salmon", "ME49_host"), ("kallisto", "ME49_host")]


def test_a_drain_marker_stops_new_runs(tmp_path, monkeypatch):
    calls = []
    _fake_tools(tmp_path, monkeypatch, calls)
    Path(Q.DRAIN_FILE).write_text("")
    assert Q._quantify_pass([RUN], None, 0, Q.CHEAP) is None
    assert calls == []


def test_the_raw_cache_fetches_host_pending_runs_last(tmp_path, monkeypatch):
    import nas_overflow_cache as C
    monkeypatch.setattr(C.Q, "NAS", str(tmp_path))
    staged = {"study": "GSE1", "run_accession": "SRR1"}
    fresh = {"study": "GSE1", "run_accession": "SRR2"}
    (tmp_path / "rnaseq" / "GSE1" / "SRR1").mkdir(parents=True)
    (tmp_path / "rnaseq" / "GSE1" / "SRR1" / Q.STAGES_FILE).write_text("{}")
    assert sorted([staged, fresh], key=C.host_salmon_pending) == [fresh, staged]
