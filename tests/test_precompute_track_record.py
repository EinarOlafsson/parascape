"""Resumable track-record partitions (instruction 64.25): reuse, invalidation, failure, cancellation, budget."""
from __future__ import annotations

import pandas as pd
import pytest

from parascape import precompute_track_record as P


def _parts(settings_b='{"k": 1}', inputs="n1"):
    return [P.Partition("Tg", "label", "a", inputs, '{"k": 1}', "c1"),
            P.Partition("Tg", "label", "b", inputs, settings_b, "c1"),
            P.Partition("Tg", "label", "c", inputs, '{"k": 1}', "c1")]


class Builder:
    def __init__(self, fail=()):
        self.calls, self.fail = [], set(fail)

    def __call__(self, p):
        self.calls.append(p.strategy)
        if p.strategy in self.fail:
            raise RuntimeError("boom")
        return pd.DataFrame({"strategy": [p.strategy], "gene_id": ["g1"], "correct": [True]})


def test_interrupted_run_resumes_only_unfinished_partitions(tmp_path):
    b = Builder()
    first = P.run(_parts(), tmp_path, P.Budget(max_jobs=1), builder=b)
    assert list(first.values()) == ["built", "pending", "pending"]
    second = P.run(_parts(), tmp_path, builder=b)
    assert second == {p.name: s for p, s in zip(_parts(), ["reused", "built", "built"])}
    assert b.calls == ["a", "b", "c"]
    assert P.run(_parts(), tmp_path, builder=b) and b.calls == ["a", "b", "c"]      # full replay builds nothing
    assert len(P.assemble(tmp_path, _parts())) == 3


def test_a_changed_setting_invalidates_only_that_method_and_a_changed_table_everything(tmp_path):
    b = Builder()
    P.run(_parts(), tmp_path, builder=b)
    b.calls.clear()
    P.run(_parts(settings_b='{"k": 2}'), tmp_path, builder=b)
    assert b.calls == ["b"]
    b.calls.clear()
    P.run(_parts(inputs="n2"), tmp_path, builder=b)
    assert b.calls == ["a", "b", "c"]


def test_a_failure_is_recorded_and_the_other_partitions_continue(tmp_path):
    states = P.run(_parts(), tmp_path, builder=Builder(fail={"b"}))
    assert list(states.values()) == ["built", "failed", "built"]
    report = P.report(tmp_path).set_index("partition")
    assert "boom" in report.loc[_parts()[1].name, "error"]
    assert len(P.assemble(tmp_path, _parts())) == 2          # a failed partition never enters the table
    assert P.run(_parts(), tmp_path, builder=Builder())[_parts()[1].name] == "built"


def test_cancellation_keeps_earlier_partitions_and_edited_bytes_are_rebuilt(tmp_path):
    calls = {"n": 0}

    def cancel():
        calls["n"] += 1
        return calls["n"] > 1
    states = P.run(_parts(), tmp_path, builder=Builder(), cancelled=cancel)
    assert list(states.values()) == ["built", "pending", "pending"]
    path = tmp_path / "partitions" / f"{_parts()[0].name}.parquet"
    path.write_bytes(path.read_bytes() + b"x")
    assert not P.valid(tmp_path, _parts()[0])
    assert P.run(_parts(), tmp_path, builder=Builder())[_parts()[0].name] == "built"


def test_a_corrupt_journal_is_refused(tmp_path):
    P.run(_parts(), tmp_path, builder=Builder())
    journal = tmp_path / "journal.json"
    journal.write_text(journal.read_text().replace('"ok"', '"OK"', 1))
    with pytest.raises(ValueError, match="corrupt"):
        P.run(_parts(), tmp_path, builder=Builder())
