"""One export bundle for a query (64.39): labelled tables, methods text, a manifest a reload verifies."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from parascape import datasets as D
from parascape import export_bundle as E


def _column():
    d = next(d for d in D.REGISTRY if d.citation and d.columns and d.organism == "Tg" and not d.derived_from)
    return d.columns[0]


def test_a_bundle_round_trips_and_labels_every_table(tmp_path):
    answer = pd.DataFrame({"gene_id": ["g1", "g2"], "label": ["x", "y"]})
    path = E.write(tmp_path, organism="Tg", query={"question": "q?"}, tables={"answer": answer},
                   kinds={"answer": "inferred"}, columns=[_column(), "unregistered_thing"],
                   strategy="map neighbourhoods")
    manifest = E.load(tmp_path)
    assert manifest["query"] == {"question": "q?"}
    assert manifest["files"]["answer.csv"]["kind"] == "inferred"
    assert manifest["files"]["methods.md"]["unregistered_columns"] == ["unregistered_thing"]
    assert (tmp_path / "methods.md").read_text().startswith("Analyses were run in ParaScape")
    assert json.loads(path.read_text())["files"]["answer.csv"]["rows"] == 2


def test_a_table_not_marked_measured_or_inferred_is_refused(tmp_path):
    with pytest.raises(ValueError, match="measured, inferred or computed"):
        E.write(tmp_path, organism="Tg", query={}, tables={"t": pd.DataFrame({"a": [1]})}, kinds={})


def test_an_edited_bundle_is_caught_on_reload(tmp_path):
    E.write(tmp_path, organism="Tg", query={}, tables={"t": pd.DataFrame({"a": [1]})},
            kinds={"t": "measured"}, columns=[])
    (tmp_path / "t.csv").write_text("a\n2\n")
    with pytest.raises(ValueError, match="changed or missing"):
        E.load(tmp_path)
