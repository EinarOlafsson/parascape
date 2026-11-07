"""Your own kept inferences (instruction 78): kept, shown as kept, never facts, never training data,
superseded loudly when measured."""
import pathlib
import re

import pandas as pd
import pytest

from parascape import gene_reading as GR
from parascape import kept_claims as KC


@pytest.fixture(autouse=True)
def store(tmp_path, monkeypatch):
    monkeypatch.setenv("PARASCAPE_KEPT_CLAIMS", str(tmp_path / "kept.jsonl"))
    monkeypatch.setattr(KC, "release_id", lambda: "0.0.10.0")
    return tmp_path / "kept.jsonl"


class Map:
    def __init__(self, measured=None):
        self.measured = measured

    def gene_facts(self, gene):
        rows = [dict(column="compartment_hyperlopit", tier=1, table="information map",
                     description="niche", observed=True, value="dense granules",
                     probability=0.9, source_tier="observed-curated")]
        if self.measured is not None:
            rows.append(dict(column="host_facing", tier=1, table="information map",
                             description="faces the host", observed=True, value=self.measured,
                             probability=0.95, source_tier="observed-direct"))
        return pd.DataFrame(rows)

    def gene_evidence(self, gene):
        return pd.DataFrame(columns=["gene_id", "trait", "value", "probability", "source_tier",
                                     "reference", "dataset"])

    def gene_partners(self, gene):
        return pd.DataFrame(columns=["partner", "side", "probability", "source_tier", "reference"])


def _claim(**kw):
    base = dict(gene="TGME49_000001", trait="host_facing", value="yes", probability=0.83,
                method="map neighbours", method_version="k=10", accepted_score=0.83)
    base.update(kw)
    return KC.KeptClaim(**base)


def test_a_claim_is_kept_with_its_provenance():
    c = KC.keep(_claim())
    got = KC.claims("TGME49_000001")
    assert len(got) == 1 and got.iloc[0].release == "0.0.10.0" and got.iloc[0].kept_at
    assert not c.cleared_bar, "0.83 is below the 90 % bar"
    assert KC.keep(_claim(accepted_score=0.95)).cleared_bar


def test_a_claim_without_its_method_is_refused():
    with pytest.raises(ValueError):
        KC.keep(_claim(method=""))


def test_a_kept_claim_never_appears_among_the_facts():
    KC.keep(_claim())
    reading = GR.read(Map(), "TGME49_000001")
    every_fact = pd.concat(list(reading.facts.values()))
    assert "host_facing" not in set(every_fact["column"])
    assert "yes" not in " ".join(s for _c, s in GR.synthesis(reading))


def test_nothing_that_builds_training_data_reads_the_store():
    """Structural: only the screens that show kept claims import the module."""
    root = pathlib.Path(KC.__file__).parent
    importers = {p.name for p in root.glob("*.py")
                 if re.search(r"\bkept_claims\b", p.read_text()) and p.name != "kept_claims.py"}
    assert importers <= {"gene_module.py", "infer_module.py"}, importers


def test_a_claim_survives_a_release_upgrade_with_its_provenance(monkeypatch):
    KC.keep(_claim())
    monkeypatch.setattr(KC, "release_id", lambda: "0.0.11.0")
    got = KC.claims("TGME49_000001")
    assert got.iloc[0].release == "0.0.10.0", "the claim says which release it was made against"


def test_a_measured_cell_supersedes_the_claim_loudly():
    KC.keep(_claim(value="yes"))
    agree = KC.check_against(Map(measured=True))
    assert len(agree) == 1 and bool(agree.iloc[0].agrees)
    disagree = KC.check_against(Map(measured=False))
    assert not bool(disagree.iloc[0].agrees)
    assert KC.check_against(Map()).empty, "unmeasured: nothing to supersede"


def test_forget_one_or_every_claim_of_a_method_version():
    a = KC.keep(_claim())
    KC.keep(_claim(gene="TGME49_000002"))
    KC.keep(_claim(gene="TGME49_000003", method="infer:host_facing", method_version="q1"))
    assert KC.revert(a.claim_id)
    assert len(KC.claims()) == 2
    assert KC.revert_method("map neighbours", "k=10") == 1
    assert list(KC.claims().method) == ["infer:host_facing"]
    assert len(KC.claims(include_reverted=True)) == 3, "reverted claims stay in the history"


def test_export_writes_every_claim_with_provenance(tmp_path):
    KC.keep(_claim())
    out = tmp_path / "kept.csv"
    assert KC.export(str(out)) == 1
    frame = pd.read_csv(out)
    assert {"release", "method", "method_version", "question", "accepted_score"} <= set(frame.columns)
