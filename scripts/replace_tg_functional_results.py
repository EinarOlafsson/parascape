#!/usr/bin/env python3
"""Replace the Toxoplasma frozen functional results after the node table changed (instruction 64.11a).

The GSE108740 log fix changed the Toxoplasma node table, so the two frozen Toxoplasma recoveries (FN-EC-01 EC
profiles, FN-DOM Pfam profiles), bound to the old node bytes, are refused by `functional_results.require_context`.
Einar (2026-11-07) chose to re-run them on a fresh split: EC seed 24 (was 23), Pfam seed 20261107 (was 20261008).
This installs the two new packaged bundles in place of the old entries, keeps every Plasmodium entry exactly, and
moves the reader pin. Each step is checked; nothing is fitted here.

    python scripts/replace_tg_functional_results.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BUNDLES = (ROOT / "results/functional_ec_ui_bundle_2026_11_07", ROOT / "results/functional_profile_ui_bundle_2026_11_07")


def main() -> int:
    """Swap the Toxoplasma entries and re-pin; print the before/after summary."""
    from parascape import functional_results as F, organisms as O
    shipped = ROOT / "parascape/data/functional_results.json"
    reader = Path(F.__file__)
    old_pin = F.FUNCTIONAL_RESULTS_SHA256
    F.load(shipped, expected_sha256=old_pin)
    document = json.loads(shipped.read_bytes())
    incoming = []
    for bundle in BUNDLES:
        receipt = json.loads((bundle / "receipt.json").read_text())
        path = bundle / "functional_results.json"
        (entry,) = F.load(path, expected_sha256=receipt["sha256"])
        assert entry.organism == O.TOXOPLASMA
        incoming.extend(json.loads(path.read_bytes())["benchmarks"])
    new_ids = {b["metadata"]["benchmark_id"] for b in incoming}
    kept = [b for b in document["benchmarks"] if b["metadata"]["organism"] != O.TOXOPLASMA]
    dropped = [b["metadata"]["benchmark_id"] for b in document["benchmarks"] if b["metadata"]["organism"] == O.TOXOPLASMA]
    assert set(dropped) == new_ids, (dropped, new_ids)
    document["benchmarks"] = incoming + kept
    merged = (F._canonical(document) + "\n").encode()
    digest = hashlib.sha256(merged).hexdigest()
    shipped.write_bytes(merged)
    code = reader.read_text()
    assert code.count(old_pin) == 1
    reader.write_text(code.replace(old_pin, digest))
    loaded = F.load(shipped, expected_sha256=digest)
    summary = {"old_sha256": old_pin, "new_sha256": digest, "replaced": sorted(dropped),
               "entries": [(b.organism, b.metadata["benchmark_id"], len(b.rows)) for b in loaded],
               "plasmodium_entries_kept_exactly": [json.dumps(b, sort_keys=True) for b in kept]
               == [json.dumps(b, sort_keys=True) for b in json.loads(merged)["benchmarks"][len(incoming):]]}
    print(json.dumps({k: v for k, v in summary.items()}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
