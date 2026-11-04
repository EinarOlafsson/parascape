# Host-facing truth set, reviewed by Einar (2026-11-04)

Built by `scripts/build_host_facing_truth.py` from his decisions on every one of the 140 candidate claims
(https://claude.ai/artifact/TbAj2n9cnBPxZ7Leb5duCh). 140 claims over 98 genes.

* **66 genes host-facing**, **32 not**.
* Rule: strongest evidence wins (direct > proximity > inferred). 3 genes had
  disagreeing claims; GRA60, GRA44 and GRA47 are explicit overrides to "not host-facing".
* Probability: a gene decided by a direct experimental claim is `observed-direct` 0.95, one decided by a
  proximity or inferred claim is `observed-curated` 0.90.
* What makes this set valuable: 32 negatives are **vacuole-resident** proteins, the hard
  kind. The earlier map had 9 such negatives and several hundred markers of compartments inside the parasite,
  which answer an easier question. Einar rejected 22 dense granule proteins that the PVM proximity screen had
  called positive, which is the distinction the first model failed to learn.

## Files
* `host_facing_truth.tsv` -- one row per gene: decision, the claim that decided it, tier, probability.
* `reviewed_claims.tsv` -- all 140 reviewed claims as recorded.
* `evidence.parquet` -- evidence rows for the information map.
