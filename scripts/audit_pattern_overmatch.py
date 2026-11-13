"""Instruction 65.04: slot patterns that select more columns than they were meant to.

Both biological eligibility defects found on 2026-11-13 are the same mistake wearing two coats.

  * `Tg_transcription_bradyzoite_tissue_cyst` carries the pattern `rna108740_Day`, which selects 14
    samples the submitter declares as cat enteroepithelial stages.
  * `Tg_protein_abundance_tachyzoite` carries the pattern `proteome_`, which selects the untreated
    tachyzoite columns AND nine columns of AP2XII-1-depleted parasites being reprogrammed to
    pre-sexual stages AND three log2 fold-changes.

A slot pattern is a PREFIX. It was written with particular columns in mind, and it keeps matching
whatever is added to the deposit later, so a slot silently acquires columns nobody chose for it.
Neither defect is visible in the slot declaration, in the registry, or in any test -- the pattern
looks right and the columns it names were never listed.

This finds them by looking at what each pattern ACTUALLY selects and asking whether those columns
describe one condition. It is offline and cheap, so it can run beside the suite.

**It catches only the half that the column names give away, and the other half is the reason
`resolve_stage_from_deposit.py` exists.** The `proteome_` defect is visible here because the columns
say `UT` against `T_24h` and `abundance` against `lfc`. The GSE108740 defect is NOT: those columns
are named `rna108740_Day3_CS4_FPKM`, and "Day" is a timepoint, not a stage. Nothing in that name
says cat, or enteroepithelial, or anything a word list could match -- only the submitter knew, and
only asking the deposit found it. The two scripts are complementary and neither replaces the other.

Writes ``results/dataset_selection_2026_11_08/pattern_overmatch.tsv``.
"""
from __future__ import annotations

import os
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import paths, slots as S  # noqa: E402

OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")

# Tokens that mark two columns as DIFFERENT conditions. A pattern selecting columns from more than
# one group on the same axis is selecting across a contrast, which is what a slot is meant to hold
# constant.
AXES = {
    "treatment": ([r"(^|_)ut(_|$)", r"untreated", r"vehicle", r"parent", r"\bwt\b", r"control"],
                  [r"(^|_)t_\d+h", r"treated", r"(^|_)ko(_|$)", r"knockdown", r"depleted",
                   r"actinomycin"]),
    "quantification": ([r"abundance", r"fpkm", r"tpm", r"ibaq", r"counts", r"cpm"],
                       [r"lfc", r"log2_fold", r"fold_change", r"_padj", r"_fdr", r"differential"]),
    "life stage": ([r"tachy"], [r"brady", r"cyst", r"oocyst", r"sporo", r"merozoite", r"day\d"]),
}


def which(column: str, axis: str) -> str:
    """Which side of a contrast a column sits on, or '' when it does not say."""
    left, right = AXES[axis]
    low = column.lower()
    on_left = any(re.search(p, low) for p in left)
    on_right = any(re.search(p, low) for p in right)
    if on_left and not on_right:
        return "a"
    if on_right and not on_left:
        return "b"
    return ""


def main() -> int:
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    columns = list(nodes.columns)
    rows = []
    for slot in S.all_slots():
        if slot.role in ("never", "metadata"):
            continue
        for pattern in slot.patterns:
            hit = [c for c in columns if c.startswith(pattern)]
            if len(hit) < 2:
                continue
            for axis in AXES:
                sides = {c: which(c, axis) for c in hit}
                a = sorted(c for c, s in sides.items() if s == "a")
                b = sorted(c for c, s in sides.items() if s == "b")
                if a and b:
                    rows.append({
                        "slot": slot.key, "axis": slot.axis, "policy": slot.policy,
                        "pattern": pattern, "contrast": axis,
                        "columns_selected": len(hit),
                        "on_one_side": len(a), "on_the_other": len(b),
                        "example_one": a[0], "example_other": b[0],
                        "note": (f"the pattern '{pattern}' selects {len(hit)} columns that differ "
                                 f"by {axis}; a slot is meant to hold that constant"),
                    })

    # The same question one level up. `Tg_transcription_bradyzoite_tissue_cyst` holds TWO patterns,
    # `rna108740_Tissue_cysts` and `rna108740_Day`, and each is internally consistent -- the
    # contrast is between them. Checking patterns one at a time cannot see it, so the columns of
    # the whole slot are checked together as well.
    for slot in S.all_slots():
        if slot.role in ("never", "metadata") or len(slot.patterns) < 2:
            continue
        hit = sorted({c for p in slot.patterns for c in columns if c.startswith(p)})
        if len(hit) < 2:
            continue
        for axis in AXES:
            sides = {c: which(c, axis) for c in hit}
            a = sorted(c for c, s in sides.items() if s == "a")
            b = sorted(c for c, s in sides.items() if s == "b")
            if a and b and not any(r["slot"] == slot.key and r["contrast"] == axis for r in rows):
                rows.append({
                    "slot": slot.key, "axis": slot.axis, "policy": slot.policy,
                    "pattern": "(across the patterns of the slot)", "contrast": axis,
                    "columns_selected": len(hit), "on_one_side": len(a), "on_the_other": len(b),
                    "example_one": a[0], "example_other": b[0],
                    "note": (f"the patterns of this slot together select {len(hit)} columns that "
                             f"differ by {axis}; each pattern is consistent, the slot is not"),
                })

    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "pattern_overmatch.tsv"), sep="\t", index=False)
    print(f"patterns selecting across a contrast: {len(frame)}")
    if len(frame):
        print()
        print(frame.groupby(["contrast"]).size().to_string())
        print()
        for r in frame.itertuples():
            print(f"  {r.slot}")
            print(f"    pattern '{r.pattern}' selects {r.columns_selected} columns, "
                  f"{r.on_one_side} vs {r.on_the_other} across {r.contrast}")
            print(f"    e.g. {r.example_one}")
            print(f"     and {r.example_other}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
