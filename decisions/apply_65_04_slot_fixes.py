"""The two slot decisions from instruction 65.04, ready to apply or to refuse.

Both were found by the biological eligibility read on 2026-11-13 and left UNAPPLIED, because each
changes the feature values every downstream model sees for that slot. The evidence and the measured
cost are in ``results/dataset_selection_2026_11_08/README_biological_eligibility.md``; this file is
only the hands.

Run it to apply both::

    python decisions/apply_65_04_slot_fixes.py            # show what would change
    python decisions/apply_65_04_slot_fixes.py --apply    # change parascape/data/slots.json
    python decisions/apply_65_04_slot_fixes.py --apply --only A

It is idempotent: applying twice changes nothing the second time, and it refuses rather than guesses
if a slot no longer looks the way it did when the evidence was gathered.

----------------------------------------------------------------------------------------------
**A. `Tg_transcription_bradyzoite_tissue_cyst` drops the pattern `rna108740_Day`.**

That pattern selects 14 samples of GSE108740 which the submitter declares `developmental stage:
Enteropithelial stages`, host `Felis catus` -- the stages in the gut of the cat, not bradyzoites in
a tissue cyst in a mouse. The slot averages 12 columns and **eight are those cat samples**.

Cost of NOT fixing it, over 7,865 genes: the slot disagrees with its own clean version at Spearman
0.79, with 3,105 genes (39 per cent) moving more than 10 percentile points and 973 (12 per cent)
moving more than 25. No rationale for the inclusion is recorded anywhere, and `expr_cyst` -- the
curated summary of the same deposit -- is documented as "the mean TISSUE-CYST FPKM", so the
hand-built column already excludes what the pattern pulls in.

**If you would rather keep the cat stages as evidence**, they deserve a slot of their own: no slot
currently asks about the enteroepithelial stages at all. Refusing this fix and adding that slot is
a coherent choice; leaving them inside a bradyzoite slot is not.

----------------------------------------------------------------------------------------------
**B. `Tg_protein_abundance_tachyzoite` narrows `proteome_` to the untreated columns.**

Its facets say `tachyzoite` and `unperturbed`, and `proteome_` selects 15 columns: 3 untreated, 9 of
AP2XII-1-depleted parasites that PRIDE titles "In vitro production of cat-restricted Toxoplasma
pre-sexual stages", and 3 log2 FOLD CHANGES averaged in with log2 abundances -- which the registry
note for that source already calls distinct quantification types. Each pattern is collapsed by row
mean before the policy runs, so all 15 become one number.

**The measured cost is small and that is part of the decision**: Spearman 0.987 against the
untreated columns alone, 126 genes (4 per cent) moving more than 10 percentile points, because the
columns are normalized and centred near zero. A real defect in what the slot claims to be, a minor
one in what it returns. Fix B is tidiness; fix A is substance.

Note that B does not address the other half: 102 genes take this slot from `protein_ibaq_log2`, the
Pru immunoprecipitation column, whose own note says "Enrichment, NOT a deep proteome". Removing that
too is a third decision and is deliberately not bundled here.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLOTS = os.path.join(ROOT, "parascape", "data", "slots.json")


def targets():
    """The repository copy ONLY, and the reason that is not an oversight.

    `slots.CATALOG` is `data_dir()/slots.json` -- the DOWNLOADED data release, not this repository.
    Editing only the repository file changes what ships at the next release and changes nothing that
    runs now, including the tests, which read the downloaded copy. A first version of this script
    edited only the repo file, re-ran 141 tests, watched them pass, and had proved nothing: the
    slots under test were the unmodified ones.

    The obvious correction -- write both -- is WRONG, and the suite says so. Hand-editing the
    downloaded copy fails `test_packaged_layout_matches_the_shared_matrix_and_gene_order`, because
    the data release pins the checksum of what it packaged. The release is meant to be rebuilt, not
    patched in place.

    So this writes the repository source, and the change reaches a running ParaScape when the data
    release is rebuilt and reinstalled. Both facts were established by trying it.
    """
    return [SLOTS]


DROP_PATTERN = "rna108740_Day"
WIDE = "proteome_"
NARROW = "proteome_log2_normalized_and_imputated_abundances_UT_"

# Fix A orphans eight columns, and instruction 41 forbids a column claimed by nobody as firmly as
# one claimed twice. They are not waste: they are the only transcriptome of the cat enteroepithelial
# stages in the map, and no slot asks about that stage at all. So fix A gives them a slot rather
# than dropping them on the floor. The facets follow the bradyzoite slot they are leaving, with the
# stage changed -- `in a host` and `gut` because these were taken from the intestine of a cat.
ENTERO_SLOT = {
    "organism": "Tg",
    "name": "transcription \u00b7 enteroepithelial stages",
    "axis": "transcription",
    "context": "enteroepithelial, cat intestine",
    "unit": "gene",
    "patterns": ["rna108740_Day"],
    "policy": "average",
    "evidence_path": ["molecular measurements", "RNA", "transcript abundance",
                      "enteroepithelial"],
    "biology_path": ["gene expression", "RNA abundance", "enteroepithelial", "unperturbed"],
    "context_path": ["in a host", "cat", "gut", "enteroepithelial", "transcription"],
    "target_family": "transcription_enteroepithelial",
    "target_columns": [],
    "role": "feature",
}

# Fix B orphans twelve columns for the same reason, so it carries the same obligation. They go to
# TWO slots rather than one, because putting an abundance and a log2 fold change in one slot is the
# defect fix B exists to remove: each pattern is collapsed by row mean, and a mean over two
# quantification types is not a quantity.
DEPLETION_SLOTS = [
    {
        "organism": "Tg",
        "name": "protein abundance \u00b7 under AP2XII-1 depletion",
        "axis": "protein abundance",
        "context": "AP2XII-1/AP2XI-2 depleted, 24-48 h, pre-sexual reprogramming",
        "unit": "gene",
        "patterns": ["proteome_log2_normalized_and_imputated_abundances_T_"],
        "policy": "separate",
        "evidence_path": ["molecular measurements", "protein", "protein abundance",
                          "merozoite"],
        "biology_path": ["gene expression", "protein abundance", "merozoite",
                         "genetic background"],
        "context_path": ["in vitro or assay-defined", "no host", "tissue-unspecified",
                         "merozoite", "protein abundance"],
        "target_family": "protein_abundance_under_ap2xii1_depletion",
        "target_columns": [],
        "role": "feature",
    },
    {
        "organism": "Tg",
        "name": "protein abundance \u00b7 AP2XII-1 depletion, fold change",
        "axis": "protein abundance",
        "context": "untreated against AP2XII-1 depleted, log2 fold change",
        "unit": "gene",
        "patterns": ["proteome_lfc_"],
        "policy": "separate",
        "evidence_path": ["molecular measurements", "protein", "protein abundance",
                          "merozoite"],
        "biology_path": ["gene expression", "protein abundance", "merozoite",
                         "genetic background"],
        "context_path": ["in vitro or assay-defined", "no host", "tissue-unspecified",
                         "merozoite", "protein abundance"],
        "target_family": "protein_abundance_ap2xii1_depletion_fold_change",
        "target_columns": [],
        "role": "feature",
    },
]


def load():
    with open(SLOTS) as fh:
        return json.load(fh)


def rows(doc):
    return doc if isinstance(doc, list) else doc["slots"]


def plan(doc, only: str = ""):
    """What each fix would do, without doing it."""
    out = []
    for r in rows(doc):
        if not isinstance(r, dict):
            continue
        pat = list(r.get("patterns") or ())
        if only in ("", "A") and r.get("context") == "bradyzoite" and DROP_PATTERN in pat:
            out.append(("A", r, [p for p in pat if p != DROP_PATTERN]))
        if only in ("", "B") and r.get("axis") == "protein abundance" and WIDE in pat:
            out.append(("B", r, [NARROW if p == WIDE else p for p in pat]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--only", default="", choices=["", "A", "B"])
    args = ap.parse_args()

    doc = load()
    todo = plan(doc, args.only)
    if not todo:
        print("nothing to do: neither slot still carries the pattern this would change.")
        print("either the fixes are applied already, or the slots have moved on and the")
        print("evidence in README_biological_eligibility.md should be re-read before acting.")
        return 0

    for label, slot, after in todo:
        print(f"[{label}] {slot['organism']}_{slot['name']}")
        print(f"      before: {slot['patterns']}")
        print(f"      after : {after}")
    if not args.apply:
        print("\nnothing written. Re-run with --apply to change parascape/data/slots.json.")
        return 0

    for _label, slot, after in todo:
        slot["patterns"] = after
    if any(label == "A" for label, _s, _a in todo):
        existing = {r.get("target_family") for r in rows(doc) if isinstance(r, dict)}
        if ENTERO_SLOT["target_family"] not in existing:
            rows(doc).append(dict(ENTERO_SLOT))
            print(f"[A] added slot: Tg_{ENTERO_SLOT['name']} "
                  f"<- {ENTERO_SLOT['patterns']}, so no column is left claimed by nobody")
    if any(label == "B" for label, _s, _a in todo):
        existing = {r.get("target_family") for r in rows(doc) if isinstance(r, dict)}
        for extra in DEPLETION_SLOTS:
            if extra["target_family"] not in existing:
                rows(doc).append(dict(extra))
                print(f"[B] added slot: Tg_{extra['name']} <- {extra['patterns']}")
    written = targets()
    for path in written:
        with open(path, "w") as fh:
            json.dump(doc, fh, indent=1)
            fh.write("\n")
    print(f"\napplied {len(todo)} change(s) to {len(written)} file(s):")
    for path in written:
        print(f"   {path}")
    print()
    print("This changes the REPOSITORY source only. It reaches a running ParaScape when the data")
    print("release is rebuilt and reinstalled -- the downloaded copy must not be patched in place.")
    if any(label == "A" for label, _s, _a in todo):
        print()
        print("Fix A moved rna108740_Day into a slot of its own rather than dropping it, because")
        print("leaving those eight columns claimed by nobody fails the column partition -- found")
        print("by applying it, not by reading. The cat enteroepithelial stages now have the only")
        print("slot in the map that asks about them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
