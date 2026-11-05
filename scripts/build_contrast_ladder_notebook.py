"""Build the contrast ladder notebook: how each study is read against its own baselines.

Instruction 79. The keyword rule it replaces read `WT` in `T. gondii VEG WT, Switch Conditions` and marked
the alkaline-stressed arm of GSE218921 as a baseline, which deleted the best alkaline contrast in the
collection. This notebook shows the parse, the pairing, the refusals and the recovered contrast, executed
against the live sample registry so the numbers in the prose are the numbers in the cells.

    python scripts/build_contrast_ladder_notebook.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())

from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_05"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"contrast_ladder_{DATE}")
NB = os.path.join(REPO, "notebooks", f"contrast_ladder_{DATE}.ipynb")
SAMPLES = os.path.join(REPO, "results", "rna_sample_quality_2026_11_04", "sample_quality.tsv")


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    nb = ExecutedNotebook("The contrast ladder: every comparison a study supports, "
                          "each against its own baseline")

    nb.md(
        "## What went wrong first",
        "Control pairing used to ask one question of each sample group: *does this label look like a "
        "control?* So it read `WT` in `T. gondii VEG WT, Switch Conditions` and marked the "
        "**alkaline-stressed** arm of GSE218921 as a baseline.",
        "`WT` is a level on the **strain** axis. The baseline for a stress question lives on the "
        "**treatment** axis. Having claimed the treated group as the control, nothing was left to pair it "
        "with, so the contrast was never formed and no column for it entered the information map -- absent "
        "rather than wrong.",
        "The cost was the cleanest bradyzoite-switch signal in the collection: BFD1 +1.36, BAG1 +4.75, "
        "ENO1 +4.91, SAG1 -1.34. The study fell back to internal normalisation against a centre built from "
        "all six groups mixed, switch samples included, so part of the stress response sat in the baseline "
        "and cancelled -- which is the mechanism behind cross-study agreement on alkaline stress of 0.236.",
        "The fix is not a longer keyword list. It is to parse the design and let the baselines nest.",
    )

    nb.code(
        "import importlib.util, csv, collections",
        f"spec = importlib.util.spec_from_file_location('cl', {os.path.join(REPO, 'scripts', 'build_contrast_ladder.py')!r})",
        "CL = importlib.util.module_from_spec(spec); spec.loader.exec_module(CL)",
        f"rows = list(csv.DictReader(open({SAMPLES!r}), delimiter='\\t'))",
        "len(rows), len({r['study'] for r in rows})",
    )

    nb.md(
        "## Step 1: the design is already written down",
        "GEO submitters record the structured part as `axis: level` pairs separated by `;`. The split is "
        "never on `,` -- a treatment level legitimately contains commas, as `high pH, low C02, low serum "
        "conditions` does.",
    )
    nb.code(
        "groups = collections.defaultdict(set)",
        "for r in rows: groups[r['study']].add(r['condition_group'])",
        "axes = collections.Counter()",
        "structured = 0",
        "for gs in groups.values():",
        "    for g in gs:",
        "        f, _ = CL.parse_factors(g)",
        "        structured += bool(f)",
        "        axes.update(f.keys())",
        "print('groups with a parseable design: %d of %d' % (structured, sum(len(v) for v in groups.values())))",
        "axes.most_common()",
    )

    nb.md(
        "## Step 2: the axis a depositor left in free text",
        "GSE218921 is the case that matters. All six tachyzoite groups carry an **identical** structured "
        "head -- `strain: VEG; treatment: ...; tissue: ...`. The genotype distinguishing WT from ROCY1 KO "
        "and BFD1 KO exists only in the sample name. Without recovering it the three strains collapse into "
        "one group and no genotype axis is seen at all.",
    )
    nb.code(
        "for g in sorted(groups['GSE218921']):",
        "    f, tail = CL.parse_factors(g)",
        "    if 'genotype' not in f:",
        "        f['genotype'] = CL.genotype_from_tail(tail) or 'unstated'",
        "    print('%-12s %-34s %s' % (f['genotype'], f['treatment'][:34], f['tissue'][:40]))",
    )

    nb.md(
        "## Step 3: pair on exactly one axis, and orient by the reference level",
        "Two groups are comparable when they differ on **one** axis. The group sitting at that axis "
        "reference level -- untreated, wild type, earliest time -- is the baseline, and every row records "
        "the word that decided it.",
        "Strain, tissue and stage never get an automatic reference. Choosing a strain is choosing an "
        "organism, not choosing a control, and silently making ME49 the baseline is how a strain "
        "difference gets reported as a treatment effect.",
    )
    nb.code(
        "contrasts, unresolved, refused = CL.build(rows, 'condition_group', 'study')",
        "print('contrasts: %d across %d studies' % (len(contrasts), len({c['study'] for c in contrasts})))",
        "collections.Counter(c['axis'] for c in contrasts).most_common()",
    )

    nb.md(
        "## The ladder for GSE218921",
        "`depth` counts how many **other** axes sit off their own reference. Depth 0 is the clean effect on "
        "an otherwise-baseline background; depth 1 is the same effect inside a mutant or a stressed "
        "background. Nothing is averaged across depths -- the treatment effect inside a knockout is a "
        "different measurement from the effect in wild type.",
        "The first treatment row is the contrast the keyword rule deleted. The depth-1 treatment rows are "
        "the nesting the design calls for: each knockout measured against its own untreated cells. The "
        "depth-1 genotype rows are the interaction terms, which the construction gives for free.",
    )
    nb.code(
        "for c in sorted((c for c in contrasts if c['study'] == 'GSE218921'),",
        "                key=lambda c: (c['axis'], c['depth'])):",
        "    print('d%s %-9s %-26s -> %-22s [%s]' % (c['depth'], c['axis'], c['baseline_level'][:26],",
        "                                           c['treated_level'][:22], c['background'][:52]))",
    )

    nb.md(
        "## Refusing rather than guessing",
        "A pair is formed only when the text no parsed factor explains is the same on both sides. If it "
        "differs, the study varies on an axis nobody parsed, and pairing across it compares two things at "
        "once. GSE132238 declares genotype and timepoint while `stress`/`unstress` hides in the sample "
        "name; without this rule `WT_stress_24h` is compared to `KO_unstress_24h` and a stress response is "
        "reported as a genotype effect.",
        "Leftover words are attributed to an axis by asking the study rather than a word list: a word "
        "belongs to axis A when its presence is constant inside every level of A. That recognises `Switch "
        "Conditions` as a restatement of `treatment: high pH` while leaving `stress` unexplained and "
        "therefore disqualifying.",
    )
    nb.code(
        "print('pairs refused: %d (real axis %d, plate or batch label %d)' % (",
        "    len(refused),",
        "    sum(1 for r in refused if r['looks_like'] == 'a real condition axis'),",
        "    sum(1 for r in refused if r['looks_like'] == 'plate or batch identifier')))",
        "collections.Counter(r['unparsed_difference'] for r in refused",
        "                    if r['looks_like'] == 'a real condition axis').most_common(8)",
    )

    nb.md(
        "## Two measurements must never print as one row",
        "Where an unparsed axis is the only thing separating two contrasts, it is named in the background "
        "(`unparsed=stress`) rather than dropped. That keeps the rows distinct and tells the reviewer "
        "exactly what needs naming. The count below must stay zero: identical rows get averaged together "
        "downstream.",
    )
    nb.code(
        "sig = collections.Counter((c['study'], c['axis'], c['baseline_level'], c['treated_level'],",
        "                           c['background']) for c in contrasts)",
        "sum(v for v in sig.values() if v > 1)",
    )

    nb.md(
        "## What is left for a person",
        "Two kinds: axes where no level matches the reference vocabulary (`etoh | shield1` -- Shield-1 "
        "stabilisation against an ethanol vehicle), and axes that exist only in free text (`bystander | "
        "infected`, `gra18 | pru`).",
        "This is the review item 79 described, an order of magnitude smaller: naming the reference level of "
        "these axes rather than naming the control for 450 runs.",
    )
    nb.code(
        "for u in unresolved[:12]:",
        "    print('%-12s %-22s %s' % (u['study'], u['axis'][:22], u['levels'][:54]))",
        "len(unresolved)",
    )

    nb.md(
        "## Gate",
        "Half met. The gate in instruction 79 is that cross-study agreement on alkaline stress rises from "
        "0.236 **and** the GSE218921 contrast exists as a column. The contrast is built and correctly "
        "oriented, which is the second half. The first half needs the expression columns rebuilt from this "
        "ladder, which is the next step and is not in this notebook.",
    )

    os.makedirs(os.path.dirname(NB), exist_ok=True)
    nb.write(NB)
    print("wrote", NB)


if __name__ == "__main__":
    main()
