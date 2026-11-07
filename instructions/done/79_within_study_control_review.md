# 79 · Within-study control review: unlock the relative expression columns
**Status: done 2026-11-06.** Gate met: wild-type alkaline stress agrees across studies at 0.439 (was 0.236), and the GSE218921 contrast is a column. Requested 2026-11-05. The highest-value curation left on the RNA side, identified while
building the expression columns.

## Why
`control_link` in the sample registry carries **no pairing at all**: every one of the 450 kept runs reads "No
verified pairing; choose the matching control within this study", and the snapshot states
`control_pairing: unverified`. Only **17 of 36 studies** yield a usable control, so 19 fall back to weaker
internal normalisation. The cost is concrete: in GSE218921 a keyword rule marked the **wild-type alkaline-switch
arm as a baseline** because its label contains "WT", so the best alkaline contrast in the snapshot never became
a column. Computed directly it is textbook -- BFD1 +1.36, BAG1 +4.75, ENO1 +4.91, SAG1 -1.34. It also explains
the one weak cross-study result: alkaline stress agrees at only 0.24 between studies, both columns being
internally normalised against centres made of different condition mixtures.

## The plan
A review of the same shape as the host-facing claim review that worked: 36 studies, one row per contrast, naming
the control group. Pre-filled from the registry where the marking is unambiguous, flagged where the keyword rule
is the only evidence, with the study title, the sample titles and the group means shown so the answer is
visible. Reviewed decisions rebuild the expression columns, converting 19 studies from internal reference to
relative-to-control.

## Gate
Cross-study agreement on alkaline stress rises from 0.236, and the GSE218921 contrast exists as a column.
Measured, not assumed.


## Automated 2026-11-05: the contrast ladder

Einar's correction, which reframed the item: *"the WT and treatment failure just means that if a study has
several levels of baselines this needs to be taken into account. In this case the untreated knockouts should
be related to WT then the treatments should be related to untreated. Can you do that? And do it
automatically?"*

Yes, and it is built: `scripts/build_contrast_ladder.py`, with
`results/contrast_ladder_2026_11_05/` and 26 tests in `tests/test_contrast_ladder.py`. The rule is no longer
"find the control" but "parse the design, then pair groups that differ on exactly one axis" -- which is what
makes the nesting fall out rather than having to be specified per study. Depth records how many other axes
sit off their reference, so the treatment effect in wild type and the same effect inside a knockout stay
separate measurements.

**It refuses rather than guesses.** A pair is formed only when the text no parsed factor explains is the same
on both sides; otherwise the study varies on an axis nobody parsed. 168 pairs refused on that basis, 36 of
them on a real condition axis. Leftover words are attributed to axes by testing the study -- a word belongs
to an axis when its presence is constant within every level of it -- which separates a restatement of a
declared factor from a genuinely hidden one, and is run inside strata of tissue, strain and stage because
those separate different experiments sharing one accession.

**85 contrasts across 15 of 31 studies, 0 rows that print identically.** The GSE218921 alkaline contrast
exists, at depth 0, correctly oriented from `control media` -- the half of the gate that was in reach.

**What is left for Einar shrank by an order of magnitude**: 30 axes needing a named reference level, not 450
runs needing a named control. Listed in `reference_needs_review.tsv` with the study, the levels and why the
vocabulary could not decide.

**Still open**: rebuilding the expression columns from the ladder, which is what the rest of the gate
(alkaline agreement above 0.236) measures. The ladder is the input to that, not the result of it.

## Done 2026-11-06: reviewed references applied, columns rebuilt, gate measured

**Reviewed decisions.** `parascape/data/contrast_reference_decisions.tsv` holds the 30 axes that needed a named
reference, each decided from the GEO text and the paper, with the evidence quoted. `build_contrast_ladder.py
--decisions` reads it; batch tokens are stripped before groups merge (168 pairs had printed identically
because of them), and a study whose only factor is `condition` now gets an axis. Result
(`results/contrast_ladder_2026_11_06/`): **178 contrasts across 19 studies, 0 indistinguishable, 0 needing a
reference.** 31 tests in `tests/test_contrast_ladder.py`.

**Columns rebuilt** (`scripts/build_expression_traits_notebook.py`,
`notebooks/expression_traits_2026_11_06.ipynb`, `results/expression_traits_2026_11_06/`). Every ladder contrast
becomes an `expr_lfc_{study}_{slug}` column, normalized `relative_to_control`, with `ladder_axis`,
`ladder_depth` and `oriented_by` in the registry. Columns with no value for any gene are no longer published
(8 dropped: GSE245775 ribosome profiling, whose runs carry zero quality weight, and GSE98227); the drop is a
row in `exclusions.tsv`.

**The gate** (`alkaline_gate.tsv`, computed in the notebook). Wild-type alkaline against normal medium:

| pair | Spearman |
|---|---|
| GSE218921 against GSE223869, both relative to control | **0.439** |
| GSE223869 against GSE208647 (the baseline pair, now with a real control on one side) | 0.373 |
| GSE218921 against GSE208647 | 0.021 |

Median 0.373 against the old 0.236; two of three pairs clear it. The low pair involves GSE208647, which the
ladder still cannot pair (its pH arms also differ in the declared stage, "Tachyzoite" against
"Bradyzoite"), so it is compared through the difference of its two internally referenced arms. GSE200962 is
not an alkaline contrast (both arms alkaline, the contrast is auxin) and is left out.

**Not in this item:** the app information map still reads the 2026-11-05 expression table
(`information_categories.EXPRESSION_COLUMNS` is pinned to it). The 2026-11-06 columns enter at the next
information-map release.
