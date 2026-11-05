# 79 · Within-study control review: unlock the relative expression columns
**Status: open (requested 2026-11-05).** The highest-value curation left on the RNA side, identified while
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
