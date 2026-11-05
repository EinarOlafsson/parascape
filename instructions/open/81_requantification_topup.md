# 81 · Re-quantify the earlier runs: trimming, equivalence classes and Gibbs samples
**Status: open (agreed 2026-11-03, scheduled after the new downloads).**

## Why three things accumulated
1. **Adapter trimming** arrived on 2026-11-04. Runs verified before it carry no `trimming` key in `done.json`,
   which is how an untrimmed result is recognised. Normal libraries lose almost nothing (0.0-0.1 % of reads), but
   40 runs of GSE243206, GSE243207 and GSE245775 mapped 0.00 % until trimming and are excluded from the
   expression columns until this runs.
2. **Equivalence classes and 30 Gibbs samples** (`--dumpEqWeights --numGibbsSamples 30`) arrived the same week,
   so runs can be re-quantified against an updated annotation and carry per-transcript uncertainty without their
   reads. Earlier runs lack both.
3. The reads are deleted after quantification, so this is a re-download.

## The plan
Re-download and re-quantify the runs verified before each change, appending the missing outputs **without
altering a verified result**: a new quantification is written beside the old one and the old numbers stay until
the new ones are checked against them. Einar's instruction stands: do this only after the new data is
downloaded. The 40 zero-parasite runs come first, because they currently contribute nothing at all.

## Gate
The 40 runs produce non-zero parasite mapping; every re-quantified run agrees with its previous values within
the expected tolerance, and any that does not is reported rather than overwritten.
