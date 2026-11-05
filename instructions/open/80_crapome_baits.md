# 80 · Identify the pull-down baits so the background becomes usable
**Status: open (requested 2026-11-05).** The blocker found by the CRAPome first pass.

## Why
The background frequency column is **not usable as a filter**: the known complexes sit inside the frequent-flyer
range (MyoA glideosome 0.38-0.88, AMA1-RON 0.76-0.85 against actin 0.97, HSP70/90 0.76-0.94), separation AUC is
only 0.743, and at a 50 % cutoff **12 of 13 complex members would be deleted**. The cause is biological: 34
experiments over **11 distinct baits**, nearly all from whole parasite lysate, so the abundant invasion
machinery is in every list. The human CRAPome had over 400 experiments. **Identifying more baits is worth more
than finishing the remaining 140 projects.**

## The plan
1. **Resolve the 10 ambiguous baits** whose symbol appears in no ToxoDB product name (SPARKEL, Cdt1, MyoL, VMP1,
   CoqFMO, PCNA, TgSIR2A and the rest) from the depositing paper rather than by guessing.
2. **Classify the 8 `other` projects** whose titles name no method; several look like genuine pull-downs
   (PXD011692 guanylate-cyclase platform, PXD044588 apical annuli, PXD063615 anti-TgPDIA3).
3. **Handle the non-protein baits** (succinate dehydrogenase, mETC supercomplex, mitochondrial ribosome, conoid,
   glycan beads) as complex-level baits rather than discarding them.
4. **Weight by independence, not count**: SPARK/SPARKEL appear in five deposits and CDPK1 in three, so 34
   experiments are fewer than 34 independent observations and the report must say so.
5. **Keep the specificity score as the usable product** in the meantime: the tagged protein is in the top 5 % of
   its own hit list in 8 of 11 affinity purifications (median rank 14 of ~580 against 28 by abundance alone).

## Gate
The frequency filter becomes usable only when a cutoff exists that removes half the frequent flyers while
keeping 85 % of the known complexes. Until then the README says so and the score is used for ranking only.
