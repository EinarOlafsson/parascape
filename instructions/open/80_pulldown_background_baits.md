# 80 · Identify the pull-down baits so the background becomes usable
**Status: steps 1-5 done 2026-11-06; the gate waits on more projects.** Baits resolved from the depositing
papers (Europe PMC full text) or ToxoDB gene names into `parascape/data/pulldown_baits_curated.tsv`, each row
with its source and the identifying sentence: SPARKEL TGME49_291010, Cdt1/TgiRD1 TGME49_247040, Cyc4
TGME49_249880, GSK3 TGME49_265330, CoqFMO TGME49_232540, MyoL TGME49_291020, PCNA1, PDIA3, PP2A-A2/C2,
AP2XI-2, TgH1-like; four 'other' projects reclassified as affinity purifications; SIR2A kept as a symbol match
(paper unreadable) and PXD011692 left unassigned rather than guessed. Rebuilt notebook
`notebooks/pulldown_background_2026_11_06.ipynb`: **22 distinct confident baits over 41 experiments (was
11 over 34)**; bait in its own top 5 % in 17 of 22 affinity purifications, median rank 6 of ~778 (28 by
abundance). The frequency filter is still not usable (AUC 0.710; 11 of 13 complex members lost at 50 %):
the gate needs baits from more compartments, which come with the PRIDE projects still being processed
(instruction 88). The specificity score stays the usable product.

Requested 2026-11-05. The blocker found by the first pass of the background table.

## Why
The background frequency column is **not usable as a filter**: the known complexes sit inside the frequent-flyer
range (MyoA glideosome 0.38-0.88, AMA1-RON 0.76-0.85 against actin 0.97, HSP70/90 0.76-0.94), separation AUC is
only 0.743, and at a 50 % cutoff **12 of 13 complex members would be deleted**. The cause is biological: 34
experiments over **11 distinct baits**, nearly all from whole parasite lysate, so the abundant invasion
machinery is in every list. The human CRAPome (Mellacheruvu et al. 2013), the published resource ours imitates, had over 400 experiments. **Identifying more baits is worth more
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

## 2026-11-07: rerun on 77 finished projects (was 73)

`results/pulldown_background_2026_11_07/`: 74 classified. Separation of frequent flyers from complex members
AUC 0.718 (was 0.710), p 0.006; at a 50 % cutoff 11 of 13 complex members would still be discarded. The gate is
still unmet; 126 of 203 projects have not finished the queue. Rerun the builder as the queue fills.

## 2026-11-08: rerun on 97 finished projects (was 77)

`results/pulldown_background_2026_11_08/`: 94 classified; 50 background experiments (affinity 40 + proximity 10),
29 with a confident bait, 25 distinct baits (was 22). Separation of frequent flyers from complex members AUC 0.689
(was 0.718), p 0.015; at a 50 % cutoff 12 of 13 complex members would still be discarded, and no cutoff meets the
gate (half the flyers removed, at most 15 % of complex members lost). More projects of the same kind (whole-parasite
lysate) do not move it; baits from other compartments are what is needed. The specificity score holds: bait in
its own top 5 % in 20 of 25 affinity purifications, median rank 8 of ~808 (31 by abundance). 18 projects are still
class `other` and 5 affinity purifications have no confident bait -- the next bait-resolution pass. 106 of 203
projects remain in the queue.
