# 80 · Identify the pull-down baits so the background becomes usable
**Status: steps 1-5 done 2026-11-06. 2026-11-17: waiting for more projects will NOT open the gate, so the gate as written cannot be reached that way.** Finished projects went from about 73 to 152 and the separation AUC went DOWN, 0.710 to 0.689 to 0.684. More whole-lysate projects do not help. What is missing is baits from OTHER COMPARTMENTS, so no number of queued projects is the answer and the 51 still in the queue are unlikely to change it. The specificity score remains the usable product: the bait is in its own top 5 per cent in 27 of 34 affinity purifications, median rank 11 of a median 825 proteins against 24 by abundance alone.** Baits resolved from the depositing
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

## A compartment-matched background would NOT open the gate either (tested 2026-11-18)

The obvious next idea, and the one this project reaches for everywhere else: the baits are
compartment-biased -- of the 25 pulldowns whose bait has an assigned compartment, 12 are nuclear
and 8 apical secretory -- so compare a protein against pulldowns of baits from a DIFFERENT
compartment rather than against all of them.

It does not work, and the reason is that the premise is wrong. If complex members looked frequent
because the baits over-represent their compartment, the 13 known complex members sitting in
nucleus or apical secretory would be more frequent than the ones elsewhere. They are not: median
frequency 0.676 inside against 0.704 outside, p = 0.53. And the 82 expected frequent flyers are
almost all somewhere else entirely (81 of 82), at a median frequency of 0.817 -- they are abundant
housekeeping proteins such as GAPDH and HSP70/90, which a whole-lysate pulldown recovers whatever
the bait is and wherever it sits.

So frequency cannot be rescued by matching on compartment, because complex members and true
background are BOTH frequent for the same reason: whole lysate. That leaves the specificity score
-- abundance relative to the background rather than presence in it -- as the usable product, which
is what step 5 already concluded. Recorded here so the matched-background route is not built.

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

**Bait-resolution pass, same day.** The 23 unresolved or `other` projects were read from their PRIDE records and
papers (`results/pulldown_background_2026_11_08/bait_resolution_candidates.tsv`, each row with a quote and its
source). 18 new rows were added to `parascape/data/pulldown_baits_curated.tsv`:
- Confident baits: RON4, GRA57, CDPK1, TgHCF101, PDIA3 (TurboID), AP2X-1/MORC/HDAC3, BCC4/APP1.
- Not confident:
  - Nd10 and MAF1a/b: symbol only, no gene id.
  - Human GBP1: a host bait.
- Nine projects reclassified as not pull-downs: PTM enrichment, whole proteome, crosslinking.

A gene id written in the PRIDE record now counts as confident; the change is recorded in the builder.

Results:
- **36 experiments with a confident bait, 30 distinct baits (was 25).**
- Specificity: the bait is in its own top 5 % in 25 of 32 affinity purifications, median rank 10 of about 808.
- The frequency filter is unchanged (AUC 0.689; 12 of 13 complex members lost at 50 %): more baits do not change
  which proteins appear in every experiment. Gate unmet.

Still open: PXD011692 (paper not open access). Next step: the queue (106 projects).

## 2026-11-17: rerun on 152 finished projects (was 97); gate verdict: NOT open

`results/pulldown_background_2026_11_17/` (builder rerun unchanged apart from the date; 149 classified). Counted from
`done.json` receipts on the NAS: 152 of 203 selected projects finished (was 97 on 8 Nov, 73 at the first gate
check), 51 still in the queue, 12 unsupported, 2 HDX-MS excluded.
- Pull-down experiments in the background: **71** (affinity 58 + proximity 13; was 50). **38** have a confident
  bait, **32 distinct baits** (was 36 and 30). 22 projects are still class `other`.
- Frequency filter: separation AUC **0.684** (was 0.689, p 0.017). At a 50 % cutoff 12 of 13 complex members are
  still lost. The best cutoff that keeps at least 85 % of complex members (0.85) removes only 27 % of the frequent
  flyers (22 of 82); the cutoff that removes 62 % of flyers (0.80) loses 31 % of complex members (4 of 13). The
  gate (half the flyers removed, at most 15 % of complex members lost) is **not met, and the filter is not usable**.
- Specificity score still holds: bait in its own top 5 % in 27 of 34 affinity purifications, median rank 11 of a
  median 825 proteins (24 by abundance alone).
- Why more projects do not help: 55 more finished projects moved the AUC from 0.710 to 0.684, i.e. not at all.
  The gate needs baits from other compartments, not more whole-lysate projects; the remaining 51 queued projects
  are unlikely to change that. The recorded conclusion stands: use the specificity score, not the frequency filter.
