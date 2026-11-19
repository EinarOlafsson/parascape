# Figure 8b: why the AF3 batch stopped at 2 of 24

Paused 2026-11-19 on Einar's instruction, after consulting the AF3 folding session
(`repo-1f`, `(work) AF3 large parasite proteins`), which I should have done BEFORE starting.

## The decisive number, verified here and not taken on trust

`toxoplasma_projects/af3_folding/reports/bait_folded_fraction.json`:

| | |
|---|---:|
| MYR1 (TGME49_254470) folded fraction | **0.0** |
| median across their 208 baits | 0.59 |
| baits at exactly 0.0 | 13 of 208 (6 %) |

Zero of 810 residues sit above AFDB pLDDT 70. Their screen uses a 30 % foldability cut, so
**MYR1 was excluded from it**, and the only piece ever co-folded is the fragment 638-810.

## Why that ends the panel as designed

Over ~400 PDB complexes with true DockQ they established that AF3 maps the folded,
domain-mediated fraction well (median DockQ 0.92) and is **at chance for disordered,
motif-mediated binding** (iLIS AUC 0.515). The first result here, ipTM 0.19, is what a fully
disordered bait looks like; 22 more would have looked the same, and the two arms would have
differed by whatever noise dominated rather than by binding.

The obvious rescue is already closed: `pep_TGME49_254470_*` and `decoypep_TGME49_254470_*` in
`af_output` are peptide tiling of MYR1 against human partners with shuffled decoys. The real
motif won 27 % of the time and the decoys made significantly MORE contacts.

## Two things I had wrong beyond the bait

* **ipTM is the wrong score to rank on.** Spearman against true DockQ is **-0.014** (p = 0.78);
  it saturates, with 353 of 408 pairs inside 0.88-0.96. Rank on interface pLDDT (0.502) or
  interface PAE (-0.585). Never on contact count, which runs **-0.519** against DockQ.
* **No MSA to reuse.** Their cache holds only the 638-810 fragment for MYR1, not full length.

## What was NOT wrong

Their hardest-won trap is length confounding: their positives averaged 800 aa against negatives
at 442, and length alone separated the arms at AUC 0.976. Checked here: host-length
**AUC 0.583, Mann-Whitney p = 0.51**. Not that confound -- but with 12 per arm that is
underpowered agreement, not the matched 0.500 their method requires.

## Both folds scored, and they are below the WRONG-POSE distribution

The AF3 session ran its 44-feature scorer over the two finished folds
(`af3_session_scored_2folds.csv` in the output directory):

| | Q99584 | Q06210 | high-accuracy poses | wrong poses |
|---|---:|---:|---:|---:|
| interface pLDDT | 53.8 | 57.2 | median 98.1 | median 69.7 |
| interface PAE (A) | 20.8 | 12.0 | median 1.0 | median 15.8 |
| iLIS / cLIS | 0.000 | 0.000 | | |
| pDockQ2 | 0.009 | 0.009 | | |

Both sit BELOW the median of wrong poses. iLIS and cLIS at exactly 0.000 mean essentially no
interface residue pairs fall under the 12 A PAE cutoff: not a weak interface, an absent one.

## The length confound is visible INSIDE my own two folds

This matters more than the scores. Q06210 looks better than Q99584 on ipTM (0.33 vs 0.19),
contacts (219 vs 60) and BSA (7718 vs 2599) -- and its partner is 699 aa against 98, seven times
longer. The length-INSENSITIVE measures barely move: interface pLDDT 53.8 to 57.2, iLIS 0.000 to
0.000, pDockQ2 identical. Every metric that appeared to separate the two pairs was tracking
partner length. That is the 0.583 host-length AUC showing its teeth at n = 2, and it confirms the
reading that 0.583 at n = 12 is underpowered agreement, not demonstrated balance -- detecting an
imbalance of the kind that bit them (0.976) would need roughly 25-30 per arm.

## Where this should land

The defensible statement is not a weak ipTM comparison. It is that **MYR1 lies in the fraction
AF3 cannot see** -- no ordered residues, below the foldability cut, peptide tiling already
failed. That is a result about the limits of the method, and a better sentence than the panel
I set out to draw.

**KEEP THE CLAIM ABOUT THE METHOD, NOT THE PROTEIN** (their caution, and they are right). These
data say AF3 cannot model MYR1 interfaces and that MYR1 has no confidently ordered domain in
AFDB. They do NOT say MYR1 lacks structure in vivo or lacks partners: MYR1 is a translocon
component, and its IP-MS partners may be real but mediated by regions or conditions AF3 cannot
represent, or indirect through the translocon. Nor does ipTM 0.19 make Q99584 a non-binder -- the
honest statement is that the prediction is uninformative, which is a different thing.

Five independent lines converge on the method limit: folded fraction 0.0 against a median of
0.59; below the 30 per cent cut, so excluded from a 40,675-pair screen by rule; peptide tiling
with shuffled decoys already failed (real motif won 27 per cent, decoys made MORE contacts);
phosphorylating empirical ToxoDB sites did not raise local pLDDT (+1.76, p = 0.42); and now two
co-folds below the wrong-pose distribution.

Still open, both suggested by the AF3 session: whether any of the 733 partners pairs with a MYR1
region that IS ordered (none is above the cut, so probably not, but the 638-810 fragment they
hold makes it a one-line check); and whether the IP-MS enrichment is simply the better panel,
since 733 partners with an imputation-aware class split stand on their own with no structure
prediction at all. MSA reuse is marginal here: 19 of my 651 partner accessions are in their
cache of 72, and their entries are FRAGMENTS keyed to exact residue ranges.

The two completed folds are kept: `TGME49_254470__Q99584` (ipTM 0.19) and
`TGME49_254470__Q06210`. The other 22 inputs remain in
`/media/carruthers/mnt1/af3/af_input/host_cofolds_2026_11_19/` and run again whenever wanted.
