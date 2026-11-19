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

## Where this should land

The defensible statement is not a weak ipTM comparison. It is that **MYR1 lies in the fraction
AF3 cannot see** -- no ordered residues, below the foldability cut, peptide tiling already
failed. That is a result about the limits of the method, and a better sentence than the panel
I set out to draw. Awaiting Einar.

The two completed folds are kept: `TGME49_254470__Q99584` (ipTM 0.19) and
`TGME49_254470__Q06210`. The other 22 inputs remain in
`/media/carruthers/mnt1/af3/af_input/host_cofolds_2026_11_19/` and run again whenever wanted.
