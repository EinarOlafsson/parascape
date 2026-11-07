# 75 · Condition fingerprints: annotate every dataset with metadata nobody published
**Status: v1 2026-11-06 -- 2 of 8 axes shipped (cell cycle, bradyzoite), each validated on held-out studies; the other 6 open.** Requested 2026-11-05. Moved out of instruction 71 product 5 so it has its own item.

## The request (Einar)
Build a fingerprint from a dataset that samples a condition axis in distinct samples (cell-cycle stages, for
example), then score every other dataset against it. "That way we not only compare, we annotate new metadata not
known before to the datasets." On both the RNA and the protein level.

## The plan
1. **Output proportions, not labels.** A bulk sample is a mixture: an asynchronous culture holds every
   cell-cycle stage, an alkaline-stressed culture is partly converted. "68 % tachyzoite, 32 % bradyzoite" is both
   more honest and more useful, and estimated composition becomes a covariate that explains between-study
   differences currently attributed to biology.
2. **Score by within-sample ranks** (ssGSEA-style), never absolute values, or a fingerprint learned in one study
   scores that study's batch. Rank scoring is also what makes one fingerprint usable across platforms and
   between RNA and protein.
3. **Build one axis at a time**, in this order, by how much data backs it: cell cycle (28 columns, a
   synchronised 12-hour series with async and blocked controls, ~5,377 periodic genes); bradyzoite
   differentiation (33 columns over several studies plus 18 knockdown columns, so it can be validated *across*
   studies); translational state (46 Ribo-seq and efficiency columns, orthogonal to transcript level); oocyst,
   merozoite and enteroepithelial (16); extracellular versus intracellular (12); culture conditions -- pH,
   serum, lipid, glucose, glutamine, host density (42); strain type (GT1/ME49/VEG, already quantified per run);
   protein level (12-column proteome time course, 8-column iTRAQ oocyst, plus the Sage tables as they finish).
4. **Validate on held-out studies that state their own condition.** Build on one study, score another that
   reports its stages, measure recovery per axis. Expect RNA -> protein transfer to be weaker than RNA -> RNA and
   report the two separately rather than assuming.
5. **An annotation is inferred, never observed.** It enters the map as `filled:<fingerprint>@<version>` or stays
   exploratory, and its source study is declared in `depends_on`. Otherwise the failure is obvious once stated:
   annotate a dataset with an inferred stage, then use that stage as evidence for a trait, and the model learns
   from its own guess.

## Gate
Each axis ships only when a held-out study that reports its own condition is recovered, with the number stated.
An axis that cannot be validated that way is published as exploratory or not at all.

## Two directions, one word (2026-11-05)
This item fingerprints a **condition** and annotates **samples**. Instruction **82** fingerprints a **gene set** from a literature meta-analysis and annotates **genes**. Same word, opposite axis of the matrix, different validation: a held-out study that states its condition here, held-out members of the set there. Keep them separate in the registry so a scorecard cannot be read across them.

## v1, 2026-11-06: cell cycle and bradyzoite conversion

`scripts/build_condition_fingerprints_notebook.py` -> `notebooks/condition_fingerprints_2026_11_06.ipynb`,
`results/condition_fingerprints_2026_11_06/`. 668 bulk runs of the 735-run snapshot scored by within-sample
percentile ranks (37 runs have no values in the snapshot and are dropped: GSE243206, GSE243207, GSE245775,
GSE98227).

| axis | built on | tested on (held out) | result |
|---|---|---|---|
| cell cycle | Xue 2020 single-cell phases (873 genes) | synchronised series `cellcycle19092`, hours 1-12 | **all 11 hourly steps advance the phase** (shuffled hours 50 %, p 0.001); 1.44 cycles in 12 h |
| bradyzoite (alkaline response) | GSE223869 wild type, alkaline vs standard | 7 studies that state the condition | pooled AUC **0.934**, 6 of 7 studies above 0.7 |
| bradyzoite (conversion) | the same, and without BFD2 | the same | pooled AUC **0.922**, 6 of 7 above 0.7 |

**Converter removed under stress** (held out): BFD1 knocked out (GSE132238) or depleted (GSE223819) drops both
fingerprints in every run (AUC 1.0), from about 1.0 to 0.51-0.65. The conversion-specific version falls only
slightly further than the general one, so most of the alkaline response already depends on BFD1; it does not
reach zero, so part of what both read is stress that needs no converter.

**New metadata, as asked:** in GSE223819 two tagged lines have a genotype field and a sample title naming
different genes (208020 against 311100 = BFD2). The line titled 311100 converts poorly even on vehicle (0.29);
the fingerprint sides with the titles. Recorded as evidence in `gse223819_label_evidence.tsv`, not applied.

Outputs are proportions (`cellcycle_<phase>` mixture and angle; `bradyzoite_*_fraction` between the build
study own tachyzoite and converted means), each with `filled:<fingerprint>@v1` and `depends_on`.

**Open, with what each needs to be validated:** translational state (needs a second Ribo-seq study with
usable counts -- GSE245775 has none in this snapshot); oocyst / merozoite / enteroepithelial; extracellular
against intracellular; culture conditions; strain type (already known per run, so it is a positive control
rather than new metadata); protein level (RNA -> protein transfer, reported separately).

**Searched for the next axes, 2026-11-06** (condition text of the 717 verified runs): no run states
extracellular, intracellular, egress, oocyst or feline enteroepithelial conditions. "Merozoites" appears in
GSE224934 and GSE228580, but both are AP2XII-1 knockdowns from one system, and the label is applied to the
0 h wild-type control as well, so it names the experiment, not the stage -- not an independent stated
condition to validate against. Sporozoites appear in one study (GSE218921). These axes wait for data that
states them, or for a curated stage label per run.
