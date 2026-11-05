# 75 · Condition fingerprints: annotate every dataset with metadata nobody published
**Status: open (requested 2026-11-05).** Moved out of instruction 71 product 5 so it has its own item.

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
