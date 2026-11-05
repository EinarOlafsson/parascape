# Contrast ladder — every comparison a study supports, each against its own baseline

Built 2026-11-05 by `scripts/build_contrast_ladder.py` from
`results/rna_sample_quality_2026_11_04/sample_quality.tsv`. Replaces the keyword rule that assigned controls
by matching words in a sample label.

## What the keyword rule got wrong

It asked one question of each sample group — "does this label look like a control?" — and so it read `WT` in
`T. gondii VEG WT, Switch Conditions` and marked the **alkaline-stressed** arm of GSE218921 as a baseline.
`WT` is a level on the *strain* axis; the baseline for a stress question lives on the *treatment* axis. Having
claimed the treated group as the control, nothing was left to pair it with, so the contrast was never formed
and **no column for it entered the information map** — absent, not wrong.

The cost was the cleanest bradyzoite-switch signal in the collection: BFD1 +1.36, BAG1 +4.75, ENO1 +4.91,
SAG1 −1.34. The study instead fell back to internal normalisation against a centre built from all six groups
mixed, switch samples included, so part of the stress response sat in the baseline and cancelled. That is the
mechanism behind cross-study agreement on alkaline stress of only 0.236.

## The rule now

A study with several levels of baseline gets several levels of contrast. Untreated knockouts are compared to
untreated wild type; treated samples are compared to the untreated samples **of their own genotype**.

1. Parse each group label into factors — treatment, genotype, strain, timepoint, tissue. 178 of 209 groups
   already carry explicit `axis: level` structure, and GEO writes it as `;`-separated pairs.
2. Keep the axes that vary inside that study.
3. Pair groups differing on **exactly one** axis.
4. Orient the pair: the group at the axis reference level — untreated, wild type, earliest time — is the
   baseline. Every row records the word that decided it in `oriented_by`.
5. `depth` counts how many *other* axes sit off their own reference. Depth 0 is the clean effect on an
   otherwise-baseline background; depth 1 is the same effect inside a mutant or a stressed background.

Nothing is averaged across depths. The treatment effect inside a knockout is a different measurement from the
effect in wild type, and the ladder keeps them apart.

## GSE218921, the case that drove this

| depth | axis | baseline → treated | background |
|---|---|---|---|
| 0 | treatment | control media → high pH, low CO₂, low serum | genotype=wild type; tachyzoites |
| 0 | genotype | wild type → bfd1 ko | untreated; tachyzoites |
| 0 | genotype | wild type → rocy1 ko | untreated; tachyzoites |
| 1 | treatment | control media → high pH… | genotype=bfd1 ko |
| 1 | treatment | control media → high pH… | genotype=rocy1 ko |
| 1 | genotype | wild type → bfd1 ko | under stress |
| 1 | genotype | wild type → rocy1 ko | under stress |
| 1 | treatment | control media → high pH… | sporozoites, genotype unstated |

Row 1 is the contrast the keyword rule deleted. Rows 4–5 are the nesting Einar asked for. Rows 6–7 are the
interaction terms, which the construction gives for free.

## Refusing rather than guessing

A pair is formed only when the text left over after removing every parsed factor is the same on both sides.
If it differs, the study varies on an axis nobody parsed, and pairing across it compares two things at once.
GSE132238 declares genotype and timepoint while `stress`/`unstress` hides in the sample name; without this
rule `WT_stress_24h` is compared to `KO_unstress_24h` and a stress response is reported as a genotype effect.

Leftover words are attributed to an axis by asking the study, not a word list: a word belongs to axis A when
its presence is constant inside every level of A. That recognises `Switch Conditions` as a restatement of
`treatment: high pH` while leaving `stress` unexplained and therefore disqualifying. The test runs **inside
strata** of tissue, strain and stage, because those axes separate different experiments sharing one
accession — the GSE218921 sporozoite arm words the same stress as `bradyzoite induction`, and study-wide
that made `switch` look unexplained and refused the alkaline contrast a second time.

## Numbers

| | |
|---|---|
| studies read | 31 |
| studies yielding contrasts | 15 |
| contrasts | 85 (treatment 27, genotype 37, time 21) |
| depth 0 / depth 1 | 46 / 39 |
| pairs refused for an unparsed axis | 168 |
| …of those, on a real condition axis | 36 |
| …on what looks like a plate or batch label | 132 |
| axes needing a named reference | 30 |
| rows that print identically | **0** |

## Three quieter mistakes this also fixes

- **An axis some groups declare and others omit** is skipped silently. GSE136123 has two different untreated
  baselines, `Pru MORC KD_UT` and `Pru_US`; only the first parses a genotype, so both read as the same
  baseline and the knockdown was measured against the parental line as though it were the same strain. A
  missing axis now gets an explicit `unstated` level, which matches no reference vocabulary, so the resulting
  contrast goes to review instead of being oriented on a guess.
- **Replicates deposited as separate groups.** GSE132238 submits `..._array1` and `..._array2`; left apart
  they emit the same contrast twice, each with half the replicates it really has. Groups with identical
  factors and identical leftover text are merged.
- **Two measurements printing as one row.** `indistinguishable_rows` must stay 0. Where an unparsed axis is
  the only thing separating two contrasts, it is named in the background (`unparsed=stress`) rather than
  dropped — which both keeps the rows distinct and tells the reviewer exactly what needs naming.

## Files

| file | contents |
|---|---|
| `contrasts.tsv` | one row per contrast: axis, baseline and treated level, depth, background, `oriented_by`, replicate counts, both group labels |
| `refused_unparsed_axis.tsv` | every refused pair, the hidden words, and whether they look like a real axis or a plate label |
| `reference_needs_review.tsv` | axes whose reference level could not be named automatically — the human review list |
| `summary.json` | the table above |

## What still needs Einar

`reference_needs_review.tsv`, 30 rows. Two kinds: axes where no level matches the reference vocabulary
(`etoh | shield1` — Shield-1 stabilisation against an ethanol vehicle), and axes that exist only in free text
(`bystander | infected`, `gra18 | pru`). Strain, tissue and stage never get an automatic reference by design,
because choosing a strain is choosing an organism rather than a control, and silently making ME49 the
baseline is how a strain difference gets reported as a treatment effect.

This is far smaller than the review item 79 originally described: naming the reference level of ~30 axes
rather than naming the control for 450 runs.

## Gate (instruction 79)

Not yet met. The gate is that cross-study agreement on alkaline stress rises from 0.236 **and** the
GSE218921 contrast exists as a column. The second half is done — the contrast is built and oriented. The
first half requires rebuilding the expression columns from this ladder, which is the next step and is not in
this folder.
