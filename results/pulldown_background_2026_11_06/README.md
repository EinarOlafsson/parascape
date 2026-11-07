# Toxoplasma background frequency, first pass (2026-11-06)

Instruction 71 product 4: how often each parasite protein is detected in an affinity-purification or proximity-labelling experiment regardless of what was pulled down, plus a specificity score for a new hit list.

Notebook: `notebooks/pulldown_background_2026_11_06.ipynb` (executed). Builder: `scripts/build_pulldown_background_notebook.py`. Re-running the builder rescans the uniform search output and picks up whatever has finished since.

## Coverage

| | n |
|---|---|
| projects selected for the uniform re-search | 203 |
| project directories on disk | 86 |
| finished (done.json receipt) | 73 |
| excluded: unsupported.json | 11 |
| excluded: NOT_COMPARABLE.json (HDX-MS) | 2 |
| finished and classified | 71 |
| still waiting in the queue | 130 |
| experiments in the background (affinity + proximity) | 41 |
| of those with a confident bait gene | 26 |
| distinct confident bait genes | 22 |
| runs behind those experiments | 274 |
| PSM rows read | 10757620 |
| PSMs at 1% FDR, target | 2682373 |
| parasite proteins in the background table | 4528 |
| class: affinity_purification | 32 |
| class: fractionation | 2 |
| class: other | 7 |
| class: proximity_labelling | 9 |
| class: ptm_enrichment | 10 |
| class: whole_proteome | 11 |

## Classification rules, verbatim

Applied in this order to the PRIDE title first and the deposited run file names second; first match wins. The matched text is kept in `project_classification.tsv` (`cue`, `decided_by`).

```
hdx                    \bhdx\b|hdxms|hydrogen[- ]deuterium
                         (HDX-MS: masses are deuterium-shifted, not abundance (instruction 69))
proximity_labelling    bioid|bira|turboid|turbo-id|\bapex2?\b|proxisome|proximity[- ]?(dependent|label)
                         (BioID / TurboID / APEX: biotin ligase fused to the bait labels neighbours)
affinity_purification  immunoprecipit|immuno-precipit|\bco-?ip\b|\bip-?ms\b|interactome|interacting prote|pull[- ]?down|precipitation of|affinity purif|anti-ha ms|associated complex|\binteractors\b|binding with|associates with|interact with|co-?purif|proximal interacting
                         (IP-MS / co-IP / tagged-bait affinity purification)
ptm_enrichment         phospho|acetylat|methylscan|methylarginine|lactylation|malonylation|nitrosylation|ubiquitylat|ubiquitinat|succinylat|crotonylat|glycosylat|\bptm\b|post-?translational modification
                         (the sample is an enriched modified subset, not the proteome)
fractionation          hyperlopit|\blopit\b|subcellular|fractionat|spatial proteom|organellar|sucrose gradient|\bbn-?page\b
                         (organelle or gradient fractionation: detection depends on the fraction)
whole_proteome         proteome|proteomic|protein abundance|expression|abundances following|antigen|immunogenic|proteogenomic
                         (whole-cell or expression proteome)

file-name cues (weaker, used only when the title gave none):
proximity_labelling    bioid|bira|turboid|_bio[_.-]|biotin
affinity_purification  (^|[_.-])ip[_.-]|_ip\d|ha[_.-]?ip|flag[_.-]?ip|beads|pulldown|pull[_-]down|compip
hdx                    hdx
ptm_enrichment         phos|_pho\d|psn|_sn[_.]|ptm
```

Counts: affinity_purification 32, fractionation 2, other 7, proximity_labelling 9, ptm_enrichment 10, whole_proteome 11.

**35 of 71 classifications hinge on a judgement call** (confidence low, none or medium, or hand-checked): either the title gave no keyword and the call rests on deposited file names, or the title matched cues for two classes. They are listed in `project_classification.tsv`; every one should be checked by a human before the table is relied on.

Hand overrides, each with its reason:

* `PXD039434` -> proximity_labelling: run names are BIO vs VEH (biotin against vehicle), i.e. proximity labelling, not the IP arm of the same paper
* `PXD068413` -> fractionation: runs are Top/Middle/Bottom gradient slices of a conoid preparation, not a bait pull-down
* `PXD068228` -> fractionation: runs are Top/Middle/Bottom gradient slices, same preparation style
* `PXD079900` -> affinity_purification: glycan-bead capture: a real affinity purification but the bait is a glycan, so it has no bait gene
* `PXD035658` -> whole_proteome: run names carry SILAC labelling at two drug concentrations (100uM_L / 10uM_H), i.e. a whole-proteome drug comparison, not a pull-down
* `PXD064338` -> whole_proteome: 32 TMT runs of infected host cells (a host immunity study): a whole-proteome experiment

## Exclusions

* **HDX-MS, `NOT_COMPARABLE.json`**: PXD023791, PXD031816. Hydrogen-deuterium exchange peptides carry an unknown number of deuteriums, so their masses are shifted by amounts the search does not model; any protein-level quantity reflects deuteration, not abundance (instruction 69 measured ~1/10 the PSMs per run of an ordinary bottom-up project). Spectra are kept; tables are not used.
* **`unsupported.json`**: PXD009604, PXD011816, PXD011828, PXD011829, PXD011830, PXD011831, PXD034607, PXD039335, PXD052664, PXD054994, PXD063409. The pipeline could not search these at all — SCIEX `.wiff` (needs a converter that does not run here), deposits with no spectra, and the two HDX projects stopped before searching (PXD034607, PXD039335).
* Projects classified as whole proteome, PTM enrichment, fractionation or other are **classified but not counted in the background**: a whole-proteome run detects nearly everything and would make every protein a frequent flyer, and a PTM enrichment detects a biased subset.
* 130 of the 203 selected projects have not finished the queue yet.

## Detection rule and its sensitivity

> A protein is **detected** in an experiment when at least **2 peptide-spectrum matches at 1 % spectrum-level FDR** (Sage `spectrum_q <= 0.01`, target hits only) map to that gene **and to no other parasite gene**.

1 % FDR is Sage's own number and the field standard. Two matches rather than one because at 1 % FDR the single-spectrum identifications are the ones most likely to be wrong, and for a background table a false detection is the costlier error: it inflates a protein's frequency and gets someone's real hit filtered out. Unique mapping because a peptide shared between paralogues would otherwise inflate whole families; a protein seen only through shared peptides is not counted, which makes frequencies slightly too low rather than too high.

| rule | proteins | frequency >= 0.5 | frequency >= 0.8 | seen once | median frequency | Spearman vs default |
|---|---|---|---|---|---|---|
| loose (>=1 match) | 5426 | 601 | 113 | 944 | 0.1220 | 0.9539 |
| default (>=2 matches) | 4528 | 437 | 81 | 927 | 0.0976 | 1.0000 |
| strict (>=3 matches, >=2 peptides) | 3922 | 329 | 55 | 871 | 0.0976 | 0.9560 |

`background.tsv` carries `frequency_strict` and `frequency_loose` beside the default, so any conclusion can be repeated under either.

**A known consequence of the unique-mapping rule.** Two nearly identical paralogues share most of their tryptic peptides, so both lose the shared ones and both look rarer than they are. Enolase 1 and enolase 2 are the clearest case in this table: they are among the most abundant proteins in a tachyzoite yet come out at a low frequency. Treat the frequency of a member of a close paralogue pair as a lower bound, and check `total_psm` beside it.

**Experiments, not runs.** One project is one experiment; a project's several runs are usually the bait and its control in replicate, which is one look at the background. Run-level counts are kept in `runs_detected` for comparison.

**Abundance.** Each protein is expressed as parts per million of its own project's total parasite signal, so projects of different depth can be averaged. `median_abundance` / `max_abundance` are MS1 label-free quantities (from `lfq.tsv`); `median_psm_ppm` / `max_psm_ppm` are summed MS2 fragment intensities and exist for every project, including the timsTOF ones that have no MS1 LFQ.

## Baits

41 affinity + proximity experiments; **26 have a confident bait gene** -- a symbol in the title that resolves to exactly one ME49 gene through ToxoDB product names, or a bait identified in the depositing paper or by its ToxoDB gene name (`parascape/data/pulldown_baits_curated.tsv`, each row with its source and the sentence that names the tagged gene) --, covering **22 distinct bait genes**.

* confident: 26
* a complex or a non-protein bait, not one tagged gene: 6
* none: 5
* ambiguous (unresolved): 1
* host bait: 1
* ambiguous (paper names the gene; the pull-down itself not verified): 1
* ambiguous (symbol match only (paper not readable)): 1

`interaction_studies.parquet` (97 author-reported studies, 31 with a bait gene) is a **separate** resource and was used only as a comparison: it is keyed by PMID and the PRIDE selection carries no PMID, so the two cannot be joined. Every bait in `project_classification.tsv` comes from the PRIDE title or run names, never from that table.

## Specificity score

```
specificity = log2( (ppm_here + 1) / (median background ppm + 1) )
            + log2( 1 / ((experiments_detected + 0.5) / (n_experiments + 0.5)) )
```

The first term asks whether there is more of the protein here than usual; the second asks whether it is unusual to see it at all. Both are reported separately in `specificity_example.tsv` (`log2_fold`, `log2_rarity`), because they fail differently: the fold is unreliable for a protein with few background measurements and the rarity term is unreliable for everything while `n_experiments` is small. A protein seen in every pull-down at its usual level scores about 0.

Worked example: **PXD016383**, the MYR1 immunoprecipitation (363 proteins in the hit list). Known MYR1-complex members and where the score put them:

| protein | gene | rank | percentile | specificity | background frequency |
|---|---|---|---|---|---|
| MYR1 | TGME49_254470 | 1 | 100.0 | 10.67 | 0.49 |
| MYR4 | TGME49_211460 | 2 | 99.7 | 8.84 | 0.25 |
| MYR3 | TGME49_237230 | 38 | 89.8 | 5.17 | 0.30 |
| GRA44 | TGME49_228170 | 52 | 86.0 | 4.82 | 0.71 |
| GRA45 | TGME49_316250 | 70 | 81.0 | 4.64 | 0.45 |
| ROP17 | TGME49_258580 | 183 | 49.9 | 1.78 | 0.83 |

## Validation

| | n | median frequency | max frequency |
|---|---|---|---|
| known complex members | 13 | 0.683 | 0.854 |
| expected frequent flyers | 82 | 0.805 | 0.951 |
| all proteins in the table | 4528 | 0.098 | 0.951 |

Separation: AUC 0.710, Mann-Whitney p = 0.00765. At a 50 % frequency cutoff, 11 of 13 complex members would be discarded as background and 73 of 82 flyers correctly caught.

The direction is right — the frequent flyers are above the complex members — but the two distributions overlap almost completely, and that is what decides whether a filter can be built. The cutoff sweep:

| cutoff | frequent flyers removed | complex members lost | proteins removed from the table |
|---|---|---|---|
| 0.25 | 78 / 82 (95 %) | 13 / 13 (100 %) | 1129 |
| 0.50 | 73 / 82 (89 %) | 11 / 13 (85 %) | 437 |
| 0.60 | 70 / 82 (85 %) | 8 / 13 (62 %) | 286 |
| 0.70 | 67 / 82 (82 %) | 5 / 13 (38 %) | 167 |
| 0.75 | 60 / 82 (73 %) | 5 / 13 (38 %) | 124 |
| 0.80 | 47 / 82 (57 %) | 4 / 13 (31 %) | 81 |
| 0.85 | 24 / 82 (29 %) | 1 / 13 (8 %) | 38 |
| 0.90 | 9 / 82 (11 %) | 0 / 13 (0 %) | 15 |
| 0.95 | 2 / 82 (2 %) | 0 / 13 (0 %) | 3 |
| 1.00 | 0 / 82 (0 %) | 0 / 13 (0 %) | 0 |

**No cutoff removes at least half the frequent flyers while losing at most 15 % of the known complex members, so the frequency column must not be used as a standalone filter.** The reason is biological rather than a bug: the glideosome, AMA1-RON, MIC2-M2AP and the MORC complex are among the most abundant proteins in a tachyzoite, so a pull-down from whole parasite lysate picks them up whatever the bait. Their frequencies (0.38-0.88) sit inside the frequent-flyer range, and a cutoff anywhere in it deletes exactly the complexes a *Toxoplasma* biologist cares about.

### What the full score does, which the frequency alone does not

The product is the specificity score, not the frequency column. Tested on the one case where the answer is certain — the tagged protein in its own pull-down — across every experiment with a confident bait:

| project | class | bait | gene | hits | bait rank | percentile | specificity | background frequency |
|---|---|---|---|---|---|---|---|---|
| PXD012103 | affinity purification | HCE1 | TGME49_239010 | 334 | 1 | 100.0 | 9.32 | 0.35 |
| PXD016383 | affinity purification | MYR1 | TGME49_254470 | 363 | 1 | 100.0 | 10.67 | 0.49 |
| PXD038848 | affinity purification | TgBDP1 | TGME49_263580 | 244 | 1 | 100.0 | 10.84 | 0.40 |
| PXD043808 | affinity purification | AP2XI-2 | TGME49_310900 | 710 | 1 | 100.0 | 10.35 | 0.25 |
| PXD043994 | affinity purification | Cyclin-4 (TgCyc4) | TGME49_249880 | 1085 | 1 | 100.0 | 11.38 | 0.08 |
| PXD067653 | affinity purification | GRA35 | TGME49_226380 | 33 | 1 | 100.0 | 12.49 | 0.64 |
| PXD051495 | affinity purification | PP2A-C2 | TGME49_215170 | 704 | 1 | 100.0 | 10.25 | 0.16 |
| PXD043924 | affinity purification | Cdt1 (TgiRD1) | TGME49_247040 | 1774 | 1 | 100.0 | 8.68 | 0.13 |
| PXD083501 | affinity purification | MyoL | TGME49_291020 | 1829 | 2 | 99.9 | 9.18 | 0.11 |
| PXD039896 | affinity purification | C3 | TGME49_286470 | 1167 | 6 | 99.6 | 8.03 | 0.23 |
| PXD059579 | proximity labelling | GSK3 | TGME49_265330 | 1748 | 9 | 99.5 | 6.75 | 0.42 |
| PXD016846 | affinity purification | MORC | TGME49_305340 | 808 | 6 | 99.4 | 10.57 | 0.69 |
| PXD063615 | affinity purification | PDIA3 | TGME49_211680 | 698 | 6 | 99.3 | 6.86 | 0.69 |
| PXD071065 | affinity purification | TgCoqFMO | TGME49_232540 | 2084 | 17 | 99.2 | 7.29 | 0.11 |
| PXD032360 | affinity purification | GRA28 | TGME49_231960 | 428 | 7 | 98.6 | 7.02 | 0.49 |
| PXD039919 | affinity purification | SPARK | TGME49_268210 | 749 | 12 | 98.5 | 6.21 | 0.16 |
| PXD044091 | affinity purification | CLAMP | TGME49_212270 | 731 | 13 | 98.4 | 4.99 | 0.16 |
| PXD037318 | affinity purification | TgH1-like | TGME49_315570 | 65 | 4 | 95.4 | 7.59 | 0.25 |
| PXD051495 | affinity purification | PP2A-A2 | TGME49_229350 | 704 | 37 | 94.9 | 5.05 | 0.06 |
| PXD039979 | affinity purification | SPARKEL | TGME49_291010 | 1141 | 233 | 79.7 | 2.59 | 0.13 |
| PXD039983 | proximity labelling | SPARK | TGME49_268210 | 1092 | 251 | 77.1 | 3.15 | 0.16 |
| PXD039986 | proximity labelling | SPARKEL | TGME49_291010 | 1398 | 323 | 77.0 | 2.92 | 0.13 |
| PXD032319 | affinity purification | ROP1 | TGME49_309590 | 38 | 13 | 68.4 | 4.16 | 0.66 |
| PXD039985 | proximity labelling | SPARKEL | TGME49_291010 | 974 | 502 | 48.6 | 0.34 | 0.13 |
| PXD039922 | affinity purification | SPARK | TGME49_268210 | 2261 | 1298 | 42.6 | 2.55 | 0.16 |
| PXD044080 | affinity purification | TGME49_306920 | TGME49_306920 | 1229 | 868 | 29.5 | 1.23 | 0.08 |

Bait in the top 5 % of its own hit list: **18 of 26** overall and **17 of 22** of the affinity purifications; median bait rank 6 of a median 778 proteins, against median rank 28 if the list is ranked by raw abundance instead. `bait_recovery.tsv`.

The split by class matters. The score is built for affinity purification, where the tagged protein is pulled out intact and should dominate its own eluate. In proximity labelling the bait is a biotin-ligase fusion whose job is to label its *neighbours*, so there is no reason for it to top the list, and in a degron experiment (`PXD039922`, SPARK-mNG-mAID) the bait is deliberately degraded, so a low rank there is the correct answer rather than a failure. Those cases are kept in the table rather than dropped, because dropping them would be choosing the evidence.

No experiment here baits the glideosome, AMA1-RON, MIC2-M2AP or MORC/HDAC3, so none of them has a pull-down in which it *should* be specific. The best specificity each member reaches in any of the 41 experiments is in `complex_best_specificity.tsv`; the highest per complex:

| complex | protein | best in | rank | percentile | specificity |
|---|---|---|---|---|---|
| AMA1-RON | AMA1 | PXD011692 | 22 | 65.0 | 4.84 |
| AMA1-RON | RON2 | PXD043994 | 132 | 87.9 | 3.73 |
| MIC2-M2AP | MIC2 | PXD032360 | 91 | 79.0 | 4.72 |
| MIC2-M2AP | M2AP | PXD043924 | 281 | 84.2 | 3.45 |
| MORC/HDAC3 | MORC | PXD016846 | 6 | 99.4 | 10.57 |
| MORC/HDAC3 | HDAC3 | PXD016846 | 45 | 94.6 | 5.21 |
| MyoA glideosome | GAP50 | PXD012103 | 6 | 98.5 | 5.68 |
| MyoA glideosome | GAP45 | PXD012103 | 7 | 98.2 | 5.58 |

Per complex:

| complex | protein | gene | experiments detected | frequency |
|---|---|---|---|---|
| AMA1-RON | AMA1 | TGME49_255260 | 31 / 41 | 0.756 |
| AMA1-RON | RON2 | TGME49_300100 | 33 / 41 | 0.805 |
| AMA1-RON | RON4 | TGME49_229010 | 35 / 41 | 0.854 |
| MIC2-M2AP | M2AP | TGME49_214940 | 24 / 41 | 0.585 |
| MIC2-M2AP | MIC2 | TGME49_201780 | 23 / 41 | 0.561 |
| MORC/HDAC3 | HDAC3 | TGME49_227290 | 20 / 41 | 0.488 |
| MORC/HDAC3 | MORC | TGME49_305340 | 28 / 41 | 0.683 |
| MyoA glideosome | ELC1 | TGME49_269442 | 14 / 41 | 0.341 |
| MyoA glideosome | GAP40 | TGME49_249850 | 24 / 41 | 0.585 |
| MyoA glideosome | GAP45 | TGME49_223940 | 27 / 41 | 0.659 |
| MyoA glideosome | GAP50 | TGME49_219320 | 33 / 41 | 0.805 |
| MyoA glideosome | MLC1 | TGME49_257680 | 28 / 41 | 0.683 |
| MyoA glideosome | MyoA | TGME49_235470 | 34 / 41 | 0.829 |

Per frequent-flyer family (curated sets; the cytosolic ribosome is the hyperLOPIT 40S/60S marker set, so the apicoplast and mitochondrial ribosomes are not mixed in):

| family | proteins | median frequency | max frequency |
|---|---|---|---|
| GAPDH | 2 | 0.780 | 0.829 |
| HSP70/HSP90 | 3 | 0.878 | 0.927 |
| actin | 1 | 0.951 | 0.951 |
| cytosolic ribosome (hyperLOPIT marker) | 68 | 0.805 | 0.927 |
| enolase | 2 | 0.244 | 0.488 |
| tubulin | 6 | 0.598 | 0.951 |

The same families as a keyword would define them, for comparison only — these sweep in the tubulin-tyrosine ligases and the organellar ribosomes, which are neither abundant nor sticky:

| family (keyword) | proteins | median frequency | max frequency |
|---|---|---|---|
| actin (any) | 12 | 0.244 | 0.951 |
| heat shock (any) | 23 | 0.195 | 0.927 |
| ribosomal protein (any) | 158 | 0.244 | 0.927 |
| tubulin (any) | 23 | 0.000 | 0.951 |

## Is it usable yet

* **The frequency column as a standalone filter: not yet, and this is the headline.** The frequent flyers are above the complex members (median 0.80 against 0.68, AUC 0.710, p = 0.00765) but the overlap is almost total: at the 50 % cutoff, 11 of the 13 known complex members are deleted. Do not filter a hit list on frequency.
* **The specificity score: usable for ranking, with the caveats below.** It puts the tagged bait in the top 5 % of its own hit list in 17 of the 22 affinity purifications with a confident bait, and on the MYR1 worked example it puts MYR1 first, MYR4 third and five of six known complex members in the top fifth. The work is done by the abundance-fold term, not the frequency term.
* **The limiting resource is independent baits (22 distinct, over 41 experiments), not proteins or spectra.** Nearly every experiment is a pull-down from whole parasite lysate, so the abundant invasion and motility machinery is in every list. A frequency becomes informative only when there are enough baits from enough different compartments that "present in most pull-downs" means stickiness rather than abundance. For comparison, the human background frequency (Mellacheruvu et al. 2013) was built from more than 400 experiments.
* **Practical advice for now:** rank a hit list by `specificity`, or by its `log2_fold` component alone, and read `frequency` as context — "this protein has been seen in 35 of 41 other pull-downs" — rather than as a filter. The `median_abundance` column is the useful part of the table today: it says how much of a protein a pull-down normally gives, which is what makes an enrichment judgeable.

## What will change as the queue fills

* **The frequency resolution.** With 41 experiments the smallest step a frequency can take is 1/41 = 0.024. Every frequency in this table is that coarse.
* **Bait independence.** 22 distinct confident baits is the binding constraint, not the protein count. Several of the finished projects come from the same laboratories and the same complexes (SPARK/SPARKEL, CDPK1, the mitochondrial complexes), so the experiments are not as independent as their number suggests. A protein that happens to sit in two of the few baited complexes looks like a frequent flyer for a reason that has nothing to do with stickiness.
* **Classification.** The 29 uncertain calls should be reviewed; a misclassified whole-proteome project in the background would raise every frequency at once.
* **Nothing here is hand-patched.** Re-run the builder and the whole table, the validation and this README are rebuilt from whatever has finished.

## Files

* `project_classification.tsv` — every finished project: class, how it was decided, the matched cue, confidence, bait and bait confidence, run and PSM counts, title.
* `background.tsv` — gene_id, experiments_detected, n_experiments, frequency, median_abundance, max_abundance, plus the strict and loose frequencies, the PSM-intensity abundances, run counts, product name and hyperLOPIT compartment.
* `specificity_example.tsv` — the MYR1 hit list scored, with both score components.
* `validation.tsv` — the known complexes and the frequent flyers (curated and keyword-defined) with their frequencies.
* `bait_recovery.tsv` — for every experiment with a confident bait, where the score put the bait in its own hit list.
* `complex_best_specificity.tsv` — the best specificity each known complex member reaches in any experiment, and where.
* `coverage.tsv` — the counts above, so a later run can be compared with this one.
