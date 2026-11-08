# Toxoplasma background frequency, first pass (2026-11-08)

Instruction 71 product 4: how often each parasite protein is detected in an affinity-purification or proximity-labelling experiment regardless of what was pulled down, plus a specificity score for a new hit list.

Notebook: `notebooks/pulldown_background_2026_11_08.ipynb` (executed). Builder: `scripts/build_pulldown_background_notebook.py`. Re-running the builder rescans the uniform search output and picks up whatever has finished since.

## Coverage

| | n |
|---|---|
| projects selected for the uniform re-search | 203 |
| project directories on disk | 111 |
| finished (done.json receipt) | 97 |
| excluded: unsupported.json | 12 |
| excluded: NOT_COMPARABLE.json (HDX-MS) | 2 |
| finished and classified | 94 |
| still waiting in the queue | 106 |
| experiments in the background (affinity + proximity) | 50 |
| of those with a confident bait gene | 29 |
| distinct confident bait genes | 25 |
| runs behind those experiments | 440 |
| PSM rows read | 15956399 |
| PSMs at 1% FDR, target | 3975917 |
| parasite proteins in the background table | 4952 |
| class: affinity_purification | 40 |
| class: fractionation | 2 |
| class: other | 18 |
| class: proximity_labelling | 10 |
| class: ptm_enrichment | 11 |
| class: whole_proteome | 13 |

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

Counts: affinity_purification 40, fractionation 2, other 18, proximity_labelling 10, ptm_enrichment 11, whole_proteome 13.

**50 of 94 classifications hinge on a judgement call** (confidence low, none or medium, or hand-checked): either the title gave no keyword and the call rests on deposited file names, or the title matched cues for two classes. They are listed in `project_classification.tsv`; every one should be checked by a human before the table is relied on.

Hand overrides, each with its reason:

* `PXD039434` -> proximity_labelling: run names are BIO vs VEH (biotin against vehicle), i.e. proximity labelling, not the IP arm of the same paper
* `PXD068413` -> fractionation: runs are Top/Middle/Bottom gradient slices of a conoid preparation, not a bait pull-down
* `PXD068228` -> fractionation: runs are Top/Middle/Bottom gradient slices, same preparation style
* `PXD079900` -> affinity_purification: glycan-bead capture: a real affinity purification but the bait is a glycan, so it has no bait gene
* `PXD035658` -> whole_proteome: run names carry SILAC labelling at two drug concentrations (100uM_L / 10uM_H), i.e. a whole-proteome drug comparison, not a pull-down
* `PXD064338` -> whole_proteome: 32 TMT runs of infected host cells (a host immunity study): a whole-proteome experiment

## Exclusions

* **HDX-MS, `NOT_COMPARABLE.json`**: PXD023791, PXD031816. Hydrogen-deuterium exchange peptides carry an unknown number of deuteriums, so their masses are shifted by amounts the search does not model; any protein-level quantity reflects deuteration, not abundance (instruction 69 measured ~1/10 the PSMs per run of an ordinary bottom-up project). Spectra are kept; tables are not used.
* **`unsupported.json`**: PXD003399, PXD009604, PXD011816, PXD011828, PXD011829, PXD011830, PXD011831, PXD034607, PXD039335, PXD052664, PXD054994, PXD063409. The pipeline could not search these at all — SCIEX `.wiff` (needs a converter that does not run here), deposits with no spectra, and the two HDX projects stopped before searching (PXD034607, PXD039335).
* Projects classified as whole proteome, PTM enrichment, fractionation or other are **classified but not counted in the background**: a whole-proteome run detects nearly everything and would make every protein a frequent flyer, and a PTM enrichment detects a biased subset.
* 106 of the 203 selected projects have not finished the queue yet.

## Detection rule and its sensitivity

> A protein is **detected** in an experiment when at least **2 peptide-spectrum matches at 1 % spectrum-level FDR** (Sage `spectrum_q <= 0.01`, target hits only) map to that gene **and to no other parasite gene**.

1 % FDR is Sage's own number and the field standard. Two matches rather than one because at 1 % FDR the single-spectrum identifications are the ones most likely to be wrong, and for a background table a false detection is the costlier error: it inflates a protein's frequency and gets someone's real hit filtered out. Unique mapping because a peptide shared between paralogues would otherwise inflate whole families; a protein seen only through shared peptides is not counted, which makes frequencies slightly too low rather than too high.

| rule | proteins | frequency >= 0.5 | frequency >= 0.8 | seen once | median frequency | Spearman vs default |
|---|---|---|---|---|---|---|
| loose (>=1 match) | 5775 | 683 | 127 | 807 | 0.1400 | 0.9630 |
| default (>=2 matches) | 4952 | 495 | 90 | 832 | 0.1000 | 1.0000 |
| strict (>=3 matches, >=2 peptides) | 4290 | 385 | 62 | 774 | 0.1000 | 0.9635 |

`background.tsv` carries `frequency_strict` and `frequency_loose` beside the default, so any conclusion can be repeated under either.

**A known consequence of the unique-mapping rule.** Two nearly identical paralogues share most of their tryptic peptides, so both lose the shared ones and both look rarer than they are. Enolase 1 and enolase 2 are the clearest case in this table: they are among the most abundant proteins in a tachyzoite yet come out at a low frequency. Treat the frequency of a member of a close paralogue pair as a lower bound, and check `total_psm` beside it.

**Experiments, not runs.** One project is one experiment; a project's several runs are usually the bait and its control in replicate, which is one look at the background. Run-level counts are kept in `runs_detected` for comparison.

**Abundance.** Each protein is expressed as parts per million of its own project's total parasite signal, so projects of different depth can be averaged. `median_abundance` / `max_abundance` are MS1 label-free quantities (from `lfq.tsv`); `median_psm_ppm` / `max_psm_ppm` are summed MS2 fragment intensities and exist for every project, including the timsTOF ones that have no MS1 LFQ.

## Baits

50 affinity + proximity experiments; **29 have a confident bait gene** -- a symbol in the title that resolves to exactly one ME49 gene through ToxoDB product names, or a bait identified in the depositing paper or by its ToxoDB gene name (`parascape/data/pulldown_baits_curated.tsv`, each row with its source and the sentence that names the tagged gene) -- covering **25 distinct bait genes**.

* confident: 29
* none: 8
* a complex or a non-protein bait, not one tagged gene: 6
* ambiguous (symbol not in any ToxoDB product name): 2
* ambiguous (unresolved): 1
* host bait: 1
* ambiguous (symbol matches several genes): 1
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
| MYR1 | TGME49_254470 | 1 | 100.0 | 10.83 | 0.52 |
| MYR4 | TGME49_211460 | 2 | 99.7 | 9.19 | 0.27 |
| MYR3 | TGME49_237230 | 19 | 95.0 | 5.98 | 0.33 |
| GRA44 | TGME49_228170 | 69 | 81.3 | 4.81 | 0.72 |
| GRA45 | TGME49_316250 | 73 | 80.2 | 4.64 | 0.45 |
| ROP17 | TGME49_258580 | 188 | 48.5 | 1.80 | 0.82 |

## Validation

| | n | median frequency | max frequency |
|---|---|---|---|
| known complex members | 13 | 0.700 | 0.860 |
| expected frequent flyers | 82 | 0.800 | 0.960 |
| all proteins in the table | 4952 | 0.100 | 0.960 |

Separation: AUC 0.689, Mann-Whitney p = 0.0145. At a 50 % frequency cutoff, 12 of 13 complex members would be discarded as background and 74 of 82 flyers correctly caught.

The direction is right — the frequent flyers are above the complex members — but the two distributions overlap almost completely, and that is what decides whether a filter can be built. The cutoff sweep:

| cutoff | frequent flyers removed | complex members lost | proteins removed from the table |
|---|---|---|---|
| 0.25 | 79 / 82 (96 %) | 13 / 13 (100 %) | 1299 |
| 0.50 | 74 / 82 (90 %) | 12 / 13 (92 %) | 495 |
| 0.60 | 70 / 82 (85 %) | 8 / 13 (62 %) | 286 |
| 0.70 | 63 / 82 (77 %) | 6 / 13 (46 %) | 156 |
| 0.75 | 58 / 82 (71 %) | 5 / 13 (38 %) | 122 |
| 0.80 | 50 / 82 (61 %) | 4 / 13 (31 %) | 90 |
| 0.85 | 19 / 82 (23 %) | 1 / 13 (8 %) | 33 |
| 0.90 | 8 / 82 (10 %) | 0 / 13 (0 %) | 13 |
| 0.95 | 2 / 82 (2 %) | 0 / 13 (0 %) | 2 |
| 1.00 | 0 / 82 (0 %) | 0 / 13 (0 %) | 0 |

**No cutoff removes at least half the frequent flyers while losing at most 15 % of the known complex members, so the frequency column must not be used as a standalone filter.** The reason is biological rather than a bug: the glideosome, AMA1-RON, MIC2-M2AP and the MORC complex are among the most abundant proteins in a tachyzoite, so a pull-down from whole parasite lysate picks them up whatever the bait. Their frequencies (0.38-0.88) sit inside the frequent-flyer range, and a cutoff anywhere in it deletes exactly the complexes a *Toxoplasma* biologist cares about.

### What the full score does, which the frequency alone does not

The product is the specificity score, not the frequency column. Tested on the one case where the answer is certain — the tagged protein in its own pull-down — across every experiment with a confident bait:

| project | class | bait | gene | hits | bait rank | percentile | specificity | background frequency |
|---|---|---|---|---|---|---|---|---|
| PXD012103 | affinity purification | HCE1 | TGME49_239010 | 334 | 1 | 100.0 | 9.61 | 0.35 |
| PXD016383 | affinity purification | MYR1 | TGME49_254470 | 363 | 1 | 100.0 | 10.83 | 0.52 |
| PXD083501 | affinity purification | MyoL | TGME49_291020 | 1829 | 1 | 100.0 | 9.26 | 0.11 |
| PXD043808 | affinity purification | AP2XI-2 | TGME49_310900 | 710 | 1 | 100.0 | 9.98 | 0.31 |
| PXD043924 | affinity purification | Cdt1 (TgiRD1) | TGME49_247040 | 1774 | 1 | 100.0 | 8.52 | 0.15 |
| PXD051495 | affinity purification | PP2A-A2 | TGME49_229350 | 704 | 1 | 100.0 | 12.90 | 0.07 |
| PXD067653 | affinity purification | GRA35 | TGME49_226380 | 33 | 1 | 100.0 | 12.44 | 0.66 |
| PXD071065 | affinity purification | TgCoqFMO | TGME49_232540 | 2084 | 1 | 100.0 | 9.70 | 0.11 |
| PXD051495 | affinity purification | PP2A-C2 | TGME49_215170 | 704 | 2 | 99.9 | 10.10 | 0.23 |
| PXD038848 | affinity purification | TgBDP1 | TGME49_263580 | 244 | 2 | 99.6 | 10.45 | 0.41 |
| PXD059579 | proximity labelling | GSK3 | TGME49_265330 | 1748 | 11 | 99.4 | 6.63 | 0.43 |
| PXD032362 | affinity purification | BCC7 | TGME49_311230 | 491 | 4 | 99.4 | 6.94 | 0.58 |
| PXD043994 | affinity purification | Cyclin-4 (TgCyc4) | TGME49_249880 | 1085 | 8 | 99.4 | 6.97 | 0.09 |
| PXD032360 | affinity purification | GRA28 | TGME49_231960 | 428 | 4 | 99.3 | 7.28 | 0.52 |
| PXD039896 | affinity purification | C3 | TGME49_286470 | 1167 | 11 | 99.1 | 7.80 | 0.27 |
| PXD063615 | affinity purification | PDIA3 | TGME49_211680 | 698 | 7 | 99.1 | 6.90 | 0.70 |
| PXD016846 | affinity purification | MORC | TGME49_305340 | 808 | 9 | 99.0 | 10.48 | 0.70 |
| PXD066820 | affinity purification | TGME49_320300 | TGME49_320300 | 1688 | 27 | 98.5 | 6.84 | 0.27 |
| PXD039919 | affinity purification | SPARK | TGME49_268210 | 749 | 15 | 98.1 | 6.49 | 0.17 |
| PXD044592 | affinity purification | Crk4 | TGME49_256070 | 1016 | 33 | 96.9 | 4.71 | 0.09 |
| PXD037318 | affinity purification | TgH1-like | TGME49_315570 | 65 | 4 | 95.4 | 7.08 | 0.25 |
| PXD044091 | affinity purification | CLAMP | TGME49_212270 | 731 | 42 | 94.4 | 4.89 | 0.17 |
| PXD039986 | proximity labelling | SPARKEL | TGME49_291010 | 1398 | 271 | 80.7 | 3.20 | 0.11 |
| PXD039979 | affinity purification | SPARKEL | TGME49_291010 | 1141 | 222 | 80.6 | 2.88 | 0.11 |
| PXD032319 | affinity purification | ROP1 | TGME49_309590 | 38 | 9 | 78.9 | 5.14 | 0.66 |
| PXD039983 | proximity labelling | SPARK | TGME49_268210 | 1092 | 237 | 78.4 | 3.43 | 0.17 |
| PXD044080 | affinity purification | TGME49_306920 | TGME49_306920 | 1229 | 526 | 57.3 | 3.20 | 0.11 |
| PXD039985 | proximity labelling | SPARKEL | TGME49_291010 | 974 | 449 | 54.0 | 0.62 | 0.11 |
| PXD039922 | affinity purification | SPARK | TGME49_268210 | 2261 | 1204 | 46.8 | 2.83 | 0.17 |

Bait in the top 5 % of its own hit list: **21 of 29** overall and **20 of 25** of the affinity purifications; median bait rank 8 of a median 808 proteins, against median rank 31 if the list is ranked by raw abundance instead. `bait_recovery.tsv`.

The split by class matters. The score is built for affinity purification, where the tagged protein is pulled out intact and should dominate its own eluate. In proximity labelling the bait is a biotin-ligase fusion whose job is to label its *neighbours*, so there is no reason for it to top the list, and in a degron experiment (`PXD039922`, SPARK-mNG-mAID) the bait is deliberately degraded, so a low rank there is the correct answer rather than a failure. Those cases are kept in the table rather than dropped, because dropping them would be choosing the evidence.

No experiment here baits the glideosome, AMA1-RON, MIC2-M2AP or MORC/HDAC3, so none of them has a pull-down in which it *should* be specific. The best specificity each member reaches in any of the 50 experiments is in `complex_best_specificity.tsv`; the highest per complex:

| complex | protein | best in | rank | percentile | specificity |
|---|---|---|---|---|---|
| AMA1-RON | AMA1 | PXD011692 | 26 | 58.3 | 4.64 |
| AMA1-RON | RON2 | PXD043994 | 127 | 88.4 | 3.81 |
| MIC2-M2AP | MIC2 | PXD032360 | 72 | 83.4 | 4.92 |
| MIC2-M2AP | M2AP | PXD043924 | 212 | 88.1 | 3.45 |
| MORC/HDAC3 | MORC | PXD016846 | 9 | 99.0 | 10.48 |
| MORC/HDAC3 | HDAC3 | PXD016846 | 47 | 94.3 | 5.28 |
| MyoA glideosome | GAP50 | PXD012103 | 9 | 97.6 | 5.75 |
| MyoA glideosome | GAP45 | PXD012103 | 11 | 97.0 | 5.58 |

Per complex:

| complex | protein | gene | experiments detected | frequency |
|---|---|---|---|---|
| AMA1-RON | AMA1 | TGME49_255260 | 38 / 50 | 0.760 |
| AMA1-RON | RON2 | TGME49_300100 | 41 / 50 | 0.820 |
| AMA1-RON | RON4 | TGME49_229010 | 42 / 50 | 0.840 |
| MIC2-M2AP | M2AP | TGME49_214940 | 29 / 50 | 0.580 |
| MIC2-M2AP | MIC2 | TGME49_201780 | 30 / 50 | 0.600 |
| MORC/HDAC3 | HDAC3 | TGME49_227290 | 25 / 50 | 0.500 |
| MORC/HDAC3 | MORC | TGME49_305340 | 35 / 50 | 0.700 |
| MyoA glideosome | ELC1 | TGME49_269442 | 20 / 50 | 0.400 |
| MyoA glideosome | GAP40 | TGME49_249850 | 30 / 50 | 0.600 |
| MyoA glideosome | GAP45 | TGME49_223940 | 33 / 50 | 0.660 |
| MyoA glideosome | GAP50 | TGME49_219320 | 40 / 50 | 0.800 |
| MyoA glideosome | MLC1 | TGME49_257680 | 36 / 50 | 0.720 |
| MyoA glideosome | MyoA | TGME49_235470 | 43 / 50 | 0.860 |

Per frequent-flyer family (curated sets; the cytosolic ribosome is the hyperLOPIT 40S/60S marker set, so the apicoplast and mitochondrial ribosomes are not mixed in):

| family | proteins | median frequency | max frequency |
|---|---|---|---|
| GAPDH | 2 | 0.810 | 0.840 |
| HSP70/HSP90 | 3 | 0.880 | 0.940 |
| actin | 1 | 0.960 | 0.960 |
| cytosolic ribosome (hyperLOPIT marker) | 68 | 0.800 | 0.940 |
| enolase | 2 | 0.250 | 0.500 |
| tubulin | 6 | 0.610 | 0.960 |

The same families as a keyword would define them, for comparison only — these sweep in the tubulin-tyrosine ligases and the organellar ribosomes, which are neither abundant nor sticky:

| family (keyword) | proteins | median frequency | max frequency |
|---|---|---|---|
| actin (any) | 12 | 0.290 | 0.960 |
| heat shock (any) | 23 | 0.240 | 0.940 |
| ribosomal protein (any) | 158 | 0.230 | 0.940 |
| tubulin (any) | 23 | 0.020 | 0.960 |

## Is it usable yet

* **The frequency column as a standalone filter: not yet, and this is the headline.** The frequent flyers are above the complex members (median 0.80 against 0.70, AUC 0.689, p = 0.0145) but the overlap is almost total: at the 50 % cutoff, 12 of the 13 known complex members are deleted. Do not filter a hit list on frequency.
* **The specificity score: usable for ranking, with the caveats below.** It puts the tagged bait in the top 5 % of its own hit list in 20 of the 25 affinity purifications with a confident bait, and on the MYR1 worked example it puts MYR1 first, MYR4 third and five of six known complex members in the top fifth. The work is done by the abundance-fold term, not the frequency term.
* **The limiting resource is independent baits (25 distinct, over 50 experiments), not proteins or spectra.** Nearly every experiment is a pull-down from whole parasite lysate, so the abundant invasion and motility machinery is in every list. A frequency becomes informative only when there are enough baits from enough different compartments that "present in most pull-downs" means stickiness rather than abundance. For comparison, the human background frequency (Mellacheruvu et al. 2013) was built from more than 400 experiments.
* **Practical advice for now:** rank a hit list by `specificity`, or by its `log2_fold` component alone, and read `frequency` as context — "this protein has been seen in 43 of 50 other pull-downs" — rather than as a filter. The `median_abundance` column is the useful part of the table today: it says how much of a protein a pull-down normally gives, which is what makes an enrichment judgeable.

## What will change as the queue fills

* **The frequency resolution.** With 50 experiments the smallest step a frequency can take is 1/50 = 0.020. Every frequency in this table is that coarse.
* **Bait independence.** 25 distinct confident baits is the binding constraint, not the protein count. Several of the finished projects come from the same laboratories and the same complexes (SPARK/SPARKEL, CDPK1, the mitochondrial complexes), so the experiments are not as independent as their number suggests. A protein that happens to sit in two of the few baited complexes looks like a frequent flyer for a reason that has nothing to do with stickiness.
* **Classification.** The 44 uncertain calls should be reviewed; a misclassified whole-proteome project in the background would raise every frequency at once.
* **Nothing here is hand-patched.** Re-run the builder and the whole table, the validation and this README are rebuilt from whatever has finished.

## Files

* `project_classification.tsv` — every finished project: class, how it was decided, the matched cue, confidence, bait and bait confidence, run and PSM counts, title.
* `background.tsv` — gene_id, experiments_detected, n_experiments, frequency, median_abundance, max_abundance, plus the strict and loose frequencies, the PSM-intensity abundances, run counts, product name and hyperLOPIT compartment.
* `specificity_example.tsv` — the MYR1 hit list scored, with both score components.
* `validation.tsv` — the known complexes and the frequent flyers (curated and keyword-defined) with their frequencies.
* `bait_recovery.tsv` — for every experiment with a confident bait, where the score put the bait in its own hit list.
* `complex_best_specificity.tsv` — the best specificity each known complex member reaches in any experiment, and where.
* `coverage.tsv` — the counts above, so a later run can be compared with this one.
