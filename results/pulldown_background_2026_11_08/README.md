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
| experiments in the background (affinity + proximity) | 60 |
| of those with a confident bait gene | 36 |
| distinct confident bait genes | 30 |
| runs behind those experiments | 585 |
| PSM rows read | 15956399 |
| PSMs at 1% FDR, target | 3975917 |
| parasite proteins in the background table | 5125 |
| class: affinity_purification | 48 |
| class: crosslinking | 1 |
| class: fractionation | 2 |
| class: other | 3 |
| class: proximity_labelling | 12 |
| class: ptm_enrichment | 13 |
| class: whole_proteome | 15 |

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

Counts: affinity_purification 48, crosslinking 1, fractionation 2, other 3, proximity_labelling 12, ptm_enrichment 13, whole_proteome 15.

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
| loose (>=1 match) | 6034 | 666 | 116 | 765 | 0.1333 | 0.9680 |
| default (>=2 matches) | 5125 | 495 | 82 | 716 | 0.1167 | 1.0000 |
| strict (>=3 matches, >=2 peptides) | 4460 | 379 | 59 | 680 | 0.1000 | 0.9661 |

`background.tsv` carries `frequency_strict` and `frequency_loose` beside the default, so any conclusion can be repeated under either.

**A known consequence of the unique-mapping rule.** Two nearly identical paralogues share most of their tryptic peptides, so both lose the shared ones and both look rarer than they are. Enolase 1 and enolase 2 are the clearest case in this table: they are among the most abundant proteins in a tachyzoite yet come out at a low frequency. Treat the frequency of a member of a close paralogue pair as a lower bound, and check `total_psm` beside it.

**Experiments, not runs.** One project is one experiment; a project's several runs are usually the bait and its control in replicate, which is one look at the background. Run-level counts are kept in `runs_detected` for comparison.

**Abundance.** Each protein is expressed as parts per million of its own project's total parasite signal, so projects of different depth can be averaged. `median_abundance` / `max_abundance` are MS1 label-free quantities (from `lfq.tsv`); `median_psm_ppm` / `max_psm_ppm` are summed MS2 fragment intensities and exist for every project, including the timsTOF ones that have no MS1 LFQ.

## Baits

60 affinity + proximity experiments; **36 have a confident bait gene** -- a symbol in the title that resolves to exactly one ME49 gene through ToxoDB product names, or a bait identified in the depositing paper or by its ToxoDB gene name (`parascape/data/pulldown_baits_curated.tsv`, each row with its source and the sentence that names the tagged gene) -- covering **30 distinct bait genes**.

* confident: 36
* none: 8
* a complex or a non-protein bait, not one tagged gene: 6
* ambiguous (symbol not in any ToxoDB product name): 2
* ambiguous (symbol only, no gene): 2
* ambiguous (unresolved): 1
* ambiguous (host bait (human GBP1)): 1
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
| MYR1 | TGME49_254470 | 1 | 100.0 | 10.74 | 0.55 |
| MYR4 | TGME49_211460 | 4 | 99.2 | 8.59 | 0.31 |
| MYR3 | TGME49_237230 | 30 | 92.0 | 5.69 | 0.32 |
| GRA44 | TGME49_228170 | 67 | 81.8 | 4.76 | 0.74 |
| GRA45 | TGME49_316250 | 80 | 78.2 | 4.48 | 0.47 |
| ROP17 | TGME49_258580 | 188 | 48.5 | 1.84 | 0.80 |

## Validation

| | n | median frequency | max frequency |
|---|---|---|---|
| known complex members | 13 | 0.683 | 0.850 |
| expected frequent flyers | 82 | 0.800 | 0.950 |
| all proteins in the table | 5125 | 0.117 | 0.950 |

Separation: AUC 0.689, Mann-Whitney p = 0.0144. At a 50 % frequency cutoff, 12 of 13 complex members would be discarded as background and 74 of 82 flyers correctly caught.

The direction is right — the frequent flyers are above the complex members — but the two distributions overlap almost completely, and that is what decides whether a filter can be built. The cutoff sweep:

| cutoff | frequent flyers removed | complex members lost | proteins removed from the table |
|---|---|---|---|
| 0.25 | 78 / 82 (95 %) | 13 / 13 (100 %) | 1369 |
| 0.50 | 74 / 82 (90 %) | 12 / 13 (92 %) | 495 |
| 0.60 | 71 / 82 (87 %) | 8 / 13 (62 %) | 284 |
| 0.70 | 61 / 82 (74 %) | 6 / 13 (46 %) | 160 |
| 0.75 | 56 / 82 (68 %) | 4 / 13 (31 %) | 116 |
| 0.80 | 47 / 82 (57 %) | 4 / 13 (31 %) | 82 |
| 0.85 | 17 / 82 (21 %) | 0 / 13 (0 %) | 27 |
| 0.90 | 8 / 82 (10 %) | 0 / 13 (0 %) | 13 |
| 0.95 | 0 / 82 (0 %) | 0 / 13 (0 %) | 0 |
| 1.00 | 0 / 82 (0 %) | 0 / 13 (0 %) | 0 |

**No cutoff removes at least half the frequent flyers while losing at most 15 % of the known complex members, so the frequency column must not be used as a standalone filter.** The reason is biological rather than a bug: the glideosome, AMA1-RON, MIC2-M2AP and the MORC complex are among the most abundant proteins in a tachyzoite, so a pull-down from whole parasite lysate picks them up whatever the bait. Their frequencies (0.38-0.88) sit inside the frequent-flyer range, and a cutoff anywhere in it deletes exactly the complexes a *Toxoplasma* biologist cares about.

### What the full score does, which the frequency alone does not

The product is the specificity score, not the frequency column. Tested on the one case where the answer is certain — the tagged protein in its own pull-down — across every experiment with a confident bait:

| project | class | bait | gene | hits | bait rank | percentile | specificity | background frequency |
|---|---|---|---|---|---|---|---|---|
| PXD012103 | affinity purification | HCE1 | TGME49_239010 | 334 | 1 | 100.0 | 9.76 | 0.34 |
| PXD016383 | affinity purification | MYR1 | TGME49_254470 | 363 | 1 | 100.0 | 10.74 | 0.55 |
| PXD043808 | affinity purification | AP2XI-2 | TGME49_310900 | 710 | 1 | 100.0 | 9.93 | 0.34 |
| PXD043924 | affinity purification | Cdt1 (TgiRD1) | TGME49_247040 | 1774 | 1 | 100.0 | 8.67 | 0.14 |
| PXD038848 | affinity purification | TgBDP1 | TGME49_263580 | 244 | 1 | 100.0 | 10.68 | 0.39 |
| PXD083501 | affinity purification | MyoL | TGME49_291020 | 1829 | 1 | 100.0 | 9.52 | 0.09 |
| PXD067653 | affinity purification | GRA35 | TGME49_226380 | 33 | 1 | 100.0 | 12.42 | 0.67 |
| PXD051551 | affinity purification | TgHCF101 | TGME49_318590 | 2171 | 1 | 100.0 | 9.56 | 0.09 |
| PXD051495 | affinity purification | PP2A-A2 | TGME49_229350 | 704 | 1 | 100.0 | 13.46 | 0.07 |
| PXD071065 | affinity purification | TgCoqFMO | TGME49_232540 | 2084 | 2 | 100.0 | 9.21 | 0.11 |
| PXD051495 | affinity purification | PP2A-C2 | TGME49_215170 | 704 | 2 | 99.9 | 10.75 | 0.21 |
| PXD032362 | affinity purification | BCC7 | TGME49_311230 | 491 | 2 | 99.8 | 7.73 | 0.60 |
| PXD032360 | affinity purification | GRA28 | TGME49_231960 | 428 | 3 | 99.5 | 7.63 | 0.52 |
| PXD065585 | affinity purification | MORC | TGME49_305340 | 735 | 6 | 99.3 | 9.23 | 0.69 |
| PXD059579 | proximity labelling | GSK3 | TGME49_265330 | 1748 | 13 | 99.3 | 6.40 | 0.44 |
| PXD041352 | affinity purification | GRA57 | TGME49_217680 | 502 | 5 | 99.2 | 7.84 | 0.70 |
| PXD043994 | affinity purification | Cyclin-4 (TgCyc4) | TGME49_249880 | 1085 | 10 | 99.2 | 7.23 | 0.07 |
| PXD063615 | affinity purification | PDIA3 | TGME49_211680 | 698 | 7 | 99.1 | 6.77 | 0.69 |
| PXD016846 | affinity purification | MORC | TGME49_305340 | 808 | 8 | 99.1 | 10.46 | 0.69 |
| PXD039919 | affinity purification | SPARK | TGME49_268210 | 749 | 10 | 98.8 | 6.91 | 0.16 |
| PXD065585 | affinity purification | AP2X-1 | TGME49_227900 | 735 | 10 | 98.8 | 8.26 | 0.14 |
| PXD039896 | affinity purification | C3 | TGME49_286470 | 1167 | 19 | 98.5 | 7.67 | 0.31 |
| PXD065585 | affinity purification | HDAC3 | TGME49_227290 | 735 | 13 | 98.4 | 7.88 | 0.50 |
| PXD066820 | affinity purification | TGME49_320300 | TGME49_320300 | 1688 | 37 | 97.9 | 6.96 | 0.26 |
| PXD044592 | affinity purification | Crk4 | TGME49_256070 | 1016 | 31 | 97.0 | 4.97 | 0.07 |
| PXD037318 | affinity purification | TgH1-like | TGME49_315570 | 65 | 4 | 95.4 | 7.33 | 0.24 |
| PXD044091 | affinity purification | CLAMP | TGME49_212270 | 731 | 66 | 91.1 | 4.45 | 0.19 |
| PXD039979 | affinity purification | SPARKEL | TGME49_291010 | 1141 | 213 | 81.4 | 3.14 | 0.09 |
| PXD039986 | proximity labelling | SPARKEL | TGME49_291010 | 1398 | 267 | 81.0 | 3.46 | 0.09 |
| PXD039983 | proximity labelling | SPARK | TGME49_268210 | 1092 | 220 | 79.9 | 3.85 | 0.16 |
| PXD063619 | proximity labelling | TgPDIA3 | TGME49_211680 | 825 | 199 | 76.0 | 3.09 | 0.69 |
| PXD032319 | affinity purification | ROP1 | TGME49_309590 | 38 | 14 | 65.8 | 4.46 | 0.69 |
| PXD044080 | affinity purification | TGME49_306920 | TGME49_306920 | 1229 | 480 | 61.0 | 3.46 | 0.09 |
| PXD039985 | proximity labelling | SPARKEL | TGME49_291010 | 974 | 397 | 59.3 | 0.89 | 0.09 |
| PXD039922 | affinity purification | SPARK | TGME49_268210 | 2261 | 1110 | 51.0 | 3.25 | 0.16 |
| PXD021516 | affinity purification | RON4 | TGME49_229010 | 1478 | 801 | 45.9 | 2.03 | 0.83 |
| PXD044081 | affinity purification | CDPK1 | TGME49_301440 | 3592 | 2021 | 43.8 | 2.23 | 0.62 |

Bait in the top 5 % of its own hit list: **26 of 37** overall and **25 of 32** of the affinity purifications; median bait rank 10 of a median 808 proteins, against median rank 26 if the list is ranked by raw abundance instead. `bait_recovery.tsv`.

The split by class matters. The score is built for affinity purification, where the tagged protein is pulled out intact and should dominate its own eluate. In proximity labelling the bait is a biotin-ligase fusion whose job is to label its *neighbours*, so there is no reason for it to top the list, and in a degron experiment (`PXD039922`, SPARK-mNG-mAID) the bait is deliberately degraded, so a low rank there is the correct answer rather than a failure. Those cases are kept in the table rather than dropped, because dropping them would be choosing the evidence.

No experiment here baits the glideosome, AMA1-RON, MIC2-M2AP or MORC/HDAC3, so none of them has a pull-down in which it *should* be specific. The best specificity each member reaches in any of the 60 experiments is in `complex_best_specificity.tsv`; the highest per complex:

| complex | protein | best in | rank | percentile | specificity |
|---|---|---|---|---|---|
| AMA1-RON | AMA1 | PXD011692 | 25 | 60.0 | 4.58 |
| AMA1-RON | RON2 | PXD043994 | 145 | 86.7 | 3.73 |
| MIC2-M2AP | MIC2 | PXD032360 | 143 | 66.8 | 3.85 |
| MIC2-M2AP | M2AP | PXD043924 | 175 | 90.2 | 3.28 |
| MORC/HDAC3 | MORC | PXD016846 | 8 | 99.1 | 10.46 |
| MORC/HDAC3 | HDAC3 | PXD065585 | 13 | 98.4 | 7.88 |
| MyoA glideosome | GAP45 | PXD041352 | 18 | 96.6 | 5.99 |
| MyoA glideosome | GAP50 | PXD012103 | 12 | 96.7 | 5.43 |

Per complex:

| complex | protein | gene | experiments detected | frequency |
|---|---|---|---|---|
| AMA1-RON | AMA1 | TGME49_255260 | 45 / 60 | 0.750 |
| AMA1-RON | RON2 | TGME49_300100 | 49 / 60 | 0.817 |
| AMA1-RON | RON4 | TGME49_229010 | 50 / 60 | 0.833 |
| MIC2-M2AP | M2AP | TGME49_214940 | 36 / 60 | 0.600 |
| MIC2-M2AP | MIC2 | TGME49_201780 | 35 / 60 | 0.583 |
| MORC/HDAC3 | HDAC3 | TGME49_227290 | 30 / 60 | 0.500 |
| MORC/HDAC3 | MORC | TGME49_305340 | 41 / 60 | 0.683 |
| MyoA glideosome | ELC1 | TGME49_269442 | 25 / 60 | 0.417 |
| MyoA glideosome | GAP40 | TGME49_249850 | 36 / 60 | 0.600 |
| MyoA glideosome | GAP45 | TGME49_223940 | 39 / 60 | 0.650 |
| MyoA glideosome | GAP50 | TGME49_219320 | 48 / 60 | 0.800 |
| MyoA glideosome | MLC1 | TGME49_257680 | 44 / 60 | 0.733 |
| MyoA glideosome | MyoA | TGME49_235470 | 51 / 60 | 0.850 |

Per frequent-flyer family (curated sets; the cytosolic ribosome is the hyperLOPIT 40S/60S marker set, so the apicoplast and mitochondrial ribosomes are not mixed in):

| family | proteins | median frequency | max frequency |
|---|---|---|---|
| GAPDH | 2 | 0.808 | 0.833 |
| HSP70/HSP90 | 3 | 0.883 | 0.933 |
| actin | 1 | 0.950 | 0.950 |
| cytosolic ribosome (hyperLOPIT marker) | 68 | 0.800 | 0.933 |
| enolase | 2 | 0.275 | 0.550 |
| tubulin | 6 | 0.550 | 0.950 |

The same families as a keyword would define them, for comparison only — these sweep in the tubulin-tyrosine ligases and the organellar ribosomes, which are neither abundant nor sticky:

| family (keyword) | proteins | median frequency | max frequency |
|---|---|---|---|
| actin (any) | 12 | 0.300 | 0.950 |
| heat shock (any) | 23 | 0.250 | 0.933 |
| ribosomal protein (any) | 158 | 0.208 | 0.933 |
| tubulin (any) | 23 | 0.033 | 0.950 |

## Is it usable yet

* **The frequency column as a standalone filter: not yet, and this is the headline.** The frequent flyers are above the complex members (median 0.80 against 0.68, AUC 0.689, p = 0.0144) but the overlap is almost total: at the 50 % cutoff, 12 of the 13 known complex members are deleted. Do not filter a hit list on frequency.
* **The specificity score: usable for ranking, with the caveats below.** It puts the tagged bait in the top 5 % of its own hit list in 25 of the 32 affinity purifications with a confident bait, and on the MYR1 worked example it puts MYR1 first, MYR4 third and five of six known complex members in the top fifth. The work is done by the abundance-fold term, not the frequency term.
* **The limiting resource is independent baits (30 distinct, over 60 experiments), not proteins or spectra.** Nearly every experiment is a pull-down from whole parasite lysate, so the abundant invasion and motility machinery is in every list. A frequency becomes informative only when there are enough baits from enough different compartments that "present in most pull-downs" means stickiness rather than abundance. For comparison, the human background frequency (Mellacheruvu et al. 2013) was built from more than 400 experiments.
* **Practical advice for now:** rank a hit list by `specificity`, or by its `log2_fold` component alone, and read `frequency` as context — "this protein has been seen in 51 of 60 other pull-downs" — rather than as a filter. The `median_abundance` column is the useful part of the table today: it says how much of a protein a pull-down normally gives, which is what makes an enrichment judgeable.

## What will change as the queue fills

* **The frequency resolution.** With 60 experiments the smallest step a frequency can take is 1/60 = 0.017. Every frequency in this table is that coarse.
* **Bait independence.** 30 distinct confident baits is the binding constraint, not the protein count. Several of the finished projects come from the same laboratories and the same complexes (SPARK/SPARKEL, CDPK1, the mitochondrial complexes), so the experiments are not as independent as their number suggests. A protein that happens to sit in two of the few baited complexes looks like a frequent flyer for a reason that has nothing to do with stickiness.
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
