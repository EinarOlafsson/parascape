# ParaScape progress

Three sections, as Einar asked on 2026-11-18: **the paper**, **the program**, **the data**.

Every number here was measured on the date given, by the command named beside it, and not carried
over from an older note. Instruction 69 recorded "17 of 1,011" for a month after it stopped being
true, which is why the command is written down next to the count.

Last measured: **2026-11-18**.

---

## 1. The paper

Instruction [101](open/101_parascape_paper.md). Twelve figures, each titled by the conclusion it
reaches rather than the data it uses. Figure 2 gates the rest.

| # | Figure — the claim it ends on | State | % |
|---|---|---|---|
| 1 | What ParaScape is: imports, modules, one gene space | parts exist, never assembled | 20 |
| 2 | **Which studies can be compared at all** — control clustering, cell-cycle normalisation, the ladder | ladder built and validated; the clustering is not written | 25 |
| 3 | RNA and protein co-expression, clustered → what the modules are | components built | 40 |
| 4 | RNA response to treatment → which conditions group, and why | 201 contrasts on one scale | 35 |
| 5 | Protein response to treatment; bradyzoite drug reach → subtype C is reached by 4 of 201 contrasts | Reach module built, swept, verified | 60 |
| 6 | Which measurement explains which other → wave-built organelles are placeable from transcription alone | localisation panel built; 4 panels to go | 25 |
| 7 | The network: traits on nodes, derived relationships on edges | edge layers exist, the graph does not | 10 |
| 8 | Strains → what differs, and whether it is the effectors | data loaded, nothing drawn | 5 |
| 9 | Species, with measured transfer per layer → when one species licenses a claim about another | one layer measured (Cp 56.7 % against 21.9 %) | 15 |
| 10 | The host: response, vacuole, co-folding, dual perturb-seq | atlas built; the per-host-gene import is missing | 20 |
| 11 | Filling the map, scored on held-out facts | scoreboard exists, the figure does not | 15 |
| 12 | Measured absence → which datasets the field most needs | limits machinery exists, never assembled | 5 |
| — | **The manuscript**, each figure in a two-cell table with the legend below | not started | 0 |

**Named next jobs:** the figure 2 control-comparability analysis; importing the GSE229505 host
differential-expression table per host gene.

**Decisions taken (Einar, 2026-11-18):** twelve figures, no merging · the network is not
co-expression drawn as a graph, so the edge layer becomes an object rather than a view · the
co-folding panel reports scale and states that it is uncontrolled · a heatmap was rejected for the
localisation panel, because colour cannot separate 0.94 from 0.93.

---

## 2. The program

24 modules, 8 built on their own shell, 32 named views, **7,709 tests passing** (`pytest tests/ -q
-m "not slow"`, 2026-11-18, 61 minutes, 0 failed, 60 skipped).

| Item | State | % |
|---|---|---|
| Every analysis in the paper runnable inside ParaScape | 1 of 12 figures has a module (Reach) | 8 |
| Each module saves its result with what it was made from | Reach does; nothing else | 10 |
| Reach: the cluster-against-treatment module | **done** — `300cac69`, 12 tests in a real window | 100 |
| Interface on the spaCR pattern (instruction 98) | done | 100 |
| Organism spaces (instruction [53](open/53_organism_spaces.md)) | Hs, Mm, Cp built; packs, cross-space bridges, per-space calibration open | 60 |
| Hugging Face hosting (instruction [99](open/99_data_on_hugging_face.md)) | core data and 1,205 components published; the instruction-53 packs and the release-process upload step remain | 75 |
| Edge layer as a first-class object (figure 7) | **blocked on figure 7** | 0 |
| Modules reorganised to match the figure structure | **blocked on the figures** | 0 |
| Modules no figure needs retired or merged | **blocked on the figures** | 0 |

The last three are deliberately not started. A module built for a figure nobody has drawn is a
module built for a guess.

---

## 3. The data

### Slots filled
`python scripts/generate_slot_table.py`, 2026-11-18: **362 slots — A 127, B 58, C 58, empty 119**,
so 243 of 362 carry something (67.1 %).

| Kind | Filled | % |
|---|---|---|
| gene | 194 / 278 | 69.8 |
| pair | 19 / 28 | 67.9 |
| host gene | 27 / 50 | 54.0 |
| metabolite | 3 / 6 | 50.0 |

177 datasets: Tg 101, Pf 55, Hs 15, Mm 4, Ag 2.

### Re-quantification (instructions [69](open/69_perturbation_omics_and_contrasts.md), [81](open/81_requantification_topup.md))
`find $PARASCAPE_OMICS/rnaseq -maxdepth 3 -name done.json | wc -l`, 2026-11-18: **1,538 runs
finished**. The service reports **671 new runs still to quantify** (unit `parascape-requant`,
active, hourly log). So 1,538 of about 2,209, **69.6 %** — and the denominator is still growing,
because downloads keep arriving. The 2026-11-17 note of "1,446 of 1,470" was true for the runs that
existed then; it is not a completion figure.

### Downloads (unit `globus-downloads-monitor`, active)
One Funk transfer running at 4.77 MB/s over the last 24 hours. Queued behind it: S-BIAD3096 (four
plate directories), then S-BIAD877.

> **Disk warning, measured not predicted.** Free space on the download target fell 3,477 → 2,690 GB
> in the 24 hours to 2026-11-18, about 790 GB a day. At that rate it runs out in roughly **three to
> four days**. This needs a decision before it becomes an incident.

### Open data instructions
| # | What is left |
|---|---|
| [41](open/41_fill_the_slots.md) | 119 empty slots; two Cryptosporidium slots filled 2026-11-18 (`2834be6c`) |
| [65](open/65_dataset_selection_audit.md) | A grades overstate coverage — bone marrow 9.8 %, dermis 22.5 %, salivary gland 19.5 % |
| [69](open/69_perturbation_omics_and_contrasts.md) | proteomics, metabolomics, cross-study profiles, the expression pack |
| [80](open/80_pulldown_background_baits.md) | 22 baits; the gate AUC fell 0.710 → 0.684 as projects doubled, so waiting does not reach it |
| [81](open/81_requantification_topup.md) | queued behind the new downloads; 40 zero-parasite runs first |
| [99](open/99_data_on_hugging_face.md) | the instruction-53 packs, and a token to publish them |

---

## Waiting on Einar
* Which figure to build first — figure 2 (control comparability) is the recommendation.
* The download disk, above.
* The Plasmodium hold (`deposits.HELD_OUT_OF_NODES`) stays until he says otherwise.
* A Hugging Face token, to publish the instruction-53 packs.
* Whether the 40 zero-parasite runs jump the re-quantification queue.
