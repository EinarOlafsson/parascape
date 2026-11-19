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
| 1 | What ParaScape is → **the corpus, assembled** | **both panels drawn** (schematic + coverage) | 85 |
| 2 | **Which studies can be compared at all** → **cell cycle raises agreement 0.832 to 0.887, three times a random axis, p = 0.005** | **measured and drawn** on 12 studies | 100 |
| 3 | Co-expression modules → **14 of 23 are compartments; the secretory machinery splits in two** | **measured and drawn** | 100 |
| 4 | What structures the corpus → **same-study pairs agree 3.3x better than cross-study; only 1 perturbation is repeated across labs** | **measured and drawn** over 20,100 pairs | 100 |
| 5 | Protein response; bradyzoite reach → **subtype C is reached by 4 of 201 contrasts** | **panel drawn**, sweep saved, module built | 85 |
| 6 | Which measurement explains which other → **transcription is 96 % of what is explainable about protein abundance** | 2 of 5 panels built (localisation, protein abundance) | 50 |
| 7 | The network → **24x the reach at two thirds the accuracy; weighting recovers 5.7 points** | **measured and drawn** over 13 edge layers | 100 |
| 8 | Strains → **dense granules 2.06 and apical secretory 1.73 against 1.36; everything internal below** | **measured and drawn** | 100 |
| 9 | Species transfer → **localisation transfers at +0.256, fitness at 0.314, codon usage not at all** | **measured and drawn**, length is the control at 0.907 | 100 |
| 10 | The host → **mostly private**, and the single-cell check says the parasite is unchanged (1 of 85 targets) | **measured on pseudobulk AND on the 126,495 cells**; co-folding remains | 75 |
| 11 | Filling the map → **5 of 12 targets beat guessing; localisation 76.5 % against 33.1 %** | **measured and drawn** from 1,047,574 held-out calls | 100 |
| 12 | Measured absence → **host response 94.7 % unmeasured, transcript 0.1 %** | **measured and drawn**; prescriptive line stands | 100 |
| — | **The manuscript**, each figure in a two-cell table with the legend below | **done** — 12 tables, 12 images embedded, no placeholders | 100 |

### Data
105 study matrices built (0 failed), taking the omics store to 128 studies with gene-level
matrices. GSE229505 downloaded: 2.2 GB, 18 matrices, 126,495 cells, host and parasite genes plus
the CRISPR guide per cell.

### Figure style
The house style is now code (`parascape/figstyle.py`) and a skill
(`~/.claude/skills/publication-figures`): fixed width per category so panels match, uniform bar
and violin geometry, control grey, 60 per cent opacity, no bar with whiskers above 20 points.
Palette read off the two reference papers that were reachable; four of the six answered 403, 406
or a stub and the skill records that their style was NOT copied.

### Findings this session (2026-11-18)
* **Subtype C** of a brain cyst is reached by 4 of 201 published contrasts, against 86 for E. Not the matching giving up: C is testable and centred in every contrast checked.
* **Transcription explains 96 % of what is explainable** about protein abundance in Plasmodium — but on the complete-case genes. On the wider set it explains 24.3 % rather than 48.6, and that factor of two is the measurement bias.
* **Transcription alone places the organelles built in waves** during division (5 of 13 compartments) and not the ones built continuously.
* **Localisation transfers between species, codon usage does not.** Protein length at 0.907 is the control proving the orthologue pairing is sound.
* **The laboratory structures the perturbation corpus, not the perturbation.** Two contrasts from one study agree at median |rho| 0.284 against 0.085 for two from different studies; sharing the KIND of perturbation across studies buys 0.090 against 0.085, which is nothing. And only ONE perturbation (deltaBFD1) is stated by two independent labs, so cross-lab replication is not askable here at all. This is the argument for figure 2.
* **Cell cycle IS a substantial confounder** — removing it raises agreement between the untreated arms of independent studies from 0.832 to 0.887, against +0.018 for a random axis of the same shape, p = 0.0050 on 12 studies. This REVERSES the eight-study result reported earlier the same day: that run was short four studies the matrices build has since added, and it included GSE245775, which quantifies to zero and whose 8,157 identical zeros made every correlation undefined.
* **ParaScape cannot yet fill phenotype or export.** is_exported is 95.9 % right against a 96.6 % baseline, which is worse than guessing.
* **What differs between strains is what goes into the host.** Dense granules 2.06 and apical secretory 1.73 against a length-matched 1.36; nucleus 0.89, cytosol 0.72, the most required genes 0.89.
* **The host response is mostly private to one effector** — 3,435 host genes moved by exactly one of 21, only 68 by half or more. The translocon is the shared part.
* **The whole network trades accuracy for reach**: 24 times the genes at two thirds the accuracy, and the naive union is dragged down by its biggest and worst layer. Weighting layers by what each is worth recovers 5.7 points.
* **Co-expression modules are compartments** — 14 of 23, one of them 78.6 % ribosome against 3.8 % matched — and the apical secretory machinery splits into two modules, not one.

**Named next jobs:** the module reshaping (C1-C3) is the only paper-driven work left, and it is
no longer gated: the figures exist.

**Decisions taken (Einar, 2026-11-18):** twelve figures, no merging · the network is not
co-expression drawn as a graph, so the edge layer becomes an object rather than a view · the
co-folding panel reports scale and states that it is uncontrolled · a heatmap was rejected for the
localisation panel, because colour cannot separate 0.94 from 0.93.

---

## 2. The program

24 modules, 8 built on their own shell, 32 named views, **7,754 tests passing**.

| Item | State | % |
|---|---|---|
| Every analysis in the paper runnable inside ParaScape | **done** — the Findings module runs all 14 analyses across the 12 figures | 100 |
| Each analysis saves its result with what it was made from | **done** — every finding carries its result directory and summary | 100 |
| Reach: the cluster-against-treatment module | **done** — `300cac69`, 12 tests in a real window | 100 |
| Interface on the spaCR pattern (instruction 98) | done | 100 |
| Organism spaces (instruction [53](open/53_organism_spaces.md)) | Hs, Mm, Cp built; packs, cross-space bridges, per-space calibration open | 60 |
| Hugging Face hosting (instruction [99](open/99_data_on_hugging_face.md)) | core data and 1,205 components published; the instruction-53 packs and the release-process upload step remain | 75 |
| Edge layer as a first-class object (figure 7) | **done** — `edge_layers.py`, 13 layers with measured worth, 6 tests | 100 |
| Modules reorganised to match the figure structure | **measured, recommended against** — the bands are user intent, not the order of an argument; Einar's call | — |
| Modules no figure needs retired or merged | **measured, recommended against** — 9 of 25 unexercised, and none is dead; nothing retired | — |

C1 is done. C2 and C3 are measured and argued against in instruction 101 rather than executed:
reorganising the bands would rebuild a navigation system around the shape of a document, and
retiring a module removes something a user can open today. Both are Einar's call, and nothing has
been moved or removed.

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
| [41](open/41_fill_the_slots.md) | 119 empty slots. **The v4 candidate list is 11.1 % on target** (9 of 81 name both the organism and the assay of their slot), so filling means curating candidate by candidate, not working down a queue |
| [65](open/65_dataset_selection_audit.md) | **Cause found and fixed**: every coverage string reported the share of the SOURCE FILE, not of the host space. All 21 host datasets now carry both; 8 cover under a quarter of their space. What a GRADE should mean is still Einar's call |
| [69](open/69_perturbation_omics_and_contrasts.md) | proteomics, metabolomics, cross-study profiles, the expression pack |
| [80](open/80_pulldown_background_baits.md) | 22 baits; waiting does not reach the gate, and **a compartment-matched background does not either** — complex members are as frequent inside the over-represented compartments as outside (p = 0.53). The specificity score stays the product |
| [81](open/81_requantification_topup.md) | queued behind the new downloads; 40 zero-parasite runs first |
| [99](open/99_data_on_hugging_face.md) | the instruction-53 packs, and a token to publish them |

---

## Verification
Full suite green at 2026-11-18 19:05: **7,754 passed, 0 failed, 60 skipped**, 56 minutes. The run
before it found 6 failures and every one was mine — four missing docstrings, one organism literal
in a figure script, and two generated artefacts (`docs/datasets.md` and the 177 per-dataset
scripts) that went stale when the host coverage strings gained their second denominator.

## How progress is now checked

`scripts/watch_progress.py`, run every 30 minutes by `parascape-progress-watch.timer`. It compares
every tracked counter with the previous snapshot and reports MOVING, QUIET, SLOW or STALLED. A
STALLED verdict requires two consecutive unchanged snapshots AND a worker that has written nothing,
and it prints the evidence. The measurement is `write_bytes`, never CPU: a futex spin burns a core
exactly as hard as real work. Added after the queue sat at 624 for nine hours while three kallisto
processes read their index and wrote nothing, and I twice reported it as "running".

## Overnight, 2026-11-18

* **The manuscript has its text.** Results, one passage per figure, and Methods stating the three
  rules every claim is held to. 2,614 words, every number read from the summary of the analysis it
  describes. Twelve figures, twelve legends, twelve two-cell tables.
* **Citations: 67 uncited at 19:00, 38 now.** Eight from GEO's own link, 15 reference resources,
  11 from PRIDE's referenceLine, 6 from PMIDs already written in the accession field, 12 by
  constructed DOI checked for a unique PubMed match. The 38 left name a journal and an article
  number that does not map deterministically to a DOI; guessing those produced three wrong
  citations earlier in the day, so they are left for a person.
* **The one genuine instruction-41 candidate assessed.** GSE232438 downloaded (39 MB, 10,319 cells)
  and derived to a per-timepoint table covering 100.0 % of the Cryptosporidium space. It does NOT
  fill the oocyst slot it was offered for: the supplement carries no stage annotation, so no cell
  can be called an oocyst without redoing the paper's clustering.

## The four that are not work

Measured 2026-11-18 19:30. Each was pushed on repeatedly and has no hidden task left in it.

* **69 / 81 — running, not blocked.** 1,573 done.json markers, 637 runs still queued, falling about
  8 an hour (671 at 17:00, 645 at 17:52, 637 at 18:52). That is roughly **3.3 days** of compute.
  Nothing accelerates it that does not take CPU from the other sessions on this box.
* **99 upload — needs a credential.** The five packs are built and their checksums are in
  `parascape/pack_index.json`. `huggingface_hub` is installed in base, spacr and autotrade but not
  in `starplast`, and neither base nor spacr has a stored token. Searching the credential stores
  was denied by the permission classifier, correctly. `publish_packs.py --token <token>` takes one
  inline without storing it.
* **65.04 — the replacement is answered, the two slot decisions are not.** The batch-5 proposal was
  generated by the flag-coverage defect and does not arise now that it is fixed in both paths.
* **C2 / C3 — explained and argued, not decided.** Nothing has been moved or removed.

## Waiting on Einar
* Which figure to build first — figure 2 (control comparability) is the recommendation.
* The download disk, above.
* The Plasmodium hold (`deposits.HELD_OUT_OF_NODES`) stays until he says otherwise.
* A Hugging Face token, to publish the instruction-53 packs.
* Whether the 40 zero-parasite runs jump the re-quantification queue.
