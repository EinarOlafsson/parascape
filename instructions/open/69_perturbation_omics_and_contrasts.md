# 69 · Perturbation omics: every gene-knockout, drug and strain study, compared condition by condition

**Status: open (requested 2026-10-17).**

## The request (Einar, 2026-10-17)
* "Users should be able to also compare expression profiles between toxoplasma strains and strains under
  different conditions like thyroxin knock out" -- clarified: a **thioredoxin** knockout; **Type I/II/III**
  strains; the output is **per-gene differential expression** (condition A vs B: log2 fold change, a
  replicate-based test, FDR, a ranked table and a volcano plot).
* "I'd like as many gene ko and drug studies as possible with RNA seq."
* "And any other omics readouts" -- proteomics, metabolomics, ribosome profiling, chromatin.

## Plan
1. **Catalogue** (notebook, APIs only -- no accession typed from memory): every *T. gondii* series in GEO
   (NCBI E-utilities), PRIDE (proteomics), MetaboLights and Metabolomics Workbench (metabolomics).
   Classify each: gene perturbation (knockout, knockdown, auxin/AID, CRISPR), drug or compound, strain
   comparison, stage or time course, host-side only, other. Record readout type, sample count, replicate
   structure and whether replicate-level quantities (counts, intensities) are downloadable.
2. **Contrast engine** (`parascape/expression_contrast.py`): replicate-level matrices with sample
   annotations (strain, condition, perturbation, drug, dose, time, replicate); limma-trend (log-CPM,
   empirical-Bayes moderated variance; Law et al. 2014) in numpy/scipy, since neither R nor DESeq2 is
   available. Validated on simulations (null p-values uniform, FDR controlled) and against a published
   dataset's reported results before any is shipped.
3. **Load studies** in priority order: thioredoxin knockout; Type I/II/III tachyzoites; then gene
   perturbations and drugs by replicate quality; then other omics with readout-appropriate tests.
   Gene identifiers through the identity layer (GT1/ME49/VEG map by numeric suffix).
4. **UI**: pick a study and two groups -> ranked differential table, volcano, and the gene list handed to
   the strategies, maps, star map and track record.

## Rules that apply
* Never type a PMID, DOI or accession from memory; record every query.
* VEuPathDB is login-only (401); do not work around it. GEO/SRA/PRIDE/MetaboLights are fine.
* Leakage: a contrast built from a dataset is that dataset's family; register its leakage family.
* RAM: the ledger at `/media/carruthers/mnt3/claude/ram_leases.md` (a plaque-assay priority asks that jobs
  >= 8 GB wait).

## Decisions (Einar, 2026-10-17)
* **Quantify every documented convention and let the user choose** -- quantifier (salmon with genome
  decoys and bias correction; kallisto), reference (ME49 for all; GT1 / VEG strain-matched references for
  strain studies), and at contrast time the normalisation (CPM, TMM, median-of-ratios, upper quartile,
  TPM) and the authors' own processed values beside the reprocessing.
* **Datasets of one type must be comparable to each other**: offer normalisation within a study (each
  sample against its own study's control: log2 fold change to the matched wild type / vehicle) as well as
  across studies.
* **Storage: the NAS**, `/nas_mnt/carruthers/Einar/parascape_omics` (references, per-run
  quantifications, run metadata, logs, the plan). The NAS is shared and was 95% full (4.1 TB free); the 58
  RNA-seq candidate studies hold ~2.5 TB of compressed reads (791 runs). Decision: **stream** -- download a
  run from ENA, quantify, keep the results, delete the reads (re-downloadable). Scratch:
  `/media/carruthers/mnt3/claude/parascape_omics_scratch` (mnt3 was 96% full, 149 GB free; the script
  waits for space).
* **A separate conda env `parascape-quant`** (salmon 2.8.0, kallisto 0.52.0, sra-tools 3.4.1), created
  with Einar's approval; no other env touched.

## Found so far (`notebooks/omics_catalogue_2026_10_17.ipynb`)
* 606 studies: GEO 371, PRIDE 209, SRA-only 19, Metabolomics Workbench 7.
* The thioredoxin knockout: **PRJNA1008424, "RNA-Seq of Toxoplasma gondii Trx1 mutant"** (SRA only, 8
  paired runs, 32 GB). It has an empty organism field, so the organism-tagged BioProject search misses it;
  the catalogue needs an all-fields search too. Related: Trx4 (PMID 38576040, no RNA-seq found yet);
  TgPRO, a redox regulator (RNA-seq GSE329845, TMT proteome PXD077563).
* Type I/II/III: **GSE60305**, RH, PLK, CTG as tachyzoites and bradyzoites (18 runs, single-end).

## Progress (2026-10-18)
* Quantification restarted LOW PROFILE at Einar's request (user unit `parascape-quant`: one run at a time,
  3 threads, nice 19, idle IO, 12 GB cap -- a 67M-pair run against ME49_host was OOM-killed at 10 GB on
  2026-10-17). 17 of 1,011 runs done. The shared scratch disk (mnt3) filled to 0 GB overnight; the job
  waits for space rather than competing. At this pace the plan takes weeks.
* First studies built into matrices: GSE103113 (GRA18, one sample per condition: profiles only, no test;
  also mis-flagged as a strain study), GSE116539 (TgCRT knockout, 2 vs 2, purified parasites, host share
  0.05%), GSE125120 (Myr1 KO in HFF; "bystander" samples are uninfected cells with ~0.2M parasite reads --
  only the infected samples support a parasite comparison).
* **Validation against a paper** (Thornton et al. 2019, PLoS Pathog, PMID 31170269, found via Europe PMC
  from GSE116539): the knockout lowers endolysosomal proteases. Reproduced: SUB1 -0.31 log2 (FDR 0.008),
  CPC2 -0.26, ASP1 -0.79, on TMM, median-of-ratios and CPM alike.
* **Trap:** TgCRT itself reads ~10x HIGHER in its own knockout (salmon and kallisto agree; the labels are
  right, since the protease result matches). Likely the construct leaves CRT sequence transcribed. A
  gene's own expression in its knockout is not evidence the knockout worked; the UI should say so when the
  knockout target is known.
* Fixed: replicate stripping turned "gra18" into "gra"; a number is now a replicate only when set off by a
  separator or named as one.
