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

## Download takeover (2026-10-20, America/Detroit)

User requested taking over the Claude `(work) Starplast` session and downloading as much as possible
without affecting the active spaCR tutorial. Scope here is transport and resource management of the
existing plan. The broader inference/action list remains paused.

- At takeover: 92/1,011 runs complete, 919 pending. Existing processing was already running at five
  threads, nice 19, idle IO and a 14 GiB memory cap. Previously claimed CPUQuota and IO bandwidth
  cgroup limits are **not enforced** on this user session; do not describe them as effective caps.
- `scripts/omics_download_cache.py` maintains a resumable cache on the root filesystem at
  `/home/carruthers/.cache/parascape-downloads`, capped at 40 GiB including partials, retaining at least
  90 GiB free. One download at a time, maximum 4 MiB/s, one CPU affinity (15), nice 19, idle IO,
  384 MiB memory cap. It yields when CPU pressure reaches 5%, IO pressure 8%, memory pressure 0.5%,
  or available RAM falls below 32 GiB. These are conservative controls, not a guarantee of zero
  effect on other jobs. No tutorial process, file, setting or GPU allocation is changed.
- Existing sequential processing consumes checksum-verified cache entries and still runs its SRA
  integrity validation; cached raw archives are deleted after successful expansion. Direct-download
  fallback is also capped at 4 MiB/s. No raw reads are retained on the shared NAS.
- User units: `parascape-downloads.service` (cache), `parascape-quant-cached.service` (consumer).
  The old transient `parascape-quant.service` is replaced only after its current run finishes.
  `parascape-download-handover.service` waits for that boundary; its durable record is
  `~/.cache/parascape-downloads/handover.json`. Do not start a second consumer while it is waiting.
  The consumer stays at five threads and has effective process CPU affinity 10–14, plus its original
  14 GiB memory cap. User services use a fixed code copy in
  `~/.local/share/parascape-download-worker/`; later source changes require an intentional redeploy.
- Progress: cache `status.json`, `journalctl --user -u parascape-downloads`, and the existing
  NAS `logs/omics_quant.log`. To stop these new jobs use `systemctl --user stop
  parascape-download-handover parascape-downloads parascape-quant-cached`; stopping a service does not
  disable its next-boot start. All cached partials are resumable.
- Transport validation: synthetic local HTTP fixtures cover HTTP range resume, servers ignoring
  Range, interrupted transfers, disk reserve/capacity, checksum corruption, producer/consumer
  locking, throttling/pressure waits and consumption without duplicate downloads. Run them without
  Qt or GPU initialization: `python -m unittest discover -s tests -p test_omics_download_cache.py -v`.

Instruction 69 remains open while the queue and its other requested work are unfinished.

## Faster processing on request (2026-10-21)

At 09:45 EDT: 114/1,011 runs complete, 54 archives (81.53 GiB) downloaded since takeover,
32 archives (37.94 GiB) waiting in the cache. The handover completed at 22:03 EDT on October 20;
no download/processing failure had been logged since takeover. Processing, rather than the
network, was limiting throughput. One current run had spent over four hours in its final tool;
its input position was observed advancing, so it was preserved.

The user then authorized faster work. The consumer memory cap increases from 14 to 20 GiB
(applied live), and its configured thread count from five to ten at the next completed-run
boundary. `parascape-speedup.service` waits for valid completion metadata and cleared scratch;
`~/.cache/parascape-downloads/speedup.json` records whether the switch has actually occurred.
A configured environment value alone does not prove the running process has adopted it.

`scripts/omics_resource_guard.py` controls only the consumer's existing process/thread CPU
affinities. `parascape-resource-guard.service` checks every five seconds, returns to five
CPUs when a tutorial capture/validation worker or resource pressure is detected, and allows
ten CPUs after 30 seconds of sustained quiet. Its stop handler restores the conservative
five-CPU affinity. The consumer also starts with that conservative affinity so a missing
guard does not enable extra CPU capacity. This uses enforced scheduler affinity, not the
unavailable user CPU cgroup controller. Ten requested threads may therefore run on five CPUs
while tutorial work needs priority. Nice 19, idle IO, the download rate, cache bound and free
disk reserve remain unchanged; no GPU or other session's process settings are changed.

State is in `~/.cache/parascape-downloads/resource-status.json`. Five synthetic guard tests
check tutorial detection, immediate yielding, delayed boosting, all-thread coverage, exclusion
of unrelated processes and exiting-worker races: `python -m unittest discover -s tests
-p test_omics_resource_guard.py -v`.

## Token-independent continuation (2026-10-21)

User requested that token exhaustion not stop the process and that it run to completion.
The ten-thread switch has taken effect (10:25 EDT); the processing service continues independently
of this chat. A persistent `parascape-queue.timer` and deterministic supervisor now recover missing
owned workers, repeat incomplete passes even when the process exits successfully, and retry with
bounded backoff without an attempt limit. User lingering is enabled; worker units are started by
the supervisor after logout/reboot rather than independently enabled at boot.

The current source plan is frozen with its hash under
`~/.cache/parascape-downloads/automation/`. It contains 1,011 study/run entries for 999 distinct
raw accessions; twelve runs belong to two studies. Both entries retain their requested outputs.
Completion requires matching provenance and SHA-256/gzip-integrity verification of all expected
compressed outputs for all entries. Invalid receipts are preserved under renamed filenames for
review and reprocessing. Active run output is left alone. Unchanged verified files reuse evidence.

The supervisor writes `completion.json` only after every entry verifies, then stops its owned
workers and timer. Temporary resource/source/storage waits stay pending. This automation uses no
model tokens or credentials and does not resume the separate paused inference/action list.
See [the unattended queue operations guide](../../docs/omics-queue.md) for status files, service
commands, verification and explicit operator pause/resume. Transport/guard/supervisor validation:
27 focused tests passed, and both user-systemd unit definitions passed verification.

### CPU guard correction (2026-10-21 afternoon)

A progress check still found 115 verified entries. Ten runnable tool threads restricted
to five CPUs were producing roughly 50% global CPU pressure despite about 75% of host CPU
capacity being idle. The old PSI-only rule therefore kept the worker permanently reduced
and also paused downloads. Resource checks now sample `/proc/stat` deltas: CPU PSI alone
does not prevent restoring the allowed ten CPUs when at least 25% host idle time remains.
Independent tutorial detection and IO/memory pressure still yield immediately. The worker's
current run is preserved; only the guard and resumable downloader need restarting to load
this transport/resource helper change. All 29 focused checks passed, including regression
cases for high PSI on an idle host and actual host CPU contention.

### Full-machine setting (2026-10-21 evening)

User stated tutorials have stopped and explicitly authorized all available CPU/RAM/GPU
capacity with no additional model-token use. Applied all-CPU scheduler affinity (0–31)
and a verified live 112 GiB memory cap. A deterministic boundary restart changes tool
threads to 32 and normal CPU/IO priority, preserving the current run. Future scratch is
on mnt1 (about 900 GiB free). Downloads have no application bandwidth cap; cache is
80 GiB with a 25 GiB disk reserve. `PARASCAPE_FULL_RESOURCES=1` bypasses CPU/IO pressure
pauses for downloads while retaining a 4 GiB available-RAM floor. The guard's explicit
`--full-resources` option keeps all CPUs available. Installed Salmon/kallisto help has
CPU thread controls and no GPU mode. No additional AI process or token usage is created.
The existing supervisor/retry/completion automation remains in control. Compatibility:
29 focused checks passed, with actual kernel affinity and memory cap verified.

### Eight-core limit (2026-10-22 evening)

User requested reducing CPU use to eight cores and a progress report. Applied CPU
affinity 0–7 to all existing queue threads and persisted it in worker/downloader/guard/
supervisor units. Both guard CPU sets are fixed to 0–7. Configured eight tool threads
with a completed-run boundary switch; current processing is preserved. Verified
128/1,011 entries, 883 remaining, 35 cached archives (77.4 GiB); automation active.

### Resume and processing-record check (2026-10-27)

User explicitly resumed downloads and requested automatic continuation and an explanation
of processing/provenance. Found Claude's separate original worker already running; the
earlier stopped-status check had missed that service. Supervisor now adopts its active
state and excludes its active output from verification, preventing duplicate writers.
A deterministic completed-run handover will restore the cache-consuming worker.

Root had 18 GiB free. Relocated only the queue cache (22 archives, 78.04 GiB) to mnt1,
verified all archive SHA-256s and metadata copies, then preserved the original cache path
as a symlink. Root free space increased to 96 GiB. PAUSE removed, persistent timer enabled;
downloads observed transferring. Supervisor verified 295/1,011 entries, 716 remaining.
Eight-CPU limit retained; no model/agent continuation started. Thirty focused tests passed.

Inspected actual completed outputs and reference records. Source/sample metadata and
versions are in done.json; tool invocations/diagnostics in Salmon cmd_info/meta_info and
kallisto run_info; source/reference checksums, commands and output hashes are retained in
their respective records. Added a guide in docs/omics-queue.md. Quantification and integrity
verification are automated; downstream differential-expression/comparison work is separate.


## 2026-10-29 offline snapshot release

User requested integration into the PyPI package and easy per-gene comparison.
Prepared version 0.0.2.0 with a separate descriptive RNA expression browser and CLI.
Snapshot scope: 296 verified entries / 23 studies from the 1,011-entry frozen queue;
competitive ME49_host results from Salmon and kallisto. Other references remain external.
No new prediction features, graph inputs or calibration claims are introduced.
Full queue and independent published-contrast validation remain unfinished.
Builder: `scripts/build_rna_snapshot.py`; guide: `docs/rna-expression.md`.
Release publication and verification evidence will be recorded after the checks.


### User scope update: GitHub data instead of bundled datasets

User requested every completed run, source-evidenced replicate/control/repeat
organization, Trx1 priority, then all application data hosted on GitHub with a
startup download prompt. The prepared wheel is ~1.6 MB; the versioned data asset
is ~221 MB / 394 checksum-pinned files. All available reference variants for
the frozen 296 runs are included; existing scientific tables retain their original
values and encodings. New data consumers use paths.data_dir rather than package
relative dataset paths. Icons and the small descriptor remain in the package.

655 release-gate checks and clean installed-wheel popup/cancel/install/offline
checks passed. Main-startup tests were updated for the new consent dialog; full
regressions and hosted publication verification are still in progress.
The live queue reached297 completed, and priority Trx1 run SRR25737844 is active;
the snapshot remains an explicit immutable296-entry cut.

## 2026-11-01: full catalogue queued, storage limits, named-expression tables
* The coverage ledger showed 109 catalogue RNA studies that were never queued. A read-only ENA review
  (`scripts/review_unplanned_rna_studies.py`) resolved all 109; `scripts/queue_unplanned_rna_studies.py`
  appended 675 Toxoplasma-named bulk runs (3.44 TB) and 620 human/mouse host-labelled runs (1.38 TB,
  `run_label` tagged) under the supervisor lock. Queue: 1011 -> 2306 associations.
* Kept out deliberately: per-cell single-cell GSE145800 (11,419 runs) and GSE145080 (5,499 runs), whose
  deposited AnnData matrices are already prepared; hosts without transcripts in ME49_host
  (cow, hamster, pig, cat); blank or non-Toxoplasma organisms; non-RNA libraries.
* The NAS raw tier is now hard-capped (~500 GB, Einar) after growing to 1.58 TB.
* All 105 GEO named-feature summary tables (13 gene / 7 transcript quantity families kept separate,
  Expression value unit unknown) were reconciled into 1,149 packs, 0 reviews.
* Normalization stance: comparable RNA values come from the uniform reprocessing (salmon + kallisto,
  ME49 and ME49+host, same settings for every run). Deposited tables keep their original units and are
  shown beside, never merged with, the reprocessed values.

## 2026-11-02: sample/context registry and Einar's comparison decisions
* `scripts/build_sample_registry.py` keeps every GEO sample annotation verbatim (371 series, 6,621 samples);
  `scripts/propose_context_mapping.py` maps tag/value pairs onto context fields through named rules plus a
  first-pass review by Claude (97% of samples have context; 149 pairs left for review).
* **Time (Einar: "check the associated publication"):** the clock's start was resolved for all 24 studies with
  unstated times from GEO text or the paper (`STUDY_TIME_REFERENCE`, evidence cited per study). Only times whose
  reference is infection are compared as infection time; others carry their own event (stress medium, drug,
  auxin, host stimulation, culture).
* **Scope (Einar: "keep everything"):** samples from other pathogens in shared studies (e.g. the ImmGenT atlas)
  stay in the registry as host context.
* **Independence (Einar: conservative):** samples sharing a mouse, donor or biological-replicate identifier count
  as one unit; without identifiers, a study's samples are not assumed independent beyond its stated design.

## 2026-11-03: what each run keeps, and the planned top-up
* Each run keeps transcript-level estimated counts, effective lengths and TPM from salmon and kallisto per
  reference, plus salmon QC, bias and command records. These are pre-normalization values: every
  normalization, differential-expression, batch-correction and meta-analysis method is applied later from
  them, so raw reads are not needed for any of those.
* From 2026-11-03 (Einar) salmon also writes equivalence classes with weights (`--dumpEqWeights`) and 30 Gibbs
  samples (`--numGibbsSamples 30`), so runs can be re-quantified against an updated annotation and carry
  per-transcript uncertainty without their reads (smoke-tested on salmon 2.8.0).
* **Queued for after all new data is downloaded (Einar):** re-download the runs finished before this change
  (370 verified at the time) and append the same extra outputs, without changing their verified results.


## 2026-11-03: uniform proteomics re-search (parascape-ms-uniform)
* Every Toxoplasma PRIDE project in `processing_audit/pride_ms_uniform_selection.json` (201) is streamed, re-searched
  with Sage 0.14.6 against one combined FASTA (ME49 + human/mouse longest isoforms + contaminants) and its RAW
  files deleted. Each run keeps centroided MS1+MS2 mzML.gz, instrument metadata and checksums; each project keeps
  Sage PSM, peptide/protein and LFQ tables with the exact config (`scripts/ms_uniform_search.py`).
* Formats: Thermo `.raw` (ThermoRawFileParser `-f 1 -g`), author mzML, zip/tar/rar/7z archives (7-Zip, same
  outside-the-folder check), Bruker timsTOF `.d` (zipped `.d` kept as the stored spectra, Sage searches it
  directly at 15 ppm precursor / 20 ppm fragment, MS2 identifications only, no MS1 LFQ).
* Recorded gaps (`unsupported.json`): SCIEX `.wiff` (4 projects; ProteoWizard needs Windows, Wine or Docker,
  none installable here) and non-spectra deposits such as PXD009604 (XML only).
* Downloads use aria2c, 8 ranged connections per file (one PRIDE connection is about 1 MB/s from here).

## 2026-11-03: Gubbels/Zarringhalam multiomics preprint (unreleased)
* Wu et al. 2026, bioRxiv 10.64898/2026.10.29.762658: 22 10x scRNA-seq libraries (WT RH + 21 transcriptional or
  epigenetic factor knockdowns, 211,052 cells on a common cell-cycle clock), scATAC-seq, CUT&RUN, pull-down MS.
  No accessions: "The data can be explored in the ToxoHub interactive application".
* Einar chose to harvest the ToxoHub v4 JSON API (`scripts/harvest_toxohub_preprint.py`, 1 request/s) into
  `parascape_omics/preprint_results/toxohub_v4`, every file flagged unreleased. **Never put it in a release.**
  When GEO/SRA has the reads, queue them in the uniform pipeline and treat the harvest as superseded.

## 2026-11-03: share spare CPUs with autotrade
* Einar: ParaScape and autotrade share the machine; when autotrade is not using the cores, ParaScape may.
  `omics_resource_guard.py` gained a `wide` mode (`--wide-cpus`, drop-in `60-wide-when-idle.conf`): after 3 min
  with non-ParaScape work under 6 busy cores and RAM under 80 GB it adds CPUs 8-14 and 24-30 and sets
  `extra_workers: 1`; it yields at once above 10 cores or 80 GB. `omics_quant.py` starts
  `PARASCAPE_EXTRA_WORKERS` (1) further workers that claim a run only while the guard reports spare workers.
* Measured queue (automation plan, not the older catalogue plan): 2,316 runs, 361 verified, 7.3 TB of reads
  left; ~135 GB/day at 16 CPUs.
* Harvest finished 2026-11-03 19:24: 5,605 API responses (829 MB; app 4.1.0, data bundle 2026-10-29): 21 regulators
  with whole-cycle and 24-window DE plus targets, 108 communities, 210 regulator pairs, 4,713 gene pages.

## 2026-11-04: single-cell collections out of the bulk queue
* GSE229505 (dual perturb-seq, 97 runs) failed with "no reads written": 10x spots hold 1 cDNA read plus 3
  technical reads. With it, four other per-cell collections still pending (GSE346252, PRJNA1047766, PRJDB12314,
  PRJNA1053644; 125 runs in all) were removed from the queue (`scripts/exclude_single_cell_studies.py`, backup in
  `automation/before-single-cell-exclusion-*`); per Einar they use deposited single-cell matrices instead.
  GSE132238 (16 runs, already verified) keeps its results. Queue now 2,227 runs, 391 verified.
* Correction 2026-11-04 06:30: the running workers kept the old queue in memory and spent ~6 h on GSE229505;
  four of its runs got bulk results through the `fastq-dump --split-3` NAS fallback. Those are not valid bulk
  measurements and were moved to `parascape_omics/rnaseq_excluded_single_cell/` (README there). Workers now
  re-read the plan every 10 min once it changes and skip runs taken out of it.

## 2026-11-04: adapter trimming before quantification (Einar)
* The per-sample quality score found 44 runs in GSE243206, GSE243207 and GSE245775 where salmon mapped
  **0.00 %** while kallisto mapped 10-32 % of the same reads. Cause: the deposits are untrimmed short-insert
  libraries. In SRR26060420, 99.5 % of reads carry the Illumina TruSeq adapter and the 150 nt reads hold a
  ~44 nt insert (mean read length after trimming 49 nt), so salmon's selective alignment rejects everything.
  `--allowDovetail` and single-end mapping both still give 0.00 %; trimming gives 5.10 % against the
  parasite-only reference (the rest are host reads).
* Einar chose: **always trim** with fastp 1.4.0 (`--detect_adapter_for_pe`, `-l 20`), and re-quantify the
  already-verified runs as part of the planned top-up re-download (eq-classes and Gibbs samples), so trimming
  costs no extra download. `omics_quant._trim` records the fastp version, read counts, mean read length before
  and after, and the adapter-trimmed share in `done.json` under `trimming`.
* Runs verified before this change carry no `trimming` key; that is how an untrimmed result is recognised.

## 2026-11-04 evening: ParaScape capped to 8 threads on CPUs 0-7 for plaque v8
The tutorials session relayed a maintainer request: plaque v8 training (unit `plaque-v8-ws-20261104`, GPU, a
single-threaded loop) was running ~7x slower than at Home because WS load was ~33 (three salmon runs, Sage and
an autotrade backtest). ParaScape was reduced without pausing anything:
* `parascape-.service.d/55-share-with-plaque-v8.conf`: `CPUAffinity=0-7` (was 0-7,16-23), `Nice=19`,
  `IOSchedulingClass=idle`. Covers every parascape-* unit, so salmon, kallisto, fastp, Sage and the download
  workers inherit it. Verified in the cgroups.
* The guard's **wide mode is off** (`60-wide-when-idle.conf` rewritten without `--wide-cpus`): it used to claim
  the idle cores whenever the machine looked quiet, which is what would starve v8 again.
* Quantifier: `PARASCAPE_WORKERS=1`, `PARASCAPE_THREADS=4`, `PARASCAPE_EXTRA_WORKERS=0`; Sage keeps 4 threads.
  Eight ParaScape threads in total.
* Cost: RNA quantification goes from 3 concurrent runs on 16 cores to 1 on 4, roughly a quarter of the
  throughput (~88 runs/day to ~20; the remaining 1,746 runs from ~22 days to ~80). Proteomics is download-bound
  and barely affected. **To restore, delete `55-share-with-plaque-v8.conf`, put `--wide-cpus
  8,9,10,11,12,13,14,24,25,26,27,28,29,30` back in `60-wide-when-idle.conf`, and set WORKERS=2, THREADS=16,
  EXTRA_WORKERS=1.**

## 2026-11-04: HDX-MS projects excluded from the uniform re-search
PXD039335 failed twice on refused downloads. The file turned out to be reachable (PRIDE is intermittently
unavailable rather than blocking us), but its filenames revealed the real problem: it is **hydrogen-deuterium
exchange MS**. Four of the 203 selected projects are HDX: **PXD034607** ("Hydrogen Deuterium Exchange of MAF1
and RalGAP", 11.3 GB), **PXD031816**, **PXD023791** and **PXD039335** (102 RAW files in total).
* **Why they must not go through this pipeline:** HDX peptides carry a variable, unknown number of deuteriums,
  so their masses are shifted by amounts the search does not model. Searching them against an undeuterated
  FASTA identifies a fraction of the spectra, and any protein-level quantity would reflect deuteration rather
  than abundance.
* **Measured, not assumed.** The two HDX projects already searched yield **1,272 and 1,883 PSMs per run**
  against **10,717 / 13,316 / 45,527** for the standard bottom-up projects in the same pipeline -- about a tenth.
* PXD034607 and PXD039335 now carry `unsupported.json` with the reason, so the queue skips them rather than
  retrying forever. PXD031816 and PXD023791, already searched, carry `NOT_COMPARABLE.json`: the spectra are kept
  but their PSM, peptide and LFQ tables must not be used as expression or co-purification evidence, and they
  must be left out of the background table.
* **Detection rule for future selections:** HDX is visible in RAW filenames (`HDXMS`, `HDX`) more reliably than
  in the project title -- only PXD034607 says so in its title, while all four say so in their filenames.
* **Confirmed working, and expected to last ~4.5 days (2026-11-04 23:00).** The tutorials session reports plaque
  v8 went from 0.8 to 6.3 windows/s and GPU 13 % to 93 % -- the same as at Home -- so v8 was purely CPU-starved
  and ParaScape was the main cause. At ~2.3 h/epoch it needs ~48 more epochs to its epoch-100 minimum. The cap
  **stays until that session pings**; do not revert it early. RNA quantification runs at ~20 runs/day meanwhile
  (~80 days for the remaining 1,746 instead of ~22). The proteomics queue is download-bound and unaffected.
  A middle ground was offered if the maintainer wants RNA back sooner: 4 threads on CPUs 0-3, which leaves v8
  28 free cores and should keep it near 6 windows/s.
* **Widened to CPUs 0-23 (2026-11-05 07:11), maintainer decision via the tutorials session:** v8 needs only
  1-2 free cores. `55-share-with-plaque-v8.conf` now sets `CPUAffinity=0-23` (Nice=19, idle IO unchanged); the
  guard has fast CPUs 0-23, protected 0-7, and **wide mode stays OFF so cores 24-31 can never be handed over** --
  that is the mechanism that would otherwise creep back onto v8. Quantifier: 3 workers x 6 threads = 18, plus
  Sage 4, so 22 threads inside a 24-core range, slack left on purpose. Verified by walking every thread of both
  ParaScape cgroups: none is permitted on 24-31, and v8 was at 98.6 % CPU immediately after.
* **Correction to the cost estimate:** the 0-7 cap was predicted to cost ~68 runs/day; measured overnight it was
  18 runs in 6 h (~72/day), i.e. ~24 days to finish the queue rather than the ~80 claimed. One run on 4 threads
  is far more efficient than three competing for 16, so stepping back to 0-7 is much cheaper than first stated.
* **Confirmed by measurement (2026-11-05 07:35):** with ParaScape on 0-23, plaque v8 holds 6.26 windows/s,
  2.3 h/epoch, GPU 97 % -- unchanged from its Home speed. **Keep the 0-23 share and wide mode OFF until the
  tutorials session reports v8 finished**, expected around 9-10 November (4.5 days from 06:50 on 2026-11-05),
  possibly longer if its plateau rule continues. If v8 slows at any point, narrow to CPUs 0-7 first and
  diagnose afterwards.

## 2026-11-05: PRIDE flaps, so partial downloads are now resumed
Three projects failed on `ConnectionRefusedError` despite a retry ladder spanning 38 minutes. Diagnosed rather
than assumed: the endpoint is **intermittently** unavailable on a timescale of seconds -- the same directory URL
answered HTTP 206 and then refused within a minute, while a small file from the same project succeeded and the
API stayed up throughout. It is not a wrong URL (ours matches the API's `publicFileLocations` date path), not a
missing file, and not a wholesale block on our address.
* **Fix: aria2c now runs with `--continue=true` and a transport failure keeps the partial file**, so each
  attempt and each later queue pass adds to it rather than starting again. Deleting the partial on failure (as
  this code did until now) threw away exactly the progress that a flapping host allows. A partial file is still
  deleted when the failure is *not* transport -- a wrong size or a truncated plain download -- because then it
  cannot be trusted.
* Aspera would be the robust answer (PRIDE publishes `prd_ascp@fasp.ebi.ac.uk` for every file, a different host
  and protocol) but no client is installed and it is not on the usual channels. Recorded as the next option if
  the resume behaviour proves insufficient.

## 2026-11-05: PRIDE over FTP, nine times faster and reliable — no Aspera key needed
Aspera turned out to be unnecessary. `~/.aspera/connect/bin/ascp` is installed (Connect 4.2.19) and
`fasp.ebi.ac.uk:33001` is reachable, but EBI's `prd_ascp` account needs the openly published
`asperaweb_id_dsa.openssh`, which Connect 4.x no longer bundles; neither the shipped
`aspera_tokenauth_id_rsa` nor the default authenticates. **Before chasing that key, the FTP endpoint PRIDE
advertises beside HTTPS for every file was tested, and it is simply better.**

Measured on the same file at the same moment (`PXD027029/OE03569.raw`, the one that had been failing):

| Transport | Speed |
|---|---|
| HTTPS, single ranged request | 165 kB/s |
| FTP, single ranged request | 560 kB/s |
| **FTP via aria2c, 6 connections** | **9.2 MB/s** (231 MB in 25 s) |

Resume also works over FTP: a second aria2c run continued the partial file rather than restarting it.
`ms_uniform_search.fetch` now derives the `ftp://` form of each PRIDE URL and tries it first, keeping HTTPS as
the fallback for the last attempts and for any URL that is not a PRIDE one. That is roughly a nine-fold
speed-up on the slowest part of the pipeline, and it removes the connection refusals, which came from the HTTPS
endpoint specifically. **Nothing is needed from Einar.**

## 2026-11-09: the plaque v8 reservation is lifted

The tutorials session reported v8 finished at 20:19, so cores 24-31 are no longer reserved and the share set up
on 2026-11-05 is reverted:

* `55-share-with-plaque-v8.conf` **deleted**, so the standing `50-cpu-limit-16.conf` governs again:
  `CPUAffinity=0-7 16-23`, Einar's half-the-machine share of 2026-11-02.
* `60-wide-when-idle.conf`: **wide mode back ON**, fast CPUs 0-7 and 16-23 with 0-7 protected, spare cores
  `--wide-cpus 8-14,24-30` borrowed only while nothing else wants them and given back as soon as other work
  appears; `--quiet-seconds 180` as before the plaque period.
* **A mistake caught on the way, recorded because it would have been invisible:** the first revert set
  `CPUAffinity=0-31`, reading "restore your normal share" as "take the whole machine". The normal share is
  HALF -- 16 logical CPUs -- and the memory file saying so was read only after the change was applied. Reverting
  to a remembered state is not the same as reverting to the documented one.
* Guard restarted and active; it reports `protected` and "waiting for sustained quiet", which is correct with a
  plaque v9 GPU evaluation running on 4 threads -- it will widen on its own once the machine is quiet.

A plaque v9 GPU evaluation runs now and other GPU jobs follow; the tutorials session says none of them needs a
CPU reservation. If any of them slows, narrow to CPUs 0-7 first and diagnose afterwards, as before.
