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
