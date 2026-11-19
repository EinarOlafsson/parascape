# 81 · Re-quantify the earlier runs: trimming, equivalence classes and Gibbs samples

**2026-11-19: the queue was stalled for nine hours by three DEADLOCKED kallisto processes.**
The hourly log read `624 new runs still to quantify` from 21:52 to 06:52 without moving. The cause
was not slowness: three `kallisto quant` processes against the `ME49_host` index, on
PRJDB16831/DRR512193, DRR512199 and DRR512194, had each been running 9.5 hours at ~595 per cent CPU
having **read 2.11 GB (the 2.04 GB index) and written zero bytes**, with their output directories
empty and their kernel state `futex_wait_queue` -- threads spinning on a lock. Eighteen cores were
busy producing nothing and the worker was waiting on them.

Einar authorised clearing them. A plain SIGTERM was enough for all three; the worker picked up new
runs within a minute. Nothing was lost, because nothing had been written.

**2026-11-19, the cause found: kallisto deadlocks on the COMBINED index in `--single` mode.**
Clearing the first three did not help -- the worker took the next three runs of the same study
(PRJDB16831 DRR512206, DRR512209, DRR512212) and they hung identically within minutes. The study
has six runs, none with a `done.json`, so the queue could never get past it.

What distinguishes it, measured:

* every study that has completed is **PAIRED**; PRJDB16831 is **single-end, 36 bp**;
* the same 100,000 reads, the same `--single -l 200 -s 30`, the same single thread:
  against `references/ME49/kallisto.idx` (0.02 GB) kallisto finishes in under a minute --
  100,000 reads processed, 6,518 pseudoaligned, EM in 289 rounds, all three outputs written;
  against `references/ME49_host/kallisto.idx` (2.04 GB) it hangs with zero bytes written;
* it is NOT a threading bug: `-t 1` hangs exactly as `-t 6` does.

So the trigger is the combined ME49_host index in single-end mode, not the data, the parameters or
the thread count. Note that only 6.5 per cent of these reads are parasite, which is why these
samples need the combined index at all -- the ME49-only result is not a substitute, it is the
control that isolates the fault.

**What this needs.** Either the worker skips kallisto against ME49_host for single-end libraries
and keeps salmon, which was running normally throughout, or PRJDB16831 is quarantined so six runs
cannot hold 624. Einar's call; nothing has been changed in the worker.

**What to watch for.** The symptom is a flat queue count with high CPU, which looks like a slow run
and is not. The check that distinguishes them is `write_bytes` in `/proc/<pid>/io`: a working quant
writes steadily, a deadlocked one has read its index and written nothing. If those three runs hang
again on retry, they should be quarantined rather than retried indefinitely -- the queue cannot
make progress while they hold the worker.

**2026-11-18, found while building figure 2 of the paper: GSE245775 quantifies to ZERO.**
"Differentiation ribosome profiling (eIF1.2)", 24 runs, and every run sums to zero parasite TPM
against BOTH `ME49` and `ME49_host`. A ribosome-profiling study of Toxoplasma with no parasite
reads is a pipeline result, not a study property, and the shape matches the adapter problem this
instruction is about: ribosome-footprint reads are ~28-30 nt and still carry adapter, so an
untrimmed read maps nowhere. It is the same failure mode as GSE60305 (instruction 89), different
cause. It should go in the re-quantification queue, and until it does, anything reading its
matrices gets 8,157 zeros rather than an error.

**Status: ready, waiting on the downloads (2026-11-06).** `scripts/omics_quant.py requantify` (deployed):
waits hourly until every planned run is quantified, then re-quantifies the 456 runs verified before
trimming -- the zero-mapped ones (GSE243206, GSE243207, GSE245775) first -- into `requant/` beside each
verified result, and writes `requant/compare.json` (Spearman of TPM per reference; "agrees" at >= 0.9,
"recovered" for a run that mapped nothing before, "DIFFERS -- review" otherwise; nothing is overwritten).
The service file is `scripts/systemd/parascape-requant.service`; **installing it is Einar's step** (a new
persistent service was not installed without his say-so). On 2026-11-17 the markers say 1,446 of 1,470 runs are done (98.4 per cent); the figure below is the one recorded on 2026-11-06, kept for the record: 723 of 999 runs done, 1,474
priority-plan runs pending. Tests: `tests/test_omics_requant.py`.

## Why three things accumulated
1. **Adapter trimming** arrived on 2026-11-04. Runs verified before it carry no `trimming` key in `done.json`,
   which is how an untrimmed result is recognised. Normal libraries lose almost nothing (0.0-0.1 % of reads), but
   40 runs of GSE243206, GSE243207 and GSE245775 mapped 0.00 % until trimming and are excluded from the
   expression columns until this runs.
2. **Equivalence classes and 30 Gibbs samples** (`--dumpEqWeights --numGibbsSamples 30`) arrived the same week,
   so runs can be re-quantified against an updated annotation and carry per-transcript uncertainty without their
   reads. Earlier runs lack both.
3. The reads are deleted after quantification, so this is a re-download.

## The plan
Re-download and re-quantify the runs verified before each change, appending the missing outputs **without
altering a verified result**: a new quantification is written beside the old one and the old numbers stay until
the new ones are checked against them. Einar's instruction stands: do this only after the new data is
downloaded. The 40 zero-parasite runs come first, because they currently contribute nothing at all.

## Gate check, measured 2026-11-17 23:30

**The gate is NOT open, and the reason is scheduling rather than anything wrong with the method.**

`parascape-requant.service` has been up for six days and has re-quantified **none** of the three
short-insert series. By design it waits for the new-run queue first, and that queue is draining
slowly: the unit logged 773 runs still to quantify at 17:50 and 748 at 22:51, about five an hour.
At that rate the re-quantification does not begin for roughly six days, so the 40 zero-parasite runs
that this plan says "come first" are in fact behind 748 others. Nothing is broken; the priority in
the plan is simply not expressed in the queue.

**The premise of the plan is confirmed.** Measured over all 48 runs of GSE243206, GSE243207 and
GSE245775, against ME49:

| | these 48 runs | a random 120 runs elsewhere |
|---|---|---|
| salmon mapping | 0.00 per cent in 48 of 48 | zero in 3 of 120, median 58.0 per cent |
| kallisto pseudoalignment | median 6.95 per cent, range 0.2 to 49.2 | median 66.9 per cent |

Salmon is not broken: it maps normally everywhere else, including on untrimmed runs (median 57.0 per
cent over 37 untrimmed runs sampled). What is specific to these three series is heavily adapter-laden
reads, which this file already records as adapter in 99.5 per cent of reads. Kallisto still maps a
little because pseudoalignment is permissive; the stricter selective alignment of salmon rejects them
outright, which is why the rate is exactly zero rather than merely low. Trimming is the right fix.

**Ready for when the queue reaches them:**
- `results/requant_gate_baseline_2026_11_17.csv` - every rate for all 48 runs across all four
  references, recorded BEFORE re-quantification, which is the only time a baseline can be taken.
- `scripts/check_requant_gate.py` - reads only, never writes into a run directory. It reports a run
  whose rate FALLS by more than the tolerance rather than overwriting it, as the gate requires.
  It currently prints `GATE: NOT OPEN, 48 of 48 runs have not been re-quantified`.

**A decision for Einar:** leave the queue to reach these runs in about six days, or let the
re-quantification take the 48 first as the plan says. Nothing was changed in the running service.

## Gate
The 40 runs produce non-zero parasite mapping; every re-quantified run agrees with its previous values within
the expected tolerance, and any that does not is reported rather than overwritten.

## Installed (2026-11-08, Einar: "Yes, install it")

`parascape-requant.service` is installed as a user unit and running (`systemctl --user status parascape-requant`),
with MemoryMax lowered to 24G and 8 threads (the shared machine RAM rules: long jobs together under ~80 GB); the
RAM ledger has its line. As designed it waits for the new downloads first: at start, "1313 new runs still to
quantify". It re-checks hourly, then re-quantifies beside each verified result (the 40 zero-parasite runs first) and
writes `requant/compare.json`; nothing is overwritten. The gate (non-zero parasite mapping for the 40, every run
agreeing within tolerance or reported) is checked when it has run.
