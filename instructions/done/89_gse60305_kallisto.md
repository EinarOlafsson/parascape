# 89 · GSE60305: kallisto fails on every run

**Status: resolved 2026-11-06 (moved to done).** The cause was not kallisto: GSE60305 was sequenced on an
**AB SOLiD 4**, and its reads are **colour space** (`T13302220001300...`). No k-mer quantifier can read
them -- kallisto exited 1 because nothing pseudoaligned, and salmon "succeeded" with **0 of 46,921,896**
reads mapped (SRR1542935), which the pipeline counted as success. Fixed in `scripts/omics_quant.py`
(deployed to the worker; active at its next restart): colour-space instruments are recorded in
`colour_space.json` before any download and never retried; salmon mapping 0 reads is a failure; a kallisto
failure after salmon mapped reads keeps the run on salmon, recorded in `done.json` (`kallisto_failed`). The
18 GSE60305 runs are recorded. The study's own processed table, `GSE60305_readCount.txt.gz` (8,072 ToxoDB
genes x 18 samples: RH, PRU, CTG at pH 7 and 8, three replicates), is already downloaded as a
`processed_table_candidate` and enters through processed-table admission. Tests:
`tests/test_omics_quant_colour_space.py`.

Found 2026-11-06. All 18 runs of GSE60305 fail at `kallisto quant` (single-end 76 nt,
`-l 200 -s 30`), each retried about three times -- 49 of one day's attempts. Salmon succeeds on every one
and its output is on disk; the pipeline marks a run failed when either quantifier fails, so the study is
absent from the release. Retries stopped on their own.

## The plan
1. Capture kallisto's own error, which the log does not keep.
2. Fix the cause if it is ours; otherwise accept the salmon quantification for this study, recorded as such
   per run.
3. Make the pipeline keep a run whose salmon succeeded when kallisto fails, with the failure recorded, so one
   tool cannot lose a whole study again.
