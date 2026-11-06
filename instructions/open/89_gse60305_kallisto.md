# 89 · GSE60305: kallisto fails on every run

**Status: open (found 2026-11-06).** All 18 runs of GSE60305 fail at `kallisto quant` (single-end 76 nt,
`-l 200 -s 30`), each retried about three times -- 49 of one day's attempts. Salmon succeeds on every one
and its output is on disk; the pipeline marks a run failed when either quantifier fails, so the study is
absent from the release. Retries stopped on their own.

## The plan
1. Capture kallisto's own error, which the log does not keep.
2. Fix the cause if it is ours; otherwise accept the salmon quantification for this study, recorded as such
   per run.
3. Make the pipeline keep a run whose salmon succeeded when kallisto fails, with the failure recorded, so one
   tool cannot lose a whole study again.
