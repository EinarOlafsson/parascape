# 88 · Proteomics downloads have slowed to a crawl

**Status: open (measured 2026-11-06).** The uniform proteomics re-search fell from ~17 projects a day
(2026-11-05) to ~4.

## Measured
- Per project, the search is not the cost: PXD030932 took 4,687 s in all and Sage 221 s of it.
- PRIDE FTP is now the bottleneck: **0.45 MB/s** straight to local disk (50 MB in 117 s), against
  **9.2 MB/s** measured the day before. Writing to the NAS is not the limit -- 247 MB/s.
- One file of PXD004426 sat at 0 KB/s for over 30 s with aria2c still running.
- `status.json` has read `done=73` since 2026-11-05 noon while prepared projects went 79 -> 83: the counter
  is stale and the progress report has to count folders instead.

## The plan
1. Download several files of a project at once (each aria2c already splits one file six ways; files are
   fetched one after another).
2. Detect a stalled transfer and restart it rather than wait out the timeout.
3. Try PRIDE's other transports again (Aspera, HTTPS) and keep whichever is fastest on the day.
4. Fix the `done` counter.
