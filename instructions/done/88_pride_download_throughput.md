# 88 · Proteomics downloads have slowed to a crawl

**Status: done 2026-11-06 (moved to done).** `scripts/ms_uniform_search.py` (deployed; used from the next
project): two files of a project at once (2 x 6 connections, within PRIDE's tolerance); aria2c gives up on a
transfer below 30 KB/s (`--lowest-speed-limit`, `--timeout=60`) and the next attempt resumes the kept
partial; before each project an 8 MB probe of FTP and HTTPS keeps the faster -- measured at once: **FTP 0.14
MB/s, HTTPS 0.71 MB/s**, the reverse of 2026-11-05. `scripts/ms_uniform_queue.py`: the done counter is
counted at every status write (deployed; active at the service's next restart). Tests:
`tests/test_ms_uniform_transport.py`.

Measured 2026-11-06. The uniform proteomics re-search fell from ~17 projects a day
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
