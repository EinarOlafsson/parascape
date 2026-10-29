# Unattended download and processing queue

**Resumed October 27:** the system disk had only 18 GiB free. The queue's own
78 GiB cache was copied to `/media/carruthers/mnt1/parascape_omics_download_cache`,
all 22 ready archives were checked against their recorded SHA-256 digests, and
`~/.cache/parascape-downloads` now links to that directory. Existing commands and
state paths remain usable. The original disk copy was removed after verification.

Claude had separately resumed the original `parascape-quant.service`. The supervisor
recognizes that active service and avoids starting a second writer. Status now names
the actual `processing_service`; active legacy outputs are excluded from audits
until complete. `parascape-legacy-handover.service` preserves the current run, then
switches to the cache-consuming worker at its completed-run boundary. If the original
worker exits earlier, the supervisor recovers through the normal cached worker.
The latest CPU limit remains eight. The queue timer has been re-enabled, PAUSE removed,
and downloads observed transferring. Read the current status rather than old handoff
snapshots, which intentionally retain their historical stopped state.

## What is processed and where the records are

The fixed RNA-seq queue both downloads and processes each planned study/run entry.
It validates downloaded SRA archives (or checks ENA FASTQ MD5s), expands reads when
needed, and runs Salmon and kallisto for every reference named in that entry's plan.
The stored results contain estimated counts, TPM and effective transcript lengths.
Raw reads are deleted after processing; downloading ahead of the worker is a bounded
cache rather than permanent raw-data archiving. Differential-expression tests and
cross-study comparisons are separate downstream work, not the queue's completion gate.

Results live under `/nas_mnt/carruthers/Einar/parascape_omics/rnaseq/<study>/<run>/`:

- `done.json` records source/sample metadata, source URLs and ENA checksums from the
  frozen plan, requested references, tool versions and the finish time.
- `salmon_<reference>/quant.sf.gz` and `kallisto_<reference>/abundance.tsv.gz` are the
  quantification tables; every requested tool/reference combination is retained.
- Salmon's `cmd_info.json` and `aux_info/meta_info.json` record invocation settings,
  index sequence/name hashes, bias correction, mapping metrics, timing and diagnostics.
- kallisto's `run_info.json` records its complete command, tool/index version and
  read-processing statistics. Single-end runs currently use the explicit fragment
  assumptions `-l 200 -s 30`; these are visible in the command records.

Reference downloads, source URLs and checksums are in
`/nas_mnt/carruthers/Einar/parascape_omics/references/<reference>/reference.json`.
The chronological command log is `logs/omics_quant.log` under the same NAS root.
The supervisor's `automation/verified.json` records hashes of completed output
tables and receipts. These are distributed processing records, not one consolidated
per-dataset report. Integrity verification is separate from biological quality review.

**Latest setting, October 22 evening:** user requested an eight-core limit. All queue
services now use CPU affinity 0–7, including the live worker and its child threads.
The guard keeps both CPU sets at 0–7, preventing an automatic increase. Tool threads
change to eight at the next completed-run boundary; the current run is preserved.
Other storage, memory and automatic recovery settings retain their previous values.
This supersedes the earlier full-CPU setting below.

**Current resource setting, October 21 evening:** user ended tutorial priority and
explicitly requested full-machine processing without added model-token use. The
worker has access to all 32 logical CPUs and a live 112 GiB memory cap. The next
completed-run restart uses 32 threads, normal CPU/IO priority and scratch on
`/media/carruthers/mnt1/parascape_omics_scratch`. The deterministic boundary service
`parascape-full-capacity` preserves the current run; its state is in
`~/.cache/parascape-downloads/full-capacity.json`. Full-resource guard mode keeps
all CPUs available. Downloads have no application bandwidth limit, a 200 GiB cache
and a 25 GiB disk reserve, yielding only below 4 GiB available RAM. These settings
supersede the earlier conservative limits below. Installed processing tools expose
CPU threading and no GPU execution option. No model or agent is started by this mode.

The existing queue runs through local user-systemd services. It does not use a chat
session, model tokens, an API key, or another AI agent. The source plan and worker
code are frozen under the user's home directory. Source selection and analytical
methods are not expanded by the supervisor.

The frozen plan contains **1,011 study/run entries and 999 distinct raw run
accessions**. Twelve run accessions appear under two study accessions. Each
study/run entry needs its own planned outputs and provenance; these entries must
not be collapsed when checking completion.

## Lifecycle

`parascape-queue.timer` runs the supervisor every two minutes, with persistent
scheduling after downtime. User lingering is already enabled on this machine, so
logging out does not terminate the user service manager. When the machine and
storage return after a reboot, the supervisor resumes unfinished work.

The supervisor starts missing download, processing and resource-guard services.
Worker crashes retain their existing systemd restart policy. A successful process
exit with unfinished entries also triggers another pass. Repeated short failed
passes use durable backoff, capped at one hour between starts, with no maximum
retry count. Missing storage, unavailable sources or low disk space remain waiting
states and cannot produce a completion receipt. Neither inactivity nor a long
CPU-active run is treated as proof that a worker should be killed.

The supervisor verifies each completed entry against the frozen plan, both tool
versions and every requested reference. Every expected compressed table must
contain its correct header and at least one row, pass full gzip integrity checks,
and receive a SHA-256 checksum. Unchanged verified files reuse cached evidence;
file size, modification time or change time invalidates it. The supervisor waits
for active runs to finish their cleanup before counting them as verified.

A malformed completion marker or a marker whose expected outputs are missing or
corrupt is renamed to `done.invalid-<timestamp>.json`. The original marker and
outputs are preserved, and the existing worker can retry that entry. New worker
completion markers are published atomically. Shared cached archives are retained
until all study entries requiring them are complete, or consumed by the existing
worker's normal expansion/cleanup path.

Only after **all 1,011 entries** pass verification does the supervisor write
`completion.json`, with the frozen-plan hash and verified-manifest hash. It then
stops its three owned workers and its timer. It never controls tutorial processes,
other sessions' services, GPU allocation, or the separate inference action list.

The existing ten-thread setting, adaptive five-/ten-CPU affinity guard, memory
cap, download bandwidth cap, cache limit and free-space reserve still apply.
CPU pressure is paired with measured host CPU headroom: when at least 25% of CPU
time is idle, constrained worker threads cannot keep their own affinity limit
permanently reduced through high global PSI alone. Tutorial detection, IO and
memory pressure still independently require yielding.
The supervisor itself has a 256 MiB cap, low CPU/IO priority and CPU affinity 15.
Initial verification reads the existing completed outputs; later checks reuse
unchanged evidence.

## Inspect progress

State directory: `~/.cache/parascape-downloads/automation/`.

- `status.json`: verified count, unfinished entries and their latest audit reasons.
- `recovery.json`: restart/backoff state.
- `verified.json`: file signatures and verified output hashes, keyed by study/run.
- `completion.json`: created only for the verified complete plan.
- `plan.jsonl` and `config.json`: frozen source plan and its checksum.

```sh
systemctl --user status parascape-queue.timer parascape-queue.service
cat ~/.cache/parascape-downloads/automation/status.json
journalctl --user -u parascape-queue.service -n 20
```

Worker units are started by the queue supervisor rather than independently
enabled at boot. This lets the supervisor recognize a finished queue without
starting processing again. Runtime code lives in
`~/.local/share/parascape-download-worker/`; changing repository files alone does
not upgrade an already-running worker. Preserve its current run when redeploying.

## Operator pause and resume

To pause this automation explicitly, create its pause marker and stop only these
owned units:

```sh
touch ~/.cache/parascape-downloads/automation/PAUSE
systemctl --user stop parascape-queue.timer parascape-queue.service \
  parascape-downloads.service parascape-quant-cached.service \
  parascape-resource-guard.service
```

The marker prevents an accidentally triggered supervisor from restarting workers.
To resume:

```sh
rm ~/.cache/parascape-downloads/automation/PAUSE
systemctl --user start parascape-queue.timer parascape-queue.service
```

## Validation

Run the transport, scheduler and orchestration tests without Qt or GPU setup:

```sh
python -m unittest discover -s tests -p 'test_omics_*.py' -v
systemd-analyze --user verify scripts/systemd/parascape-queue.service \
  scripts/systemd/parascape-queue.timer
```

The synthetic checks cover range resume, capacity and corruption, conservative
CPU restoration, retries after successful incomplete exits, bounded backoff with
unlimited attempts, active-run protection, study aliases, verified completion and
shutdown, post-completion reboot behavior, and explicit operator pauses.


As of 2026-10-29, the resnet_1 relocation is checksum-verified and complete;
the download cap is 200 GiB, with the 25 GiB free-space reserve and eight-core
policy retained. OS signal-paused workers are now reported as `signal_paused`;
the supervisor does not override external operator pauses.

The PyPI expression snapshot is separate from the live archive. Its build and
coverage are documented in [RNA expression](rna-expression.md).


The 2026-10-29 priority override uses `automation/priority-plan.jsonl` for execution
only, with the same 1,011 entries reordered to put PRJNA1008424 first. The signed
audit plan stays unchanged. `parascape-priority-handover.service` switches the
worker after its current completed run; downloader priority is already updated.
