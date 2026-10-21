# Unattended download and processing queue

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
