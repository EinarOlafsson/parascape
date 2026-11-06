# 97 · Errors become ParaScape GitHub issues, and a routine fixes them

**Status: built 2026-11-06.** Worker deployed (`https://parascape-relay.einar-olafsson.workers.dev`,
D1 `parascape-support`); GitHub App `parascape-first-responder` (App ID 5647965, installation 179382040)
created by Einar, key set by Einar. Gate met: a test report filed issue #1 as app/parascape-first-responder
(closed). Client: `parascape/report_dialog.py` (asks once per error, shows the report word for word,
Send / Don't send / always; Help ▸ Report a problem; Help ▸ Error reports ask/always/never) installed from
`app.main`. Routine: Einar chose a check in the working session -- every 4 h (session cron, 7-day expiry,
re-created when it lapses). Einar: "like in spaCR we need to set up a mechanism where errors
get sent as ParaScape GitHub issues, and you routinely check the GitHub issues for ParaScape and fix them."

## What spaCR does
`spacr/crashreport.py` builds an anonymised report of an error (traceback with paths and user names
removed, module, versions) and posts it to Einar's Cloudflare worker `spacr-first-responder`
(`services/spacr-first-responder`): `POST /report` files it as an issue on `EinarOlafsson/spacr` through
the `spacr-first-responder` GitHub App, de-duplicated by an error signature (a repeat comments on the
open issue), with per-install and daily caps; `POST /diagnose` is the free hosted AI. The App private key
lives only in the worker's secret store.

## The plan
1. **One ParaScape worker** (`services/parascape-relay`, from spaCR's): `/report` (files on
   `EinarOlafsson/parascape`), `/diagnose` and `/ask` (the free ParaScape AI of instruction 94), the same
   quotas and anonymisation. Deployed on Einar's Cloudflare account (approved 2026-11-06 for the AI).
2. **GitHub access -- needs Einar**: either install the existing `spacr-first-responder` GitHub App on the
   `parascape` repository, or create a `parascape-first-responder` App; then set its private key as the
   worker secret himself (`wrangler secret put GITHUB_APP_PRIVATE_KEY`). Claude does not handle the key.
   Until then `/report` answers "not configured" and the client keeps the report locally.
3. **The client**: an uncaught exception (and the console's ERROR lines) opens a small dialog -- what
   happened in plain words, the exact anonymised text that would be sent, "Send report" / "Don't send" /
   "Always send". Nothing leaves the machine without that consent. Help ▸ Report a problem sends one by
   hand.
4. **The routine**: a scheduled Claude routine reads open issues labelled `error-report` on
   `EinarOlafsson/parascape`, reproduces each, fixes it on nightly with a test, comments on and closes the
   issue (or labels it `needs-einar` when it cannot be reproduced or needs a decision).

## Gate
A forced error in ParaScape files one issue (a second identical error comments on it), and the routine
closes a real one with a fix and a test.
