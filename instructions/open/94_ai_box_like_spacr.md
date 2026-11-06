# 94 · The AI box with spaCR's choices, and a free ParaScape AI

**Status: built 2026-11-06.** The AI box lists spaCR's four choices (Claude, ChatGPT, Gemini -- shown as
"not installed" with install and sign-in steps when absent -- and ParaScape AI (free)), selecting the first
that can answer. The free one asks `POST /diagnose` on the deployed `parascape-relay` worker
(`services/parascape-relay`, Gemma 4 26B on Workers AI, 4000-neuron daily cap so spaCR keeps most of the
account allowance) and is told what is on screen. Client: `parascape/relay.py`, `parascape/chat.py`. Einar: "the AI box should have the options that spaCR does. And one
should be the free ParaScape AI."

## What spaCR offers (`spacr/qt/ai/providers.py`)
Claude (via Claude Code), ChatGPT (via Codex CLI), Gemini (via Gemini CLI), and **spaCR AI (free)** -- a hosted
provider served by Einar's Cloudflare worker (`services/spacr-first-responder`, Workers AI, daily quotas).
Every provider is listed whether installed or not; an uninstalled one offers its install and sign-in.

## The plan
1. ParaScape lists all four, not only those on PATH; choosing an absent CLI shows how to install / sign in.
2. **ParaScape AI (free)**: a worker of the same shape (`services/parascape-ai`, Workers AI, per-install and
   daily caps, no account needed) and a client provider. Deploying it publishes a service on Einar's
   Cloudflare account, so it is deployed only when Einar says so; until then the entry says it is not yet
   switched on.

## Gate
The AI box offers the same four kinds of choice as spaCR, and the free one answers.
