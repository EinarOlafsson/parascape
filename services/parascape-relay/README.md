# parascape-relay

A Cloudflare Worker that lets ParaScape users without a GitHub or AI account

- send an **anonymised error report**, filed as a GitHub issue on `EinarOlafsson/parascape` by the
  `parascape-first-responder` GitHub App (App ID 5647965, installation 179382040), de-duplicated by an
  error signature (a repeat adds a comment to the open issue), and
- use the **free ParaScape AI**: `POST /diagnose` answers an error report, or a plain question when the
  report carries no exception.

Adapted from spaCR's `services/spacr-first-responder` (same routes, quotas, anonymisation and tests);
the prompts and the module catalog are ParaScape's. `src/navigation_catalog.json` is generated from
`parascape/modules.py` -- regenerate it when modules change.

Deployed at `https://parascape-relay.einar-olafsson.workers.dev` on Einar's Cloudflare account, D1
database `parascape-support`. Both relays draw on the account's one free Workers AI allowance;
`AI_DAILY_NEURON_CAP` is 4000 here so spaCR keeps most of it.

## Routes

`GET /health`, `GET /allowance?install_id=<uuid>`, `POST /report`, `POST /diagnose`. Request and
response shapes are spaCR's (see its README); the version field is `parascape_version`.

## Owner setup

Done by Claude on 2026-11-06: `wrangler d1 create parascape-support`, the schema, `wrangler deploy`,
`IP_HASH_SALT`. Done by Einar: the GitHub App. **Einar only** (the key never passes through Claude):

```sh
cd services/parascape-relay
openssl pkcs8 -topk8 -nocrypt -in ~/Downloads/parascape-first-responder.*.private-key.pem \
  -out /tmp/parascape-key.pkcs8.pem
npx wrangler secret put GITHUB_APP_PRIVATE_KEY < /tmp/parascape-key.pkcs8.pem
rm /tmp/parascape-key.pkcs8.pem
```

Redeploy after a change: `npx wrangler deploy` (secrets are kept). Tests: `node --test test/*.mjs`.
