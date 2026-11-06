import assert from "node:assert/strict";
import test from "node:test";

import catalog from "../src/navigation_catalog.json" with { type: "json" };
import { SUMMARY_PROMPT, SYSTEM_PROMPT } from "../src/prompts.js";

// The catalog is generated from parascape/modules.py; every module the assistant may name must be
// in it, with what it leads to, and the prompt must carry each row.
test("every ParaScape module is in the catalog with its next steps", () => {
  const keys = new Set(catalog.modules.map((m) => m.key));
  for (const required of ["gene", "map", "network", "infer", "drop", "questions"]) {
    assert.ok(keys.has(required), required);
  }
  for (const m of catalog.modules) {
    for (const next of m.leads_to) assert.ok(keys.has(next), `${m.key} leads to unknown ${next}`);
    assert.ok(SYSTEM_PROMPT.includes(`${m.key} | ${m.name} | `), m.key);
  }
});

test("module links use the parascape scheme and nothing is spaCR", () => {
  assert.match(SYSTEM_PROMPT, /parascape:\/\/module\/<key>/);
  assert.doesNotMatch(SYSTEM_PROMPT, /spacr/i);
  assert.doesNotMatch(SUMMARY_PROMPT, /spacr/i);
});

test("an error answer has the four sections and ends with a cause line", () => {
  for (const heading of ["What happened", "Likely cause", "How to fix it", "Is this a ParaScape bug?"]) {
    assert.ok(SYSTEM_PROMPT.includes(`**${heading}**`), heading);
  }
  assert.match(SYSTEM_PROMPT, /CAUSE: <kebab-case-slug>/);
});
