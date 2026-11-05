# 82 · Fingerprints from a literature meta-analysis: gene-set signatures

**Status: open (requested 2026-11-05).** The second kind of fingerprint, and the one that bears directly on the
host-facing failure.

## The request (Einar)
"Fingerprints can be based on a meta-study analysis, like parasitophorous membrane host-facing parasite proteins
and parasite proteins in host cytosol. A meta-analysis of the literature can make those annotations and a
fingerprint can then be developed."

## Yes — and it runs along the other axis of the matrix
Instruction 75 fingerprints a **condition** and annotates **samples**: built from columns that differ in that
condition, used to say where a dataset sits on an axis. This one fingerprints a **gene set** and annotates
**genes**: built from the data profiles of proteins known to share an annotation, used to say how well an
unannotated protein matches that group. Same word, opposite direction, and they validate differently:

| | Built from | Annotates | Validated by |
|---|---|---|---|
| Condition fingerprint (75) | samples differing in a condition | datasets and samples | a held-out study that states its own condition |
| **Set fingerprint (82)** | genes known to share an annotation | **genes** | **held-out members of the set** (the harness, unchanged) |

Mechanically a set fingerprint **is** the lookalike strategy of instruction 77: rank the genome by resemblance
to a seed set. What this item adds is the seed sets, and they have to come from the literature because no single
experiment defines them.

## Why this is the highest-value version of the idea
Host-facing failed for a specific reason: the evidence could not separate "stays in the vacuole" from "exported
into the host", and the available negatives were markers of compartments *inside* the parasite, which answer an
easier question. The sets Einar names are exactly that missing distinction. A meta-analysis should produce a
**compartment ladder**, each rung its own label with its own fingerprint:

1. injected at the moving junction (rhoptry contents delivered on entry)
2. PVM, **vacuole-facing** side
3. PVM, **host-facing** side
4. vacuole lumen and the intravacuolar network — **the hard negatives**
5. exported to the host cytosol
6. host nucleus
7. host organelle or host plasma membrane

Rungs 2 and 3 are the pair that no current assay distinguishes, and rung 4 is the negative set whose absence
capped host-facing precision at 0.822.

## The plan
1. **Mine and grade, do not pool.** Every assignment carries its evidence grade (direct microscopy or
   fractionation / proximity labelling / inferred from phenotype) and a verified PMID, exactly as the 140-claim
   review did. The review of those 140 is the template, and its lesson is the reason for the ladder: 22 proteins
   the PVM proximity screen called host-facing are vacuole-resident, and only reading the papers separated them.
2. **One fingerprint per rung**, built from the data profiles of its direct-evidence members only, scored by
   held-out membership through the existing harness, with the rungs' confusion with each other reported — the
   useful output is "this protein resembles rung 3 and not rung 4", which is the question that matters.
3. **Report confusion, not just accuracy.** A fingerprint that cannot separate rungs 2 and 3 should say so;
   that is a finding about the evidence, and it is what the host-facing run already showed numerically.
4. **Expect ranking, not filling.** A literature meta-analysis realistically yields a few hundred annotated
   proteins. The harness threshold for a frozen test is 400 labels, so most rungs will be `no_frozen_test`:
   they may rank and prioritise, and may not write cells. Stated up front so the result is not a surprise.

## Gate
Each rung ships when held-out members are recovered above the best baseline, **and** the confusion with the
adjacent rung is reported. The ladder succeeds as a whole if rung 4 gains enough members to lift host-facing's
attainable precision above 0.822 — the number that currently makes that question unfillable at any label count.
