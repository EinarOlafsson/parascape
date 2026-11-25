# Reading rules: the free Toxoplasma Affinage run, Stage 1

These are adapted from Di Bernardo and Cheeseman 2026 (arXiv 2607.02217), Stage 1. You read each abstract ONCE and record findings for each listed ParaScape gene that the paper names.

## What counts as a finding

Keep a finding ONLY if THIS paper reports a direct experiment establishing it for THIS gene product.

**Admit:**
- **Localisation** by imaging (IFA, an endogenous tag, immuno-EM, expansion microscopy) or by fractionation.
- **Partners** by co-IP, pulldown, BioID or TurboID, cross-linking, yeast two-hybrid, or reconstitution. If a partner is named, put its ParaScape gene id in `partner_gene_ids` when it is one of the listed genes; otherwise give its name in `partner_names`.
- **Complexes** shown biochemically or structurally.
- **Activity**: an enzymatic activity from an in-vitro assay or active-site mutagenesis, or a molecular function shown directly (binding, transport, phosphorylation of a named substrate).
- **Phenotype of loss or depletion**: knockout, knockdown, AID, conditional KO or mutant. The phenotype must be measured; name the step affected (invasion, egress, replication, motility, gliding, differentiation, virulence in mice, …).
- **Secretion and host-facing**: the protein is secreted into the PV, the PVM or the host cell, or acts on a named host protein or pathway, shown by experiment.
- **Structure**: an experimental structure (crystal, cryo-EM, NMR).

**Exclude:**
- expression or transcript-level correlation only;
- association or omics-list membership only;
- a prediction (domain, AlphaFold, homology, "putative");
- a claim that the paper only raises, proposes or discusses;
- findings about a different gene that only shares a name;
- findings in another organism, unless the abstract shows it for the Toxoplasma protein.

**Status:**
- `established`: the experiment supports it.
- `refuted`: the paper rules it out. Keep refutations; they are informative.

**No findings:** if a paper has no admissible finding for a listed gene, write NOTHING for that gene. Never pad.

**The listed gene is the subject.** Every `finding` sentence has the listed gene as its subject. For example, write "ASP5 is required to export GRA84 to the host nucleus", not "GRA84 export depends on ASP5". A result about a protein that is NOT on the paper's gene list is not written under another gene; put it in the doubts file instead (see below).

**Named in the abstract.** The abstract itself must name the protein: by a name, an id, or an unambiguous description. Never credit a listed gene with a finding because its product name fits a complex or pathway the abstract describes in general terms. For example, "we purified ATP synthase" does NOT give findings to every ATP synthase subunit on the list (2026-11-24: seven such findings were excluded).

**Same gene, other ids.** `TGGT1_`, `TGME49_` and `TGVEG_` ids with the same number are the same gene, so `TGGT1_226250` is `TGME49_226250`. Use the product description to match a name that the abstract uses for a listed gene. Write a finding only when the match is clear.

**Ectopic expression.** Expressing a parasite protein in host or cancer cells, rather than having the parasite deliver it, does NOT establish host_facing. Do not write it.

**Supplement batches (`supp_NNN.json`).** A paper there was read before for the genes in `already_read`; its `genes` list holds only the genes that the first pass left out. Write findings ONLY for the genes in `genes`. A partner may be one of the `already_read` genes: give its id in `partner_gene_ids`. Do not repeat findings for the `already_read` genes.

**Doubts file.** Write one line per doubt to `findings/<batch>.doubts.jsonl`:
```
{"pmid": "...", "gene_id": "... or \"\"", "doubt": "wrong_gene_match|unlisted_protein|ambiguous_match|other", "note": "one sentence"}
```
Use it when:
- a listed gene looks matched to the wrong protein (`wrong_gene_match`);
- the paper establishes something about a protein with no listed gene (`unlisted_protein`, naming the protein);
- you cannot tell which listed gene a name refers to (`ambiguous_match`).


## Output (one JSON object per line, in the batch output file)

```
{"pmid": "...", "gene_id": "TGME49_...", "kind": "localisation|partner|complex|activity|phenotype|host_facing|structure",
 "finding": "one plain declarative sentence, as the paper establishes it",
 "method": "e.g. IFA of endogenous HA tag; co-IP; AID knockdown plaque assay",
 "method_tier": 1|2|3,
 "status": "established|refuted",
 "compartment": "one of the ParaScape compartments below, or \"\"",
 "partner_gene_ids": [], "partner_names": [],
 "stage": "tachyzoite|bradyzoite|sexual|sporozoite|\"\"",
 "open_question": "what the abstract says is still unknown, or \"\""}
```

**Method tier:**
- 1 = genetic loss-of-function with a measured phenotype, an in-vitro activity with mutants, or a structure;
- 2 = imaging, co-IP, pulldown or proximity labelling;
- 3 = anything else admitted.

**ParaScape compartments** (use the closest; use "" if none fits):
19S proteasome, 20S proteasome, 40S ribosome, 60S ribosome, ER, Golgi, IMC, PM - integral, PM - peripheral, apical, apicoplast, cytosol, dense granules, endomembrane vesicles, micronemes, mitochondrion - membranes, mitochondrion - soluble, nucleolus, nucleus - chromatin, nucleus - non-chromatin, rhoptries, tubulin cytoskeleton, plus the extra terms PV lumen, PVM, host cytosol, host nucleus, conoid, basal complex, apical annuli, residual body, cyst wall.
