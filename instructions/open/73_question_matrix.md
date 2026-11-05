# 73 · The question matrix: 100–1000 trait questions, each with an answer, a probability and an uncertainty

**Status: open (requested 2026-11-05).** Supersedes the framing of instruction 72, which treated host-facing,
function and partners as *the* three traits. They are three questions among hundreds.

## The request (Einar, 2026-11-05)
* "I do not only care about host-facing parasite proteins. I care about systematically generating as much
  knowledge with the information available today on Toxoplasma from the Toxoplasma literature. Whether a protein
  is host-facing or not is just one trait. There are hundreds or thousands of information columns I would like to
  systematically predict, with assigned uncertainty to each prediction."
* His examples: importance for nutrient acquisition, and under low nutrient conditions; cell division; host-cell
  invasion; egress; switching to bradyzoite and from sporozoite; dissemination; survival in the liver;
  bradyzoite survival; a favourable immune response in mice, and in humans; the cell cycle; translation; protein
  scavenging; transcription; autophagy; methylation; acetylation; phosphorylation; lipid biogenesis and
  scavenging; the most likely binding partners; what a protein is expressed with; what it is expressed in
  response to. "A very small sample."
* **Answer types differ**: some boolean, some a category, some a ranked list (partners). Each needs a
  probability *and* an uncertainty.
* "The models will keep getting refined the more data comes in, but the actual architecture we should be able to
  set up already."

## The finding that reshapes the problem
**286 of the 411 numeric columns in `nodes.parquet` are measured for 5,000 or more genes.** Most of his questions
are therefore not prediction problems for most genes — the answer is already measured. The questions divide
into four classes, and only one of them is hard:

| Class | Examples from his list | Labels available | What the architecture must do |
|---|---|---|---|
| **A · Measured genome-wide** | nutrient acquisition (`fit_no_glucose`, `fit_no_glutamine`, 7,393 genes); low-nutrient (serum and lipid restriction, 7,395); survival in liver, lung, spleen, brain, heart, peritoneum (`fit_invivo_*`, 7,460 each); dissemination (the cross-tissue pattern of the same set); favourable immune response **in mice** (`fit_ifng` vs `fit_naive_bmdm`, 7,402); cell division (in vitro growth, 7,325) | ~7,400 genes each | **Carry, do not predict.** The work is turning a continuous score into a calibrated answer with an uncertainty — which needs each screen's published FDR, the thing found missing for the Wang screen. Prediction is needed only for the ~700 genes a screen missed |
| **B · Measured for a subset — a well-posed prediction problem** | egress, F-actin, apicoplast, replication (`screen_*`, 318 genes); localisation (716 markers); function class (242); phospho/acetyl/methyl/ubiquitin sites (717) | 242–717 | Standard scoring, and **far better posed than host-facing was with 65** |
| **C · Few or no labels** | a favourable immune response **in humans** (no human screen exists — identified 2026-11-04); autophagy; protein scavenging; sporozoite-to-tachyzoite switching | < 100, sometimes 0 | Curation (the host-facing review is the template) or a new experiment. Say so rather than predicting into a void |
| **D · Not boolean** | most likely binding partners; what a protein is expressed with; what it is expressed in response to | crosslinks 1,630 genes; co-expression 8,052; 266 expression columns | **Ranked lists with a probability per item**, not a single cell |

## Architecture
Five parts. Only the second is per-question work.

1. **A question registry.** One row per question: `id`, the question in Einar's own words, `answer_type`
   (boolean / category / ranked_list / number), `label_source`, the tier of those labels, `evidence_allowed`
   (families), `evidence_forbidden` (`depends_on`, so a question is never answered from its own label source),
   `class` (A–D above) and `status`. This generalises the `Trait` registry in `parascape/information_map.py`,
   which already carries tier, family, normalization, role and depends_on.
2. **A label loader per question** — the only bespoke code. For class A it reads a screen column and a
   threshold; for B a curated table; for D a pair or neighbour table.
3. **One scoring harness, written once.** Orthogroup-grouped folds, the stated baselines, isotonic calibration,
   the Clopper–Pearson lower bound, the fill rule, the no-knowledge split, permutation importance by family.
   This has now been written three times (host-facing, function, partners); factor it out rather than a fourth.
4. **Two numbers per cell, not one.** `probability` = calibrated P(answer), and `uncertainty` = the width of its
   interval (a credible interval for a boolean, the spread of the fold predictions, the standard error for a
   measured score). A confident 0.5 and an ignorant 0.5 are different cells and must look different.
5. **Question generation, not question writing.** 100–1000 questions are not typed by hand; they are derived
   from the label sources: 26 hyperLOPIT niches, 15 function classes, ~45 fitness and condition columns, every
   GO term with at least 10 experimental genes, 266 expression columns as "responds to X", 313 domain families
   as "has domain X" (measured, not predicted). That reaches several hundred immediately, and each arrives with
   its label source already attached.

## The honest caveat that must ship with the matrix
**A thousand questions is not a thousand independent pieces of knowledge.** Many will be variations of a few
underlying axes: the seven in-vivo tissue columns are largely one dissemination signal, and the serum, lipid and
density screens share a growth axis. The matrix must therefore publish a **redundancy report** — questions
clustered by the correlation of their *labels*, stating how many independent axes the set really contains.
Without it the result is an impressive-looking matrix carrying five facts. The same rule as everywhere else
applies: report each question against the hardest honest comparison, and "cannot be answered yet" is a result.

## Order of work
1. The question registry and the refactored scoring harness (the architecture itself).
2. Class A done properly: every genome-wide screen carried as a calibrated answer with uncertainty, which alone
   answers a dozen of his questions for ~90 % of genes. Needs the published per-gene statistics for each screen.
3. Class B scored through the harness: egress and the other imaging phenotypes, localisation niches, function.
4. Class D as ranked lists: partners, co-expression neighbours, expression response.
5. Class C listed as gaps with what each would need, so the missing experiments are visible rather than implied.
