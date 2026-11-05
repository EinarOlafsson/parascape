# 74 · A navigator and a home screen: organise ParaScape by what the user wants, not by what it contains

**Status: open (requested 2026-11-05).** Implements the navigator half of instruction 70 and supersedes its
sketch. Related: 73 (the question matrix the navigator must expose).

## The request (Einar, 2026-11-05)
* Bring the spaCR navigator into ParaScape with ParaScape logic. The same AI as spaCR, with a ParaScape prompt,
  so the user can query any question or anything else in ParaScape.
* Search for strategies, questions, genes. A result is a **clickable node**. A strategy node opens the strategy
  and shows **what is downstream** of it, as in spaCR. A gene search opens the gene **at the centre of a
  network**, neighbours clickable, the network fading to the sides and upwards, zoomable; clicking a neighbour
  **shifts the view** to it. Clicking a gene offers nodes for learning about it, one of which is the information
  map, which opens the network full screen.
* "The current layout for ParaScape is confusing." He wants a home screen where each path is a node, nodes
  nested in categories, and asks what structure would guide a user intuitively.

## What already exists (so this is navigation, not a rewrite)
* `parascape/chat.py` — the assistant, already the spaCR console's shape: it shells out to whichever vendor CLI
  the user is logged into (`claude`, `codex`, `gemini`), so there is no API key, and it already builds a briefing
  from the live state (selected gene, evidence, colouring, filters, standing caveats). **It needs a ParaScape
  prompt that knows the question registry, not a new assistant.**
* `parascape/help_search.py` — the search box beside Help, built deliberately like spaCR's. **It needs to return
  nodes rather than jump straight to a destination.**
* `parascape/star_map.py` (1,641 lines) — the existing network view. **It needs gene-centred recentring, zoom
  and the fade, not a new graph widget.**
* `parascape/guided_panel.py` — a question-at-a-time "start here" flow. **It becomes one node on the home
  screen, not the entry point.**
* `parascape/information_map_panel.py` — shipped 2026-11-05; the full-screen destination he describes.

## The home screen: organise by intention, and keep strategies behind questions
**Recommendation, and it disagrees with the suggestion of making "inference strategies" the first category.**
Nobody arrives wanting a random forest. They arrive with a gene, a question, or a list of hits, and the reason
the current layout loses people is that it presents machinery before intent -- the same complaint behind
instruction 70 ("as it is now, I get lost right away"). Strategies stay reachable, but behind the question they
answer, plus one explicit "all strategies" node for whoever does want the machinery.

Six categories, each a row of nodes:

| Category | Nodes | Existing panel behind it |
|---|---|---|
| **1 · One gene** | What is known about this gene · Its network · Its partners · Its expression · Its structure · What has been written about it | information map panel, star map, pairs, RNA panels, AF3/3Di, mentions |
| **2 · One question, all genes** | The question catalogue (host-facing, needed without glucose, involved in egress, located where, ...) · Fill a gap for my gene list · *All strategies (advanced)* | questions registry + harness (73), strategy panel |
| **3 · Compare conditions** | Two conditions · A perturbation series · Strains · Stages · Genes that move together | compare, coexpression, RNA design panels |
| **4 · My own data** | Paste a gene list · Score my pull-down against the background · Compare my screen | compare panel, CRAPome specificity, importer |
| **5 · Where it came from** | Datasets and provenance · References · Downloads | dataset coverage, deposited values, component downloads |

### Trust is integrated, not a category (Einar, 2026-11-05)
A sixth category, "How much can this be trusted?", was proposed and **rejected**: if someone has to visit a
separate room to learn whether a number is reliable, the default experience is untrustworthy. His correction is
the better design, and it becomes a rule the panels must follow:
* a node that shows a claim shows its **probability, its uncertainty and its source tier inline**;
* a **question node** carries its own scoreboard -- the cross-validated estimate and the frozen-test result side
  by side -- and how many genes it can answer at the 90 % bar, rather than linking elsewhere;
* **coverage** belongs on the question or trait itself ("716 of 8,140 genes have an observed fact; the rest are
  not measured, which is not a negative");
* **"cannot be answered yet"** is a state on the question catalogue, not a page: too few labels renders as
  unanswerable with the reason;
* **effective independent axes** sits in the catalogue header, qualifying the whole list;
* provenance and references stay reachable from each claim; category 5 holds only dataset-level browsing.

## Node behaviour
* **Every search result is a node**, with a kind (gene, question, strategy, dataset, panel, guide page) and a
  one-line statement of what clicking it does. `help_index` already ranks these; it gains a node kind.
* **A strategy node** shows the strategy, its inputs, and **what is downstream**: which questions it answers,
  which evidence families it may see, which it is forbidden (73), and its last scoreboard.
* **A gene node** opens the gene-centred network: the gene at the centre, neighbours by evidence kind
  (crosslink, co-expression, co-fitness, shared compartment), edge thickness by probability, **click to
  recentre**, zoom, and a fade toward the edges so depth reads as distance. Each neighbour carries the same node
  menu, one entry of which opens the information map full screen for that gene.
* **The assistant is reachable from every node** with that node's context in its briefing, so "what does this
  mean" has a referent. The ParaScape prompt must state what the map does and does not contain -- that nothing
  is filled unless its frozen-test bound cleared 90 %, and that an empty cell is not a negative.

## Phases
1. **Home screen with categorised nodes** and the node model; existing panels become destinations. Nothing is
   removed, so no workflow breaks.
2. **Search returns nodes**; strategy nodes gain the downstream view.
3. **Gene-centred network**: recentring, zoom, fade, node menu, full-screen information map.
4. **ParaScape assistant prompt** aware of the question registry, reachable from every node.
5. Retire the docks that the home screen replaces, only once each has a node.

## The node catalogue: every node, and what it shows when opened
Nesting is three deep and no deeper: **category → node → destination view**. A category is a section with a
title, one line of purpose, and a row of cards. A card shows its title and a one-line statement of what clicking
it does. Strategies nest one level further, under the question they answer. Every claim shown anywhere carries
its probability, uncertainty and source tier inline, per the rule above, and the assistant is reachable from
every node with that node's context in its briefing.

### 1 · One gene — "everything known about one protein"
Entry: the identifier box, a search result, or a gene picked from a dropped list.
| Node | When opened |
|---|---|
| **What is known about this gene** | The information map full screen: one row per trait, tier 1 conclusions first then tier 2 measurements, each with value, probability, uncertainty, source tier and the reference that supports it. A measured negative is amber; **not measured is grey italic** and says the cell is empty, not "no" |
| **Its network** | The gene at the centre; neighbours coloured by evidence kind (crosslink, co-expression, co-fitness, shared compartment), edge width by probability, fading toward the edges; click a neighbour to recentre, scroll to zoom; each neighbour carries this same node menu |
| **Its partners** | The pair table for this gene: partner, side (parasite or host), probability, source tier, reference, and the standing line that a partner not listed has not been measured |
| **Its expression** | Its row across the 266 expression columns, relative-to-control first, then internal reference, then absolute, each with its normalization, effective sample count and the mean sample quality behind it |
| **Its structure** | AF3 confidence, 3Di coverage and which source it came from (measured, stitched or aligned, with identity and coverage), transmembrane topology, signal peptide, and the sequence regions |
| **What has been written about it** | Publication and mention counts with the attention caveat stated, because study attention is a predictor of nothing biological |

### 2 · One question, all genes — "ask one thing of every protein"
| Node | When opened |
|---|---|
| **The question catalogue** | Every question as a row: the question in plain words, answer type, class, labels available, the cross-validated estimate **and** the frozen-test result, genes answerable at the 90 % bar, and a state: answerable, not enough labels, or no data. The header carries **effective independent axes**, so the list cannot read as more knowledge than it holds |
| **One question** (opened from the catalogue) | Its scoreboard inline, then the ranked gene list with probability and uncertainty per gene, filterable to filled versus exploratory, with the evidence families it used and those it was forbidden to see |
| **Answer a question for my gene list** | Takes a dropped or pasted list and returns that question's answer per gene, each marked filled or exploratory |
| **All strategies** *(advanced, rendered last)* | One card per strategy; opening one shows what it is, its inputs, **what is downstream** — which questions it answers, which evidence families it may and may not see — and its last scoreboard |

### 3 · Compare conditions — "what changes between two situations"
Two conditions · A perturbation series · Strains (GT1, ME49, VEG) · Stages (tachyzoite, bradyzoite, oocyst,
merozoite) · Genes that move together. Each opens the existing comparison or co-expression view with the
contrast preselected, the fold change shown with its standard error, and the number of studies supporting it.

### 4 · My own data — "bring a list and get it annotated"
| Node | When opened |
|---|---|
| **Drop or paste a list** (first card; the window also accepts a drag-and-drop anywhere) | The file name, sheet, detected separator and chosen identifier column, all changeable; an identifier resolution report with the unresolved ones **listed, never silently dropped**; any extra per-gene columns with what each would be treated as, left for the user to confirm. Reader: `parascape/importer.py`, which accepts `.csv .tsv .txt .tab .xlsx .xls .parquet`, sniffs the separator rather than trusting the extension, and obeys its own rule that **the range decides, not the filename** |
| **Annotate my list from the map** | Each gene with its map facts, probabilities and sources; exportable |
| **Score my pull-down against the background** | The CRAPome specificity ranking for the hit list, with frequency shown as context only and the warning that frequency alone must not be used as a filter |
| **Compare my screen** | Joins the list to existing evidence, keeping unresolved identifiers visible |

### 5 · Where it came from — "which measurement, which paper"
Datasets and provenance · References · Downloads and components. Dataset-level browsing only; per-claim
provenance stays on the claim.
