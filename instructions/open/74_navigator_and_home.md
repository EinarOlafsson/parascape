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
| **5 · How much can this be trusted?** | Scoreboards (cross-validated and frozen test) · Coverage · Effective independent axes · What cannot be answered yet | harness results, scoreboard, track record |
| **6 · Where it came from** | Datasets and provenance · References · Downloads | dataset coverage, deposited values, component downloads |

Category 5 on the home screen is deliberate and unusual: this project's value rests on knowing which cells are
measured, which are inferred and which are unanswerable, so the honesty surface is a destination rather than a
footnote.

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
