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
| **4 · My own data** | Paste a gene list · Score my pull-down against the background · Compare my screen | compare panel, background frequency specificity, importer |
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
| **Score my pull-down against the background** | The background frequency specificity ranking for the hit list, with frequency shown as context only and the warning that frequency alone must not be used as a filter |
| **Compare my screen** | Joins the list to existing evidence, keeping unresolved identifiers visible |

### 5 · Where it came from — "which measurement, which paper"
Datasets and provenance · References · Downloads and components. Dataset-level browsing only; per-claim
provenance stays on the claim.

## Final structure (Einar, 2026-11-05): the query unit decides level 1
He replaced the five intention categories with three, and the reason is better than the original: **the map is a
matrix, so a query is a row, a block of rows, or a column.** Level 1 is therefore the unit being asked about,
not the user's mood. `Sources` is added as a fourth because provenance browsing belongs nowhere else; per-claim
provenance still travels with the claim.

**Level 1: One gene · Gene list · Label · Sources.** Search and the assistant are present at every level.

### Level 2 — One gene (one row of the matrix)
| Node | Opened |
|---|---|
| Everything known | The information map for that gene: trait, value, probability, uncertainty, source tier, reference; tier 1 then tier 2; measured negative amber, **not measured grey italic** |
| Network | The gene centred, neighbours by evidence kind, edge width by probability, fade outward, click to recentre, scroll to zoom; every neighbour carries this same level-2 menu |
| Partners | Partner, side, probability, tier, reference, and the line that an unlisted partner is unmeasured |
| Expression | Its row over the 266 columns, relative-to-control first, each with normalization, effective n and mean sample quality |
| Structure | AF3 confidence, 3Di coverage and source (measured, stitched, aligned with identity), topology, signal peptide, sequence regions |
| Literature | Mention and publication counts, with the attention caveat stated |

### Level 2 — Gene list (a block of rows)
| Node | Opened |
|---|---|
| Drop or paste | File name, sheet, sniffed separator, identifier column, all changeable; resolution report with unresolved identifiers **listed, not discarded**; extra per-gene columns shown with what each would be treated as. Reader `importer.py`: `.csv .tsv .txt .tab .xlsx .xls .parquet`, and **the range decides, not the filename** |
| Annotate from the map | Every gene with its facts, probabilities and sources; exportable |
| Answer a label for the list | A chosen label answered per gene, each marked filled or exploratory |
| Score as a pull-down | background frequency specificity ranking; frequency as context only, with the warning it must not be a filter |
| Compare as a screen | Joins the list to existing evidence, unresolved identifiers kept visible |

### Level 2 — Label (one column of the matrix)
| Node | Opened |
|---|---|
| The label catalogue | Every label: plain-words question, answer type, labels available, cross-validated **and** frozen-test numbers, genes answerable at 90 %, state (answerable / too few labels / no data). Header: **effective independent axes** |
| One label | Its scoreboard inline, then the ranked genes with probability and uncertainty, filterable filled versus exploratory, with the evidence families used and those forbidden |
| How this label is answered | The strategy behind it: inputs, what is downstream, which families it may and may not see, last scoreboard |
| Compare two labels | Two columns side by side — which genes agree, which disagree; for expression labels this is the condition comparison, with standard errors and study counts |

### Level 2 — Sources
Datasets and provenance · References · Downloads and components. Dataset-level only.

## Node names: one word where possible (Einar, 2026-11-05)
"Preferentially each node should be one or at least fewer than three but preferentially one or two words."
So the **title carries the name and the one-line summary carries the meaning** — the card already has both, and
a title that explains itself is a summary in the wrong place. Final names, all one word:

| Level 1 | Level 2 |
|---|---|
| **Gene** | Facts · Network · Partners · Expression · Structure · Papers |
| **List** | Drop · Annotate · Answer · Pulldown · Screen · Export |
| **Label** | Catalogue · Strategy · Compare · Fill |
| **Sources** | Datasets · References · Downloads |

Level 1 is four single words, which also makes the four parallel: `Gene` is a row, `List` is a block of rows,
`Label` is a column, `Sources` is where they came from. The summaries do the explaining, for example:
* `Gene → Facts` — "every trait recorded for this gene, with its probability and source"
* `Gene → Papers` — "what has been published about it, and how much of that is attention rather than biology"
* `List → Drop` — "a csv, tsv, txt, xls, xlsx or parquet of genes, with any data per gene"
* `List → Answer` — "answer one label for every gene in the list"
* `List → Pulldown` — "rank a hit list against the pull-down background"
* `Label → Catalogue` — "every question, how well it is answered, and how many genes it reaches"
* `Label → Strategy` — "how this label is answered, and what is downstream"
* `Label → Fill` — "write the answers that clear the 90 % bar; nothing else"
A single-word title must never be the only thing a user reads: the summary is mandatory, not decorative, and a
node without one fails its test.

## Corrected node design (Einar, 2026-11-05): the nodes must expose the inference, not browse the map
His objection, and it was right: the previous design described a **browser for the information map** and left
out what ParaScape actually is -- the 39 strategies in 9 mechanism families (`strategy_catalog`), the UMAP and
embedding search, label inference, discoveries and claims. "That is basically the current application, so are
you throwing away the entire application?" **Nothing is thrown away.** Every existing panel is a destination
below, and the three nodes that were missing are the point of the program: infer a label, characterise a list,
explore the space.

**Level 1: Gene · List · Label · Space · Sources.** `Space` is added because the embedding is neither a row, a
block nor a column -- it is the geometry all rows live in, and the UMAP gallery is a first-class feature.

### Gene — one row
| Node | Opened |
|---|---|
| **Card** | The summary page: the AlphaFold structure viewable, the datasets this gene appears in, the papers it is referenced in, and **two tables** — a standardised one with the most relevant data for any gene, and a second showing **how this gene is unusual** (rare or absent domain, lineage-specific, outlier expression, extreme fitness, unique modification, structure unlike anything). Replaces the separate literature node, which is folded in here |
| **Facts** | Every trait recorded, with probability, uncertainty, source tier and reference |
| **Predicted** | What the strategies infer for this gene that nobody has measured: the value, which strategy produced it, its held-out score, and whether it cleared the 90 % bar. Marked exploratory when it did not |
| **Similar** | The genes most like this one, by a chosen kind of evidence — the "borrow from neighbours" family and the embedding neighbours |
| **Network** | The edge view, with the three levels of edge control described below |
| **Partners** | Measured interactions only, with their own probabilities |
| **Expression** | Its row over the expression columns, relative-to-control first |
| **Structure** | The numbers behind the Card's viewer: AF3 confidence, 3Di source and coverage, topology, signal peptide, regions |

### List — a block of rows
Family 5 of the strategy catalogue already is "a list in, a ranking and a profile out"; these nodes expose it.
| Node | Opened |
|---|---|
| **Drop** | Import a csv, tsv, txt, tab, xls, xlsx or parquet of genes with any data per gene; sniffed separator, chosen identifier column, and unresolved identifiers **listed, not discarded** |
| **Profile** | **What characterises this list**: which labels, compartments, domains, expression patterns and fitness signatures are over-represented against a matched background, each with an effect size and a p-value, and the warning when the list is too small to say |
| **Mechanism** | **What the list points to**: the shared pathway, complex or process that best explains it, the strongest common evidence, and the competing explanations it does not distinguish |
| **Rank** | Rank the whole genome by similarity to the list — the lookalike strategy — returning candidate members with probabilities |
| **Answer** | Answer one chosen label for every gene in the list |
| **Pulldown** | Rank the list against the pull-down background: specificity per protein, frequency as context only |
| **Compare** | Two lists against each other: shared genes, and what differs in their profiles |
| **Export** | The list with everything attached, and the provenance of each value |

### Label — one column
This is where the inference lives.
| Node | Opened |
|---|---|
| **Catalogue** | Every label: the question in plain words, labels available, cross-validated and frozen-test numbers, genes answerable at the bar, and its state. Header carries effective independent axes |
| **Infer** | **Run strategies to predict this label** for the genes nobody has measured: pick strategies or let the scorecard pick, see the held-out score, the calibration and the ranked predictions. The machine-learning core of the program |
| **Strategies** | The 39 strategies in their 9 mechanism families, filtered to those that can address this label, each with its self-test, scorecard task and techniques |
| **Scorecard** | Per strategy for this label: cross-validated estimate, frozen-test result, baselines, calibration error, which evidence families it used and which it was forbidden |
| **Fill** | Write the answers that clear the 90 % frozen-test bar, and nothing else; everything below stays exploratory |
| **Compare** | Two labels side by side: where they agree, where they disagree, and whether they are really one axis |

### Space — the geometry of all rows
| Node | Opened |
|---|---|
| **UMAP** | The embedding gallery, with a choice of which evidence block defines the space |
| **Search** | A gene's neighbourhood in that space |
| **Clusters** | The regions, and what each is enriched for |
| **Project** | Drop a list and see where it lands |

### Sources
Datasets · References · Downloads.

## Network edges: three levels of control
| Level | What the user sets |
|---|---|
| **1 · General** | One control, "connect by", with presets: **Interaction** (measured only) · **Expression** (expressed together) · **Behaviour** (behaves alike across perturbations) · **Place** (same compartment) · **Structure** (similar fold) · **Everything** |
| **2 · Layers** | Per-evidence toggles with one slider each: crosslinks, pull-downs, co-expression, co-fitness, shared compartment, shared domain, 3Di similarity, orthology — so a user can ask for "crosslinks plus co-expression, nothing else" |
| **3 · Granular** | Per layer: probability floor, minimum number of supporting studies, correlation cutoff, top-k rank cutoff, whether inferred edges are admitted at all, maximum neighbours per node, and the layout controls (fade distance, zoom, recentre depth) |
Level 1 is the default and the only thing visible until asked for. An edge always states which layer produced it
and its probability, so a network cannot imply more certainty than its weakest edge.

## Final: three level-1 nodes, provenance and the embedding baked in (Einar, 2026-11-05)
He removed `Sources` and asked for it to be baked in. The same argument retires `Space`, which I had added a
moment earlier: a UMAP is a *view of the rows*, not a fourth kind of query, so it belongs inside the three.
**Level 1 is exactly the three shapes of the matrix: Gene (a row) · List (a block of rows) · Label (a column).**

**Where `Sources` went.** Provenance is not a place, it is an attribute of every number:
* every value shows its source tier and reference **inline**, and a value can be opened to its dataset;
* `Gene → Card` carries the datasets the gene appears in and the papers it is referenced in, which is what the
  old `Datasets` and `References` nodes were for;
* downloads happen **on demand**: a node that needs absent data offers to fetch it then, through the existing
  component download dialog, instead of a Downloads page nobody visits before they need it.

### Gene — one protein
| Node | Description shown on the card |
|---|---|
| **Card** | Its structure, the datasets it appears in, the papers that mention it, and two tables: the standard facts for any gene, and what makes this one unusual |
| **Facts** | Every trait ever recorded for it, with how certain each is and who measured it |
| **Predicted** | What the strategies think is true but nobody has measured, and whether that claim is strong enough to keep |
| **Similar** | The genes that most resemble it, by whichever evidence you choose to judge resemblance |
| **Network** | Its neighbours as a map you can walk, with you deciding what counts as a connection |
| **Partners** | The proteins it has been shown to bind, and how good that evidence is |
| **Expression** | Where and when it is switched on, preferring measurements made against a control |
| **Structure** | The numbers behind the fold: confidence, topology, signal peptide, and the regions either side of the membrane |

### List — many proteins
| Node | Description shown on the card |
|---|---|
| **Drop** | Bring your own genes from a spreadsheet or a paste, with any data you have for them |
| **Profile** | What this set has in common that a random set would not: places, domains, conditions, phenotypes |
| **Mechanism** | The process or complex that best explains the set, and the explanations it cannot rule out |
| **Rank** | The rest of the genome ordered by how well each gene would fit this set |
| **Answer** | One question answered for every gene in the set, each with its own certainty |
| **Pulldown** | Which of your hits are specific and which turn up in everybody's pull-downs |
| **Project** | Where the set sits in the map of all genes, and what surrounds it |
| **Compare** | Two sets against each other: what they share, and what separates them |
| **Export** | Everything attached, with the origin of every value |

### Label — one question across all proteins
| Node | Description shown on the card |
|---|---|
| **Catalogue** | Every question the data can be asked, how well each is currently answered, and for how many genes |
| **Infer** | Run the strategies on one question and predict it for the genes nobody has measured |
| **Strategies** | The thirty-nine ways to infer, grouped by mechanism, each stating what it may and may not use |
| **Scorecard** | How a strategy did on answers it was not allowed to see, against the simple baselines |
| **Clusters** | The regions of the map this question lights up, and what else lives there |
| **Fill** | Write only the answers strong enough to keep; everything weaker stays marked exploratory |
| **Compare** | Two questions side by side, including whether they are secretly the same question |

## Scope: an overhaul of the application, with the goal unchanged (Einar, 2026-11-05)
"This is an overhaul of the application. The goal remains the same. It was just not intuitive enough to get
through the application before, that's why we are doing this overhaul; additionally we are adding a lot of new
data and the app needs to accommodate it."

So two drivers, and both are about reaching the same goal rather than changing it:
1. **Navigability.** The capability was there and unreachable. Hence three level-1 nodes shaped like the data,
   one-word names with the meaning in the summary, strategies behind the questions they answer, and provenance
   attached to values instead of parked in a page.
2. **Scale.** The data grew and will keep growing, so the app must accommodate it by design:
   * **Nothing loads wholesale.** The shipped map is already 53 MB and the catalogue will run to hundreds of
     labels and thousands of columns; a panel reads the gene, the label or the list in front of it, never every
     table. `load_published` already returns what exists and names what is missing -- that pattern generalises.
   * **The catalogue is searched, not scrolled.** With hundreds of questions a list is unusable: ranked search
     (`help_index`) is the primary way in, with filters on state, class, label count and answerability.
   * **Columns are found by meaning, not by name.** Tier, family and normalization are registry fields already;
     the UI must group by them, because `expr_lfc_GSE223877_DD_BFD1_Shield_1a` is unreadable and there are 266
     such columns today.
   * **Data arrives on demand.** A node that needs an absent table offers to fetch it, so the install stays
     small as the corpus grows.
   * **Every release is dated and frozen**, so a prediction can be tested against what was published later --
     the prospective check in instruction 73 depends on it.
3. **Nothing is discarded.** Every panel that exists becomes a destination behind a node. The overhaul changes
   how the application is entered, not what it can do -- and the new data is what the new nodes exist to reach.

## Phase 1 built (2026-11-05): the node model and the home screen
Added: `parascape/nodes.py`, `parascape/home_panel.py`, `parascape/gene_drop.py`,
`tests/test_nodes.py`, `tests/test_home_panel.py`, `tests/test_gene_drop.py`.
Changed: `parascape/app.py` (`_home`, installed last so it is the tab that is up),
`parascape/help_index.py` (a `node` kind, first in `KIND_ORDER`, and `node_entries`),
`parascape/help_search.py` (`open_node`), `parascape/importer.py` (the separator sniff).

* **Three level-1 categories** (`Gene`, `List`, `Label`) and **26 level-2 cards**: Gene 9, List 9,
  Label 8, plus **7 level-3 question nodes** under `Label ▸ Catalogue`, one per registered question
  except the retired `has_domain_ipr000719`. Nesting is category → node → view and no deeper, which
  the suite checks. `Label` is drawn in the settled order (Catalogue · Infer · Strategies ·
  Scorecard · Map · Runs · Fill · Compare); `Strategies` is no longer "advanced", because losing
  nothing of the previous ParaScape means the 39 strategies stay first-class. The advanced-last
  mechanism is kept and tested, for whichever node needs it next.
* **One search bar is the navigator** (`nodes.search`): genes by accession or product, questions,
  strategies, datasets, panels and guide pages, ranked and interleaved by `help_index.search` over
  all of them at once, every result drawn as a card in a `Found` section above the home sections.
  Clicking a gene result selects that gene and opens its card, not a generic page; Down from the box
  moves into the results and Return opens one, as the Help field behaves. Two new destination kinds
  carry it: `gene` (select the gene, raise its card) and `entry` (hand the row to
  `help_search.open_entry`, which already knows where each kind lives). A node made by the search
  is marked `from_search`, which is the only exemption from the one-or-two-word title rule -- a
  result is named by the thing it found.
* **The layout is spaCR's.** `theme.RADIUS` gained `tile: 16` and `theme.SPACING` is spaCR's scale
  (4/8/12/16/24/32); the cards follow its `QPushButton#AppTile` rule -- no resting rim, an accent
  rim when focused or hovered, because a rim on every tile made the pane read as a ruled table --
  and an unavailable card keeps a dashed one so it cannot be mistaken for a working tile.
* **Nothing was removed.** `_home` is called after the other installs, so every dock is still built
  and still tabbed; `test_nothing_was_removed_to_make_room_for_the_home_screen` names the thirteen
  and raises each one. The home dock is the tab that is up on startup.
* **`Map` is one node under all three categories** (`nodes.MAP_LENS`), the same `umap_gallery_panel`
  entered with three contexts -- a gene marked, a set projected, a column coloured -- rather than a
  fourth category; `Project` and `Clusters` are gone. The embedding strategies stay under
  `Label ▸ Strategies` and `Label ▸ Infer`.
* **The user's own layer has its hook, not its implementation.** `Node.produces_inference` marks the
  four nodes that show an inference (`Gene ▸ Predicted`, `List ▸ Answer`, `Label ▸ Infer`,
  `Label ▸ Fill`), which is where **Keep** will attach once the harness produces claims to keep. Two
  structural guarantees already hold: a node that shows an inference must also show its probability,
  uncertainty and source tier, and `Gene ▸ Facts` **raises** if it is ever marked as showing an
  inference -- measured and inferred stay separate surfaces.
* **Unavailable is a rendered state with a reason**, which is how the missing work is now visible on
  the screen rather than discovered by clicking: `Gene ▸ Similar` (needs the phase-3 neighbour
  space), `List ▸ Profile`, `List ▸ Mechanism`, `List ▸ Rank` (no set-level view exists),
  `List ▸ Pulldown` (no background table), `Label ▸ Runs` (the harness writes
  `results/harness_holdout_log.jsonl`, but no panel reads it back),
  `List ▸ Answer/Map/Export` (until a list is
  dropped), `Label ▸ Fill` (no question has cleared the 90 % frozen-test bound), and the two
  declared gaps, which render as "cannot be answered yet" rather than as a page.
* **Trust is inline and structural.** `Node.__post_init__` refuses a claim-bearing node that does
  not declare probability, uncertainty and source tier, so the rule cannot be dropped by whoever
  draws the next card. A question node carries its scoreboard -- cross-validated beside frozen test
  -- its coverage sentence ("7,393 of 8,140 genes ... the rest are not measured, which is not a
  negative") and its state; the catalogue header carries effective independent axes, reported as an
  **upper bound** because the label-correlation clustering of instruction 73 has not been computed.
* **A dropped file is a gene list.** `gene_drop` wraps `importer`, adds the resolution report in
  which an unresolved identifier is listed rather than dropped (and an ambiguous one is never
  guessed), and refuses an unknown extension with `importer.READABLE`. The window accepts a drop
  anywhere through an event filter and routes it to `List ▸ Drop`.
* **A defect the drop path found, and one it could not fix.** `importer.read_any` sniffed its
  separator with `csv.Sniffer` and no candidate list, so on a one-column list of gene identifiers it
  chose a LETTER: `gene\nTGME49_208830` came back as three columns split on "e". The sniff is now
  restricted to `,`, tab, `;` and `|`, and a file with none of them is read as one column. Separate
  and NOT fixed here: `app.apply_import` calls `identity.GeneIndex.load().resolve`, and
  `GeneIndex` has neither `load` nor `resolve` -- the call is inside a `try/except` that silently
  falls back to `resolve=None`, so no import has ever resolved a previous or strain accession.
* **Still phase 1 only.** A question node opens its row in the catalogue rather than a per-question
  panel (phase 2), `Gene ▸ Network` opens the star map as it is, without recentring or the fade
  (phase 3), and the assistant is not yet reachable from a node (phase 4).

## Two additions (Einar, 2026-11-05): a user's own layer, and the embedding as a lens
### The embedding was undersold as `Clusters`
He asked whether the UMAP strategies were represented by `Label → Clusters`, and doubted it. Correctly: the
embedding is strategy **family 1 of nine** ("search the map space -- maps, walks and held-out labels", 01-06),
plus the UMAP gallery panel, and "the regions this question lights up" is one use of it.

**Fix: `Map` is one node name appearing under all three level-1 items, entering the same view with different
context** -- a lens, not a fourth category:
* `Gene → Map` — the space with this gene marked and its neighbourhood around it;
* `List → Map` — the space with the set projected onto it, and what surrounds it (replaces `Project`);
* `Label → Map` — the space coloured by the answer, showing where it concentrates (replaces `Clusters`).
The embedding **strategies** stay where they belong, under `Label → Strategies` and `Label → Infer`: families
1 and 2 are ways to infer, not ways to browse, and the scorecard judges them like any other.

### A user's own layer, so an inference can be kept
"There should be mechanisms that connect the one protein and labels, where if inference is done the user can
choose to add it to their own ParaScape dataset and then it should turn up in `Predicted`."

* **Accepting.** Any node that produces an inference — `Label → Infer`, `List → Answer`, `Gene → Predicted` —
  offers **Keep**. Keeping writes the claim to the user's own layer, never to the shipped map.
* **What is recorded with it**: the value, its probability and uncertainty, which strategy and version produced
  it, the question, the score it was accepted at, whether it cleared the 90 % frozen-test bar, and the date.
  A kept claim is reproducible from its record or it is not kept.
* **Where it appears.** `Gene → Predicted` shows kept claims beside the unaccepted ones, visibly distinct from
  both measured facts and from predictions nobody has accepted. `Facts` **never** shows a kept inference:
  measured and inferred stay separate surfaces, which is the whole point of the probability and source columns.
* **`Label → Fill` is the bulk form of the same action**: write the answers that clear the bar to the user's
  layer, with the option to keep weaker ones explicitly marked exploratory. Nothing is ever written to the
  shipped map by the application.
* **Rules.** A kept inference is **never training data** — the no-self-training rule from instruction 72, now
  enforceable because the layer is separate. It can be reverted, individually or by strategy version. It is
  exportable with its provenance, so it can be shared or cited. A released map upgrade that *measures* a kept
  claim supersedes it and says so, rather than silently agreeing with it.

## Label, designed against Einar's four requirements (2026-11-05)
He accepts `Gene` (exploratory) and `List` (compare your data to ParaScape), and names `Label` as the hard part:
"this is where inference happens and where the data is actually used to generate new knowledge." His four
requirements, and what each forces:

| Requirement | What it forces in the Label nodes |
|---|---|
| Lose nothing from the previous ParaScape | all 39 strategies in their 9 mechanism families remain, reachable from `Strategies`, each keeping its own self-test and scorecard task; the UMAP gallery, discoveries, claims and strategy panels are the destinations |
| All the new strategies implemented, accessible, intuitive, **controlled** | the new ones (positive-unlabelled weighting, propensity weighting, matrix completion, per-family stacking, lookalike ranking, selective prediction, network propagation, multi-task, shrinkage, condition fingerprints) enter the **same** catalogue with the **same** quality control: nested cross-validation, a frozen partition, baselines, calibration. An untested strategy is listed as untested, never offered as if it were equal |
| Easy to add new data as it arrives | a column cannot enter without an information category and a family; the picker groups by those, so new data appears in the right place without a UI change |
| Easy to save inference runs with full metadata | `Runs` is a node, and a run is reproducible or it is not saved |

### Label — level 2
| Node | Description on the card |
|---|---|
| **Catalogue** | Every question the data can be asked, how well each is answered, and for how many genes |
| **Infer** | Choose a question, choose which data may be used, choose a strategy, and predict it for the genes nobody has measured |
| **Strategies** | Every way to infer, grouped by mechanism, each stating what it may use and how it was tested |
| **Scorecard** | How a strategy did on answers it was not allowed to see, against the simple baselines |
| **Map** | The space of all genes coloured by this question, and where the answer concentrates |
| **Runs** | Every inference you have run, what data it used, and the result — re-runnable |
| **Fill** | Keep only the answers strong enough; everything weaker stays marked exploratory |
| **Compare** | Two questions side by side, including whether they are secretly the same question |

### `Infer` is a three-step builder, because that is the intuitive form of "some or all data"
1. **Question** — pick it, or write a new one from a label source; the state says whether it can be answered.
2. **Data** — pick *some or all* by **information category**, family and tier. The panel states what is in,
   what was **excluded automatically to prevent circularity** (with the categories named and counted), and what
   a force-include would cost. This step is what makes a UMAP-plus-clustering run honest rather than circular.
3. **Strategy** — pick one, several, or let the scorecard choose; each shows what it needs and whether it has
   ever been tested on this kind of question.
Then it runs, and returns the cross-validated and frozen-test numbers together, the baselines, the calibration,
and a **Keep** action that writes the accepted claims to the user's own layer.
