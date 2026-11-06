# 91 · Navigator: one click shows what comes next, a second click opens

**Status: built 2026-11-06, waiting on Einar.** One click on a module shows it with what it leads to
(`modules.DOWNSTREAM`, at most four per step, edges drawn, the clicked path highlighted); a click on a next
step shows what follows it; a second click on a highlighted module opens it; Back / Escape return to the
results. "map" now gives 5 cards (Map, 7 genes, a question, 4 methods, 7 settings and panels) instead of 24;
a group opens to its members. The Navigator is emptied every time Home is shown. Tests in
`tests/test_navigator.py`.

Requested 2026-11-06. Einar: "I want the network if possible like in spaCR of what is
downstream of map. And a second click to open the module. Also there are way too many options when I type
map. I'm sure some of them can be consolidated for navigator ... navigator should be reset every time I go
back to it."

## What spaCR does (`spacr/qt/widgets/home.py`, `_EasyWorkflowPage` / `_EasyWorkflowView`)
spaCR's own hint: "Explore: Click once to open the module network. Click twice to open the module." One click
on a module card centres it and grows its direct downstream modules beside it (up to four, the rest on
request), joined by edges; clicking a downstream node expands it in turn; the node clicked twice opens.

## The plan
1. **First click explores, second click opens**, for module results. Genes, questions and datasets keep one
   click (they have no downstream).
2. **A downstream graph for ParaScape's modules**: Map → Clusters, Network, Gallery, Export; Gene → Partners,
   Expression, Structure, Network; Drop → List, Screen, Infer; Infer → Scorecard, Fill, Runs ... declared in
   `parascape/modules.py` beside each module, so it is one list to keep right.
3. **Consolidate the results**: one card per destination (a module and a help entry that opens the same
   panel are one card); help/guide entries folded under their module; at most one row of cards.
4. **Reset on return**: every time Home is shown the field is cleared, cards and graph removed, the hint back.

## Gate
Typing "map" shows a handful of cards; one click on Map shows its downstream modules; a second click opens it;
going Home shows an empty Navigator.
