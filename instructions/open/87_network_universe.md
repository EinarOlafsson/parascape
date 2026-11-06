# 87 · The network as a navigable universe

**Status: open (instruction 74, goal 8; requested 2026-11-05).** Einar: "a more universe feel where the
current nodes are in focus on a map where the nodes farther away from the center of the screen become more
obscure, but you can click any of the visible nodes and click and pull your way to them. There should be
some options to the right that determine the edges and the organization of the network."

**Status: built 2026-11-06 (first version), waiting on Einar.** `parascape/network_module.py` +
`parascape/network_reading.py`, through the parascape-module skill: every gene drawn (as on the map, or
arranged by the links -- spectral embedding of the largest connected part, then UMAP, in the background,
kept per choice of links); the kinds of link ticked (13 measured layers, plain names and colours) above a
chosen strength percentile; the focus gene sharp (its links bright, its partners faint, the strongest named)
and the rest fading by link steps; wheel zooms about the pointer, drag pans, click travels to a gene; the
reading says what it is linked to by kind and how often a gene's links name its label against chance. Not
yet: links inferred by strategies as a layer, true 3D depth of field. Tests in
`tests/test_network_reading.py`, `tests/test_module_shell.py`.

**Restated 2026-11-06.** Einar: "the network still looks nothing like I want. Now you just get a gene with its
connected genes and no context. I want to see the entire network with the gene's region zoomed in, with the
ability to zoom the network and change the edge definitions of the network."

## The plan
1. Depth of field: nodes sharp near the centre, progressively faded and blurred with distance from it.
2. Click any visible node to travel to it; click-and-drag to pull through the space.
3. A right-hand panel: what defines an edge (co-mention, co-fitness, interaction, co-expression ...), the
   threshold, and the layout.
4. No spaCR equivalent exists, so this is built rather than ported, in spaCR's visual language.
