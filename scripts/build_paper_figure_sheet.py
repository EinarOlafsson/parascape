"""The twelve figures of the ParaScape paper, real panels where they exist and placeholders where
they do not.

Einar, 2026-11-19: "we will make the paper figure with real data and with placeholders. then we
will fill in the data for the paper. then we will publish the paper and parascape needs to allign
with the paper analasys ... this is a developing target. but we do need to develope the target!"

So this is the target, written down. Every panel the twelve-figure plan calls for is listed here
once. A panel whose drawing exists is copied in; a panel whose measurement has not been made is
drawn as a marked placeholder at the SAME size, so the page does not reflow when the real one
arrives and nobody mistakes an empty panel for a finding.

It is deliberately a declaration and not a search. A script that globbed for whatever panels happen
to exist could never report a MISSING one, and the missing ones are the point: they are the list of
what the paper still needs, and `summary.json` is that list in a form the Findings module and the
progress report can both read.

Writes `results/paper_figure_sheet_<date>/`: every panel as `figNN<letter>.png`, plus a manifest.
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import sys
from dataclasses import dataclass
from datetime import date

import matplotlib
matplotlib.use("Agg")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import figstyle as FS      # noqa: E402


@dataclass(frozen=True)
class Panel:
    """One panel of one figure, and where its drawing comes from.

    `source` is a glob relative to `results/`, newest match wins, or None when nothing draws it
    yet. `why` says what is missing, and is printed on the placeholder.
    """
    figure: str
    letter: str
    title: str
    source: str | None
    why: str = ""


#: The twelve figures. Order is the order of the argument, not the order they were measured in.
PANELS = (
    Panel("01", "a", "What ParaScape holds", "figure1_schematic_*/panel_schematic.png"),
    Panel("01", "b", "What the corpus reaches", "corpus_inventory_*/panel_coverage.png"),
    Panel("02", "a", "Which studies can be compared",
          "paper_panels_*/figure02_comparability.png"),
    Panel("03", "a", "The co-expression modules", "paper_panels_*/figure03_modules.png"),
    Panel("04", "a", "What structures a contrast", "paper_panels_*/figure04_what_structures.png"),
    Panel("05", "a", "Which subtype a treatment reaches",
          "bradyzoite_reach_*/panel_reach_by_subtype.png"),
    Panel("05", "b", "The evidence behind those counts",
          "paper_panels_*/figure05b_reach_violin.png"),
    Panel("06", "a", "Which evidence places a protein",
          "coexpression_components_localization_*/panel_which_evidence_places_a_protein.png"),
    Panel("06", "b", "What sets protein abundance",
          "paper_panels_*/figure06b_protein_abundance.png"),
    Panel("06", "c", "Which data explains which",
          "paper_panels_*/figure06c_explains_which.png"),
    Panel("07", "a", "The whole network against one layer",
          "paper_panels_*/figure07_network_layers.png"),
    Panel("08", "a", "What differs between strains",
          "paper_panels_*/figure08_strain_variation.png"),
    Panel("08", "b", "Host and parasite proteins co-folded", None,
          "AF3 co-folding: Einar asked for it in this figure; not run"),
    Panel("09", "a", "What transfers between species",
          "paper_panels_*/figure09_species_transfer.png"),
    Panel("10", "a", "The host response per effector", "paper_panels_*/figure10_host_shared.png"),
    Panel("10", "b", "The parasite in the same cells",
          "paper_panels_*/figure10b_single_cell.png"),
    Panel("11", "a", "What the map can fill", "paper_panels_*/figure11_filling_the_map.png"),
    Panel("12", "a", "What nothing measures",
          "paper_panels_*/figure12a_measured_absence.png"),
    Panel("12", "b", "How much is known per gene",
          "paper_panels_*/figure12b_categories_per_gene.png"),
)


def newest(pattern: str) -> str | None:
    """The newest drawing matching a glob, or None when nothing does."""
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    return found[-1] if found else None


def main() -> None:
    out = os.path.join(ROOT, "results", f"paper_figure_sheet_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)

    rows, real, placed = [], 0, 0
    for panel in PANELS:
        name = f"fig{panel.figure}{panel.letter}"
        found = newest(panel.source) if panel.source else None
        if found:
            shutil.copyfile(found, os.path.join(out, f"{name}.png"))
            source = os.path.relpath(found, ROOT)
            pdf = found[:-4] + ".pdf"
            if os.path.exists(pdf):
                shutil.copyfile(pdf, os.path.join(out, f"{name}.pdf"))
            real += 1
            print(f"  {name}  {panel.title[:44]:44s} <- {source}")
        else:
            why = panel.why or "not measured"
            figure = FS.placeholder(panel.title, why, letter_text=panel.letter)
            for extension in ("png", "pdf"):
                figure.savefig(os.path.join(out, f"{name}.{extension}"),
                               dpi=600 if extension == "png" else None,
                               bbox_inches="tight", pad_inches=0.02)
            matplotlib.pyplot.close(figure)
            source = None
            placed += 1
            print(f"  {name}  {panel.title[:44]:44s} .. PLACEHOLDER ({why})")
        rows.append({"figure": panel.figure, "letter": panel.letter, "title": panel.title,
                     "real": bool(found), "source": source,
                     "why_missing": None if found else (panel.why or "not measured")})

    figures = sorted({p.figure for p in PANELS})
    whole = [f for f in figures
             if all(r["real"] for r in rows if r["figure"] == f)]
    print(f"\n{real} of {len(PANELS)} panels drawn from data; {placed} are placeholders.")
    print(f"{len(whole)} of {len(figures)} figures are complete: {', '.join(whole)}")
    print("\nStill to measure or draw:")
    for row in rows:
        if not row["real"]:
            print(f"  figure {row['figure']}{row['letter']}  {row['title']}"
                  f"\n      {row['why_missing']}")

    json.dump({"panels": len(PANELS), "real": real, "placeholders": placed,
               "figures": len(figures), "complete_figures": whole, "rows": rows},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
