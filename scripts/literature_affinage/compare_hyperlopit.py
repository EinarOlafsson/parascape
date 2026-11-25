"""Agreement of the literature compartment with hyperLOPIT, for genes that both place.

A secreted destination counts as its organelle of origin (PV lumen and PVM from dense granules, host
cytosol and host nucleus from rhoptries or dense granules), and the numbered hyperLOPIT sub-clusters
("rhoptries 1", "ER 2") as their organelle. Writes `tables/compartment_vs_hyperlopit.tsv` in the work
folder and prints the rate.

    python compare_hyperlopit.py WORK_FOLDER NODES_PARQUET
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

ORIGIN = {"PV lumen": {"dense granules"}, "PVM": {"dense granules", "rhoptries"},
          "host cytosol": {"dense granules", "rhoptries"}, "host nucleus": {"dense granules", "rhoptries"},
          "cyst wall": {"dense granules"}, "conoid": {"apical"}, "apical annuli": {"IMC", "apical"},
          "basal complex": {"IMC"}}


def coarse(label: str) -> str:
    """'rhoptries 1' -> 'rhoptries'."""
    return re.sub(r"\s+\d+$", "", str(label)).strip()


def agrees(lit: str, lopit: str) -> bool:
    """Whether a literature compartment and a hyperLOPIT one name the same place."""
    lopit = coarse(lopit)
    return lit == lopit or lopit in ORIGIN.get(lit, set())


#: The organelle each fraction belongs to. Imaging rarely tells chromatin from the rest of the nucleus,
#: or the matrix from the membranes, so the organelle level is the fair comparison for a paper.
ORGANELLE = {"nucleus - chromatin": "nucleus", "nucleus - non-chromatin": "nucleus", "nucleolus": "nucleus",
             "mitochondrion - membranes": "mitochondrion", "mitochondrion - soluble": "mitochondrion",
             "PM - integral": "plasma membrane", "PM - peripheral": "plasma membrane",
             "40S ribosome": "cytosol", "60S ribosome": "cytosol", "19S proteasome": "cytosol",
             "20S proteasome": "cytosol"}


def same_organelle(lit: str, lopit: str) -> bool:
    """Agreement at the organelle level: fractions of one organelle count as one place."""
    if agrees(lit, lopit):
        return True
    a, b = ORGANELLE.get(lit, lit), ORGANELLE.get(coarse(lopit), coarse(lopit))
    return a == b


def main(argv=None) -> int:
    """Compare and write the table."""
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        raise SystemExit(__doc__)
    work, nodes = Path(args[0]), pd.read_parquet(args[1], columns=["gene_id", "compartment"])
    labels = pd.read_parquet(work / "tables" / "literature_gene_labels.parquet")
    t = labels[labels["lit_compartment"].astype(str).ne("")][["gene_id", "lit_compartment"]].merge(
        nodes.rename(columns={"compartment": "hyperlopit"}), on="gene_id", how="left")
    t["hyperlopit"] = t["hyperlopit"].fillna("unassigned").astype(str)
    both = t[t["hyperlopit"].ne("unassigned")].copy()
    both["agree"] = [agrees(a, b) for a, b in zip(both["lit_compartment"], both["hyperlopit"])]
    both["same_organelle"] = [same_organelle(a, b) for a, b in zip(both["lit_compartment"], both["hyperlopit"])]
    t = t.merge(both[["gene_id", "agree", "same_organelle"]], on="gene_id", how="left")
    t.to_csv(work / "tables" / "compartment_vs_hyperlopit.tsv", sep="\t", index=False)
    off = both[~both["same_organelle"]]
    pairs = (off.assign(h=off["hyperlopit"].map(coarse)).groupby(["lit_compartment", "h"]).size()
             .sort_values(ascending=False).head(8))
    print(json.dumps({"placed_by_papers": len(t), "both_place": len(both),
                      "agree": int(both["agree"].sum()),
                      "agree_percent": round(100 * both["agree"].mean(), 1) if len(both) else None,
                      "same_organelle": int(both["same_organelle"].sum()),
                      "same_organelle_percent": round(100 * both["same_organelle"].mean(), 1) if len(both) else None,
                      "papers_only": int((t["hyperlopit"] == "unassigned").sum()),
                      "commonest_disagreements": {f"{a} | {b}": int(n) for (a, b), n in pairs.items()}},
                     indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
