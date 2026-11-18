"""Which host gene each Toxoplasma effector moves, per host gene.

Instruction 101, item A2. `hosttx_effectors` is the dual perturb-seq study (Butterworth and Sanchez
et al., Cell Host Microbe 2023; PMID 37827122; GSE229505): pooled CRISPR knockout of parasite
effectors read out by single-cell RNA-seq of both the parasite and the infected host cell.

WHAT WAS ALREADY LOADED, AND WHY IT IS NOT ENOUGH. The release carries the per-effector SUMMARY --
a Hotelling T2, and the 737,726-row host differential-expression table reduced to 20 principal
components with its norm and a count of substantially changed genes. That answers which effectors
change the host and by how much. It CANNOT answer which host gene any effector moves, because a
principal component is a direction through the host transcriptome and not a measurement of a gene.
Figure 10 needs the gene.

So this reads the differential-expression table itself and resolves both ends: the parasite end
from an effector name to a ToxoDB ME49 accession, the host end from a symbol to the human gene
identifier the host space uses. Both resolutions are reported as coverage rather than assumed, and
an unresolved row is counted and kept out rather than dropped silently.

Writes `results/hosttx_per_host_gene_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
from datetime import date

import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import perturbation_atlas as PA   # noqa: E402

SOURCE = os.path.join("/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/DNA",
                      "CRISPR_screen", "37827122", "hosttx_effectors_DE_host_genes.csv")

#: Corrected probability below which an effect on a host gene is kept as an effect.
MAX_Q = 0.05

#: Effector names the product annotation does not carry, and what to look for instead. Everything
#: else resolves by finding its own name as a word in the `product` column of the node table.
#: Resolving from the table rather than from a written-out list of accessions is not fastidiousness:
#: the first version of this file carried 21 accessions typed from memory and SIX were wrong,
#: including GRA70 pointing at TGME49_222100, which is SOS1. A wrong accession here attributes one
#: effector's host effects to a different gene and nothing downstream could detect it.
ALIASES = {
    "HCE1": "E2F4-associated EZH2-inducing gene regulator",   # HCE1 is TEEGR
    "IST": "inhibitor of STAT1 transcription",
    "222100": "effector protein SOS1",                        # the effector this study discovered
}


def resolve_effectors(targets, cache: str) -> tuple:
    """Effector name -> accession, found in the node table, and the names that did not resolve.

    An unresolved effector is REPORTED and its rows are kept out. It is never guessed.
    """
    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"), columns=["gene_id", "product"])
    nodes = nodes.dropna(subset=["product"])
    product = nodes["product"].astype(str)
    resolved, ambiguous = {}, {}
    for target in targets:
        name = str(target)
        if name in ALIASES:
            hits = nodes[product.str.contains(ALIASES[name], case=False, regex=False, na=False)]
        else:
            hits = nodes[product.str.upper().str.contains(rf"\b{re.escape(name.upper())}\b",
                                                          regex=True, na=False)]
        if len(hits) == 1:
            resolved[name] = str(hits.iloc[0]["gene_id"])
        elif len(hits) > 1:
            ambiguous[name] = [str(g) for g in hits["gene_id"].head(4)]
    unresolved = sorted(set(map(str, targets)) - set(resolved))
    return resolved, unresolved, ambiguous


def resolve_hosts(symbols) -> tuple:
    """Host symbol -> human gene identifier, and how many symbols did not resolve.

    `host_symbol_map` drops a symbol that names more than one gene rather than guessing which, so
    a symbol missing here is either unknown to the human gene space or ambiguous within it.
    """
    table = PA.host_symbol_map()
    human = table[table["species"] == "human"].drop_duplicates("symbol")
    mapping = dict(zip(human["symbol"], human["gene_id"]))
    upper = {s: str(s).upper() for s in symbols}
    resolved = {s: mapping[u] for s, u in upper.items() if u in mapping}
    return resolved, sorted(set(symbols) - set(resolved))


def main() -> None:
    if not os.path.exists(SOURCE):
        raise SystemExit(f"source not found: {SOURCE}")
    raw = pd.read_csv(SOURCE)
    print(f"{len(raw):,} rows, {raw.target.nunique()} effectors, {raw.gene.nunique():,} host genes")

    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")
    effectors, unresolved_effectors, ambiguous = resolve_effectors(raw["target"].unique(), cache)
    if ambiguous:
        print("ambiguous, left out:", {k: v for k, v in ambiguous.items()})
    hosts, unresolved_hosts = resolve_hosts(raw["gene"].unique())
    print(f"effectors resolved: {len(effectors)} of {raw.target.nunique()}"
          + (f"  (unresolved: {', '.join(unresolved_effectors)})" if unresolved_effectors else ""))
    print(f"host symbols resolved: {len(hosts):,} of {raw.gene.nunique():,} "
          f"({len(hosts) / raw.gene.nunique():.1%})")

    frame = raw.assign(
        gene_id=raw["target"].map(effectors),
        host_id=raw["gene"].map(hosts),
    ).dropna(subset=["gene_id", "host_id"])
    frame = frame.rename(columns={"avg_log2FC": "hosttx_log2fc", "p_val_bh": "hosttx_q",
                                  "gene": "host_symbol", "target": "effector"})
    frame = frame[["gene_id", "host_id", "host_symbol", "effector", "hosttx_log2fc", "hosttx_q",
                   "pct.1", "pct.2"]]
    kept = frame[frame["hosttx_q"] < MAX_Q].copy()
    kept["bridge"] = "host"

    print(f"\n{len(frame):,} rows with both ends resolved; "
          f"{len(kept):,} at q < {MAX_Q} over {kept.gene_id.nunique()} effectors "
          f"and {kept.host_id.nunique():,} host genes")
    per = kept.groupby("effector").size().sort_values(ascending=False)
    print("\nhost genes moved, per effector:")
    print(per.head(8).to_string())
    print(f"median across effectors: {int(per.median())}")

    # DOES IT RECOVER THE BIOLOGY OF ITS OWN PAPER? The study reports that SOS1 (TGME49_222100) is
    # necessary for sustained host STAT6 signalling. If the import is keyed correctly, knocking SOS1
    # out should LOWER the canonical targets of that pathway, and it does.
    stat = ["SOCS1", "CISH", "CCL17", "CCL22", "CCL24", "ALOX15", "MRC1", "CD200R1", "STAT6"]
    sos1 = kept[kept["effector"] == "222100"]
    found = sos1[sos1["host_symbol"].isin(stat)][["host_symbol", "hosttx_log2fc", "hosttx_q"]]
    print(f"\nSOS1 knockout, canonical STAT targets among the {len(sos1):,} host genes it moves:")
    print(found.to_string(index=False) if len(found) else "  none")
    if len(found) and (found["hosttx_log2fc"] < 0).all():
        print("  all of them DOWN, which is what losing sustained STAT6 signalling predicts")

    out = os.path.join(ROOT, "results", f"hosttx_per_host_gene_{date.today():%Y_%m_%d}")
    os.makedirs(out, exist_ok=True)
    kept.to_parquet(os.path.join(out, "hosttx_host_gene_effects.parquet"), compression="zstd")
    per.rename("host_genes_moved").to_csv(os.path.join(out, "per_effector.tsv"), sep="\t")
    pd.Series(unresolved_hosts, name="symbol").to_csv(
        os.path.join(out, "unresolved_host_symbols.tsv"), sep="\t", index=False)
    json.dump({"source": SOURCE, "pmid": "37827122", "accession": "GSE229505",
               "rows_in": int(len(raw)), "effectors_in": int(raw.target.nunique()),
               "host_symbols_in": int(raw.gene.nunique()),
               "effectors_resolved": len(effectors),
               "unresolved_effectors": unresolved_effectors,
               "ambiguous_effectors": {k: v for k, v in ambiguous.items()},
               "resolved_effectors": effectors,
               "host_symbols_resolved": len(hosts),
               "rows_both_ends_resolved": int(len(frame)),
               "max_q": MAX_Q, "rows_kept": int(len(kept)),
               "effectors_kept": int(kept.gene_id.nunique()),
               "host_genes_kept": int(kept.host_id.nunique()),
               "median_host_genes_per_effector": int(per.median()),
               "sos1_stat_targets": {r.host_symbol: float(r.hosttx_log2fc)
                                     for r in found.itertuples()}},
              open(os.path.join(out, "summary.json"), "w"), indent=2)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
