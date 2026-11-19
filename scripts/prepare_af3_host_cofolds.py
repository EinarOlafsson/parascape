"""Figure 8b: prepare the AF3 co-folds of parasite effectors with their host partners.

Einar asked for the AF3 co-folding work in this figure. This script does everything that does NOT
need the GPU: it chooses the pairs, fetches the sequences, writes one AF3 input per complex and
writes the runner. The GPU step is left to be started deliberately, because the documented queue on
this machine is committed elsewhere and a large AF3 job would cut in front of it.

WHICH PAIRS, AND WHY A CONTROL. The only parasite-host PHYSICAL evidence in the corpus is
`myr1_host_ip`: one bait, 733 human proteins, with an IP enrichment per partner. Co-folding the
apparent best partners alone would produce scores with nothing to compare them against -- AF3
returns an interface score for any two chains, including two that never meet. So a matched control
arm is prepared from the same experiment: proteins that went through the same pulldown, the same
mass spectrometry and the same identification, and differ in the one thing that is supposed to
matter.

AND THE ARMS ARE CLASSES, NOT A RANKING. The first version of this script took the twelve most
enriched and the twelve most depleted. Both were wrong. 80 of the 720 partners sit above +20 log2,
which is not a ratio of two measurements but "not detected in the control" with the control
imputed at a detection limit -- so the top twelve would have been the twelve largest imputations.
And only 6 partners fall below -20, so a "most depleted" arm could not be filled, and would in any
case compare one imputation against another. The arms are now: seen ONLY IN THE BAIT, against seen
IN BOTH AND LEVEL, sampled from each class with a fixed seed.

A CAVEAT THAT BELONGS IN THE FIGURE, not just here: MYR1 is a translocon component. It threads
other effectors across the vacuole membrane, so a host protein enriched in a MYR1 pulldown may be
bound by something MYR1 carried rather than by MYR1. A high interface score is evidence the pair
CAN be modelled together; it is not evidence that MYR1 is the binder.

THE MSA STEP IS CAPPED. Stock jackhmmer has no `--seq_limit` and has been measured on this machine
at 118 GB, which OOM-killed the box. The runner this writes caps the MSA subprocess, and refuses to
start if the cap is not in place.

Writes `results/af3_host_cofolds_<date>/` and AF3 inputs under `af_input/host_cofolds_<date>/`.
"""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import date

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

#: Where AF3 lives on this machine, and where its inputs go.
AF3_ROOT = os.environ.get("PARASCAPE_AF3_ROOT", "/media/carruthers/mnt1/af3")

#: The bait, and the table of its host partners.
BRIDGES = os.path.join(ROOT, "parascape", "data", "host_bridges.parquet")

#: How many enriched pairs, and the same number of depleted ones as the matched control. Twelve
#: each is a day of GPU rather than a week, and twelve against twelve separates a shifted
#: distribution from a tie without pretending to resolve a small difference.
N_EACH = 12

#: Above this, an "enrichment" is not a ratio of two measurements. 80 of the 720 partners sit above
#: +20 log2 -- a 1,000,000-fold enrichment -- with 75 distinct values among them. That is the
#: signature of ONLY IN THE BAIT: the protein was not detected in the control, the control
#: intensity was imputed at a detection limit, and the ratio is the imputation, not the biology.
#: Ranking those 80 and taking the top twelve would be ranking by imputation noise, which is the
#: same error as reading a count of zero as a measurement.
ONLY_IN_BAIT_ABOVE = 20.0

#: The control arm is partners DETECTED IN BOTH and not enriched -- the proteins that were in the
#: same lysate, went through the same pulldown and the same mass spectrometry, and came out level.
#: The obvious alternative, the most DEPLETED partners, is wrong twice over: only 6 of 720 are
#: below -20, so the arm could not be filled from them, and "only in the control" is just as
#: imputed as "only in the bait" -- it would compare one imputation against another.
CONTROL_WITHIN = 0.5

#: UniProt, for the host sequences. The parasite ones are already in the node table.
UNIPROT = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
PAUSE_S = 0.34          # polite, and well inside what the service asks for


def host_sequence(accession: str) -> str | None:
    """One host sequence from UniProt, or None when the accession no longer resolves."""
    try:
        with urllib.request.urlopen(UNIPROT.format(accession=accession), timeout=30) as response:
            body = response.read().decode()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError):
        return None
    lines = [line.strip() for line in body.splitlines() if line and not line.startswith(">")]
    return "".join(lines) or None


def choose(bridges: pd.DataFrame) -> pd.DataFrame:
    """Partners seen ONLY in the bait, against partners seen in both and not enriched.

    Neither arm is a ranking, because the quantity cannot carry one at either end: see
    `ONLY_IN_BAIT_ABOVE` and `CONTROL_WITHIN`. Both arms are sampled with a fixed seed from their
    class, so the choice is reproducible and is not a quiet ranking by imputation noise.
    """
    usable = bridges.dropna(subset=["host_ip_enrichment_log2"]).copy()
    value = usable["host_ip_enrichment_log2"]
    only_in_bait = usable[value > ONLY_IN_BAIT_ABOVE]
    level = usable[value.abs() <= CONTROL_WITHIN]
    if len(only_in_bait) < N_EACH or len(level) < N_EACH:
        raise SystemExit(f"not enough pairs: {len(only_in_bait)} only-in-bait, {len(level)} level")
    seed = 0
    enriched = only_in_bait.sample(n=N_EACH, random_state=seed).assign(arm="only in the bait")
    control = level.sample(n=N_EACH, random_state=seed).assign(arm="in both, not enriched")
    print(f"  {len(only_in_bait)} partners are only in the bait (> {ONLY_IN_BAIT_ABOVE:g} log2), "
          f"{len(level)} are in both and level (within {CONTROL_WITHIN:g});")
    print(f"  {N_EACH} sampled from each with seed {seed}, because neither class carries a rank.")
    return pd.concat([enriched, control], ignore_index=True)


def main() -> None:
    cache = os.environ.get("PARASCAPE_CACHE", "")
    if not cache:
        raise SystemExit("set PARASCAPE_CACHE")

    bridges = pd.read_parquet(BRIDGES)
    pairs = choose(bridges)
    print(f"{len(bridges):,} host partners recorded; {len(pairs)} pairs chosen")

    nodes = pd.read_parquet(os.path.join(cache, "nodes.parquet"))
    nodes = nodes.set_index(nodes["gene_id"].astype(str))

    stamp = f"{date.today():%Y_%m_%d}"
    inputs = os.path.join(AF3_ROOT, "af_input", f"host_cofolds_{stamp}")
    out = os.path.join(ROOT, "results", f"af3_host_cofolds_{stamp}")
    os.makedirs(inputs, exist_ok=True)
    os.makedirs(out, exist_ok=True)

    rows, written, missing = [], 0, []
    for pair in pairs.itertuples():
        parasite = str(pair.gene_id)
        if parasite not in nodes.index:
            missing.append((parasite, "parasite gene not in the node table"))
            continue
        parasite_seq = nodes.at[parasite, "sequence"]
        if not isinstance(parasite_seq, str) or not parasite_seq:
            missing.append((parasite, "no parasite sequence"))
            continue
        host_seq = host_sequence(pair.host_id)
        time.sleep(PAUSE_S)
        if not host_seq:
            missing.append((pair.host_id, "UniProt did not serve a sequence"))
            continue
        name = f"{parasite}__{pair.host_id}"
        job = {"name": name, "modelSeeds": [1],
               "sequences": [{"protein": {"id": "A", "sequence": parasite_seq}},
                             {"protein": {"id": "B", "sequence": host_seq}}]}
        with open(os.path.join(inputs, f"{name}.json"), "w") as handle:
            json.dump(job, handle)
        written += 1
        rows.append({"name": name, "parasite": parasite, "host_id": pair.host_id,
                     "host_name": pair.host_name, "arm": pair.arm,
                     "host_ip_enrichment_log2": float(pair.host_ip_enrichment_log2),
                     "parasite_residues": len(parasite_seq), "host_residues": len(host_seq)})
        print(f"  {name:34s} {pair.arm:19s} enrichment {pair.host_ip_enrichment_log2:7.2f}  "
              f"{len(parasite_seq):5d} + {len(host_seq):5d} aa")

    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(out, "pairs.tsv"), sep="\t", index=False)
    json.dump({"bait": "MYR1 (translocon component -- see the caveat in the module docstring)",
               "partners_recorded": int(len(bridges)), "pairs_written": written,
               "n_each_arm": N_EACH, "only_in_bait_above": ONLY_IN_BAIT_ABOVE,
               "control_within": CONTROL_WITHIN,
               "af3_inputs": inputs, "not_prepared": missing,
               "gpu": "NOT started: the documented GPU queue is committed elsewhere"},
              open(os.path.join(out, "summary.json"), "w"), indent=2)

    print(f"\n{written} AF3 inputs written to {inputs}")
    if missing:
        print(f"{len(missing)} pairs not prepared:")
        for what, why in missing:
            print(f"  {what}: {why}")
    print(f"\nPairs table: {out}/pairs.tsv")
    print("\nThe GPU step is NOT started. Nothing here has run AlphaFold.")


if __name__ == "__main__":
    main()
