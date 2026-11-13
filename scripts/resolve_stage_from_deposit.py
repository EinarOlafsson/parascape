"""Instruction 65.04: ask the DEPOSIT which life stage a column was measured in.

The eligibility read (`verify_biological_eligibility.py`) can only compare declarations, and for 27
candidates the column name declares no stage. This script asks the authority that does know: the
submitter, through the per-sample characteristics recorded in GEO.

It is the same principle the origin-paper check settled on -- ask the deposit before the paper --
applied to biology rather than provenance. The deposit is better evidence than a column name and
better evidence than my reading, because it is the person who did the experiment saying what the
sample was.

**What it does with the answer matters more than the fetch.** Three outcomes:

  * the series declares ONE stage across every sample, so any column from it carries that stage;
  * the series declares SEVERAL, and the column maps to particular samples by name, so the column
    carries the stage of those samples -- and if they disagree with the slot, that is a finding;
  * the series declares several and the column cannot be mapped, which stays unknown.

Writes ``results/dataset_selection_2026_11_08/stage_from_deposit.tsv``.
"""
from __future__ import annotations

import os
import re
import sys
import time
import urllib.request

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

sys.path.insert(0, os.path.join(ROOT, "scripts"))

import verify_biological_eligibility as E  # noqa: E402

from parascape import datasets, slots as S  # noqa: E402

OUT = os.path.join(ROOT, "results", "dataset_selection_2026_11_08")
SOFT = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        "?acc={acc}&targ=gsm&form=text&view=brief")

# What the submitter writes, mapped to the stage vocabulary the slots use. The left-hand side is
# matched against the whole characteristics line, so it catches `developmental stage: Tachyzoites`
# and `cell type: tachyzoites` alike.
DECLARED = [
    (r"enterop?ithelial|enteroepithelial", "enteroepithelial"),
    (r"tissue cyst", "bradyzoite"),
    (r"pre-?bradyzoite", "pre-bradyzoite"),
    (r"bradyzoite", "bradyzoite"),
    (r"tachyzoite", "tachyzoite"),
    (r"oocyst", "oocyst"),
    (r"sporozoite", "sporozoite"),
    (r"merozoite", "merozoite"),
    (r"gametocyte", "gametocyte"),
    (r"schizont", "schizont"),
    (r"trophozoite", "trophozoite"),
    (r"\bring\b", "ring"),
]


def fetch(accession: str) -> str:
    """The brief SOFT record for every sample of a series, or an empty string."""
    for attempt in range(3):
        try:
            with urllib.request.urlopen(SOFT.format(acc=accession), timeout=90) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(3 * (attempt + 1))
    return ""


def samples(text: str) -> list:
    """(title, declared stages) for each sample, in the order the series lists them."""
    out, title, said = [], None, set()
    for line in text.splitlines():
        if line.startswith("!Sample_title"):
            if title is not None:
                out.append((title, said))
            title, said = line.split("=", 1)[1].strip(), set()
        elif line.startswith("!Sample_characteristics") and title is not None:
            value = line.split("=", 1)[1].strip().lower()
            for pattern, stage in DECLARED:
                if re.search(pattern, value):
                    said.add(stage)
                    break
    if title is not None:
        out.append((title, said))
    return out


def tokens(pattern: str, key: str) -> list:
    """The part of a column pattern that could name a sample.

    `rna108740_Day` -> ['day']; `rna108740_Tissue_cysts` -> ['tissue','cysts']. The accession-derived
    prefix is dropped because every sample of the series shares it and it discriminates nothing.
    """
    stem = re.sub(r"^[a-z]+\d+_?", "", str(pattern or "").lower())
    return [t for t in re.split(r"[^a-z0-9]+", stem) if t and not t.isdigit()]


def matching(pattern: str, key: str, records: list) -> list:
    """The samples whose title the column pattern selects."""
    toks = tokens(pattern, key)
    if not toks:
        return []
    return [(t, s) for t, s in records if all(tok in t.lower() for tok in toks)]


def main() -> int:
    by_key = {e.key: e for e in datasets.REGISTRY}
    by_slot = {s.key: s for s in S.all_slots()}
    elig = pd.read_csv(os.path.join(OUT, "biological_eligibility.tsv"), sep="\t")
    need = elig[elig.stage_verdict == "unknown"]

    wanted = sorted({r.source_id for r in need.itertuples()
                     if str(getattr(by_key.get(r.source_id), "accession", "") or "").startswith("GSE")})
    series = {}
    for key in wanted:
        acc = by_key[key].accession
        print(f"asking the deposit: {acc} ({key})")
        series[key] = samples(fetch(acc))
        time.sleep(0.5)

    rows = []
    for r in need.itertuples():
        entry, slot = by_key.get(r.source_id), by_slot.get(r.slot)
        records = series.get(r.source_id)
        # The stages the SLOT asks about, read from its own facets. Parsing them out of the
        # evidence sentence instead reported `Tg_transcription_oocyst_sporozoite` as asking only
        # for sporozoite, so oocyst data -- which is exactly what that slot is for -- came out as a
        # disagreement. A slot that covers two stages is satisfied by either.
        asked = sorted(E.stage_words(slot)) if slot is not None else []
        row = {"slot": r.slot, "dataset_id": r.dataset_id, "source_id": r.source_id,
               "accession": getattr(entry, "accession", ""), "slot_context": slot.context,
               "slot_asks": ", ".join(asked)}
        if not records:
            row.update(outcome="no deposit to ask",
                       evidence="the source records no GEO series, so the deposit cannot be asked")
            rows.append(row)
            continue

        declared = sorted({s for _, said in records for s in said})
        hit = matching(r.pattern if hasattr(r, "pattern") else "", r.source_id, records)
        if not hit:
            pat = str(r.dataset_id).split(":", 1)[-1]
            hit = matching(pat, r.source_id, records)
        hit_stages = sorted({s for _, said in hit for s in said})

        row["series_declares"] = ", ".join(declared) or "no stage in any sample"
        row["samples_matched"] = len(hit)
        row["column_stage"] = ", ".join(hit_stages)

        if not declared:
            row.update(outcome="the deposit declares no stage",
                       evidence="no sample of the series names a life stage; the paper is the only "
                                "authority left")
        elif len(declared) == 1:
            row.update(outcome="resolved", column_stage=declared[0],
                       evidence=f"every sample of {entry.accession} declares '{declared[0]}'")
        elif hit_stages and len(hit_stages) == 1:
            row.update(outcome="resolved by sample names",
                       evidence=f"the column selects {len(hit)} sample(s), all declared "
                                f"'{hit_stages[0]}'")
        else:
            row.update(outcome="several stages, column unmapped",
                       evidence=f"{entry.accession} declares {declared}; the column selects "
                                f"{len(hit)} sample(s)")
        rows.append(row)

    frame = pd.DataFrame(rows)
    frame.to_csv(os.path.join(OUT, "stage_from_deposit.tsv"), sep="\t", index=False)

    print()
    print(frame["outcome"].value_counts().to_string())

    # The point of the exercise: a column whose samples the submitter assigned to a stage the slot
    # does not ask for.
    # Only rows where the deposit actually ANSWERED. A row with no deposit to ask has an empty
    # column stage, and listing those as disagreements buried the one real finding in 23 rows that
    # say nothing.
    answered = frame[frame.column_stage.notna() & (frame.column_stage.astype(str) != "")
                     & (frame.slot_asks.astype(str) != "")]
    bad = answered[~answered.apply(lambda x: any(a and a in str(x.column_stage)
                                                 for a in str(x.slot_asks).split(", ")), axis=1)]
    if len(bad):
        print()
        print("THE SLOT ASKS FOR ONE STAGE AND THE DEPOSIT DECLARES ANOTHER:")
        print(bad[["slot", "dataset_id", "slot_asks", "column_stage", "evidence"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
