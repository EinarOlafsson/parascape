"""Distributable packs for the five further spaces (instruction 53, WP1 and WP4-6).

A pack is what lets a space leave `results/` and be installed by someone else, checked against a published
SHA256. `packs.build` refuses any column whose licence is not verified open (`O`) or share-alike (`SA`), which
is the whole point: a distributed pack must not carry data we have no right to distribute.

**What that excludes, and why it is left out rather than guessed at.** The Ensembl gene builds, Bgee expression
calls and Ensembl Compara orthologs are open and go in. The measurement columns derived from GEO series --
`expr_pv_liver` and the four `expr_midgut_*` -- do NOT: a GEO submission carries no stated licence of its own,
and deciding that a depositor's data may be redistributed is a judgement about someone else's rights, not a
technical detail. Those columns stay in `results/` for local use and the pack says which are missing.

    python scripts/build_space_packs.py
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
SOURCE = os.path.join(ROOT, "results", "host_gene_spaces_2026_11_09")
OUT = os.path.join(ROOT, "results", "space_packs_2026_11_09")
VERSION = _dt.date.today().strftime("%Y.%m.%d")

#: Column -> licence code, by what produced it. Only `O` and `SA` may be distributed.
OPEN_SOURCES = {
    "gene_id": ("O", "Ensembl gene build, no restrictions on use"),
    "gene_name": ("O", "Ensembl gene build, no restrictions on use"),
    "biotype": ("O", "Ensembl gene build, no restrictions on use"),
    "seqid": ("O", "Ensembl gene build, no restrictions on use"),
    "start": ("O", "Ensembl gene build, no restrictions on use"),
    "end": ("O", "Ensembl gene build, no restrictions on use"),
    "strand": ("O", "Ensembl gene build, no restrictions on use"),
    "bgee_top_score": ("O", "Bgee, CC0"),
    "bgee_tissues_present": ("O", "Bgee, CC0"),
    "human_ortholog": ("O", "Ensembl Compara, no restrictions on use"),
    "human_ortholog_high_confidence": ("O", "Ensembl Compara, no restrictions on use"),
}
#: Columns held back from every pack, with the reason a reader can check.
WITHHELD = {
    "expr_pv_liver": "derived from GEO GSE220039, which states no licence; redistribution is the depositor's "
                     "call and not ours to assume",
    "expr_midgut_36h": "derived from GEO GSE284537, which states no licence",
    "expr_midgut_48h": "derived from GEO GSE284537, which states no licence",
    "expr_midgut_4d": "derived from GEO GSE284537, which states no licence",
    "expr_midgut_7d": "derived from GEO GSE284537, which states no licence",
}


def _space(code: str):
    from parascape import space_catalogue
    return space_catalogue._host_space(code)


def main() -> int:
    """Build one pack per space, and say what each carries and what it leaves behind."""
    import pandas as pd

    from parascape.spaces import build_pack

    os.makedirs(OUT, exist_ok=True)
    report = {}
    for code in ("Rn", "Fc", "Ag", "As", "Pv"):
        space = _space(code)
        source = os.path.join(SOURCE, space.nodes)
        if not os.path.exists(source):
            print(f"{code}: no table at {source}")
            continue
        table = pd.read_parquet(source)
        held = [c for c in table.columns if c in WITHHELD]
        unknown = [c for c in table.columns if c not in OPEN_SOURCES and c not in WITHHELD]
        if unknown:
            raise SystemExit(f"{code}: no licence declared for {unknown}; declare it or withhold it")
        packed = table.drop(columns=held)
        staging = os.path.join(OUT, f"_stage_{code}")
        shutil.rmtree(staging, ignore_errors=True)
        os.makedirs(staging)
        packed.to_parquet(os.path.join(staging, space.nodes), index=False)
        licences = {c: OPEN_SOURCES[c][0] for c in packed.columns}
        notice = [f"ParaScape {space.species} gene space ({code}), version {VERSION}", "",
                  "Sources and licences, per column:"]
        notice += [f"  {c}: {OPEN_SOURCES[c][1]}" for c in packed.columns]
        if held:
            notice += ["", "Columns present in the build but NOT distributed here:"]
            notice += [f"  {c}: {WITHHELD[c]}" for c in held]
        notice_path = os.path.join(staging, "NOTICE.txt")
        open(notice_path, "w").write("\n".join(notice) + "\n")
        archive = os.path.join(OUT, f"parascape_{code}_{VERSION}.zip")
        sha = build_pack(space, staging, VERSION, licences, archive,
                         extra={"NOTICE.txt": notice_path})
        shutil.rmtree(staging, ignore_errors=True)
        report[code] = dict(species=space.species, kind=space.kind, genes=len(packed),
                            columns=list(packed.columns), withheld=held,
                            archive=os.path.basename(archive), sha256=sha,
                            bytes=os.path.getsize(archive))
        print(f"{code}: {len(packed):,} genes, {len(packed.columns)} columns"
              f"{', withheld ' + ', '.join(held) if held else ''} -> "
              f"{os.path.basename(archive)} ({os.path.getsize(archive) // 1024} KiB)")
    with open(os.path.join(OUT, "SHA256"), "w") as handle:
        for code, row in sorted(report.items()):
            handle.write(f"{row['sha256']}  {row['archive']}\n")
    json.dump(report, open(os.path.join(OUT, "packs.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
