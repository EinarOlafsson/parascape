"""Build Toxoplasma gene x raw-file protein tables from every finished uniform Sage re-search (parascape.ms_genes).

Reads prepared_ms_uniform/<project>/sage_*/lfq.tsv for every project with done.json; writes
prepared_ms_genes/<project>/<sage folder>/genes.parquet (log2, median-centred Top3) and peptides.parquet (unique
peptides per gene), and prepared_ms_genes/summary.csv (one row per project and instrument folder). Read-only on
the search outputs; re-run safely (outputs are replaced).

    python scripts/ms/build_ms_gene_tables.py [--root /nas_mnt/carruthers/Einar/parascape_omics]
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from parascape import ms_genes as M  # noqa: E402

LABELLED = re.compile(r"itraq|tmt|tandem.?mass.?tag", re.I)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default="/nas_mnt/carruthers/Einar/parascape_omics")
    args = ap.parse_args(argv)
    src = os.path.join(args.root, "prepared_ms_uniform")
    out = os.path.join(args.root, "prepared_ms_genes")
    os.makedirs(out, exist_ok=True)
    rows = []
    for done in sorted(glob.glob(os.path.join(src, "*", "done.json"))):
        project = os.path.basename(os.path.dirname(done))
        for folder in sorted(glob.glob(os.path.join(os.path.dirname(done), "sage_*"))):
            path = os.path.join(folder, "lfq.tsv")
            if not os.path.exists(path):
                rows.append({"project": project, "instrument": os.path.basename(folder), "note": "no lfq.tsv"})
                continue
            lfq = pd.read_csv(path, sep="\t")
            G, peptides = M.gene_table(lfq)
            dest = os.path.join(out, project, os.path.basename(folder))
            os.makedirs(dest, exist_ok=True)
            G.to_parquet(os.path.join(dest, "genes.parquet"))
            peptides.rename("unique_peptides").to_frame().to_parquet(os.path.join(dest, "peptides.parquet"))
            labelled = bool(LABELLED.search(" ".join(G.columns) + " " + project))
            rows.append({**M.summary(project, lfq, G, labelled), "instrument": os.path.basename(folder),
                         "note": ""})
            print(rows[-1], flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(out, "summary.csv"), index=False)
    print(f"{len(S)} project folders; {int((S.get('genes_quantified', 0) >= 100).sum())} quantify >= 100 genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
