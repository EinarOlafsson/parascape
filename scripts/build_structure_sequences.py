"""Genome-wide structure sequences for the information map (Einar, 2026-11-20).

For every ME49 AlphaFold DB model: the gene (from the CIF `_ma_target_ref_db_details.gene_name`),
the amino-acid sequence, a DSSP-style secondary-structure string (H/E/C, Kabsch-Sander hydrogen
bonds; the same code as toxoplasma_projects/af3_folding/analysis_05) and per-residue pLDDT. The
Foldseek 3Di string comes from `foldseek createdb` + `convert2fasta` on the same models.

Writes /media/carruthers/mnt2/parascape/structure_seqs/{secondary_structure.parquet, 3di.fasta}.
"""
import collections
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

STRUCT = "/home/carruthers/foldseek/me49_1/structures"
OUT = "/media/carruthers/mnt2/parascape/structure_seqs"
sys.path.insert(0, "/media/carruthers/mnt3/claude/toxoplasma_projects/af3_folding")
AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def dssp_like(Nx, CAx, Cx, Ox):
    n = len(CAx)
    ss = np.full(n, "C", dtype="<U1")
    if n < 6:
        return ss
    unit = lambda v: v / np.maximum(np.linalg.norm(v, axis=-1, keepdims=True), 1e-9)
    H = np.full((n, 3), np.nan)
    H[1:] = Nx[1:] + 1.01 * unit(unit(Nx[1:] - Cx[:-1]) + unit(Nx[1:] - CAx[1:]))
    d = lambda a, b: np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        E = 0.084 * 332 * (1 / d(Ox, Nx) + 1 / d(Cx, H) - 1 / d(Ox, H) - 1 / d(Cx, Nx))
    E = np.nan_to_num(E, nan=0.0, posinf=0.0, neginf=0.0)
    ii, jj = np.indices(E.shape)
    E[np.abs(ii - jj) < 2] = 0.0
    HB = E < -0.5
    t4 = np.zeros(n, bool)
    t4[:n - 4] = HB[np.arange(n - 4), np.arange(4, n)]
    for i in np.flatnonzero(t4[:-1] & t4[1:]):
        ss[i + 1:i + 5] = "H"
    R = np.zeros_like(HB); R[1:, :-1] = HB[:-1, 1:]
    U = np.zeros_like(HB); U[1:, :] = HB[:-1, :]
    W = np.zeros_like(HB); W[:, :-1] = HB[:, 1:]
    bridge = (HB & HB.T) | (R & R.T) | (U & W.T) | (U.T & W)
    bridge &= np.abs(ii - jj) >= 3
    br = bridge.any(1) | bridge.any(0)
    ss[br & (ss != "H")] = "E"
    return ss


def one(fname):
    path = os.path.join(STRUCT, fname)
    cols, started, gene, acc = [], False, None, None
    at = collections.defaultdict(dict); res = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("_ma_target_ref_db_details.gene_name"):
                gene = ln.split()[-1]
            elif ln.startswith("_ma_target_ref_db_details.db_accession"):
                acc = ln.split()[-1]
            elif ln.startswith("_atom_site."):
                cols.append(ln.strip().split(".")[1]); started = True
            elif started and ln.startswith("ATOM"):
                f = ln.split()
                if len(f) < len(cols):
                    continue
                dd = dict(zip(cols, f))
                a = dd.get("label_atom_id")
                k = int(dd["label_seq_id"])
                res[k] = dd.get("label_comp_id")
                if a in ("N", "CA", "C", "O"):
                    at[k][a] = (float(dd["Cartn_x"]), float(dd["Cartn_y"]), float(dd["Cartn_z"]),
                                float(dd["B_iso_or_equiv"]))
    resi = sorted(k for k, v in at.items() if {"N", "CA", "C", "O"} <= set(v))
    if not resi or len(resi) > 3000:          # very long models: skip the O(n^2) step
        return None
    g = lambda a: np.array([at[k][a][:3] for k in resi])
    ss = dssp_like(g("N"), g("CA"), g("C"), g("O"))
    return {"gene_id": gene, "uniprot": acc, "length": len(resi),
            "sequence": "".join(AA.get(res[k], "X") for k in resi), "ss": "".join(ss),
            "plddt": np.round([at[k]["CA"][3] for k in resi], 1).tolist()}


if __name__ == "__main__":
    files = sorted(f for f in os.listdir(STRUCT) if f.endswith(".cif"))
    with ProcessPoolExecutor(max_workers=8) as ex:
        rows = [r for r in ex.map(one, files, chunksize=20) if r]
    df = pd.DataFrame(rows)
    df.to_parquet(os.path.join(OUT, "secondary_structure.parquet"), index=False)
    print(f"{len(df)} models, {df.gene_id.nunique()} genes -> secondary_structure.parquet")
