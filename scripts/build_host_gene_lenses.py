"""Instruction 53, WP9 route 1: the host GENE spaces get measurements and the same two lenses.

Route 2 (the host PROTEIN spaces, keyed on UniProt accessions) shipped on 2026-11-13. This is
route 1: the Ensembl gene spaces `HUMAN_SPACE` (58,988 ENSG) and `MOUSE_SPACE` (18,428 ENSMUSG),
laid out the same two ways as the protein lenses, so a reader comparing a protein window with a
gene window is comparing two layouts that make the same two claims.

* **`tissue`** -- two genes are similar when they are expressed in the same tissues and cell types.
* **`parasite_facing`** -- two genes are similar in how they respond to infection, and in what the
  parasite does to them.

**Why this is not the protein lens with the identifiers swapped.** The protein lenses read the host
protein tables, whose columns were already summarised per UniProt accession. The gene tables carried
almost nothing: 2 GTEx columns and 7 Human Protein Atlas categories for human, Bgee top score and
human-ortholog fields for mouse. So the measurements are built here from the raw files already on
disk under `$PARASCAPE_DATA/spaces/Hs` and `spaces/Mm`. Nothing is downloaded.

**The mouse parasite-facing lens was refused for the protein table and is built here.** The protein
table carried ONE parasite-facing measurement (the macrophage surface, 147 of 16,473 proteins), which
is why `build_host_layouts.py` refused that lens. The mouse gene table now carries ten, from six
Toxoplasma, Cryptosporidium and Plasmodium series, so the lens is real rather than imputation wearing
the shape of a map. The floor that refused it is kept and still applies.

**Every value is measured on the species whose table it enters.** No column is carried across the
protein-to-gene link here: a gene lens built from protein measurements spread over a many-to-many
mapping would be a map of the mapping. The link is a separate feature (`parascape/host_links.py`).

    python scripts/build_host_gene_lenses.py --probe     # coverage only, no layout
    python scripts/build_host_gene_lenses.py             # both species, both lenses
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import zipfile

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from parascape import paths  # noqa: E402
from parascape.embedding import EmbeddingSpec, embed as project  # noqa: E402

#: A lens must measure at least this share of the genes, or it is refused rather than drawn. Carried
#: over from `build_host_layouts.py`, where it caught a one-column mouse lens that produced 16,473
#: distinct positions out of 147 measured genes.
MINIMUM_MEASURED = 0.10

#: Where the raw deposits live. `PARASCAPE_DATA` by default.
SPACES = os.path.join(paths.dataset_root(), "spaces")

#: The gene packs the measurements are added to. Their node tables are the baseline: nothing in them
#: may be lost, which `spaces.validate_nodes(nodes, previous)` enforces.
PACKS = {"Hs": os.path.join(ROOT, "results", "human_pack_2026_11_07",
                            "parascape_Hs_2026.11.07.zip"),
         "Mm": os.path.join(ROOT, "results", "mouse_gene_space_2026_11_07",
                            "parascape_Mm_2026.11.07.zip")}
NODES = {"Hs": "hs_gene_nodes.parquet", "Mm": "mm_gene_nodes.parquet"}
LINKS = {"Hs": "hs_protein_gene_links.parquet", "Mm": "mm_protein_gene_links.parquet"}


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _space(*parts) -> str:
    return os.path.join(SPACES, *parts)


# --------------------------------------------------------------------------- identifier mapping
def _unversioned(ids) -> pd.Series:
    """ENSG00000223972.5 -> ENSG00000223972. The version is the annotation's, not the gene's."""
    return pd.Series(ids, dtype="string").str.split(".", n=1).str[0]


def symbol_map(code: str, universe: set) -> dict:
    """lowercased symbol or synonym -> gene id, for the symbols that name exactly one gene.

    A symbol naming several genes in the universe is DROPPED rather than resolved to the first one.
    Picking one would attach a real measurement to an arbitrary gene, and a reader has no way to see
    that happened. The counts are reported, so what the ambiguity costs is visible.
    """
    pairs: list = []
    if code == "Hs":
        hgnc = pd.read_csv(_space("Hs", "reference", "HGNC_2026-09-25", "hgnc_complete_set.txt"),
                           sep="\t", dtype=str, low_memory=False)
        hgnc = hgnc[hgnc.ensembl_gene_id.notna()]
        for column in ("symbol", "alias_symbol", "prev_symbol"):
            if column not in hgnc.columns:
                continue
            part = hgnc[[column, "ensembl_gene_id"]].dropna()
            part = part.assign(**{column: part[column].str.split("|")}).explode(column)
            pairs.append(part.rename(columns={column: "name", "ensembl_gene_id": "gene_id"}))
    else:
        mgi = pd.read_csv(_space("Mm", "reference", "MGI_2026-09-21", "MRK_ENSEMBL.rpt"),
                          sep="\t", header=None, dtype=str, low_memory=False)
        # MGI accession (0), symbol (1), name (2), cM (3), chromosome (4), Ensembl gene id (5).
        part = mgi[[1, 5]].dropna()
        part.columns = ["name", "gene_id"]
        pairs.append(part)
        synonyms = pd.read_csv(_space("Mm", "reference", "MGI_2026-09-21", "MRK_List2.rpt.gz"),
                               sep="\t", dtype=str, low_memory=False)
        column = next((c for c in synonyms.columns if "Synonym" in c), None)
        symbol = next((c for c in synonyms.columns if c.strip() == "Marker Symbol"), None)
        if column and symbol:
            syn = synonyms[[column, symbol]].dropna()
            syn = syn.assign(**{column: syn[column].str.split("|")}).explode(column)
            byname = dict(zip(part.name.str.lower(), part.gene_id))
            syn["gene_id"] = syn[symbol].str.lower().map(byname)
            pairs.append(syn[[column, "gene_id"]].dropna().rename(columns={column: "name"}))
    table = pd.concat(pairs, ignore_index=True)
    table["name"] = table.name.str.strip().str.lower()
    table = table[table.gene_id.isin(universe) & (table.name != "")]
    table = table.drop_duplicates(["name", "gene_id"])
    counts = table.name.value_counts()
    unique = table[table.name.map(counts) == 1]
    print(f"  {code}: {len(unique):,} symbols name one gene; "
          f"{int((counts > 1).sum()):,} name several and are left out")
    return dict(zip(unique.name, unique.gene_id))


def by_symbol(values: pd.Series, mapping: dict, universe) -> pd.Series:
    """A symbol-keyed series re-keyed on gene id, keeping the mean where two symbols meet one gene."""
    frame = pd.DataFrame({"name": pd.Series(values.index, dtype="string").str.strip().str.lower(),
                          "value": pd.to_numeric(pd.Series(values.to_numpy()), errors="coerce")})
    frame["gene_id"] = frame.name.map(mapping)
    frame = frame.dropna(subset=["gene_id", "value"])
    return frame.groupby("gene_id").value.mean().reindex(universe)


def by_gene(values: pd.Series, universe) -> pd.Series:
    """An ENSG/ENSMUSG-keyed series aligned to the universe, versions stripped, duplicates averaged."""
    frame = pd.DataFrame({"gene_id": _unversioned(values.index).to_numpy(),
                          "value": pd.to_numeric(pd.Series(values.to_numpy()), errors="coerce")})
    frame = frame.dropna()
    return frame.groupby("gene_id").value.mean().reindex(universe)


# --------------------------------------------------------------------------- counts to a contrast
def _cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts per million per column. The libraries differ in depth; a raw ratio reads that depth."""
    total = counts.sum(axis=0).replace(0, np.nan)
    return counts.divide(total, axis=1) * 1e6


def log2_ratio(frame: pd.DataFrame, case: list, control: list, logged: bool = False) -> pd.Series:
    """log2 of the mean of `case` over the mean of `control`, per row.

    `logged` says the values already are on a log scale, in which case the difference of the means
    IS the log ratio and taking a ratio again would be wrong. Both appear among these deposits, and
    getting it the wrong way round is silent: the numbers stay plausible.

    Only the named columns are read, and each is coerced to a number. Several of these deposits put
    a gene name, a biotype or a date beside their counts, and a library size summed over those is
    not a library size.
    """
    if not case or not control:
        raise ValueError("a contrast needs samples on both sides")
    counts = frame[list(case) + list(control)].apply(pd.to_numeric, errors="coerce")
    if logged:
        return counts[case].mean(axis=1) - counts[control].mean(axis=1)
    scaled = _cpm(counts)
    return np.log2((scaled[case].mean(axis=1) + 1.0) / (scaled[control].mean(axis=1) + 1.0))


def _tar_table(path: str, member: str, key: str, value: str, comment=None) -> pd.Series:
    with tarfile.open(path) as tar:
        with tar.extractfile(member) as raw, gzip.open(raw, "rt") as fh:
            table = pd.read_csv(fh, sep="\t", dtype=str, comment=comment, low_memory=False)
    column = value if value in table.columns else table.columns[-1]
    return pd.Series(pd.to_numeric(table[column], errors="coerce").to_numpy(),
                     index=table[key].to_numpy())


def _tar_matrix(path: str, pick, key: str, value: str, comment=None) -> pd.DataFrame:
    """One column per archive member `pick` accepts, named by what `pick` returns."""
    out = {}
    with tarfile.open(path) as tar:
        for info in tar.getmembers():
            name = pick(info.name)
            if not name:
                continue
            with tar.extractfile(info) as raw, gzip.open(raw, "rt") as fh:
                table = pd.read_csv(fh, sep="\t", dtype=str, comment=comment, low_memory=False)
            column = value if value in table.columns else table.columns[-1]
            out[name] = pd.Series(pd.to_numeric(table[column], errors="coerce").to_numpy(),
                                  index=table[key].to_numpy())
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- the tissue lens
def human_tissue(universe, cells: int = 40) -> tuple:
    """GTEx v10 tissue medians and Human Protein Atlas single-cell types, both keyed on ENSG.

    Expression is used as log1p of the value, not raw: a map built on raw TPM is a map of the few
    thousand genes with a large number, and every other gene sits in one corner of it.
    """
    out, records = {}, {}
    gtex = _space("Hs", "expression", "GTEx_v10",
                  "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
    table = pd.read_csv(gtex, sep="\t", skiprows=2, dtype=str, low_memory=False)
    table.index = _unversioned(table.Name).to_numpy()
    # The two GTEx columns the pack already carries, under their existing names: adding them again
    # under a new name would give the same measurement two votes in the layout.
    already = {"Cells_Cultured_fibroblasts": "gtex_fibroblast_median_tpm",
               "Liver_Hepatocyte": "gtex_lcm_hepatocyte_median_tpm"}
    for tissue in table.columns:
        if tissue in ("Name", "Description") or tissue in already:
            continue
        name = "gtex_" + re.sub(r"[^a-z0-9]+", "_", tissue.lower()).strip("_") + "_median_tpm"
        out[name] = np.log1p(by_gene(pd.Series(table[tissue].to_numpy(), index=table.index),
                                     universe))
        records[name] = dict(source="GTEx v10 median TPM (log1p)", context=tissue,
                             accession="phs000424", license="O")

    archive = _space("Hs", "expression", "HPA_v25.1_single_cell_type", "rna_single_cell_type.tsv.zip")
    with zipfile.ZipFile(archive) as z, z.open("rna_single_cell_type.tsv") as raw:
        hpa = pd.read_csv(io.TextIOWrapper(raw, "utf-8"), sep="\t", dtype=str, low_memory=False)
    hpa["nCPM"] = pd.to_numeric(hpa["nCPM"], errors="coerce")
    wide = hpa.pivot_table(index="Gene", columns="Cell type", values="nCPM", aggfunc="mean")
    # Capped at the cell types measured in the most genes, by the same coverage rule the mouse lens
    # uses, so the two species are cut the same way. Every GTEx tissue is kept: there are 66 of them
    # and they are the bulk half of the claim.
    wide = wide[wide.notna().sum().sort_values(ascending=False).head(cells).index]
    for cell in wide.columns:
        name = "hpa_sc_" + re.sub(r"[^a-z0-9]+", "_", str(cell).lower()).strip("_") + "_ncpm"
        out[name] = np.log1p(by_gene(wide[cell], universe))
        records[name] = dict(source="Human Protein Atlas v25.1 single cell type nCPM (log1p)",
                             context=str(cell), accession="", license="SA")
    return out, records


def mouse_tissue(universe, entities: int = 40) -> tuple:
    """Bgee 15.2 expression score per anatomical entity, keyed on ENSMUSG.

    The entities are chosen by COVERAGE -- the forty Bgee tested in the most genes -- rather than by
    taste. A curated tissue list would be a second claim on top of the lens, and nothing in the file
    says which forty a reader would have picked.

    A gene Bgee tested and called absent keeps its measured score; a gene Bgee never tested in that
    entity stays missing. That distinction is the one the mouse table was built on.
    """
    path = _space("Mm", "expression", "Bgee_15.2", "Mus_musculus_expr_simple.tsv.gz")
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False,
                        usecols=["Gene ID", "Anatomical entity name", "Expression score"])
    frame.columns = ["gene_id", "entity", "score"]
    frame["score"] = pd.to_numeric(frame.score, errors="coerce")
    frame = frame[frame.gene_id.isin(universe)]
    chosen = frame.entity.value_counts().head(entities).index
    frame = frame[frame.entity.isin(chosen)]
    wide = frame.pivot_table(index="gene_id", columns="entity", values="score", aggfunc="mean")
    out, records = {}, {}
    for entity in wide.columns:
        name = "bgee_" + re.sub(r"[^a-z0-9]+", "_", str(entity).lower()).strip("_") + "_score"
        out[name] = wide[entity].reindex(universe)
        records[name] = dict(source="Bgee 15.2 expression score", context=str(entity),
                             accession="", license="O")
    return out, records


# --------------------------------------------------------------------------- the parasite lens
#: Human series read for the parasite-facing lens. Each is a published deposit; each column names its
#: accession and its paper, which is what Einar's 2026-11-13 decision requires before it may ship.
HUMAN_SERIES = {
    "tg_mueller_glia_ghpt_log2fc": dict(
        accession="GSE145836", parasite="Toxoplasma gondii",
        citation="Smith JR et al. Ophthalmology Science 2025;5:100971 (PMID 41693810), "
                 "GEO GSE145836: human retinal Mueller glial cells infected with T. gondii GHPT",
        context="retinal Mueller glia, T. gondii GHPT over non-infected"),
    "tg_mueller_glia_gt1_log2fc": dict(
        accession="GSE145836", parasite="Toxoplasma gondii",
        citation="Smith JR et al. Ophthalmology Science 2025;5:100971 (PMID 41693810), "
                 "GEO GSE145836: human retinal Mueller glial cells infected with T. gondii GT1",
        context="retinal Mueller glia, T. gondii GT1 over non-infected"),
    "tg_cytotrophoblast_log2fc": dict(
        accession="GSE263615", parasite="Toxoplasma gondii",
        citation="Ander SE et al. mSphere 2024 (PMID 38771057), GEO GSE263615: cytotrophoblasts "
                 "from trophoblast stem cells, mock-infected or infected with T. gondii",
        context="cytotrophoblast, infected over mock"),
    "tg_syncytiotrophoblast_log2fc": dict(
        accession="GSE263615", parasite="Toxoplasma gondii",
        citation="Ander SE et al. mSphere 2024 (PMID 38771057), GEO GSE263615: "
                 "syncytiotrophoblasts, mock-infected or infected with T. gondii",
        context="syncytiotrophoblast, infected over mock"),
    "pf_erythroblast_basophilic_log2fc": dict(
        accession="GSE233739", parasite="Plasmodium falciparum",
        citation="Neveu G et al. 2023 (PMID 37493969), GEO GSE233739: transcriptional profiling of "
                 "primary human erythroblasts infected with P. falciparum",
        context="basophilic erythroblast, infected over uninfected"),
    "pf_erythroblast_orthochromatic_log2fc": dict(
        accession="GSE233739", parasite="Plasmodium falciparum",
        citation="Neveu G et al. 2023 (PMID 37493969), GEO GSE233739: transcriptional profiling of "
                 "primary human erythroblasts infected with P. falciparum",
        context="orthochromatic erythroblast, infected over uninfected"),
    "ch_enterocyte_log2fc": dict(
        accession="GSE241570", parasite="Cryptosporidium hominis",
        citation="GEO GSE241570 (BioProject of the human air-liquid-interface enterocyte series): "
                 "C. hominis infection of human small intestinal epithelial cells",
        context="air-liquid-interface enterocyte, C. hominis over non-infected"),
    "ch_enterocyte_second_log2fc": dict(
        accession="GSE241571", parasite="Cryptosporidium hominis",
        citation="GEO GSE241571 (companion series to GSE241570): C. hominis infection of human "
                 "small intestinal epithelial cells",
        context="air-liquid-interface enterocyte, second series, C. hominis over non-infected"),
    "cp_enterocyte_10h_log2fc": dict(
        accession="GSE185247", parasite="Cryptosporidium parvum",
        citation="Deng S et al. 2022 (GEO GSE185247): a genetic screen identifies a protective "
                 "type III interferon response to Cryptosporidium requiring TLR3",
        context="human enterocyte 10 h, infected over bystander"),
    "cp_enterocyte_48h_log2fc": dict(
        accession="GSE185247", parasite="Cryptosporidium parvum",
        citation="Deng S et al. 2022 (GEO GSE185247): a genetic screen identifies a protective "
                 "type III interferon response to Cryptosporidium requiring TLR3",
        context="human enterocyte 48 h, wild-type-infected over uninfected control"),
    "pv_blood_symptomatic_log2fc": dict(
        accession="GSE273483", parasite="Plasmodium vivax",
        citation="GEO GSE273483: systems approach identifies monocyte imbalance in symptomatic and "
                 "asymptomatic P. vivax malaria",
        context="whole blood, symptomatic malaria over healthy control"),
    "pv_blood_asymptomatic_log2fc": dict(
        accession="GSE273483", parasite="Plasmodium vivax",
        citation="GEO GSE273483: systems approach identifies monocyte imbalance in symptomatic and "
                 "asymptomatic P. vivax malaria",
        context="whole blood, asymptomatic malaria over healthy control"),
}

MOUSE_SERIES = {
    "tg_lung_log2fc": dict(
        accession="GSE119855", parasite="Toxoplasma gondii",
        citation="Singhania A et al. Nature Communications 2019 (GEO GSE119855): global "
                 "transcriptional profiling of the interferon network, T. gondii ME49, four tissues",
        context="whole lung, wild-type infected over uninfected"),
    "tg_blood_log2fc": dict(
        accession="GSE119855", parasite="Toxoplasma gondii",
        citation="Singhania A et al. Nature Communications 2019 (GEO GSE119855): global "
                 "transcriptional profiling of the interferon network, T. gondii ME49, four tissues",
        context="whole blood, wild-type infected over uninfected"),
    "tg_spleen_log2fc": dict(
        accession="GSE119855", parasite="Toxoplasma gondii",
        citation="Singhania A et al. Nature Communications 2019 (GEO GSE119855): global "
                 "transcriptional profiling of the interferon network, T. gondii ME49, four tissues",
        context="spleen, wild-type infected over uninfected"),
    "tg_liver_log2fc": dict(
        accession="GSE119855", parasite="Toxoplasma gondii",
        citation="Singhania A et al. Nature Communications 2019 (GEO GSE119855): global "
                 "transcriptional profiling of the interferon network, T. gondii ME49, four tissues",
        context="liver, wild-type infected over uninfected"),
    "tg_bmdm_log2fc": dict(
        accession="GSE103113", parasite="Toxoplasma gondii",
        citation="Nyonda MA et al. 2018 (GEO GSE103113): bone-marrow-derived macrophage "
                 "transcriptomes under T. gondii Pru infection and the GRA18 mutants",
        context="bone-marrow-derived macrophage, T. gondii Pru over uninfected"),
    "tg_chronic_28d_log2fc": dict(
        accession="GSE117504", parasite="Toxoplasma gondii",
        citation="GEO GSE117504: host transcriptional responses during chronic T. gondii infection "
                 "reveal sex-specific differences",
        context="chronic infection, day 28, infected over uninfected"),
    "cp_ileum_a_log2fc": dict(
        accession="GSE168680", parasite="Cryptosporidium parvum",
        citation="Gullicksrud JA et al. Mucosal Immunology 2022 (GEO GSE168680): enterocyte-ILC1 "
                 "crosstalk drives innate interferon-gamma-mediated control of Cryptosporidium",
        context="intestinal epithelium, experiment A, infected over uninfected"),
    "cp_ileum_b_log2fc": dict(
        accession="GSE168680", parasite="Cryptosporidium parvum",
        citation="Gullicksrud JA et al. Mucosal Immunology 2022 (GEO GSE168680): enterocyte-ILC1 "
                 "crosstalk drives innate interferon-gamma-mediated control of Cryptosporidium",
        context="intestinal epithelium, experiment B, infected over uninfected"),
    "pb_brain_cerebral_malaria_log2fc": dict(
        accession="GSE162535", parasite="Plasmodium berghei",
        citation="GEO GSE162535: RNA-seq of the mouse brain reveals inflammatory genes associated "
                 "with early artesunate treatment of experimental cerebral malaria",
        context="brain, experimental cerebral malaria over control"),
    "py_liver_log2fc": dict(
        accession="GSE44831", parasite="Plasmodium yoelii",
        citation="Liehl P et al. Nature Medicine 2014 (GEO GSE44831): interferon-mediated innate "
                 "immune responses against malaria parasite liver stages",
        context="liver, infected over mock, BALB/c"),
}


def human_parasite(universe, symbols: dict) -> dict:
    """The twelve human infection contrasts, each computed from its own deposit on disk."""
    out = {}

    glia = _tar_matrix(_space("Hs", "infection", "GSE145836", "GSE145836_RAW.tar"),
                       lambda n: (re.search(r"RS-(GHPTD|GT1D|TD)\d", n) or [None])[0],
                       "Geneid", "", comment="#")
    groups = {p: [c for c in glia.columns if c.startswith("RS-" + p)] for p in ("GHPTD", "GT1D", "TD")}
    for key, prefix in (("tg_mueller_glia_ghpt_log2fc", "GHPTD"),
                        ("tg_mueller_glia_gt1_log2fc", "GT1D")):
        out[key] = by_gene(log2_ratio(glia, groups[prefix], groups["TD"]), universe)

    # Keyed on "Gene ID", not on "ENSEMBL": that column holds the record URL rather than the
    # accession, and keying on it silently produced a column with no value for any gene.
    tropho = _tar_matrix(_space("Hs", "infection", "GSE263615", "GSE263615_RAW.tar"),
                         lambda n: (re.search(r"(CYT|SYN)_(mock|infected)_\d", n) or [None])[0],
                         "Gene ID", "Total gene reads")
    for key, cell in (("tg_cytotrophoblast_log2fc", "CYT"),
                      ("tg_syncytiotrophoblast_log2fc", "SYN")):
        case = [c for c in tropho.columns if c.startswith(cell + "_infected")]
        control = [c for c in tropho.columns if c.startswith(cell + "_mock")]
        out[key] = by_gene(log2_ratio(tropho, case, control), universe)

    ery = pd.read_csv(_space("Hs", "infection", "GSE233739", "GSE233739_RawGeneCounts.txt.gz"),
                      sep="\t", index_col=0, low_memory=False)
    for key, stage in (("pf_erythroblast_basophilic_log2fc", "BasoE"),
                       ("pf_erythroblast_orthochromatic_log2fc", "OrthoE")):
        case = [c for c in ery.columns if c.startswith(stage + "_Infected")]
        control = [c for c in ery.columns if c.startswith(stage + "_Uninfected")]
        out[key] = by_gene(log2_ratio(ery, case, control), universe)

    for key, accession in (("ch_enterocyte_log2fc", "GSE241570"),
                           ("ch_enterocyte_second_log2fc", "GSE241571")):
        # Named .csv.gz and written with tabs, and the second one carries an unnamed row number as
        # its first column. Read by separator and by column name rather than by position.
        table = pd.read_csv(_space("Hs", "infection", accession,
                                   f"{accession}_normalized_filtered_data.csv.gz"),
                            sep="\t", low_memory=False)
        out[key] = by_symbol(pd.Series(table["LogFC"].to_numpy(), index=table["geneID"].to_numpy()),
                             symbols, universe)

    # The two arms of GSE185247 are NOT the same experiment at two times, which the file names
    # suggest and the columns deny: the 10 h file compares infected with bystander cells, the 48 h
    # file compares wild-type-infected with uninfected controls. Each is read by its own columns.
    for key, hours, case_mark, control_mark in (
            ("cp_enterocyte_10h_log2fc", "10h", "infected", "bystander"),
            ("cp_enterocyte_48h_log2fc", "48h", "wt_", "control_")):
        table = pd.read_csv(_space("Hs", "infection", "GSE185247",
                                   f"GSE185247_ARG_{hours}_GEO_processed.txt.gz"),
                            sep="\t", index_col=0, low_memory=False)
        case = [c for c in table.columns if c.startswith(case_mark) or case_mark in c]
        control = [c for c in table.columns if c.startswith(control_mark) or control_mark in c]
        if not case or not control:
            raise ValueError(f"GSE185247 {hours}: no {case_mark}/{control_mark} columns")
        out[key] = by_symbol(log2_ratio(table, case, control, logged=True), symbols, universe)

    blood = pd.read_csv(_space("Hs", "infection", "GSE273483",
                               "GSE273483_Pvivax_RNAseq_counts.txt.gz"),
                        sep="\t", index_col=0, low_memory=False)
    control = [c for c in blood.columns if c.startswith("HC_")]
    for key, prefix in (("pv_blood_symptomatic_log2fc", "SM_"),
                        ("pv_blood_asymptomatic_log2fc", "AM_")):
        case = [c for c in blood.columns if c.startswith(prefix)]
        out[key] = by_symbol(log2_ratio(blood, case, control), symbols, universe)
    return out


def mouse_parasite(universe, symbols: dict) -> dict:
    """The ten mouse infection contrasts: four Toxoplasma tissues, macrophages, chronic infection,
    two Cryptosporidium experiments, cerebral malaria and the malaria liver stage."""
    out = {}

    path = _space("Mm", "infection", "GSE119855",
                  "GSE119855_Raw_counts_Toxoplasma_WT_IFNR_KO_4tissues.xlsx")
    counts = pd.read_excel(path, sheet_name=0).set_index("Ensembl_ID")
    titles = _series_titles(_space("Mm", "infection", "GSE119855", "GSE119855_series_matrix.txt.gz"))
    # "Lung_WT_Case [MS0691]" -> sample MS0691 is a wild-type infected lung.
    arms: dict = {}
    for title in titles:
        match = re.match(r"(\w+)_WT_(Case|Control) \[(MS\d+)\]", title)
        if match and match.group(3) in counts.columns:
            arms.setdefault((match.group(1), match.group(2)), []).append(match.group(3))
    for key, tissue in (("tg_lung_log2fc", "Lung"), ("tg_blood_log2fc", "Blood"),
                        ("tg_spleen_log2fc", "Spleen"), ("tg_liver_log2fc", "Liver")):
        out[key] = by_gene(log2_ratio(counts, arms[(tissue, "Case")], arms[(tissue, "Control")]),
                           universe)

    bmdm = pd.read_excel(_space("Mm", "infection", "GSE103113", "GSE103113_processed_data.xlsx"),
                         sheet_name="BMDMs transcripts")
    logged = bmdm.set_index("Name")[["Pru - log₂ total RPKM", "ui - log₂ total RPKM"]]
    out["tg_bmdm_log2fc"] = by_symbol(
        log2_ratio(logged, [logged.columns[0]], [logged.columns[1]], logged=True), symbols, universe)

    chronic = _tar_matrix(_space("Mm", "infection", "GSE117504", "GSE117504_RAW.tar"),
                          lambda n: (re.search(r"[IU][FM]_28d_\d+", n) or [None])[0],
                          "gene_id", "expected_count")
    chronic.index = pd.Series(chronic.index, dtype="string").str.split("_", n=1).str[0].to_numpy()
    out["tg_chronic_28d_log2fc"] = by_gene(
        log2_ratio(chronic, [c for c in chronic.columns if c[0] == "I"],
                   [c for c in chronic.columns if c[0] == "U"]), universe)

    crypto = pd.read_csv(_space("Mm", "infection", "GSE168680",
                                "GSE168680_Gullicksrud_GEO_processed.txt.gz"),
                         sep="\t", index_col=0, low_memory=False)
    for key, arm in (("cp_ileum_a_log2fc", "A_"), ("cp_ileum_b_log2fc", "B_")):
        case = [c for c in crypto.columns if c.startswith(arm + "inf")]
        control = [c for c in crypto.columns if c.startswith(arm + "uninf")]
        out[key] = by_symbol(log2_ratio(crypto, case, control, logged=True), symbols, universe)

    brain = pd.read_csv(_space("Mm", "infection", "GSE162535", "GSE162535_All.HTSeq.counts.txt.gz"),
                        sep="\t", index_col=0, low_memory=False)
    out["pb_brain_cerebral_malaria_log2fc"] = by_symbol(
        log2_ratio(brain, [c for c in brain.columns if c.startswith("MB_")],
                   [c for c in brain.columns if c.startswith("CB_")]), symbols, universe)

    liver = pd.read_csv(_space("Mm", "infection", "GSE44831", "GSE44831_matrix_BALB_C57BL.txt.gz"),
                        sep="\t", low_memory=False)
    out["py_liver_log2fc"] = by_symbol(
        pd.Series(pd.to_numeric(liver.log2_fold_change, errors="coerce").to_numpy(),
                  index=liver.gene_id.to_numpy()), symbols, universe)
    return out


def _series_titles(path: str) -> list:
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                return [p.strip().strip('"') for p in line.rstrip("\n").split("\t")[1:]]
    return []


# --------------------------------------------------------------------------- putting it together
BUILDERS = {"Hs": (human_tissue, human_parasite, HUMAN_SERIES),
            "Mm": (mouse_tissue, mouse_parasite, MOUSE_SERIES)}


def read_pack(code: str, out_dir: str) -> tuple:
    """The node table and the protein-gene link table out of the shipped pack, plus their names."""
    with zipfile.ZipFile(PACKS[code]) as z:
        z.extract(NODES[code], out_dir)
        z.extract(LINKS[code], out_dir)
        metadata = json.loads(z.read("column_metadata.json"))
        manifest = json.loads(z.read("manifest.json"))
    return (pd.read_parquet(os.path.join(out_dir, NODES[code])), metadata,
            manifest.get("column_licenses", {}))


def build(code: str, out_dir: str, probe: bool) -> dict:
    print(f"\n{code}")
    os.makedirs(out_dir, exist_ok=True)
    baseline, metadata, licences = read_pack(code, os.path.join(out_dir, "previous"))
    universe = pd.Index(baseline.gene_id.astype(str))
    print(f"  {len(universe):,} genes, {len(baseline.columns)} columns in the 2026.11.07 pack")
    symbols = symbol_map(code, set(universe))

    tissue_builder, parasite_builder, series = BUILDERS[code]
    tissue, tissue_records = tissue_builder(universe)
    parasite = parasite_builder(universe, symbols)

    columns = dict(tissue)
    columns.update(parasite)
    frame = baseline.set_index(baseline.gene_id.astype(str))
    for name, values in columns.items():
        if name in frame.columns:
            raise ValueError(f"{code}: {name} would overwrite a column the pack already carries")
        # float32: these are log expression and log ratios, where the seventh decimal place is not a
        # measurement, and the pack is a download. float64 made the human table 33 MB.
        #
        # An INFINITE value becomes missing, not a large number. GSE44831 is cuffdiff output and
        # writes `inf` and `-inf` for a fold change whose other arm was zero: the direction is known
        # and the size is not, which is exactly what "unknown" means here. Left in, it reached the
        # layout as a float32 infinity, UMAP refused the matrix, the harness fell back to PCA and
        # PCA failed on the same values -- three steps away from the one column that caused it.
        column = pd.to_numeric(values.reindex(universe), errors="coerce").to_numpy(dtype=np.float32)
        frame[name] = np.where(np.isfinite(column), column, np.nan).astype(np.float32)
    nodes = frame.copy().reset_index(drop=True)

    lenses = {
        "tissue": [c for c in (list(tissue) + [c for c in ("gtex_fibroblast_median_tpm",
                                                           "gtex_lcm_hepatocyte_median_tpm",
                                                           "bgee_top_score")
                                               if c in nodes.columns])],
        "parasite_facing": list(parasite),
    }
    for lens, features in lenses.items():
        measured = int(nodes[features].notna().any(axis=1).sum())
        print(f"  {lens}: {len(features)} features, {measured:,} of {len(nodes):,} genes with at "
              f"least one value ({100 * measured / len(nodes):.1f} per cent)")
        for name in features:
            have = int(nodes[name].notna().sum())
            if have == 0:
                raise ValueError(f"{code}/{lens}: {name} has no value for any gene")
            if lens == "parasite_facing":
                print(f"    {name}: {have:,} genes ({100 * have / len(nodes):.1f} per cent)")

    summary = dict(space=code, genes=len(nodes), baseline_columns=list(baseline.columns),
                   added_columns=list(columns), lenses={})
    if probe:
        for lens, features in lenses.items():
            summary["lenses"][lens] = dict(features=features, probe=True)
        return summary

    nodes.to_parquet(os.path.join(out_dir, NODES[code]), index=False)
    node_sha = _sha256(os.path.join(out_dir, NODES[code]))
    for lens, features in lenses.items():
        summary["lenses"][lens] = layout(code, lens, nodes, features, node_sha, out_dir)

    # Per-column provenance, licence and -- for every GEO-derived column -- the accession and the
    # paper, which is what `packs.build` now demands before it will carry one.
    for name, record in tissue_records.items():
        metadata[name] = record
    sources = {}
    for name, record in series.items():
        metadata[name] = dict(source=f"{record['accession']} ({record['parasite']})",
                              context=record["context"], accession=record["accession"],
                              license="G", citation=record["citation"])
        sources[name] = dict(accession=record["accession"], citation=record["citation"])
    for name, record in tissue_records.items():
        licences[name] = record["license"]
    for name in series:
        licences[name] = "G"
    summary.update(column_licenses=licences, column_sources=sources,
                   node_file_sha256=node_sha)
    with open(os.path.join(out_dir, f"{code.lower()}_column_metadata.json"), "w") as fh:
        json.dump(metadata, fh, indent=1, sort_keys=True)
    return summary


def layout(code: str, lens: str, nodes, features: list, node_sha: str, out_dir: str) -> dict:
    """One lens layout for one gene space, or a recorded refusal when too little is measured."""
    measured = int(nodes[features].notna().any(axis=1).sum())
    share = measured / len(nodes)
    if share < MINIMUM_MEASURED:
        print(f"  {code}/{lens}: REFUSED -- {100 * share:.1f} per cent of genes measured")
        return dict(lens=lens, refused=True, features=features, genes=len(nodes),
                    genes_with_any_value=measured,
                    reason=f"only {100 * share:.1f} per cent of genes carry any value")
    # NO `random_state`, on purpose, and it is not a loss of reproducibility we had.
    #
    # umap-learn overrides `n_jobs` to 1 whenever a seed is set, so a seeded run optimises the
    # layout on ONE core. The human tissue lens -- 58,988 genes over 106 features -- was still in
    # that single-threaded phase after two and a half hours on this machine. And the seed was never
    # buying the thing it costs: instruction 53 records that two rebuilds of the SAME lens from the
    # SAME table differ by up to 97 units even inside one process, so umap-learn here does not
    # honour the seed as a guarantee. Paying hours for a promise the library does not keep is the
    # worse of the two trades. What ties a layout to its data is the SHA-256 pair recorded in the
    # file, and the shipped file is the artefact.
    spec = EmbeddingSpec(name=f"{code}_gene_{lens}", blocks=(), extra_columns=tuple(features),
                         na_policy="median", max_missing=1.0, method="umap", n_components=3,
                         random_state=None)
    coords, used, rows, meta = project(nodes[["gene_id"] + features], spec,
                                       log=lambda m: print("   ", m), return_metadata=True)
    if not np.all(rows):
        raise ValueError(f"{code}/{lens}: a host gene layout must keep every gene")
    xyz = np.asarray(coords, dtype=np.float32)
    record = dict(meta)
    # `ordered_gene_ids` is every gene id again -- 58,988 of them -- and the file already carries
    # them as its own array. Kept, it made the summary 156,000 lines of identifiers.
    record.pop("ordered_gene_ids", None)
    record.update(space=code, lens=lens, features=list(used), requested_columns=features,
                  genes=len(nodes), genes_with_any_value=measured, node_file=NODES[code],
                  node_file_sha256=node_sha,
                  coordinates_sha256=hashlib.sha256(xyz.tobytes()).hexdigest())
    path = os.path.join(out_dir, f"{code.lower()}_gene_graph_{lens}.npz")
    np.savez_compressed(path, xyz=xyz, gene_ids=nodes.gene_id.to_numpy(dtype=str),
                        layout_metadata=np.array(json.dumps(record)))
    print(f"  {code}/{lens}: wrote {os.path.basename(path)}")
    return record


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--space", choices=sorted(BUILDERS))
    ap.add_argument("--probe", action="store_true", help="report coverage; build no layout")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "host_gene_lenses_2026_11_13"))
    args = ap.parse_args(argv)

    summary = {}
    for code in ([args.space] if args.space else sorted(BUILDERS)):
        summary[code] = build(code, args.out, args.probe)
    if not args.probe:
        with open(os.path.join(args.out, "summary.json"), "w") as fh:
            json.dump(summary, fh, indent=1, default=str)
    print(f"\n{args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
