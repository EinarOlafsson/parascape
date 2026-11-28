"""The label catalogue: every label the Infer benchmark scores, where it comes from and how it was made.

Einar, 2026-11-25: every label carries metadata saying where its column comes from, how it was
treated and generated, and what classes it contains; and a label NAME says what kind of evidence
it is when that is not obvious (observed in Toxoplasma, or copied from another organism).

`ENTRIES` is written by hand, once per label: the evidence type, the source datasets (resolved
through `parascape.datasets`, so citations are never typed twice), the steps that turn the source
into the label, and what the benchmark leaves out of the features to avoid circularity.
`catalogue(L, exclusions)` adds what can only be measured on the data: genes labelled, the classes
and their sizes, the columns actually left out. Every benchmark run writes it next to its results.

Evidence types:

- **observed** - measured in Toxoplasma (an experiment, or a curator reading one);
- **copied by orthology** - measured in another organism and carried over through orthologs: a
  claim about the ortholog, not about the Toxoplasma protein; never mixed with observed labels;
- **read from papers** - experimental findings read from Toxoplasma papers by the literature run.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

OBSERVED = "observed"
COPIED = "copied by orthology"
PAPERS = "read from papers"
EVIDENCE = (OBSERVED, COPIED, PAPERS)

ENTRIES: dict = {
    "effect on host": dict(
        evidence=OBSERVED, sources=("pvm_proximity", "hosttx_effectors", "organelle_surface_turboid"),
        column="information map trait `host_facing`",
        made=("`host_facing` of the information map (a protein, or part of it, found beyond the parasite "
              "plasma membrane): yes for genes with that evidence, 'other' for every other gene "
              "(positive-unlabelled). The first version used hyperLOPIT internal markers as negatives and "
              "copied the niche (Cramer V 0.99), so it was dropped."),
        left_out=("the trait's own assays (PVM proximity, host-transcription effectors, organelle-surface "
                  "proximity); none of their columns is in any view, so nothing has to be dropped")),
    "enzymatic function (EC class)": dict(
        evidence=OBSERVED, sources=("toxodb_ec_numbers",), column="ec_number",
        made=("First digit of the ToxoDB EC number (1 oxidoreductase ... 7 translocase). ToxoDB EC numbers "
              "mix curated and sequence-based assignments; the evidence code is not in the column."),
        left_out="ec_number"),
    "hyperLOPIT niche": dict(
        evidence=OBSERVED, sources=("lopit_tgon",), column="information map trait `compartment_hyperlopit`",
        made="Curated hyperLOPIT marker niche of Barylyuk 2020; non-markers unlabelled.",
        left_out="the whole place view and every lopit column"),
    "phenotype score, fibroblasts": dict(
        evidence=OBSERVED, sources=("crispr_invitro",), column="fit_invitro_hff",
        made="Fitness score of the genome-wide CRISPR screen in human foreskin fibroblasts (continuous).",
        left_out="fit_invitro_hff; NOTE the other fibroblast screens stay in, so this is close to a "
                 "repeat measurement, not an inference"),
    "phenotype score, IFN-gamma macrophages": dict(
        evidence=OBSERVED, sources=("crispr_macrophage",), column="fit_ifng",
        made="Wang 2020 IFN-gamma macrophage screen; sign flipped so negative = loss costs growth.",
        left_out="every column of the same screen"),
    "phenotype score, in vivo (brain)": dict(
        evidence=OBSERVED, sources=("crispr_invivo_composite",), column="fit_invivo_brain",
        made="Giuliano 2024 in vivo composite score, brain (continuous).",
        left_out="every column of the same screen"),
    "cell-cycle phase": dict(
        evidence=OBSERVED, sources=("xue_singlecell",), column="cellcycle_phase",
        made="Peak cell-cycle phase from the single-parasite atlas of Xue 2020.",
        left_out="every cellcycle column"),
    "nutrient dependence (glucose)": dict(
        evidence=OBSERVED, sources=("crispr_glucose_limitation",), column="fit_glucose_dependence",
        made="Fitness without glucose relative to complete medium (continuous).",
        left_out="every column of the same screen"),
    "imaging phenotype": dict(
        evidence=OBSERVED, sources=("splitcas9_imaging_screen",), column="screen_any_phenotype",
        made=("Any visible defect (replication, apicoplast, F-actin or egress) in the arrayed splitCas9 "
              "imaging screen, scored by eye; 320 genes screened, 99 with a defect; unscreened genes "
              "unlabelled, not 'normal'."),
        left_out="screen_any_phenotype"),
    "KEGG pathway (copied by orthology)": dict(
        evidence=COPIED, sources=("kegg_tgo",), column="parascape/data/interaction_labels/kegg_pathway_label.tsv",
        made=("KEGG T. gondii ME49 (tgo) gene-pathway links, NCBI gene ids mapped to TGME49 through NCBI E-utilities "
              "(857 of 881 genes); global and overview maps (tgo011xx/012xx) left out; each gene takes its most "
              "specific pathway with >= 10 genes (763 genes, 35 pathways)."),
        left_out="nothing (no feature is built from KEGG)"),
    "hyperLOPIT niche (copied by orthology)": dict(
        evidence=COPIED, sources=("ortholopit",), column="ortholopit_label",
        made="The hyperLOPIT niche of the gene's orthologs in other apicomplexan hyperLOPIT maps (localization.py).",
        left_out="the whole place view and every lopit column"),
    "compartment (papers)": dict(
        evidence=PAPERS, sources=("literature_affinage",), column="lit_compartment",
        made=("Compartment found experimentally in single-gene papers (place findings only); genes with "
              "several compartments are unlabelled; classes under 10 genes dropped."),
        left_out="the whole place view and every lopit column"),
    "host-facing (papers)": dict(
        evidence=PAPERS, sources=("literature_affinage",), column="lit_host_facing",
        made="Host-facing by a paper finding; 'not host-facing' only for genes some paper studied.",
        left_out="nothing (no feature is read from papers)"),
    "interaction module (STRING, measured in Toxoplasma)": dict(
        evidence=OBSERVED, sources=("string_tgon_measured",), column="node column (interaction_labels.COLUMNS)",
        made=("STRING v12.0 T. gondii VEG pairs with Toxoplasma's own experimental or curated-database "
              "evidence (channels `experiments`, `database`; nothing transferred), joined to ME49 by locus "
              "number; Louvain communities (resolution 1, seed 0); communities under 10 genes dropped; "
              "named by their best-connected gene."),
        left_out="nothing (no feature is built from STRING)"),
    "interaction module (STRING, copied by orthology)": dict(
        evidence=COPIED, sources=("string_tgon_transferred",), column="node column (interaction_labels.COLUMNS)",
        made=("STRING v12.0 T. gondii VEG pairs whose `experiments_transferred` or `database_transferred` "
              "score is at least 700, i.e. evidence measured in other organisms; co-expression and text "
              "mining NOT used (co-expression is circular with the transcript view). Louvain communities "
              "as above."),
        left_out="nothing; families (orthogroups) keep paralogs in one fold"),
    "complex (CORUM human, copied by orthology)": dict(
        evidence=COPIED, sources=("corum_human_by_orthology",), column="node column (interaction_labels.COLUMNS)",
        made=("CORUM human complexes; subunit UniProt -> Ensembl gene (HGNC) -> OrthoMCL 6.21 orthogroup -> "
              "Toxoplasma genes of that orthogroup (orthogroups with more than 3 Toxoplasma genes not used); "
              "a complex needs 2 mapped orthogroups; each gene takes the most specific complex that keeps "
              "10 Toxoplasma genes. Overlapping CORUM entries (spliceosome A/B/C) stay separate classes."),
        left_out="nothing; families (orthogroups) keep paralogs in one fold"),
}

# Sources that are not datasets of the node table are described here (none at present: STRING and CORUM
# joined the dataset registry for release 0.0.14.0).
EXTRA_SOURCES: dict = {
    "kegg_tgo": dict(name="KEGG PATHWAY, Toxoplasma gondii ME49 (tgo)", citation="Kanehisa M et al. Nucleic Acids Res "
                     "2023;51:D587-D592", pmid="36300620", url="https://rest.kegg.jp/link/pathway/tgo"),
    "ortholopit": dict(name="hyperLOPIT niches transferred from orthologs (ParaScape localization.py)",
                       citation="", pmid="", url=""),
}


def source_record(key: str) -> dict:
    """Name, citation, PubMed id and URL of one source, from the dataset registry or EXTRA_SOURCES."""
    if key in EXTRA_SOURCES:
        return {"key": key, **EXTRA_SOURCES[key]}
    from parascape import datasets
    try:
        d = datasets.get(key)
    except Exception:  # noqa: BLE001 - an unknown key is reported, not fatal
        return {"key": key, "name": "UNKNOWN SOURCE", "citation": "", "pmid": "", "url": ""}
    return {"key": key, "name": d.name, "citation": d.citation or "", "pmid": d.pmid or "", "url": d.url or ""}


def catalogue(L: dict, excluded: dict | None = None) -> pd.DataFrame:
    """One row per label: its metadata plus what is measured on the data (genes, classes, exclusions).

    `L` is the benchmark label dict (name -> (values, kind, target columns)); `excluded` maps a label
    to the feature columns actually left out for it. A label without an entry raises: a label nobody
    has described must not reach a result.
    """
    rows = []
    for name, (y, kind, _t) in L.items():
        if name not in ENTRIES:
            raise KeyError(f"label {name!r} has no catalogue entry in parascape/label_catalogue.py")
        e = ENTRIES[name]
        srcs = [source_record(k) for k in e["sources"]]
        lab = pd.notna(y)
        if kind == "continuous":
            v = np.asarray(y, float)[lab]
            classes = f"continuous: median {np.median(v):.3g}, range {v.min():.3g} to {v.max():.3g}" if len(v) else ""
        else:
            counts = pd.Series(np.asarray(y, object)[lab]).value_counts()
            classes = "; ".join(f"{c} ({n})" for c, n in counts.items())
        ex = sorted(excluded.get(name, [])) if excluded else []
        rows.append({"label": name, "evidence": e["evidence"], "kind": kind,
                     "genes_labelled": int(lab.sum()), "classes": classes,
                     "source_column": e["column"], "how_made": e["made"],
                     "left_out_rule": e["left_out"], "left_out_columns": len(ex),
                     "left_out_list": ", ".join(ex),
                     "sources": "; ".join(f"{s['name']} ({s['citation'] or 'citation to confirm'}"
                                          f"{', PMID ' + s['pmid'] if s['pmid'] else ''})" for s in srcs),
                     "source_urls": " ".join(s["url"] for s in srcs if s["url"])})
    return pd.DataFrame(rows)


def write(L: dict, out_dir: str, excluded: dict | None = None) -> pd.DataFrame:
    """Write `label_catalogue.tsv` and `.json` into a run folder and return the table."""
    t = catalogue(L, excluded)
    t.to_csv(os.path.join(out_dir, "label_catalogue.tsv"), sep="\t", index=False)
    json.dump(t.to_dict("records"), open(os.path.join(out_dir, "label_catalogue.json"), "w"), indent=1)
    return t
