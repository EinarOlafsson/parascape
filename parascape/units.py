"""The unit of a numeric column, with where that statement comes from (instruction 64.11).

A value scorecard reports MAE and RMSE in the unit of the column it predicts, so the card must say what that unit
is. Units here are curated per dataset from the source publication or the code that derives the column, and each
carries its basis; columns matched only by a name pattern say so; anything else is "unit not recorded" rather than
guessed.
"""
from __future__ import annotations

import re

_INVIVO = ("fit_invivo_PE", "fit_invivo_lung", "fit_invivo_liver", "fit_invivo_spleen", "fit_invivo_heart",
           "fit_invivo_brain")
_PF_STAGES = ("expr_ring", "expr_early_trophozoite", "expr_late_trophozoite", "expr_schizont",
              "expr_gametocyte_ii", "expr_gametocyte_v", "expr_ookinete")

#: column -> (unit, basis). Curated per dataset; the basis names the source that states it.
CURATED = {
    "fit_invitro_hff": ("phenotype score: log2 fold change of guide abundance after growth in fibroblasts "
                        "(negative = the parasite grows worse without the gene)",
                        "Sidik et al. 2016, PMID 27594426 (dataset crispr_invitro)"),
    **{c: ("composite score: differential fitness multiplied by its significance, dimensionless "
           "(negative = depleted in the tissue)",
           "Giuliano et al. 2024, PMID 38977907 (dataset crispr_invivo_composite)") for c in _INVIVO},
    # The rna108740_* columns themselves are ln(FPKM + 1) of GEO's GSE108740_FPKM.xlsx (checked 2026-11-07); the
    # summaries below undo that log and take log2 once (instruction 64.11a).
    "expr_tachy": ("log2(FPKM + 1) of the mean tachyzoite FPKM",
                   "GSE108740 (GEO, FPKM); computed in build_graph.py, corrected 2026-11-07 (64.11a)"),
    "expr_cyst": ("log2(FPKM + 1) of the mean tissue-cyst FPKM",
                  "GSE108740 (GEO, FPKM); computed in build_graph.py, corrected 2026-11-07 (64.11a)"),
    "expr_max": ("log2(FPKM + 1) of the highest GSE108740 sample",
                 "GSE108740 (GEO, FPKM); computed in build_graph.py, corrected 2026-11-07 (64.11a)"),
    **{c: ("ln(FPKM + 1), despite the _FPKM in the name", "GSE108740 (GEO FPKM table), logged upstream; checked 2026-11-07")
       for c in ("rna108740_Tachyzoites_T2_FPKM", "rna108740_Tachyzoites_T4_FPKM", "rna108740_Tissue_cysts_A_FPKM",
                 "rna108740_Tissue_cysts_B_FPKM", "rna108740_Day3_CS4_FPKM", "rna108740_Day3_CS5_FPKM",
                 "rna108740_Day3_CS6_FPKM", "rna108740_Day5_CS2_FPKM", "rna108740_Day5_K4_14_FPKM",
                 "rna108740_Day7_CS7_FPKM", "rna108740_Day7_CS8_FPKM", "rna108740_Day7_CS9_FPKM")},
    "piggybac_mis": ("mutagenesis index score, 0 to 1 (low = the gene resists disruption, likely essential)",
                     "Zhang et al. 2018, PMID 29724925 (dataset plasmodb_pf3d7_piggybac)"),
    "piggybac_mfs": ("mutant fitness score (low = fitness cost)",
                     "Zhang et al. 2018, PMID 29724925 (dataset plasmodb_pf3d7_piggybac)"),
    **{c: ("transcript abundance as served by PlasmoDB for the Su seven-stage RNA-seq; the report does not state "
           "the normalisation", "PlasmoDB expression report (dataset plasmodb_pf3d7_expression)")
       for c in _PF_STAGES},
    "mean_plddt": ("pLDDT, 0 to 100: AlphaFold confidence averaged over residues",
                   "AlphaFold DB (Varadi et al. 2024; Jumper et al. 2021)"),
}

#: (pattern, unit) for columns not curated one by one; matched case-insensitively on the name.
PATTERNS = (
    (r"_tpm$|^tpm_|_tpm_", "TPM (transcripts per million)"),
    (r"_fpkm$|_fpkm_", "FPKM (fragments per kilobase per million)"),
    (r"lfc$|_lfc_|log2fc|_log2$|log2_fc", "log2 fold change"),
    (r"plddt", "pLDDT, 0 to 100"),
    (r"_percentile$", "percentile among genes, 0 to 100"),
)

UNRECORDED = "unit not recorded"


def unit_of(column: str) -> tuple:
    """(unit, basis) for a column: curated first, then a name pattern (said so), else unit not recorded."""
    if column in CURATED:
        return CURATED[column]
    for pattern, unit in PATTERNS:
        if re.search(pattern, str(column), flags=re.IGNORECASE):
            return unit, "read from the column name, not curated"
    return UNRECORDED, ""


def sentence(column: str) -> str:
    """One line for a card: what MAE and RMSE are measured in."""
    unit, basis = unit_of(column)
    if unit == UNRECORDED:
        return f"MAE and RMSE are in the units of {column}, which are not recorded."
    return f"MAE and RMSE are in the units of {column}: {unit} ({basis})."
