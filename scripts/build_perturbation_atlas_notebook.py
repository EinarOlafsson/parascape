"""The perturbation atlas recovery benchmark (instruction 71, product 2), as an executed notebook.

    python scripts/build_perturbation_atlas_notebook.py
"""
from __future__ import annotations

import datetime as _dt
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_runner import ExecutedNotebook  # noqa: E402

DATE = "2026_11_07"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "results", f"perturbation_atlas_{DATE}")
NB = os.path.join(REPO, "notebooks", f"perturbation_atlas_{DATE}.ipynb")

G = "expr_lfc_"
GROUPS = {
    "BFD1 stabilised by Shield-1": [
        G + "GSE139537_condition_shield1_vs_etoh_genotype_delta_bfd1_dd_bfd1_ty",
        G + "GSE223877_treatment_standard_medium_3um_shield_1_vs_standard_medium_ve",
        G + "GSE223877_treatment_standard_medium_3um_shield_1_vs_standard_medium_ve_"],
    "BFD1 lost under stress": [
        G + "GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_24_hr",
        G + "GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_48_hr",
        G + "GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_72_hr",
        G + "GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_200385",
        G + "GSE218921_genotype_bfd1_ko_vs_wild_type_treatment_high_ph_low_c02_low"],
    "Alkaline stress against normal medium, wild type": [
        G + "GSE132238_condition_stress_vs_unstress_genotype_wt_time_24_hr",
        G + "GSE132238_condition_stress_vs_unstress_time_48_hr_genotype_wt",
        G + "GSE132238_condition_stress_vs_unstress_time_72_hr_genotype_wt",
        G + "GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me",
        G + "GSE223869_treatment_alkaline_stress_media_ambient_co2_vs_standard_cult"],
    "AP2XII-1 knockdown": [
        G + "GSE224934_iKD_AP2XII_1_IAA_12h_biological_repeat",
        G + "GSE224934_iKD_AP2XII_1_IAA_24h_biological_repeat",
        G + "GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_1",
        G + "GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_2",
        G + "GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_3"],
}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    sys.path.insert(0, REPO)
    nb = ExecutedNotebook("Which known perturbation does my mutant look like?",
                          namespace={"REPO": REPO, "OUT": OUT, "GROUPS": GROUPS})
    nb.md("**What this is.** The perturbation atlas of instruction 71: every within-study contrast of the "
          "uniformly reprocessed RNA-seq is a signature (a shrunken log2 fold change against that study's own "
          "control), and a query is ranked against all of them by rank correlation.",
          "**The test.** Four perturbations measured by independent studies, curated by reading the contrast "
          "labels. Each member is queried with its own study held out entirely; it passes when another "
          "lab's experiment of the same perturbation comes first. Chance is the expected best rank of the "
          "partners among all candidates. AP2XII-1 knockdown is two deposits of one system (the same lab), "
          "so it is the weakest test of independence and is reported as such.",
          f"Built {_dt.date.today().isoformat()}.")
    nb.code("import os, pandas as pd",
            "from parascape import perturbation_atlas as P",
            "store = P.load()",
            "print(store['release'], store['lfc'].shape[1], 'signatures over', store['registry'].study.nunique(), 'studies')",
            "missing = {k: [m for m in v if m not in store['registry'].index] for k, v in GROUPS.items()}",
            "print('members not in the store:', {k: v for k, v in missing.items() if v})")
    nb.code("rec = P.recovery(store, GROUPS)",
            "rec.to_csv(os.path.join(OUT, 'recovery.tsv'), sep='\\t', index=False)",
            "summary = rec.groupby('perturbation').agg(queries=('query', 'size'), top1=('top1', 'mean'),",
            "    top5=('top5', 'mean'), median_best_rank=('best_partner_rank', 'median'),",
            "    median_chance_rank=('chance_rank', 'median')).round(2)",
            "summary.to_csv(os.path.join(OUT, 'summary.tsv'), sep='\\t')",
            "summary")
    nb.code("rec[['perturbation', 'study', 'best_partner_rank', 'chance_rank', 'candidates', 'best_partner']]")
    nb.md("### An example query, as a user would see it",
          "The BFD1 knockout under stress from GSE218921, its own study held out: the top ten matches.")
    nb.code("q = GROUPS['BFD1 lost under stress'][-1]",
            "hits = P.query(store['lfc'][q], store, exclude_study='GSE218921').head(10)",
            "hits.to_csv(os.path.join(OUT, 'example_query.tsv'), sep='\\t', index=False)",
            "hits[['rank', 'study', 'contrast', 'spearman', 'genes']]")
    nb.write(NB)


if __name__ == "__main__":
    main()
