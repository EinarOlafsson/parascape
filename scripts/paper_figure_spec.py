"""The twelve ParaScape figures, panel by panel: the contract the figure build reads.

Einar, 2026-11-19, on the figures: build the figures agreed in instruction 101, with the
spacr-graphs skill, at least 6 and preferably 12 panels per figure, and with the bradyzoite
single-cell analysis present. Panels had been drawn one at a time and drifted from the plan.

So the specification is DATA THE BUILD CONSUMES. Every panel names its graph type (from the
skill) and its data source. A panel whose measurement does not exist yet has source=None and is
drawn as PLACEHOLDER DATA IN ITS REAL FORM, marked in the legend, so the figure structure can be
judged before the data exists.

Rules every figure obeys (instruction 101 and the skill): titled by its conclusion and ending in
a sentence about the parasite; every panel the same shape of claim; never organised by data type.
"""
from __future__ import annotations

from dataclasses import dataclass

BRADY = "/media/carruthers/mnt3/claude/brady_metabolic_classes"
BRADY_H5AD = f"{BRADY}/data/processed/v150f_clustered.h5ad"
BRADY_RES = f"{BRADY}/results_v150f"
AF3_REPORTS = "/media/carruthers/mnt3/claude/toxoplasma_projects/af3_folding/reports"


@dataclass(frozen=True)
class Panel:
    letter: str
    what: str
    form: str
    source: str | None
    note: str = ""


@dataclass(frozen=True)
class Figure:
    number: int
    conclusion: str
    panels: tuple


FIGURES = (
    Figure(1, "ParaScape gathers 177 datasets into one gene space; for most Toxoplasma genes, "
              "host response and perturbation phenotype have never been measured.", (
        Panel("a", "what a user brings, the gene space, and the questions it answers",
              "schematic", None, "redesign for information content, not layout"),
        Panel("b", "datasets per organism", "bar", "corpus_inventory_*/datasets.tsv"),
        Panel("c", "datasets per level of biology", "bar", "corpus_inventory_*/datasets.tsv"),
        Panel("d", "datasets by year of publication", "line", None),
        Panel("e", "genes with any measurement, per category (categories defined in legend)",
              "bar", "corpus_inventory_*/coverage.tsv"),
        Panel("f", "categories measured per gene", "violin",
              "measured_absence_*/categories_per_gene.tsv"),
    )),
    Figure(2, "Cell cycle is the dominant confounder of the published Toxoplasma "
              "transcriptome: removing it brings disagreeing studies into agreement.", (
        Panel("a", "control-to-control agreement between studies", "heatmap", None),
        Panel("b", "agreement before vs after removing cell cycle, per study", "scatter",
              "control_comparability_*/comparability.tsv"),
        Panel("c", "gain from cell cycle vs random axes of equal rank", "violin", None),
        Panel("d", "agreement gain per study", "bar", None),
        Panel("e", "contrasts made comparable, by ladder depth", "bar", None),
        Panel("f", "cell-cycle phase composition of each study's controls", "stacked bar", None),
    )),
    Figure(3, "Co-expression modules are, for the most part, single organelles built on a "
              "shared schedule.", (
        Panel("a", "named modules: share in one compartment vs matched baseline", "bar",
              "coexpression_modules_*/modules.tsv"),
        Panel("b", "every module against its own baseline", "scatter",
              "coexpression_modules_*/modules.tsv"),
        Panel("c", "module sizes", "box", "coexpression_modules_*/modules.tsv"),
        Panel("d", "distance above the detection limit, per module", "bar",
              "coexpression_modules_*/modules.tsv", "split at 12 conditions per panel"),
        Panel("e", "module expression across the cell cycle", "line", None),
        Panel("f", "RNA module vs protein module agreement", "scatter", None),
    )),
    Figure(4, "Treatment responses group by laboratory before biology.", (
        Panel("a", "contrast similarity: same study, same kind, other kind", "box",
              "what_structures_contrasts_*/pairs.tsv", "y label: Spearman (rho)"),
        Panel("b", "contrasts clustered by response", "heatmap", None),
        Panel("c", "variance apportioned to lab, strain, medium, platform, perturbation",
              "stacked bar", None),
        Panel("d", "within-lab vs across-lab agreement per perturbation class", "box", None),
        Panel("e", "one perturbation in two laboratories", "scatter", None),
        Panel("f", "conditions grouping once lab is removed", "heatmap", None),
    )),
    Figure(5, "A tissue cyst holds several bradyzoite states, and published treatments reach "
              "them unequally: one subpopulation is reached by almost nothing.", (
        Panel("a", "15,177 bradyzoites coloured by region", "umap", BRADY_H5AD),
        Panel("b", "the same cells coloured by published subtype A-E", "umap", BRADY_H5AD),
        Panel("c", "cell-cycle S/M score across the cyst", "umap", BRADY_H5AD),
        Panel("d", "marker genes per region", "dot plot", f"{BRADY_RES}/cluster_markers.tsv"),
        Panel("e", "published contrasts reaching each subtype", "bar",
              "bradyzoite_reach_*/contrasts.tsv"),
        Panel("f", "q of every contrast against every subtype", "violin",
              "bradyzoite_reach_*/contrasts.tsv"),
        Panel("g", "treatment by region", "heatmap", f"{BRADY_RES}/drug_by_cluster.tsv"),
        Panel("h", "drug-target expression by region (CDPK1, PKG, DHFR-TS)", "violin",
              BRADY_H5AD),
        Panel("i", "druggability by region", "bar", f"{BRADY_RES}/druggability.tsv"),
        Panel("j", "in vivo fitness cost by region", "bar",
              f"{BRADY_RES}/gse229548_invivo_screen_by_region.tsv"),
        Panel("k", "cyst-wall genes by region", "bar",
              f"{BRADY_RES}/cyst_wall_and_bradyzoite_sources.tsv"),
        Panel("l", "differentiation time course mapped onto the cyst", "bar",
              f"{BRADY_RES}/differentiation_timecourse_mapping.tsv"),
    )),
    Figure(6, "Organelles built in waves during division can be placed from transcription "
              "alone; those made continuously cannot.", (
        Panel("a", "compartment placed from RNA alone vs all evidence", "scatter",
              "coexpression_components_localization_*/scores.tsv"),
        Panel("b", "protein abundance explained by each layer", "bar",
              "what_sets_protein_abundance_*/layers.tsv"),
        Panel("c", "which measurement explains which", "bar",
              "which_data_explains_which_*/pairs.tsv"),
        Panel("d", "translation efficiency vs protein abundance", "scatter", None),
        Panel("e", "mRNA decay vs protein abundance", "scatter", None),
        Panel("f", "modification vs protein abundance", "scatter", None),
    )),
    Figure(7, "A graph whose nodes carry traits and whose edges carry derived relationships "
              "places proteins that no single relationship can.", (
        Panel("a", "nodes with traits, edges with relationships", "schematic", None),
        Panel("b", "genes answered vs accuracy, per edge layer", "scatter",
              "network_layers_*/layers.tsv"),
        Panel("c", "accuracy per edge layer", "bar", "network_layers_*/layers.tsv"),
        Panel("d", "genes answered per edge layer", "bar", "network_layers_*/layers.tsv"),
        Panel("e", "union vs best single layer, per compartment", "scatter", None),
        Panel("f", "overlap between edge layers", "heatmap", None),
    )),
    Figure(8, "Secreted effectors vary between strains far more than housekeeping genes, and "
              "AlphaFold cannot model the interfaces of the most disordered of them.", (
        Panel("a", "nonsynonymous / synonymous SNPs per gene, by class", "violin",
              "strain_variation_by_class_*/per_gene.tsv"),
        Panel("b", "fraction of residues with pLDDT >= 70, host-facing genes", "violin",
              f"{AF3_REPORTS}/bait_folded_fraction.json"),
        Panel("c", "strain variation vs ordered fraction, per effector", "scatter", None),
        Panel("d", "strain variation by secretory class", "box", None),
        Panel("e", "SNP density along an effector", "line", None),
        Panel("f", "the most variable effectors", "bar", None),
    )),
    Figure(9, "Where a protein is transfers between apicomplexan species; how its message is "
              "written does not.", (
        Panel("a", "chance-corrected agreement between orthologues, per layer", "bar",
              "species_transfer_*/layers.tsv"),
        Panel("b", "localisation transfer per compartment", "bar", None),
        Panel("c", "fitness: Toxoplasma vs Plasmodium orthologues", "scatter", None),
        Panel("d", "codon usage: Toxoplasma vs Plasmodium orthologues", "scatter", None),
        Panel("e", "transfer vs orthologue sequence identity", "scatter", None),
        Panel("f", "Toxoplasma localisation predicting Cryptosporidium", "bar", None),
    )),
    Figure(10, "Most of the host response is shared across effectors, and the parasite itself "
               "changes in the same cell for almost none of them.", (
        Panel("a", "host genes moved per effector: own vs shared", "stacked bar",
              "host_response_shared_*/per_effector.tsv"),
        Panel("b", "parasite vs host evidence per knockout, single cell", "scatter",
              "dual_perturb_single_cell_*/*_targets.tsv"),
        Panel("c", "host genes shared by number of effectors", "bar",
              "host_response_shared_*/effectors_per_host_gene.tsv"),
        Panel("d", "tissue-specific vs shared host response", "bar", None),
        Panel("e", "host-parasite co-fold coverage and confidence", "scatter", None,
              "report coverage, state plainly it is uncontrolled"),
        Panel("f", "cells per knockout in the dual perturb-seq", "bar",
              "dual_perturb_single_cell_*/*_targets.tsv"),
    )),
    Figure(11, "The map fills in localisation and cell-cycle phase at a validated rate, and "
               "phenotype barely above guessing.", (
        Panel("a", "accuracy above baseline, best strategy per target", "bar",
              "filling_the_map_*/best_per_target.tsv"),
        Panel("b", "every strategy per target", "box",
              "filling_the_map_*/by_target_and_strategy.tsv"),
        Panel("c", "accuracy vs answer rate", "scatter", "filling_the_map_*/best_per_target.tsv"),
        Panel("d", "calibration of predicted probabilities", "line", None),
        Panel("e", "new assignments per target at the validated rate", "bar", None),
        Panel("f", "held-out accuracy by evidence category removed", "bar", None),
    )),
    Figure(12, "The most valuable new dataset for the field is a host-response screen: almost "
               "no Toxoplasma gene has been measured that way.", (
        Panel("a", "smallest detectable effect per dataset", "box", None),
        Panel("b", "genes measured per category, with detection limits", "bar", None),
        Panel("c", "absence as a claim: genes with no effect above the detection limit", "bar",
              None),
        Panel("d", "the least-measured genes and what is known of them", "dot plot", None),
        Panel("e", "genes gained per category by one new dataset of each kind", "bar", None),
        Panel("f", "conditions no dataset has tested", "heatmap", None),
    )),
)


def summary() -> str:
    lines = []
    for f in FIGURES:
        real = sum(1 for p in f.panels if p.source)
        lines.append(f"Figure {f.number:2d}: {len(f.panels):2d} panels "
                     f"({real} real, {len(f.panels) - real} placeholder)")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
    total = sum(len(f.panels) for f in FIGURES)
    real = sum(1 for f in FIGURES for p in f.panels if p.source)
    print(f"\n{total} panels, {real} from data, {total - real} placeholder")
