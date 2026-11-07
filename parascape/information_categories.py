"""Information categories: what KIND of knowledge each data column embodies, so inference cannot be circular.

Instruction 73, "Information categories: the mechanism that stops circular inference" (Einar,
2026-11-05). The rule, verbatim in effect: **every column in the map carries an information
category. Inferring a label of category X excludes every column of category X by default, whatever
its name, family or study.** Including one is possible but needs an explicit reason that is recorded
in the run and shown on the scorecard, so a circular run is visible rather than impossible.

Why a per-column category rather than the per-question `evidence_forbidden` list this generalises:
that list already failed once. On the harness's first real run the family importance exposed three
further pooled CRISPR fitness screens whose family names contain no "crispr" -- `gra17_synthlethal`
was top of the table, with `invivo_platform` and `gra12` behind it -- each answering a fitness
question with the same kind of measurement as its own labels. A human reading an importance table
caught it; a category catches it by construction. The same hole is wider for the UMAP strategies
Einar asked about: an embedding built from "all the data" and scored against held-out localisation
labels is circular if any input column carries `place`, and the circularity is in the *input
selection*, not in the model, so nothing downstream can see it.

**Categories are about the kind of knowledge, not the assay.** Two columns from completely different
instruments share a category when they tell you the same kind of thing (hyperLOPIT and an
ortholog-transferred TAGM label are both `place`), and two columns from one deposit do not when they
do not (the RNA arm of a Ribo-seq study is `transcript`, its RPF and TE arms are `translation`).

**Annotation is derived from provenance, not from names.** `datasets.provenance(column, organism)`
gives each node column its source dataset, `level` and `kind` (transcription/RNAseq,
DNA/CRISPR_screen, post_translation/LOPIT...), and that pair is the primary rule
(`LEVEL_KIND_CATEGORIES`). A dataset whose level and kind are true but not decisive gets a record of
its own (`DATASET_CATEGORIES`) justified from the registry entry -- a pooled CRISPR screen read out
as host transcription is `host_response`, not `fitness`. A single column that differs from its
dataset gets an exact override (`COLUMN_CATEGORIES`). Only columns with no provenance at all --
the trait tables built in `results/`, and fifteen computed node columns -- are resolved by name
pattern, and those are flagged as judged so a reviewer can argue with them. `derivation(column)`
reports which of those routes was used for any column, and `derivations()` is the whole reviewable
table.

**Nothing is left silently unannotated.** `unannotated()` lists every shipped column with no
category at all; it is empty, and `tests/test_information_categories.py` is the gate that keeps it
empty. A column whose category genuinely cannot be told is to be left out of the mapping so that the
gate fails loudly, which is a great deal cheaper than a wrong category passing quietly.

    categories_of("fit_no_glucose")        # ('fitness',)
    categories_of("pvm_proximity_positive")  # ('place', 'interaction')
    excluded_for("place")                  # every column that carries localisation knowledge
    unannotated()                          # () -- the gate
"""
from __future__ import annotations

import re

import pandas as pd

from . import datasets, organisms

#: The seventeen categories of instruction 73, each defined by the KIND of knowledge it is rather than
#: by the instrument that produced it. The vocabulary is fixed: `tests/test_information_categories`
#: checks it against the instruction, because a category invented in passing is a category nobody
#: excludes on.
CATEGORIES = {
    "place": "Where in the cell, or at which interface, the protein is -- localisation knowledge.",
    "fitness": "How much the parasite needs the gene to survive or grow in some condition.",
    "transcript": "How much RNA of the gene there is, in some state or contrast.",
    "translation": "How much of that RNA is being translated -- ribosome occupancy and efficiency.",
    "protein_abundance": "How much of the protein there is.",
    "modification": "Which covalent modifications the molecule carries (PTM sites, m6A).",
    "interaction": "Which other molecules it is physically with.",
    "sequence": "What the primary sequence is and what can be read directly off it.",
    "structure": "What the folded shape is, and how confidently it is known.",
    "domain": "Which known families and domains the protein contains.",
    "conservation": "How the gene is distributed and constrained across species and strains.",
    "chromatin": "The state of the locus itself: accessibility, histone marks, factor occupancy.",
    "phenotype": "What is observed when the gene is perturbed -- imaging, drug sensitivity, differentiation.",
    "function": "What the protein is held to DO: GO terms, EC numbers, coarse functional class.",
    "host_response": "How the host cell or host immune system responds to the parasite or protein.",
    "attention": "Bookkeeping: publication counts, study effort -- a predictor of nothing biological.",
    "identity": "Names and accessions -- never a feature.",
}

#: Categories whose columns are never offered as model features at all, whatever the question. An
#: accession is not evidence, and a gene name that contains "rhoptry" is the answer written into the
#: question. `identity` is the only one: everything else is excluded per question, not globally.
NEVER_A_FEATURE = ("identity",)

# --------------------------------------------------------------------------- rule 4: level and kind
#: The primary rule: `(level, kind)` from `datasets.provenance` -> categories. The on-disk level
#: already states what kind of measurement a dataset is (`datasets` module docstring), so this is
#: provenance doing the work rather than a name heuristic. Deliberately incomplete where a pair is
#: genuinely ambiguous -- `post_translation/proteomics` carries PTM searches, melting curves,
#: proximity labelling and fraction partitions -- so an unmapped dataset of that kind surfaces as
#: unannotated rather than inheriting a wrong category.
LEVEL_KIND_CATEGORIES = {
    ("DNA", "CRISPR_screen"): ("fitness",),
    ("DNA", "imaging_screen"): ("phenotype",),
    ("DNA", "ATACseq"): ("chromatin",),
    ("DNA", "CUTandTAG"): ("chromatin",),
    ("DNA", "ChIPchip"): ("chromatin",),
    ("transcription", "RNAseq"): ("transcript",),
    ("transcription", "microarray"): ("transcript",),
    ("transcription", "scRNAseq"): ("transcript",),
    ("transcription", "LongRead"): ("transcript",),
    ("transcription", "RIPseq"): ("transcript",),
    ("transcription", "MeRIP"): ("transcript", "modification"),
    ("translation", "Ribo-seq"): ("translation",),
    ("translation", "RiboSeq"): ("translation",),
    ("translation", "proteomics"): ("protein_abundance",),
    ("post_translation", "LOPIT"): ("place",),
    ("post_translation", "BioID"): ("interaction",),
    ("post_translation", "proximity_labelling"): ("interaction",),
    ("post_translation", "IPMS"): ("interaction",),
    ("post_translation", "IP-MS"): ("interaction",),
    ("post_translation", "XLMS"): ("interaction",),
    ("post_translation", "phosphoproteomics"): ("modification",),
    ("post_translation", "structure"): ("structure",),
    ("post_translation", "computed structure"): ("structure",),
    ("reference", "domains"): ("domain",),
    ("reference", "identity"): ("identity",),
    ("reference", "orthology"): ("conservation",),
    ("reference", "sequence"): ("sequence",),
    ("reference", "structure"): ("structure",),
    ("reference", "variation"): ("conservation",),
    ("reference", "crosslink_MS"): ("interaction",),
    ("reference", "annotation"): ("sequence",),
}

# --------------------------------------------------------------------------- rule 3: the dataset
#: Datasets whose `(level, kind)` is true but not decisive, each with the reason read off its own
#: registry entry. This is where the three screens that bit us are nailed down, and where the levels
#: that bundle several kinds of knowledge (`post_translation/proteomics`) are resolved one deposit at
#: a time. Every entry is `key -> (categories, why)`.
DATASET_CATEGORIES = {
    # Pooled CRISPR screens whose family name says nothing about fitness. All three of the ones the
    # harness's importance table surfaced are here, and they are `fitness` by the level/kind rule
    # already -- recorded explicitly so the check is a statement and not a coincidence.
    "gra17_synthlethal": (("fitness",), "pooled CRISPR screen: RH vs RH-dgra17 depletion by passage"),
    "invivo_platform": (("fitness",), "pooled in-vivo CRISPR screen: mean log fold-change of guides"),
    "gra12": (("fitness",), "pooled CRISPR screens: median L2FC in vitro and in vivo, DISCO score"),
    "gse253884_5": (("fitness",),
                    "pooled in-vivo CRISPR screen of two targeted libraries. The column names say "
                    "'hyperlopit' only because the LIBRARY was chosen from hyperLOPIT-unassigned "
                    "proteins; the measurement is guide depletion, so this is fitness and NOT place. "
                    "A name-based rule gets this one exactly wrong."),
    "second_background_fitness": (("fitness",), "guide depletion over eight passages in a second strain"),
    "oxidative_stress_screen": (("fitness",), "pooled CRISPR screening score under oxidative challenge"),
    "curated_enteric_fitness": (("fitness",),
                                "literature-curated oocyst yield of gene disruptions through the "
                                "feline stage: fitness, though its level is reference/literature"),
    "hosttx_effectors": (("host_response", "place"),
                         "the host-transcription effector screen. Its level is DNA/CRISPR_screen but "
                         "the readout is the HOST transcriptional response (Hotelling T2 and the "
                         "per-effector signature), not parasite depletion. `place` as well since "
                         "2026-11-05: an effector that changes host transcription has very nearly "
                         "been shown to reach the host, so for a localisation question it is close "
                         "to the conclusion, and the category rule must be what excludes it rather "
                         "than a hand-written family list"),
    "differentiation_screen": (("phenotype", "fitness"),
                               "pooled CRISPR screen read out as guide enrichment in "
                               "reporter-positive parasites over the bulk. The ratio largely cancels "
                               "depletion, so the knowledge is a differentiation phenotype; "
                               "`fitness` is kept because it is a pooled screen and the exclusion "
                               "errs towards safety"),
    "spacr_escrt_screen": (("phenotype", "host_response", "fitness"),
                           "pooled image-based CRISPR screen of host ESCRT recruitment to the "
                           "vacuole: a perturbation phenotype that is measured on a host protein"),
    # Localisation assays, whatever the instrument.
    "lopit_tgon": (("place",), "hyperLOPIT / TAGM marker niches -- the canonical place columns"),
    "organelle_surface_turboid": (("interaction", "place"),
                                  "TurboID at an organelle surface: proximity (interaction) that is "
                                  "read as which surface the protein sits on (place)"),
    "pvm_proximity": (("place", "interaction"),
                      "PVM proximity labelling: the study's statement that the protein is AT the "
                      "parasitophorous vacuole membrane, from a proximity assay"),
    "secretome_partition": (("place",), "how a secreted protein splits between soluble and vesicular fractions"),
    # post_translation/proteomics, one deposit at a time.
    "pride_acetylation": (("modification",), "acetylation sites per gene"),
    "pride_glycosylation": (("modification",), "O-fucosylated peptides per gene"),
    "pride_lactylation": (("modification",), "lactylation sites per gene"),
    "pride_nitrosylation": (("modification",), "nitrosylation sites per gene"),
    "pride_ubiquitination": (("modification",), "GlyGly (ubiquitin / SUMO) sites per gene"),
    "toxodb_arginine_methylation": (("modification",), "arginine methylation sites per gene"),
    "toxodb_palmitome": (("modification",), "palmitoylation enrichment ratios"),
    "myristoylome": (("modification",), "N-myristoylation confidence"),
    "cdpk1_substrates": (("modification",), "thiophosphorylated peptides: a phosphorylation readout"),
    "pride_proximity": (("interaction",), "proximity partners per gene"),
    "melting_temperature": (("structure",), "a protein's melting temperature: a property of the fold"),
    "thermal_shift_cetsa": (("structure",),
                            "how far the melting curve moves with calcium: fold stability and "
                            "ligand engagement, not abundance"),
    # Levels that are true but too coarse.
    "gse302108_utr5": (("sequence",),
                       "5' UTR length, upstream AUGs, uORFs and Kozak score. The deposit is RNAseq, "
                       "but these columns are features of the transcript's SEQUENCE and carry no "
                       "abundance; a learner using them is not reading expression"),
    "toxodb_ec_numbers": (("function",),
                          "EC numbers state what the enzyme does. Filed under `phenotype` until "
                          "2026-11-05 for want of anywhere better; `function` now exists because a "
                          "conclusion about what a protein does is not a perturbation observation"),
    "iedb_bcell": (("host_response",),
                   "antibody epitopes: what the host immune system recognises on this protein"),
    "toxodb_epitopes": (("host_response",), "IEDB epitopes mapped to the gene, as above"),
    "curated_drug_sensitivity": (("phenotype",), "curated sensitivity shifts of mutants to compounds"),
    "curated_resistance_alleles": (("phenotype", "sequence"),
                                   "resistance substitutions: a drug phenotype carried by named "
                                   "sequence variants"),
    "pubmed": (("attention",), "publication and paper counts: study effort"),
    "pmc_oa": (("attention",), "full-text counts: study effort"),
    "toxodb_identity": (("identity",), "gene identifiers and product names"),
    "orthomcl": (("identity", "conservation"),
                 "the orthogroup accession. It is the protected GROUP key of every split, so it is "
                 "never a feature (`identity`), and it is also homology knowledge (`conservation`)"),
    "foldseek_struct": (("structure",), "structural neighbours at TM >= 0.7 over AlphaFold models"),
}

# --------------------------------------------------------------------------- rule 1: single columns
#: Exact per-column overrides, for the handful of columns that differ from their own dataset. Each
#: carries its reason, because every one of these is a place where the dataset-level answer is wrong.
COLUMN_CATEGORIES = {
    "iron_depletion_rna_log2fc": (("transcript",),
                                  "the RNA arm of a proteomics deposit: transcript, not abundance"),
    "best_model_agreement": (("interaction", "structure"),
                             "how well the measured crosslinks agree with the predicted complex: "
                             "an interaction measurement used as a structural check"),
    "n_holes": (("interaction", "attention", "transcript", "fitness"),
                "a degree count in the graph of pairs that co-expression AND co-fitness support "
                "while the literature does not co-mention them (`build_graph.structural_holes`). It "
                "embodies interaction knowledge, is defined by the ABSENCE of literature "
                "(`attention`), and is computed from co-expression and co-fitness, so all four are "
                "declared. Composite columns are where a single category would be a lie"),
    "attention_depth": (("attention",), "focal / substantive / incidental: how hard the gene was studied"),
    "lit_tier": (("attention",), "literature tier of the gene: study effort, no provenance recorded"),
    "sequence": (("sequence",), "the amino-acid sequence itself"),
    "dtm_class": (("sequence",), "DeepTMHMM class of the protein, read off the sequence"),
    "tm_kd_mean": (("sequence",), "hydrophobicity of the predicted TM helices"),
    "tm_kd_max": (("sequence",), "hydrophobicity of the predicted TM helices"),
    "tm_len_mean": (("sequence",), "length of the predicted TM helices"),
    "has_pf_ortholog": (("conservation",), "whether a Plasmodium falciparum ortholog exists"),
    "has_cp_ortholog": (("conservation",), "whether a Cryptosporidium ortholog exists"),
    "lineage_specific": (("conservation",), "no ortholog outside the lineage"),
    "has_structure": (("structure",), "whether a model of the fold exists at all"),
    "ortholopit_label": (("place",),
                         "hyperLOPIT niche transferred from an ortholog. TAGM-derived, so `place`: "
                         "the trap the trait name hides, since no dataset entry claims it"),
    "ortholopit_donors": (("place",), "which orthologs donated the transferred hyperLOPIT niche"),
    "ortholopit_accuracy": (("place",), "measured accuracy of the transferred hyperLOPIT niche"),
    "ortholopit_accepted": (("place",), "whether the transferred hyperLOPIT niche was accepted"),
    # The information map's own tier-1 conclusion columns (`values.parquet`).
    "host_facing": (("place",),
                    "the conclusion that the protein reaches beyond the parasite plasma membrane: a "
                    "statement about where it is"),
    "compartment_hyperlopit": (("place",), "hyperLOPIT marker niche"),
    "compartment_go_cc": (("place",), "experimental GO cellular component"),
    "compartment_go_cc_n": (("place",), "number of GO cellular-component terms recorded"),
    "partner_best": (("interaction",), "the strongest recorded partner"),
    "partner_best_prob": (("interaction",), "probability of the strongest recorded partner"),
    "partner_n_ge95": (("interaction",), "partners recorded at probability >= 0.95"),
    "partner_n_ge50": (("interaction",), "partners recorded at probability >= 0.50"),
    "function_class": (("function",),
                       "coarse functional class: a conclusion about what the protein does, which is "
                       "the thing a function question predicts. Its own category since 2026-11-05"),
    "function_class_n": (("function",), "as function_class"),
    "function_go_bp": (("function",), "experimental GO biological process; as function_class"),
    "function_go_bp_n": (("function",), "as function_go_bp"),
    "function_go_mf": (("function",), "experimental GO molecular function; as function_class"),
    "function_go_mf_n": (("function",), "as function_go_mf"),
    "ec_number": (("function",), "the EC number itself; as toxodb_ec_numbers"),
    "has_ec": (("function",), "whether an EC number is recorded; as toxodb_ec_numbers"),
    # Domain traits that are not purely domain content.
    "domclass_tm_rich": (("domain", "sequence"),
                         "three or more TM helices: a DeepTMHMM count, not a domain match"),
    "domclass_rhoptry_kinase_ropk": (("domain", "place"),
                                     "a domain match OR a protein kinase whose PRODUCT NAME matches "
                                     "rhoptry|ROP\\d+. The second arm reads a localisation word out "
                                     "of the annotation, so this column carries `place`"),
    # Sequence-trait columns that are structural rather than sequence.
    "tdi_source": (("structure",), "which structure the 3Di string came from"),
    "tdi_coverage": (("structure",), "fraction of residues with a real 3Di letter"),
    "tdi_n_fragments": (("structure",), "how many models were merged into the 3Di string"),
    "tdi_align_identity": (("structure",), "identity of the alignment that transferred the 3Di string"),
    "tdi_align_coverage": (("structure",), "coverage of the alignment that transferred the 3Di string"),
}

# --------------------------------------------------------------------------- rule 2: name overrides
#: Declared name patterns that BEAT the dataset default, each with its reason. There are only two
#: situations: Ribo-seq deposits ship an RNA arm beside the RPF and TE arms, and the information map
#: ships the imaging screen's phenotype columns. Both would otherwise inherit the deposit's category,
#: which is the error instruction 73 names explicitly ("translational-efficiency and Ribo-seq
#: columns are translation, not transcript" -- and the converse holds for their RNA arm).
OVERRIDE_PATTERNS = (
    (re.compile(r"^rna\d+_"), ("transcript",),
     "the RNA arm of a Ribo-seq or RNAseq deposit: RNA abundance"),
    (re.compile(r"^rpf\d+_"), ("translation",), "ribosome-protected fragments: translation"),
    (re.compile(r"^te\d+_"), ("translation",), "translational efficiency: translation, never transcript"),
    (re.compile(r"^screen_[a-z_]+_phenotype$|^screen_any_phenotype$|^screen_scorers_agree$"),
     ("phenotype",), "the split-Cas9 imaging screen's phenotype calls"),
)

# --------------------------------------------------------------------------- rule 5: last resort
#: Name patterns for columns with NO provenance at all: the trait tables built under `results/`, and
#: whatever a build script computed without registering a dataset. Everything resolved here is
#: reported as `judged`, so the review table shows exactly how much of the annotation rests on a name.
#: Each entry is `(pattern, categories, route, reason)`. `release_registry` means the derivation
#: comes from that release's own per-column registry -- `domain_lookup.tsv`, `column_registry.tsv`
#: or the documented build rules of `sequence_traits` -- which is evidence of the same kind as the
#: dataset registry, just kept beside the table instead of in `datasets.py`. `name_pattern` means
#: the name alone decided, and nothing shipped resolves that way.
PATTERNS = (
    (re.compile(r"^dom_(IPR|PF)\d+$"), ("domain",), "release_registry",
     "domain_lookup.tsv gives this column's InterPro or Pfam accession and description"),
    (re.compile(r"^domclass_"), ("domain",), "release_registry",
     "a named domain class built from the InterPro and Pfam lists quoted in the release README"),
    (re.compile(r"^(n_interpro_domains|n_pfam_domains|has_any_domain|domain_novel)$"),
     ("domain",), "release_registry", "counts and flags over the same domain table"),
    (re.compile(r"^tdi_"), ("structure",), "release_registry",
     "sequence_traits README: a Foldseek 3Di string, the structural alphabet and not the sequence"),
    (re.compile(r"^seq_|_len$|^n_tm$|^has_signal_peptide$|^flanks_identical$|_side$"),
     ("sequence",), "release_registry",
     "sequence_traits README: the sequence, its regions, their lengths and the DeepTMHMM topology"),
    (re.compile(r"(?i)^expr_(rel|abs|lfc|lfcse)_.*riboseq"), ("translation",), "release_registry",
     "column_registry.tsv gives this column's contrast as a RIBOseq library, so it is translation"),
    (re.compile(r"^expr_(rel|abs|lfc|lfcse)_"), ("transcript",), "release_registry",
     "column_registry.tsv gives the study, contrast and normalization of an RNA-abundance column"),
)

# --------------------------------------------------------------------------- the columns
# Generated by scripts/build_information_categories_review.py from the shipped tables; the
# manifest is checked against those files by tests/test_information_categories.py.

#: Every column of `parascape/data/nodes.parquet` (443: 411 numeric or boolean, 32 of text), the
#: table the map is built on.
NODE_COLUMNS = (
    "gene_id", "length", "sequence", "dtm_class",
    "n_tm", "has_signal_peptide", "is_tm", "tm_kd_mean",
    "tm_kd_max", "tm_len_mean", "lopit_confident", "interpro_id",
    "interpro_desc", "pfam_id", "n_interpro", "has_domain",
    "orthogroup", "product", "paralog_number", "has_pf_ortholog",
    "has_cp_ortholog", "lineage_specific", "mean_plddt", "has_structure",
    "rna108740_Tachyzoites_T2_FPKM", "rna108740_Tachyzoites_T4_FPKM", "rna108740_Tissue_cysts_A_FPKM", "rna108740_Tissue_cysts_B_FPKM",
    "rna108740_Day3_CS4_FPKM", "rna108740_Day3_CS5_FPKM", "rna108740_Day3_CS6_FPKM", "rna108740_Day5_CS2_FPKM",
    "rna108740_Day5_K4_14_FPKM", "rna108740_Day7_CS7_FPKM", "rna108740_Day7_CS8_FPKM", "rna108740_Day7_CS9_FPKM",
    "rna206344_Unsporulated_R1", "rna206344_Unsporulated_R2", "rna206344_Sporulating_R1", "rna206344_Sporulating_R2",
    "rna206344_Sporulated_R1", "rna206344_Sporulated_R2", "n_phosphosites", "has_phospho",
    "fit_invitro_hff", "fit_invivo_PE", "fit_invivo_lung", "fit_invivo_liver",
    "fit_invivo_spleen", "fit_naive_bmdm", "fit_ifng", "fit_invivo_young2019",
    "lopit_map", "lopit_mcmc", "lopit_prob_map", "lopit_prob_mcmc",
    "lopit_methods_agree", "compartment", "lopit_unified", "ortholopit_label",
    "ortholopit_donors", "ortholopit_accuracy", "ortholopit_accepted", "compartment_source",
    "compartment_best", "expr_tachy", "expr_cyst", "expr_max",
    "expr_sporulated", "crispr_gra17ko_phenotype", "crispr_gra17_synthlethal_delta", "crispr_gra17_candidate",
    "crispr_gra12s1_l2fc_invitro", "crispr_gra12s1_l2fc_invivo", "crispr_gra12s1_disco", "crispr_gra12s2_l2fc_invitro",
    "crispr_gra12s2_l2fc_invivo", "crispr_gra12s2_disco", "crispr_invivo_platform_lfc", "fit_hyperlopit_unassigned_invivo_lib1",
    "fit_hyperlopit_unassigned_invivo_lib2", "hosttx_T2", "hosttx_padj", "protein_ibaq_log2",
    "hosttx_signature_pc01", "hosttx_signature_pc02", "hosttx_signature_pc03", "hosttx_signature_pc04",
    "hosttx_signature_pc05", "hosttx_signature_pc06", "hosttx_signature_pc07", "hosttx_signature_pc08",
    "hosttx_signature_pc09", "hosttx_signature_pc10", "hosttx_signature_pc11", "hosttx_signature_pc12",
    "hosttx_signature_pc13", "hosttx_signature_pc14", "hosttx_signature_pc15", "hosttx_signature_pc16",
    "hosttx_signature_pc17", "hosttx_signature_pc18", "hosttx_signature_pc19", "hosttx_signature_pc20",
    "hosttx_signature_norm", "hosttx_signature_n_de", "invivo_TZ_1", "invivo_TZ_2",
    "invivo_WholeBrain_Acute_1", "invivo_WholeBrain_Acute_2", "invivo_WholeBrain_Acute_3", "invivo_WholeBrain_Chronic_1",
    "invivo_WholeBrain_Chronic_2", "invivo_WholeBrain_Chronic_3", "invivo_BZ_28DPI_1", "invivo_BZ_28DPI_3",
    "invivo_BZ_90DPI_1", "invivo_BZ_90DPI_2", "invivo_BZ_120DPI_2", "invivo_BZ_120DPI_3",
    "stress_5-1_unstress_red_24", "stress_5-2_unstress_red_24", "stress_5-24_stress_green_48", "stress_5-27_stress_green_48",
    "stress_6-1_unstress_red_24", "stress_6-4_unstress_red_24", "stress_6-16_stress_green_48", "stress_6-18_stress_green_48",
    "stress_7-1_unstress_red_24", "stress_7-18_stress_green_48", "rna22258_tachyzoite", "rna22258_bradyzoite",
    "brain168465_1d_base_mean", "brain168465_1d_lfc", "brain168465_2d_base_mean", "brain168465_2d_lfc",
    "brain168465_4d_base_mean", "brain168465_4d_lfc", "brain168465_7d_base_mean", "brain168465_7d_lfc",
    "brain168465_14d_base_mean", "brain168465_14d_lfc", "rpf99395_extracellular_r1", "rpf99395_extracellular_r2",
    "rpf99395_intracellular_r1", "rpf99395_intracellular_r2", "rna99395_extracellular_r1", "rna99395_extracellular_r2",
    "rna99395_intracellular_r1", "rna99395_intracellular_r2", "te99395_extracellular_r1", "te99395_extracellular_r2",
    "te99395_intracellular_r1", "te99395_intracellular_r2", "rpf129869_confluent_r1", "rpf129869_confluent_r2",
    "rpf129869_confluent_r3", "rna129869_confluent_r1", "rna129869_confluent_r2", "rna129869_confluent_r3",
    "rpf129869_subconfluent_r1", "rpf129869_subconfluent_r2", "rpf129869_subconfluent_r3", "rna129869_subconfluent_r1",
    "rna129869_subconfluent_r2", "rna129869_subconfluent_r3", "te129869_confluent_r1", "te129869_confluent_r2",
    "te129869_confluent_r3", "te129869_subconfluent_r1", "te129869_subconfluent_r2", "te129869_subconfluent_r3",
    "rpf245775_parent_tachy_r1", "rpf245775_parent_tachy_r2", "rpf245775_parent_tachy_r3", "rna245775_parent_tachy_r1",
    "rna245775_parent_tachy_r2", "rna245775_parent_tachy_r3", "rpf245775_parent_prebrady_r1", "rpf245775_parent_prebrady_r2",
    "rpf245775_parent_prebrady_r3", "rna245775_parent_prebrady_r1", "rna245775_parent_prebrady_r2", "rna245775_parent_prebrady_r3",
    "rpf245775_eif12ko_tachy_r1", "rpf245775_eif12ko_tachy_r2", "rpf245775_eif12ko_tachy_r3", "rna245775_eif12ko_tachy_r1",
    "rna245775_eif12ko_tachy_r2", "rna245775_eif12ko_tachy_r3", "rpf245775_eif12ko_prebrady_r1", "rpf245775_eif12ko_prebrady_r2",
    "rpf245775_eif12ko_prebrady_r3", "rna245775_eif12ko_prebrady_r1", "rna245775_eif12ko_prebrady_r2", "rna245775_eif12ko_prebrady_r3",
    "te245775_parent_tachy_r1", "te245775_parent_tachy_r2", "te245775_parent_tachy_r3", "te245775_parent_prebrady_r1",
    "te245775_parent_prebrady_r2", "te245775_parent_prebrady_r3", "te245775_eif12ko_tachy_r1", "te245775_eif12ko_tachy_r2",
    "te245775_eif12ko_tachy_r3", "te245775_eif12ko_prebrady_r1", "te245775_eif12ko_prebrady_r2", "te245775_eif12ko_prebrady_r3",
    "bfd2_rip_log2_ip_over_input", "cellcycle19092_async_r1", "cellcycle19092_async_r2", "cellcycle19092_blocked_r1",
    "cellcycle19092_blocked_r2", "cellcycle19092_1h_r1", "cellcycle19092_1h_r2", "cellcycle19092_2h_r1",
    "cellcycle19092_2h_r2", "cellcycle19092_3h_r1", "cellcycle19092_3h_r2", "cellcycle19092_4h_r1",
    "cellcycle19092_4h_r2", "cellcycle19092_5h_r1", "cellcycle19092_5h_r2", "cellcycle19092_6h_r1",
    "cellcycle19092_6h_r2", "cellcycle19092_7h_r1", "cellcycle19092_7h_r2", "cellcycle19092_8h_r1",
    "cellcycle19092_8h_r2", "cellcycle19092_9h_r1", "cellcycle19092_9h_r2", "cellcycle19092_10h_r1",
    "cellcycle19092_10h_r2", "cellcycle19092_11h_r1", "cellcycle19092_11h_r2", "cellcycle19092_12h_r1",
    "cellcycle19092_12h_r2", "rna51780_tachy_r1", "rna51780_tachy_r2", "rna51780_mero_r3",
    "rna51780_mero_r4", "rna51780_mero_r5", "rna51780_mero_r6", "cpsf4rna168155_ut_1",
    "cpsf4rna168155_ut_2", "cpsf4rna168155_iaa_7h_1", "cpsf4rna168155_iaa_7h_2", "cpsf4rna168155_iaa_24h_1",
    "cpsf4rna168155_iaa_24h_2", "cpsf4rna168155_iaa_48h_1", "cpsf4rna168155_iaa_48h_2", "restriction200962_7th_gt_2_kk31_kim_s1_l001",
    "restriction200962_7th_gt_20_kk27_kim_s14_l001", "restriction200962_7th_gt_22_kk29_kim_s15_l001", "restriction200962_7th_gt_4_kk32_kim_s2_l001", "restriction200962_kk17_p2wo_ph8_2_s5",
    "restriction200962_kk18_p2wo_ph8_2_s1", "restriction200962_kk19_p2w_ph8_2_s22", "restriction200962_kk20_p2w_ph8_2_s16", "restriction200962_kk21_p5wo_ph8_2_s13",
    "restriction200962_kk22_p5wo_ph8_2_s15", "restriction200962_kk23_p5w_ph8_2_s12", "restriction200962_kk24_p5w_ph8_2_s11", "restriction200962_kk5_p2wo_ph7_4_s6",
    "restriction200962_kk6_p2wo_ph7_4_s19", "restriction200962_kk7_p2w_ph7_4_s14", "restriction200962_kk8_p2w_ph7_4_s23", "morc_MORC_UT1",
    "morc_MORC_UT2", "morc_MORC_UT3", "morc_MORC_IAA1", "morc_MORC_IAA2",
    "morc_MORC_IAA3", "morc_MORC-KD-BFD1-KO_UT1", "morc_MORC-KD-BFD1-KO_UT2", "morc_MORC-KD-BFD1-KO_UT3",
    "morc_MORC-KD-BFD1-KO_IAA1", "morc_MORC-KD-BFD1-KO_IAA2", "morc_MORC-KD-BFD1-KO_IAA3", "morc_DD-BFD1-Ty_UT1",
    "morc_DD-BFD1-Ty_UT2", "morc_DD-BFD1-Ty_UT3", "morc_DD-BFD1-Ty_Shield1", "morc_DD-BFD1-Ty_Shield2",
    "morc_DD-BFD1-Ty_Shield3", "proteome_log2_normalized_and_imputated_abundances_UT_R1", "proteome_log2_normalized_and_imputated_abundances_UT_R2", "proteome_log2_normalized_and_imputated_abundances_UT_R3",
    "proteome_log2_normalized_and_imputated_abundances_T_24h_R1", "proteome_log2_normalized_and_imputated_abundances_T_24h_R2", "proteome_log2_normalized_and_imputated_abundances_T_24h_R3", "proteome_log2_normalized_and_imputated_abundances_T_32h_R1",
    "proteome_log2_normalized_and_imputated_abundances_T_32h_R2", "proteome_log2_normalized_and_imputated_abundances_T_32h_R3", "proteome_log2_normalized_and_imputated_abundances_T_48h_R1", "proteome_log2_normalized_and_imputated_abundances_T_48h_R2",
    "proteome_log2_normalized_and_imputated_abundances_T_48h_R3", "proteome_lfc_UT_Vs_T_24h_log2_fold_change", "proteome_lfc_UT_Vs_T_32h_log2_fold_change_1", "proteome_lfc_UT_Vs_T_48h_log2_fold_change_2",
    "phospho_up_sites", "phospho_up_ratio", "phospho_down_sites", "phospho_down_ratio",
    "phospho_sites_measured", "oocyst_itraq_115_113", "oocyst_itraq_116_113", "oocyst_itraq_115_114",
    "oocyst_itraq_116_114", "oocyst_itraq_115_113_1", "oocyst_itraq_116_113_1", "oocyst_itraq_115_114_1",
    "oocyst_itraq_116_114_1", "crispr_reporter_strain_p8_log2", "diff_reporter_log2_mNG_over_bulk", "cetsa_calcium_ed_score",
    "cyst_wall_max_spectral", "cyst_wall_n_baits", "oxidative_stress_screen_score", "n_acetylation_sites",
    "n_proximity_partners", "n_nitrosylation_sites", "n_ubiquitination_sites", "n_lactylation_sites",
    "n_o_fucosyl_peptides", "snp_total_all_strains", "snp_nonsynonymous", "snp_synonymous",
    "snp_noncoding", "snp_stop_codon", "codon_enc", "codon_gc3",
    "codon_cai_ribosomal", "atac_promoter_ut", "cuttag_hdac3_promoter_ut", "palmitome_odya_vs_hydroxylamine_log2",
    "palmitome_odya_vs_palmitate_log2", "iedb_epitope_count", "h4_acetylation_chip_score", "macrophage_expression_percentile",
    "escrt_recruitment_xgboost", "escrt_recruitment_maxvit", "n_m6a_peaks", "sexual_stage_8dpi_log2fc",
    "secretome_soluble_over_vesicle_log2", "antisense_expression_percentile", "melting_temperature_tm", "n_crosslink_partners",
    "pvm_proximity_positive", "cdpk1_thiophospho_peptides", "mrna_remaining_5h_actinomycin", "myristoylation_confidence",
    "n_arginine_methylation_sites", "novel_transcript_models", "ees_vs_tachyzoite_log2", "ec_number",
    "has_ec", "n_bcell_epitopes", "screen_actin_phenotype", "screen_apicoplast_phenotype",
    "screen_egress_phenotype", "screen_replication_phenotype", "screen_any_phenotype", "screen_scorers_agree",
    "resistance_allele_count", "resistance_compound_count", "resistance_substitutions", "resistance_compounds",
    "drug_compounds_tested", "drug_sensitivity_shifts", "drug_sensitivity_directions", "enteric_oocyst_yield",
    "enteric_sporulation", "enteric_measurements", "cellcycle_phase", "cellcycle_pseudotime",
    "stage_enriched_derived", "stage_margin_derived", "n_bioid_studies", "n_ipms_studies",
    "n_publications", "n_fulltext", "lit_tier", "n_papers_focal",
    "n_papers_substantive", "n_papers_incidental", "attention_depth", "n_xlink_partners",
    "n_struct_similar", "n_ipms_partners", "best_model_agreement", "n_host_targets",
    "n_holes", "af3_sequence_coverage", "af3_mean_plddt", "af3_plddt_q25",
    "af3_confident_sequence_fraction", "af3_very_confident_sequence_fraction", "af3_low_confidence_modelled_fraction", "af3_confident_rg_angstrom",
    "af3_confident_contacts_per_residue", "af3_ptm", "fit_invivo_heart", "fit_invivo_brain",
    "fit_lipid_rich_p8", "fit_lipid_limited_p8", "fit_lipid_rich_p4p5", "fit_lipid_limited_p4p5",
    "fit_serum_differential_p8", "fit_serum_differential_p4p5", "te302107_tachy_r1", "te302107_tachy_r2",
    "mrna_log2_remaining_4h_actinomycin", "fit_complete_medium_2025", "fit_no_glutamine", "fit_no_glucose",
    "fit_glucose_dependence", "fit_glucose_dependence_fdr", "utr5_length", "utr5_n_uaugs",
    "utr5_n_uorfs", "utr5_n_oorfs", "utr5_n_inframe_ext", "utr5_kozak_score",
    "bzsub_A_expr", "bzsub_B_expr", "bzsub_C_expr", "bzsub_D_expr",
    "bzsub_E_expr", "iron_depletion_protein_log2fc", "iron_depletion_protein_padj", "iron_depletion_rna_log2fc",
    "surface_apicoplast_log2fc", "surface_apicoplast_stringent", "surface_mitochondrion_log2fc", "surface_mitochondrion_stringent",
    "surface_er_log2fc", "surface_er_stringent", "fit_density_low", "fit_density_high",
    "fit_density_dependence", "fit_density_dependence_log10padj", "fit_density_dim",
)

#: `results/sequence_traits_2026_11_04/sequence_traits.parquet` (22).
SEQUENCE_TRAIT_COLUMNS = (
    "seq_full", "seq_mature", "seq_nterm_flank", "seq_cterm_flank",
    "tdi_full", "tdi_mature", "tdi_nterm_flank", "tdi_cterm_flank",
    "n_tm", "has_signal_peptide", "nterm_side", "cterm_side",
    "seq_len", "mature_len", "nterm_len", "cterm_len",
    "tdi_source", "tdi_coverage", "tdi_n_fragments", "tdi_align_identity",
    "tdi_align_coverage", "flanks_identical",
)

#: The 299 InterPro and Pfam families of `results/domain_traits_2026_11_04/domain_traits.parquet`
#: that occur in at least 10 genes, each shipped as a `dom_<accession>` flag.
DOMAIN_FAMILY_ACCESSIONS = (
    "IPR000008", "IPR000048", "IPR000152", "IPR000177", "IPR000182", "IPR000195", "IPR000209",
    "IPR000210", "IPR000222", "IPR000225", "IPR000253", "IPR000330", "IPR000387", "IPR000403",
    "IPR000408", "IPR000504", "IPR000571", "IPR000608", "IPR000626", "IPR000629", "IPR000717",
    "IPR000719", "IPR000727", "IPR000742", "IPR000795", "IPR000873", "IPR000884", "IPR001005",
    "IPR001163", "IPR001214", "IPR001251", "IPR001353", "IPR001394", "IPR001412", "IPR001471",
    "IPR001487", "IPR001594", "IPR001609", "IPR001611", "IPR001623", "IPR001650", "IPR001680",
    "IPR001752", "IPR001757", "IPR001806", "IPR001841", "IPR001849", "IPR001881", "IPR001932",
    "IPR001965", "IPR002048", "IPR002073", "IPR002110", "IPR002130", "IPR002347", "IPR002423",
    "IPR002553", "IPR003029", "IPR003323", "IPR003409", "IPR003439", "IPR003593", "IPR003607",
    "IPR003609", "IPR003959", "IPR003960", "IPR004000", "IPR004087", "IPR004088", "IPR004161",
    "IPR004273", "IPR004843", "IPR005225", "IPR006073", "IPR006186", "IPR006195", "IPR006387",
    "IPR007052", "IPR007125", "IPR007226", "IPR007863", "IPR008250", "IPR008271", "IPR008927",
    "IPR008978", "IPR008979", "IPR008984", "IPR008991", "IPR009000", "IPR009030", "IPR009057",
    "IPR009060", "IPR009072", "IPR009080", "IPR009091", "IPR010920", "IPR011009", "IPR011011",
    "IPR011012", "IPR011047", "IPR011053", "IPR011249", "IPR011333", "IPR011527", "IPR011545",
    "IPR011701", "IPR011990", "IPR011992", "IPR012337", "IPR012340", "IPR013087", "IPR013320",
    "IPR013584", "IPR013602", "IPR013766", "IPR014001", "IPR014014", "IPR014756", "IPR015424",
    "IPR015500", "IPR015655", "IPR015915", "IPR015940", "IPR015947", "IPR016024", "IPR016035",
    "IPR016135", "IPR016181", "IPR017441", "IPR017853", "IPR017871", "IPR017907", "IPR018108",
    "IPR018200", "IPR018247", "IPR018253", "IPR018303", "IPR019734", "IPR019775", "IPR019787",
    "IPR019821", "IPR020103", "IPR020472", "IPR020568", "IPR020845", "IPR020846", "IPR022086",
    "IPR022742", "IPR023088", "IPR023174", "IPR023298", "IPR023299", "IPR023313", "IPR023395",
    "IPR023753", "IPR023828", "IPR024317", "IPR024743", "IPR026983", "IPR027409", "IPR027410",
    "IPR027413", "IPR027417", "IPR027640", "IPR027916", "IPR028078", "IPR028352", "IPR028889",
    "IPR029000", "IPR029021", "IPR029033", "IPR029044", "IPR029045", "IPR029052", "IPR029055",
    "IPR029058", "IPR029063", "IPR029071", "IPR034204", "IPR035699", "IPR035706", "IPR035892",
    "IPR035969", "IPR035979", "IPR036047", "IPR036188", "IPR036236", "IPR036249", "IPR036259",
    "IPR036291", "IPR036322", "IPR036390", "IPR036412", "IPR036427", "IPR036457", "IPR036465",
    "IPR036612", "IPR036640", "IPR036691", "IPR036755", "IPR036770", "IPR036852", "IPR036855",
    "IPR036865", "IPR036866", "IPR036869", "IPR036890", "IPR036915", "IPR038765", "IPR041228",
    "IPR041466", "IPR041569", "IPR041658", "IPR043129", "IPR043502", "IPR043519", "IPR044492",
    "IPR045864", "IPR045877", "IPR046341", "IPR047575", "PF00004", "PF00005", "PF00009",
    "PF00022", "PF00024", "PF00063", "PF00069", "PF00071", "PF00076", "PF00082",
    "PF00085", "PF00106", "PF00118", "PF00122", "PF00125", "PF00149", "PF00153",
    "PF00160", "PF00176", "PF00179", "PF00225", "PF00226", "PF00227", "PF00233",
    "PF00240", "PF00270", "PF00271", "PF00400", "PF00415", "PF00439", "PF00443",
    "PF00481", "PF00498", "PF00501", "PF00566", "PF00583", "PF00642", "PF00650",
    "PF00664", "PF00702", "PF00847", "PF00856", "PF01399", "PF01423", "PF01529",
    "PF01602", "PF01926", "PF02338", "PF02493", "PF03028", "PF03144", "PF04092",
    "PF04969", "PF05193", "PF07690", "PF07992", "PF08373", "PF08393", "PF09717",
    "PF12146", "PF12314", "PF12774", "PF12775", "PF12777", "PF12780", "PF12781",
    "PF12796", "PF13499", "PF13540", "PF13639", "PF14295", "PF14531", "PF14580",
    "PF14733", "PF17852", "PF17862", "PF18198", "PF18199",
)

#: The named classes, counts and flags of the same table (14).
DOMAIN_SUMMARY_COLUMNS = (
    "domclass_protein_kinase", "domclass_ap2_apiap2", "domclass_srs_sag",
    "domclass_ef_hand", "domclass_wd40", "domclass_ankyrin",
    "domclass_rna_binding", "domclass_proteasome_ubiquitin", "domclass_rhoptry_kinase_ropk",
    "domclass_tm_rich", "n_interpro_domains", "n_pfam_domains",
    "has_any_domain", "domain_novel",
)

DOMAIN_TRAIT_COLUMNS = tuple(f"dom_{a}" for a in DOMAIN_FAMILY_ACCESSIONS) + DOMAIN_SUMMARY_COLUMNS

#: `results/expression_traits_2026_11_06/expression_traits.parquet` (571). Study, contrast and
#: normalization of each are in that release's `column_registry.tsv`. A Ribo-seq library column is
#: `translation` rather than `transcript`; there is none in this release (GSE245775, the one Ribo-seq
#: study, has no values in the snapshot and is not published), but the rule stays.
EXPRESSION_COLUMNS = (
    "expr_lfc_GSE103113_condition_gra18_vs_pru_treatment_none", "expr_lfcse_GSE103113_condition_gra18_vs_pru_treatment_none", "expr_lfc_GSE103113_condition_gra18_vs_pru_treatment_lps",
    "expr_lfcse_GSE103113_condition_gra18_vs_pru_treatment_lps", "expr_lfc_GSE103113_treatment_lps_vs_none_condition_pru", "expr_lfcse_GSE103113_treatment_lps_vs_none_condition_pru",
    "expr_lfc_GSE103113_treatment_lps_vs_none_condition_gra18", "expr_lfcse_GSE103113_treatment_lps_vs_none_condition_gra18", "expr_abs_GSE103113",
    "expr_rel_GSE116539_crt_ko", "expr_rel_GSE116539_wildtype", "expr_abs_GSE116539",
    "expr_lfc_GSE125120_condition_infected_vs_bystander_genotype_wild_type", "expr_lfcse_GSE125120_condition_infected_vs_bystander_genotype_wild_type", "expr_lfc_GSE125120_condition_infected_vs_bystander_genotype_myr1_ko",
    "expr_lfcse_GSE125120_condition_infected_vs_bystander_genotype_myr1_ko", "expr_lfc_GSE125120_genotype_myr1_ko_vs_wild_type_condition_bystander", "expr_lfcse_GSE125120_genotype_myr1_ko_vs_wild_type_condition_bystander",
    "expr_lfc_GSE125120_genotype_myr1_ko_vs_wild_type_condition_infected", "expr_lfcse_GSE125120_genotype_myr1_ko_vs_wild_type_condition_infected", "expr_abs_GSE125120",
    "expr_lfc_GSE132238_condition_stress_vs_unstress_genotype_wt_time_24_hr", "expr_lfcse_GSE132238_condition_stress_vs_unstress_genotype_wt_time_24_hr", "expr_lfc_GSE132238_condition_stress_vs_unstress_time_48_hr_genotype_wt",
    "expr_lfcse_GSE132238_condition_stress_vs_unstress_time_48_hr_genotype_wt", "expr_lfc_GSE132238_condition_stress_vs_unstress_time_72_hr_genotype_wt", "expr_lfcse_GSE132238_condition_stress_vs_unstress_time_72_hr_genotype_wt",
    "expr_lfc_GSE132238_condition_stress_vs_unstress_genotype_delta_bfd1_time_24_hr", "expr_lfcse_GSE132238_condition_stress_vs_unstress_genotype_delta_bfd1_time_24_hr", "expr_lfc_GSE132238_condition_stress_vs_unstress_genotype_delta_bfd1_time_48_hr",
    "expr_lfcse_GSE132238_condition_stress_vs_unstress_genotype_delta_bfd1_time_48_hr", "expr_lfc_GSE132238_condition_stress_vs_unstress_genotype_delta_bfd1_time_72_hr", "expr_lfcse_GSE132238_condition_stress_vs_unstress_genotype_delta_bfd1_time_72_hr",
    "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_condition_unstress_time_24_hr", "expr_lfcse_GSE132238_genotype_delta_bfd1_vs_wt_condition_unstress_time_24_hr", "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_24_hr",
    "expr_lfcse_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_24_hr", "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_time_48_hr_condition_unstress", "expr_lfcse_GSE132238_genotype_delta_bfd1_vs_wt_time_48_hr_condition_unstress",
    "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_time_72_hr_condition_unstress", "expr_lfcse_GSE132238_genotype_delta_bfd1_vs_wt_time_72_hr_condition_unstress", "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_48_hr",
    "expr_lfcse_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_48_hr", "expr_lfc_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_72_hr", "expr_lfcse_GSE132238_genotype_delta_bfd1_vs_wt_condition_stress_time_72_hr",
    "expr_lfc_GSE132238_time_48_hr_vs_24_hr_condition_unstress_genotype_wt", "expr_lfcse_GSE132238_time_48_hr_vs_24_hr_condition_unstress_genotype_wt", "expr_lfc_GSE132238_time_72_hr_vs_24_hr_condition_unstress_genotype_wt",
    "expr_lfcse_GSE132238_time_72_hr_vs_24_hr_condition_unstress_genotype_wt", "expr_lfc_GSE132238_time_72_hr_vs_48_hr_condition_unstress_genotype_wt", "expr_lfcse_GSE132238_time_72_hr_vs_48_hr_condition_unstress_genotype_wt",
    "expr_lfc_GSE132238_time_48_hr_vs_24_hr_condition_stress_genotype_wt", "expr_lfcse_GSE132238_time_48_hr_vs_24_hr_condition_stress_genotype_wt", "expr_lfc_GSE132238_time_48_hr_vs_24_hr_genotype_delta_bfd1_condition_unstress",
    "expr_lfcse_GSE132238_time_48_hr_vs_24_hr_genotype_delta_bfd1_condition_unstress", "expr_lfc_GSE132238_time_72_hr_vs_24_hr_condition_stress_genotype_wt", "expr_lfcse_GSE132238_time_72_hr_vs_24_hr_condition_stress_genotype_wt",
    "expr_lfc_GSE132238_time_72_hr_vs_48_hr_condition_stress_genotype_wt", "expr_lfcse_GSE132238_time_72_hr_vs_48_hr_condition_stress_genotype_wt", "expr_lfc_GSE132238_time_72_hr_vs_24_hr_genotype_delta_bfd1_condition_unstress",
    "expr_lfcse_GSE132238_time_72_hr_vs_24_hr_genotype_delta_bfd1_condition_unstress", "expr_lfc_GSE132238_time_72_hr_vs_48_hr_genotype_delta_bfd1_condition_unstress", "expr_lfcse_GSE132238_time_72_hr_vs_48_hr_genotype_delta_bfd1_condition_unstress",
    "expr_lfc_GSE132238_time_48_hr_vs_24_hr_condition_stress_genotype_delta_bfd1", "expr_lfcse_GSE132238_time_48_hr_vs_24_hr_condition_stress_genotype_delta_bfd1", "expr_lfc_GSE132238_time_72_hr_vs_24_hr_condition_stress_genotype_delta_bfd1",
    "expr_lfcse_GSE132238_time_72_hr_vs_24_hr_condition_stress_genotype_delta_bfd1", "expr_lfc_GSE132238_time_72_hr_vs_48_hr_condition_stress_genotype_delta_bfd1", "expr_lfcse_GSE132238_time_72_hr_vs_48_hr_condition_stress_genotype_delta_bfd1",
    "expr_abs_GSE132238", "expr_rel_GSE136123_Pru_FR_18h", "expr_rel_GSE136123_Pru_FR_7h",
    "expr_rel_GSE136123_Pru_MORC_KD_IAA_24h", "expr_rel_GSE136123_RH_MORC_KD_IAA_24h", "expr_rel_GSE136123_Pru_MORC_KD_IAA_48h",
    "expr_rel_GSE136123_RH_MORC_KD_IAA_7h", "expr_rel_GSE136123_Pru_MORC_KD_UT", "expr_rel_GSE136123_Pru_US",
    "expr_rel_GSE136123_RH_MORC_KD_UT", "expr_abs_GSE136123", "expr_lfc_GSE139537_condition_shield1_vs_etoh_genotype_wt_me49_delta_ku80",
    "expr_lfcse_GSE139537_condition_shield1_vs_etoh_genotype_wt_me49_delta_ku80", "expr_lfc_GSE139537_condition_shield1_vs_etoh_genotype_delta_bfd1_dd_bfd1_ty", "expr_lfcse_GSE139537_condition_shield1_vs_etoh_genotype_delta_bfd1_dd_bfd1_ty",
    "expr_lfc_GSE139537_genotype_delta_bfd1_dd_bfd1_ty_vs_wt_me49_delta_ku80_conditi", "expr_lfcse_GSE139537_genotype_delta_bfd1_dd_bfd1_ty_vs_wt_me49_delta_ku80_conditi", "expr_lfc_GSE139537_genotype_delta_bfd1_vs_wt_me49_delta_ku80_condition_shield1",
    "expr_lfcse_GSE139537_genotype_delta_bfd1_vs_wt_me49_delta_ku80_condition_shield1", "expr_lfc_GSE139537_genotype_delta_bfd1_dd_bfd1_ty_vs_wt_me49_delta_ku80_conditi_2", "expr_lfcse_GSE139537_genotype_delta_bfd1_dd_bfd1_ty_vs_wt_me49_delta_ku80_conditi_2",
    "expr_abs_GSE139537", "expr_lfc_GSE140197_time_12hr_vs_4hr_treatment_dmso_vehicle", "expr_lfcse_GSE140197_time_12hr_vs_4hr_treatment_dmso_vehicle",
    "expr_lfc_GSE140197_time_12hr_vs_8hr_treatment_dmso_vehicle", "expr_lfcse_GSE140197_time_12hr_vs_8hr_treatment_dmso_vehicle", "expr_lfc_GSE140197_time_8hr_vs_4hr_treatment_dmso_vehicle",
    "expr_lfcse_GSE140197_time_8hr_vs_4hr_treatment_dmso_vehicle", "expr_lfc_GSE140197_time_12hr_vs_4hr_treatment_24_2_nm_trte", "expr_lfcse_GSE140197_time_12hr_vs_4hr_treatment_24_2_nm_trte",
    "expr_lfc_GSE140197_time_12hr_vs_8hr_treatment_24_2_nm_trte", "expr_lfcse_GSE140197_time_12hr_vs_8hr_treatment_24_2_nm_trte", "expr_lfc_GSE140197_time_12hr_vs_4hr_treatment_6_nm_trte",
    "expr_lfcse_GSE140197_time_12hr_vs_4hr_treatment_6_nm_trte", "expr_lfc_GSE140197_time_12hr_vs_8hr_treatment_6_nm_trte", "expr_lfcse_GSE140197_time_12hr_vs_8hr_treatment_6_nm_trte",
    "expr_lfc_GSE140197_time_8hr_vs_4hr_treatment_24_2_nm_trte", "expr_lfcse_GSE140197_time_8hr_vs_4hr_treatment_24_2_nm_trte", "expr_lfc_GSE140197_time_8hr_vs_4hr_treatment_6_nm_trte",
    "expr_lfcse_GSE140197_time_8hr_vs_4hr_treatment_6_nm_trte", "expr_lfc_GSE140197_treatment_24_2_nm_trte_vs_dmso_vehicle_time_4hr", "expr_lfcse_GSE140197_treatment_24_2_nm_trte_vs_dmso_vehicle_time_4hr",
    "expr_lfc_GSE140197_treatment_6_nm_trte_vs_dmso_vehicle_time_4hr", "expr_lfcse_GSE140197_treatment_6_nm_trte_vs_dmso_vehicle_time_4hr", "expr_lfc_GSE140197_treatment_24_2_nm_trte_vs_dmso_vehicle_time_12hr",
    "expr_lfcse_GSE140197_treatment_24_2_nm_trte_vs_dmso_vehicle_time_12hr", "expr_lfc_GSE140197_treatment_24_2_nm_trte_vs_dmso_vehicle_time_8hr", "expr_lfcse_GSE140197_treatment_24_2_nm_trte_vs_dmso_vehicle_time_8hr",
    "expr_lfc_GSE140197_treatment_6_nm_trte_vs_dmso_vehicle_time_12hr", "expr_lfcse_GSE140197_treatment_6_nm_trte_vs_dmso_vehicle_time_12hr", "expr_lfc_GSE140197_treatment_6_nm_trte_vs_dmso_vehicle_time_8hr",
    "expr_lfcse_GSE140197_treatment_6_nm_trte_vs_dmso_vehicle_time_8hr", "expr_abs_GSE140197", "expr_rel_GSE150334_iKD_AP2IX_5_Auxin",
    "expr_rel_GSE150334_iKD_AP2IX_5_Auxin_2", "expr_abs_GSE150334", "expr_lfc_GSE156685_genotype_unstated_vs_wild_type_reference_background",
    "expr_lfcse_GSE156685_genotype_unstated_vs_wild_type_reference_background", "expr_abs_GSE156685", "expr_lfc_GSE158231_genotype_tgif2k_b_knock_out_vs_wild_type_treatment_untreated",
    "expr_lfcse_GSE158231_genotype_tgif2k_b_knock_out_vs_wild_type_treatment_untreated", "expr_lfc_GSE158231_genotype_tgif2k_b_knock_out_vs_wild_type_treatment_5um_sodiu", "expr_lfcse_GSE158231_genotype_tgif2k_b_knock_out_vs_wild_type_treatment_5um_sodiu",
    "expr_lfc_GSE158231_treatment_5um_sodium_arsenite_vs_untreated_genotype_wild_typ", "expr_lfcse_GSE158231_treatment_5um_sodium_arsenite_vs_untreated_genotype_wild_typ", "expr_lfc_GSE158231_treatment_5um_sodium_arsenite_vs_untreated_genotype_tgif2k_b",
    "expr_lfcse_GSE158231_treatment_5um_sodium_arsenite_vs_untreated_genotype_tgif2k_b", "expr_abs_GSE158231", "expr_lfc_GSE159515_stage_bradyzoites_vs_tachyzoites_strain_me49",
    "expr_lfcse_GSE159515_stage_bradyzoites_vs_tachyzoites_strain_me49", "expr_lfc_GSE159515_strain_rh_vs_me49_stage_tachyzoites", "expr_lfcse_GSE159515_strain_rh_vs_me49_stage_tachyzoites",
    "expr_abs_GSE159515", "expr_lfc_GSE175919_treatment_6_hours_of_treatment_with_mc1742_0_5_m_vs_6_hours", "expr_lfcse_GSE175919_treatment_6_hours_of_treatment_with_mc1742_0_5_m_vs_6_hours",
    "expr_abs_GSE175919", "expr_rel_GSE178355_METTL3_16h", "expr_rel_GSE178355_METTL3_4h",
    "expr_rel_GSE178355_METTL3_0h", "expr_rel_GSE178355_WTAP_16h", "expr_rel_GSE178355_WTAP_4h",
    "expr_rel_GSE178355_WTAP_0h", "expr_rel_GSE178355_YTH1_16h", "expr_rel_GSE178355_YTH1_4h",
    "expr_rel_GSE178355_YTH1_0h", "expr_abs_GSE178355", "expr_rel_GSE189677_RH",
    "expr_rel_GSE189677_RHpit", "expr_rel_GSE189677_iKD_PT2_ATc", "expr_rel_GSE189677_iKD_PT2_ATc_2",
    "expr_abs_GSE189677", "expr_lfc_GSE193614_genotype_unstated_vs_wild_type_reference_background", "expr_lfcse_GSE193614_genotype_unstated_vs_wild_type_reference_background",
    "expr_abs_GSE193614", "expr_lfc_GSE193779_genotype_ap2x_4_knock_out_vs_wild_type_reference_background", "expr_lfcse_GSE193779_genotype_ap2x_4_knock_out_vs_wild_type_reference_background",
    "expr_abs_GSE193779", "expr_lfc_GSE198001_CaM_indirect_knockout_with_Atc_rep1_iCaM_ATC_1", "expr_lfcse_GSE198001_CaM_indirect_knockout_with_Atc_rep1_iCaM_ATC_1",
    "expr_lfc_GSE198001_CaM_indirect_knockout_with_Atc_rep2_iCaM_ATC_2", "expr_lfcse_GSE198001_CaM_indirect_knockout_with_Atc_rep2_iCaM_ATC_2", "expr_lfc_GSE198001_CaM_indirect_knockout_with_Atc_rep3_iCaM_ATC_3",
    "expr_lfcse_GSE198001_CaM_indirect_knockout_with_Atc_rep3_iCaM_ATC_3", "expr_abs_GSE198001", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcyc5_bz",
    "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcyc5_bz", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcyc5_bz_2", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcyc5_bz_2",
    "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_bz", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_bz", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_bz_2",
    "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_bz_2", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_tz", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_tz",
    "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_tz_2", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_24h_with_auxin_vs_me49_tgcycp2_tz_2", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcyc5_bz_aux_re",
    "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcyc5_bz_aux_re", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcyc5_bz_aux_re_2", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcyc5_bz_aux_re_2",
    "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_bz_aux_r", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_bz_aux_r", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_bz_aux_r_2",
    "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_bz_aux_r_2", "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_tz_aux_r", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_tz_aux_r",
    "expr_lfc_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_tz_aux_r_2", "expr_lfcse_GSE200962_treatment_48h_ph_7_4_5_co2_no_auxin_vs_me49_tgcycp2_tz_aux_r_2", "expr_lfc_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy",
    "expr_lfcse_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy", "expr_lfc_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2", "expr_lfcse_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2",
    "expr_lfc_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3", "expr_lfcse_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3", "expr_lfc_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3_4",
    "expr_lfcse_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3_4", "expr_lfc_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3_4_5", "expr_lfcse_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3_4_5",
    "expr_lfc_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3_4_5_6", "expr_lfcse_GSE200962_treatment_72h_ph_8_2_ambient_air_24h_with_auxin_vs_me49_tgcy_2_3_4_5_6", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep1",
    "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep1", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep2", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep2",
    "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_rep", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_rep", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_rep_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_rep_2", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_rep", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_rep",
    "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_rep_2", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_rep_2", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep1",
    "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep1", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep2", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep2",
    "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_rep", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_rep", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_rep_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_rep_2", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_rep", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_rep",
    "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_rep_2", "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_rep_2", "expr_lfc_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep1_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcyc5_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep1_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_re",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_re", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_re", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_re", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_re_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcycp2_tz_aux_re_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_re",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_re", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_re", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_re", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_re_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_tz_aux_re_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3_4", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3_4", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3", "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3_4", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3_4",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2_3", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_bz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2_3", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_tz_aux_re", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_tz_aux_re",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_tz_aux_re_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_tz_aux_re_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_tz_aux_re", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_tz_aux_re",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_tz_aux_re_2", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_tz_aux_re_2", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3_4", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcyc5_bz_aux_rep_2_3_4",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2_3", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2_3", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2_3_4",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep1_vs_me49_tgcycp2_bz_aux_re_2_3_4", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3_4", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcyc5_bz_aux_rep_2_3_4", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2_3",
    "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2_3", "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2_3_4", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_bz_aux_re_2_3_4",
    "expr_lfc_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_tz_aux_re_2_3", "expr_lfcse_GSE200962_treatment_me49_tgcycp2_tz_aux_rep2_vs_me49_tgcycp2_tz_aux_re_2_3", "expr_abs_GSE200962",
    "expr_rel_GSE208647_ME49EW_cyst_30_days", "expr_rel_GSE208647_ME49_Tir1_strain_pH7_0_48hrs", "expr_rel_GSE208647_ME49_Tir1_strain_pH8_2_48hrs",
    "expr_abs_GSE208647", "expr_lfc_GSE218921_genotype_bfd1_ko_vs_wild_type_tissue_replicating_tachyzoites", "expr_lfcse_GSE218921_genotype_bfd1_ko_vs_wild_type_tissue_replicating_tachyzoites",
    "expr_lfc_GSE218921_genotype_rocy1_ko_vs_wild_type_tissue_replicating_tachyzoite", "expr_lfcse_GSE218921_genotype_rocy1_ko_vs_wild_type_tissue_replicating_tachyzoite", "expr_lfc_GSE218921_genotype_bfd1_ko_vs_wild_type_treatment_high_ph_low_c02_low",
    "expr_lfcse_GSE218921_genotype_bfd1_ko_vs_wild_type_treatment_high_ph_low_c02_low", "expr_lfc_GSE218921_genotype_rocy1_ko_vs_wild_type_treatment_high_ph_low_c02_low", "expr_lfcse_GSE218921_genotype_rocy1_ko_vs_wild_type_treatment_high_ph_low_c02_low",
    "expr_lfc_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me", "expr_lfcse_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me", "expr_lfc_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me_2",
    "expr_lfcse_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me_2", "expr_lfc_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me_2_3", "expr_lfcse_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me_2_3",
    "expr_lfc_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me_2_3_4", "expr_lfcse_GSE218921_treatment_high_ph_low_c02_low_serum_conditions_vs_control_me_2_3_4", "expr_abs_GSE218921",
    "expr_lfc_GSE223621_genotype_pru_ova_tdtomato_delta_bfd1_mneongreen_vs_pru_ova_t", "expr_lfcse_GSE223621_genotype_pru_ova_tdtomato_delta_bfd1_mneongreen_vs_pru_ova_t", "expr_lfc_GSE223621_genotype_pru_ova_tdtomato_delta_bfd2_hxgprt_vs_pru_ova_tdtom",
    "expr_lfcse_GSE223621_genotype_pru_ova_tdtomato_delta_bfd2_hxgprt_vs_pru_ova_tdtom", "expr_abs_GSE223621", "expr_lfc_GSE223819_genotype_me49_tir1_200385_mng_maid_vs_me49_tir1_treatment_ve",
    "expr_lfcse_GSE223819_genotype_me49_tir1_200385_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfc_GSE223819_genotype_me49_tir1_208020_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfcse_GSE223819_genotype_me49_tir1_208020_mng_maid_vs_me49_tir1_treatment_ve",
    "expr_lfc_GSE223819_genotype_me49_tir1_224630_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfcse_GSE223819_genotype_me49_tir1_224630_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfc_GSE223819_genotype_me49_tir1_253790_mng_maid_vs_me49_tir1_treatment_ve",
    "expr_lfcse_GSE223819_genotype_me49_tir1_253790_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfc_GSE223819_genotype_me49_tir1_306620_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfcse_GSE223819_genotype_me49_tir1_306620_mng_maid_vs_me49_tir1_treatment_ve",
    "expr_lfc_GSE223819_genotype_me49_tir1_311100_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfcse_GSE223819_genotype_me49_tir1_311100_mng_maid_vs_me49_tir1_treatment_ve", "expr_lfc_GSE223819_genotype_me49_tir1_200385_mng_maid_vs_me49_tir1_treatment_50",
    "expr_lfcse_GSE223819_genotype_me49_tir1_200385_mng_maid_vs_me49_tir1_treatment_50", "expr_lfc_GSE223819_genotype_me49_tir1_208020_mng_maid_vs_me49_tir1_treatment_50", "expr_lfcse_GSE223819_genotype_me49_tir1_208020_mng_maid_vs_me49_tir1_treatment_50",
    "expr_lfc_GSE223819_genotype_me49_tir1_224630_mng_maid_vs_me49_tir1_treatment_50", "expr_lfcse_GSE223819_genotype_me49_tir1_224630_mng_maid_vs_me49_tir1_treatment_50", "expr_lfc_GSE223819_genotype_me49_tir1_253790_mng_maid_vs_me49_tir1_treatment_50",
    "expr_lfcse_GSE223819_genotype_me49_tir1_253790_mng_maid_vs_me49_tir1_treatment_50", "expr_lfc_GSE223819_genotype_me49_tir1_306620_mng_maid_vs_me49_tir1_treatment_50", "expr_lfcse_GSE223819_genotype_me49_tir1_306620_mng_maid_vs_me49_tir1_treatment_50",
    "expr_lfc_GSE223819_genotype_me49_tir1_311100_mng_maid_vs_me49_tir1_treatment_50", "expr_lfcse_GSE223819_genotype_me49_tir1_311100_mng_maid_vs_me49_tir1_treatment_50", "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1",
    "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1", "expr_lfc_GSE223819_treatment_5um_iaa_vs_vehicle_pbs_genotype_me49_tir1", "expr_lfcse_GSE223819_treatment_5um_iaa_vs_vehicle_pbs_genotype_me49_tir1",
    "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_200385", "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_200385", "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_208020",
    "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_208020", "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_224630", "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_224630",
    "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_253790", "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_253790", "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_306620",
    "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_306620", "expr_lfc_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_311100", "expr_lfcse_GSE223819_treatment_50um_iaa_vs_vehicle_pbs_genotype_me49_tir1_311100",
    "expr_abs_GSE223819", "expr_lfc_GSE223869_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_bfd2_vs_me49d", "expr_lfcse_GSE223869_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_bfd2_vs_me49d",
    "expr_lfc_GSE223869_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_bfd2_vs_me49d_2", "expr_lfcse_GSE223869_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_bfd2_vs_me49d_2", "expr_lfc_GSE223869_treatment_alkaline_stress_media_ambient_co2_vs_standard_cult",
    "expr_lfcse_GSE223869_treatment_alkaline_stress_media_ambient_co2_vs_standard_cult", "expr_lfc_GSE223869_treatment_alkaline_stress_media_ambient_co2_vs_standard_cult_2", "expr_lfcse_GSE223869_treatment_alkaline_stress_media_ambient_co2_vs_standard_cult_2",
    "expr_abs_GSE223869", "expr_lfc_GSE223877_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_hxgprt_ptub1", "expr_lfcse_GSE223877_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_hxgprt_ptub1",
    "expr_lfc_GSE223877_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_hxgprt_ptub1_2", "expr_lfcse_GSE223877_genotype_me49delta_ku80delta_bfd1_bfd1_tydelta_hxgprt_ptub1_2", "expr_lfc_GSE223877_treatment_standard_medium_3um_shield_1_vs_standard_medium_ve",
    "expr_lfcse_GSE223877_treatment_standard_medium_3um_shield_1_vs_standard_medium_ve", "expr_lfc_GSE223877_treatment_standard_medium_3um_shield_1_vs_standard_medium_ve_2", "expr_lfcse_GSE223877_treatment_standard_medium_3um_shield_1_vs_standard_medium_ve_2",
    "expr_abs_GSE223877", "expr_lfc_GSE224934_iKD_AP2XII_1_IAA_12h_biological_repeat", "expr_lfcse_GSE224934_iKD_AP2XII_1_IAA_12h_biological_repeat",
    "expr_lfc_GSE224934_iKD_AP2XII_1_IAA_24h_biological_repeat", "expr_lfcse_GSE224934_iKD_AP2XII_1_IAA_24h_biological_repeat", "expr_abs_GSE224934",
    "expr_lfc_GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_1", "expr_lfcse_GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_1", "expr_lfc_GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_2",
    "expr_lfcse_GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_2", "expr_lfc_GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_3", "expr_lfcse_GSE228580_ME49_iKD_AP2XII_1_IAA_60h_biological_repeat_3",
    "expr_abs_GSE228580", "expr_rel_GSE240648_cKd_ATc1", "expr_rel_GSE240648_cKd_ATc2",
    "expr_rel_GSE240648_cKd_ATc3", "expr_rel_GSE240648_cKd_ATc4", "expr_rel_GSE240648_cKd",
    "expr_abs_GSE240648", "expr_rel_GSE268652_UT", "expr_rel_GSE268652_IAA",
    "expr_abs_GSE268652", "expr_rel_GSE268750_UT", "expr_rel_GSE268750_IAA",
    "expr_abs_GSE268750", "expr_rel_GSE268892_iSNF2L_IAA_0h_RNA_seq", "expr_rel_GSE268892_iSNF2L_IAA_12h_RNA_seq",
    "expr_rel_GSE268892_iSNF2L_IAA_24h_RNA_seq", "expr_abs_GSE268892", "expr_lfc_GSE271902_Illumina_RFTS_KD_IAA",
    "expr_lfcse_GSE271902_Illumina_RFTS_KD_IAA", "expr_lfc_GSE271902_Illumina_SNF2h_KD_IAA", "expr_lfcse_GSE271902_Illumina_SNF2h_KD_IAA",
    "expr_lfc_GSE271902_Nanopore_SNF2h_KD_IAA", "expr_lfcse_GSE271902_Nanopore_SNF2h_KD_IAA", "expr_abs_GSE271902",
    "expr_lfc_GSE285297_KO", "expr_lfcse_GSE285297_KO", "expr_abs_GSE285297",
    "expr_rel_GSE285353_RH_Tir1_IAA", "expr_rel_GSE285353_RH_iAP2XI_3_IAA", "expr_abs_GSE285353",
    "expr_rel_GSE286090_RNA_Seq_G5A_HA_tz1", "expr_rel_GSE286090_RNA_Seq_G5A_HA_tz2", "expr_rel_GSE286090_RNA_Seq_G5A_HA_tz3",
    "expr_rel_GSE286090_RNA_Seq_G5A_KO_tz1", "expr_rel_GSE286090_RNA_Seq_G5A_KO_tz2", "expr_rel_GSE286090_RNA_Seq_G5A_KO_tz3",
    "expr_rel_GSE286090_RNA_Seq_G5A_HA_bz1", "expr_rel_GSE286090_RNA_Seq_G5A_HA_bz2", "expr_rel_GSE286090_RNA_Seq_G5A_HA_bz3",
    "expr_rel_GSE286090_RNA_Seq_G5A_KO_bz1", "expr_rel_GSE286090_RNA_Seq_G5A_KO_bz2", "expr_rel_GSE286090_RNA_Seq_G5A_KO_bz3",
    "expr_abs_GSE286090", "expr_lfc_GSE286338_A", "expr_lfcse_GSE286338_A",
    "expr_lfc_GSE286338_B", "expr_lfcse_GSE286338_B", "expr_lfc_GSE286338_C",
    "expr_lfcse_GSE286338_C", "expr_lfc_GSE286338_D", "expr_lfcse_GSE286338_D",
    "expr_lfc_GSE286338_E", "expr_lfcse_GSE286338_E", "expr_lfc_GSE286338_F",
    "expr_lfcse_GSE286338_F", "expr_lfc_GSE286338_G", "expr_lfcse_GSE286338_G",
    "expr_lfc_GSE286338_H", "expr_lfcse_GSE286338_H", "expr_abs_GSE286338",
    "expr_rel_GSE303275_iAlg11ATc", "expr_rel_GSE303275_iAlg11", "expr_abs_GSE303275",
    "expr_rel_GSE304633_RH_iAP2XI_3_IAA_2h", "expr_rel_GSE304633_RH_iAP2XI_3_IAA_8h", "expr_rel_GSE304633_RH_iAP2XI_3_IAA_0h",
    "expr_abs_GSE304633", "expr_lfc_GSE313582_Illumina_GCN5b_KD_IAA24h", "expr_lfcse_GSE313582_Illumina_GCN5b_KD_IAA24h",
    "expr_lfc_GSE313582_Illumina_PHD1_KD_IAA24h", "expr_lfcse_GSE313582_Illumina_PHD1_KD_IAA24h", "expr_lfc_GSE313582_Illumina_GCN5b_KD_IAA48h",
    "expr_lfcse_GSE313582_Illumina_GCN5b_KD_IAA48h", "expr_lfc_GSE313582_Illumina_PHD1_KD_IAA48h", "expr_lfcse_GSE313582_Illumina_PHD1_KD_IAA48h",
    "expr_abs_GSE313582", "expr_rel_GSE314937_iALG11_ATc1", "expr_rel_GSE314937_iALG11_ATc2",
    "expr_rel_GSE314937_iALG11_ATc3", "expr_rel_GSE314937_iALG11", "expr_abs_GSE314937",
    "expr_rel_GSE334992_Toxoplasma_gondii_RH_strain_TGGT1_245740_knock", "expr_rel_GSE334992_Toxoplasma_gondii_RH_strain_TGGT1_245740_knock_2", "expr_rel_GSE334992_Toxoplasma_gondii_RH_strain_TGGT1_245740_knock_3",
    "expr_rel_GSE334992_Toxoplasma_gondii_RH_strain_wild_type_tachyzoi", "expr_abs_GSE334992", "expr_lfc_GSE343877_Toxoplasma_gondii_ME49_bradyzoites_eIF5A_K49Q",
    "expr_lfcse_GSE343877_Toxoplasma_gondii_ME49_bradyzoites_eIF5A_K49Q", "expr_lfc_GSE343877_Toxoplasma_gondii_ME49_bradyzoites_eIF5A_K49R", "expr_lfcse_GSE343877_Toxoplasma_gondii_ME49_bradyzoites_eIF5A_K49R",
    "expr_lfc_GSE343877_Toxoplasma_gondii_ME49_tachyzoites_eIF5A_K49Q", "expr_lfcse_GSE343877_Toxoplasma_gondii_ME49_tachyzoites_eIF5A_K49Q", "expr_lfc_GSE343877_Toxoplasma_gondii_ME49_tachyzoites_eIF5A_K49R",
    "expr_lfcse_GSE343877_Toxoplasma_gondii_ME49_tachyzoites_eIF5A_K49R", "expr_abs_GSE343877", "expr_rel_GSE98227_B4_P3_fastq",
    "expr_abs_GSE98227", "expr_lfc_PRJNA1008424_genotype_trx1_mutant_vs_wild_type_reference_background", "expr_lfcse_PRJNA1008424_genotype_trx1_mutant_vs_wild_type_reference_background",
    "expr_abs_PRJNA1008424", "expr_abs_PRJNA1064470", "expr_abs_PRJNA1121513",
    "expr_lfc_PRJNA1133933_Treatment", "expr_lfcse_PRJNA1133933_Treatment", "expr_abs_PRJNA1133933",
    "expr_rel_PRJNA1148922_AP2XII_8_C", "expr_rel_PRJNA1148922_AP2XII_8_IAA12", "expr_rel_PRJNA1148922_AP2XII_8_IAA24",
    "expr_abs_PRJNA1148922", "expr_abs_PRJNA1164455", "expr_abs_PRJNA1165311",
    "expr_abs_PRJNA1165326", "expr_abs_PRJNA1233816", "expr_rel_PRJNA1252425_CTL",
    "expr_rel_PRJNA1252425_PIC", "expr_abs_PRJNA1252425", "expr_rel_PRJNA1282873_D1",
    "expr_rel_PRJNA1282873_D2", "expr_rel_PRJNA1282873_D3", "expr_rel_PRJNA1282873_V1",
    "expr_rel_PRJNA1282873_V2", "expr_rel_PRJNA1282873_V3", "expr_abs_PRJNA1282873",
    "expr_rel_PRJNA1309712_ROP64_T1", "expr_rel_PRJNA1309712_ROP64_T2", "expr_rel_PRJNA1309712_ROP64_T3",
    "expr_rel_PRJNA1309712_ROP65_T1", "expr_rel_PRJNA1309712_ROP65_T2", "expr_rel_PRJNA1309712_ROP65_T3",
    "expr_rel_PRJNA1309712_RH_T1", "expr_rel_PRJNA1309712_RH_T2", "expr_rel_PRJNA1309712_RH_T3",
    "expr_abs_PRJNA1309712", "expr_rel_PRJNA1446183_gra47_1", "expr_rel_PRJNA1446183_gra47_2",
    "expr_rel_PRJNA1446183_gra47_3", "expr_rel_PRJNA1446183_RH_1", "expr_rel_PRJNA1446183_RH_2",
    "expr_rel_PRJNA1446183_RH_3", "expr_abs_PRJNA1446183", "expr_abs_PRJNA1522152",
    "expr_abs_PRJNA1531573", "expr_abs_PRJNA252680", "expr_lfc_PRJNA705300_Model_organism_or_animal_sample_from_Toxoplasm",
    "expr_lfcse_PRJNA705300_Model_organism_or_animal_sample_from_Toxoplasm", "expr_lfc_PRJNA705300_Model_organism_or_animal_sample_from_Toxoplasm_2", "expr_lfcse_PRJNA705300_Model_organism_or_animal_sample_from_Toxoplasm_2",
    "expr_abs_PRJNA705300", "expr_abs_PRJNA725560", "expr_rel_PRJNA725840_Model_organism_or_animal_sample_from_Toxoplasm",
    "expr_rel_PRJNA725840_Model_organism_or_animal_sample_from_Toxoplasm_2", "expr_abs_PRJNA725840", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU",
    "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_2", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_3", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_4",
    "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_5", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_6", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_7",
    "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_8", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_9", "expr_rel_PRJNA738301_Microbe_sample_from_Toxoplasma_gondii_PRU_10",
    "expr_abs_PRJNA738301",
)

#: `results/information_map_v0_2026_11_04/values.parquet` (20): the map's own trait columns.
INFORMATION_MAP_COLUMNS = (
    "host_facing", "compartment_hyperlopit", "compartment_go_cc",
    "compartment_go_cc_n", "partner_best", "partner_best_prob",
    "partner_n_ge95", "partner_n_ge50", "function_class",
    "function_class_n", "function_go_bp", "function_go_bp_n",
    "function_go_mf", "function_go_mf_n", "screen_actin_phenotype",
    "screen_apicoplast_phenotype", "screen_egress_phenotype", "screen_replication_phenotype",
    "screen_any_phenotype", "screen_scorers_agree",
)

#: The shipped tables, each as the exact column list of the file named. Enumerated rather than read
#: from disk so that `unannotated()` is answerable in an installed copy with no `results/` directory,
#: and so a column added to a table without a category shows up as a diff in this file.
TABLES = {
    "nodes.parquet": NODE_COLUMNS,
    "sequence_traits.parquet": SEQUENCE_TRAIT_COLUMNS,
    "domain_traits.parquet": DOMAIN_TRAIT_COLUMNS,
    "expression_traits.parquet": EXPRESSION_COLUMNS,
    "information_map/values.parquet": INFORMATION_MAP_COLUMNS,
}

#: Every column of every shipped table, each named once, in table order.
SHIPPED_COLUMNS = tuple(dict.fromkeys(column for columns in TABLES.values() for column in columns))

#: How a column's categories were arrived at, strongest evidence first. `provenance:dataset` is as
#: hard as `provenance:level_kind` -- it reads the registry entry for that deposit -- while
#: `name_pattern` is the only route that rests on a name alone.
DERIVED_FROM = ("column_override", "name_override", "provenance:dataset", "provenance:level_kind",
                "release_registry", "name_pattern", "none")
#: Which routes count as derived from recorded provenance rather than judged here, for the honesty
#: count. `release_registry` is included: the per-column registry shipped beside a trait table states
#: the study, contrast or accession of every column in it, which is the same kind of evidence as
#: `datasets.py`. Everything else -- an exact override, a name that beats its dataset, a bare name
#: pattern -- is a judgement and is listed as one in the review table and the README.
FROM_PROVENANCE = ("provenance:dataset", "provenance:level_kind", "release_registry")


def _provenance(column: str, organism: str):
    """The dataset behind a node column, preferring the organism's own table over a shared name."""
    found = datasets.provenance(column, organism)
    return found if found is not None else datasets.provenance(column)


def _resolve(column: str, organism: str) -> tuple:
    """(categories, derived_from, reason, provenance_key, level_kind) for one column, by the rule order.

    The order is the point: an exact override, then a declared name override, then the dataset, then
    its level and kind, then -- only for columns with no provenance -- a name pattern. Anything that
    falls through comes back empty, which is what `unannotated` reports and the gate test refuses.
    """
    found = _provenance(column, organism)
    level_kind = f"{found.level}/{found.kind}" if found is not None else ""
    default = LEVEL_KIND_CATEGORIES.get((found.level, found.kind), ()) if found is not None else ()
    recorded = DATASET_CATEGORIES[found.key][0] if found is not None and found.key in DATASET_CATEGORIES \
        else default
    if column in COLUMN_CATEGORIES:
        categories, why = COLUMN_CATEGORIES[column]
        # An override that says what the registry already says has overridden nothing, and must not
        # be counted as a judgement: the honesty count is only useful if it is exact.
        if categories == recorded:
            return (categories, "provenance:dataset" if found.key in DATASET_CATEGORIES
                    else "provenance:level_kind",
                    f"{level_kind}, and the column entry agrees: {why}", found.key, level_kind)
        return categories, "column_override", why, found.key if found else "", level_kind
    for pattern, categories, why in OVERRIDE_PATTERNS:
        if pattern.search(column):
            # A declared name rule that agrees with the dataset's own level and kind has not
            # overridden anything, so it is reported as what it is: provenance.
            if categories == default:
                return (categories, "provenance:level_kind",
                        f"level/kind {level_kind}, and the name agrees: {why}", found.key, level_kind)
            return categories, "name_override", why, found.key if found else "", level_kind
    if found is not None:
        if found.key in DATASET_CATEGORIES:
            categories, why = DATASET_CATEGORIES[found.key]
            return categories, "provenance:dataset", why, found.key, level_kind
        if (found.level, found.kind) in LEVEL_KIND_CATEGORIES:
            return (LEVEL_KIND_CATEGORIES[(found.level, found.kind)], "provenance:level_kind",
                    f"level/kind {level_kind}: {found.provides}", found.key, level_kind)
        return (), "none", f"dataset {found.key} has an unmapped level/kind {level_kind}", found.key, level_kind
    for pattern, categories, route, why in PATTERNS:
        if pattern.search(column):
            return categories, route, why, "", ""
    return (), "none", "no provenance and no rule matches: left unannotated on purpose", "", ""


def categories_of(column: str, organism: str = organisms.TOXOPLASMA) -> tuple:
    """The information categories one column carries, or `()` when it has none.

    Empty is a real answer and means the column is unannotated: callers must treat it as unusable
    rather than as uncategorised-therefore-safe, which is the failure mode this module exists to end.
    """
    return _resolve(str(column), organism)[0]


def derivation(column: str, organism: str = organisms.TOXOPLASMA) -> dict:
    """One reviewable row: the column, its categories, how they were derived and how confident that is.

    `confidence` is `high` when the categories come from recorded provenance -- the dataset registry's
    level and kind, its own entry, or the per-column registry shipped with a trait release --
    `judged` when this module decided against or without that record, and `none` when nothing
    matched at all.
    """
    categories, route, why, key, level_kind = _resolve(str(column), organism)
    confidence = "high" if route in FROM_PROVENANCE else ("judged" if categories else "none")
    return {"column": str(column), "categories": ",".join(categories), "derived_from": route,
            "reason": why, "provenance_key": key, "provenance_level_kind": level_kind,
            "confidence": confidence}


def derivations(columns=None, organism: str = organisms.TOXOPLASMA) -> list:
    """`derivation` for every shipped column (or the given ones), with the table each came from."""
    if columns is None:
        rows = []
        for table, members in TABLES.items():
            for column in members:
                row = derivation(column, organism)
                row["table"] = table
                rows.append(row)
        return rows
    return [derivation(column, organism) for column in columns]


def categories_table(columns=None, organism: str = organisms.TOXOPLASMA) -> pd.DataFrame:
    """The review table of instruction 73: one row per column, with its derivation and confidence."""
    frame = pd.DataFrame(derivations(columns, organism))
    order = [c for c in ("column", "table", "categories", "derived_from", "reason",
                         "provenance_key", "provenance_level_kind", "confidence") if c in frame.columns]
    return frame[order]


def unannotated(columns=None, organism: str = organisms.TOXOPLASMA) -> tuple:
    """Columns with no information category at all -- the gate the whole mechanism depends on.

    Empty for the shipped tables, and kept empty by `tests/test_information_categories.py`. A column
    here cannot be used as evidence by any question, because nothing can be said about whether using
    it would be circular.
    """
    members = SHIPPED_COLUMNS if columns is None else [str(c) for c in columns]
    return tuple(c for c in members if not categories_of(c, organism))


def excluded_for(label_category, columns=None, organism: str = organisms.TOXOPLASMA) -> tuple:
    """Every column that must be excluded when inferring a label of this category (or categories).

    The rule of instruction 73: a column is excluded if ANY of its categories matches the label's.
    `identity` columns are in the answer for every category, because they are never features at all.
    Unknown category names raise rather than quietly excluding nothing, which would be the silent
    failure the mechanism replaces.
    """
    wanted = (label_category,) if isinstance(label_category, str) else tuple(label_category)
    unknown = sorted(set(wanted) - set(CATEGORIES))
    if unknown:
        raise ValueError(f"not information categories: {unknown}; the vocabulary is {sorted(CATEGORIES)}")
    members = SHIPPED_COLUMNS if columns is None else [str(c) for c in columns]
    out = []
    for column in members:
        carried = set(categories_of(column, organism))
        if carried & set(wanted) or carried & set(NEVER_A_FEATURE):
            out.append(column)
    return tuple(out)


def is_feature(column: str, organism: str = organisms.TOXOPLASMA) -> bool:
    """Whether a column may be offered as a feature at all: annotated, and not an identity column."""
    carried = categories_of(column, organism)
    return bool(carried) and not set(carried) & set(NEVER_A_FEATURE)


def counts(columns=None, organism: str = organisms.TOXOPLASMA) -> pd.DataFrame:
    """How many columns carry each category, and how many of those were derived from provenance.

    A column with two categories counts once under each, so the total exceeds the column count; that
    is the honest way to report it, since the exclusion fires on any match. A column shipped in two
    tables (the imaging screen's phenotypes are in both `nodes` and the map) is counted once.
    """
    seen, rows = set(), []
    for row in derivations(columns, organism):
        if row["column"] in seen:
            continue
        seen.add(row["column"])
        rows.append(row)
    out = []
    for name in sorted(CATEGORIES):
        members = [r for r in rows if name in (r["categories"].split(",") if r["categories"] else [])]
        out.append({"category": name, "columns": len(members),
                    "from_provenance": sum(r["derived_from"] in FROM_PROVENANCE for r in members),
                    "judged": sum(r["derived_from"] not in FROM_PROVENANCE for r in members)})
    return pd.DataFrame(out)


def label_categories(question) -> tuple:
    """The categories of a question's own LABELS -- what `excluded_for` is then asked about.

    Taken from the question's declared `label_category` where it has one, and otherwise read off the
    column named in `label_source` (`'nodes.parquet:fit_no_glucose <= -2.0'` -> `fit_no_glucose`), so
    a question that forgot to declare one still gets the exclusion by construction. `()` when neither
    route resolves, in which case only the hand-written `evidence_forbidden` list applies and the run
    record says so.
    """
    declared = getattr(question, "label_category", "") or ""
    if declared:
        return tuple(part.strip() for part in declared.split(",") if part.strip())
    source = str(getattr(question, "label_source", "") or "")
    for token in re.split(r"[\s:,<>=()]+", source):
        if token and categories_of(token):
            return categories_of(token)
    return ()
