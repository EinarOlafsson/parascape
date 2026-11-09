"""Recorded source refusals scoped to a question, not an entire publication.

A publication can contain both useful evidence and an unsuitable table for one particular analysis. Recording the refusal at its original question and entity scope preserves the reason without removing unrelated valid evidence or pretending that unavailable observations are negative experimental results.
"""

REFUSED_CANDIDATES = {
    # Eight GEO series proposed for the IFN-gamma macrophage slot, opened 2026-08-19. Not one is a
    # transcriptome of the parasite inside an activated macrophage: the term that matched is in the
    # summaries, not in the samples. The slot's verdict stands.
    ("Tg_transcription · in IFN-gamma macrophage", "GSE313582"):
        "GCN5b/PHD1 knockdown transcriptome in RH Ku80 tachyzoites, HFF host, no macrophage",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE313273"):
        "ChIP-seq of GCN5b, PHD1, MORC, HDAC3 and histone marks; not a transcriptome at all",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE313048"):
        "ATAC-seq of GCN5b knockdown tachyzoites in HFF monolayers",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE334992"):
        "RNA-seq of purified EXTRACELLULAR tachyzoites, wild type against one knockout, no host",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE329845"):
        "TgPRO knockout against wild type under actinomycin D; a redox-adaptation experiment",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE300509"):
        "ME49 grown on different host cell LINES; the host axis is cell line, not activation state",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE275112"):
        "AP2XII-8 knockdown and CUT&Tag; a cell-cycle transcription factor in HFF",
    ("Tg_transcription · in IFN-gamma macrophage", "GSE266204"):
        "AP2XII-9 knockdown and CUT&Tag; likewise",
    # The six ribosome-profiling series proposed for the cell-cycle translation slot are the same
    # nine that the recorded sweep already opened one at a time. None is synchronised or sorted.
    ("Tg_translation · per cell-cycle phase", "GSE302108"):
        "5'UTR MPRA -- reporter constructs, not the endogenous translatome",
    ("Tg_translation · per cell-cycle phase", "GSE302107"):
        "ribosome profiling of unsynchronised tachyzoites; no cell-cycle axis",
    ("Tg_translation · per cell-cycle phase", "GSE243206"):
        "eIF4E1 depletion driving bradyzoite formation; a stage axis, not a cycle axis",
    ("Tg_translation · per cell-cycle phase", "GSE129869"):
        "host-context ribosome profiling, already shipped as its own dataset",
    ("Tg_translation · per cell-cycle phase", "GSE99395"):
        "intracellular against extracellular ribosome profiling, already shipped",
    ("Tg_translation · per cell-cycle phase", "GSE43722"):
        "the unfolded protein response; a stress axis",
    # Opened 2026-08-19 for `Pf_glycosylation`, and refused on three counts at once: the modification
    # is measured on RECOMBINANT protein, the substrate list is a domain prediction, and the biology
    # is the mosquito stages rather than the blood stage the slot names.
    ("Pf_glycosylation · asexual blood stage", "PXD033470"):
        "tryptophan C-mannosylation shown on recombinant SPATR and MTRAP, with the substrate set "
        "defined by TSR-domain content rather than measured in the parasite, and in the "
        "transmission stages rather than the asexual blood stage",
    # And the one proposed for protein turnover is the study that already fills a different slot.
    ("Tg_protein turnover", "PXD033642"):
        "the Ca2+-responsive thermal-shift proteome (PMID 35976251), which already fills "
        "`Tg_target engagement / thermal shift`; melting behaviour is not a degradation rate",
    # --- Verified 2026-11-09 against the GEO and PRIDE records themselves
    # (`results/slot_candidate_verification_2026_11_09.tsv`: 160 slot-candidate pairs over 26 slots and
    # 97 distinct accessions, every verdict carrying a verbatim quote from the record). 144 are the wrong
    # assay, organism or stage for the slot that proposed them, and are refused here so the literature
    # index cannot propose them again. The 8 admissible ones are NOT refused: see instruction 41.
    ("Pf_RNA-binding protein targets", "GSE267361"):
        "wrong assay: deposit is RNA-seq of PfMD3 knockout gametocytes; slot needs CLIP-seq or RIP-seq of an RNA-binding protein",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE136008"):
        "wrong assay: deposit is microarray of gametocytes +/- BIX-01294; slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE159516"):
        "wrong assay: deposit is transcript microarray time course, sub-lethal fosmidomycin; slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE208865"):
        "wrong assay: deposit is RNA-seq of PfNIF4 knockdown vs wild type; slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE212052"):
        "wrong assay: deposit is AP2-LT transcription-factor ChIP-seq; slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE215426"):
        "wrong assay: deposit is histone-PTM ChIP-seq, schizont/merozoite/ring; slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE215427"):
        "wrong assay: deposit is RNA-seq of schizonts vs merozoites; slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE24416"):
        "wrong assay: deposit is microarray transcript time course (Dd2 IDC); slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_S-nitrosylation · asexual blood stage", "GSE304394"):
        "wrong assay: deposit is CUT&RUN for protein lactylation (panKla, H4K12la); slot needs S-nitrosylation proteomics (e.g. biotin switch MS)",
    ("Pf_TF binding · per factor", "GSE295105"):
        "wrong assay: deposit is RNA-seq of pfap2-g reporter line (design text copied from the ChIP series); slot needs ChIP-seq or CUT&RUN of a transcription factor",
    ("Pf_TF binding · per factor", "GSE298024"):
        "wrong assay: deposit is H3K9me3 ChIP-seq, stress conditions; slot needs ChIP-seq or CUT&RUN of a transcription factor",
    ("Pf_TF binding · per factor", "GSE298061"):
        "wrong assay: deposit is RNA-seq, choline depletion / DHA / heat shock; slot needs ChIP-seq or CUT&RUN of a transcription factor",
    ("Pf_TF binding · per factor", "GSE299944"):
        "wrong assay: PfHP1 is a heterochromatin protein and the rest are histone marks, not a TF",
    ("Pf_antigenic variation family expression", "GSE186984"):
        "wrong assay: deposit is ChIP-seq of bromodomain proteins PfBDP7/PfBDP1; slot needs var/rif/stevor transcript expression per gene",
    ("Pf_antigenic variation family expression", "GSE296188"):
        "wrong assay: deposit is long-read DNA sequencing of subtelomeric deletions; slot needs var/rif/stevor transcript expression per gene",
    ("Pf_chromatin accessibility", "GSE239393"):
        "wrong assay: deposit is ChIP-seq of PfMORC; slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_chromatin accessibility", "GSE246114"):
        "wrong assay: deposit is ChIP-seq of PfMYST and H4 acetylation marks; slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_chromatin accessibility", "GSE288146"):
        "wrong assay: deposit is scRNA-seq of PfAP2-V knockdown; slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_chromatin accessibility", "GSE298024"):
        "wrong assay: deposit is H3K9me3 ChIP-seq, stress conditions; slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_chromatin accessibility", "GSE298061"):
        "wrong assay: deposit is RNA-seq, choline depletion / DHA / heat shock; slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_chromatin accessibility", "GSE304394"):
        "wrong assay: deposit is CUT&RUN for protein lactylation (panKla, H4K12la); slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_chromatin accessibility", "GSE337518"):
        "wrong assay: deposit is ChIP-seq of PfDBSP, a novel dsDNA-binding protein (1868 bound gene elements); slot needs ATAC-seq, DNase-seq, MNase-seq or FAIRE",
    ("Pf_co-fitness", "GSE103042"):
        "wrong assay: deposit is host+parasite RNA-seq, macaque bone marrow; slot needs co-fitness pairs from several genetic screens; organism also not accepted (P. coatneyi in M. mulatta)",
    ("Pf_co-fitness", "GSE103259"):
        "wrong assay: deposit is host+parasite RNA-seq, macaque whole blood; slot needs co-fitness pairs from several genetic screens; organism also not accepted (P. coatneyi in M. mulatta)",
    ("Pf_co-fitness", "GSE103507"):
        "wrong assay: deposit is host+parasite RNA-seq, macaque whole blood; slot needs co-fitness pairs from several genetic screens; organism also not accepted (P. cynomolgi in M. mulatta)",
    ("Pf_co-fitness", "GSE104223"):
        "wrong assay: deposit is host+parasite RNA-seq, macaque whole blood; slot needs co-fitness pairs from several genetic screens; organism also not accepted (P. cynomolgi in M. mulatta)",
    ("Pf_co-fitness", "GSE128115"):
        "wrong assay: deposit is host+parasite RNA-seq from macaque blood (MaHPIC); slot needs co-fitness pairs from several genetic screens",
    ("Pf_co-fitness", "GSE85478"):
        "wrong assay: deposit is nascent-RNA-seq and Pol II ChIP-seq, intraerythrocytic cycle; slot needs co-fitness pairs from several genetic screens",
    ("Pf_co-fitness", "GSE94273"):
        "wrong assay: deposit is host+parasite RNA-seq, macaque bone marrow; slot needs co-fitness pairs from several genetic screens; organism also not accepted (P. cynomolgi in M. mulatta)",
    ("Pf_co-fitness", "GSE99486"):
        "wrong assay: deposit is host+parasite RNA-seq, macaque whole blood; slot needs co-fitness pairs from several genetic screens; organism also not accepted (P. cynomolgi in M. mulatta)",
    ("Pf_essentiality in a second background", "GSE119516"):
        "wrong assay: deposit is microarray time course NF54 vs K13-promoter piggyBac mutant; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE188542"):
        "wrong assay: deposit is RNA-seq of one piggyBac mutant +/- DHA or bortezomib; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE205012"):
        "wrong assay: deposit is RNA-seq of a single piggyBac mutant (KIC5) vs NF54; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE214535"):
        "wrong assay: deposit is H3R2me2s CUT&Tag, schizont; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE223497"):
        "wrong assay: deposit is RNA-seq of a single piggyBac mutant (PF3D7_1136600) vs NF54; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE228947"):
        "wrong assay: deposit is MNase-seq nucleosome profiling of PfSnf2L-KO; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE228948"):
        "wrong assay: deposit is RNA-seq of PfSnf2L knockout / inhibitor, single-gene perturbation; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_essentiality in a second background", "GSE284527"):
        "wrong assay: deposit is H3K9ac / GCN5 ChIP-seq and CUT&Tag, schizont; slot needs knockout/piggyBac essentiality in a second strain",
    ("Pf_fitness · gametocyte", "GSE149772"):
        "wrong assay: deposit is RNA-seq of AP2-G5 / AP2-G4 mutant lines; slot needs pooled mutant or CRISPR/piggyBac fitness screen in gametocytes",
    ("Pf_fitness · gametocyte", "GSE149773"):
        "wrong assay: deposit is ChIP-seq of AP2-G5-GFP / AP2-G-GFP; slot needs pooled mutant or CRISPR/piggyBac fitness screen in gametocytes",
    ("Pf_fitness · gametocyte", "GSE157709"):
        "wrong assay: deposit is histone ChIP-seq (H3K9ac, H3K9me3) in AP2-G5 truncation line; slot needs pooled mutant or CRISPR/piggyBac fitness screen in gametocytes",
    ("Pf_fitness · gametocyte", "GSE175931"):
        "wrong assay: deposit is whole-genome sequencing of AP2-G5 mutant lines; slot needs pooled mutant or CRISPR/piggyBac fitness screen in gametocytes",
    ("Pf_fitness · gametocyte", "GSE228947"):
        "wrong assay: deposit is MNase-seq nucleosome profiling of PfSnf2L-KO; slot needs pooled mutant or CRISPR/piggyBac fitness screen in gametocytes",
    ("Pf_fitness · gametocyte", "GSE228948"):
        "wrong assay: deposit is RNA-seq of PfSnf2L knockout / inhibitor, single-gene perturbation; slot needs pooled mutant or CRISPR/piggyBac fitness screen in gametocytes",
    ("Pf_fitness · transmission", "GSE149772"):
        "wrong assay: deposit is RNA-seq of AP2-G5 / AP2-G4 mutant lines; slot needs mutant screen of mosquito transmission fitness",
    ("Pf_fitness · transmission", "GSE149773"):
        "wrong assay: deposit is ChIP-seq of AP2-G5-GFP / AP2-G-GFP; slot needs mutant screen of mosquito transmission fitness",
    ("Pf_fitness · transmission", "GSE157709"):
        "wrong assay: deposit is histone ChIP-seq (H3K9ac, H3K9me3) in AP2-G5 truncation line; slot needs mutant screen of mosquito transmission fitness",
    ("Pf_fitness · transmission", "GSE175931"):
        "wrong assay: deposit is whole-genome sequencing of AP2-G5 mutant lines; slot needs mutant screen of mosquito transmission fitness",
    ("Pf_glycosylation · asexual blood stage", "GSE136008"):
        "wrong assay: deposit is microarray of gametocytes +/- BIX-01294; slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE159516"):
        "wrong assay: deposit is transcript microarray time course, sub-lethal fosmidomycin; slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE208865"):
        "wrong assay: deposit is RNA-seq of PfNIF4 knockdown vs wild type; slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE212052"):
        "wrong assay: deposit is AP2-LT transcription-factor ChIP-seq; slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE215426"):
        "wrong assay: deposit is histone-PTM ChIP-seq, schizont/merozoite/ring; slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE215427"):
        "wrong assay: deposit is RNA-seq of schizonts vs merozoites; slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE24416"):
        "wrong assay: deposit is microarray transcript time course (Dd2 IDC); slot needs glycoproteomics / glycan-site MS",
    ("Pf_glycosylation · asexual blood stage", "GSE304394"):
        "wrong assay: deposit is CUT&RUN for protein lactylation (panKla, H4K12la); slot needs glycoproteomics / glycan-site MS",
    ("Pf_host receptor binding", "GSE159124"):
        "wrong assay: deposit is RNA-seq of schizont and gametocyte stages; slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE159126"):
        "wrong assay: deposit is RNA-BisSeq (m5C); slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE166258"):
        "wrong assay: deposit is transcript microarray of CHMI-derived parasites; slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE166390"):
        "wrong assay: deposit is H3K9me3 ChIP-seq; slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE180985"):
        "wrong assay: deposit is RNA-seq of sorted male/female gametocytes; slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE186820"):
        "wrong assay: deposit is RNA-seq of parasite isolates from cerebral vs uncomplicated malaria; slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE274138"):
        "wrong assay: deposit is scRNA-seq of P. falciparum 3D7 +/- beta-hydroxybutyrate; slot needs host receptor binding assay per parasite gene",
    ("Pf_host receptor binding", "GSE301430"):
        "wrong assay: deposit is RNA-seq of human monocytes exposed to parasites; slot needs host receptor binding assay per parasite gene; organism also not accepted (P. falciparum with human monocytes)",
    ("Pf_host transcriptional effect per effector", "GSE78931"):
        "wrong assay: deposit is dual RNA-seq of infected vs uninfected hepatoma cells, no per-effector perturbation; slot needs dual perturb-seq: host transcript change per knocked-out effector",
    ("Pf_interaction degree · crosslink MS", "GSE208757"):
        "wrong assay: deposit is ORC1 ChIP-seq (replication origins); slot needs crosslinking mass spectrometry",
    ("Pf_interaction degree · crosslink MS", "GSE278140"):
        "wrong assay: deposit is ChIP-seq of AP2-P and MORC (Micro-C companion); slot needs crosslinking mass spectrometry",
    ("Pf_interaction degree · crosslink MS", "GSE85478"):
        "wrong assay: deposit is nascent-RNA-seq and Pol II ChIP-seq, intraerythrocytic cycle; slot needs crosslinking mass spectrometry",
    ("Pf_interaction · proximity labelling", "GSE208865"):
        "wrong assay: deposit is RNA-seq of PfNIF4 knockdown vs wild type; slot needs BioID / TurboID proximity proteomics",
    ("Pf_interaction · proximity labelling", "GSE231838"):
        "wrong assay: deposit is RNA-seq of female gametocytes vs zygotes (TurboID proteomics is not in this GEO record); slot needs BioID / TurboID proximity proteomics",
    ("Pf_interaction · proximity labelling", "GSE267413"):
        "wrong assay: deposit is bulk RNA-seq of var-gene SLI lines (not proximity labelling); slot needs BioID / TurboID proximity proteomics",
    ("Pf_invasion and egress phenotype", "GSE164175"):
        "wrong assay: deposit is RNA-seq of PP1 knockdown gametocytes/schizonts; slot needs high-content imaging of invasion/egress per gene",
    ("Pf_invasion and egress phenotype", "GSE180257"):
        "wrong assay: deposit is scRNA-seq of P. cynomolgi liver stages; slot needs high-content imaging of invasion/egress per gene; organism also not accepted (P. cynomolgi (+ mouse))",
    ("Pf_invasion and egress phenotype", "GSE180258"):
        "wrong assay: deposit is SuperSeries (ChIP-seq + scRNA-seq) of P. cynomolgi hypnozoites; slot needs high-content imaging of invasion/egress per gene; organism also not accepted (P. cynomolgi (+ mouse))",
    ("Pf_invasion and egress phenotype", "GSE291233"):
        "wrong assay: deposit is RNA-seq after chloroquine / MMV compounds; slot needs high-content imaging of invasion/egress per gene",
    ("Pf_invasion and egress phenotype", "GSE317694"):
        "wrong assay: deposit is HP1 CUT&Tag after inhibitor treatment (microscopy only used to pick compounds); slot needs high-content imaging of invasion/egress per gene",
    ("Pf_lipid composition", "GSE159516"):
        "wrong assay: deposit is transcript microarray time course, sub-lethal fosmidomycin; slot needs lipidomics",
    ("Pf_lipid composition", "GSE176469"):
        "wrong assay: deposit is transcript microarray time course, sub-lethal pyrazoleamide; slot needs lipidomics",
    ("Pf_lipid composition", "GSE197916"):
        "wrong assay: deposit is CUT&RUN for H3K9me3 and H3K4me3 under media conditions; slot needs lipidomics",
    ("Pf_lipid composition", "GSE218998"):
        "wrong assay: deposit is transcript microarray time course, KAE609 / pyrazoleamide; slot needs lipidomics",
    ("Pf_metabolic flux", "GSE13578"):
        "wrong assay: deposit is spotted oligo microarray, PfAdoMetDC/ODC co-inhibition; slot needs isotope-labelled precursor flux measurements",
    ("Pf_metabolic flux", "GSE136169"):
        "wrong assay: deposit is transcript microarray time course of apicoplast-negative PfMev parasites; slot needs isotope-labelled precursor flux measurements",
    ("Pf_metabolic flux", "GSE176469"):
        "wrong assay: deposit is transcript microarray time course, sub-lethal pyrazoleamide; slot needs isotope-labelled precursor flux measurements",
    ("Pf_metabolite levels", "GSE123495"):
        "wrong assay: deposit is host+parasite RNA-seq from macaque pilot infection; slot needs metabolomics (LC/GC-MS) of metabolite levels",
    ("Pf_metabolite levels", "GSE127079"):
        "wrong assay: deposit is host+parasite RNA-seq from macaque infection; slot needs metabolomics (LC/GC-MS) of metabolite levels",
    ("Pf_metabolite levels", "GSE128115"):
        "wrong assay: deposit is host+parasite RNA-seq from macaque blood (MaHPIC); slot needs metabolomics (LC/GC-MS) of metabolite levels",
    ("Pf_metabolite levels", "GSE94274"):
        "wrong assay: deposit is SuperSeries of MaHPIC macaque RNA-seq; slot needs metabolomics (LC/GC-MS) of metabolite levels; organism also not accepted (mixed macaque/Plasmodium (knowlesi, coatneyi, cynomolgi))",
    ("Pf_protein abundance · other life stages", "GSE111333"):
        "wrong assay: deposit is RNA-seq of P. berghei AP2 SNP mutant vs wild type; slot needs proteomics (LC-MS) of non-asexual stages",
    ("Pf_protein abundance · other life stages", "GSE116219"):
        "wrong assay: deposit is ChIP-seq of SMC3; chromatin proteomics data not in GEO record",
    ("Pf_protein abundance · other life stages", "GSE123839"):
        "wrong assay: deposit is m6A-seq + RNA-seq, 12/24/36 hpi; slot needs proteomics (LC-MS) of non-asexual stages",
    ("Pf_protein abundance · other life stages", "GSE128115"):
        "wrong assay: deposit is host+parasite RNA-seq from macaque blood (MaHPIC); slot needs proteomics (LC-MS) of non-asexual stages",
    ("Pf_protein abundance · other life stages", "GSE152608"):
        "wrong context: GEO deposit is RNA-seq only and stage is asexual; record says: RNA sequencing was performed on circulating ring stage parasites and LC-MS/MS on the corresponding sequestered mature forms.",
    ("Pf_protein abundance · other life stages", "GSE164175"):
        "wrong assay: deposit is RNA-seq of PP1 knockdown gametocytes/schizonts; slot needs proteomics (LC-MS) of non-asexual stages",
    ("Pf_protein abundance · other life stages", "GSE208865"):
        "wrong assay: deposit is RNA-seq of PfNIF4 knockdown vs wild type; slot needs proteomics (LC-MS) of non-asexual stages",
    ("Pf_protein abundance · other life stages", "GSE267361"):
        "wrong assay: deposit is RNA-seq of PfMD3 knockout gametocytes; slot needs proteomics (LC-MS) of non-asexual stages",
    ("Pf_protein abundance · other life stages", "PXD046874"):
        "wrong assay: phosphosite meta-analysis, not abundance; record says: Phosphorylation in the Plasmodium falciparum proteome: A meta-analysis of publicly available data sets",
    ("Pf_protein abundance · other life stages", "PXD048354"):
        "wrong context: real Pf proteome but ring stage under DHA/DMSO, i.e. asexual not other stages",
    ("Pf_protein abundance · other life stages", "PXD060218"):
        "wrong context: measures plasma proteins, not parasite proteome in a parasite life stage",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE119516"):
        "wrong assay: deposit is microarray time course NF54 vs K13-promoter piggyBac mutant; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE157709"):
        "wrong assay: deposit is histone ChIP-seq (H3K9ac, H3K9me3) in AP2-G5 truncation line; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE175931"):
        "wrong assay: deposit is whole-genome sequencing of AP2-G5 mutant lines; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE188542"):
        "wrong assay: deposit is RNA-seq of one piggyBac mutant +/- DHA or bortezomib; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE205012"):
        "wrong assay: deposit is RNA-seq of a single piggyBac mutant (KIC5) vs NF54; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE223497"):
        "wrong assay: deposit is RNA-seq of a single piggyBac mutant (PF3D7_1136600) vs NF54; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE228947"):
        "wrong assay: deposit is MNase-seq nucleosome profiling of PfSnf2L-KO; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · asexual blood stage", "GSE228948"):
        "wrong assay: deposit is RNA-seq of PfSnf2L knockout / inhibitor, single-gene perturbation; slot needs per-gene sexual-conversion rate from a mutant screen",
    ("Pf_stage-conversion phenotype · gametocyte", "GSE149772"):
        "wrong assay: deposit is RNA-seq of AP2-G5 / AP2-G4 mutant lines; slot needs per-gene gametocyte-conversion phenotype",
    ("Pf_stage-conversion phenotype · gametocyte", "GSE149773"):
        "wrong assay: deposit is ChIP-seq of AP2-G5-GFP / AP2-G-GFP; slot needs per-gene gametocyte-conversion phenotype",
    ("Pf_stage-conversion phenotype · gametocyte", "GSE157709"):
        "wrong assay: deposit is histone ChIP-seq (H3K9ac, H3K9me3) in AP2-G5 truncation line; slot needs per-gene gametocyte-conversion phenotype",
    ("Pf_stage-conversion phenotype · gametocyte", "GSE175931"):
        "wrong assay: deposit is whole-genome sequencing of AP2-G5 mutant lines; slot needs per-gene gametocyte-conversion phenotype",
    ("Pf_stage-conversion phenotype · gametocyte", "GSE228947"):
        "wrong assay: deposit is MNase-seq nucleosome profiling of PfSnf2L-KO; slot needs per-gene gametocyte-conversion phenotype",
    ("Pf_stage-conversion phenotype · gametocyte", "GSE228948"):
        "wrong assay: deposit is RNA-seq of PfSnf2L knockout / inhibitor, single-gene perturbation; slot needs per-gene gametocyte-conversion phenotype",
    ("Pf_stage-conversion phenotype · liver stage", "GSE119516"):
        "wrong assay: deposit is microarray time course NF54 vs K13-promoter piggyBac mutant; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE188542"):
        "wrong assay: deposit is RNA-seq of one piggyBac mutant +/- DHA or bortezomib; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE205012"):
        "wrong assay: deposit is RNA-seq of a single piggyBac mutant (KIC5) vs NF54; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE214535"):
        "wrong assay: deposit is H3R2me2s CUT&Tag, schizont; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE223497"):
        "wrong assay: deposit is RNA-seq of a single piggyBac mutant (PF3D7_1136600) vs NF54; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE228947"):
        "wrong assay: deposit is MNase-seq nucleosome profiling of PfSnf2L-KO; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE228948"):
        "wrong assay: deposit is RNA-seq of PfSnf2L knockout / inhibitor, single-gene perturbation; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_stage-conversion phenotype · liver stage", "GSE284527"):
        "wrong assay: deposit is H3K9ac / GCN5 ChIP-seq and CUT&Tag, schizont; slot needs per-gene stage-conversion phenotype in liver stage",
    ("Pf_transcription · liver stage", "GSE180256"):
        "wrong organism: ChIP-seq and P. cynomolgi; record says: ChIP-seq analysis of P. cynomolgy sprozoites and blood stage parasites for one heterochromatic (H3K9me2) and one euchromatic (H3K9/14ac) mark.",
    ("Pf_transcription · liver stage", "GSE180257"):
        "wrong organism: P. cynomolgi is not in the accepted Plasmodium list; record says: scRNA-seq",
    ("Pf_transcription · liver stage", "GSE180258"):
        "wrong organism: P. cynomolgi is not in the accepted Plasmodium list; record says: This SuperSeries is composed of the SubSeries listed below.",
    ("Pf_transcription · liver stage", "GSE215426"):
        "wrong assay: ChIP-seq, and mostly blood stage; record says: Chromatin immunoprecipitation DNA-sequencing (ChIP-seq) for histone PTMs in Plasmodium falciparum schizonts, merozoites, and rings",
    ("Pf_transcription · liver stage", "GSE215429"):
        "wrong context: blood stages; only P. berghei hepatic merozoites ChIP are liver",
    ("Pf_transcription · mosquito stages", "GSE180256"):
        "wrong organism: P. cynomolgi and ChIP-seq; record says: ChIP-seq analysis of P. cynomolgy sprozoites and blood stage parasites for one heterochromatic (H3K9me2) and one euchromatic (H3K9/14ac) mark.",
    ("Pf_transcription · mosquito stages", "GSE180257"):
        "wrong organism: P. cynomolgi; record says: scRNA-seq",
    ("Pf_transcription · mosquito stages", "GSE180258"):
        "wrong organism: P. cynomolgi and not mosquito stages; record says: This SuperSeries is composed of the SubSeries listed below.",
    ("Pf_transcription · mosquito stages", "GSE269994"):
        "wrong context: P. berghei knockdown gametocytes; record says: About 108 purified activated gametocytes from PbIMC1gcKD [+]/[-] TMP parasites were lysed",
    ("Pf_transcription · mosquito stages", "GSE271105"):
        "wrong context: P. berghei knockout activated gametocytes; record says: and the activated gametocytes (25℃ for 15 min) were collected",
    ("Pf_transcription · mosquito stages", "GSE284425"):
        "wrong organism: mosquito transcriptome record; parasite genes not described; record says: mRNA-seq circadian profiles generated for female Anopheles stephensi mosquitoes kept under 2 different light conditions.",
    ("Pf_translation efficiency · asexual blood stage", "GSE1124"):
        "wrong assay: deposit is human whole-blood microarray of malaria patients; slot needs Ribo-seq plus matched RNA; organism also not accepted (H. sapiens)",
    ("Pf_translation efficiency · asexual blood stage", "GSE117613"):
        "wrong assay: deposit is human whole-blood microarray, CM vs SMA; slot needs Ribo-seq plus matched RNA; organism also not accepted (H. sapiens)",
    ("Pf_translation efficiency · asexual blood stage", "GSE34404"):
        "wrong assay: deposit is human whole-blood expression + genotype eQTL; slot needs Ribo-seq plus matched RNA; organism also not accepted (H. sapiens)",
    ("Pf_translation efficiency · asexual blood stage", "GSE35858"):
        "wrong assay: deposit is human whole-blood microarray of Vietnamese patients; slot needs Ribo-seq plus matched RNA; organism also not accepted (H. sapiens)",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE136008"):
        "wrong assay: deposit is microarray of gametocytes +/- BIX-01294; slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE159516"):
        "wrong assay: deposit is transcript microarray time course, sub-lethal fosmidomycin; slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE208865"):
        "wrong assay: deposit is RNA-seq of PfNIF4 knockdown vs wild type; slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE212052"):
        "wrong assay: deposit is AP2-LT transcription-factor ChIP-seq; slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE215426"):
        "wrong assay: deposit is histone-PTM ChIP-seq, schizont/merozoite/ring; slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE215427"):
        "wrong assay: deposit is RNA-seq of schizonts vs merozoites; slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE24416"):
        "wrong assay: deposit is microarray transcript time course (Dd2 IDC); slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "GSE304394"):
        "wrong assay: deposit is CUT&RUN for protein lactylation (panKla, H4K12la); slot needs ubiquitin/SUMO-enrichment proteomics",
    ("Pf_ubiquitination / SUMOylation · asexual blood stage", "PXD042503"):
        "wrong assay: deposit is global LC-MS of clinical severe-malaria isolates, no PTM enrichment stated; slot needs ubiquitin/SUMO-enrichment proteomics",
    # --- 2026-11-09, refused on the DEPOSITED DATA rather than on the record (instruction 41). The
    # GEO record reads exactly like a usable ChIP-seq of a transcription factor, and the verification
    # pass judged it admissible on that basis. The peak files say otherwise: replicate 1 holds 25
    # peaks and replicate 2 holds 1,727 (19 against 1,680 in the choline-depleted arm), so the
    # replicates do not agree on what was bound. Assigning peaks to genes and keeping only those
    # reproducible in BOTH replicates leaves 18 genes with choline, 12 without, 21 in all -- a short
    # candidate list, not a per-gene binding map, and for the other 5,699 genes "not bound" would be
    # unmeasured rather than measured-as-negative. Gene coordinates from Ensembl Protists r63.
    ("Pf_TF binding \u00b7 per factor", "GSE339382"):
        "replicates disagree: 25 peaks against 1,727; only 21 genes reproduce in both, so this is a candidate list and not a per-gene binding map",
}


def records(catalog=None):
    """Return the existing scoped refusals with explicit organism and entity unit."""
    from . import slots
    catalog = slots.all_slots() if catalog is None else catalog
    questions = {slot.organism + '_' + slot.name: slot for slot in catalog}
    result = []
    for (question, accession), reason in sorted(REFUSED_CANDIDATES.items()):
        if question not in questions:
            continue
        slot = questions[question]
        result.append({'source_id': 'candidate:' + accession, 'organism': slot.organism,
            'unit': slot.unit, 'question': question, 'reason': reason})
    return result
