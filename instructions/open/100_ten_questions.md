================================================================================
TEN QUESTIONS THE DATA ALREADY IN PARASCAPE CAN ANSWER
================================================================================

STATE
    The bradyzoite analysis is the worked example and it is now a module: five
    single-cell subtypes of a mouse brain cyst (`brady_subtypes`) met with the
    201 stored within-study contrasts, asking which subtype a treatment reaches.
    Engine `parascape/cluster_axis.py`, module `parascape/reach_module.py`,
    commits 055155ea, 5c852a76, 300cac69.

    What made that question worth asking was not the single-cell data and not
    the drug data. It was that the two were never put together, and that the
    combination answers something the field needs: how do we treat an infection
    we cannot currently clear.

    The ten below are chosen the same way. Each one combines datasets that are
    in the release and are not currently read together, each has an answer a
    reviewer could argue with, and each has a stated reason it could come out
    empty. Ranked by what the answer would be worth, not by how easy it is.

    The columns named are real. Dataset keys are `parascape/datasets.py`.

    THE WORKED EXAMPLE, CHECKED (2026-11-18). The five subtypes were run against
    all 201 stored contrasts, 1,000 matched sets each, 150 markers per subtype.
    135 contrasts reach at least one subtype, and the hits land where they should
    if the mechanism works: 56 of them in GSE200962 and 17 in GSE132238, the
    bradyzoite differentiation studies, with single hits in drug contrasts such as
    the MC1742 treatment of GSE175919.

    The result worth following: subtype **C is reached by 4 of the 201 contrasts**,
    against 83 for B and 86 for E. That is not the matching giving up -- C is
    testable and its matched sets are centred in all 40 contrasts checked for it --
    so C is the subtype the differentiation and drug contrasts in the literature do
    not move. If a compound has to reach every subtype to clear a cyst, C is the
    one nothing in this catalogue reaches.

WHY IT MATTERS
    ParaScape holds 177 datasets and the overwhelming majority are read one at a
    time. A dataset read alone answers the question its authors already answered.
    The value of putting them in one gene space is the question nobody could ask
    with either half, and that value is only realised if somebody writes the
    questions down.

--------------------------------------------------------------------------------
1.  WHICH HOST CELL CAN AN EFFECTOR ACTUALLY WORK IN?
--------------------------------------------------------------------------------
    An exported effector that binds a host protein the host cell does not make
    cannot do anything in that cell. Toxoplasma infects neurons, macrophages,
    fibroblasts and muscle, and the effector repertoire is one repertoire.

    Combine   `starpath_xlms` and `pf_host_bridge_xlms` (parasite-host crosslinks,
              residue level) and `hosttx_effectors`, against host expression per
              cell type: `host_immune_cell_rna`, `host_brain_cell_rna`,
              `host_tissue_proteome`, `host_monocyte_proteome`,
              `host_fibroblast_proteome`, `host_macrophage_surfaceome`.
    Answer    an effector x host-cell-type table saying, for each effector, which
              cell types make its partner and which do not. Then: is any effector
              class restricted to one cell type, and does that match where the
              parasite forms cysts?
    Could be  empty if every partner is a housekeeping protein made everywhere,
              which would itself be the finding -- the effectors target what is
              always there, and tissue tropism is not explained by the partners.
    Honest    a protein absent from a proteome is a detection limit, not a zero.
              Coverage is 9.8% for bone marrow and 22.5% for dermis, so most
              cells of this table must read n.m. and the A grades on those
              datasets overstate what they cover (instruction 65).
    Today     the data is loaded. The join is not written.

--------------------------------------------------------------------------------
2.  IS A STAGE-RESTRICTED DRUG EXPLAINED BY WHEN ITS TARGET IS MADE?
--------------------------------------------------------------------------------
    Nothing clears a cyst. The usual explanation is that bradyzoites are dormant,
    which is a statement about the parasite. The testable alternative is about
    the drug: its target is not made in that stage.

    Combine   `pf_target_engagement`, `pf_meltome` and `pf_hsp90_chemoproteome`
              (which protein a compound actually engages) against stage
              expression: `plasmodb_pf3d7_expression`, `pf_derived_stage_labels`,
              `pf_idc_timing`, and on the Toxoplasma side `stage_enriched` with
              `brady_subtypes`.
    Answer    for each compound with an engaged target, the stage profile of that
              target, and whether it predicts the stage the compound is known to
              act in. Calibrate on compounds whose stage specificity is known,
              then read the predictions for the ones where it is not.
    Could be  empty if engagement is promiscuous, in which case the question
              becomes which compounds have a single engaged target at all.
    Honest    an engaged target is not necessarily the killing target. The claim
              is association, and the scoreboard holds it to that.
    Today     needs one new mechanism: compound -> engaged protein -> stage.

--------------------------------------------------------------------------------
3.  THE CONDITIONAL ESSENTIALOME: WHAT IS ONLY ESSENTIAL SOMEWHERE?
--------------------------------------------------------------------------------
    A gene dispensable in a dish and essential inside a macrophage is a drug
    target the standard screen threw away. Seven screens in the release differ in
    exactly one condition each, which is the design this question needs.

    Combine   `crispr_invitro`, `crispr_invivo_composite`, `crispr_macrophage`,
              `crispr_serum_restriction`, `crispr_glucose_limitation`,
              `crispr_parasite_density`, `crispr_young2019`.
    Answer    per gene, which conditions it is required in, as a profile rather
              than a single fitness score; then the genes whose profile has one
              entry, which is the shortlist.
    Could be  dominated by noise: the difference of two noisy scores is noisier
              than either. This is the question that most needs the matched null
              and the smallest-detectable-difference machinery already built, and
              it must report how small a difference each pair of screens could
              have resolved.
    Honest    screens differ in library, passage and laboratory as well as in
              condition. A difference between two screens is a difference between
              two experiments until a shared control says otherwise.
    Today     the `reach` engine does most of it: each condition is an axis, and
              a gene set is a cluster. A profile view is the missing piece.

--------------------------------------------------------------------------------
4.  DOES THE PARASITE MOVE ITS PROTEINS, OR ONLY ITS TRANSCRIPTS?
--------------------------------------------------------------------------------
    Every stage conclusion in the field rests on transcription. A protein can
    change compartment without changing abundance, and that would be invisible to
    all of it.

    Combine   `lopit_tgon` (hyperLOPIT, protein compartment), `lopit_pfal`,
              `lopit_cpar` against the stage and stress transcriptomes
              (`stage_enriched`, `gse132248_stress`, `morc_depletion`,
              `brady_subtypes`) and the proteomes (`proteome_total`,
              `proteome_pru`, `iron_depletion_proteome`, `oocyst_itraq`).
    Answer    which compartments change their membership between stages, and
              which genes carry that change.
    Could be  blocked rather than empty, which is the real finding here: hyperLOPIT
              is tachyzoite-only, so for most genes this question cannot be
              answered at all. Writing down which genes those are is a declared
              gap and the strongest case in the release for a new experiment.
    Honest    a compartment call has a confidence; a reassignment within the
              confidence of both calls is not a move.
    Today     the data is loaded; the comparison is not written.

--------------------------------------------------------------------------------
5.  WHICH ENZYME HAS NO HOST COUNTERPART AND IS MADE WHEN WE NEED TO KILL?
--------------------------------------------------------------------------------
    The classical drug-target argument, done once across every dataset instead of
    one pathway at a time.

    Combine   `toxodb_ec_numbers` and `pf_enzyme_classification` (what reaction),
              `orthomcl` (is there a human orthologue), `host_gtex_transcriptome`
              (is the human counterpart expressed in the infected tissue),
              `foldseek_struct` and `pf_foldseek_struct` (is the fold itself
              distinctive), stage expression, and `crispr_invitro` plus
              `crispr_invivo_composite` (is it required), with `metabolome_iron`
              and `lipidome_vesicles` for whether the pathway runs.
    Answer    a ranked shortlist with each criterion shown separately, so a
              reader can see which criterion each candidate fails.
    Could be  short, and that is useful: a shortlist of four is a result.
    Honest    no orthologue is not no counterpart -- a different fold can catalyse
              the same reaction. The Foldseek screens already showed that
              host-facing domains are structurally NOVEL rather than mimics, so
              fold distinctiveness is evidence about druggability, not about
              mimicry.
    Today     every column is present. Nothing joins them.

--------------------------------------------------------------------------------
6.  DOES A RESISTANCE MUTATION SIT ON AN INTERFACE?
--------------------------------------------------------------------------------
    If resistance works by breaking a contact, the mutated residues should land
    where the protein touches something. Crosslinks give residues, so this is
    answerable at residue level rather than at gene level.

    Combine   `pf_resistome` and `plasmodb_pf3d7_snps` (which residues change
              under drug pressure in the field), `pf_ip_ms`, `pf_host_bridge_xlms`
              and `starpath_xlms` (which residues touch what),
              `pf_alphafold_confidence` and `local_af3` (where those residues are
              in the fold), `pf_cdpk1_dependent_sites` and `phosphosites` (whether
              the residue is modified).
    Answer    per resistance locus, whether its mutated residues are enriched at
              crosslinked interfaces or at modified sites, against a null of
              residues matched for solvent exposure and conservation.
    Could be  empty, and the matched null is what makes that informative: surface
              residues are both more mutable and more crosslinkable, so an
              unmatched comparison would find an effect whatever the truth.
    Honest    the PTM coordinate trap applies -- ToxoDB `modification_site` is the
              PEPTIDE start, not the residue. Naive use scores below chance.
    Today     needs a residue-level join, which is new; the gene-level one exists.

--------------------------------------------------------------------------------
7.  WHICH EXPORTED PROTEINS ARE BOTH SEEN BY IMMUNITY AND INVARIANT?
--------------------------------------------------------------------------------
    A vaccine target has to be visible and has to not change. Both halves are
    measured and they are measured in different datasets.

    Combine   `pf_iedb_tcell`, `pf_iedb_bcell`, `iedb_bcell`, `toxodb_epitopes`
              (what immunity sees), `plasmodb_pf3d7_snps` and `pf_resistome`
              (what varies in the field), `plasmodb_pf3d7_exportpred` and
              `pf_secretome` (what is exported), `pf_gametocyte_proteome` and
              stage expression (when it is made).
    Answer    the exported, epitope-dense, field-invariant set, with the usual
              three outcomes: in it, out of it, and not assayed for epitopes.
    Could be  nearly empty, because selection is the reason epitope-dense proteins
              vary. An empty answer here is a real statement about vaccine design.
    Honest    epitope coverage is heavily biased towards proteins somebody already
              studied. Absence of epitopes is absence of assay, and the matched
              null has to control for how much a protein has been looked at --
              `pf_literature` and `pubmed` are in the release for exactly this.
    Today     joinable now; the literature-matched null is the work.

--------------------------------------------------------------------------------
8.  IS THE STAGE SWITCH TRANSLATIONAL?
--------------------------------------------------------------------------------
    If protein output changes where mRNA does not, the switch is not
    transcriptional, and every transcriptome-based stage marker is a proxy for
    something else.

    Combine   `pf_riboseq`, `gse302107_riboseq`, `pf_translation_efficiency`,
              `pf_mrna_dynamics` and `gse99395`/`gse129869` against the matched
              transcript levels, with `pf_m6a_nanopore` (is the message marked),
              `pf_rna_dependence` and `gse223620_bfd2_rip` (is it bound),
              `toxodb_codon_usage` and `pf_codon_usage` (is it decoded slowly).
    Answer    the genes whose translation efficiency changes with stage while
              their transcript does not, and whether m6A or RNA binding marks
              them more than matched genes.
    Could be  empty if translation efficiency mostly tracks transcript level,
              which would close the question cleanly.
    Honest    ribosome profiling and RNA-seq come from different experiments; a
              ratio of two studies is a ratio of two batches. Only within-study
              contrasts count, which is the rule the contrast store already keeps.
    Today     the contrast store holds the within-study pairs. The join is new.

--------------------------------------------------------------------------------
9.  DOES THE PARASITE REWIRE THE SAME HOST PROGRAMME IN EVERY TISSUE?
--------------------------------------------------------------------------------
    The host side of question 1, and the same shape as the bradyzoite analysis
    with the roles swapped: the clusters are host cell programmes and the axis is
    infection.

    Combine   `host_hff_tg_infection` (fibroblast, Toxoplasma),
              `host_hepatocyte_pf_infection` (hepatocyte, Plasmodium),
              `host_bmdm_baseline` and `host_macrophage_surfaceome`,
              `host_k562_rhoptry_screen`, `pv_host_uptake`, with
              `host_gtex_transcriptome` and `host_immune_cell_rna` for what each
              cell type can express at all.
    Answer    which host programmes move on infection in one cell type and not in
              another, each with the smallest effect that comparison could have
              resolved.
    Could be  confounded by species and parasite as well as by tissue; the
              comparison is only honest where the same parasite is read in two
              host cells, and the write-up has to say which pairs those are.
    Honest    a host programme is a gene set somebody defined. The matched null
              must be matched on expression in that cell type, or the answer is
              about how much each cell type expresses.
    Today     runnable with the `reach` engine as built, once host programmes are
              available as gene sets.

--------------------------------------------------------------------------------
10. WHAT DOES AN UNANNOTATED GENE INHERIT FROM THE COMPLEX IT IS IN?
--------------------------------------------------------------------------------
    Roughly half the genome has no function recorded. Physical membership of a
    complex whose other members agree on a function is the strongest evidence
    available for the ones that have nothing.

    Combine   `starpath_xlms`, `ipms_corpus`, `bioid_corpus`, `pulldown_corpus`,
              `cotranslation_edges` (translated together), `pf_rna_dependence`,
              `lopit_tgon` (same compartment), with `interpro` and `orthomcl` for
              what is already known.
    Answer    for each unannotated gene, the function its partners agree on and
              how much they agree, scored on genes whose annotation was held out.
    Could be  bounded by the pulldown background problem: a bait that pulls
              everything agrees with everything. Instruction 80 is this, and its
              AUC fell from 0.710 to 0.684 as the project count doubled, so the
              honest version reports the gate and the trend.
    Honest    clustering must never see the genes being scored; co-translation and
              crosslinking are different evidence and must not be summed into one
              score as though they were replicates.
    Today     the parts exist across `infer`, `predict` and the scoreboard; what
              is missing is the one answer per gene with its uncertainty.

--------------------------------------------------------------------------------
WHAT TO DO
--------------------------------------------------------------------------------
    Nothing, until Einar picks. Each of the ten is between a day and a week, and
    building the wrong three is worse than building one.

    If a default is wanted: 3 and 9 are nearly free because the `reach` engine
    already does the shape of them, and 1 and 5 are the ones whose answer another
    laboratory would use.

HOW TO KNOW IT WORKED
    A question is answered when its module reports all three outcomes -- found,
    looked for and not found with the limit of detection, and not measured --
    and when a held-out check says how often an answer of that shape is right.
    A ranked list with no null behind it is not an answer.

TRAPS
    Every one of the ten can be made to produce a result by dropping the matched
    null, and every one of those results would be about expression level, assay
    coverage or category size. The null is the work.
