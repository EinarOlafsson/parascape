================================================================================
THE PARASCAPE PAPER: TWELVE FIGURES, AND SHAPING THE PROGRAM AROUND THEM
================================================================================

STATE
    Einar decided on 2026-11-18 to plan a ParaScape paper first, then change
    ParaScape to fit it. Three pieces of work, in this order and not any other:

        1. define the figure structure        -- DONE, twelve figures below
        2. produce the figures and the text   -- open, figure by figure
        3. reshape the program around them    -- open, and conditional on 2

    The third is the one that is easy to start too early. A module built for a
    figure that has not been drawn is a module built for a guess.

    Instruction 100 (ten questions) is SUPERSEDED by this file. What survived of
    it is recorded under "What the triage of instruction 100 established".

WHY IT MATTERS
    ParaScape is 177 datasets, 362 slots and twelve shell modules with no single
    statement of what it is FOR. A paper forces that statement, and the figure
    list is the statement: if a figure cannot be drawn from the program, either
    the figure is wrong or the program is missing something, and either way it
    is cheaper to find out on paper.

    The second reason is Einar's own: every analysis in the paper should be
    reproducible inside ParaScape by a user who did not write it. That is the
    same requirement spaCR carries, and it is what turns a paper into a tool.

--------------------------------------------------------------------------------
THE TWELVE FIGURES
--------------------------------------------------------------------------------
    Each one is titled by the CONCLUSION it reaches, not by the data it uses.
    A figure organised by data type is a chapter of a database description.

    State: [ ] not started · [~] partly built · [x] drawable today

    1.  [~] WHAT PARASCAPE IS. Imports, modules, the one gene space.
            Ends on: the corpus, assembled.

    2.  [ ] WHICH STUDIES CAN BE COMPARED AT ALL.
            Einar 2026-11-18: cluster the studies by how well their CONTROLS
            (untreated, wild type) agree; then add normalisation steps that make
            biological sense -- cell cycle stage first -- and ask which studies
            move into agreement. Conditions inside an agreeing cluster are a
            fair comparison; conditions across clusters are not.
            Ends on: which biological variable, removed, makes the literature
            agree -- the claim being that cell cycle is the dominant confounder
            in the published Toxoplasma transcriptome.
            Already working in miniature: the contrast ladder took two alkaline
            studies from agreeing at 0.236 to agreeing well once the controls
            were read from the design, and 59 of the 201 stored contrasts sit
            above ladder depth 0, i.e. are comparisons no study made directly.
            This figure earns figures 3 to 11. It goes early for that reason.

    3.  [~] RNA AND PROTEIN CO-EXPRESSION, CLUSTERED.
            Ends on: what the modules are.

    4.  [~] RNA RESPONSE TO TREATMENT, CLUSTERED.
            Ends on: which conditions group, and why.

    5.  [~] PROTEIN RESPONSE TO TREATMENT; THE BRADYZOITE DRUG REACH.
            The Reach module is the worked example and is built (`300cac69`).
            Ends on: subtype C is reached by 4 of 201 published contrasts
            against 86 for E, and it is not the matching giving up.

    6.  [~] WHICH MEASUREMENT EXPLAINS WHICH OTHER MEASUREMENT.
            Every panel the same shape: one measurement predicting another,
            cross-validated, mechanism in the caption. Panels: how much of
            protein abundance transcription explains, then translation
            efficiency, decay, modification; and the localisation panel, which
            is BUILT (`dc5baf51`).
            Ends on: wave-built organelles are placeable from transcription
            alone (5 of 13 compartments; rhoptries 0.975, micronemes 0.936) and
            the compartments made continuously are not (plasma membrane 0.67).

    7.  [ ] THE NETWORK.
            Einar's correction, 2026-11-18, and it changes what gets built: this
            is NOT co-expression drawn as a graph. Nodes carry TRAITS (a gene is
            transmembrane, its hyperLOPIT compartment) and edges carry DERIVED
            RELATIONSHIPS -- "x and y are both required for dissemination to the
            kidney but not the lung", "compound x moves x and not y".
            Co-expression is one edge type among many, not the graph.
            Ends on: what the combined graph shows that no single edge type does.
            Consequence for the program: the edge layer is first-class, not a
            view of the map. `star_edges`, `edge:cotranslation` and `bridge:host`
            are the start of it.

    8.  [ ] STRAINS. Ends on: what differs, and whether it is the effectors.

    9.  [ ] SPECIES, WITH MEASURED TRANSFER PER LAYER.
            Ends on: when a result in one species licenses a claim in another.
            Precedent: Toxoplasma orthologue locations predict Cryptosporidium
            location at 56.7 % against a 21.9 % baseline. Localisation clearly
            transfers; whether essentiality does is unmeasured.

    10. [ ] THE HOST. Response, vacuole proteins, AlphaFold co-folding, dual
            perturb-seq.
            Ends on: the shared host response and the tissue-specific part.
            Co-folding wording, decided by Einar 2026-11-18: **the strength is
            SCALE, and it cannot be controlled, and that is fine.** So the panel
            reports coverage -- how many parasite-host pairs folded, at what
            confidence -- and says plainly that it is uncontrolled. Where the
            Plasmodium crosslinks or MYR1 overlap a prediction, say what they
            show as an observation, never as a hit rate: a rate computed on that
            few pairs is a number with no precision pretending to be reassurance.

    11. [ ] FILLING THE MAP, SCORED ON HELD-OUT FACTS.
            Ends on: new host-facing, partner and function assignments at a
            validated rate. Without this the paper is a database description
            rather than the inference engine ParaScape is for.

    12. [ ] MEASURED ABSENCE.
            Genes and conditions no dataset touches, each with the smallest
            effect that would have been detected, so absence becomes a positive
            claim rather than a blank.
            Ends on: which datasets would be most valuable to the field -- which
            is how Einar framed the close, and it is prescriptive, not a caveat.

--------------------------------------------------------------------------------
THE TEST EVERY FIGURE MUST PASS
--------------------------------------------------------------------------------
    Einar, 2026-11-18: *"none of it has yet a biological question endpoint, and
    some of yours explicitly do not have a biological question at the end and
    rather have an experimental one."*

    So: **a figure ends in a sentence about the parasite.** "Which experiments
    are redundant" is an endpoint for a funding panel. "Transcription alone
    places the organelles built in waves during division, and not the ones built
    continuously" is a figure.

    And Einar's prior criterion, which is the filter applied FIRST:

        1.  INFORMATION GAP. It must not have been answered before -- including
            by us. Instruction 71 product 3 had already answered one of the ten.

    The structures by which combining datasets creates knowledge, as agreed:

        2a. META-ANALYSIS. Normalise two datasets until comparable, then map the
            information in one onto the other. The bradyzoite shape.
        2b. POOLING FOR POWER. Normalise several so they can be added, and gain
            the power for inference no single study could support.
        2c. COMPARING THE UNCOMPARED. Normalise so that conditions never set
            against each other can be. Einar: this is the heart of the current
            expression work and has not yet produced new knowledge. 59 of 201
            contrasts are already of this kind.
        2d. MEASURED ABSENCE. One dataset's absence is a detection limit; the
            union of many, each with its own limit, makes absence a claim.
        2e. DISCORDANCE AS THE SIGNAL. Where two methods of the same quantity
            disagree, the disagreement localises an artefact or a real subtlety.
            hyperLOPIT MAP and MCMC disagree for 26.0 % of 3,827 proteins, and
            that discordance has only ever been used defensively.
        2f. VARIANCE APPORTIONMENT. How much of an effect is biology and how
            much is lab, strain, medium, platform.
        2g. COMPOSITION. Chain datasets so the product is a rate or a lag no
            dataset contains.
        2h. PREDICTION AS A STATEMENT ABOUT THE DATA. Which evidence type
            predicts which other. Its BIOLOGICAL form is figure 6.
        2i. TRANSFER ACROSS A BOUNDARY. Fill a sparse system from a dense one
            and treat how much transfers as the result.

--------------------------------------------------------------------------------
WHAT THE TRIAGE OF INSTRUCTION 100 ESTABLISHED (2026-11-18)
--------------------------------------------------------------------------------
    Ten questions were proposed and Einar challenged them one at a time. The
    result is recorded because the FAILURE MODE is the reusable part: the ten
    were generated by pairing dataset NAMES and asking what the pairing would
    mean, without checking that both ends had data or that anybody had already
    answered it.

    WITHDRAWN (4):
      * effector x host-cell-type reachability -- the Toxoplasma crosslink set
        is parasite-only. The only parasite-to-host bridges in the release are
        `myr1_host_ip` (ONE bait) and `pf_host_bridge_xlms` (Plasmodium to
        erythrocyte). One bait is not a repertoire.
      * drug stage specificity from target engagement -- engagement data is
        Plasmodium, which has no cyst, so the motivation does not survive.
      * hyperLOPIT across stages -- tachyzoite only, 47.0 % of genes, MAP and
        MCMC disagreeing on 26.0 %. There is no second stage to compare with.
        It is a case for an experiment, not an analysis.
      * the stage switch being translational -- done by the source paper in
        Plasmodium (`pf_riboseq` is footprints and matched mRNA at five points);
        impossible in Toxoplasma (`gse302107_riboseq` is one stage).
      * host rewiring across tissues -- answered 2026-11-07 by instruction 71
        product 3, and the answer was that the data is too thin: one perturbation
        is stated by two independent studies.

    DOWNGRADED TO PLUMBING (2): the conditional essentialome (each paper already
    reports its own condition; the real job is putting seven screens on one
    scale, and `crispr_macrophage` has an INVERTED sign convention); epitopes
    against invariance (concrete and runnable, but the novelty is the
    study-effort-matched null rather than the question).

    GATED (1): guilt by complex, behind instruction 80 -- a bait that pulls
    everything agrees with everything, and the gate AUC fell 0.710 to 0.684 as
    the project count doubled.

--------------------------------------------------------------------------------
WHAT TO DO
--------------------------------------------------------------------------------
    A. THE FIGURES AND THE TEXT
       A1. [ ] Figure 2: the control-comparability analysis. Cluster studies on
               control agreement; add cell-cycle normalisation; measure which
               studies move. Einar named this; it gates the rest.
       A2. [ ] Import the GSE229505 host differential-expression table PER HOST
               GENE. What is loaded today is the per-effector summary: 20 PCA
               coordinates, a norm and a count, 252 screened with 22 full
               signatures. It answers which effectors change the host and by how
               much, and CANNOT answer which host gene an effector moves.
       A3. [ ] Figures 3-5, 8, 9: draw from the existing releases.
       A4. [x] Figure 6, localisation panel -- built, `dc5baf51`.
       A5. [ ] Figure 6, remaining panels: transcription, translation
               efficiency, decay and modification against protein abundance.
       A6. [ ] Figure 7: the trait-and-relationship network.
       A7. [ ] Figures 10-12.
       A8. [ ] The manuscript. Each figure in a two-cell table with the legend
               in the bottom row, as for the bradyzoite paper.

    B. EVERY ANALYSIS RUNNABLE INSIDE PARASCAPE
       B1. [x] The bradyzoite drug reach -- the Reach module, `300cac69`.
       B2. [ ] One module or view per figure whose analysis a user can re-run.
               A figure with no route through the program is the gap this names.
       B3. [ ] Each saves its result WITH what it was made from, as Reach does.

    C. RESHAPING THE PROGRAM (after the figures exist, not before)
       C1. [x] The edge layer is an object: `parascape/edge_layers.py`, thirteen
               kinds of relationship each carrying the worth figure 7 measured,
               and the layer built from the compartment labels refused by name
               for recovering a compartment. Six tests.

       C2. [!] REORGANISE THE BANDS TO MATCH THE FIGURE STRUCTURE -- recommended
               AGAINST, with the reason, rather than done.
               The bands are Look up, Explore, Ask, Bring your own and Data, and
               instruction 98 chose them deliberately: they are what the person
               came to DO. The figure structure is a different axis entirely --
               it is the order of an argument about the corpus. Rebuilding a
               navigation system around the shape of a document would make the
               program harder to use for everyone who is not reading the paper,
               which is everyone. The one place the paper DID change the
               catalogue is where it found something missing: the Findings
               module, in Data, because the analyses are corpus-level.
               Einar's call. Nothing has been moved.

       C3. [!] RETIRE OR MERGE MODULES NO FIGURE NEEDS -- measured, not done,
               because retiring a module removes something a user can open today
               and that is not a change to make unasked.
               Nine of 25 modules are not exercised by any figure: structure,
               cell, gallery, runs, fill, drop, screen, jobs, console.
               But look at WHAT they are. `drop` and `screen` are how a user
               brings their own data, which no figure does because the paper
               uses published data only. `jobs` and `console` are how anything
               long-running is watched. `runs` and `fill` are how an inference is
               kept. `gallery` and `cell` are how a map is compared and drawn.
               Not one of them is dead; they are the parts of the program the
               PAPER does not need, which is a different thing. The honest
               conclusion is that "no figure needs it" is the wrong test for
               whether a module should exist, and C3 as written should be
               dropped. The right test is whether a module is reachable, honest
               about what it cannot do, and used -- which is instruction 98 and
               is already held by the module contract.

HOW TO KNOW IT WORKED
    A figure is done when (a) it ends in a sentence about the parasite, (b) the
    analysis behind it runs inside ParaScape, and (c) a reader could repeat it
    from the saved record of what it was made from.

    The paper is done when every figure passes all three. "The script ran" is
    not a check.

TRAPS
    * ORGANISING BY DATA TYPE. A figure called "co-expression" or "perturbation"
      is a database chapter, and a reviewer reads twelve of them as twelve
      clusterings of one table.
    * A GRAB BAG FIGURE. "Four more questions like this" does not survive. Every
      panel of a figure is the same SHAPE of claim.
    * BUILDING THE MODULE BEFORE THE FIGURE. See section C.
    * THE HEATMAP REFLEX. Einar asked whether a heatmap suited the localisation
      panel; it did not, because colour cannot separate 0.94 from 0.93, and once
      the values are printed into the cells the colour is decoration on a table
      whose crossings are invisible. Pick the form from what the data has to
      show. Four forms were drawn before that panel was right.
    * PANEL LAYOUT. Labels are placed by MEASUREMENT, against candidate
      positions, rejecting any that overlaps a label, covers a point or leaves
      the page. Measure on `matplotlib.text.Text`, never on the annotation,
      which includes the leader and measures several times too large. Only
      Regular, Italic and Bold of Open Sans are installed on this machine.
