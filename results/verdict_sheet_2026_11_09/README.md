# Verdict sheet: the biological questions awaiting Einar's reading (instruction 59)

Built 2026-11-09 from `results/questions_2026_11_06/summary.csv`. **27 of 42 questions have no verdict.**

`grade` is the software's own statement about whether the method worked. **The verdict is the biology** -- is the answer true, useful, or an artefact -- and only you can give it. Write PASS, FAIL or INCONCLUSIVE against each qid and I will record it.

## Q51 · Tg · positive_unlabeled · graded *reliable*

**Which uncharacterised proteins behave like known rhoptry proteins?**

* Anchor: PMID:41975467
* What the run says: 123 positives against the rest of the proteome, 15 bags; the 300 genes that most resemble your list are listed, scored only by models that did not train on them. Nothing on the list is assumed about any other gene: unlabelled is treated as unknown, not as negative.
* It called: Q51_positive_unlabeled_candidates.csv: TGME49_285870 (SAG-related sequence SRS20A); TGME49_275670 (alveolin domain containing intermediate filament); TGME49_293430 (hypothetical protein); TGME49_233030 (gliding-associated protein GAP70); TGME49_312950 (RAP-1 family protein); TGME49_316260 (hypothetical protein)

**Verdict:** ______

## Q52 · Tg · positive_unlabeled · graded *reliable*

**Which proteins behave like the known PV-membrane nutrient channels?**

* Anchor: PMID:42555613
* What the run says: 193 positives against the rest of the proteome, 15 bags; the 300 genes that most resemble your list are listed, scored only by models that did not train on them. Nothing on the list is assumed about any other gene: unlabelled is treated as unknown, not as negative.
* It called: Q52_positive_unlabeled_candidates.csv: TGME49_220890 (dense granule protein GRA12A); TGME49_247195 (rhoptry discharge factor RDF1); TGME49_270700 (membrane protein MYR2); TGME49_311230 (basal complex component BCC7); TGME49_316710 (hypothetical protein); TGME49_309590 (rhoptry protein ROP1)

**Verdict:** ______

## Q53 · Tg · link_prediction · graded *reliable*

**Which genes are synthetic-lethal-like partners of GRA17 in the data?**

* Anchor: PMID:37498952
* What the run says: ValueError: the cofitness layer is the correlation structure of its source measurements, so completing it from its own visible edges recomputes a correlation (AUROC 0.99 however the non-pairs are drawn) rather than predicting a missed contact; strategy 34 predicts it from the other evidence instead
* It called: no per-gene table written

**Verdict:** ______

## Q54 · Tg · trait_regression · graded *reliable*

**Which effectors change host transcription, judged from the host-effect columns?**

* Anchor: PMID:37827122
* What the run says: ValueError: window shape cannot be larger than input array shape
* It called: no per-gene table written

**Verdict:** ______

## Q55 · Tg · conjunctions · graded *reliable*

**Which nuclear-targeted effector candidates are there beyond the known four?**

* Anchor: PMID:39292011
* What the run says: A map built without compartment or screen_any_phenotype; 0 clusters are enriched for a COMBINATION of the two beyond what either label alone explains (interaction ratio >= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted AND fitness-conferring in vivo.
* It called: no per-gene table written

**Verdict:** ______

## Q56 · Tg · condition_shift · graded *weak*

**Which in-vivo fitness hits are likely paracrine-rescued rather than cell-autonomous?**

* Anchor: PMID:39654402
* What the run says: The part of fit_invivo_PE that fit_invitro_hff does not explain -- the condition-specific effect -- for 7,324 genes. Other measurements predict it at out-of-fold rank correlation 0.076: little of it is predictable from the rest of the table. Negative shift: needed more in the condition than the baseline predicts.
* It called: Q56_condition_shift_predicted.csv: TGME49_277240 (NTPase I); TGME49_310110 (hypothetical protein); TGME49_500287 (phosphatidylinositol 4-kinase, putative); TGME49_215055 (hypothetical protein); TGME49_250930 (hypothetical protein); TGME49_500055 (arginine and glutamate-rich protein, putative)

**Verdict:** ______

## Q57 · Tg · supervised_classifier · graded *reliable*

**Which intrinsically disordered proteins look like exported effectors?**

* Anchor: PMID:38747635
* What the run says: A class-balanced logistic regression on 364 permitted measurements, trained on 3,804 labelled genes, calls 2,000 unlabelled genes (probability >= 0.00). 'what drives each class' lists the measurements with the largest weights -- the model's reasons, to be read as hypotheses about what defines each class.
* It called: Q57_supervised_classifier_calls.csv: TGME49_224760 (SAG-related sequence SRS40E); TGME49_284010 (5'-3' exonuclease, N-terminal resolvase family d); TGME49_251740 (AP2 domain transcription factor AP2XII-9); TGME49_278390 (Toxoplasma gondii family A protein); TGME49_320250 (SAG-related sequence SRS15A); TGME49_252200 (palmitoyltransferase DHHC7)

**Verdict:** ______

## Q58 · Tg · paralog_divergence · graded *reliable*

**Which ROP-family proteins look like tethers rather than enzymes?**

* Anchor: PMID:41291278
* What the run says: 3,452 paralog pairs scored by how differently they behave across the permitted measurements (one minus the correlation of their profiles). The top of the list is where gene duplication was followed by a change of job -- sub- or neo-functionalisation -- and the bottom is redundancy.
* It called: no per-gene table written

**Verdict:** ______

## Q61 · Tg · cluster_guilt · graded *reliable*

**Which conoid and apical-complex proteins are still unassigned to a position?**

* Anchor: PMID:41366525
* What the run says: 10 clusters are enriched for one compartment label at q < 0.05 and lift >= 2 (hypergeometric, BH-corrected over every cluster x label); their 321 unlabelled members are called with it.
* It called: Q61_cluster_guilt_calls.csv: TGME49_287960 (hypothetical protein); TGME49_250720 (hypothetical protein); TGME49_291830 (hypothetical protein, conserved); TGME49_318500 (cpw-wpc domain-containing protein); TGME49_223030 (tetratricopeptide repeat-containing protein); TGME49_262890 (hypothetical protein)

**Verdict:** ______

## Q62 · Tg · multiplex_modules · graded *weak*

**Which micronemal proteins partition into a secretion module and a junction module?**

* Anchor: PMID:42470170
* What the run says: 6 communities that the layers coexpression, cofitness, cotranslation, domain, ip_ms, struct, xlms agree on; 7,822 genes placed in one. 2,000 unlabelled members are called by their community's majority compartment.
* It called: Q62_multiplex_modules_calls.csv: TGME49_313630 (VPS13 domain-containing protein); TGME49_209000 (HECT-domain (ubiquitin-transferase) domain-conta); TGME49_225745 (hypothetical protein); TGME49_294820 (type I fatty acid synthase, putative); TGME49_203665 (hypothetical protein); TGME49_242625 (ATPase family associated with various cellular a)

**Verdict:** ______

## Q64 · Tg · condition_shift · graded *weak*

**Which apicoplast proteins become dispensable when the organelle is bypassed?**

* Anchor: PMID:40025025
* What the run says: The part of fit_complete_medium_2025 that fit_invitro_hff does not explain -- the condition-specific effect -- for 7,325 genes. Other measurements predict it at out-of-fold rank correlation 0.082: little of it is predictable from the rest of the table. Negative shift: needed more in the condition than the baseline predicts.
* It called: Q64_condition_shift_predicted.csv: TGME49_327600 (exoribonuclease, putative); TGME49_325000 (hypothetical protein); TGME49_325300 (proteasome subunit alpha type, putative); TGME49_326000 (RNA binding protein, putative); TGME49_327000 (U5 small nuclear ribonucleoprotein component (U5); TGME49_302009 (ORF B)

**Verdict:** ______

## Q65 · Tg · understudied_first · graded *reliable*

**Which predicted transporters have no publication at all?**

* Anchor: PMID:41582196
* What the run says: 7,390 genes have no focal or substantive publication; 1,723 of those without a compartment label receive an agreed call, ranked by agreement times novelty. These are the cheapest discoveries: a hypothesis the data supports about a gene nobody has written about.
* It called: Q65_understudied_first_candidates.csv: TGME49_209000 (HECT-domain (ubiquitin-transferase) domain-conta); TGME49_225745 (hypothetical protein); TGME49_266010 (phosphatidylinositol 3- and 4-kinase); TGME49_286270 (hypothetical protein); TGME49_271145 (hypothetical protein); TGME49_253870 (hypothetical protein)

**Verdict:** ______

## Q66 · Tg · conjunctions · graded *reliable*

**Which genes are induced in bradyzoites and localize to the cyst wall?**

* Anchor: PMID:37081202
* What the run says: A map built without compartment or stage_enriched_derived; 2 clusters are enriched for a COMBINATION of the two beyond what either label alone explains (interaction ratio >= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted AND fitness-conferring in vivo.
* It called: no per-gene table written

**Verdict:** ______

## Q67 · Tg · condition_shift · graded *weak*

**Which differentiation factors are translationally rather than transcriptionally controlled?**

* Anchor: PMID:38782906
* What the run says: The part of te245775_parent_prebrady_r1 that rna245775_parent_prebrady_r1 does not explain -- the condition-specific effect -- for 7,880 genes. Other measurements predict it at out-of-fold rank correlation 0.596: the condition-specific need has a signature, and 'predicted' extends it to genes never screened in the condition. Negative shift: needed more in the condition than the baseline predicts.
* It called: Q67_condition_shift_predicted.csv: TGME49_500314 (protein phosphatase PPM11C); TGME49_500435 (HECT-domain (ubiquitin-transferase) -containing ); TGME49_500337 (ubiquitin-like protein, putative); TGME49_500146 (guanine-nucleotide exchange factor, putative); TGME49_500048 (hypothetical protein, conserved); TGME49_500081 (hypothetical protein, conserved)

**Verdict:** ______

## Q70 · Tg · condition_shift · graded *weak*

**Which uncharacterised dense-granule proteins have a stage-specific phenotype?**

* Anchor: PMID:40941387
* What the run says: The part of fit_invivo_brain that fit_invitro_hff does not explain -- the condition-specific effect -- for 7,324 genes. Other measurements predict it at out-of-fold rank correlation 0.183: the condition-specific need has a signature, and 'predicted' extends it to genes never screened in the condition. Negative shift: needed more in the condition than the baseline predicts.
* It called: Q70_condition_shift_predicted.csv: TGME49_302053 (ribosomal protein S7); TGME49_303450 (coatomer gamma subunit); TGME49_328000 (hypothetical protein); TGME49_322900 (hypothetical protein); TGME49_329620 (hypothetical protein); TGME49_302057 (ribosomal protein S7)

**Verdict:** ______

## Q73 · Tg · conjunctions · graded *reliable*

**Which AP2 factors are cell-cycle rather than stage regulators?**

* Anchor: PMID:40920096
* What the run says: A map built without cellcycle_phase or stage_enriched_derived; 0 clusters are enriched for a COMBINATION of the two beyond what either label alone explains (interaction ratio >= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted AND fitness-conferring in vivo.
* It called: no per-gene table written

**Verdict:** ______

## Q75 · Tg · supervised_classifier · graded *reliable*

**Which palmitoylated proteins are IMC or plasma-membrane residents?**

* Anchor: PMID:39690890
* What the run says: A class-balanced logistic regression on 364 permitted measurements, trained on 3,804 labelled genes, calls 2,000 unlabelled genes (probability >= 0.00). 'what drives each class' lists the measurements with the largest weights -- the model's reasons, to be read as hypotheses about what defines each class.
* It called: Q75_supervised_classifier_calls.csv: TGME49_224760 (SAG-related sequence SRS40E); TGME49_284010 (5'-3' exonuclease, N-terminal resolvase family d); TGME49_251740 (AP2 domain transcription factor AP2XII-9); TGME49_278390 (Toxoplasma gondii family A protein); TGME49_320250 (SAG-related sequence SRS15A); TGME49_252200 (palmitoyltransferase DHHC7)

**Verdict:** ______

## Q80 · Tg · conjunctions · graded *reliable*

**Which hypothetical genes are stage-specifically induced?**

* Anchor: PMID:41761089
* What the run says: A map built without stage_enriched_derived or compartment; 2 clusters are enriched for a COMBINATION of the two beyond what either label alone explains (interaction ratio >= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted AND fitness-conferring in vivo.
* It called: no per-gene table written

**Verdict:** ______

## Q82 · Tg · paralog_divergence · graded *reliable*

**Which alveolin-family paralogues divided the subpellicular network between them?**

* Anchor: PMID:42113865
* What the run says: 3,452 paralog pairs scored by how differently they behave across the permitted measurements (one minus the correlation of their profiles). The top of the list is where gene duplication was followed by a change of job -- sub- or neo-functionalisation -- and the bottom is redundancy.
* It called: no per-gene table written

**Verdict:** ______

## Q83 · Tg · conjunctions · graded *reliable*

**Which apical-ring proteins are assembly factors rather than structural?**

* Anchor: PMID:41315181
* What the run says: A map built without compartment or cellcycle_phase; 2 clusters are enriched for a COMBINATION of the two beyond what either label alone explains (interaction ratio >= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted AND fitness-conferring in vivo.
* It called: no per-gene table written

**Verdict:** ______

## Q84 · Tg · cluster_guilt · graded *reliable*

**Which microtubule inner proteins are still hypotheticals?**

* Anchor: PMID:41626803
* What the run says: 10 clusters are enriched for one compartment label at q < 0.05 and lift >= 2 (hypergeometric, BH-corrected over every cluster x label); their 321 unlabelled members are called with it.
* It called: Q84_cluster_guilt_calls.csv: TGME49_287960 (hypothetical protein); TGME49_250720 (hypothetical protein); TGME49_291830 (hypothetical protein, conserved); TGME49_318500 (cpw-wpc domain-containing protein); TGME49_223030 (tetratricopeptide repeat-containing protein); TGME49_262890 (hypothetical protein)

**Verdict:** ______

## Q86 · Tg · condition_shift · graded *weak*

**Which genes need glucose specifically, beyond general fitness?**

* Anchor: PMID:40025025
* What the run says: The part of fit_no_glucose that fit_complete_medium_2025 does not explain -- the condition-specific effect -- for 7,393 genes. Other measurements predict it at out-of-fold rank correlation 0.325: the condition-specific need has a signature, and 'predicted' extends it to genes never screened in the condition. Negative shift: needed more in the condition than the baseline predicts.
* It called: Q86_condition_shift_predicted.csv: TGME49_249180 (bifunctional dihydrofolate reductase-thymidylate); TGME49_300690 (RNA polymerase C2, putative); TGME49_223150 (START domain-containing protein); TGME49_232940 (heat shock protein HSP20); TGME49_500174 (KRUF family protein); TGME49_500183 (glyoxalase, putative)

**Verdict:** ______

## Q87 · Tg · supervised_classifier · graded *reliable*

**Which iron-regulated proteins are apicoplast or mitochondrial?**

* Anchor: PMID:39348385
* What the run says: A class-balanced logistic regression on 364 permitted measurements, trained on 3,804 labelled genes, calls 2,000 unlabelled genes (probability >= 0.00). 'what drives each class' lists the measurements with the largest weights -- the model's reasons, to be read as hypotheses about what defines each class.
* It called: Q87_supervised_classifier_calls.csv: TGME49_224760 (SAG-related sequence SRS40E); TGME49_284010 (5'-3' exonuclease, N-terminal resolvase family d); TGME49_251740 (AP2 domain transcription factor AP2XII-9); TGME49_278390 (Toxoplasma gondii family A protein); TGME49_320250 (SAG-related sequence SRS15A); TGME49_252200 (palmitoyltransferase DHHC7)

**Verdict:** ______

## Q89 · Pf · understudied_first · graded *works when tuned*

**Which ApiAP2 factors still have no target set in the data?**

* Anchor: PMID:39419713
* What the run says: 4,540 genes have no focal or substantive publication; 2,000 of those without a lopit_pf_location label receive an agreed call, ranked by agreement times novelty. These are the cheapest discoveries: a hypothesis the data supports about a gene nobody has written about.
* It called: Q89_understudied_first_candidates.csv: PF3D7_0100100 (erythrocyte membrane protein 1, PfEMP1); PF3D7_0100400 (rifin); PF3D7_0100700 (Plasmodium exported protein, unknown function, f); PF3D7_0101100 (exported protein family 4); PF3D7_0101600 (rifin); PF3D7_0101800 (stevor)

**Verdict:** ______

## Q93 · Pf · understudied_first · graded *works when tuned*

**Which essential genes are structurally druggable but unpublished?**

* Anchor: PMID:40066064
* What the run says: 4,540 genes have no focal or substantive publication; 2,000 of those without a lopit_pf_location label receive an agreed call, ranked by agreement times novelty. These are the cheapest discoveries: a hypothesis the data supports about a gene nobody has written about.
* It called: Q93_understudied_first_candidates.csv: PF3D7_0100100 (erythrocyte membrane protein 1, PfEMP1); PF3D7_0100400 (rifin); PF3D7_0100700 (Plasmodium exported protein, unknown function, f); PF3D7_0101100 (exported protein family 4); PF3D7_0101600 (rifin); PF3D7_0101800 (stevor)

**Verdict:** ______

## Q96 · Pf · conjunctions · graded *untestable*

**Which febrile-temperature responders are also phospho-regulated?**

* Anchor: PMID:41960867
* What the run says: A map built without stage_enriched_derived or has_phospho; 4 clusters are enriched for a COMBINATION of the two beyond what either label alone explains (interaction ratio >= 1.25). Each is a kind of gene defined by both: nuclear AND S-phase, secreted AND fitness-conferring in vivo.
* It called: no per-gene table written

**Verdict:** ______

## Q100 · Pf · condition_shift · graded *reliable*

**Which surface proteins are reversibly regulated with population state?**

* Anchor: PMID:40354414
* What the run says: The part of expr_gametocyte_v that expr_gametocyte_ii does not explain -- the condition-specific effect -- for 5,720 genes. Other measurements predict it at out-of-fold rank correlation 0.486: the condition-specific need has a signature, and 'predicted' extends it to genes never screened in the condition. Negative shift: needed more in the condition than the baseline predicts.
* It called: Q100_condition_shift_shifted_genes.csv: PF3D7_0707000 (tRNA Lysine); PF3D7_0706800 (tRNA Threonine); PF3D7_1103200 (tRNA Glycine); PF3D7_1438300 (tRNA Methionine); PF3D7_0312600 (tRNA Valine); PF3D7_1438200 (tRNA Selenocysteine)

**Verdict:** ______

