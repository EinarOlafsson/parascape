# Function class rules (2026-11-04), verbatim

Applied by `scripts/build_function_truth_notebook.py`. Two rule kinds per class.

* **GO rule**: a gene matches the class if one of its *experimental* GO terms is the seed term or a descendant of it over `is_a`/`part_of` in the `go-basic.obo` snapshot of 2026-11-04. A class reached this way is tier `observed-direct`, probability 0.95.
* **Text rule**: a Python regular expression (`re.I`) searched in `"<recommended protein name>. <FUNCTION comment>"` of a reviewed UniProt entry. A class reached only this way is tier `observed-curated`, probability 0.90.
* **Screen rule**: a True phenotype column of the arrayed imaging screen supports one process class and nothing else. Tier `observed-screen`, probability 0.80.

Never used to assign a class: InterPro/Pfam domains, ToxoDB product names, IEA/ISS/ISO GO annotations, UniProt keywords and domain features, and CRISPR fitness scores.

## `host_cell_modulation`

GO seeds:

* `GO:0035821` modulation of process of another organism (biological_process, 291 descendants, structural IC 6.36)
* `GO:0044003` symbiont-mediated perturbation of host process (biological_process, 190 descendants, structural IC 6.97)
* `GO:0052200` symbiont-mediated response to host defenses (biological_process, 31 descendants, structural IC 9.55)
* `GO:0035325` Toll-like receptor binding (molecular_function, 2 descendants, structural IC 11.71)

Text patterns:

* `suppress\w*[^.]{0,60}host`
* `modulat\w*[^.]{0,60}host`
* `host (innate )?immun`
* `host (interferon|cytokine|transcription|gene expression)`
* `immunity-related GTPase`
* `\bIRG[AB]?\d*\b`
* `\bIFN-gamma\b`
* `export\w*[^.]{0,40}(into|to) the host`
* `tether\w*[^.]{0,60}host`
* `host (mitochondri|nucleus|nuclear|cytosol|endocytic|lysosom)`
* `\bvirulence (factor|in mice)\b`
* `NF-kappa-?B`
* `\bSTAT[136]\b`
* `effector[^.]{0,40}(export|translocat)`
* `host cell cycle`

Screen columns: none

## `invasion_motility`

GO seeds:

* `GO:0044409` symbiont entry into host (biological_process, 43 descendants, structural IC 9.09)
* `GO:0044406` adhesion of symbiont to host (biological_process, 28 descendants, structural IC 9.69)
* `GO:0046812` host cell surface binding (molecular_function, 1 descendants, structural IC 12.29)
* `GO:2000145` regulation of cell motility (biological_process, 231 descendants, structural IC 6.69)
* `GO:0060352` cell adhesion molecule production (biological_process, 0 descendants, structural IC 14.55)
* `GO:0098609` cell-cell adhesion (biological_process, 55 descendants, structural IC 8.74)

Text patterns:

* `invasion`
* `invad\w+`
* `adhesin`
* `attachment (of|to)`
* `moving junction`
* `motility`
* `glid\w+`
* `conoid`
* `apical complex`
* `binding to host cells?`
* `microneme secretion`
* `rhoptry secretion`
* `\bMJ complex\b`

Screen columns: none

## `egress`

GO seeds:

* `GO:0035891` exit from host cell (biological_process, 18 descendants, structural IC 10.30)

Text patterns:

* `egress`
* `exit from (the )?host`
* `lysis of the host cell`
* `permeabiliz\w+[^.]{0,60}vacuol`
* `parasite exit`

Screen columns: `screen_egress_phenotype`

## `cell_cycle_division`

GO seeds:

* `GO:0007049` cell cycle (biological_process, 413 descendants, structural IC 5.86)
* `GO:0051301` cell division (biological_process, 120 descendants, structural IC 7.63)
* `GO:0007059` chromosome segregation (biological_process, 86 descendants, structural IC 8.11)
* `GO:0042127` regulation of cell population proliferation (biological_process, 258 descendants, structural IC 6.53)
* `GO:0006260` DNA replication (biological_process, 81 descendants, structural IC 8.19)

Text patterns:

* `cell (division|cycle)`
* `daughter (cell|parasite)`
* `centrosome`
* `cytokinesis`
* `\bbudding\b`
* `mitosis`
* `mitotic`
* `G[12]/[MS]`
* `intracellular (multiplication|replication)`
* `parasite (division|replication|multiplication)`
* `chromosome`
* `\bcell polarity\b`

Screen columns: `screen_replication_phenotype`

## `protein_trafficking_secretion`

GO seeds:

* `GO:0008104` intracellular protein localization (biological_process, 348 descendants, structural IC 6.10)
* `GO:0045055` regulated exocytosis (biological_process, 36 descendants, structural IC 9.34)
* `GO:0016192` vesicle-mediated transport (biological_process, 210 descendants, structural IC 6.83)
* `GO:1903305` regulation of regulated secretory pathway (biological_process, 44 descendants, structural IC 9.06)
* `GO:0099500` vesicle fusion to plasma membrane (biological_process, 2 descendants, structural IC 12.96)
* `GO:0009306` protein secretion (biological_process, 22 descendants, structural IC 10.03)
* `GO:0007030` Golgi organization (biological_process, 15 descendants, structural IC 10.55)

Text patterns:

* `traffick\w+`
* `secretion`
* `secretory`
* `exocytosis`
* `export of`
* `import of[^.]{0,60}into the`
* `vesicul\w+`
* `localization to[^.]{0,40}membrane`
* `escorter`
* `sort\w+ (of|into)`
* `\bGolgi\b`
* `endoplasmic reticulum`
* `dense granule (transport|exocytosis)`

Screen columns: none

## `signalling`

GO seeds:

* `GO:0004672` protein kinase activity (molecular_function, 115 descendants, structural IC 6.44)
* `GO:0006468` protein phosphorylation (biological_process, 19 descendants, structural IC 10.23)
* `GO:0007165` signal transduction (biological_process, 569 descendants, structural IC 5.39)
* `GO:0009190` cyclic nucleotide biosynthetic process (biological_process, 2 descendants, structural IC 12.96)
* `GO:0030551` cyclic nucleotide binding (molecular_function, 7 descendants, structural IC 10.29)
* `GO:0004721` phosphoprotein phosphatase activity (molecular_function, 36 descendants, structural IC 8.08)
* `GO:0070887` cellular response to chemical stimulus (biological_process, 669 descendants, structural IC 5.16)

Text patterns:

* `\bkinase\b`
* `\bcyclase\b`
* `phosphorylat\w+`
* `phosphatase`
* `second messenger`
* `\bcGMP\b`
* `\bcAMP\b`
* `signal transduc\w+`
* `calcium[- ]dependent`
* `\bpseudokinase\b`
* `senses? `
* `detects and transduces`

Screen columns: none

## `transport_homeostasis`

GO seeds:

* `GO:0022857` transmembrane transporter activity (molecular_function, 990 descendants, structural IC 3.34)
* `GO:0055085` transmembrane transport (biological_process, 476 descendants, structural IC 5.65)
* `GO:0015267` channel activity (molecular_function, 171 descendants, structural IC 5.87)
* `GO:0050801` monoatomic ion homeostasis (biological_process, 73 descendants, structural IC 8.34)
* `GO:0019725` cellular homeostasis (biological_process, 103 descendants, structural IC 7.85)
* `GO:0140115` export across plasma membrane (biological_process, 18 descendants, structural IC 10.30)

Text patterns:

* `transporter`
* `\bchannel\b`
* `permeable to`
* `homeostasis`
* `import across`
* `uptake of`
* `\befflux\b`
* `extrusion of`
* `\btransports?\b[^.]{0,40}?\b(across|into|out of)\b`
* `\btransports?\b[^.]{0,30}?\b(ions?|water|glycerol|urea|pyruvate|glucose|hexose|amino acids?|lysine|arginine|methionine|nucleotides?|sugars?|lipids?|cholesterol)\b`

Screen columns: none

## `metabolism`

GO seeds:

* `GO:0044281` small molecule metabolic process (biological_process, 1476 descendants, structural IC 4.02)
* `GO:0005975` carbohydrate metabolic process (biological_process, 315 descendants, structural IC 6.25)
* `GO:0055086` nucleobase-containing small molecule metabolic process (biological_process, 555 descendants, structural IC 5.43)
* `GO:0006520` amino acid metabolic process (biological_process, 174 descendants, structural IC 7.10)
* `GO:0006091` generation of precursor metabolites and energy (biological_process, 96 descendants, structural IC 7.95)
* `GO:0140640` catalytic activity, acting on a nucleic acid (molecular_function, 444 descendants, structural IC 4.50)

Text patterns:

* `catalyz\w+`
* `biosynthe\w+`
* `metaboli\w+`
* `glycolysis`
* `glycolytic`
* `salvage`
* `\bnucleotide`
* `\bpurine\b`
* `\bpyrimidine\b`
* `tricarboxylic acid`
* `pentose phosphate`
* `\bfolate\b`
* `\bdehydrogenase\b`
* `\bsynthase\b`
* `\btransferase\b(?! )`
* `amino acid (biosynthesis|metabolism)`

Screen columns: none

## `lipid_membrane_biology`

GO seeds:

* `GO:0006629` lipid metabolic process (biological_process, 578 descendants, structural IC 5.37)
* `GO:0008289` lipid binding (molecular_function, 94 descendants, structural IC 6.72)
* `GO:0061024` membrane organization (biological_process, 234 descendants, structural IC 6.67)
* `GO:0006869` lipid transport (biological_process, 107 descendants, structural IC 7.79)
* `GO:0008654` phospholipid biosynthetic process (biological_process, 43 descendants, structural IC 9.09)

Text patterns:

* `\blipid`
* `phospholipid`
* `phosphatidyl`
* `fatty acid`
* `\bsterol\b`
* `cholesterol`
* `membrane (fusion|organization|remodel)`
* `binds lipid`
* `\bflippase\b`
* `\bpore[- ]forming\b`
* `(membranous|intravacuolar) network`

Screen columns: none

## `apicoplast_mitochondrion`

GO seeds:

* none -- this class is text-derived only (the matching GO terms are cellular-component terms, which are a different trait of the map)

Text patterns:

* `apicoplast`
* `(?<!host )mitochondri\w+`
* `\bplastid\b`
* `isoprenoid`
* `methylerythritol`
* `\bMEP pathway\b`
* `\bFASII\b`
* `respiratory chain`

Screen columns: `screen_apicoplast_phenotype`

## `translation`

GO seeds:

* `GO:0006412` translation (biological_process, 84 descendants, structural IC 8.14)
* `GO:0042254` ribosome biogenesis (biological_process, 76 descendants, structural IC 8.28)
* `GO:0006399` tRNA metabolic process (biological_process, 121 descendants, structural IC 7.62)
* `GO:0003735` structural constituent of ribosome (molecular_function, 0 descendants, structural IC 13.29)
* `GO:0045182` translation regulator activity (molecular_function, 5 descendants, structural IC 10.71)

Text patterns:

* `\btranslation (of|initiation|elongation|termination)\b`
* `translational (control|regulation|repression)`
* `ribosom\w+`
* `\btRNA\b`
* `aminoacyl`
* `\belongation factor\b`
* `\binitiation factor\b`

Screen columns: none

## `transcription_chromatin`

GO seeds:

* `GO:0003700` DNA-binding transcription factor activity (molecular_function, 16 descendants, structural IC 9.21)
* `GO:0006355` regulation of DNA-templated transcription (biological_process, 100 descendants, structural IC 7.89)
* `GO:0006325` chromatin organization (biological_process, 63 descendants, structural IC 8.55)
* `GO:0006397` mRNA processing (biological_process, 35 descendants, structural IC 9.38)
* `GO:0010468` regulation of gene expression (biological_process, 822 descendants, structural IC 4.86)
* `GO:0003677` DNA binding (molecular_function, 130 descendants, structural IC 6.26)
* `GO:0003723` RNA binding (molecular_function, 163 descendants, structural IC 5.94)

Text patterns:

* `(?<!host )transcription`
* `chromatin`
* `histone`
* `promoter`
* `\bmRNA\b`
* `\bsplic\w+`
* `(?<!host )gene expression`
* `DNA[- ]binding`
* `\bRNA polymerase\b`
* `\bcapping\b`

Screen columns: none

## `stress_differentiation`

GO seeds:

* `GO:0006950` response to stress (biological_process, 584 descendants, structural IC 5.36)
* `GO:0048869` cellular developmental process (biological_process, 1366 descendants, structural IC 4.13)
* `GO:0022611` dormancy process (biological_process, 26 descendants, structural IC 9.79)
* `GO:0045165` cell fate commitment (biological_process, 245 descendants, structural IC 6.61)
* `GO:0045595` regulation of cell differentiation (biological_process, 697 descendants, structural IC 5.10)

Text patterns:

* `bradyzoite`
* `differentiat\w+`
* `\bcyst\b`
* `cyst wall`
* `dormanc\w+`
* `\bstress\b`
* `heat shock`
* `oxidative`
* `starvation`
* `chronic[- ]stage`
* `\bpersistence\b`
* `reactivation`

Screen columns: none

## `proteostasis`

GO seeds:

* `GO:0006508` proteolysis (biological_process, 45 descendants, structural IC 9.03)
* `GO:0006457` protein folding (biological_process, 5 descendants, structural IC 11.96)
* `GO:0006914` autophagy (biological_process, 45 descendants, structural IC 9.03)
* `GO:0008233` peptidase activity (molecular_function, 81 descendants, structural IC 6.94)
* `GO:0051604` protein maturation (biological_process, 72 descendants, structural IC 8.36)
* `GO:0070647` protein modification by small protein conjugation or removal (biological_process, 38 descendants, structural IC 9.26)

Text patterns:

* `\bprotease\b`
* `peptidase`
* `proteoly\w+`
* `(?<!host )autophag\w+`
* `chaperon\w+`
* `\bfolding\b`
* `ubiquitin`
* `proteasom\w+`
* `degradation of`
* `processing[- ]mediated maturation`
* `\bcleav\w+\b`
* `\bisomerase\b`

Screen columns: none

## `cytoskeleton`

GO seeds:

* `GO:0007010` cytoskeleton organization (biological_process, 216 descendants, structural IC 6.79)
* `GO:0030029` actin filament-based process (biological_process, 100 descendants, structural IC 7.89)
* `GO:0003779` actin binding (molecular_function, 3 descendants, structural IC 11.29)
* `GO:0007017` microtubule-based process (biological_process, 204 descendants, structural IC 6.87)
* `GO:0005200` structural constituent of cytoskeleton (molecular_function, 3 descendants, structural IC 11.29)
* `GO:0003774` cytoskeletal motor activity (molecular_function, 10 descendants, structural IC 9.83)
* `GO:0003789` actin filament severing activity (molecular_function, 0 descendants, structural IC 13.29)
* `GO:0003785` actin monomer binding (molecular_function, 0 descendants, structural IC 13.29)

Text patterns:

* `\bactin\b`
* `tubulin`
* `microtubule`
* `\bmyosin\b`
* `cytoskelet\w+`
* `inner membrane complex`
* `\bIMC\b`
* `\bspectrin\b`
* `\bseptin\b`
* `peripheral annuli`
* `\bF-actin\b`

Screen columns: `screen_actin_phenotype`

