# Domain traits (InterPro and Pfam), tier 2 of the information map

Built by scripts/build_domain_traits_notebook.py (notebook notebooks/domain_traits_2026_11_04.ipynb). One row per gene, the 8,140 genes of parascape/data/nodes.parquet.

## Files
* domain_traits.parquet: 313 columns (index gene_id): 299 family flags dom_<accession>, 10 named classes domclass_*, 4 counts/flags.
* domain_lookup.tsv: column, accession, description, gene count, source for every family column.
* coverage.tsv: how many genes have what.

## Frequency cutoff
A family must occur in at least 10 genes. Columns at each cutoff:
* 5: 514 InterPro / 226 Pfam / 740 total
* 10: 221 InterPro / 78 Pfam / 299 total
* 20: 82 InterPro / 21 Pfam / 103 total
* 50: 26 InterPro / 6 Pfam / 32 total
Chosen: 10, giving 221 InterPro + 78 Pfam = 299 family columns.

## Matching rules for the named classes (verbatim)
* **protein_kinase**: InterPro ['IPR000719', 'IPR001245', 'IPR020635', 'IPR027916']; Pfam ['PF00069', 'PF07714', 'PF14531']
* **ap2_apiap2**: InterPro ['IPR001471']; Pfam ['PF00847']
* **srs_sag**: InterPro ['IPR007226', 'IPR036755', 'IPR028352']; Pfam ['PF04092']
* **ef_hand**: InterPro ['IPR002048', 'IPR011992']; Pfam ['PF00036', 'PF13499', 'PF13202', 'PF13405', 'PF13833']
* **wd40**: InterPro ['IPR001680', 'IPR036322', 'IPR019775', 'IPR020472']; Pfam ['PF00400']
* **ankyrin**: InterPro ['IPR002110', 'IPR036770']; Pfam ['PF00023', 'PF12796', 'PF13606', 'PF13637', 'PF13857']
* **rna_binding**: InterPro ['IPR000504', 'IPR003954', 'IPR001313', 'IPR033133', 'IPR004087', 'IPR004088', 'IPR036612']; Pfam ['PF00076', 'PF00806', 'PF00013']
* **proteasome_ubiquitin**: 136 InterPro and 59 Pfam entries whose description matches the regex below.
* proteasome_ubiquitin regex (case-insensitive, on the InterPro and Pfam descriptions): proteasome|ubiquitin|HECT|F-box|cullin|SUMO|NEDD8|neddylation
* rhoptry_kinase_ropk: IPR027916 or PF14531 (Apicomplexa protein kinase-like), OR protein_kinase AND product name matches regex rhoptry|\bROP\d+ (case-insensitive). The second arm uses the product annotation, not a domain.
* tm_rich: DeepTMHMM n_tm >= 3.
* Broad superfamily entries (for example IPR011009 Protein kinase-like superfamily) are deliberately excluded.

Gene counts per class:
* domclass_protein_kinase: 173 genes
* domclass_ap2_apiap2: 45 genes
* domclass_srs_sag: 110 genes
* domclass_ef_hand: 82 genes
* domclass_wd40: 129 genes
* domclass_ankyrin: 30 genes
* domclass_rna_binding: 104 genes
* domclass_proteasome_ubiquitin: 199 genes
* domclass_rhoptry_kinase_ropk: 50 genes
* domclass_tm_rich: 668 genes

Known-gene spot checks (CDPK1, ROP18, AP2XI-5, SAG1), missed: none

## Coverage
4967 of 8140 genes have at least one domain; 3173 (39.0%) have none (domain_novel). Check against the raw file: 3 InterPro and 3 Pfam lists differ from the raw file; 1 gene(s) at 20 IDs.
Source: ToxoDB gene attributes (datasets/interpro_tgon.csv, registry key interpro; the registry holds no PMID or accession).

## Host-facing cross-check
domain_novel in 46 of 65 host-facing genes (70.8%) versus 3173 of 8140 genes overall (39.0%); Fisher exact, host-facing vs all other genes: odds ratio 3.83, p = 2.55e-07. Only 65 labelled genes: a hint, not a finding.

## Caveats
* An InterPro or Pfam assignment is itself a prediction from the amino-acid sequence, not an experiment. False means not recognised, not proven absent.
* Domain annotation partly records study effort, not only conserved architecture (registry note).
* Missing is not zero: domain_novel genes may carry real domains that no database knows.
