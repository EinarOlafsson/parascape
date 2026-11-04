# Host-facing Toxoplasma gondii proteins — candidate truth table (2026-11-04)

Candidate rows for Einar to review BEFORE any of this is used to score predictions.
Nothing here is a finished label: every row carries an empty `einar_decision` column.

Files:

- `candidates.tsv` — 140 rows, 98 distinct ME49 gene IDs (plus 3 rows whose ID is unresolved).
- this README.

## Definition used

Host-facing = the protein, or a cleaved part of it, is found beyond the parasite plasma
membrane facing the host: the PV membrane (PVM), the PV lumen **only** when the protein is
exported further or engages host factors, the host cytosol, host nucleus, host plasma
membrane, a host organelle (e.g. host mitochondria), or injected rhoptry/RON proteins in the
host or at the moving junction. Residing in a dense granule or a rhoptry is **not** enough.

`host_compartment` is one of: `PVM`, `host_cytosol`, `host_nucleus`, `host_plasma_membrane`,
`host_organelle`, `moving_junction`, `PV_lumen_host_interacting`,
`organelle_only_not_exported`.

## Rules followed (all of them, verbatim in effect)

1. **Every row needs a PMID verified to exist and to match the claim.** All 56 distinct PMIDs
   were checked with NCBI E-utilities `esummary.fcgi` (db=pubmed, retmode=json); title,
   journal, year and first author were read back and the `year`/`first_author` columns in the
   table are the values returned by NCBI, not from memory. One PMID (42791346, Valleau et al.,
   EMBO J 2026) was verified with `efetch.fcgi` instead because it was fetched in the same call
   as its abstract; its title/journal/year were read back the same way. No PMID was written
   from memory. Where a claim could not be tied to a verified PMID, the row was not written.
2. **Mapped to ME49 gene IDs** from
   `/media/carruthers/mnt3/claude/repo/starplast/parascape/data/toxodb_identity.tsv` (8,843
   rows; `gene_id`, `gene_name`, `previous_ids`, `product`). Papers that used TGGT1_ IDs are
   noted in `aliases`; the ME49 equivalent was taken from the identity table or from a UniProt
   ToxoDB cross-reference, never guessed. RH88 numbering (TGRH88_) was treated as
   non-convertible. Unresolved IDs are marked `id_unresolved` (see below).
3. **Evidence is specific**: `method` names the technique (IFA in infected host cells, IEM,
   topology/selective permeabilization, host-side proximity labelling, export reporter,
   fractionation, co-IP with a host protein, GO IDA/EXP), `host_compartment` names the
   compartment reached, and `evidence_text` is a <=25-word quote or close paraphrase of the
   abstract.
4. **Strength is not upgraded**: `direct` = shown in the host or at the PVM by microscopy, IEM,
   topology or fractionation; `proximity` = host-side labelling only; `inferred` = host
   phenotype or host-protein interaction without a localization in the host compartment.
5. **Unlabeled is not negative.** No row was written as a negative from absence of evidence.
   `organelle_only_not_exported` is used 0 times: no paper read here explicitly tested export
   of a given protein and reported none for that protein. (Coffey et al. 2018, PMID 30377279,
   state that ASP5 substrates as a class "none appear to be exported" — that is a statement
   about a class, not a per-protein export test, so it was not turned into rows.)
6. **Conservative over large.** Dense granule / IVN-only proteins (GRA1, GRA2, GRA4, GRA9,
   GRA12), the parasite surface antigen SAG1, cyst wall proteins (CST1, MAG2), parasite-side
   kinases (WNG1, PV lumen), ASP5 (Golgi-resident) and GRA42/GRA43 (PV lumen; they position
   GRA35 but were not themselves shown host-facing) were left out of the literature rows.
   GRA83 was left out: it is PV-localized and its host effect (reduced host p65 translocation)
   is a phenotype, not a localization.

## What was searched

Date of all queries: 2026-11-04. Python 3 standard library from
`/home/carruthers/anaconda3/envs/starplast/bin/python`, requests paced at <= ~2.5/s with
tool/email identification.

### 1. Literature (PubMed E-utilities)

- `esearch.fcgi` (db=pubmed, sort=relevance, retmax=8) over **126 queries** in three rounds:
  - round 1 (72 queries): one per seed effector or machinery component — GRA15, GRA16, GRA24,
    GRA28, TgIST, GRA18, GRA6, GRA7, MAF1/MAF1b, GRA17/GRA23, ROP16, ROP18, ROP17, ROP5,
    ROP38, RON2/4/5/8, TgWIP, TEEGR, HCE1, GRA3, GRA5, GRA14, MYR1-4, ASP5, GRA44/45, PPM3C,
    GRA2/4/8/12/25/35/42/43/55/57/70/71/72, MAF1b, ROP1/2/4/39/40/47/54, toxofilin, TgPL1 —
    plus screen-shaped queries ("host nucleus effector BioID screen", "host cytosol export
    screen HA tag", "TurboID host interface proximity labelling", "MYR1-dependent exported
    effectors proteomics", "host mitochondria PVM association effector", "host plasma membrane
    secreted protein", review queries for the Hakimi, Saeij, Boothroyd, Bradley, Coppens
    labs).
  - round 2 (36 queries): gaps — GRA45, GRA44, GRA54, GRA61-63, GRA64, GRA80-83, GRA47/GRA60,
    WNG1, NSM, MAG1, host ESCRT/lysosome/Golgi/autophagy at the PVM, hyperLOPIT spatial
    proteomics, "exportome mass spectrometry host cytosol", "host nucleus proteomics infection
    parasite protein".
  - round 3 (18 queries): topology and mechanism gaps — GRA6/CAMLG, ROP17 host-cytosolic face,
    GRA15 TRAF, GRA7 Irga6, GRA16 c-Myc, Neospora GRA host nucleus, SPOTs/TOM70.
- `efetch.fcgi` (retmode=xml) for all **234 unique PMIDs** returned; abstracts were digested by
  regex for host-facing sentences, and every row traces to a sentence in that digest or in a
  full abstract read in full.
- `efetch.fcgi` (db=pmc) for three papers to recover gene IDs from the full text:
  PMC6495377 (HCE1 = TGGT1_239010, i.e. the same gene as TEEGR, TGME49_239010),
  PMC12367318 (GRA84 = TGME49_294970), PMC8071938 (MIC1/MIC4 paper, no Toxo gene IDs given).
- UniProt REST (`rest.uniprot.org/uniprotkb/<acc>.json`) for the nine accessions returned by
  QuickGO, to map them to ToxoDB IDs.

### 2. GO annotations (QuickGO REST)

```
https://www.ebi.ac.uk/QuickGO/services/annotation/search
  ?taxonId=5811&taxonUsage=descendants&goUsage=descendants
  &goId=GO:0033643,GO:0020005,GO:0030430,GO:0042025,GO:0020002,GO:0043657,GO:0020003,GO:0085036
  &evidenceCode=ECO:0000269,ECO:0000314,ECO:0000353,ECO:0000315,ECO:0000316,ECO:0000270,
                ECO:0007005,ECO:0007001
  &evidenceCodeUsage=descendants&limit=200&page=1
```

Returned **14 annotations only** (1 page), covering 9 UniProt accessions:

| accession | symbol | GO | evidence | reference |
|---|---|---|---|---|
| A0A140H546 | MAF1b1 | GO:0020005 | EXP, EXP, IDA | PMID:26920761, 28567444, 24781109 |
| A0A193AUK9, S8GJB7 | MAF1a1 | GO:0020005 | EXP | PMID:26920761 |
| B6KV60 | RON2 | GO:0020002 (host cell plasma membrane) | EXP x3 | PMID:16244709, 19247437, 21347343 |
| F6N7K6 | PSD1pv | GO:0020003 (symbiont-containing vacuole) | EXP x2 | PMID:22563079, 24429285 |
| O00834 / Q9XZH7 | MIC1 / MIC4 | GO:0044174 (host cell endosome) | IDA | PMID:33912181 |
| Q5Y808 | ROP4 | GO:0020005 | IDA | PMID:15470260 |
| V5BCL0 | PLP1 | GO:0020005 | IDA | PMID:19095897 |

Used in the table: MAF1a1 (PVM), RON2 (host plasma membrane), ROP4 (PVM, `id_unresolved`),
PLP1 (PVM, `needs_check`), MIC1/MIC4 (host plasma membrane, `inferred`, `needs_check`).
Not used: PSD1pv — PV lumen only, with no export or host interaction shown, so it does not
meet the definition. MAF1b1's GO row is covered by the literature rows for MAF1b.

### 3. Local deposits (read-only)

Inspected: `deposit_host_k562_rhoptry_screen.tsv`,
`deposit_host_k562_rhoptry_gene_evidence.parquet`, `deposit_organelle_surface_turboid.tsv`,
`deposit_host_hff_tg_infection.tsv`, `deposit_host_bmdm_baseline.tsv`,
`deposit_host_hepatocyte_pf_infection.tsv`, `secretome_partition.tsv`, `pvm_proximity.tsv`,
and their `parascape/datasets.py` entries.

- **`pvm_proximity.tsv`** (73 positives / 1,201 likely negatives; Cygan et al. 2021, mBio,
  PMID 34749525, verified) — host-cytosolic-face PVM miniTurbo labelling in human fibroblasts.
  All 73 positives are included as `PVM` / `proximity` rows. Those 73 overlap the literature
  rows for GRA7, GRA16, GRA24, GRA28, TgIST, TEEGR, GRA14, GRA17, GRA23, GRA35, GRA44, GRA45,
  GRA47, GRA54, GRA55, GRA60, MYR1-3, ROP16, ROP17, ROP54, MAF1a, MAG1 and others — the
  duplicate rows differ in `source` and `strength` on purpose, so the proximity set can be
  dropped or kept independently of the direct set. The 1,201 likely-negatives were **not**
  written as negatives.
- **`deposit_host_k562_rhoptry_screen.tsv`** / `..._gene_evidence.parquet` — the screen rows
  are *human* genes, so they give no Toxoplasma rows. Its paper (Valleau et al., EMBO J 2026,
  PMID 42791346, verified) does support a host-facing claim for the parasite MIC1/4/6 complex
  binding host N-glycans at the host plasma membrane; MIC1, MIC4 and MIC6 are recorded from it.
- **`deposit_organelle_surface_turboid.tsv`** — not used. Its baits face the *parasite* cytosol
  at the apicoplast/mitochondrion/ER outer membranes (bioRxiv 10.64898/2026.08.05.743015,
  preprint, no PMID), which is not host-facing.
- `deposit_host_hff_tg_infection.tsv`, `deposit_host_bmdm_baseline.tsv`,
  `deposit_host_hepatocyte_pf_infection.tsv` — host transcriptome/proteome responses, no
  parasite protein localization. Not used.
- `secretome_partition.tsv` — soluble-vs-vesicle partition of secreted material, keyed by
  TGGT1_; its own dataset note says it is not a list of what the parasite exports. Not used.

## Counts

By compartment x strength (rows, not unique genes):

| compartment | direct | proximity | inferred | total |
|---|---|---|---|---|
| PVM | 30 | 73 | 1 | 104 |
| host_nucleus | 11 | 0 | 0 | 11 |
| host_cytosol | 5 | 0 | 0 | 5 |
| host_organelle | 5 | 0 | 0 | 5 |
| host_plasma_membrane | 4 | 0 | 2 | 6 |
| moving_junction | 4 | 0 | 0 | 4 |
| PV_lumen_host_interacting | 0 | 0 | 5 | 5 |
| organelle_only_not_exported | 0 | 0 | 0 | 0 |
| **total** | **59** | **73** | **8** | **140** |

Status: 130 `ok`, **5 `needs_check`**, **5 `id_unresolved`**.

- `needs_check`: GRA66 (PMID 42818749 is a 2026 bioRxiv preprint, not peer reviewed);
  PLP1 (GO IDA to the PVM rests on an egress-permeabilization study, PMID 19095897, where the
  protein itself is micronemal); MIC1 and MIC4 (GO host-cell-endosome IDA comes from
  recombinant protein added to macrophages, PMID 33912181, not from infection);
  GRA25 (knockout host-immune phenotype only, PMID 24711568, no host-compartment localization).
- `id_unresolved`: three MAF1b/TgMAF1 rows (UniProt MAF1b1 is keyed to TGRH88_000130 and its
  ToxoDB cross-references list both TGME49_220950 and TGME49_279100, so the HMA-competent ME49
  paralog is not pinned; note that TGME49_220950 is also one of the 73 PVM proximity hits);
  ROP2 (the ROP2/ROP8 family — the identity table gives ROP2A = TGME49_215785 and ROP8 =
  TGME49_215775, while Sinai & Joiner 2001 studied "ROP2"); ROP4 (the identity table gives
  ROP4 = TGME49_295125 but the UniProt ROP4 entry cross-references TGME49_295110, which the
  identity table calls ROP7).

## Known gaps

- **GO contributes almost nothing**: 14 experimental annotations in host-related terms for all
  of taxon 5811. Host cell cytoplasm (GO:0030430) and host cell nucleus (GO:0042025) returned
  **zero** experimental annotations for Toxoplasma, although a dozen effectors are published in
  those compartments. The table is therefore literature-driven.
- **No unbiased Toxoplasma host-side exportome screen exists in what was found.** The only
  TurboID exportome paper the searches surfaced is *Theileria annulata* (PMID 38747635), which
  is a different organism and is not in the table. The nearest Toxoplasma equivalents are the
  PVM-face miniTurbo study (PMID 34749525, used) and the in-vacuole BirA* GRA screens
  (PMID 27486190 tachyzoites, PMID 32374791 bradyzoites) — those identify *vacuolar* GRAs, not
  host-compartment localization, so they were not converted into host-facing rows.
- Henry et al. 2024 (PMID 39292011, verified) name TgIST, GRA16, GRA24 and GRA28 as the four
  nuclear-targeted effectors blocking IFN-gamma; all four already have their own primary rows.
- MYR1-dependent-export membership (MYR4, GRA44, GRA45, PPM3C) is recorded as
  `PV_lumen_host_interacting` / `inferred`. These are machinery, not effectors; if the scoring
  task is "is this protein itself beyond the parasite plasma membrane facing the host", Einar
  may want to drop these five rows or re-tag them.
- Several well-known names were deliberately not given rows for lack of a verified
  host-facing primary claim in the abstracts read: ROP1, ROP38, ROP40, GRA8, GRA19-22, GRA55
  (except via the proximity deposit), GRA57/70/71, GRA83, WNG1, TgPL1, SAG1, AMA1, CST1.
- Full texts were only read for three papers; everything else rests on the abstract. A row
  whose `evidence_text` is thinner than Einar expects is a candidate for full-text checking,
  not a settled fact.
- Strain-of-origin is recorded only in `aliases`/`method` prose, not as its own column; several
  PVM/host claims were made in RH or GT1 and transferred to the ME49 ID by gene identity.
