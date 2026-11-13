# 53 · One space per organism, then links between pathogen, vector and host

**Status: open (designed 2026-09-26). 2026-11-08: WP0, WP1, WP2 and WP8 done, WP3 (Hs, Mm) done to its acceptance, WP6 has Cp and Pb, WP7/WP9 partly; rat, cat and vectors (WP4-5) not started -- they need data acquisition, much of it behind the VEuPathDB login (item 38).**

Crash recovery: all implementation commits through `05650fd` survived and were pushed. The 45
focused recovery checks passed. Documentation CI passed; test CI had 4,138 passes and one stale
packaging assertion that disallowed subpackages. The assertion now checks that every ParaScape
source package is included in the single distribution, including `parascape.spaces`.
All 50 publishing/release checks pass with this correction. Follow-up remote CI remains pending.

### Continued implementation

The user authorized sustained work on 2026-09-27, prioritizing token efficiency. See
`WORK_LOG.md` for completed milestones, pending dependencies and actual blockers. Current scope:
WP8's lossless per-space calibration publishing and WP1's versioned pack/build framework.
Acceptance: publishing a subset preserves untouched results/provenance; corrupt existing files
are not overwritten; a pack round-trips hashes, rejects tampering/unsafe paths, installs atomically,
and refuses to build over lost identifiers, columns or measured values.

Pf display fix implemented: per-space display recipes resolve 118 Pf features for all 5,720 genes.
`inference_spec` explicitly retains the recipe used by the shipped calibration. The executed
`notebooks/pf_display_layout_2026_09_27.ipynb` compares 17 target-specific/no-target strategy groupings
and both species' compatibility input matrices against pre-change fingerprints: all unchanged.
Only Pf coordinates and layout metadata changed; every edge array, both node tables, Tg layout,
slots and calibration were preserved. Focused checks: 156 passed. UI/docstring checks: 305 passed,
1 skipped offscreen; the real shared-OpenGL test then passed under Xvfb. UI/optimizer checks:
219 passed. Search/leakage/graph/strategy checks: 405 passed, 2 skipped. A final fingerprint check
confirmed all 17 groupings and both compatibility matrices unchanged. The clean 53.4 MB wheel
passes validation and package modules parse under Python 3.10 syntax. The earlier full-suite
GUI/GPU baseline failures remain documented below; no full-suite rerun or release was performed.

WP8 publication implemented: updates merge by organism/strategy, retain each result's original
sweep provenance and atomically replace the file. Corrupt existing files fail before writing.
Generated pages include all published spaces. Calibration/registry checks: 37 passed; shipped
calibration values were not changed. Remaining WP8 strategy applicability work is still open.

WP1 framework implemented: `packs.py` verifies complete ZIP manifests and SHA256 values, rejects
unsafe members and unreviewed/local-only license codes, downloads with bounded size and timeout,
and activates immutable versions atomically. `spaces.build_pack` and `scripts/build_space.py`
validate canonical unique IDs, graph order/indices and per-cell preservation against the previous
table. Registry paths resolve active packs without downloading or silently falling back. CLI and
usage are in `docs/space-packs.md`. 380 relevant checks passed and a clean 53.4 MB wheel includes
the new subpackage. Source-specific builders, acquisition and published download URLs remain open.

### Implementation pass, 2026-09-27

R3 implemented and validated with the baseline test failures recorded below. All 161 datasets now declare their organism;
the four host deposits declare Hs/Mm. The mixed protein cache is replaced by
`hs_host_proteins.parquet` (20,989 rows) and `mm_host_proteins.parquet` (15,590 rows). No identifier,
name, measurement or bridge endpoint was lost. The 54 identity-only rows were resolved from local
UniProt mapping files and human bridge sources, not guessed from names. The executed migration is
`notebooks/split_host_tables_2026_09_27.ipynb`; its manifest records checksums and coverage. The
script refuses unassigned/conflicting rows, unknown columns and overwriting later acquisitions.
Readers and builders select one host species at a time, including the slot viewer and generated
atlas. `slots.json` and the generated atlas remain byte-identical. Full host gene spaces, packs and
the Space menu remain subsequent work packages.

**Correction to the Pf follow-up:** the handoff's claim that strategies are unaffected is false.
`Context.blocks()` calls `embedding.default_spec()`. A proposed per-species recipe expands the
display from 20 source columns to 118 resolved features, but changes inference grouping from 50
to 62 blocks without a target, and changes it for every shipped calibration holdout. The experiment
is recorded in `notebooks/pf_layout_dependency_2026_09_27.ipynb`; production recipes/layouts were not
changed during that first audit. The subsequent display fix above explicitly separates display
and calibrated inference recipes. Do not reuse the old scores if the inference groupings change.

Validation (Python 3.12, pandas 3, offscreen Qt, local GPU environment):

* Focused host/deposit/dataset/organism/slot checks: 246 passed. The new migration and reader
  regressions also pass, including identity-only rows, zero/false measurements, ambiguous species,
  and refusing to overwrite later acquisitions. Leakage/search/source checks: 178 passed, 2 skipped.
* Full suite: **4,145 passed, 9 skipped, 7 failed** in 19 minutes. All seven failure cases reproduce
  on untouched `f02d075`: import-column counts, gated point sizes, persisted display settings,
  sprite caching, CPU UMAP fallback with cuML available, offscreen shortcut focus, and attention
  colors. These are existing GUI state/test-order and GPU-environment issues, not a green suite.
  A baseline run of `test_app_controls.py`, `test_display.py`, `test_build_graph.py`,
  `test_help_search.py`, and `test_app_smoke.py` with `--randomly-seed=42` reproduced six of them
  (11 failed, 356 passed). The seventh is reproduced by running
  `test_searching_an_exact_gene_id_selects_it` immediately before
  `test_a_gated_set_recedes_the_rest_of_the_map_without_recoloring_it` in `test_app_controls.py`
  with `-p no:randomly`: the selected gene retains its larger point size.
* `add_deposits.py --dry-run` found no loss. The migration checked every retained cell and parquet
  round trip; all 36,579 protein identifiers and all human bridge endpoints remain present.
* Regenerating the slot catalogue and atlas produced byte-identical files. Parasite node tables,
  graphs and calibration are unchanged. Both species' slot audits have zero orphan columns,
  double ownership or missing citations.
* A clean wheel builds and passes `check_wheel.py` (53.4 MB). Negative checks reject a wheel with
  the old mixed cache or a missing mouse cache. `release.py check`, Python 3.11 syntax parsing and
  `git diff --check` pass. No release/version bump is part of this milestone.

* **R0 done:** `parascape/organisms.py` declares Tg and Pf, and `tests/test_organisms.py` holds it to every literal it replaces.
* **R1 done:** these now read the registry:
  * `slots.SPECIES_TABLES`, `SPECIES_PREFIXES` and `SPECIES_BRIDGE_TABLES`;
  * `app.SPECIES` and `DEFAULT_SPECIES`, and the gene record link;
  * `strategies.Context.shipped`, its graph file, `other()` (the space's partner) and `_guess_organism` (registry detection first);
  * the calibration `TARGETS`/`NUMBERS`;
  * the slot generator's `STAGES_BY_ORGANISM` (regenerated `slots.json` byte-identical).
* **R2 started:** `test_no_file_gains_an_organism_literal` is a ratchet. It measured 156 bare "Tg"/"Pf" literals in 32 files, and any file that gains one fails. Lower the numbers as files move to the registry.
* **R3 implemented:** explicit dataset organisms and separate Hs/Mm protein tables (see above).

The user, 2026-09-26: "we should also make the host and vector datasets more comprehensive!
(plasmodium, cryptosporidium, human, mouse, feline, anophelus, rat, datasets and informationslots
need to be added) each in its own space first then we can start connecting pathogen to vector and
host."

## STATE

Measured on nightly @446e20e.

* **There is no organism registry.** About 150 `"Tg"`/`"Pf"` literals sit in 50 files, and at
  least five separate tables map species to files or ID prefixes. Three naming schemes coexist:
  `Tg/Pf`, `tgon/pfal/cpar` (`localization.py:34`), and `"Toxoplasma gondii"` (`app.py:277`).
* **"host" is two species in one table.** `host_proteins.parquet` has 36,579 rows keyed on UniProt
  accession. Its human columns (fibroblast_tpm, 19,087 rows) and mouse columns (brain_tpm 9,763,
  bmdm_tpm 15,437) never overlap. That breaks instruction 39's rule of one table per species.
* **Host slots belong to the parasites.** `slots.json` has 290 slots. Tg has 116 gene, 26 host_gene,
  10 pair and 3 metabolite slots; Pf has 99, 24, 9 and 3. So "Tg_host transcriptome · human
  fibroblast" is a Toxoplasma slot.
* **Two lookups silently default to Tg.**
  * `datasets.organism_of` (`datasets.py:2670`) returns "Tg" for every `host_*` dataset.
  * `strategies._guess_organism` treats anything that is not PF as Tg.
* **How Pf was added:** a separate builder (`plasmodium.build_all`) plus `pf_graph.py` copying Tg's
  graph construction, mirrored slots (`_pf_mirror`), a guarded rebuild (`build_plasmodium.py`), and
  ID-space tests (`test_pf_graph.py`, `test_plasmodium.py`). Every new space copies this template.
* **Public data only.** The VEuPathDB service API needs a key we do not have. The public
  `common/downloads/Current_Release` files need none; `archive.discover_veupathdb` already walks
  them.

Where each organism assumption lives (file:line):

| Concern | Where |
|---|---|
| Table/prefix maps | `slots.py:44-66`, `table_organism` :69 (first prefix match wins) |
| App | `app.py:277-310` SPECIES/load; Species menu :1498; `open_species` :2644; record link :4243; `EDGE_TYPES` :60, `COLOR_MODES` :83 |
| Strategies | `strategies.py` `shipped`, `_guess_organism`, `Context.shipped`, `family_of`, graph file, `other()` (the other one of Tg/Pf), synthetic Pf partner |
| Orthology | `strategy_catalog.py` `_source_default`, `_through_orthologs` (joins on `orthogroup`; Tg and Pf share 2,526 OG6 groups) |
| Leakage | `search.py` `SAME_QUANTITY`, `LAYER_SOURCES` (keyed on column name only); `leakage.py` |
| Slot generator | `generate_slot_table.py` `ORGANISMS`, `HOST_CONTEXTS_TG/PF`, `HOST_FAMILIES`, `HOST_COLUMNS`, `STAGES_BY_ORGANISM`, `all_slots`, `_rows` |
| Datasets / deposits | `Dataset` has no organism field; `organism_of` works from the key prefix; `Deposit.organism` is one of {Tg, Pf, host} |
| Host | `host.py` table names, `uniprot_index` (human/mouse), GTEx, FANTOM5, `TISSUE_REFERENCES` |
| Calibration | `calibrate_strategies.py` TARGETS/NUMBERS/`--organism`; `calibration.write` **overwrites the whole file** |

### Raw data acquired (2026-09-26)

5.9 GB is downloaded under `<PARASCAPE_DATA>/spaces/<code>/<kind>/<source>/`.

* **What is there:** Hs, Mm, Rn, Fc, Ag, As, Cp and Pb, plus `shared/` (gene2pubmed, gene_orthologs,
  BioGRID, BioGRID PTMs and ORCS, IntAct).
* **How each deposit is recorded:** URLS.txt holds the resolving query per file, and SHA256SUMS
  verifies it (0 mismatches).
* **Per space:** a `MANIFEST.json` with release, licence code, sanity check and target slot family.
* **Also there:** `spaces/README.md`, and the scripts and an annotated notebook under
  `spaces/_acquisition/`.
* **Structure:** mean pLDDT per model from the AlphaFold DB search API for every space. There is no
  per-proteome summary file.
* **Hs highlights:** Ensembl 116 plus UniProt, GTEx v10, HPA 25.1 (HPA is CC BY 4.0 now, not SA),
  Bgee 15.2, DepMap 24Q4 (from figshare+; the portal is behind Cloudflare), COMPARTMENTS, OpenCell,
  STRING 12.0, HuRI, BioPlex, CORUM (non-commercial), PaxDb, iPTMnet, gnomAD v4.1, GOA, and 9
  infection series (Tg, Pf, Pv, Pb, Cp).
* **Mm:** Tabula Muris Senis bulk, IMPC DR24 viability and phenotypes, ImmGen, and 9 infection
  series.
* **Rn:** BodyMap, RatGTEx.
* **Fc:** cat small-intestine epithelium transcriptome and phosphoproteome under Tg (PMIDs 37491273,
  38003154).
* **Ag:** MozAtlas, blood-meal time course, PRIDE salivary gland/hemolymph/saliva, 5 infection
  series.
* **Cp:** Walzer 2024 single-cell atlas, oocyst proteomes. No genome-wide screen exists.
* **Pb:** PlasmoGEM Bushell 2017 and Stanway 2019 (the real tables), Russell 2023, Malaria Cell
  Atlas.

Traps it found:

* **VEuPathDB is closed to anonymous downloads.** Release 71 public downloads return 404 and the
  service returns 401, so the Cp/Pb/Ag/As references come from Ensembl Genomes r63.
* **Two PlasmoGEM files in the old dataset tree are fakes.** They are 1.8 kB HTML "preparing to
  download" pages, and `reference/plasmodb/pb_transfer/` does not exist. The real tables are under
  `spaces/Pb/essentiality/`.
* **The FANTOM5 path `host/fantom5/24670764/` does not exist.** The file is now in
  `spaces/Mm/expression/FANTOM5_E-MTAB-3579/`.
* **Several id spaces are split or retired:**
  * Cat: three id spaces (retired Felis_catus_9.0, and ENSFCTG).
  * *A. gambiae*: NCBI moved to AGAMI1_ ids with no AGAP cross-reference; bridge through UniProt.
  * *A. stephensi*: four id spaces.
  * Rat: STRING, COMPARTMENTS and PaxDb use old ENSRNOP ids.
  * Cp: no open IOWA-ATCC↔cgd map.
  * Pb: some tables use old 6-digit PBANKA ids.
* **gene2pubmed barely covers the parasites** (Cp: 29 PMIDs).
* **Junk to delete by hand:** `spaces/As/abundance/PXD001647/` holds 48 kB of unrelated Hydra and
  Daphnia mzTabs, recorded as `rejected_unusable`.

## WHY IT MATTERS

Adding a species today means editing about 30 files. Every new species multiplies the hard-coded
literals and the chances for silent Tg defaults. Host and vector biology can only become inference
targets once each species has its own table, slots, leakage families and calibration.

## WHAT TO DO

### WP0, registry and refactor (serial; everything else waits for it)

Add `parascape/organisms.py` with one `Space` per organism. A space declares:

* `code` (Tg, Pf, Pb, Cp, Hs, Mm, Rn, Fc, Ag, As), `species`, `reference`, and `kind`
  (parasite, host or vector);
* `gene_regex`, `alias_prefixes`, `id_source`, `identity_file` and `symbol_prefix`;
* the file paths for `nodes`, `graph` and `mentions`;
* `orthology` and `orthomcl_abbrev`;
* `hosts`, `vectors` and `partner` (Tg↔Pf, Pb→Pf);
* `contexts` (stages or tissues), `record_url` and `pubmed_query`;
* `edge_labels` and `colour_modes`;
* `calibration` tier, `distribution` (wheel or pack) and `builder`.

Helpers: `get`, `spaces(kind)`, `detect(ids)` (a majority full-match vote), `nodes_path`,
`graph_path` and `display_name`.

Refactor in steps. The suite stays green after each one.

* **R0.** The registry for Tg and Pf only, plus `tests/test_organisms.py`. The test asserts that the
  registry reproduces every current literal.
* **R1.** Replace each constant with a derived alias that keeps its old name.
  `Context.other()` → `organisms.get(code).partner`.
* **R2.** Tests read the organism set from the registry. Add a lint test that allowlists the
  remaining `"Tg"`/`"Pf"` literals, and shrink the allowlist over time.
* **R3.** Add `Dataset.organism`. Split `Deposit.organism` "host" into "Hs" and "Mm".
* **R4.** The slot generator reads the registry. Acceptance: `slots.json` is byte-identical.

### After WP0, in parallel worktrees (owned files do not overlap)

| WP | What | Owns | Accept |
|---|---|---|---|
| 1 | Space framework and on-demand packs | `spaces/__init__.py`, `spaces/veupath.py`, `packs.py`, `paths.py`, `scripts/build_space.py`, `check_wheel.py` | A pack round-trips its sha256; a tampered pack is refused; a lost column refuses the build |
| 2 | Slot generator per space | `generate_slot_table.py`, `scripts/slot_spaces/*.py` | Tg/Pf output byte-identical; each new space emits prefixed slots |
| 3 | Human (Hs), mouse (Mm) | `spaces/vertebrate.py`, `hs.py`, `mm.py`, `dataset_registry/hosts.py` | ENSG/ENSMUSG unique; ≥95% of reviewed UniProt resolved; GTEx matches today's values |
| 4 | Rat (Rn), cat (Fc) | `spaces/rn.py`, `fc.py` | Every column Fc takes through orthology carries `derived_from`; each empty slot has a verdict |
| 5 | Anopheles gambiae (Ag), stephensi (As) | `spaces/anopheles.py`, `dataset_registry/vectors.py` | AGAP IDs; As assembly recorded; blood-meal slots `separate` |
| 6 | Cryptosporidium (Cp), P. berghei (Pb), then P. vivax (Pv) | `spaces/crypto.py`, `pberghei.py`, `dataset_registry/apicomplexa.py` | OG6 present; Cp LOPIT loaded; Pf `pb_transferred_*` equal to the projection of the Pb space |
| 7 | Leakage per space | `search.py` families, `LAYER_DERIVED_FROM`, `leakage.py` | Attack tests: GTEx→HPA, BioGRID→STRING, DepMap→ORCS, cross-space transfer |
| 8 | Strategies across N spaces | `strategies.py` (`partner`, `bridge`, `applicable`), `calibration.py` (merge per space) | Publishing one space leaves the others byte-identical; strategies that do not apply are recorded as N/A |
| 9 | UI | `app.py` Space menu (Parasites / Hosts / Vectors, a download for any space not installed) | Every installed space opens offscreen |

**Phase 2** follows once the spaces exist:

* orthology bridges (OG6, Ensembl Compara, OrthoDB; 1:n relations kept explicit);
* host–pathogen PPI re-keyed to `(space_a, id_a, space_b, id_b)`, using HPIDB and the IntAct subset;
* vector–parasite links (stage↔tissue);
* `Context.bridge(code, kind)`;
* a "Linked spaces" dock that carries a selection across windows.

### Data per space

Resolve every identifier through an API; never type one. Licence codes: **O** open, **SA**
share-alike, **V** verify the terms first, **X** local build only.

* **Hs.**
  * GTEx v10 (O). HPA tissue, cell type and subcellular (SA). DepMap summaries (O).
  * UniProt/GO, with experimental codes only for targets (O).
  * STRING, BioGRID and IntAct (O). HuRI, BioPlex, OpenCell and CORUM (V). AlphaFold (O).
  * iPTMnet (V). **PhosphoSitePlus (X)**. gnomAD constraint (V). Bgee (O). gene2pubmed (O).
* **Mm.**
  * Tabula Muris Senis (O). FANTOM5, already present. ENCODE (O). ImmGen (V).
  * IMPC (O). MGI (V). BioGRID ORCS (O). STRING. AlphaFold.
  * Tg and Pb infection series (GEO).
* **Rn.** Rat BodyMap (O). RatGTEx (V). Bgee. RGD (V). STRING. AlphaFold. Tg-infected brain (GEO).
* **Fc.** Little measured data. Use Bgee or Expression Atlas where present, AlphaFold and STRING.
  Look for feline intestinal organoid series with Tg sexual stages. Most slots will be transferred
  and declared empty, each with a verdict.
* **Ag.**
  * Tissue atlases. Blood-meal time courses.
  * Plasmodium-infected midgut and salivary gland (GEO).
  * Midgut and salivary proteomes (PRIDE). Ag1000G (V). Insecticide resistance (V). AlphaFold.
* **As.** Lab transmission data (far less).
* **Cp.**
  * CryptoDB IOWA-ATCC, with `cgd` aliases. OG6. LOPIT, already in the tree.
  * Life-cycle scRNA and time courses (GEO). Oocyst and sporozoite proteomes (PRIDE).
  * There is no known genome-wide CRISPR screen; search for one before assuming.
* **Pb.** PlasmoGEM screens get a native home instead of existing only as transfers into Pf.
* **Pf (expand).** Malaria Cell Atlas (V). MalariaGEN summaries (V). PRIDE stage proteomes.

**Size.** The wheel is 52.6 MB (PyPI limit 100 MB).

* **In the wheel:** Tg, Pf, Cp and Pb (under ~65 MB).
* **As packs** in `user_cache_dir()/spaces/<code>/<version>/` with a sha256 manifest and per-column
  licences: Hs (~15-25 MB), Mm, Rn, Fc, Ag, As and Pv.

**Calibration.**

* **Tier A (full sweep):** Tg, Pf and Hs (Hs capped at ~8k rows).
* **Tier B (3 seeds, defaults + best):** Mm, Pb, Cp and Ag.
* **Tier C (self-test only):** Rn, Fc, As and Pv.

## HOW TO KNOW IT WORKED

* WP0: the full suite passes; `slots.json` is byte-identical; the registry reproduces every old
  literal; the lint test passes.
* Each space: the ID-space tests from `test_pf_graph.py` carried over; each slot filled, or empty
  with a verdict; its leakage attack tests pass; its self-tests and calibration published without
  touching any other space's numbers.

## TRAPS

* **The Pf display previously used only 20 features.** The display follow-up now uses its own
  slots (145 source columns, 118 resolved features). `inference_spec` preserves the calibrated
  grouping separately: broad raw-column coverage never proved that grouping changes were harmless.
  The two executed notebooks above record both the dependency and the final invariance checks.

* R3 resolved the mixed human/mouse host cache. Host readers and deposits must continue to select
  one explicit organism; full host gene spaces remain pending.
* `calibration.write` now merges measured organism/strategy entries. Unswept entries keep their
  original provenance; changed strategy recipes still require a fresh calibration.
* HPA consensus integrates GTEx; BioGRID ORCS and OGEE reuse DepMap; STRING's combined channel
  contains BioGRID and IntAct. These are one family each: holding one out must hold out the others.
* The Pf host contexts in `generate_slot_table.py` already name Anopheles midgut and salivary gland.
  Those become bridge slots into the Ag space, not Pf columns.

## Progress through instruction 64 (2026-11-08)

* **WP0 (registry):** done earlier (R0-R4); the organism-literal ratchet now also holds every file touched in 64
  (literals replaced by `parascape.organisms` constants, 2026-11-07).
* **WP1 (space framework and packs):** done. `packs.py` (hash-checked install, tamper and size refusal, interrupted
  activation), `spaces.build_pack` (lost IDs/columns/values refused; `extra=` members for notices and link tables),
  `space_catalogue.py` (ready / installed / installable / unavailable, never calibrated without a track record).
* **WP3 (Hs, Mm):** done to its acceptance. Hs: 58,988 unique ENSG genes (GTEx v10, GENCODE 39); **19,268 of the
  19,400 reviewed human proteins with an Ensembl cross-reference resolve to a gene in the universe (99.3 %, bar
  95 %)**; GTEx values replay exactly from source text (foundation, 64.21); HPA annotations attached; distributable
  pack with per-column licences (O / SA) and notice. Mm: 18,428 unique ENSMUSG genes (Bgee, Ensembl Compara), human-
  ortholog locations kept as inferred, pack built (64.22). Host pilots measured on host truth (64.28/64.29).
* **WP7 (leakage per space):** the category and measured-redundancy rules apply to any space through the harness;
  per-space attack tests (GTEx->HPA etc.) not yet written.
* **WP9 (UI):** host spaces listed in the Species menu with their state and installable from there; a host gene opens
  from its pack through a parasite gene's links (64.23/64.24). A host space as a full map window is not built.
* **Phase 2:** organism-qualified bridges with relation types (`bridges.py`), counterpart resolution in its own space.
* **WP2 (slot generator per space): built 2026-11-08.** `scripts/generate_slot_table.py` no longer hardcodes two
  organisms. A further space declares its own biology in `SPACE_SLOTS` -- stages, which axes are stage-resolved, the
  default stage, and the contexts that belong to another parasite laboratory -- and `_space_mirror` asks the base
  Toxoplasma questions of it under the same rules as the Plasmodium mirror. Nothing is translated: the question
  shape carries across apicomplexans, the stage list never does.
  * **Cryptosporidium parvum: 63 slots**, graded against its own 3,923 genes (`SPACE_TABLES`), not against the
    Toxoplasma gene count. Nine carry real measurements: membrane topology 100 %, orthogroup 100 %, sequence
    basics 100 %, identity 100 %, conservation 100 %, domain content 65.6 %, measured hyperLOPIT location 28.1 %,
    single-cell stage share 37.1 %, derived stage label 37.1 %. Catalogue: 360 slots (was 297).
  * **Fitness slots are omitted for Cryptosporidium on purpose:** `fitness_stages` is empty because no pooled,
    counted genome-scale mutant library exists in it, and a slot for an experiment nobody can run is a denominator
    that makes every coverage figure read low for no reason -- the rule that already shortened the Plasmodium list.
  * **Two slots are declared rather than mirrored:** measured hyperLOPIT location (Guerin 2023) and the
    single-cell stage share (Tandel 2019). The Toxoplasma slots that would mirror onto them name their own study,
    which `STUDY_SPECIFIC` rightly refuses, so they are declared with their columns instead of smuggled past it.
  * **Coverage counts answers, not cells:** `_space_table` reads empty strings as missing, because these tables
    store "no domain" and "no location" as `""`. Left alone, domain content read 100 % covered when 1,348 of
    3,923 genes have no domain at all.
  * **A further space is declared, never registered.** Registering a third space inside the generator changed
    what every later caller saw -- two registry tests began to fail by test order alone -- so `species_of` reads
    the species from the declaration and `organisms.SPACES` keeps holding Tg and Pf. Tests
    `test_the_catalog_has_both_organisms_and_stable_prefixed_keys` and
    `test_every_slot_has_an_explicit_organism_scoped_query` now assert the property rather than a fixed pair.
  * **Plasmodium berghei gets no slots, on purpose.** The Plasmodium arm already queries berghei, vivax, knowlesi
    and yoelii by design, because the genome-scale fitness and liver-stage data are largely berghei. A parallel Pb
    table would ask every Plasmodium question twice and split the data that fills it.
* **WP8 (strategies across N spaces): built 2026-11-09.** `parascape/strategy_applicability.py` reads what a
  space holds in the words strategies use to declare what they need ("a categorical column", "two or more edge
  layers", "the other organism's table") and answers, per strategy, whether it could run and which requirement
  is missing. Nothing is run: this is applicability, never a measure of how well a strategy works, which stays
  the track record. `scripts/build_strategy_applicability_notebook.py` ->
  `notebooks/strategy_applicability_2026_11_08.ipynb`, `results/strategy_applicability_2026_11_08/`.
  * **Of 39 strategies: Tg 39, Pf 39, Cp 27.**
  * **Every one of the 12 blocked in Cryptosporidium is blocked by missing EVIDENCE, not by the organism:**
    9 wait on a graph (no edge layer, no interaction, crosslink/IP-MS, structural-similarity or co-mention
    layer), 2 on a literature-attention column, 1 on a counterpart table. So the build order is: the
    Cryptosporidium graph first, then its Toxoplasma counterpart table; most of the 12 open with the first.
  * A requirement the checker has never seen is **allowed**, not refused, so a new strategy cannot be silently
    blocked everywhere by wording the checker does not know.
  * A naming trap, found and fixed: the shipped graphs call interaction evidence `xlms` and `ip_ms`, so a check
    looking for "crosslink" alone reported Plasmodium as having no interaction layer when it has both. The
    layer-name fragments are now spelled out in `INTERACTION_LAYERS`, and a test asserts both shipped parasite
    spaces run every strategy -- a block there is a bug in the checker, not a fact about the space.
  * Tests: `tests/test_strategy_applicability.py` (5).
* **WP5 and the Pv half of WP6, built 2026-11-09** in the same pass as WP4, because the same wrong assumption
  had blocked all three. `parascape/spaces/anopheles.py` (`Ag`, **kind `vector`** -- the parasite completes its
  sexual cycle in the mosquito, so ookinete, oocyst and sporozoite happen inside this animal) and
  `parascape/spaces/vivax.py` (`Pv`).

  | Space | Kind | Genes | Build |
  |---|---|---:|---|
  | Rn | host | 23,141 | Ensembl 113 mRatBN7.2 |
  | Fc | host | 19,588 | Ensembl 113 Felis_catus_9.0 |
  | Ag | vector | 12,614 | Ensembl Metazoa 60 AgamP4 |
  | Pv | parasite | 6,523 | Ensembl Protists 63 PvP01 v2 |

  **P. vivax is kept apart from falciparum, not folded into it:** vivax forms hypnozoites and falciparum does
  not, so a slot about dormancy means something different in each and one table answering for both would be a
  false equivalence.

  **A pattern that refused real genes:** the first `Pv` id pattern was digits only, and the build names its
  organelle genes apart -- `PVP01_API05000`, `PVP01_MIT02100` -- so 33 genes (30 apicoplast, 3 mitochondrial)
  failed validation. Falciparum already handles exactly this; the vivax pattern now matches it.

  All four are listed by `space_catalogue` as **unavailable**: gene build only, no pack.

  **Both new spaces now hold measurements too, from arms of series already read for other species.** I had
  said these were "already on disk"; only the falciparum halves were, and the other halves were downloaded:
  * **Pv: liver-stage expression for all 6,523 genes** (`expr_pv_liver`), from the vivax arm of GSE220039 --
    the same humanised-liver series read for falciparum. **Controls only**: it pairs them with a PI4K-inhibitor
    arm, which is a drug experiment and a different question. 9.05 M parasite reads over three controls, and
    counts renormalised to the parasite total because the parasite is a minority of a liver biopsy.
  * **Ag: midgut expression at four timepoints** (`anopheles_midgut_36h/48h/4d/7d`) for all 12,614 genes, from the
    Anopheles arm of GSE284537 -- the mosquito's own genes responding to infection, where the parasite arm of
    the same files filled the Plasmodium mosquito-stage slot. Control (dsGFP) arms only, pseudo-bulked.

  883 tests pass including the Species menu, and 750 again after the measurements were added.

  **Renamed 2026-11-10 before it could bite.** Those columns were first called `expr_midgut_*` -- exactly the
  names the Plasmodium table carries for the PARASITE arm of the same files. One name meaning two species'
  expression is how a join goes quietly wrong, and it is the same collision the mouse and human PaxDb columns
  hit the same night (`host slot 'mouse brain' matches more than one organism`). The mosquito's are now
  `anopheles_midgut_*`; packs rebuilt, categories updated, 809 tests pass.
* **WP5 completed 2026-11-09 with the second vector:** `parascape/spaces/stephensi.py` (`As`, Ensembl Metazoa
  AsteS1, 13,336 genes). A separate space rather than a strain of gambiae: stephensi is a different species
  that carries malaria across South Asia and has spread into Africa, where it breeds in urban containers
  rather than rural pools -- a question about where a vector lives, or which insecticide reaches it, has a
  different answer in each. The same reason vivax is kept apart from falciparum.
  Gene build only; no stephensi measurement is on disk, which `what_each_space_needs.tsv` records.

**The catalogue now holds nine spaces** -- Tg, Pf, Hs, Mm, Rn, Fc, Pv (parasite or host) and Ag, As (vector) --
against two at the start of the week. 798 tests pass.

* **WP7 done 2026-11-09: the leakage rules now apply per space, and five attacks prove it.**
  `tests/test_space_leakage.py`. The rules were written once in the harness and claimed to hold in any space;
  that claim was worth only the attacks it survived, and the first attack found a real hole.

  **The hole: all five new spaces carried NO category on any column.** The closure excludes evidence BY
  category, so a column with none is not excluded by anything -- a space nobody annotated has no leakage
  protection at all, silently. Fixed with rules in `OVERRIDE_PATTERNS`, and the first test asserts every
  column of every further space is categorised, so the next space cannot arrive unprotected.

  The five attacks, each one either already possible here or a rename away:
  1. *An uncategorised column* -- caught, as above.
  2. *Coordinates and annotation as features*: `seqid`, `start`, `end`, `strand`, `biotype`, `gene_name` are
     **identity**, the one category that is never a feature. A model reading a start coordinate is reading
     gene neighbourhood, which is how gene order leaks into an answer about biology.
  3. *Annotation transferred along an ortholog, used to predict what it came from*: the mouse space copies
     locations from human orthologs. Those columns are **place**, so a question about location excludes them
     -- they ARE that location, copied. A pattern bug let `location_from_human_ortholog` through while
     catching `hpa_location_from_human_ortholog`; the test caught it before it was committed.
  4. *The ortholog link itself* is **conservation**, the same kind of evidence as an orthogroup, so a question
     about conservation excludes it.
  5. *An unknown category name* raises rather than excluding nothing, which is the silent failure the
     vocabulary replaces.

  1,059 tests pass over categories, harness, redundancy, claims, strategies, catalogue and slots.

* **Packs built 2026-11-09 for all five new spaces** (`scripts/build_space_packs.py` ->
  `results/space_packs_2026_11_09/`, with `SHA256` and `packs.json`). Each was installed from its archive and
  its contents checked, so the round trip is tested rather than assumed.

  | Space | Genes | Pack | Withheld |
  |---|---:|---|---|
  | Rn | 23,141 | 733 KiB | -- |
  | Fc | 19,588 | 625 KiB | -- |
  | Ag | 12,614 | 198 KiB | the four `expr_midgut_*` |
  | As | 13,336 | 211 KiB | -- |
  | Pv | 6,523 | 104 KiB | `expr_pv_liver` |

  **A pack may only carry columns whose licence is verified open**, which is what `packs.build` enforces. The
  Ensembl builds, Bgee calls and Compara orthologs are open and go in. The five columns derived from GEO series
  do NOT: a GEO submission states no licence of its own, and deciding that a depositor's data may be
  redistributed is a judgement about someone else's rights rather than a technical detail. Those columns stay
  in `results/` for local use and each pack's `NOTICE.txt` names what was withheld and why. **Einar to confirm**
  -- if GEO deposits are fine to redistribute it is a one-line change.

  Two faults the packing exposed:
  * Every new space **declared a graph file that does not exist**, because no edges have been computed. The
    attribute is now unset rather than naming a missing file: `build_pack` validates whatever a space declares,
    so a claimed graph is a promise the pack cannot keep.
  * The catalogue **hardcoded the three newest spaces as `unavailable`**, so three real packs stayed hidden.
    It now reads every pack space's state from disk (`PACK_SPACES`), and `install` accepts them, so the next
    pack appears without a code change. All nine spaces report their true state; 813 tests pass.

* **WP3 acceptance applied to all five new spaces, 2026-11-09** (`scripts/check_space_acceptance.py` ->
  `results/space_acceptance_2026_11_09/`). The Hs/Mm acceptance was: at least 95 % of that species' reviewed
  UniProt proteins carrying an Ensembl cross-reference resolve to a gene in the universe. It asks the only
  question that matters about a gene universe -- is anything real missing -- from a source that did not build it.

  | Space | Genes | Check used | Rate | Verdict |
  |---|---:|---|---:|---|
  | Rn | 23,141 | reviewed proteome (8,239 proteins) | 100.0 % | passes |
  | Fc | 19,588 | reviewed proteome (234 proteins, 94 resolvable) | 100.0 % | passes, on a thin sample |
  | Pv | 6,523 | independent gene list | 95.1 % | passes, barely |
  | **Ag** | 12,614 | independent gene list | **90.8 %** | **FAILS** |
  | As | 13,336 | none available | -- | cannot be checked |

  **The measure does not transfer to every species, and that is a finding.** It needs a reviewed proteome, and
  the less-studied species barely have one: As has **9** reviewed proteins and Pv **11**. A 95 % bar over nine
  proteins says something about nine proteins, not about a universe of 13,336 genes. So each space gets
  whichever check is meaningful -- the reviewed proteome where one exists, otherwise coverage of a gene list
  produced by someone else's pipeline (the genes a published series quantified).

  **Ag fails, and the acceptance test is what caught it.** The authors of GSE284537 quantified 13,845 AGAP
  genes; the Ensembl Metazoa AgamP4 build holds 12,614, so 1,271 genes other people measure are missing from
  our universe. **The Ensembl build is the wrong source for this species** -- VectorBase is where Anopheles
  annotation is maintained, and its AgamP4 gene set is the larger one. Rebuilding `Ag` from VectorBase is the
  fix; until then the space is listed as it is and its coverage is recorded here.

  `Fc` passes on 94 resolvable genes, which is a weak sample and is marked as such; `As` cannot be checked at
  all until either a reviewed proteome or someone else's gene list exists for it.

  **Diagnosed, 2026-11-09: the Anopheles gap is annotation, not a stale release, and VectorBase is behind the
  login.** Ensembl Metazoa release 63 was fetched and holds 12,627 genes against release 60's 12,614 -- thirteen
  more, so the release is not the cause. The difference is between the two annotations themselves: VectorBase's
  AgamP4 carries **1,271 AGAP genes Ensembl does not** (ordinary ids: AGAP000729, AGAP001105, AGAP001129 ...),
  and Ensembl carries 40 VectorBase does not. **VectorBase is a VEuPathDB site and its download paths 404
  without a login** (`https://vectorbase.org/common/downloads/...`), so rebuilding `Ag` from the right source
  is the first concrete, named thing the VEuPathDB API key (instruction 38) would unblock -- not a general
  request for access, but one file for one failing space.

## WP6 started: Cryptosporidium parvum (2026-11-08)

`parascape/spaces/crypto.py` (`CRYPTO_SPACE`, code `Cp`, CryptoDB `cgd` identifiers, organisms.CRYPTOSPORIDIUM) and
`scripts/build_cryptosporidium_space_notebook.py` -> `notebooks/cryptosporidium_space_2026_11_08.ipynb`,
`results/cryptosporidium_space_2026_11_08/` (gene table, column metadata, benchmark rows). From the sources already on
disk, each hashed: 3,923 Iowa II genes (18 non-`cgd` rows listed and left out) with product, OrthoMCL group, length;
measured hyperLOPIT location (Guerin 2023), measured stage (Tandel 2019 single-cell), InterPro, DeepTMHMM.

**Benchmark on Cryptosporidium truth:** the hyperLOPIT location measured on a gene's Toxoplasma orthologs, mapped to
Cryptosporidium locations on training genes only (folds over whole OrthoMCL groups), predicts the measured
Cryptosporidium location of 626 of 1,086 labelled genes **56.7 % right, against 21.9 % for the commonest class**
(20 location classes).

Not distributed: CryptoDB tables and the Guerin supplement carry licence code V until their redistribution terms
are verified; the table is local. P. berghei follows; P. vivax remains.

## WP6: Plasmodium berghei ANKA (2026-11-08)

`parascape/spaces/berghei.py` (`BERGHEI_SPACE`, code `Pb`, PBANKA identifiers) and
`scripts/build_berghei_space_notebook.py` -> `notebooks/berghei_space_2026_11_08.ipynb`,
`results/berghei_space_2026_11_08/`. 4,912 PBANKA genes from the PlasmoDB-68 annotated proteome on disk, with the
PlasmoGEM blood-stage knockout screen in its native home (Bushell et al. 2017, PMC5509546 Table S1: phenotype,
relative growth rate and interval, the named P. falciparum ortholog).

* **WP6 acceptance holds:** the P. falciparum `pb_transferred_phenotype` equals the projection of this space through
  the named orthologs for all 2,395 Pf genes it reaches (0 differ); 40 further Pf genes carry a transfer not reached
  from these named orthologs.
* **Benchmark on P. berghei truth:** the P. falciparum piggyBac mutagenesis index of the named ortholog separates
  essential from dispensable P. berghei knockouts at AUC 0.789 (1,964 genes, 1,153 essential).

Local only until the PlasmoDB annotation terms are verified (code V); PlasmoGEM supplement is PMC open access.

## WP4 started 2026-11-09: the rat and cat spaces -- and a wrong assumption corrected

**These were recorded as blocked behind the VEuPathDB login for weeks, and they never were.** A rat and a cat
are vertebrates; their gene builds are on Ensembl and nothing about them needs a parasite database. The
assumption was carried from the VEuPathDB note and applied to species nobody had checked. Checked on
2026-11-09: rat (mRatBN7.2), cat (Felis_catus_9.0) and Anopheles gambiae (AgamP4, Ensembl Metazoa) are all
public, and so is **P. vivax** (PvP01 v2, Ensembl Protists) -- so WP5 and the Pv part of WP6 are not blocked
either.

`parascape/spaces/rat.py` (`RAT_SPACE`, Rn) and `parascape/spaces/cat.py` (`CAT_SPACE`, Fc);
`scripts/build_host_gene_spaces_notebook.py` -> `notebooks/host_gene_spaces_2026_11_09.ipynb`,
`results/host_gene_spaces_2026_11_09/`. Builds hashed in the dataset taxonomy with their URLs.

| Space | Genes | Protein coding | Named |
|---|---:|---:|---:|
| Rn | 23,141 | 23,098 | 21,346 |
| Fc | 19,588 | 19,564 | 16,010 |

Both pass `spaces.validate_nodes` (unique, non-missing, canonical ids) and both are listed by
`space_catalogue` as **unavailable**, which is the honest state because no pack is built.

**Both now hold measurements as well**, the same two the mouse space carries, so the three host spaces answer
the same question the same way (Bgee expression calls, `datasets/host/bgee_expression`, hashed):

| Space | Genes tested by Bgee | Expressed somewhere | Tissues |
|---|---:|---:|---:|
| Rn | 18,134 of 23,141 | 17,802 | 22 |
| Fc | 18,862 of 19,588 | 18,596 | 12 |

`bgee_top_score` is the best expression score over the tissues a gene is called present in and
`bgee_tissues_present` how many those are. A gene Bgee tested and found absent everywhere keeps a measured
**zero**; a gene Bgee never tested stays **missing**. That distinction is the whole point: absent is a result,
untested is not.

Cat coverage turned out better than expected -- 96 % of its genes are tested, against 78 % for rat -- so the
earlier note that Bgee's cat coverage is thin was wrong.

**The Bgee source is NOT in `datasets.REGISTRY`, on purpose.** That registry describes what the application
ships, and it refuses an organism whose table is not shipped (`unknown dataset organism: Rn`). Rat and cat
tables live in `results/` until a pack is built, so the source is recorded where it belongs for now: hashed in
the notebook manifest, with its URL beside the files. It joins the registry when the packs do.

The cat matters out of proportion to its data: it is the only species where Toxoplasma completes its sexual
cycle, so the enteroepithelial stages exist in no other host, and several slots naming feline tissue have never
had a species to hang on.

What each needs next is written in `what_each_space_needs.tsv`: expression per tissue first (Bgee has rat; its
cat coverage is thin), then human-ortholog links for the rat, and the enteroepithelial series for the cat.

## 2026-11-10: NCBI is not a way round the VectorBase login, and here is why

The Anopheles blocker was re-tested rather than taken from the note. `vectorbase.org/common/downloads/...`
still answers **404** without a login, so that half stands.

The obvious workaround is NCBI, which carries the same species and keeps the VectorBase accessions
(`locus_tag=AgaP_AGAP000002`, `Dbxref=VectorBase:AGAP000002`). It is **not** a substitute, and the numbers say
why. `GCF_000005575.2` is **AgamP3**, while the space is built on **AgamP4**:

| | genes |
|---|---:|
| NCBI AgamP3, AGAP accessions | 13,233 |
| our Ag space, Ensembl Metazoa AgamP4 | 12,614 |
| in NCBI and not in ours | 2,086 |
| in ours and not in NCBI | 1,467 |

The difference runs **both ways**. These are two annotations of two assembly versions, not an old copy and a new
one, so merging them would produce a gene universe matching neither source -- precisely the silent
mis-assignment the strain-suffix rule exists to prevent (decision 2b). Taking the union would add 2,086 genes
whose coordinates belong to a different assembly.

So the conclusion stands, now with evidence rather than inheritance: **rebuilding `Ag` correctly needs the
AgamP4 annotation from VectorBase, and that needs the VEuPathDB key.** One file, for one failing space. NCBI is
recorded here as tried and rejected, with the counts, so it is not re-tried as a shortcut.

## 2026-11-10: what WP9's last piece actually needs, checked rather than estimated

"A host space as a full map window is not built" reads like a UI task with nothing blocking it. It is not, and
the difference matters for scheduling it.

**A map is a layout, and no host layout exists.** `hs_host_proteins.parquet` (27 columns) and
`mm_host_proteins.parquet` (10) carry no coordinate column of any kind, and the only layout files shipped are
`graph.npz` and `pf_graph.npz` -- Toxoplasma and Plasmodium. Opening a host space as a map therefore needs an
embedding computed over the host table first: 58,988 human genes for `Hs`, 18,428 for `Mm`.

So the remaining work is in three parts, and only the third is the window:

1. **Decide what makes two host genes similar**, which is not inherited from the parasite side. The parasite
   map is built from measurement blocks the harness knows; the human table is GTEx tissue medians plus HPA
   annotations plus ortholog-derived fields, and choosing among them IS the claim the map would make.
2. **Build and freeze the layout**, with the node-file SHA-256 recorded the way `pf_graph.npz` does, so the
   map cannot drift from the table it was built from -- the tie that broke twice this week.
3. **Open it in the window**, which is the module work this item has been describing.

Step 2 is real machine time on a shared box (a 58,988-row embedding), so this piece is gated on machine time in
the same way items 69, 80 and 81 are, not merely on effort. Recorded here so it is not picked up as a quick UI
job and then found to need an afternoon of compute.

## 2026-11-13: WP9 step 1 and step 2 -- host layouts as two selectable lenses

Einar chose **"both, as selectable lenses"** for what makes two host genes similar. Built by
`scripts/build_host_layouts.py` into `results/host_layouts_2026_11_13/`.

* **`tissue`** -- similar when expressed in the same tissues and cell types.
* **`parasite_facing`** -- similar in how they respond to infection and in what the parasite does to
  them: vacuole enrichment, the fibroblast and hepatocyte response to infection, rhoptry discharge,
  the erythrocyte surface.

Every host gene keeps a position in both: a gene no lens measures is NOT dropped, because dropping
it would make a map of what has been studied look like a map of the proteome.

| lens | genes | features | with at least one value |
|---|---:|---:|---|
| Hs tissue | 21,020 | 11 | 19,314 (91.9 per cent) |
| Hs parasite-facing | 21,020 | 9 | 20,492 (97.5 per cent) |
| Mm tissue | 16,473 | 6 | 16,417 (99.7 per cent) |
| **Mm parasite-facing** | 16,473 | 1 | **147 (0.9 per cent) -- REFUSED** |

**Two corrections to what this item previously recorded.**

**1. The compute estimate was wrong, and the row count behind it was wrong.** This file said 58,988
human genes and called the embedding "real machine time... an afternoon of compute", which is why
the piece was scheduled alongside items 69, 80 and 81. The table holds **21,020 human and 16,473
mouse** genes, and all four lenses build in about two minutes on eight cores. It was never gated on
machine time.

**2. The mouse parasite-facing lens is refused, not built.** It has one column and 147 of 16,473
genes carry a value. Built anyway -- and it was, before the guard -- it produced **16,473 distinct
positions**, because UMAP spreads points apart even when their inputs are identical: 99 per cent of
that map was the repulsion term acting on median-imputed copies of one row. It looks like a map and
is noise, which is worse than no map because nobody can see which it is. `MINIMUM_MEASURED = 0.10`
now refuses any lens below a tenth, and records the refusal in the summary rather than omitting it.

The thinness is itself the finding: **the mouse table carries one parasite-facing measurement**, the
macrophage surface. A second would make the lens real.

**Still to do (step 3):** opening these in the window, which is the module work. Each layout records
its features, its recipe and the SHA-256 of the node file it was built from, so the window can
refuse a layout that does not match the table it is shown beside.

**The host layouts are NOT bit-reproducible, and the freeze has to mean the file rather than the
recipe.** Two rebuilds of the same lens from the same table differ by up to 97 units, and a test
that ran both inside ONE process differs too -- so it is not the environment between runs, it is
umap-learn here, which does not honour `random_state` as a guarantee. The Toxoplasma layout DID
reproduce to the bit on 2026-11-13, so whatever makes that one stable does not carry over to this
path, and I have not found what it is.

What this costs: the tie recorded in each layout (`node_file_sha256` plus `coordinates_sha256`) still
says exactly which table a layout belongs to and detects any later drift, which is what WP9 step 2
asked for. What it cannot promise is that rebuilding reproduces the same coordinates. So the shipped
file is the artefact, not a recipe that regenerates it -- and anyone rebuilding gets a different,
equally valid arrangement of the same data. Worth deciding before step 3 whether that is acceptable
for a map people will compare between sessions.

## 2026-11-13: WP9 step 3 is blocked deeper than the layout, and the layout was never the blocker

Checked before building the window, and it changes what step 3 needs.

**The registered host spaces have no node table shipped.** `HUMAN_SPACE` and `MOUSE_SPACE` declare
`hs_gene_nodes.parquet` and `mm_gene_nodes.parquet` with gene patterns `ENSG\d{11}` and
`ENSMUSG\d{11}`. **Neither file exists in the data directory.** So `organisms.codes(available=True)`
cannot return them and nothing can open them, with or without a layout.

**The tables that DO exist are a different object.** `hs_host_proteins.parquet` (21,020 rows) and
`mm_host_proteins.parquet` (16,473) -- the two this item named when it scoped the work -- are keyed
on **UniProt accessions** (`A0A024R1R8`, `A0A087WPF7`), not on Ensembl gene ids. They are the host
PROTEIN tables that back the host bridges, not the host GENE spaces. `validate_nodes` would refuse
either of them for its space, correctly, because the identifiers do not match the declared pattern.

So the layouts built above are coherent and useful -- they lay out the host protein tables, which is
what this item pointed at -- but they **cannot be served as the Hs or Mm space** without one of:

1. **building the host gene tables** the spaces declare, which is the WP3 acquisition work this item
   records as "done to its acceptance" for a different meaning of done; or
2. **declaring the protein tables as their own spaces**, with a UniProt pattern, which is a design
   decision about whether ParaScape has a host gene space, a host protein space, or both; or
3. **mapping the protein layout onto gene ids** through the existing bridges, which inherits every
   ambiguity of a protein-to-gene mapping.

**This supersedes the 2026-11-10 note above.** That note said the missing piece was a layout and
that the work was gated on machine time. The layout takes two minutes and is now built; the actual
gap is that the host spaces point at tables nobody has shipped, and choosing which of the three
routes above to take is Einar's call rather than a detail to settle while building a window.

## 2026-11-13: Einar decided "both" -- host protein spaces now, host gene spaces after

Einar answered the question above on 2026-11-13: **both**. Route 2 now (declare the host protein tables as
spaces of their own), route 1 afterwards (the Ensembl gene tables, so `HUMAN_SPACE` and `MOUSE_SPACE` work as
declared). The two stay separate spaces: a protein is not a gene, and mapping one onto the other inherits every
ambiguity between them (one gene with several reviewed proteins, one protein with several genes).

### Route 2, built: the human and mouse host protein spaces

`parascape/spaces/host_proteins.py` declares `HUMAN_PROTEIN_SPACE` (code `HsP`, 21,020 proteins) and
`MOUSE_PROTEIN_SPACE` (`MmP`, 16,473), keyed on UniProt accessions (the documented UniProtKB pattern; every row
of both tables matches it). Each opens in a map window through its lenses, the WP9 layouts:

| space | lens | menu wording | state |
|---|---|---|---|
| HsP | tissue | Laid out by tissue | opens |
| HsP | parasite_facing | Laid out by what the parasite does to them | opens |
| MmP | tissue | Laid out by tissue | opens |
| MmP | parasite_facing | -- | **refused**, shown disabled with its reason (147 of 16,473 proteins, 0.9 per cent) |

* **File ▸ Species** gains one submenu per protein space, one entry per lens. A refused lens is listed and
  disabled with its reason rather than hidden. The window title names the lens.
* **A layout is served only beside the exact table it was built from**: `host_proteins.load` checks the SHA-256
  the layout records against the table on disk, the protein order, and that every coordinate is finite. A
  layout built from a different table is refused (tested with a changed table).
* **The layouts ship with the next data release.** The three files were copied from
  `results/host_layouts_2026_11_13/` into `parascape/data/` (the release builder takes every `.npz` there).
  Until that release exists the downloaded data directory has the tables and not the layouts, so
  `layout_path` looks in the data directory first and then in the source tree. That fallback is safe only
  because of the SHA-256 check above. `data_manifest.json` is regenerated at release time, not here.
* **Picking a protein shows a protein panel**: its name, its UniProt link and every column, with "unknown" for
  a missing value and the columns that placed this map marked. Before this, the parasite panel called every
  host protein "effectively uncharacterised", because it reads literature columns a host table does not have.
* **Declared, not registered in `organisms.SPACES`**, for two reasons: registering changes what every slot,
  strategy and calibration caller sees (the WP2 lesson), and a UniProt accession does not say which species
  it belongs to, so `organisms.detect` could not tell human from mouse. `organisms.by_species` finds them, so
  record links work. A protein window is never reopened at start-up.
* No edge layers: none has been built for a host table, and drawing the parasite ones would be meaningless.

### Route 1: further along than this file said

**Correction to the note above.** It said the gene tables `hs_gene_nodes.parquet` and `mm_gene_nodes.parquet`
"do not exist". They exist, in the host gene packs built on 2026-11-07 and tracked in git:
`results/human_pack_2026_11_07/parascape_Hs_2026.11.07.zip` (58,988 Ensembl genes, with
`hs_gene_graph.npz` holding coordinates for every gene) and
`results/mouse_gene_space_2026_11_07/parascape_Mm_2026.11.07.zip` (18,428 genes, `mm_gene_graph.npz`). They are
not in the data directory because these spaces are distributed as packs: `paths.space_file` resolves them from
the installed pack, never from the data directory. So "neither file exists in the data directory" was true and
was the wrong place to look.

What was built today: an installed host gene pack now **opens as a map window** ("Mus musculus genes",
"Homo sapiens genes"), from the same Species menu entry that installs it. The node table is checked with
`validate_nodes` and the layout with `validate_graph` before anything is drawn, and picking a gene shows its
Ensembl record and every value. Rat and cat have no layout, so they stay install-only.

What route 1 still needs, none of it a large download:

1. **A published download for the packs.** They are installable only from a checkout (`results/`), so a user
   of the released program sees them as unavailable. Attaching the two zips (6.3 MB and 2.3 MB) to a release,
   with their SHA256 files, and pointing `space_catalogue` at that URL is the step. Einar to choose where they
   are published.
2. **What the gene layouts claim.** The 2026-11-07 gene coordinates predate the lens decision and their recipe
   is not recorded beside them. Rebuilding them as the same two lenses (tissue; parasite-facing) needs
   measurements the gene tables do not yet carry. The raw files are already here, under
   `STARPLAST_DATA/spaces/Hs` and `spaces/Mm`: GTEx v10, HPA v25.1 tissue and single-cell, Bgee 15.2,
   PaxDb, OpenCell, and the infection series (human GSE145836, GSE241570/1 and others; mouse GSE103113,
   GSE119855 and others), plus BioGRID, IntAct, ORCS and NCBI orthologs in `spaces/shared`. Nothing more
   needs downloading for a first version. Each column must keep its licence code: the GEO-derived ones stay
   out of a distributed pack until Einar confirms GEO deposits may be redistributed (open since 2026-11-09).
3. **The link between the two kinds of space.** `hs_protein_gene_links.parquet` and
   `mm_protein_gene_links.parquet` are already inside the packs; carrying a selection from a protein window
   to its genes (and back) is the Phase 2 "linked spaces" work, and every link is many-to-many.

### Still open in the window, for both kinds of host space

* The Map module "What makes genes similar" box is empty for a host map, because it lists parasite
  measurement blocks. It should show the lens features. "Build this map" would rebuild from parasite blocks;
  until that is fixed the shipped layout is the map to read.
* The default colour is the first category column (`rbc_surface_found_senegal` for human proteins), which is
  an odd first view. A host table has few categories; binned expression would be a better default.
* The layouts are not bit-reproducible (noted above); the shipped files are the artefact.

Tests: `tests/test_host_protein_spaces.py` (13): both tables validate as UniProt spaces with every value kept;
gene ids and protein ids never pass for each other; every built lens loads beside its own table; a layout from
a changed table is refused; the refused mouse lens stays refused; protein spaces stay out of the registry and
`detect`; a protein window opens through each lens and shows the protein panel; an installed mouse gene pack
opens as its own map. Rendered offscreen in the real Window (Map module) for human proteins by both lenses and
for mouse genes, and looked at.

**The full suite found one test this change had to rewrite.** The focused run made while building this did not
cover it: `tests/test_module_shell.py::test_the_species_menu_lists_host_spaces_with_their_state` asserted that an
installed host gene pack leaves a *disabled* menu entry, which was the old behaviour -- an installed pack was a
dead entry. It now opens its map, so the entry stays enabled, and the test asserts that instead, keeping what
instruction 64.23 asked for: the state is spelled out and never reads as calibrated. A pack with no layout (rat,
cat) is still disabled once installed. Full suite on the rebased branch: 7,569 passed, 42 skipped, 1 h 6 min.

## 2026-11-13: Einar decided the two pack questions

* **GEO-derived columns may be redistributed, with citation.** "Yes, with citation." Every GEO-derived column
  ships in its pack with its GEO accession and the paper behind it named in `NOTICE.txt`. This closes the
  question open since 2026-11-09. It is not the one-line change the earlier note promised: `packs.build`
  accepts only `O` and `SA`, so GEO columns need their own licence class (public deposit, cited) that the
  builder accepts only when an accession and a citation are both recorded. Moved into the hosting work below.
* **All ParaScape data moves to Hugging Face, in the account of Einar.** "all of the data that parascape uses
  should be hosted on huggingface in my account. otherwise the space limits on github will fill quickly."
  The gene and protein packs therefore go to a Hugging Face dataset, not to GitHub release assets, and the
  downloader reads from there. The 0.0.11.x data release already being published to GitHub on this date
  finishes as planned; the move follows it.
