# Instructions: what is done and what is left

**New here? Read [START_HERE.md](START_HERE.md) first.**

One file per task. This index is the status table; the files carry the reasoning and the
"done when" condition.

## Done

| # | Task |
|---|---|
| 68 | [Broaden Discoveries beyond localization](done/68_discoveries_function_coverage.md) — **done 2026-11-08**: functional targets, benchmarks and 79,289 precomputed functional claims (EC major classes, graded; none fills) · built; Einar's review pending -- reopen on his feedback |
| 74 | [Navigator and home](done/74_navigator_and_home.md) — **done 2026-11-08**: phases delivered by 84-97; Infer step 3 lets you answer with any harness method, graded on held-out genes · built; Einar's review pending -- reopen on his feedback |
| 70 | [Inference scoreboard and navigator](done/70_inference_scoreboard_and_navigator.md) — **done 2026-11-08**: every action delivered; Tools > Question scoreboard shows all 198 questions scored on held-out genes · built; Einar's review pending -- reopen on his feedback |
| 90 | [Release 0.0.10.0 under the new licence](done/90_release_under_the_new_licence.md) — **done 2026-11-08**: approved by Einar; 0.0.10.0 on PyPI (Python 3.12+) and GitHub with the corrected data asset |
| 64 | [Explore published organism/host evidence and inspect precomputed inferences with ground-truth scorecards](done/64_information_space_and_inference_atlas.md) — **done 2026-11-08: 41/41 sub-items, released as ParaScape 0.0.10.0** · built; Einar's review pending -- reopen on his feedback |
| 84 | [The Navigator, exactly as spaCR's](done/84_navigator_exactly_like_spacr.md) — **built, waiting on Einar to compare** · a View switch, search on Return, spaCR's 224 px square cards with icon and light centred name · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 85 | [Retire the old ParaScape home screen](done/85_retire_the_old_home_screen.md) — **built, waiting on Einar** · the landing is the only home; nothing the old node-card dock did is lost · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 86 | [Every module in spaCR's module shell](done/86_spacr_module_shell.md) — **built, waiting on Einar: every available tile is a module in the shell** · settings left, console, AI, actions, live preview, hyperparameter search, collapsible · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 87 | [The network as a navigable universe](done/87_network_universe.md) — **built (first version), waiting on Einar** · whole network, gene region zoomed, zoomable, edge definitions · depth of field from the centre, click to travel, edge and layout controls on the right · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 91 | [Navigator: one click shows what comes next, a second click opens](done/91_navigator_module_network.md) — **built, waiting on Einar** · downstream module graph, fewer results, reset on return · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 92 | [The Map module preview in 3D](done/92_map_preview_in_3d.md) — **built, waiting on Einar** · a rotatable 3D preview, painted in software · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 93 | [Every module in the same spaCR theme](done/93_every_module_one_theme.md) — **built, waiting on Einar** · spaCR fields, containers at spaCR pane opacity (60 %), the gene card as a shell module · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 94 | [The AI box with spaCR's choices, and a free ParaScape AI](done/94_ai_box_like_spacr.md) — **built** · all four providers listed; the free one needs a worker Einar agrees to deploy · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 96 | [The ParaScape module mould, as a skill](done/96_parascape_module_skill.md) — **built: Map rebuilt, skill written** · get Map right with Einar, then write it down as the `parascape-module` skill and build the rest with it · built; Einar's review pending -- reopen on his feedback (2026-11-07) |
| 83 | [Move ParaScape from PyQt6 to PySide6](done/83_pyside6_migration.md) — done 2026-11-05: 6,550 passed; five latent defects found, two that would have shipped |
| 95 | [spaCR's dock](done/95_spacr_dock.md) — **built** · the left column of categories and module rows, on Home and in every module |
| 97 | [Errors become ParaScape GitHub issues, and a routine fixes them](done/97_errors_become_github_issues.md) — **built** · spaCR's relay + consent dialog; GitHub App key set by Einar; scheduled triage |
| 62 | [The hold-out track record](done/62_holdout_track_record.md) — **done 2026-11-07** · pooled record agrees with the scorecard (Spearman 0.988 over 67 cells); questions page shows it; map neighbours metric renamed |
| 60 | [What would make ParaScape most useful](done/60_usefulness_review.md) — **done 2026-11-07** · all ten changes built: methods text + BibTeX, run provenance, vector figures, what is unusual about a gene, per-label calibration, missingness/dominance warnings, imported screens reach the strategies, gallery by structure, tutorial videos, a regression guard for the empty-Strategies report |
| 77 | [The borrowed profile-inference strategies](done/77_borrowed_inference_strategies.md) — **done 2026-11-06** · six scored on frozen partitions through the harness, three untested with the reason; none beats boosted trees, positive-unlabelled bagging matches them; nothing may fill |
| 76 | [Dataset quality, including citations normalised by age](done/76_dataset_quality_and_citations.md) — **done 2026-11-06** · sample quality predicts replicate agreement (0.47, 41 studies); citations do not (−0.19), so they weight nothing |
| 79 | [Within-study control review](done/79_within_study_control_review.md) — **done 2026-11-06** · reviewed references applied; 178 contrasts, 0 indistinguishable; wild-type alkaline agrees across studies at 0.439 (was 0.236) |
| 89 | [GSE60305: colour-space SOLiD reads, not a kallisto fault](done/89_gse60305_kallisto.md) — 0 of 47 M reads could map; pipeline now records colour-space runs and treats 0 mapped as failure; processed readCount table queued |
| 88 | [Proteomics downloads: parallel files, stall restarts, the faster transport](done/88_pride_download_throughput.md) — FTP 0.14 vs HTTPS 0.71 MB/s today; counter fixed |
| 66 | [Recover missing source files and audit host candidate versions](done/66_source_recovery_and_host_candidates.md) — 4/4 complete: source recovery, RBC review, journal equivalence and verified host ambiguity correction |
| 78 | [A user's own layer: keeping an inference without contaminating the map](done/78_kept_inference_layer.md) — **done** · kept claims appear in Predicted, never in Facts, never as training data |
| 68.01 | [Browse functions, class members and inference/test coverage in Discoveries](done/68_01_discoveries_label_browser.md) — 41/30 labels, exact functional membership replay, 413 final checks |
| 01 | Nineteen gene codes resolved |
| 02 | Six expression datasets integrated |
| 03 | GUI bugs |
| 04 | `projectionMatrix` compatibility — clicking genes |
| 05 | Visual quality |
| 06 | Studies with no supplementary data |
| 07 | Verify shipped columns against GEO |
| 08 | Standalone data and auto-download |
| 09 | Cell-cycle dataset |
| 10 | Search more targets |
| 11 | README and API reference |
| 12 | Test every aspect of the library |
| 13 | Jobs, stopping, and freeing memory |
| 14 | Tooltips, docstrings, and the scoring explainer |
| 15 | The UMAP gallery — grid, scroll, incremental |
| 20 | Validation: `refit`, candidates, orthogonal evidence, per-category scores in Inference |
| 19 | Logging — opt-in, per-level console control |
| 18 | Navigate / Select modes, 2D and 3D gating |
| 16 | Automated walk — UMAP, clustering, per-category scoring |
| 17 | Annotation store |
| 23 | Left panel becomes "color by" |
| 22 | Cell diagram under the compartment list |
| 25 | American spelling in user-facing text |
| 24 | Forty black-and-white logo drafts |
| 21 | Scoring objectives in the interface |
| 26 | Save and load results, per tab and all at once |
| 27 | Import data, with the preprocessing offered |
| 25b | The identifier rename — job 2 of the spelling task |
| 28 | Every module at 100% coverage |
| 29 | Three more leaks in the circularity guard |
| 30 | Fold in datasets already present and expand missing biological axes |
| 31 | Replace file-shaped feature blocks with biological question slots |
| 32 | GPU paths for k-means, DBSCAN and t-SNE |
| 33 | Crossed-factor findings and fragmentation diagnostic |
| 34 | Repair two metric columns and the CLI finding count |
| 35 | Simplify lighting, make point modes distinct, and add volumetric ray tracing |
| 36 | GPU PBR sphere points and stable GPU density-ray tracing |
| 37 | Continuous flashlight, material lab, and ray-rendering comparison |
| 48 | The orphan alarm never looks at the host table |
| 49 | A Strategies tab: 32 ways to infer something, each testing itself |
| 50 | The September 2026 data audit: a wrong citation, 21 deposits, two leaks closed |
| 51 | How good each strategy is: calibration over settings, targets and seeds |
| 52 | The second September 2026 data audit: six Plasmodium deposits, five refusals |
| 54 | [Scorecards, techniques and five advanced strategies](done/54_scorecards_and_advanced_models.md) |
| 55 | [Strategy cards: four standard bars, the test explained, a real failure and success](done/55_strategy_cards.md) |
| 56 | [Pregenerated map gallery and label x map scores](done/56_umap_gallery_label_scores.md) |
| 56 | [Star map: gene-gene links by source, with provenance, including your own runs](done/56_star_map.md) |
| 56 | [The parasite-density screen (Giuliano et al., Cell 2026): fitness at low and high density, and what crowding needs](done/56_lourido_density_screen.md) |
| 57 | ["Start here": a guided path from what a user has to the strategies worth running](done/57_start_here_tab.md) |
| 58 | [The blank panel: a tab switched away from comes back as bare background](done/58_blank_panel_bug.md) |
| 61 | [Tutorial videos: six screen captures of the real application](done/61_tutorial_videos.md) |
| 54 | [Many more pregenerated maps, searched for structure, and a navigable maps panel](done/54_many_umap_maps_and_the_score_table.md) — shipped in 0.50.0 |
| 58 | [Star map: stop the hover twitch, show the larger network, define the connections](done/58_star_map_steady_larger_defined.md) — merged and shipped in 0.50.0 |
| 64.01 | [Define organism-qualified entities, questions and provenance-preserving alias resolution](done/64_01_entity_query_contract.md) — 472 focused checks passed |
| 64.02 | [Inventory registered evidence and reconcile installed storage coverage](done/64_02_dataset_inventory.md) — 162 sources, 180 rows, 16 scoped refusals; 525 focused checks passed |
| 64.03 | [Trace each output to source identities, transforms, mapping losses and visible gaps](done/64_03_measurement_provenance.md) — 162 sources, 628 addresses; exact 15,437-value host pilot |
| 65.01 | [Versioned citation-rate, recency and comprehensiveness preferences](done/65_01_dataset_selection_policy.md) — 389 focused checks passed; publication audit follows |
| 65.02 | [Audit all registered publication identities and annual citation rates](done/65_02_publication_identity_audit.md) — 162 sources accounted for, 95 PMID/DOI identities resolved; remaining roles/gaps explicit |
| 65.03 | [Search every slot and retain biological comparison gaps](done/65_03_slot_literature_discovery.md) — 297 slots searched, 2,032 distinct publications retained, six diagnostic traps triaged; no promotions |

| 66.02 | [Review quantitative RBC groups, donors, incumbent quantities and admission gaps](done/66_02_rbc_candidate_review.md) — 2,653 groups, 2,061 unambiguous noncolliding protein projections; retain incumbents |
| 66.03 | [Verify rhoptry journal/preprint data equivalence and upgrade citation provenance](done/66_03_rhoptry_journal_version.md) — 20,010 genes × three fields exact; 39 ambiguous symbol mappings separately queued |

| 66.01 | [Reconcile all registered source locations and recover named processed inputs](done/66_01_source_recovery.md) — 162 sources, 69 bindings, 538 container candidates; 717 checks passed |

| 64.07 | [Declare strategy/query/output meanings and technique validation roles](done/64_07_strategy_capability_contracts.md) — 39 strategies, 40 techniques, 624 routes; 855 final checks passed |

| 64.08 | [Version inference, held-out rows, model metadata and scorecard identities](done/64_08_result_artifact_contract.md) — ten round trips across six tasks/five roles; 840 checks passed |

| 64.09 | [Reconcile all individual outcomes and preserve scorecard denominators](done/64_09_scorecard_aggregation.md) — 1,038,372 rows, 4,112 cohorts, zero arithmetic mismatches; 835 checks passed |

## Open

The current execution order is [instruction 64](open/64_information_space_and_inference_atlas.md):
the user's October 7 goals, 40 controlled action cards, a ground-truth contract for all 39 current
strategies, and the full progress table. Existing instructions below retain their evidence and
remaining work; instruction 64 integrates their priorities. Repost the complete progress table
whenever an action is completed, replacing its percentage with a green tick.

| # | Task |
|---|---|
| 67.01 | [Keep native work continuation active and log actual stops](done/67_01_work_continuation_watchdog.md) — timer active, 12 historical boundaries/4 restrictions logged; 403 checks passed |
| 66.04 | [Correct ambiguous host projections and preserve all original gene evidence](done/66_04_host_symbol_projections.md) — 39 projections withheld, 20,010 gene rows retained; 575 follow-up checks passed |
| 82 | [Fingerprints from a literature meta-analysis: gene-set signatures](open/82_literature_set_fingerprints.md) — **v2 (2026-11-08): abstract + full-text mining (86 rows, 896 papers); 6 rungs pass recovery, membrane sides still not separable; host-cytosol rung unbuildable from the literature** · the compartment ladder; rung 4 is the hard-negative set whose absence capped host-facing at 0.822 |
| 81 | [Re-quantify the earlier runs: trimming, equivalence classes and Gibbs samples](open/81_requantification_topup.md) — **ready; waits on downloads and on Einar installing the service** · queued behind the new downloads; 40 zero-parasite runs first · 2026-11-08: service installed on Einar go-ahead; waits for the 1,313 new runs |
| 80 | [Identify the pull-down baits so the background table becomes usable](open/80_pulldown_background_baits.md) — **22 baits (was 11); gate waits on more projects** · 11 baits over 34 experiments is too few; the frequency filter would delete 12 of 13 complex members |
| 75 | [Condition fingerprints: annotate every dataset with metadata nobody published](open/75_condition_fingerprints.md) — **v1: cell cycle and bradyzoite shipped, validated on held-out studies; 6 axes open** · phase advances in all 11 hourly steps; bradyzoite AUC 0.92-0.93; exposes a GSE223819 label swap
| 73 | [The question matrix: 100-1000 trait questions, each with an answer, a probability and an uncertainty](open/73_question_matrix.md) — architecture requested 2026-11-05; 286 of 411 numeric columns already cover >=5,000 genes, so most questions are carried rather than predicted · 2026-11-06: host-facing scored on frozen partitions (item 77); category rule vs the 21 map questions measured, 3 differ, decision for Einar · 2026-11-08: catalogue 158 admitted (Einar: keep all); 198 questions through the harness, median frozen AUC 0.758, 7 may fill; leak check holds |
| 72 | [The information map: gene × trait values, a matching probability matrix, calibrated gap filling at 95 % precision](open/72_information_map.md) — v0 observed-facts build started 2026-11-04 |
| 71 | [What the aggregated omics are for: co-expression → localization, perturbation atlas, host response atlas, Toxoplasma background frequency](open/71_uses_of_the_aggregated_omics.md) — planning (2026-11-04); starts on verified runs while downloads continue · 2026-11-06: product 1 through the harness: co-expression alone, median AUC 0.897 on 13 niches (0.934 with all evidence) · 2026-11-08: product 3 query and cross-lab recovery (deltaBFD1 host response reproduces weakly) |
| 69 | [Perturbation omics: every gene-knockout, drug and strain study, compared condition by condition](open/69_perturbation_omics_and_contrasts.md) — 711-study catalogue; Compare tab and limma-trend engine shipped on nightly; 2,191 RNA-seq runs and 203 PRIDE projects being reprocessed uniformly; proteomics, metabolomics, cross-study profiles and the expression pack remain |
| 65 | [Audit every dataset choice and retain the citation-rate/recency/comprehensiveness rule](open/65_dataset_selection_audit.md) — 3/4 bounded items complete; batch 5 scored 2026-11-08 (1 replacement proposed, nothing installed); origin-paper verification remains 65.04 |
| 63 | [Claims: generate knowledge with a measured certainty, then test it independently](open/63_claims_generate_and_verify.md) — steps 1+2 shipped in 0.54.0; October 7 literature audit complete (no Tg candidate passes, Pf candidate adds only two claims; no promotion). Next: orthology/experimental verifier coverage, frozen prospective tests, genes outside the tested range · 2026-11-06: first frozen prospective claims (747 genes the glucose screen missed, hashed) and a grader · 2026-11-08: orthology verifier -- independent for Pf location (548 confident held-out calls at 88 %, was 90), not for Tg |
| 38 | [Build a provenance-first pan-Apicomplexan dataset archive](open/38_pan_apicomplexan_dataset_archive.md) — machinery done 2026-08-18 (248 sources); re-running discovery is blocked by the VEuPathDB login since release 71 (re-checked 2026-11-06), gene builds via Ensembl Genomes |
| 41 | Fill as many slots as possible: Toxoplasma, Plasmodium, and host |
| 59 | [One hundred biological questions this software can answer](open/59_biological_questions.md): 100 questions with entry point, strategy, settings and a resolved literature anchor; Q01-Q50 executed and graded yields/thin/no · 2026-11-06: Q51-Q100 run (42 encoded, 40 ran, 8 skipped with the reason); Q91 is agreement with ExportPred, not export |
| 53 | [One space per organism](open/53_organism_spaces.md), then links between pathogen, vector and host. R0-R3, Pf display, pack framework and lossless calibration publishing implemented. Remaining: organism builders and published packs · Hs/Mm spaces · Space menu · Cp/Pb, Ag/As, Rn/Fc spaces · per-space leakage and calibration · cross-space bridges (version numbers to be assigned; 0.47-0.54 went to other work). Small follow-ups: retry active learning with conformal uncertainty; next Tg/host data audit; `spaces/As/abundance/PXD001647/` is junk to delete. See WORK_LOG.md for active progress and blockers. · 2026-11-08: WP1 + WP3 (Hs 99.3 % reviewed proteins resolved, Mm) done; WP6 started: Cryptosporidium space, Tg ortholog locations predict Cp location 56.7 % vs 21.9 % |

## The one ordering constraint that is not negotiable — now satisfied

**17 must not ship before 20**, and 20 landed on 2026-08-12. Annotation is easy to build and easy to
build wrongly. Without validation it gives the project a mechanism for manufacturing unvalidated
claims, which is the one thing it exists to prevent — and a candidate list looks identical whether it
is 90% right or 6%.

The constraint does not disappear now that 20 exists; it becomes a requirement on 17. A saved
annotation carries the validated precision and recall for its category, the cluster's composition,
and the enrichment over prevalence — the numbers the Validation tab now computes. Saving a row
without them would be the same failure by a shorter route.

## Standing constraints on all of it

Measurement, inference, absence and annotation never read as one another.

No candidate is shown without the number that says how much to believe it. On the full proteome one
row out of 7,710 (configuration x compartment) reaches F1 0.5 — `PM - peripheral 1`, at exactly
0.500. A cluster that looks pure is a lead to check, not a result.
