# Disk budget and cleanup list

Set by the PI on 2 Nov 2026: **50 GB total for the project's added storage** (raised from 20 GB for `~/scfm_models`), "to ensure we are not shortchanging ourselves". Track what can be deleted; **clean up when the paper is ready, not before.** Nothing below is to be deleted until then.

Counted against the budget: `~/scfm_models`, the scFM virtual environments, `data/yeast`, and anything new (DREAM5, fine-tuning, later jobs). Pre-existing project data (`data/perturbseq`, `data/beeline`, M3D, etc.) and `results/` are listed for information and are not counted. Machine: one 954 GB NVMe, 311 GB free on 2 Nov (66% used).

## Counted (2 Nov 2026, 21:50)

| Location | Size | Contents | Deletable when the paper is ready? |
|---|---|---|---|
| `~/scfm_models/step3` | 5.7 GB | frozen step-3 extraction (scPRINT arms A/B, R1, scGPT, controls, K562 step 3); SHA-256 manifest in `results/scfm_step3/` | No: the scored inputs; keep until the repository tag, then recomputable (about 1.5 GPU-hours) |
| `~/scfm_models/smoke_out` | 2.7 GB | step-2 smoke-test outputs | **Yes** (superseded by step 3) |
| `~/scfm_models/scprint` | 1.5 GB | scPRINT checkpoints (medium-v1.5 222 MB; large-v1, small-v1) | Re-downloadable; keep while fine-tuning runs |
| `~/scfm_models/step4_scoring`, `step4_k562` | 1.0 GB + 0.5 GB | derived scoring exports (BEELINE, K562) | **Yes** after the final numbers are archived (re-exportable by `analysis/scfm_score_export.py`) |
| `~/scfm_models/raw` | 0.9 GB | downloaded raw count matrices (GEO) | Re-downloadable; keep until K562/BEELINE inputs are archived |
| `~/scfm_models/scgpt_human` | 0.2 GB | scGPT whole-human checkpoint | Re-downloadable |
| `~/scfm_models/step2`, `lamin_*` | 0.2 GB | gene maps, leakage audit, lamin instances | Keep (small; the audit is cited) |
| `~/venvs/scfm` | 7.1 GB | scGPT environment (torch 2.6, CUDA) | **Yes**, rebuildable from the install notes in `docs/STATUS.md` |
| `~/venvs/scfm_alt` | < 0.01 GB | scprint 2.2.0 environment (links) | Keep |
| `data/yeast` | 2.1 GB | Kemmeren deleteome and series matrices (never scored; verdict 1) | **Yes** if the yeast test is dropped; keep if DREAM5/yeast follow-up goes ahead |

Total counted: about 22 GB of 50 GB.

## Not counted, for information

| Location | Size | Note |
|---|---|---|
| `results/` | 27 GB (gitignored) | largest: `clr_logpt_test` 7.5 GB; `baselines_k562ess_truthgenes_genie3` 2.5 GB; `grnboost2_superseded` 1.9 GB (**superseded GRNBoost2 runs, deletable**); `grnboost2_check` 1.8 GB; `beeline_baselines` and `baselines_k562ess_genie3` 1.6 GB each |
| `data/perturbseq` | 11 GB | K562 Perturb-seq inputs: needed |
| `data/beeline` | 0.75 GB | BEELINE inputs: needed |
| `~/venvs/clr-baselines`, `gb2-legacy-check*` | 0.6 + 0.5 + 0.5 GB | `gb2-legacy-check-py311-failed` is a failed environment: **deletable** |
| `docs/whitepapers/figures/pb_runs/` | 0.16 GB (untracked) | old PaperBanana runs; keep only the round-4 layout reference |

## New items (append below; one row each: path, size, source, deletable yes/no)

| `data/dream5` | 0.24 GB now (login-free expression archive, zip + extracted); expected about 0.5 GB with the Synapse training files and scoring code | DREAM5 network-inference data (login-free GeneNetWeaver-site archive; Synapse `syn2787209` for TF lists, gene IDs, gold standards). Counted. Added 3 Nov 2026 (`docs/plans/dream5_design.md`) | **Yes**, re-downloadable (hashes in `docs/dream5_selections/download_log.tsv`); keep until the DREAM5 scoring is archived |
| `results/dream5` (planned) | expected 2-4 GB, ceiling 5 GB for all of DREAM5 | method score caches (PIDC edge lists about 1.5 GB are the largest), refit-replicate tables. Counted for DREAM5 even though under `results/` | **Yes** for PIDC edge lists and score caches after the final numbers are archived; keep the committed tables |
| `~/venvs/dream5` | 0.07 GB | `synapseclient` 4.14 only (downloads; holds no credentials) | **Yes**, rebuildable |

| `results/dream5/net{1,3,4}/*.scores` (step 5, 4 Nov 2026) | measured at the time of writing: 1.0 GB (expression exports 13 + 65 + 80 MB, symmetric score caches 11 / 82 / 142 MB per method on net1 / net3 / net4, TF-row caches 2-16 MB, GRNBoost2 edge lists 9-70 MB); expected final about 2.5 GB after the PIDC gzip job, but up to about 5.5 GB while the PIDC edge lists are plain tsv (Julia writes both orientations: 0.08 / 1.3 / 2.1 GB for nets 1 / 3 / 4; job 2093_d5_pidc_gzip compresses them once the score caches exist; measured 80 MB on net1), then 11 / 82 / 142 MB caches and the all-genes tree runs (GENIE3 full matrices 22 / 163 / 283 MB, GRNBoost2 edge lists) | method score caches (upper-triangle doubles, SHA-256 sidecars; manifests `results/dream5/SCORE_MANIFEST.tsv` and `results/dream5/net4/refit/REFIT_MANIFEST.tsv` are the committed record). Counted for DREAM5 | **Yes** after the final numbers are archived: PIDC edge lists (gzip them after conversion), expression exports, GRNBoost2 edge lists; keep the committed manifests |
| `results/dream5/net4/refit` | < 0.01 GB | per-replicate SHA-256 sidecars of the 100 experiment-jackknife refits (hashes only: storing 100 x 4 gene-gene matrices would take 57 GB; the truth stage recomputes the refits and checks them against these hashes) | Keep (small) |

DREAM5 projected total: well under 5 GB, so the counted total stays at about 22 to 27 GB of 50 GB.

## scFM fine-tuning

Step 3.10 (E17, E18), 3 Nov 2026. The fine-tuned state dicts (14 scPRINT at 74 MB, 4 scGPT at 207 MB, about 1.9 GB written) were deleted by the extraction jobs 1041 and 1052 after the truth-free QC gate; `~/scfm_models` is 13 GB after the deletion (15 GB while the state dicts existed). Counted total: the project stays well under 50 GB.

| Path | Size | Source | Deletable |
|---|---|---|---|
| `~/scfm_models/step5_ft/scprint_attn` | 1.8 GB | fine-tuned scPRINT attention (`.attn.npy`, `.nets.npz`, json), 14 cases; frozen by `results/scfm_ft/sha256_manifest.csv` | No while the arm is cited (the scored input); recomputable by `analysis/scfm_finetune_scprint.py` plus `scfm_extract_scprint.py --load-state-dir` (about 2 GPU-hours) |
| `~/scfm_models/step5_ft/scgpt_attn` | 0.25 GB | fine-tuned scGPT attention, 4 cases | As above (about 1 GPU-hour) |
| `~/scfm_models/step5_ft/scoring` | 0.23 GB | float32 exports of the fine-tuned networks for the R scoring (`scfm_score_export_ft.py`) | **Yes**, re-exportable from the manifest-checked attention files |
| `~/scfm_models/step5_ft/scprint_ckpt`, `scgpt_ckpt` | < 0.001 GB | per-case training json only (loss curves, SHA-256 of the deleted state dicts); the `.pt` files are gone | Keep (tiny, record of the training) |
| `results/scfm_ft`, `results/scfm_score_ft` | < 0.001 GB | manifest, QC tables, metric table, contrasts (force-added to git) | No |

| `data/dream5/truth_extracted`, `truth_extracted_tools`, `truth_extracted_universe` | about 0.1 GB incl. converted null densities (17 MB) | Truth stage 6 Nov 2026: extracted official scoring code, gold standards, Python port and converted densities (gitignored) | **Yes** after the DREAM5 scoring is archived; the zip is re-sealed and re-extractable |
| `results/dream5/_pidc_par_test` | 0.08 GB | scratch test of parallel PIDC on net1 | **Yes**, any time |
| `data/dream5/published_extracted` | about 0.26 GB (predictions 0.25 GB, results xls and CSV, symlink dir `flat_pred`, tidy table) | Scorer gate 10 Nov 2026: extraction of `Network_predictions.zip` and `DREAM5_NetworkInference_MethodScoringResults.xls` (PI's downloads, `data/dream5/synapse/published`), gitignored | **Yes** after the gate report is archived; re-extractable from the zip, `flat_pred` and the tidy CSV rebuilt by `analysis/dream5_published_table.py` |

| `results/dream5/truth_net{1,3,4}` (11 Nov 2026) | < 0.5 MB each (CSV tables, universe lists, run_info) | Truth stage primary scoring runs of our methods (step 7 part 3). Force-added to git | No (the committed record of the truth stage) |

| `results/robust_z/{m3d,k562}` (11 Nov 2026) | 36 MB (tables and CSVs about 2 MB; `cache/` rds with replicate-level results 34 MB, gitignored and not committed) | Exploratory robust row-z run (`analysis/robust_z/`). Tables, CSVs and logs force-added to git | `cache/` **yes** (recomputable by the two scheduler jobs, about 5.5 h); the rest no |

| `results/granularity_log` (11 Nov 2026) | 0.9 MB (CSV curves, bootstrap draws, two PNGs, logs; force-added to git) | Log-scale depth re-test (`analysis/granularity_log.R`, `depth_log_analysis.py`) | No (the committed record of the re-test) |

| `results/explore_ideas47/` (12 Nov 2026) | 2.3 MB measured (rds readout caches of ideas 4 and 7: 3 scopes x 100 replicates x 10 truths per readout; no matrices are kept) | exploratory ideas 4 and 7 (`docs/plans/explore_ideas47_plan.md`); caches not committed, small CSVs copied to `docs/explore_ideas47_tables/` | **Yes** (recomputable in about 3 h; keep the committed tables) |

| `results/idea5` (12 Nov 2026) | 0.02 GB kept (prep X matrix 16 MB, model run records and curves, score tables); 6.7 GB of exported logits was transient and deleted after scoring | exploratory idea 5 (`docs/plans/explore_idea5_plan.md`); small CSVs committed | **Yes** (recomputable in about 2.5 h on GPU 1 plus 1 h CPU) |

| `results/explore_ideas12/` (12 Nov 2026) | 0.52 GB kept (rds readout caches 0.36 GB, attention graphs 0.14 GB, memberships and tables); 1.9 GB of uint8 count matrices deleted after scoring | exploratory ideas 1-2 (consensus co-membership sweep, partner-less gate); regenerable from `analysis/explore_ideas12/` | **Yes** after the numbers are archived |

| `results/explore_stack/` (12 Nov 2026) | 0.8 MB measured (CSV, md, logs; rds readout caches under `cache/` about 0.3 MB, not committed; no matrices kept) | exploratory stacked a3r layers (`docs/plans/explore_stack_plan.md`); small tables force-added to git | **Yes** for `cache/` (recomputable in about 40 min); keep the committed tables |

| `data/complex_confirm_raw` (13 Nov 2026) | 1.3 MB (Complex Portal 559292 complextab, UniProt 559292 map, SHA256SUMS) | complex confirmation design (`docs/plans/complex_confirm_design.md`); gitignored. `results/complex_confirm` (planned) expected under 1 GB (three net4 score caches of about 0.14 GB plus tables). Counted | **Yes**, re-downloadable (hashes in the design); caches deletable after the numbers are archived |

| `results/scoreboard/` (13 Nov 2026) | 0.4 MB measured (CSV tables, equivalence gates, logs; no matrices kept) | scoreboard figure compute (`docs/plans/scoreboard_plan.md`, `analysis/scoreboard/`); small tables force-added to git; figures `docs/whitepapers/figures/graphical_abstract_p3{,b}.{png,pdf}` (0.9 MB) | No (the committed record of the scoreboard cells computed now) |
