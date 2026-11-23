# Global-temperature ladder and multi-scale tcaCLR scores (exploratory, descriptive)

Design and reading rules: `docs/plans/global_temperature_ladder_plan.md` (written before any score). This directory holds the code, the job files (NOT submitted) and this note. Nothing here changes the registered tcaCLR or tcaCLR2 or any width-ladder file; every truth is spent.

## What was built

| File | Role |
|---|---|
| `gt_lib.R` | the global-temperature operator (recipe-ladder V2 lines moved into functions unchanged: `gt_topmat`, `gt_w`, `gt_effv`, `gt_tau_for_eff`, `gt_make_P`), the log2 tau ladder (`GTL`), the operator list, the permutation-edge loop `gt_edges` (`wl_edges` with the operator builder swapped). Sources `analysis/width_ladder/wl_lib.R` (loaders, `wl_eval`, `wl_tca_vec`, edge machinery); no second operator or scorer. |
| `gt_probe.R` | CLR-only probe (no truth): effective-partner distribution over log2 tau -8 to +4; used to fix the ladder endpoints. |
| `gt_run.R` | `--mode gate` (reference operator + V2 point, gates 1 and 2, real truths, gates only), `--mode edge` (real-data temperature table, permutation edges of reference, V2 and the 13 taus, t table), `--mode score` (inline gates, then the ladder under R9 and RE, V2 reference, MS1-3, comparators). Smoke switches `--synthetic 1`, `--log2taus`, `--perms`, `--tag`. |
| `gt_assemble.R` | response-surface tables, best temperature with edge flags, effective-partner table, t table, ledger, gate summary, heatmaps and curves. |
| `jobs/81xx_*.sh` | edge jobs 8100-8105 (m3d, ecomac, k562, yeast, kemmeren, hepg2), score jobs 8110-8115 (each `# after:` its own edge job), assemble 8190 (`# after:` all six score jobs). Format `# mem`, `# cpus`, comma-separated `# after`. |

Ladder: log2 tau = -4, -3.5, ..., +2 (13 points; mean effective partners about 1.2 to 44 on M3D/EcoMAC), absolute and identical on every dataset; V2's own temperature of each dataset (mean effective partners 10; M3D log2 -0.886, EcoMAC -0.692, yeast -1.494, HepG2 -1.346) is a labelled reference point `V2`, scored under both depth rules, outside MS and outside the best-tau argmax. Depth rules R9 (t = 9) and RE (permutation edge at that tau, permuted operator built with the same tau). Comparators CLR, tcaCLR (keff 10, cap 50, t 9), tcaCLR2 (stored t2). MS1-3 over the 13 ladder points (MS3: log2 tau -2.5 to +1).

## Inputs

Matrices and truths exactly as the width ladder (`wl_load` / `wl_matrix`); stored edge files `results/tcaclr2*/edge_*.tsv` (reference edge check and the tcaCLR2 comparator's t2); reference tables `results/tcaclr2_scoreboard/cells_{reference_tcaclr,tcaclr2}_all.csv`, `results/direct_edge_diag/{gate_direct_edge,variants_*}.csv`, `docs/whitepapers/figures/scoreboard_data.csv` (CLR); gate 2: `results/recipe_ladder/variants_{m3d,ecomac}.csv` and `variant_info_*.csv`.

## Outputs (`results/global_temp_ladder/`)

Per dataset `<ds>/`: `tau_table.csv` (edge job: V2's tau and effective partners per tau on the real data: mean, median, 10th, 90th percentile, min, max), `tau_table_score.csv`/`tau_table_gate.csv` (recomputed by the score/gate job; the score job stops if its V2 tau or effective-partner means differ from the edge job's), `edge.tsv` (edge, t, capped flag per operator), `edge_perm.csv`, `edge_gate.csv`, `gate.csv`, `cells.csv` (long: dataset, rule, score, tau, t, truth, scope, metric, value), `run.csv`, logs. Assembled: `surface_all.csv`, `best_tau.csv`, `effective_partners.csv`, `t_by_tau.csv`, `ledger.csv`, `gate_summary.csv`, `heatmap_<ds>_<rule>.png`, `curves_<ds>_<rule>.png` (curves: dotted grey CLR, dashed blue tcaCLR, dash-dot red tcaCLR2, solid green/purple/orange MS1/MS2/MS3, brown triangle and vertical line the V2 reference). Files tagged `_smoke` are smoke output, not results; `<ds>/cells_gate.csv` and `tau_table_gate.csv` are from `--mode gate`.

## Gates (real truths, gates only; run now on M3D, EcoMAC and yeast; tolerance 1e-9)

| Dataset | gate 1 (registered operator vs stored tables), cells / max abs diff | gate 2 (V2 vs `results/recipe_ladder/variants_*.csv`), cells / max abs diff |
|---|---|---|
| M3D | 48 / 3.6e-14 | 7 / 4.4e-15 (SC direct-edge ratio 7.8177, prints 7.82) |
| EcoMAC | 48 / 4.3e-14 | 7 / 4.9e-15 (SC direct-edge ratio 6.6574, prints 6.66) |
| yeast net4 | 3 / 6.0e-14 | none stored |

K562, Kemmeren, HepG2 run the same gate 1 inline at the start of their score jobs (they stop the run on failure); not run separately now. Gate 2 cells: pooled direct-edge ratio SC and all, per-TF median AUROC SC and all, hub-controlled ratio SC and all, and V2's tau (`tau_global` of `variant_info`).

First gate attempt (`m3d/gate_FAILED.csv`, kept): gate 2 failed on the median per-TF AUROC (0.6046 against stored 0.6025) while the ratios agreed to 1e-14. Cause: a metric-definition difference, not the operator. The recipe ladder's stored median AUROC is `tf_rows` of `analysis/direct_edge_diag/lib.R` (TFs with at least 5 targets); the ladder's `tf_auroc_median` is the scoreboard `sb_per_tf_median_auroc` (at least 3 targets), which is what gate 1 and the width ladder use. Gate 2 now recomputes the recipe figure with `tf_rows` on the V2 score (metric name `tf_auroc_median_tf_rows`); the two definitions are never mixed in a table (the ladder's tables carry only `tf_auroc_median`).

## Assemble bug (width ladder) and what was changed

The width ladder's final heatmap failed with "figure margins too large". Reproduced: `wl_assemble.R` opens each heatmap with height = 120 + 46 n px at 130 dpi (n = number of metric rows) and `par(mar = c(5, 17, 3, 1))`; for the one-row datasets (yeast, Kemmeren, HepG2: n = 1, 166 px = 1.3 in) the 8 margin lines (1.6 in) exceed the device. `gt_assemble.R` sizes the device in inches from the margins (height = max(4.5, 0.2 (mar bottom + top) + 0.32 n + 1.5) in, width 11.5 in) and wraps every figure in `try`, so one failing figure cannot stop the tables. Tested on the EcoMAC smoke (12 rows) and a one-row copy. `analysis/width_ladder/wl_assemble.R` was NOT changed (the width ladder's tables, written before the heatmaps, are unaffected; its heatmaps for the failing datasets were not drawn).

## Smoke results (nice 19; limit 8 cores / 20 GB; whole-process-tree peak = sum of RSS over the command and all descendants, shared pages counted per process, so conservative)

| Run | Settings | Wall | Tree peak |
|---|---|---|---|
| edge, EcoMAC | 3 permutations x 5 operators (reference, V2, log2 tau -3, -1, +1), 3 workers x 2 threads; real-data phase 100 s, permutations 125 s | 259 s | 5.9 GB (main process 1.3 GB, worker 2.0 GB) |
| score, EcoMAC, synthetic truths | log2 tau -3, -1, +1 + V2 + reference, 8 threads; 204 cells, all six MS scores | 277 s | 3.2 GB |
| edge, HepG2 (largest) | 2 permutations x 3 operators (reference, V2, log2 tau -1), 2 workers x 4 threads; real-data phase 133 s, permutations 151 s | 315 s | 10.9 GB (worker 5.9 GB) |
| gate, M3D / EcoMAC / yeast (real truths) | reference + V2 at t = 9 (and stored t2 for the reference), 8 threads | 98 s / 127 s / 83 s | 3.2 / 3.3 / 3.9 GB |
| assemble | EcoMAC smoke and a one-row copy | 2 s each | |

First edge smoke attempt hung: `mclapply` with ONE task runs in the parent, which then fitted before forking the permutation workers (the known fork-after-fit deadlock). The real-data phase of the edge job now uses `parallel::mcparallel` (always forks). Also fixed: the ladder list is `GTL` (a first name `GT` overwrote the grid-truth list of `grid_truths.R`).

## Declared resources and full-run estimate

Score jobs fork no workers (BLAS threads only), declared at about 2x the width-ladder real peaks (M3D 4.9, EcoMAC 4.6, K562 8.4, yeast 5.2, Kemmeren 5.2, HepG2 13.1 GB main process; same arrays here). Edge jobs fork workers: declared at least 2.5 x the whole-tree peak (EcoMAC measured 5.9 GB with 3 workers, about 7.5 GB projected for 4 workers: 24 GB; HepG2 measured 10.9 GB: 30 GB; K562, yeast, Kemmeren not measured, 3 workers projected about 11 GB: 30 GB).

| Job | mem GB | cpus | est. wall |
|---|---|---|---|
| 8100 edge M3D | 24 | 8 (4 workers x 2) | 12 min |
| 8101 edge EcoMAC | 24 | 8 | 12 min |
| 8102 edge K562 | 30 | 6 (3 x 2) | 15 min |
| 8103 edge yeast | 30 | 6 | 15 min |
| 8104 edge Kemmeren (optional) | 30 | 6 | 20 min |
| 8105 edge HepG2 (optional) | 30 | 8 (2 x 4) | 30 min |
| 8110 score M3D | 10 | 8 | 14 min |
| 8111 score EcoMAC | 10 | 8 | 14 min |
| 8112 score K562 | 16 | 8 | 30 min |
| 8113 score yeast | 12 | 8 | 17 min |
| 8114 score Kemmeren (optional) | 12 | 8 | 19 min |
| 8115 score HepG2 (optional) | 28 | 8 | 55 min |
| 8190 assemble | 4 | 2 | 1 min |

Estimates: edge from the smoke (per permutation about one fit plus 15 operators at 5 to 10 s each, 10 permutations over the workers, plus the real-data fit) and from the width ladder's full edge times (M3D 377 s, EcoMAC 631, K562 839, yeast 819, Kemmeren 1199, HepG2 2118 s); score from the width ladder's full score times (M3D 851 s, EcoMAC 807, K562 1744, yeast 1038, Kemmeren 1131, HepG2 3286 s: same number of evaluations, 15 operators instead of 14, cheaper operators with fixed cap 50, depth up to 256 at the concentrated end), nothing measured at full size here. Critical path: HepG2 edge then score, about 1.5 h; the whole set within the 90 GB / 56 thread budget (memory is the binding limit: about three edge jobs at a time) about 2.5 h; without HepG2 and Kemmeren about 1.2 h. If the optional jobs are not submitted, drop them from the `# after:` list of 8190 (the assemble skips datasets without cells).

## Things to know

- Observed in the smokes (EcoMAC, 3 permutations, not a result): the RE depth at V2's temperature is t = 52 (the registered operator: 18); at log2 tau -3 the permuted operator has edge 0.99996 and t is capped at 256 and flagged `capped`; at log2 tau +1, t = 9. HepG2 smoke: V2 temperature t = 240. So RE walks run deep at the concentrated end.
- The ladder is the same absolute tau on every dataset, so the effective-partner distribution differs by dataset (V2's tau does too); `effective_partners.csv` records it per dataset and temperature.
- V2's temperature is data dependent and is computed by the edge job (real-data fit in a forked child, no truth), re-computed by the score job, which stops if it differs.
- Hub-controlled cells use one stratifier per dataset (stationary distribution of the registered operator), not one per temperature.
- Kemmeren and HepG2 have no stored CLR value, so CLR is not gated there.
- Per-TF AUROC in all ladder tables is `tf_auroc_median` (scoreboard definition, at least 3 targets), which differs from the 5-target figure stored in the recipe ladder tables; do not compare them directly.
