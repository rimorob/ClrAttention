# Width ladder and multi-scale tcaCLR scores (exploratory, descriptive)

Design and reading rules: `docs/plans/width_ladder_plan.md` (written before any score). This directory holds the code, the job files (not submitted) and this note. Nothing here changes the registered tcaCLR or tcaCLR2; every truth is spent.

## What was built

| File | Role |
|---|---|
| `wl_lib.R` | width loop helpers, multi-scale accumulation hooks, dataset/truth loaders, metric evaluator. Reuses `gl_make_P`, `gl_walk`, `gl_zmat` (tcaclr_depth/grid_lib.R), the edge machinery (`ad_fit`, `ad_permute_genes`, `ad_Psym`, `ad_eigs`, `ed_t2` of tcaclr2/edge_depth.R), `ap_fast`, `sb_per_tf_median_auroc`, `rz_hub_stratify`, the truth builders of `sb_load_ecoli` / `build_truths_ecoli` / `ideas47_lib.R` / `grid_truths.R`. No second operator or scorer. |
| `wl_run.R` | `--mode gate` (reference operator vs registered tables), `--mode edge` (permutation edges of reference + 13 widths, t table), `--mode score` (ladder, R9 and RE, MS1-3, comparators). `--synthetic 1` swaps all truths for random planted ones (smoke only). |
| `wl_assemble.R` | response-surface tables, best width per truth with edge flags, t table, ledger, heatmaps and curves (PNG). |
| `jobs/80xx_*.sh` | 6 edge jobs (8000-8005), 6 score jobs (8010-8015), assemble (8090). Format `# mem`, `# cpus`, comma-separated `# after`. |

Scores per dataset: 13 widths x {R9: t = 9, RE: t from the per-width edge}, MS1 / MS2 / MS3 per rule, comparators CLR, tcaCLR (keff 10, cap 50, t 9), tcaCLR2 (keff 10, cap 50, stored t2). Cap = min(4 keff, G - 1). The edge is computed once per permutation fit for all operators (10 permutations, registered seeds). Metrics per score: M3D / EcoMAC 12 (pooled AUPR/base raw and hub-controlled on SC, all, complexes; direct TF edges SC, all: pooled ratio, median per-TF AUROC, hub-controlled pooled ratio), K562 10 (5 truths x raw/hub), yeast / Kemmeren / HepG2 1.

Inputs: matrices and truths exactly as the registered scoreboard / stage 2 (see plan); stored edge files `results/tcaclr2*/edge_*.tsv` (t2 of tcaCLR2 comparator and the edge check); reference tables `results/tcaclr2_scoreboard/cells_{reference_tcaclr,tcaclr2}_all.csv`, `results/direct_edge_diag/{gate_direct_edge,variants_*}.csv`, `docs/whitepapers/figures/scoreboard_data.csv` (CLR).

Outputs (`results/width_ladder/<ds>/`): `edge.tsv` (edge, t, capped flag per width), `edge_perm.csv`, `edge_gate.csv`, `gate.csv`, `cells.csv` (long: dataset, rule, score, keff, cap, t, truth, scope, metric, value), `run.csv`, logs. Assembled in `results/width_ladder/`: `surface_all.csv` (+ per dataset), `best_width.csv` (best single width, `best_at_edge` flag: low edge 2 / next to low edge 3 / next to high edge 96 / high edge 128 / interior), `t_by_width.csv`, `ledger.csv`, `gate_summary.csv`, `heatmap_<ds>_<rule>.png` (log2 of score / tcaCLR; red above, blue below; dot = best single width; columns: widths, MS1-3, CLR, tcaCLR, tcaCLR2), `curves_<ds>_<rule>.png` (metric vs keff; dotted grey CLR, dashed blue tcaCLR, dash-dot red tcaCLR2, solid green / purple / orange MS1 / MS2 / MS3).

## Gate (real truths, gate only; all passed, tolerance 1e-9)

Reference operator (keff 10, cap 50) through the ladder code: tcaCLR (t 9), tcaCLR2 (stored t2), CLR where a stored value exists.

| Dataset | cells compared | max abs diff |
|---|---|---|
| M3D | 48 | 3.6e-14 |
| EcoMAC | 48 | 4.3e-14 |
| K562 | 30 | 2.0e-14 |
| yeast net4 | 3 | 6.0e-14 |
| Kemmeren | 2 (no stored CLR) | 5.0e-14 |
| HepG2 | 2 (no stored CLR) | 4.4e-15 |

M3D and EcoMAC include the pooled raw and hub-controlled cells, tf_aupr_ratio, tf_auroc_median (`cells_*`, `gate_direct_edge.csv`) and the hub-controlled direct ratio (`variants_*.csv`). Edge gate: the full M3D edge job (10 permutations, 14 operators) reproduced the stored edge exactly (0.88591846, t = 25); EcoMAC, K562, yeast, Kemmeren, HepG2 are checked inside their edge jobs (the job stops on a difference above 1e-6 or a different t2). `--mode score` refuses to run without `gate.csv`.

## Smoke (synthetic truths, real matrices, nice 19)

| Run | Settings | Wall | Peak (whole process tree) |
|---|---|---|---|
| edge, M3D | 3 permutations x 4 operators, 3 workers x 2 threads | 64 s | 5.9 GB |
| score, M3D | widths 4, 16, 64 + reference, 8 threads | 250 s (180 cells, all 6 MS scores) | 3.6 GB |
| edge, HepG2 | 2 permutations x 3 operators (reference, 2, 128), 2 workers x 4 threads | 222 s | 11.3 GB |
| score, HepG2 | widths 2, 128 + reference, 8 threads | 733 s | 9.1 GB |
| assemble | smoke M3D | under 1 min | |

Full real runs done only for the gate and for the M3D edge stage (344 s, 9.0 GB whole tree with 4 workers x 2 threads). The smoke files carry the tag `_smoke` and are not results.

## Declared resources and wall-time estimates

Score jobs fork no workers (BLAS threads only); edge jobs fork permutation workers and are declared from the measured whole-tree peak times at least 2.5 (M3D 9.0 GB -> 24; HepG2 11.3 GB -> 30; k562, yeast, kemmeren are interpolated, not measured: 3 workers about 12 GB -> 30; EcoMAC taken as M3D).

| Job | mem GB | cpus | est. wall |
|---|---|---|---|
| 8000 edge M3D | 24 | 8 | 6 min (measured) |
| 8001 edge EcoMAC | 24 | 8 | 8 min |
| 8002 edge K562 | 30 | 6 | 12 min |
| 8003 edge yeast | 30 | 6 | 13 min |
| 8004 edge Kemmeren (optional) | 30 | 6 | 16 min |
| 8005 edge HepG2 (optional) | 30 | 8 | 50 min |
| 8010 score M3D | 10 | 8 | 15 min |
| 8011 score EcoMAC | 10 | 8 | 15 min |
| 8012 score K562 | 16 | 8 | 40 min |
| 8013 score yeast | 12 | 8 | 30 min |
| 8014 score Kemmeren (optional) | 12 | 8 | 30 min |
| 8015 score HepG2 (optional) | 24 | 8 | 65 min |
| 8090 assemble | 4 | 2 | 1 min |

Estimates are extrapolated from the smoke and gate timings (main costs: dense products, `.softmax_rows` at caps up to 512, rank percentiles over all pairs, hub stratification), not measured at full size, except M3D edge. Critical path: HepG2 edge then score, about 2 h; the whole set within the 90 GB / 56 thread budget about 2.5 to 3 h; without HepG2 and Kemmeren about 1.2 h. HepG2 and Kemmeren are separate optional jobs (HepG2 is the slow one). If the optional jobs are not submitted, drop them from the `# after:` list of 8090.

## Things to know

- Width 10 is not on the ladder and the registered cap (50) differs from 4 x keff (40): the registered cells are a separate reference operator, plotted as comparators.
- t from the edge rule is large at the narrow end (M3D keff 2: edge 0.980, t = 147; HepG2 keff 2 smoke: t = 194) and 9 to 10 at keff 128. The cap t <= 256 (flagged `capped`) did not trigger in any run so far.
- Hub-controlled cells use one stratifier per dataset (stationary distribution of the registered operator), not one per width.
- Kemmeren and HepG2 have no stored CLR value, so CLR is not gated there. The smoke runs load RegulonDB tables inside `sb_load_ecoli` but score only synthetic truths.
- The M3D edge result in `results/width_ladder/m3d/edge.tsv` already exists (it is the first part of job 8000, identical when rerun).
