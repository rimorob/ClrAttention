# tcaCLR2 (exploratory)

Implementation of `docs/plans/tcaclr2_definition.md`: tcaCLR with the walk depth set by the permuted-data noise edge instead of t = 9. The definition is not changed here. Exploratory; no truth is read in computing t2.

## Scripts (run from the repo root)

| script | does |
|---|---|
| `edge_depth.R` | edge of the registered operator (keff 10, cap 50, lazy 0.5) over 10 per-gene permutations (seeds `AD$seed + 7000 + p`, the same permuted matrices as `adaptive/spectral_edge.R`), mean of the largest non-trivial eigenvalue of P_sym; `t2 = max(1, round(ln 0.05 / ln edge))`. Reuses `adaptive_lib.R` (`ad_permute_genes`, `ad_fit`, `ad_P_from_fit`, `ad_Psym`, `ad_eigs`) and `grid_lib.R` (`gl_make_P`); the only change from `spectral_edge.R` is cap 50 instead of `gl_cap(keff)` = 40. |
| `make_scores_tcaclr2.R` | `genes.txt` + packed upper-triangle `tcaCLR2.scores` (format of `tcaclr_confirm/make_scores.R`), registered `sk_a3r` on the registered sequential chain P^t2, t2 read from the edge TSV. `--t 9 --name tcaCLR_t9 --compare <score dir>` is the self-test against the stored registered tcaCLR. |
| `score_tcaclr2.R` | pooled co-membership AUPR (stage-2 truths `gl_truths`, `ap_fast` tie convention) of registered tcaCLR (t = 9, cap 50) and tcaCLR2, ratio, and the gene delete-half jackknife of the registered scorer (`tc_draw_halves`, `tc_jack_method`, `tc_interval`, seed 20261117, 500 halves) for those two cells only. |
| `collect_tcaclr2.R` | merges edge TSVs and score tables into `results/tcaclr2/summary.csv`. |
| `run_tcaclr2.sh` | `edge`, `score` or `files` stage over the five datasets (`k562_corum net4_cp_yeast armA_cp_human armB_cyc2008 m3d_tf_regulons`). |
| `jobs/7500_tcaclr2_edge_depth.sh`, `jobs/7510_tcaclr2_score.sh` | job files (7510 after 7500). `analysis/tcaclr2/jobs` is gitignored (`jobs/` rule): force-add. The `files` stage (tcaCLR2.scores for all five datasets) is not in either job file. |

## Inputs
Matrices and truths exactly as stage 2 of `analysis/tcaclr_depth` (`grid_truths.R`): K562 `results/baselines_k562ess_truthgenes_genie3/expression.tsv` with CORUM; DREAM5 net4 yeast with the frozen Complex Portal yeast truth; HepG2 Arm A and Kemmeren Arm B expression and frozen truth files (SHA-256 guard); M3D `results/baselines_m3d/expression.tsv` with RegulonDB regulons (scopes SC and all). Stored registered scores `results/tcaclr_confirm/armA|armB/scores` are read (never written) for the equality checks.

## Outputs (only `results/tcaclr2/`; the scripts refuse any output path containing `tcaclr_confirm`)
- `edge_<dataset>.tsv`: one row: dataset, G, N, edge, edge_min, edge_max, t2, tstar (unrounded), keff, cap, perms, eps, seed, n_constant_genes, top_eig_min (sanity: trivial eigenvalue); `edge_<dataset>.log`.
- `score_<dataset>/results.csv` (one row per scope): AUPR tcaCLR and tcaCLR2, `ratio`, `diff`, `diff_*` (interval of the difference: se_gene, gene_lo/hi normal, pct_lo/hi percentile, skew, n_rep, n_na) and `ratio_*` (same `tc_interval` applied to the replicate ratios; an addition of this script, not in stage 2); `checks.csv` (Arm A/B: in-memory tcaCLR vs stored file; with `--scores_dir`: in-memory tcaCLR2 vs `tcaCLR2.scores`); `score.log`.
- `summary.csv`; `scores_<dataset>/` (`genes.txt`, `tcaCLR2.scores`) from the `files` stage.

## Declared resources (measured full-size peaks plus margin)
- 7500: mem 16 GB, cpus 8. Edge uses 4 workers x 2 clr threads (Arm A: 2 workers x 4 threads). Measured worker peak 1.2 / 2.0 / 1.9 / 2.1 / 5.2 GB (M3D, K562, net4, Arm B, Arm A), so total about 5-9 GB for four workers, about 11.5 GB for Arm A.
- 7510: mem 10 GB, cpus 8 (one dataset at a time). Measured process peak 2.0 / 2.7 / 2.8 / 3.1 / 5.7 GB (M3D, K562, net4, Arm B, Arm A).

## Smoke times (full size, `nice 19`, 8 cores, one script at a time)
| dataset (G x N) | edge, wall (workers x threads) | edge | t2 | score, wall | make_scores, wall |
|---|---|---|---|---|---|
| M3D (4297 x 466) | 111 s (4 x 2) | 0.88592 | 25 | 152 s | not run |
| K562 (5671 x 707) | 259 s (4 x 2) | 0.90326 | 29 | 152 s | not run |
| net4 yeast (5950 x 536) | 220 s (4 x 2) | 0.86766 | 21 | 111 s | not run |
| Kemmeren Arm B (6169 x 1283) | 478 s (4 x 2) | 0.87625 | 23 | 138 s | 152 s (t2), 101 s (self-test t = 9) |
| HepG2 Arm A (9622 x 739) | 722 s (2 x 4) | 0.91049 | 32 | 303 s | not run |

Total: edge stage about 30 min, score stage about 14 min (plus the collector, seconds); at full size these smoke numbers are the full-size times.

## Checks
- Self-test: `make_scores_tcaclr2.R --dataset armB_cyc2008 --t 9 --name tcaCLR_t9 --compare results/tcaclr_confirm/armB/scores`: genes identical, maximum absolute difference 0 over 19,025,196 pairs.
- `score_tcaclr2.R` in-memory tcaCLR (t = 9) vs stored tcaCLR.scores on the truth pairs: Arm B 2.7e-15, Arm A 1.8e-15; in-memory tcaCLR2 vs `tcaCLR2.scores` of `make_scores_tcaclr2.R` (Arm B, t2 = 23): 3.6e-15 (squaring chain vs sequential chain, rounding only).
- Registered tcaCLR AUPR reproduces stage 2 (K562 0.105435, Arm B 0.207586374919325).

## Flags
- The edge differs a little from the cap-40 values of `spectral_edge.R` (M3D 0.8849 -> 0.8859, K562 0.9027 -> 0.9033), as expected from a wider operator; t* 24.5 -> 24.7 (M3D), 29.4 -> 29.4 (K562).
- `round()` is R's (half to even); no edge is near a half-integer boundary (unrounded t*: 21.10, 22.68, 24.73, 29.45, 31.95).
- The definition note does not say whether constant genes are removed before the edge (`ad_load` drops them in `spectral_edge.R`); none of the five matrices has a constant gene, so the question does not arise here.
- The ratio interval is the same `tc_interval` on replicate ratios; the definition note only says "the gene-jackknife interval of the existing scorer", which stage 2 applies to the difference. Both are written.
