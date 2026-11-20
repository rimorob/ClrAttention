# Adaptive (label-free) setting studies for aCLR / tcaCLR

Implements items 3, 8, 9, 11 and 12 of Revision 2 of `docs/plans/tcaclr_depth_grid_plan.md` (review: `docs/plans/fable_adaptive_review.md`). EXPLORATORY; no truth is read; smoke-tested only, nothing submitted. The operator is always the grid library's (`analysis/tcaclr_depth/grid_lib.R`: `gl_make_P`, `gl_walk`, `gl_variants`, `gl_ap_T`); there is no second operator builder. Run from the repo root.

## Scripts

| Script | What it does | Main output |
|---|---|---|
| `adaptive_lib.R` | shared helpers: data load, per-gene permutation, rank-normal transform, blocks and halves, symmetrised-by-degree operator, RSpectra eigen-solver | |
| `spectral_edge.R` | noise edge of the operator on per-gene permuted data (mean and range of the largest non-trivial eigenvalue over `--perms` permutations), observed top eigenvalues, number above the edge, t* = ln(eps)/ln(lambda_edge) for eps 0.05 and 0.01 | `spectral_edge_<tag>.tsv`, `..._spectrum.tsv`, `.png` |
| `graph_signal_lik.R` | blocked held-out Gaussian-copula graph-signal likelihood in the operator eigenbasis (t, amplitude, noise floor), surface over the 17 grid depths at each keff, best t and one-standard-error t, N_eff and G_eff | `gs_surface_<tag>.tsv`, `gs_selected_<tag>.tsv`, `gs_neff_<tag>.tsv`, `gs_pooled_<tag>.tsv`, `.png` |
| `planted_sim.R` | planted-module simulation in real null data; CLR, aCLR, tcaCLR, diffusion-PMI (exploratory) over a small keff x depth grid; AUPR per cell, oracle cell, regret of the edge rule and of the likelihood rule | `ps_cells/rules/oracle/regret/summary/datasets_<tag>.tsv` |
| `negative_control.R` | per-gene permuted copies through both rules; flags "no structure" | `nc_perm_<tag>.tsv`, `nc_summary_<tag>.tsv` |
| `run_adaptive.R` | runs the first, second and fourth on one development matrix (restartable per step) | files above in `--out` |
| `smoke_test.R` | checks (operator is the grid library's; P_sym symmetric and equals dense eigen; planting keeps marginals; AP from positives equals `ap_fast`; permuted control) | |

## Inputs and outputs

- Inputs: `results/baselines_m3d/expression.tsv` (4,297 genes x 466 experiments; blocks = the stored laboratory splits `results/lab_splits/split_assignments.csv`, 39 laboratories) and `results/baselines_k562ess_truthgenes_genie3/expression.tsv` (5,671 genes x 707 knockdowns; blocks = k-means clusters of the knockdown profiles, k = round(N_eff) = 40, seed fixed).
- Outputs: `results/tcaclr_depth/adaptive_m3d/`, `adaptive_k562/`, `adaptive_planted_m3d/` (job files below). Seeds are fixed (`AD$seed` 20261121 plus offsets).
- Jobs (header style of `analysis/tcaclr_depth/jobs/`): `jobs/7440_tcd_adaptive_m3d.sh` (mem 12, cpus 16), `jobs/7441_tcd_adaptive_k562.sh` (mem 20, cpus 16), `jobs/7450_tcd_planted_sim.sh` (mem 12, cpus 8). Run the planted simulation and the negative controls before any real-data selection is used (plan item 11).

## Declared resources and measured smoke-test numbers (8 threads, nice 19, full-size matrices, one machine shared with other work)

| Script (smoke scale) | Wall | Peak RSS |
|---|---|---|
| spectral_edge, M3D 4,297 genes, 3 keffs, 2 permutations, k_obs 100 | 108 s | 1.65 GB |
| spectral_edge, K562 5,671 genes, 2 keffs, 2 permutations | 166 s | 2.93 GB |
| graph_signal_lik, M3D, 1 split (2 directions), keffs 10 and 80, r 400 | 150 s | 1.72 GB |
| graph_signal_lik, K562, 1 split, keffs 10 and 80, r 400 | 240 s | 2.89 GB |
| graph_signal_lik, K562, 1 split, keff 10, r 1000 | 531 s | 2.49 GB |
| negative_control, M3D, 1 control (nested edge, 2 permutations), keff 10, 1 split | 127 s | 1.52 GB |
| planted_sim, M3D background 3,000 genes, 1 dataset (keffs 10, 28; 7 depths; 2 edge permutations) | 158 s | 1.89 GB |

Unit costs: one MI fit about 11 s (M3D) and 22 s (K562); one operator plus top-100 eigenvalues 7 to 10 s; top-400 eigenvectors about 20 to 25 s (M3D); graph-signal likelihood about 30 s (M3D) and 45 s (K562) per keff and direction.
Extrapolated job wall times (4 workers): 7440 about 1.5 h, 7441 about 2 h, 7450 about 1 h. The declared memory is 4 workers x measured per-worker peak plus margin. If a K562-like matrix of about 9,600 genes is used, scale memory and time by about (9,600/5,671)^2 = 2.9 (per worker about 8 GB, so use 2 workers or declare about 40 GB).

## Definitions fixed in the code

- Operator: P = 0.5 I + 0.5 A (the grid library). A is not reversible, so D^(1/2) P D^(-1/2) with D = stationary mass is not symmetric. Symmetrised form: W = (A + A')/2, d = rowSums(W), P_sym = 0.5 I + 0.5 D^(-1/2) W D^(-1/2) (symmetric, spectrum in [0, 1]). The second-largest modulus of the non-symmetric P is reported alongside (`edge_nonsym_mean`, `obs_nonsym2`).
- Edge: the largest non-trivial eigenvalue; n_above counted against the mean edge and against the largest edge over the permutations; `obs_lambda2_perm_p` = (1 + #permutations with lambda_2 at least the observed)/(P + 1).
- Graph-signal likelihood: rank-normal per gene within each half; graph and (amplitude, noise floor) fitted on one half, scored on the other half, both directions pooled; the unreachable-entry floor is the noise variance (a mode whose lambda^t vanishes keeps variance sigma^2); trailing modes beyond r = 400 get variance sigma^2 + c lambda_r^t. Score per gene and experiment, averaged within held-out blocks, then over blocks. Selection: best t, one-SE t (paired SE, and a second column with the SE of the best score itself).
- Negative control: `--seperms 0` (job setting) takes the edge from the other controls (leave-one-out); `--seperms > 0` is the nested version.
