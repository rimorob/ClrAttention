# tcaCLR global-temperature ladder and multi-scale scores: exploratory, descriptive

Written 23 Nov 2026 before any real score of the ladder exists (no truth has been scored with this operator; the probe `analysis/global_temp_ladder/gt_probe.R` reads only the MI/CLR fit of M3D and EcoMAC, no truth). EXPLORATORY and DESCRIPTIVE: every truth used here is spent, nothing here is a test, and no setting of the registered tcaCLR (keff 10, cap 50, t = 9, clipped) or tcaCLR2 changes. Any single temperature that looks best on these truths is a selection and needs a fresh test; the multi-scale scores below have no selection step. Requested by the PI after the recipe ladder (`results/recipe_ladder/stage_table.md`, variant V2) found that replacing the per-gene temperature (every row spread over exactly keff = 10 effective partners) by ONE global temperature lifts tcaCLR above CLR on direct TF-to-target edges on M3D and EcoMAC; the ladder asks how that depends on the temperature and whether it transfers to the other truths.

## Operator (reused, not rebuilt)

The operator is variant V2 of `analysis/recipe_ladder/01_ladder.R`, with its lines moved into functions unchanged (`gt_topmat`, `gt_w`, `gt_effv`, `gt_tau_for_eff`, `gt_make_P` in `analysis/global_temp_ladder/gt_lib.R`):

- candidates of gene g: its top 50 POSITIVE clipped-CLR scores (self excluded; `order(decreasing = TRUE)`, ties as in the recipe ladder); fewer if fewer are positive; a gene without a positive score attends to itself;
- weights exp((s - max s) / tau) over the candidates, row-normalised, then the registered finish step (`gl_finish_A`), lazy weight 0.5 (`P = 0.5 I + 0.5 A`), the registered walk (`gl_walk`), second CLR on log P^t over reachable pairs, clipped at 0, Stouffer, symmetric read-out (`wl_tca_vec`);
- tau is ONE number per run (units of the clipped CLR score), the same for every gene.

V2's own temperature of a dataset is the one that gives a mean over rows of exp(entropy) of exactly 10, found by the recipe ladder's bisection (log tau in [log 1e-4, log 1e4], 80 iterations; mean over rows that have at least one candidate). It is data dependent (M3D tau = 0.5412, log2 = -0.886; EcoMAC tau = 0.6188, log2 = -0.692), so it is carried as a labelled REFERENCE point (id `V2`) of each dataset next to the ladder, not as a ladder point; the ladder itself is the same absolute tau grid on every dataset.

## Settings

| Axis | Values |
|---|---|
| tau ladder (half-octave steps, absolute, the same on every dataset) | log2 tau = -4, -3.5, ..., +2 (13 points), tau = 0.0625 to 4 |
| reference point | `V2`: the dataset's own tau for mean effective partners 10 (bisection above); scored under both depth rules; not in the multi-scale scores, not in any best-tau argmax |
| candidate cap | 50 (as the registered operator and V2) |
| lazy weight | 0.5; clipS and clipO TRUE (registered clipping); Stouffer, symmetric read-out |
| depth rule R9 | fixed t = 9 (tcaCLR style) |
| depth rule RE | t from the permutation-edge rule at that temperature (tcaCLR2 style): edge = mean over 10 per-gene permutations (seeds ED seed + 7000 + p, as registered) of the largest non-trivial eigenvalue of the symmetrised lazy operator built from the PERMUTED data with the SAME tau and the same top-50 rule (tau is not re-solved on the permuted data); t = max(1, round(ln 0.05 / ln edge)); one fit per permutation, all temperatures and the registered reference operator read from it |

Endpoints (fixed here from the CLR-only probe, `results/global_temp_ladder/probe/probe_{m3d,ecomac}.csv`, mean effective partners = mean over rows of exp(entropy)):

| log2 tau | M3D mean (median; 10th-90th pct) | EcoMAC mean (median; 10th-90th pct) |
|---|---|---|
| -4.0 | 1.21 (1.00; 1.00-1.84) | 1.16 (1.00; 1.00-1.66) |
| -3.0 | 1.54 (1.05; 1.00-2.62) | 1.41 (1.02; 1.00-2.25) |
| -1.0 | 8.81 (4.55; 1.02-24.1) | 7.02 (3.07; 1.00-19.9) |
| 0.0 | 22.3 (23.5; 1.63-42.7) | 19.1 (17.4; 1.23-40.8) |
| +1.0 | 36.1 (43.3; 7.2-48.4) | 33.3 (40.9; 4.0-47.9) |
| +2.0 | 44.6 (48.6; 32.7-49.6) | 43.0 (48.1; 25.3-49.5) |
| V2 | 10.0 (5.44; 1.03-26.8) | 10.0 (4.95; 1.02-27.4) |

so the ladder runs from mean effective partners about 1.2 (very concentrated) to about 44 (nearly flat over the 50 candidates; the cap bounds it at 50), 13 points, with the V2 reference between -1.0 and -0.5. The full probe grid (log2 tau -8 to +4) is in the probe files. For every temperature and every dataset the real-data distribution of effective partners per row (mean, median, 10th and 90th percentile, min, max; rows with at least one candidate) is written to `<ds>/tau_table.csv` and assembled in `effective_partners.csv`; K562, yeast, Kemmeren and HepG2 are read on the same absolute tau grid, their effective-partner distributions differ and are recorded, not re-tuned.

Fixed here, before any score (implementation choices the request did not spell out):
1. If the permuted operator of a temperature has edge >= 1 - 1e-6 (no decay) or ln 0.05 / ln edge > 256, t is set to 256 and the point is flagged `capped` in the t table (the earlier grids' maximum depth is 256). Nothing else is altered. Very concentrated operators are expected to be disconnected or nearly so, so this flag may trigger at the low end; it is reported, not tuned.
2. Multi-scale scores use the pair scores over ALL gene pairs (upper triangle). Rank percentile = average-tie rank / number of pairs. The within-temperature z (MS2) uses the mean and sd over all pairs.
3. The hub-controlled scope uses ONE stratifier for all cells of a dataset: H = q_i + q_j with q the stationary distribution of the registered operator (keff 10, cap 50), exactly as the scoreboard, the direct-edge diagnostic and the width ladder; it is not recomputed per temperature. 10 quantile strata of H (`rz_hub_stratify`).
4. Best temperature = argmax of the metric over the 13 ladder points of a rule (first in case of an exact tie); the V2 reference is shown beside it and takes no part. It is "at an edge" when it is log2 tau = -4 or +2; "next to an edge" (-3.5, +1.5) is also listed. Depth rules are read separately.
5. A failing gate stops the run; the t table, effective-partner table, gate tables and ledger are written before any heatmap.
6. The registered per-gene operator (keff 10, cap 50) is carried as the comparators tcaCLR (t = 9) and tcaCLR2 (t = t2 of the registered edge file), plus CLR; they are not settings of the ladder.

## Multi-scale scores (no tuning, no selection), for each depth rule

| Score | Definition |
|---|---|
| MS1 | mean over the 13 ladder temperatures of each pair's rank percentile (within the temperature's score over all pairs) |
| MS2 | mean over the 13 temperatures of the pair's score standardised within the temperature (z over all pairs) |
| MS3 | as MS1 over the 8 temperatures log2 tau = -2.5, -2, ..., +1 (mean effective partners about 2 to 36 on M3D) |

Compared with each single temperature, the V2 reference, CLR, tcaCLR and tcaCLR2. Differences and ratios are shown; no interval, no p-value.

## Datasets, truths, metrics (as the width ladder; all spent; the single-cell sets are out of scope)

| Dataset | Matrix | Metrics per score |
|---|---|---|
| M3D (4297 genes) and EcoMAC (4189, ComBat) | scoreboard matrices (`sb_load_ecoli`) | pooled AUPR / base rate (raw and hub-controlled) on RegulonDB regulons SC, all and Complex Portal complexes; direct TF-to-target edges (SC, all): pooled AUPR ratio, median per-TF AUROC, hub-controlled pooled ratio (as `analysis/direct_edge_diag`) |
| K562 pseudobulk (5671 genes) | ideas47 set-up | pooled AUPR / base rate (raw and hub-controlled): CORUM, DoRothEA A-B, A-B-C, CollecTRI, KD response |
| yeast net4 (5950), Kemmeren (6169; CYC2008, Arm B), HepG2 (9622; Complex Portal, Arm A) | stage-2 matrices and frozen truth files (hash-guarded) | pooled AUPR / base rate (raw) |

All truth builders, scorers, the hub stratifier, `gl_walk`, `gl_zmat`, the edge machinery (`ad_fit`, `ad_permute_genes`, `ad_Psym`, `ad_eigs`, `ed_t2`) and the evaluator `wl_eval` are the existing ones (the evaluator, loaders and walk are sourced from `analysis/width_ladder/wl_lib.R`, which is unchanged). The new code only builds the global-temperature operator from the recipe-ladder lines and loops over temperatures.

## Gates (spent truths, gates only; tolerance 1e-9; a failure stops the run)

1. The registered per-gene operator (keff 10, cap 50) through this code, t = 9, reproduces the stored tcaCLR cells of `results/tcaclr2_scoreboard/cells_reference_tcaclr_all.csv` (every pooled raw and hub-controlled value, tf_aupr_ratio, tf_auroc_median), the direct-edge numbers of `results/direct_edge_diag/gate_direct_edge.csv`, the hub-controlled direct-edge ratios of `results/direct_edge_diag/variants_{m3d,ecomac}.csv`, and, with the stored t2, the tcaCLR2 and CLR cells (the width-ladder gate comparison, same arithmetic). The reference edge recomputed by the edge job must equal the stored edge file (edge to 1e-6, t2 identical).
2. V2's own temperature of M3D and EcoMAC (tau of `variant_info_*.csv` to 1e-9) at t = 9 reproduces `results/recipe_ladder/variants_{m3d,ecomac}.csv` stage V2_global_temp (pooled direct-edge ratio SC and all, median per-TF AUROC, hub-controlled ratio; the stored median AUROC of the recipe ladder is `tf_rows` of `analysis/direct_edge_diag/lib.R`, TFs with at least 5 targets, whereas the ladder's own metric `tf_auroc_median` is the scoreboard's `sb_per_tf_median_auroc`, TFs with at least 3 targets: gate 2 recomputes the recipe's figure with `tf_rows` on the V2 score, and the two definitions are not mixed in any table): M3D SC 7.82 and EcoMAC SC 6.66 to printed precision, and everything to 1e-9. The score job runs both gates inline on the reference and V2 operators before any ladder point is scored; `--mode gate` runs them alone. Other datasets have gate 1 only (no stored V2 value).

## Outputs (`results/global_temp_ladder/`)

Per dataset: `tau_table.csv` (real-data effective partners per temperature and V2's tau), `edge.tsv` (edge, t, capped flag per temperature), `edge_perm.csv`, `edge_gate.csv`, `gate.csv`, `cells.csv` (long: dataset, rule, score, tau, t, truth, scope, metric, value), `run.csv`, logs. Assembled: `surface_all.csv` (temperature by metric, per rule, with V2, MS and comparators), `best_tau.csv` (with edge flags), `effective_partners.csv`, `t_by_tau.csv`, `ledger.csv`, `gate_summary.csv`, heatmaps and curves (PNG).

## Ledger and multiplicity statement

Every setting scored is a row of the ledger table (dataset, rule, score, number of metric cells) written by the assemble step. Settings per dataset: 13 ladder temperatures x 2 rules = 26, the V2 reference x 2 rules = 2, and 3 multi-scale scores x 2 rules = 6, i.e. 34, plus three comparators (CLR, tcaCLR, tcaCLR2) that are not settings. Metric cells per setting: 12 on M3D and EcoMAC (6 pooled + 6 direct-edge), 10 on K562, 1 on each of yeast, Kemmeren, HepG2: 34 x (12 + 12 + 10 + 3) = 1258 setting-by-metric cells. No correction is applied because nothing is tested. The cells are correlated (same pair scores, overlapping truths), the best-temperature locations are the maximum of 13 values per cell, V2's temperature was itself chosen by looking at M3D and EcoMAC direct-edge results in the recipe ladder (a selection already seen on these two datasets), and a "best" is optimistically biased on spent truths. Nothing here may be quoted as an accuracy result for a temperature, only as the shape of the response surface; any temperature taken from it is a selection and needs a fresh test.
