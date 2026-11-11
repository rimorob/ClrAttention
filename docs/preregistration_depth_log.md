# Preregistration: depth on a logarithmic grid (log-scale re-test of "depth is not a resolution dial")

Written 2026-11-11, committed before any new compute. Analysis plan only; no new result has been computed at the time of writing.

## 1. Status of the evidence (read this first)

- This is NOT a confirmatory or blind test. The truths (E. coli M3D with RegulonDB; K562 essential compendium) and the operators were all used in the earlier linear-grid test (`analysis/granularity.R`, `results/granularity/`, Supplementary Note 12), and the linear-grid curves have been seen, including on a log2 axis (section 9). The preregistration fixes the analysis and the readings in advance; it does not make the result confirmatory.
- No truth set is added, removed or changed after the first result of the new run. No operator, grid point or statistic is changed after the first result. Anything done differently from this document is listed as a deviation in `docs/depth_log_results.md`, with the reason and whether it was decided before or after seeing results.
- Power is tiny: n = 3 truths in E. coli and n = 5 in K562. With one-sided exact permutation tests the smallest attainable P is 1/6 = 0.167 (n = 3) and 1/120 = 0.0083 (n = 5). The primary statistic is therefore descriptive. It is not a decisive test, and a P value is reported but is not part of the reading rule.

## 2. Question

The reach of a lazy random walk grows geometrically with the number of steps and relaxation times scale like a geometric grid, so depth t should be analysed on a log scale. The earlier test used a linear grid t = 1..32, and several best depths sat at the grid edge (t = 32 for E. coli TF regulons and K562 CollecTRI; 26, with the 1%-interval touching 32, for K562 KD-response modules). Question: on a geometric grid that extends to t = 256, does the depth at which a truth is best recovered increase from fine to coarse truths (log2 t_best increasing with truth coarseness), without being censored by the grid edge?

## 3. Data, truths and operators (identical to `analysis/granularity.R`)

- **E. coli**: M3D, 466 experiments (`data/E_coli_v4_Build_6`, 4,297 probes, constant rows dropped), RegulonDB extract, all evidence classes (C, S, W, ?). Truths in coarseness order: 1 operons (same-operon pairs are the positives; sets with >= 2 genes) < 2 TF regulons (single-TF, 5-500 genes, operon pairs excluded) < 3 sigmulons without sigma70 and sigmaS (>= 5 genes, operon pairs excluded).
- **K562**: Perturb-seq essential compendium (`data/perturbseq/K562_essential_normalized_bulk_01.h5ad`), network half A (same split, `set.seed(20261015)`, p_max 1e-3, as `perturbseq_bulk.R`). Truths in coarseness order: 1 CORUM complexes < 2 KD-response modules (half B) < 3 DoRothEA A-B < 4 DoRothEA A-C < 5 CollecTRI.
- **Operator** (same as before): `ClrAttention` MI estimate (bins "hg", transform "none"), calibration normal + Stouffer, soft attention with keff 10 and cap 50 ("soft10"), self-loop and lazy walk P = 0.5 I + 0.5 A (alpha = 0.5), read-out S = (P^t + t(P^t)) / 2 with zero diagonal. t = 0 is CLR itself (the calibrated CLR score matrix).
- **Score**: co-membership AUPR (`ap_fast`, same tie convention) on the same pair tables as before. Regulon coherence is also written, as before, but is not used in any reading.

## 4. Depth grid

t in {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256} plus t = 0 (CLR). 16 attention depths, spaced about half an octave apart, spanning 8 octaves. Walk powers are computed by repeated multiplication by P (P^t for all t = 1..256 in sequence); only grid points are scored. Grid points 1, 2, 4, 8, 16, 32 must reproduce the saved linear-grid values (`results/granularity/granularity_curves.csv`, AUPR, base rate and coherence) to 1e-9, and t = 0 as well; the job aborts on a failed equivalence check.

## 5. Definitions

- **t_best** (primary) = the grid point in {0} union the 16 depths with the highest AUPR. t_best = 0 means CLR beats every attention depth; it is then ranked below every positive depth (as log2 t is undefined, only the rank is used, and Spearman correlation uses ranks only). Ties among t_best values receive average ranks.
- **t_best_ge1** (sensitivity only) = the highest-AUPR grid point among t >= 1. Reported alongside the primary, never used for a reading. It matters for K562 CORUM, where CLR beat every depth in the linear grid and the argmax over t >= 1 is a depth at which attention is worse than CLR (this choice of primary convention was made after seeing that the linear-grid K562 Spearman correlation depends on it: 0.7 with CLR-best ranked lowest, 0.2 over t >= 1).
- **Truth coarseness rank** = the number in the truth name (E. coli 1-3; K562 1-5), as in `granularity.R`.
- **Within-1-percent interval**: over t >= 1, the contiguous run of grid points around t_best_ge1 whose AUPR >= 0.99 x max AUPR (over t >= 1). Width in octaves = log2(t_hi / t_lo); it is 0 when only the argmax is inside, and it is limited by the grid spacing (about 0.4-0.6 octave per step). The interval is flagged right-censored if t_hi = 256.
- **Edge censoring**: a truth is censored if its primary t_best = 256. Truths whose 1%-interval touches 256 are counted separately (descriptive).
- **Module size**: number of member genes (present in the expression matrix) of a set; median over the sets of a truth.

## 6. Primary statistic (per organism, separately)

Spearman rank correlation rho between truth coarseness rank and log2(t_best) (ranks of t_best, so t_best = 0 is lowest), with the exact permutation P value over all n! assignments (n = 3: 6; n = 5: 120), one-sided for a positive association (two-sided also reported). If t_best is identical for all truths of an organism, rho is undefined and is treated as 0.

## 7. Secondary (descriptive)

1. **Depth against module size**: OLS slope of log2 t_best on log2 median module size, over truths with t_best >= 1 (truths with t_best = 0 are excluded from the fit and listed). Reported pooled over the 8 truths of both organisms and per organism. Bootstrap CI over modules: 1,000 draws; in each draw, within each truth, modules are resampled with replacement (same number as the truth has), the pair weights are updated, t_best (same definition, grid incl. t = 0) and the median module size are recomputed, and the slope refitted. Percentile 2.5%/97.5% interval. A draw in which fewer than 3 truths (pooled) or fewer than 3 truths (per organism, n = 3 needs all three) have t_best >= 1 is dropped and counted.
   - Pair weights in a draw: each module k has multiplicity m_k (0, 1, 2, ...). A positive pair gets weight w = (sum over modules k containing the pair of m_k) / (number of modules containing the pair); positive pairs with w = 0 are dropped from that draw (neither positive nor negative); all negative pairs keep weight 1. This reproduces the original AUPR exactly when every m_k = 1 (asserted in the job). Weighted AP keeps the `ap_fast` tie convention (cumulative weights within tie blocks). The operators (expression-derived, truth-free) are held fixed across draws.
   - Seed: `set.seed(20261015)`, draws generated per truth in the order of the truth list, so that the draws are reproducible and the same for every statistic.
2. **Fraction of truths censored at the edge t = 256** (primary t_best = 256), per organism and in total; also the fraction whose 1%-interval touches 256.
3. Per-truth bootstrap percentile interval of t_best (same draws), and the distribution over draws of the per-organism Spearman rho (fraction of draws with rho >= 0.9 and with rho <= 0.5). Descriptive.
4. Sensitivity: the primary statistic recomputed with t_best_ge1 (section 5), shown next to the primary and never used for the reading.

## 8. Tertiary (descriptive)

Width of the within-1-percent interval in octaves, per truth (section 5), with its edge flag.

## 9. Disclosure of what was seen before this document

The saved linear-grid curves were re-plotted on a log2 axis before this preregistration was written (exploratory, `results/granularity_log/step1_*`). With t_best from the linear grid (CLR counted as the lowest rank when it beats all depths): E. coli best t = 4, 32, 13 (Spearman rho = 0.5); K562 best t = 0 (CLR), 26, 4, 7, 32 (rho = 0.7; 0.2 over t >= 1 only). Two of 8 best depths sit at the t = 32 edge; the 1%-interval of a third (K562 KD-response modules) touches it. Under the readings below, these linear-grid values would give "not resolution-like" (E. coli rho = 0.5 <= 0.5). The new run extends the range to t = 256 and a geometric grid; it can move the E. coli and K562 optima that were censored at 32.

## 10. Pre-specified readings (primary rho of section 6 and the censoring of section 5)

- **(a) "resolution-like on a log axis"**: rho >= 0.9 in BOTH organisms and no truth (of either organism) censored at the edge (no primary t_best = 256).
- **(b) "not resolution-like"**: rho <= 0.5 in EITHER organism.
- **(c) "inconclusive"**: otherwise (that is, neither organism has rho <= 0.5, and either an organism has 0.5 < rho < 0.9 or a truth is edge-censored).
- (a) and (b) cannot both hold. For n = 3 the attainable rho values are -1, -0.5, 0, 0.5, 1 (no ties; ties give other values), for n = 5 steps of 0.1. So "rho <= 0.5" in E. coli is met by every ordering except the perfect one.
- The reading is stated literally in the results document; no further interpretation is part of this preregistration.

## 11. Outputs and code

- `analysis/granularity_log.R` (compute: curves on the geometric grid, equivalence check, module-bootstrap draws), writing to `results/granularity_log/`: `granularity_log_curves.csv` (AUPR, base rate, coherence per truth and t), `granularity_log_best_t.csv`, bootstrap draws (`depth_log_bootstrap.csv.gz`, per draw per truth t_best and median size), `granularity_log_truth_sizes.csv`, equivalence-check report and log.
- `analysis/depth_log_analysis.R` or `.py`: the statistics of sections 6-8 exactly as written, figure (log2 axis, one panel per truth, t_best, 1%-interval, edge marked) and tables, then `docs/depth_log_results.md`. The analysis code is written after this document is committed and before any new result exists.
- The linear-grid script, its outputs and the Supplementary text are not touched.
