# Planted-module simulation: implementation choices not fixed by the memo

Implements `docs/plans/fable_theory_memo.md` section 5.4. Written and committed before any code ran. The twelve predictions R1 to R12 are quoted from the memo unchanged; this file only fixes how they are scored. Nothing here may be revised after results are seen; any later deviation goes in the results file as a labelled deviation, never as an edit here. Synthetic data only; no project truth file or real dataset is read.

## 1. Seeds and layout

- Seeds s = 1..10. The gene layout and all loadings depend on s only (`set.seed(20261117 + s)`), so the 12 cells of one seed share one planted structure. The data draw uses `set.seed(20262000 + 100*s + 10*Nidx + noiseidx)` (Nidx 1..3, noiseidx 1 Gaussian, 2 Student-t); the identity and copula marginal arms share the same latent draw (paired).
- Gene identities are a random permutation of 1..3,000, so the index carries no information.
- Round-robin loadings: module k of a size takes lambda = c(0.55, 0.70, 0.85)[((k - 1 + s - 1) %% 3) + 1]. The offset s-1 rotates the assignment across seeds, because sizes 40 (5 modules) and 100 (2 modules) cannot contain every lambda in one draw; over seeds every size sees every lambda.
- Non-member genes (1,790) are drawn without overlap into classes in this order: 60 satellites (3 per module of size >= 20, 20 modules), 50 hubs, 30 regulators; the rest (1,650) are the target pool. Hubs, satellites, regulators and targets are disjoint.
- Regulator r = 1..30 has m = c(5, 15, 40)[((r - 1) %% 3) + 1] targets and b = c(0.35, 0.50)[(((r - 1) %/% 3) %% 2) + 1], so each (m, b) combination has 5 regulators. Targets are drawn uniformly from the pool genes currently in fewer than two regulons. Regulator loads 0.6 on its own activity factor; each target loads b on it.
- Hubs: each loads 0.25 on exactly 33 factors (30% of the 80 module + 30 regulon = 110 factors) drawn at random. Satellites load 0.8 lambda on their module factor.
- Three global factors, every gene loads N(0, 0.15^2) on each (drawn once per seed in the structure stage).
- Unit variance: x_g = sum_f L_gf F_f + sqrt(max(0.05, 1 - sum L_gf^2)) e_g, then divided by its theoretical SD sqrt(sum L^2 + residual variance). For all genes except hubs this SD is 1; hubs (sum L^2 about 2) are divided down, which scales their loadings too. The "true hub loading sum" of gene g is h_g = sum_f |L_gf / sd_g| over all 113 factors (continuous because of the global loadings, so there are no ties between strata).
- Student-t arm: only the unique noise e_g is t with 4 df divided by sqrt(2); factors stay Gaussian.
- Copula arm: x'_g = exp(sigma_g x_g), sigma_g ~ U(0.3, 1.0) drawn in the structure stage (same for all cells of a seed). Applied to the standardised value of each gene, a monotone transform (a Gaussian copula in the Gaussian-noise arm).
- Factors F ~ N(0,1) iid.

## 2. Methods (package entry points, nothing tuned)

- |r|: `abs(cor(t(X)))`.
- CLR, aCLR, tcaCLR: `ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L)`; A = `$operator`; P = 0.5 I + 0.5 A (`rz_lazy_P`); P^t by repeated sparse-by-dense products. aCLR = (P^9 + P^9')/2 with zero diagonal (the project's K10 readout). tcaCLR = the reference code of `docs/plans/corum_analysis.md` section 4 (a3r) at t = 9: log P^t, unreachable (non-finite) entries NA, diagonal NA, row mean and SD (n - 1) over reachable entries, negatives and NA to 0, (Z + Z')/sqrt(2), zero diagonal.
- Wide-row arm: the same object with `build_operator(softmax_keff = 40, softmax_cap = 200L)`; aCLR_k40 and tcaCLR_k40; CLR is the same matrix as the main arm.
- Diagnostics, all at t = 9 unless stated:
  - `tcaclr_t1`: the tcaCLR code at t = 1; identity check against the direct computation (Stouffer of the row z of the CLR scores S over each row's candidate set C_i = {j: A_ij > 0, j != i}, mean and SD (n - 1) over C_i, negatives to 0, non-candidates 0). The comparison is on the one-directional Z matrix over rows with more than 10 candidates and on the final matrix over pairs both of whose rows have more than 10 candidates; the reported number is the maximum absolute difference (R6, tolerance 1e-10).
  - `sink_raw`: Sinkhorn-Knopp balancing of the lazy operator P (not of A: A lacks total support; the lazy diagonal supplies it; balancing P alone keeps the walk lazy) by alternating row and column scaling to max |sum - 1| < 1e-9 (cap 2,000 iterations, row scaling last; the residual is recorded); P_b = D_r P D_c; raw readout = (P_b^9 + P_b^9')/2 with zero diagonal. `sink_tca` = tcaCLR code on P_b^9.
  - `pi_lift`: the project's "CLR on log lift" (`analysis/corum_exploratory/corum_followup.R`): L = log(max(P^9, 1e-300)); Lf = L - log pi_j with pi from `rz_stationary` (teleport 1e-3); symmetrise (Lf + Lf')/2; row z with diagonal excluded, sample SD with the project's denominator ncol - 2, negatives to 0, Stouffer. `colmean_lift`: identical with log(colMeans(P^9))_j in place of log pi_j.
  - `infl_1.5`, `infl_2`: linear-domain CLR (same row z / clip / Stouffer, diagonal excluded) of the elementwise power (P^9)^r, r in {1.5, 2}. Row normalisation of MCL and the row geometric mean are row scalings, to which the row z is invariant, so this is the fixed-exponent control of memo section 3.2 (exp(Z) is (M/g_i)^(1/sigma_i); here the exponent is a constant r in place of 1/sigma_i).
  - Depth grid t in {3, 6, 9, 18, 36, 72}: aCLR and tcaCLR at each t from one incremental walk; scored on T1, hub share and sigma_bg only (descriptive).
  - `hubscore`: reference only, not a method: pi_i + pi_j from the main operator's stationary distribution, scored on T1 (R8's "hub-score level").

## 3. Truths, negatives, AUPR

- Evaluated pairs: all 4,498,500 pairs of the 3,000 genes (upper triangle). AUPR is the project's `ap_fast` (sourced from `analysis/regulondb.R`, not rewritten): NA scores ranked last, radix order, ties resolved at the block end (tied positives get the precision at the end of the tie block; this penalises methods with many exact zeros equally, as everywhere else in the project). "Over base rate" = AUPR / (positives / evaluated pairs).
- T1: positives = pairs in the same module (24,345). Negatives = all other pairs.
- T2, size bin b in {5, 10, 20, 40, 70, 100}: positives = same-module pairs in modules of size b. Negatives = all pairs except positives of the other size bins (those are not negatives of a size-bin truth: they are positives of the same family). This is a choice the memo leaves open.
- T3: positives = pairs both in {regulator} union {targets} of some regulon. Negatives = all other pairs (module pairs stay negatives, as the memo says). Stratified T3 for R3: by (b, m) of the regulon, with pairs positive in a different stratum (and not in this one) excluded from the stratum's universe.
- T4: positives = the 600 regulator-target pairs (counting a pair once). Negatives = all other pairs (including target-target pairs).
- Satellite-member pairs are negatives in all truths. Reported separately: the satellite false-positive share = the fraction of the top 1% of all pairs (45k, by score; ties at the cut broken by radix order) that are satellite - own-module member pairs (2,430 exist), and the probability that a satellite - own-module member pair outranks a member-member pair of the same module (ties 0.5), pooled per module size.
- Macro per-module AP: for each module, positives = its within-module pairs, negatives = its member-to-outside pairs (i in module, j not in module); AP by `ap_fast`; mean over modules, by size.
- Hub-stratified AUPR: the project's `rz_hub_stratify` rule (10 equal-count strata of H over the upper triangle, within-stratum percentile of the score), vector implementation (gated equal to `rz_hub_stratify` on a toy matrix). Two stratifiers: H = h_i + h_j (true hub loading sum) and H = rowmean(MI)_i + rowmean(MI)_j (diagonal excluded). Scored on T1 and T3.
- Hub share of the top 1%: fraction of the top 1% of pairs that contain at least one hub gene. Primary hub set: the 30 genes (1%) with the highest stationary probability under the main operator (`rz_stationary`); the same set for every method. Secondary: the 50 planted hubs. Chance level for the primary set is 2.0%.

## 4. Mechanism curves

- Per module, the mean score over its within-module pairs. Slopes: module-level regression of log(mean score) (aCLR, CLR) or of the mean score in z-units (tcaCLR) on log n with a lambda fixed effect (one point per module); primary fit on modules with n >= 20, all-size fit reported as secondary. Slope = mean over the seed's fit, then across seeds.
- sigma_bg per row: SD (n - 1) of log P^t over reachable off-diagonal entries; reported as the median across all rows, coefficient of variation across rows, and median for member rows by size, at t = 9 (and at the depth grid for the main operator).
- Retention: m_t(n) = mean over members i of module size n of sum_{j in module} P^t_ij (self included) at t = 1 and t = 9; plus the in-module share of A (sum over j in module, j != i, of A_ij).
- Unreachable share of P^9 pairs, number of rows with at most 10 candidates, per compendium.

## 5. Scoring rules for R1 to R12 (fixed now; the prediction text is unchanged)

Primary arm = Gaussian noise, identity marginals, evaluated separately at N = 150, 400, 800 (seed means over 10 seeds). A clause is "true at N" or not; the count of true (clause, N) cells decides: held = all true, did not hold = none true, partly held = otherwise. Reported with the cell counts.

- A > B on a metric: seed-mean difference > 0 and the difference positive in at least 9 of 10 seeds. A >= B fails only if B > A holds by that rule. A ~ B: |seed-mean difference| <= max(SD across seeds of the difference, 0.10 x the larger seed mean). "Intermediate": strictly between the seed means of the two neighbours.
- "Growing with n" (R1): over the four bins 20, 40, 70, 100, the seed-mean gap is larger at 100 than at 20 with at most one inversion among the 6 pairs of bins (Kendall).
- Slopes (R1): "about -1" = in [-1.25, -0.75]; "about -0.5" = in [-0.75, -0.25]; "|slope| < 0.3" as written; the n >= 20 module-level fit.
- T2 comparisons use raw per-bin AUPR (the enrichment is a monotone rescale within a bin; across bins raw AUPR is confounded with the rising base rate, which is stated wherever R2's "below their n = 20 level" is judged).
- R1 clauses: tcaCLR > CLR and CLR > aCLR on T2 in each bin 20, 40, 70, 100 (8 tests); aCLR deficit (CLR - aCLR) growing; tcaCLR - CLR gap growing; three slope conditions. The memo's failure test is also reported: tcaCLR slope as steep as CLR's (|slope_tca| >= |slope_CLR|) or no n >= 20 advantage over aCLR.
- R2 (bin 5): CLR >= tcaCLR; tcaCLR ~ aCLR; each of the three methods' bin-5 AUPR below its bin-20 AUPR; bin 10 intermediate between bins 5 and 20 for each method.
- R3: T3 pooled aCLR ~ tcaCLR; tcaCLR > CLR and aCLR > CLR (pooled); the (b = 0.35, m = 40) stratum has the largest gap (mean of aCLR and tcaCLR minus CLR, AUPR) of the six strata.
- R4: T4 pooled: CLR >= tcaCLR; tcaCLR >= aCLR; |r| competitive = |r| AUPR >= 0.8 x the best of CLR, aCLR, tcaCLR.
- R5 (heavy tail): seed SD of the CLR n >= 20 slope larger in the Student-t arm than in the Gaussian arm (identity marginals); the mean tcaCLR - CLR gap on T2 bins >= 20 larger in the t arm than in the Gaussian arm. (copula): the sign of tcaCLR - CLR on T1 equals the identity arm's sign; |r| T1 AUPR is below CLR's and below tcaCLR's. Judged at each N (Gaussian noise for the copula clauses).
- R6: maximum absolute difference <= 1e-10 in every one of the 120 compendia (keff 10 arm, rows with more than 10 candidates).
- R7: sink_raw ~ aCLR on T2 bins 20, 40, 70, 100; primary hub share of the top 1% of sink_raw < that of aCLR (9 of 10 seeds).
- R8: for aCLR and tcaCLR: T1 AUPR at t = 36 and at t = 72 each below that at t = 9; AUPR(72) closer to the `hubscore` T1 AUPR than AUPR(9) is (on the log scale); hub share at 36 and 72 above that at 9; the rise (share(72) - share(9)) of aCLR larger than that of tcaCLR.
- R9 (n >= 20, per size 20, 40, 70, 100): satellite-outranks-member probability above 0.5 for aCLR and below 0.5 for CLR and for tcaCLR (3 conditions per size and N).
- R10: median sigma_bg at k40 below that at k10 (tcaCLR operator); tcaCLR - aCLR on T1 smaller in the k40 arm than in the k10 arm. CLR is the same matrix by construction (reported, trivially true).
- R11: pi_lift has the larger relative SD (SD / mean across seeds) of T1 AUPR than each of CLR, aCLR, tcaCLR; pi_lift's primary hub share of the top 1% exceeds that of aCLR and of tcaCLR; colmean_lift's seed-mean T1 AUPR is strictly between aCLR's and tcaCLR's.
- R12: infl_1.5 and infl_2 are each not above tcaCLR on T1 (tcaCLR >= by the 9 of 10 rule counted as "does not reach" when tcaCLR > inflation holds); sink_tca >= tcaCLR on hub-stratified T1 under the true hub loading sum stratifier (primary) and under the MI row mean stratifier (reported).
- H2 bin-edge rule (primary arm only, Gaussian, identity, per N): for each bin n in {5, 10, 20, 40, 70, 100} count seeds with T2 AUPR(tcaCLR) - AUPR(aCLR) > 0 (strictly). Edge(N) = the smallest n with count >= 9. Reported also: whether every larger bin also has count >= 9 (contiguous). Overall = the largest of the three Edge(N) (the smallest size that qualifies at every N), with "not reached" reported where a N has no qualifying bin; the smallest of the three is reported beside it.

## 6. Compute and files

- One job per cell (12 cells, 10 seeds each, results cached per compendium in `results/sim_planted/units/`), threads 4 and a worker count of 4 inside scoring; at most four cell jobs run at once (four `# after` chains), 16 cpus in total. Memory declared at 1.4 times the measured pilot peak. A completed compendium is never recomputed; a job that dies resumes at the first missing compendium.
- Aggregation (`analysis/sim_planted/sim_summarize.R`) reads only the unit files and writes `docs/plans/sim_planted_results.md`.
- Code: `analysis/sim_planted/sim_lib.R` (generator, methods, metrics), `sim_run.R` (one cell), `sim_summarize.R`. Nothing under `R/` or `src/` is edited.
