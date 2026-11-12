# EXPLORATORY plan: adaptive depth by hub decile on a log grid (ideation memo, section 2, idea 3)

Written 12 Nov 2026 and committed BEFORE any compute for this experiment. Exploratory under the PI's standing rules: M3D is the tuning set, K562 is spent (CORUM and DoRothEA were used in earlier tests). No claim, no p-values, nothing here is selected for the paper or the pre-read. After this commit the plan is not changed; every deviation is logged in `docs/explore_idea3_results.md` with its timing relative to results.

Source of the idea: `docs/plans/fable_ideation.md`, section 2, idea 3 (and section 0, which records that no method with a per-gene or per-stratum depth exists in `analysis/` or `R/`). Conventions inherited from `docs/preregistration_depth_log.md` (log grid, equivalence gate, `gc()` before forking), `analysis/lab_splits.R` (the label-free held-out criterion) and `analysis/robust_z/` (operator code path, MI row mean as the method-neutral hub stratifier, hub-stratified scoring, jackknife streams).

## 1. Question

Does a depth that depends on a gene's hub degree, chosen with the label-free lab-split criterion on a log grid, improve on one global depth (t = 9), and does it improve on the controls that mix all depths with one global parameter (personalised PageRank and the heat kernel)? Is the decile-to-decile variation of the chosen depth larger than what a rewired operator produces (fitting noise)?

## 2. Operator, data, equivalence gate

- **Operator** (unchanged package path, via `analysis/robust_z/robust_z_lib.R`, V0): MI (bins "hg", transform "none"), CLR with normal calibration and Stouffer combination, softmax attention keff 10 cap 50 (package `clr:::.softmax_rows`, row normalisation and self-loop for empty rows), lazy walk P = 0.5 I + 0.5 A. Walk powers P^t by repeated multiplication, t = 1..128.
- **M3D**: 466 replicate-averaged experiments, 4,297 genes (constant rows dropped), 39 labs; truths: RegulonDB strong evidence (SC, base rate 0.0943) and all evidence (all, 0.1205), as `analysis/baselines.R` and `robust_z.R` define them (same-operon pairs removed).
- **K562** (run only if the gate in section 11 opens): essential Perturb-seq pseudobulk half A, 5,671 truth genes (`results/baselines_k562ess_truthgenes_genie3/expression.tsv`), each profile its own "lab"; truths DoRothEA A-B, CollecTRI (`collectri_sets` on the 5,671 genes), CORUM all and CORUM without the 8 RP-dominated complexes ("noRP"), constructed exactly as `robust_z.R`.
- **Equivalence gate (first, before any new quantity is reported; the job aborts on failure).**
  1. The readout function of this experiment with every row at t = 9 (S = (P^9 + P^9')/2, zero diagonal) must reproduce the pooled AUPR of `results/baselines_m3d/metrics.csv` (attention, SC and all) to 1e-9. For K562: DoRothEA A-B, CollecTRI and CORUM all from `results/baselines_k562ess_truthgenes_genie3/metrics.csv`, to 1e-9.
  2. The held-out criterion machinery (global, symmetrised, upper triangle, the `lab_splits.R` form) at t in {1, 2, 3, 4, 6, 8, 12, 16} for the first 8 lab splits x 2 directions (M3D; all 6 profile splits for K562) must reproduce the V0 curves in `results/robust_z/{m3d,k562}/heldout_curves.csv` to 1e-9.
  3. The hub-stratified enrichments at t = 9 (MI row mean and attention pi stratifiers) must reproduce the V0 att9 rows of `results/robust_z/{m3d,k562}/pooled_metrics.csv` (`enr_hub_mi`, `enr_hub_att`) to 1e-9.

## 3. Hub statistic and strata

- **Primary hub statistic**: H_i = row mean of the raw MI matrix on the full data, diagonal excluded, (rowSums(M) - diag(M)) / (G - 1). It does not depend on the attention operator (method-neutral), and uses no truth. Deciles by rank, ties broken by order of appearance (`cut(rank(H, ties.method = "first"), 10, labels = FALSE)`): equal-sized strata, decile 1 = least connected, decile 10 = most connected.
- **Sensitivity**: the same with H = pi, the stationary distribution of the full-data V0 attention walk (`rz_stationary`, teleport 1e-3). Everything in sections 4-9 is run for both stratifiers and both are reported; the MI row mean is primary.
- Strata are fixed on the full data and applied to the rows of every lab split (a gene's decile does not change between splits). The hub statistic uses all labs, including each split's held-out half; it uses no truth.

## 4. Depth grid and label-free depth per decile (form 1)

- **Grid**: t in {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128} (14 depths, about half an octave apart).
- **Splits**: M3D 16 balanced random lab splits (`lab_split()` as `lab_splits.R`; labs shuffled and assigned to the half with fewer experiments; `set.seed(20261014)`, splits generated in sequence, so the first 8 are those of `robust_z`), each used in both directions: 32 criterion curves per decile. K562: 6 random profile halves (the `robust_z` splits), both directions: 12 curves.
- **Criterion for decile d (restricted to rows in d)**: build the operator on the train half; the target is the test half's CLR (V0, upper triangle) pairs at or above the 99th percentile (top 1%, symmetric indicator T). For each depth t in the grid, score = the train row of P^t, P^t_ij, for all pairs (i, j) with i in decile d and j != i; the criterion is the average precision (`ap_fast`, the repo's tie convention) of these scores against T over the pooled pairs of the decile. No truth is used.
- **Selection**: the criterion curve of decile d is the mean over the 32 (12) split-directions at each grid depth. t_d = the argmax of the mean curve, with the tie rule: the within-1% interval is the contiguous run of grid points around the argmax whose mean criterion is at least 0.99 x the maximum; if it contains more than one grid point, t_d is the grid point nearest to the log2-centre of the run, (log2 t_lo + log2 t_hi) / 2 (equal distance: the smaller t). Reported per decile: t_d, the interval [t_lo, t_hi] and its width in octaves, log2(t_hi / t_lo) (0 when the argmax is alone), the SE of the mean curve at the maximum, and the per-split-direction argmax distribution (median, IQR) as a sensitivity that is not used for selection.
- **Trend**: sign of the Spearman correlation of decile index (1..10) with log2 t_d, reported with the number of deciles at or above 1 octave from t = 9 (that is t_d <= 4 or t_d >= 24 on the grid).
- **Readout (form 1)**: S_ij = (P^{t_i}_ij + P^{t_j}_ji) / 2 with t_i = t_d of gene i's decile, from the full-data operator; implemented by taking row i of P^{t_i} and symmetrising, zero diagonal.
- The global depth on this grid is also chosen the same way with the global (all rows, symmetrised, upper triangle) criterion, tie rule as above, and scored as a control ("global grid t_g"); t = 9 is the project's depth and is the reference for every contrast.

## 5. Second form: per-gene knee of log k_i(t)

- k_i(t) = 1 / sum_j (P^t_ij)^2 on the full-data operator, at the 14 grid depths. For gene i the curve y_i(m) = log k_i(t_m) over x_m = log2 t_m is normalised to [0, 1] on both axes (gene's own minimum and maximum over the grid) and the knee is the Kneedle point: the grid depth that maximises y_norm - x_norm (ties: the smaller t). Parameter-free.
- **Smoothing within decile (form 2)**: t_d^knee = the median over the decile's genes of log2(knee_i), snapped to the nearest grid depth in log2 (ties: the smaller t). Read-out as form 1 with t_d^knee. Sensitivity "2raw": t_i = the unsmoothed per-gene knee snapped to the grid.
- No held-out splits or labels are used for form 2.

## 6. Controls (one global parameter each, chosen by the same criterion)

Both mix all depths with fixed weights and are evaluated with the original global criterion (symmetrised S, upper triangle, `lab_splits.R`) on the same 32 (12) split-directions; walk powers t = 0..128 are accumulated, weights truncated at 128 and renormalised to sum to 1.
- **PPR**: K = sum_t (1 - a) a^t P^t, a in {0.5, 0.7, 0.8, 0.85, 0.9, 0.95} (mean depth a / (1 - a) = 1, 2.3, 4, 5.7, 9, 19).
- **Heat kernel**: K = sum_t e^{-tau} tau^t / t! P^t, tau in {2, 4, 6, 8, 12, 16, 24, 32} (mean depth tau).
- The parameter is the argmax of the mean criterion curve over the parameter grid with the same 1% tie rule (log2 of 1 / (1 - a) for PPR, log2 tau for HK). Read-out S = (K + K') / 2, zero diagonal, from the full-data operator.

## 7. Rewired-operator null (done first)

- **Rewiring**: directed degree-preserving edge swaps of the attention matrix A (non-self-loop edges): repeatedly pick two edges (i -> j) and (k -> l) and replace them by (i -> l) and (k -> j), rejecting swaps that create a self-loop or a duplicate edge; each edge keeps its weight and source row. Row weight multisets, out-degrees and in-degrees are preserved; community structure is destroyed. 10 x (number of edges) attempted swaps. Then P = 0.5 I + 0.5 A_rewired.
- **Step 0 (knee null, minutes; before any criterion curve is computed)**: on the full-data operator, the per-decile knees t_d^knee of the real operator against 50 rewirings (seeds fixed), same deciles. Statistics: SD over the 10 deciles of log2 t_d^knee, the range (octaves), and Spearman of decile index with log2 t_d^knee. Descriptive comparison of the real value against the null distribution (fraction of rewirings at or above the real SD). **Stated reading rule**: if the real between-decile SD is not above the 95th percentile of the 50 rewirings, "adaptive depth by knee is fitting row-weight and degree structure, not module structure" is reported literally.
- **Criterion null (inside the selection stage)**: for each split-direction one rewiring of the train operator (seeded per split) is walked and its decile criterion against the real test target is computed; the mean curves give t_d on the rewired operator; reported: SD over deciles of log2 t_d, range, against the real values. Also four fixed random permutations of the gene-to-decile labels (the real operator, decile sizes preserved) give four SDs of log2 t_d. **Stated reading rule**: if the real SD of log2 t_d is not larger than the rewired-operator SD, "adaptive" is reported as fitting noise.

## 8. Scoring and truths

- **Readouts scored on the full data**: t = 9 (reference); form 1 with MI deciles (primary) and with pi deciles; form 2 (MI deciles) and 2raw; global grid t_g; PPR; heat kernel; and, as a descriptive check on M3D only, form 1 with a random decile labelling (the first permutation of section 7).
- **Metrics**: pooled co-membership AUPR (`ap_fast`) on the truth's pair table, as `robust_z.R`; raw, and hub-stratified (`rz_hub_stratify`: within-stratum percentile in 10 strata of H_i + H_j) with H = MI row mean (primary stratifier) and H = pi. Paired delete-half gene jackknife, 100 replicates, seed 20261016, the stream order of `robust_z.R`, difference to t = 9 +/- 1.96 SD of the paired differences (descriptive intervals, not a test).
- **Truths**: M3D SC and all (raw and both hub-stratified); K562 (if run) DoRothEA A-B and CollecTRI (regulons), CORUM all and noRP (complexes), raw and both hub-stratified.
- Selection (sections 4-6) is written to a file before any truth is loaded by the scoring stage.

## 9. Pass / fail rules, copied literally from the memo

Memo, idea 3: "Pass: hub-stratified SC +0.005 over t = 9 and raw not below -0.005." Made operational and fixed here: differences are in pooled AUPR (absolute, not enrichment) on M3D SC, readout minus t = 9; "hub-stratified" is the hub-stratified AUPR with the primary stratifier (MI row mean) for the readout and for t = 9 (the pi stratifier is evaluated the same way and reported as a sensitivity). **PASS** for a readout iff (hub-stratified SC AUPR difference >= +0.005) AND (raw SC AUPR difference >= -0.005). Otherwise FAIL. The rule is evaluated for form 1 (primary), form 2, and, for reference, for each control. The memo's secondary requirement ("if decile depth does not beat them, it is not worth its complexity") is reported as the difference of form 1 to the better of PPR and the heat kernel on the same two quantities, with no separate pass line. The memo's collapse warning ("if decile depths come out nearly equal ... this collapses to t = 9") is reported as the number of deciles with t_d not equal to t_g.

## 10. Outputs

`results/explore_idea3/{m3d,k562}/`: knee null, held-out curves per decile and per control, selected depths with intervals, the equivalence report, scored tables, logs; rds caches are not committed. Code in `analysis/explore_idea3/`. Figure: criterion curves per decile on a log t axis (M3D and K562), primary stratifier, with the pi version and the rewired null. Document: `docs/explore_idea3_results.md`. Disk: rds caches under `results/` (small, a few hundred MB at most for K562) with a row in `docs/disk_budget.md`.

## 11. K562 gate (decided now, applied to the M3D numbers)

K562 (about 4 h, 40+ GB) is run only if EITHER (i) form 1 meets the literal pass rule of section 9 on M3D SC, OR (ii) the M3D decile-depth pattern is clearly different from t = 9, defined as BOTH: at least 3 of the 10 deciles have t_d at least 1 octave from 9 (t_d <= 4 or t_d >= 24), AND the real between-decile SD of log2 t_d exceeds the rewired-operator SD of section 7. The document states which condition opened or closed the gate.

## 12. Compute

Scheduler jobs (`clr-sched`, `jobs_pc/sched`; `# mem`, `# cpus`, `# after`): knee null (small); selection stage (workers fork per split; `gc()` before forking, `gc()` after every walk); scoring stage. Memory is declared from measurement of a small test and expected to OOM-retry if wrong; nice 19, per-job capped scopes at 1.3 x the declaration.
