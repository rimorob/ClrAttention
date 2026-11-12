# Plan: exploratory runs of ideation-memo ideas 1 and 2 (K562, in-sample)

*12 November 2026. Written and committed BEFORE any experiment compute. After this commit the plan is not changed; every deviation is logged in `docs/explore_ideas12_results.md`. Source memo: `docs/plans/fable_ideation.md` section 2, ideas 1 and 2.*

**Status: EXPLORATORY, IN-SAMPLE.** K562 truths are spent for confirmation (CORUM, DoRothEA and the KD-response modules were all scored before). No p-values are computed or presented as confirmatory; intervals and point estimates only; no variant is selected for the paper; nothing is written into a paper or pre-read. Pass/fail below is a literal bookkeeping rule copied from the memo, not a claim.

(Only compute done before this commit: a timing test of Leiden on a synthetic random graph, and reading gene-name columns of the h5ad files to count gene overlap with the 5,671 K562 truth genes: gwps 5,494 of 5,671 genes, RPE1 4,976. No truth was opened.)

## 0. Common material

- **Data.** K562 essential Perturb-seq pseudobulk half A, `results/baselines_k562ess_truthgenes_genie3/expression.tsv` (5,671 genes x 707 knockdown profiles, genes with zero range dropped), exactly as `analysis/robust_z/robust_z.R`. Each profile is its own "lab".
- **K10 operator and references (the package code path).** MI (B-spline, bins "hg", transform none) -> Stouffer CLR (bilateral, truncated at 0) -> softmax keff 10, cap 50, rows normalised -> lazy walk alpha 0.5 -> attention score S_9 = (P^9 + P^9')/2, t = 9 fixed (no depth search). "K10" = this attention at t = 9. "CLR" = the bilateral CLR matrix.
- **Truths (as `robust_z`, K562):** CORUM_all (raw), CORUM_noRP, CORUM_RP (the 8 ribosome-dominated entries), CORUM_small (3-6 genes), CORUM_large (>= 20 genes), DoRothEA A-B, KD-response (modules from half B). Added here: CollecTRI (same universe construction as `baselines.R`; its jackknife stream is the one `robust_z` draws and discards). These are the "seven K562 truths" plus CollecTRI.
- **Scoring (copied from `robust_z.R`).** Pooled AUPR via `ap_fast` (tied scores are handled as one block, no random tie-breaking). Enrichment = AUPR / base rate. Paired delete-half gene jackknife, 100 replicates, seed 20261016, same RNG stream and order as `baselines.R`/`robust_z.R` (so replicates are identical across all variants and identical to the stored `robust_z` V0 replicates). Interval for a difference = d +/- 1.96 SD of the paired replicate differences. Three scopes: raw; hub-stratified with H = pi (stationary distribution of the **V0 K10** walk, the same for every variant); hub-stratified with H = row mean of raw MI. Stratification = within-stratum percentile in 10 strata of H_i + H_j (`rz_hub_stratify`). Contrasts: each new readout minus CLR and minus K10 (t = 9). Reported in both AUPR units and enrichment units; tables carry raw numbers and no bolding.
- **Equivalence gate (before any new readout is reported).** My driver recomputes CLR and K10 (t = 9) on the full data via the package path and must reproduce, to absolute 1e-9, (a) the pooled AUPR and all 100 jackknife replicates, for all three scopes, of `results/robust_z/k562/cache/V0__pair.rds` and `V0__att9.rds` on the 7 `robust_z` truths; (b) the CLR and attention AUPR of `results/baselines_k562ess_truthgenes_genie3/metrics.csv` for DoRothEA A-B, CORUM, KD response and CollecTRI. If any check fails, STOP and report; no new readout is computed.
- **Pass/fail units.** Pass/fail differences are in pooled raw AUPR (not enrichment), against K10 at t = 9, point estimates (intervals are shown but do not enter the rule). CORUM thresholds are in raw pooled enrichment. K10 reference values (from `robust_z`): CORUM_all 3.51x, DoRothEA A-B AUPR 0.2003, CORUM_small AUPR 0.00663; CLR: CORUM_all 6.39x.
- **Order of analyses.** (1) equivalence gate; (2) idea-2 label-free diagnostic (fraction of genes with rho < 0.5; may stop idea 2); (3) idea-1 half-split operators and Leiden sweep, with the rewired-operator null and reproducibility weights (label-free, no truth); (4) idea-1 scoring; (5) idea-2 scoring if its diagnostic did not stop it; (6) the M3D lab-split version of idea 1, only if idea 1 passes. Truths are first touched by the new readouts in steps 4 and 5; every readout listed below is scored, once, with no readout added after seeing a score.
- **Compute/tools.** Leiden: `leidenalg` 0.12.0 + `python-igraph` 1.0.0 from the existing `~/venvs/scfm` (read-only use; nothing installed). Jobs through `clr-sched` (`# mem`, `# cpus`, `# after`), nice 19; declarations recalibrated from measured peaks; `gc()` before forking; OOM-retry accepted. Code in `analysis/explore_ideas12/`; results in `results/explore_ideas12/` (rds and counts caches not committed; small CSVs copied to `docs/explore_ideas12_tables/`).

## 1. Idea 1: consensus co-membership over a log-resolution sweep of the attention walk

**Graph.** W = (A + A')/2 with A the K10 softmax operator (keff 10, cap 50, row-stochastic; diagonal dropped), weighted undirected, built (i) on each half-dataset and (ii) on the full data. (The memo's alternative, P^t at the headline t, is NOT run.)

**Partitioning.** Leiden with `RBConfigurationVertexPartition` (modularity with resolution parameter gamma), edge weights W, `n_iterations = -1` (to convergence). gamma grid {2^-6, 2^-5, ..., 2^6} (13 values). Seeds s = 1..5 (random_seed 20261200 + s).

**Halves.** 6 random 50/50 splits of the 707 profiles, produced by the `robust_z` `lab_split()` under `set.seed(20261014)` (the same 6 splits as the `robust_z` depth stage); each split gives two half-datasets (353/354 profiles); each half: MI -> CLR -> A -> W (with the 1e-9 noise fix of the depth stage for zero-range genes, seed 20261014 + split). 12 half graphs x 13 gamma x 5 seeds = 780 partitions. Plus the full-data graph x 13 gamma x 5 seeds = 65 partitions (control "FULL").

**Co-membership frequency.** freq_gamma(i~j) = fraction of runs at gamma (60 half-runs for CONS; 5 runs for FULL) in which i and j share a community. Score_ij = sum_gamma w_gamma freq_gamma(i~j), sum w = 1, diagonal 0.

**Weights (all three reported; CONS-R is the primary readout for the pass rule).**
- **R (reproducibility; label-free).** For each gamma, split and seed s, nVI between the partition of half 1 and of half 2 at seed s (VI / ln G; 30 values per gamma). Rewired-operator null: each half graph is independently rewired (igraph degree-preserving edge swaps, 10 x |E| swaps, edge weights permuted across edges; rewire seed 20261300 + 100 x split + 10 x half), partitioned with the same gamma and seeds, nVI between the two rewired halves (30 values per gamma). Excess e_gamma = max(0, mean nVI_null - mean nVI_obs); w_gamma = e_gamma / sum e. If sum e = 0, uniform. (Higher reproducibility than the rewired null earns weight; no truth is used.)
- **U (uniform):** w_gamma = 1/13. The memo flags that a gain that vanishes under uniform weights is coarse-graining in disguise.
- **G (the memo's alternative):** w_gamma = gamma / sum gamma.

**Readouts scored (jackknife, three scopes, vs CLR and vs K10):** CONS-R (primary), CONS-U, CONS-G, FULL-R and FULL-U (full-data graph, same R weights / uniform; control for what halving contributes), NULL-U and NULL-R (consensus built from the rewired-operator partitions, with uniform weights and with the observed R weights; control for degree structure alone). References: CLR, K10. Diagnostic only (point estimates, raw scope, no jackknife): CONS at each single gamma (13 readouts).

**Pass/fail (copied literally from the memo).** "Pass if CORUM raw >= 6.4x (CLR) with DoRothEA within -0.005 of K10; fail if CORUM < K10 (3.51x) or DoRothEA -0.01." Operationalised: CORUM = CORUM_all pooled raw enrichment; DoRothEA = DoRothEA A-B pooled raw AUPR minus K10's. PASS: CORUM_all >= 6.4 and dDoRothEA >= -0.005. FAIL: CORUM_all < 3.51 or dDoRothEA <= -0.01. Otherwise: neither (reported as "between"). Applied to CONS-R. A PASS with CONS-U failing the same rule is flagged "disguise risk realised" (memo). Memo falsifier (i): "if an estimator that is provably n-invariant inside exchangeable blocks still loses to CLR on ribosome-dominated CORUM under the same pooled metric, the estimand story is wrong"; reported as CORUM_all and CORUM_RP of CONS-R minus CLR.

**Conditional second stage: M3D lab-split version (only if CONS-R PASSES).** Same construction on E. coli M3D (4,297 genes x 466 experiments, 39 labs), 6 lab splits with the `robust_z` `lab_split()` (seed 20261014), same grid, seeds, weights and nulls; truths RegulonDB strong-evidence (SC) and all-evidence regulons, scored as `robust_z` m3d (same jackknife stream); equivalence gate against `results/robust_z/m3d/cache` V0 to 1e-9 first. Readouts CONS-R, CONS-U, CONS-G, NULL-U; contrasts vs CLR and K10 (t = 9), raw and both hub-stratified. The memo gives no pass rule for this stage; it is reported descriptively. M3D is the tuning set: exploratory only. If CONS-R does not pass, this stage is not run.

## 2. Idea 2: cross-compendium partner-less gate

**Detector.** Compendium X = K562 essential half A (CLR on the 5,671 genes, as above). Replicate compendium Y = K562 genome-wide (gwps) half A, 1,553 profiles, defined exactly as `analysis/perturbseq_bulk.R` (energy-test p < 1e-3, non-targeting excluded, `set.seed(20261015)` shuffle, first half; genes missing in more than 1% of profiles dropped, remaining non-finite set to 0, zero-range genes dropped). Universe U = genes present in both (5,494 expected). CLR in Y is computed on the genes of U only (MI -> Stouffer CLR, as the package); X's CLR is the full 5,671-gene matrix with columns restricted to U. For gene i: T_X(i), T_Y(i) = its top 20 positive-CLR partners (columns in U, i excluded; fewer if fewer positives; deterministic order). Overlap o_i = |T_X ∩ T_Y|. **Chance rule (fixed now):** p_i = P(Hypergeometric(population |U| - 1, successes |T_X|, draws |T_Y|) >= o_i); soft weight **rho_i = min(1, ln p_i / ln 0.05)** (rho = 1 when p <= 0.05, rho = 0 when p = 1; so rho < 0.5 means p > 0.224, i.e. in practice overlap 0). Genes outside U get rho = 1 (no information to silence), counted. k = 20 and the 0.05 rule are fixed; nothing is tuned.

**Diagnostic first (label-free, before any truth).** Fraction of genes with rho < 0.5. Operationalisation of the memo's text ("near 10.1% ... if it is 40% the idea dies here"): fraction <= 0.20: proceed; >= 0.30: **stop, idea 2 is not scored, report the fraction**; between 0.20 and 0.30: proceed, flagged "not near 10.1%". Also reported: fraction with p <= 0.05 (hard gate), fraction with overlap 0, number outside U.

**Operator.** B_ij = A_ij rho_i rho_j (i != j), B_ii = 1 - sum_{j != i} B_ij, A the K10 operator of section 0 (S1's form); then lazy walk alpha 0.5, t = 9, S_9 = (P^9 + P^9')/2, exactly as K10.

**Readouts scored (jackknife, three scopes, vs CLR and vs K10), if not stopped:** G-SOFT (rho as above; primary), G-HARD (rho = 1[p <= 0.05]), HALF-SOFT (half-split detector: the same top-20 hypergeometric rho between the two halves of each of the 6 profile splits of idea 1, averaged over the 6 splits; tests the memo's claim that the half-split detector loses power), R-SOFT (RPE1 half A as Y, same rules, defined as `perturbseq_bulk.R` on the RPE1 file; universe = genes in both; tests the "shared global programs" failure mode with a different cell line), PERM-1..5 (G-SOFT rho randomly permuted across genes, seeds 20261401..20261405; control for silencing a random set of genes with the same rho distribution). Diagnostic fractions are reported for all four detectors.

**Pass/fail (copied literally from the memo).** "Pass: small complexes (3-6) >= +0.005 and DoRothEA >= -0.005." Operationalised on G-SOFT: CORUM_small pooled raw AUPR minus K10's >= +0.005 and DoRothEA A-B pooled raw AUPR minus K10's >= -0.005 (point estimates, AUPR units; note +0.005 is about x1.75 of K10's 0.0066 on that truth). Otherwise FAIL. Memo falsifier (iii) cannot be evaluated (no synthetic comparison here); the memo's expected effect (small complexes x1.3-1.8, ribosome unchanged) is reported against the table without comment.

## 3. Jobs (declared up front; recalibrate from `peaks.tsv`)

| job | what | mem | cpus |
|---|---|---|---|
| 5000_ei12_gate | equivalence gate; caches of M, CLR, K10 operator, pi, MI row mean; CLR and K10 reference scoring | 28 GB | 8 |
| 5001_ei12_i2diag | idea-2 detectors (gwps, RPE1 CLR), rho, diagnostic, stop/go flag | 16 GB | 8 |
| 5002_ei12_halves | 12 half-dataset operators, edge lists, top-20 lists for HALF-SOFT | 14 GB | 8 |
| 5003_ei12_leiden | Leiden sweep, nVI, rewired null, count matrices (python) | 6 GB | 8 |
| 5004_ei12_score1 | idea-1 scoring and summary | 28 GB | 8 |
| 5005_ei12_score2 | idea-2 scoring and summary (exits cleanly if stopped) | 28 GB | 8 |

## 4. Reporting

`docs/explore_ideas12_results.md`: method as run, equivalence, jobs and peak memory, all tables, literal pass/fail against the rules above, a deviation log, and a reading limited to the tables. `docs/STATUS.md` is not edited here.
