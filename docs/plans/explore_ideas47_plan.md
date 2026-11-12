# Plan: exploratory runs of ideation-memo ideas 4 and 7 (K562, in-sample)

*12 November 2026. Written and committed BEFORE any experiment compute. After this commit the plan is not changed; every deviation is logged in `docs/explore_ideas47_results.md`. Source memo: `docs/plans/fable_ideation.md` section 2, ideas 4 (topology-gated hybrid readout) and 7 (soft MCL inflation between walk steps). Same batch, same truths, same scoring code, same jackknife and same equivalence gate as `docs/plans/explore_ideas12_plan.md` (ideas 1 and 2), so that the four ideas are comparable.*

**Status: EXPLORATORY, IN-SAMPLE.** K562 truths are spent for confirmation. No p-values are computed or presented as confirmatory; intervals and point estimates only; no variant is selected for the paper or the pre-read; nothing is written into a paper or pre-read. Pass/fail below is a literal bookkeeping rule copied from the memo, not a claim. a3r (bilateral CLR on log P^t) was not adopted by the PI on 12 Nov as the default readout; here it is an ingredient (idea 4) and a reference column.

(Only compute done before this commit: one label-free script on the full K562 expression matrix (MI, CLR, K10 operator; no truth read, no truth file opened) that printed the distribution of four candidate definitions of the clustering coefficient c_i (quantiles, Spearman correlations between them, the density mode). It was used only to make the definition of c_i below concrete and to see that the textbook definition gives values of order 0.02, so that the product g = c_i c_j is small; that is why the percentile and soft-threshold forms are added as pre-declared variants.)

## 0. Common material (identical to the ideas 1-2 plan unless stated)

- **Data.** K562 essential Perturb-seq pseudobulk half A, `results/baselines_k562ess_truthgenes_genie3/expression.tsv`, 5,671 genes x 707 knockdown profiles (zero-range genes dropped), as `analysis/robust_z/robust_z.R`.
- **K10 and CLR (package code path).** MI (B-spline, bins "hg", transform none) -> Stouffer CLR (bilateral, truncated at 0) -> softmax keff 10, cap 50, rows normalised -> lazy walk alpha 0.5 (P = 0.5 I + 0.5 A) -> S_9 = (P^9 + P^9')/2, t = 9 fixed. "K10" = this attention at t = 9; "CLR" = the bilateral CLR matrix; "a3r" = the reference code of `docs/plans/corum_analysis.md` section 4 applied to P^9 of K10 (CORUM_all 0.10544 AUPR = 7.95x, `results/corum_expert2/k562/repro_check.csv`).
- **Truths (K562).** The eight of the ideas 1-2 plan: CORUM_all (= "CORUM raw"), CORUM_noRP, CORUM_RP, CORUM_small (3-6 genes), CORUM_large (>= 20 genes), DoRothEA A-B, KD-response (modules from half B), CollecTRI; constructed exactly as `robust_z.R` (CollecTRI as `baselines.R`). Added here (the task asks for "merged noRP" explicitly): CORUM_merged and CORUM_merged_noRP, built as `analysis/corum_expert2/corum_width_readout.R` (`merge_sets` of `analysis/clr_logpt_test.R`; raw complexes of 3-200 genes; noRP = merged sets with RP fraction < 0.3). 10 truths in total.
- **Jackknife.** Paired delete-half gene jackknife, 100 replicates, `set.seed(20261016)`, draws in the order of `robust_z.R` (DoRothEA A-B; DoRothEA A-C and CollecTRI drawn, the first discarded and the second kept as the CollecTRI stream; CORUM_all; KD_response; the 100 delete-20% draws of `robust_z`, advanced and unused; CORUM_noRP, CORUM_RP, CORUM_small, CORUM_large), then CORUM_merged and CORUM_merged_noRP appended at the end of the stream, so every replicate of the eight common truths is identical to the ideas 1-2 batch and to the stored `robust_z` V0 replicates. Interval for a difference = d +/- 1.96 SD of paired replicate differences.
- **Scopes.** Raw; hub-stratified with H = pi (stationary distribution, teleport 1e-3, of the K10 lazy walk; same for every readout); hub-stratified with H = row mean of raw MI (`robust_z`: (rowSums(M) - diag(M))/(G - 1)). Stratification = within-stratum percentile in 10 strata of H_i + H_j (`rz_hub_stratify`).
- **Metrics.** Pooled AUPR via `ap_fast` (ties as one block); enrichment = AUPR / base rate. Contrasts: each readout minus K10 and minus CLR (idea 4: also minus a3r; idea 7: also minus min(P,P') and minus the degree-scaled readout), in AUPR and enrichment units, with the jackknife interval. Per-complex (macro) AP for the raw-CORUM complexes with >= 10 pairs, by size bin 3-9 / 10-24 / 25-49 / >= 50 (point estimates, as `corum_expert2`).
- **Equivalence gate (before any new readout is computed or reported; the script stops on failure and nothing else is reported).** My driver recomputes CLR and K10 (t = 9) and must reproduce to absolute 1e-9: (a) pooled AUPR and all 100 jackknife replicates, for all three scopes, of `results/robust_z/k562/cache/V0__pair.rds` (CLR) and `V0__att9.rds` (K10) on the 7 `robust_z` truths; (b) the CLR and attention AUPR of `results/baselines_k562ess_truthgenes_genie3/metrics.csv` for DoRothEA A-B, CORUM, KD response and CollecTRI; (c) CLR, K10 and a3r AUPR on the seven `corum_expert2` truths including CORUM_merged and CORUM_merged_noRP (`results/corum_expert2/k562/t9_jackknife_EXPLORATORY.csv`), a3r to 1e-9 and its CORUM_all enrichment printed as 7.95. Idea 7 adds: (d) with r = 1 the inflation code reproduces the K10 walk S_9 to 1e-12 (max abs), and the r = 1 held-out criterion curves reproduce the V0 curves of `results/robust_z/k562/heldout_curves.csv` to 1e-9.
- **Compute and tools.** Jobs via `clr-sched` (`# mem`, `# cpus`, `# after`), nice 19, per-job capped scope; `gc()` before forking; OOM-retry accepted; declarations recalibrated from `peaks.tsv`. Code in `analysis/explore_ideas47/` (sources `analysis/robust_z/robust_z_lib.R` read-only); results in `results/explore_ideas47/` (rds caches not committed, small CSVs copied to `docs/explore_ideas47_tables/`; the disk row is appended to `docs/disk_budget.md`). A `--quick` smoke test on a gene subset is run after this commit; it reads truths in-sample and is disclosed as such.

## 1. Idea 4: topology-gated hybrid readout on a common z scale

**Gate c_i.** Graph: the K10 softmax operator A (keff 10, cap 50, row-stochastic, diagonal dropped), W = (A + A')/2, wh_ij = W_ij / max(W). N_i = the 10 genes with the largest A_ij in row i (ties by gene index). Weighted clustering coefficient of Onnela et al. (2005) of gene i's top-10 attention neighbourhood: c_i = (1 / (10 x 9)) sum over ordered pairs j != l in N_i of (wh_ij wh_il wh_jl)^(1/3). No parameter.

**Two term matrices on a common scale.** T_CLR = the bilateral CLR matrix; T_a3r = a3r of K10 P^9 (bilateral Z-score of log P^9, row statistics over reachable pairs, Stouffer combination). Each is standardised by the mean m and SD s of its upper-triangle off-diagonal entries (all G(G-1)/2 pairs): Zc = (T_CLR - m_C)/s_C, Za = (T_a3r - m_a)/s_a. Blend z_ij = g_ij Zc_ij + (1 - g_ij) Za_ij, diagonal 0. Both terms are bilateral row-standardised scores, so no rank conversion is used.

**Readouts scored (jackknife, three scopes, per-complex AP), fixed now, each scored once:**

| name | g_ij | note |
|---|---|---|
| B1 (primary) | c_i c_j | the memo's form, raw c |
| B2 | max(c_i, c_j) | the memo's "or the max" |
| B3 | s_i s_j | soft threshold: s_i = pnorm((c_i - mode(c)) / sd(c)); mode(c) = argmax of `stats::density(c)` with defaults ("soft threshold on c at the mode of its distribution") |
| B4 | max(s_i, s_j) | soft threshold, max form |
| B5 | r_i r_j | r_i = rank(c_i)/G (percentile; spread-out gate) |
| B6 (control) | 0.5 | constant blend: separates "blending two scales" from "gating by topology" |
| B7 (control) | s_i s_j with s permuted across genes | `set.seed(20261112)`; destroys the topology, keeps the g distribution |
| B8 (sensitivity) | c_i c_j | B1 with scale-only standardisation (Z = T/s, no centring) |

References: CLR, K10, a3r (all scored identically). Label-free diagnostics reported: quantiles of c, its mode, share of genes with c above the mode, Spearman of c with the MI row mean and with pi, the distribution of g for B1-B5.

**Pass/fail (copied literally from the memo).** "Pass: CORUM >= a3r (7.95x) *and* DoRothEA >= K10; fail if neither improves on a3r." Operationalised on pooled raw point estimates: CORUM = CORUM_all enrichment compared with a3r's (unrounded, 7.95x); DoRothEA = DoRothEA A-B AUPR compared with K10's. PASS: CORUM_all(B) >= CORUM_all(a3r) and DoRothEA(B) >= DoRothEA(K10). FAIL: CORUM_all(B) <= CORUM_all(a3r) and DoRothEA(B) <= DoRothEA(a3r) ("neither improves on a3r"; improvement = strictly greater). Otherwise "between". Applied to B1; reported for all of B1-B8. Intervals are shown but do not enter the rule. The memo's expectation (CORUM +0.00 to +0.02 over a3r, regulons about a3r) is reported against the table without comment.

## 2. Idea 7: soft MCL inflation between walk steps

**Operator.** P as above (lazy, alpha 0.5, K10). Walk with inflation: Q_0 = I; Q_s = rownorm( (Q_{s-1} P)^r ) (entrywise power r of the one-step expansion, then each row divided by its sum), s = 1..t. r = 1 is exactly the K10 walk Q_t = P^t. Readout S = (Q_t + Q_t')/2, diagonal 0.

**Grid.** r in {1.0 (control, = K10), 1.1, 1.2, 1.5}; t = 1..20.

**Label-free tuning (r jointly with t, the lab-split criterion of `analysis/lab_splits.R` as run by `robust_z` on K562).** Scheme: random profile halves: each knockdown profile is its own "lab" (`lab <- colnames(X)`); 6 splits made by the `robust_z` `lab_split()` (labs shuffled, then each assigned to the half with fewer profiles) under `set.seed(20261014)`; per split `set.seed(20261014 + split)` before the zero-range jitter (sd 1e-9) of the two halves; MI (bspline, hg, threads 1) and V0 CLR of each half. For train half d and test half 3 - d: the operator is built from the train-half CLR (keff 10, cap 50, lazy 0.5), walked with inflation r; criterion(r, t, split, direction) = `ap_fast` of the upper triangle of S_t against the indicator of the test half's top-1% bilateral CLR pairs (quantile 0.99 of the upper triangle). No truth is used. Selection: (i) per r: t*_r = the `robust_z` rule (per split geometric mean of the two directions' argmax t, rounded; median over the 6 splits, rounded half up); (ii) joint: the criterion averaged over 6 splits x 2 directions on the 4 x 20 grid, (r*, t*) = argmax, ties to smaller r. Curves, t*_r and (r*, t*) are written before any truth is scored.

**Readouts scored (jackknife, three scopes, per-complex AP), fixed now:** references CLR, K10 (r = 1, t = 9), K10-min = min(P^9, P^9') (diagonal 0; reciprocity, `corum_expert2` A5), K10-deg = (d_i P^9_ij + d_j P^9_ji)/2 with d_i = 1 + number of gene i's pairs in the global top-1% of CLR (`corum_expert2` degree-scaled), a3r. Inflation: IN-J = inflation at (r*, t*) (primary); IN-r{1.1, 1.2, 1.5}-own = each r at its own t*_r; IN-r{1.1, 1.2, 1.5}-t9 = each r at t = 9 (isolates the inflation from the depth choice); IN-J-min = min(Q_t, Q_t') at (r*, t*) (sensitivity: inflation plus reciprocity). Readouts that coincide with another (same r and t) are scored once and aliased in the tables.

**Pass/fail (copied literally from the memo).** "Disguise risk: high. It targets exactly what reciprocity targets, with a parameter instead of none; if it does not beat reciprocity on small complexes in-sample it should be dropped." Operationalised on pooled raw point estimates on CORUM_small (AUPR): IN-J beats reciprocity iff AUPR(IN-J) > AUPR(K10-min); the same comparison with K10-deg is reported alongside (the memo names reciprocity; the task adds the degree-scaled readout). Verdict: DROP if AUPR(IN-J) <= AUPR(K10-min); otherwise "not dropped" (flagged "beats both" or "beats min only"). The memo gives no numeric guard for DoRothEA ("as the guard"); the guard operationalised here, mirroring ideas 1-2, is dDoRothEA A-B AUPR vs K10 >= -0.005 (a breach is reported next to the verdict, not folded into it). The memo's expectation (small x1.2-1.5, large unchanged, regulons -0.005, CORUM raw may fall, satellites worse) is reported against the table without comment. The same rule is also tabulated for each other inflation readout.

## 3. Jobs (declared up front; recalibrated from `peaks.tsv`)

| job | what | mem | cpus |
|---|---|---|---|
| 5100_ei47_i4 | equivalence gate (a)-(c); caches of M, CLR, A, pi, MI row mean; c_i and diagnostics; reference + B1-B8 scoring (3 scopes, macro); idea-4 tables | 28 GB | 8 |
| 5110_ei47_i7depth | idea-7 depth stage: 6 splits x 2 directions x 4 r x 20 t, 6 forked workers; gate (d) for the r = 1 curves; t*_r and (r*, t*) written | 28 GB | 6 |
| 5120_ei47_i7score | after 5110: gate (a)-(d); idea-7 readouts at the chosen depths (3 scopes, macro); idea-7 tables | 28 GB | 8 |

## 4. Reporting

`docs/explore_ideas47_results.md`: method as run, equivalence, jobs and peak memory, all tables, literal pass/fail against the rules above, a deviation log, and a reading limited to the tables. `docs/STATUS.md` is not edited here.
