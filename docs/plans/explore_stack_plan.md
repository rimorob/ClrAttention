# EXPLORATORY plan: stacked a3r layers ("does a3r lend itself to multiple CLR passes?")

Written 12 Nov 2026 and committed BEFORE any compute for this experiment. Exploratory under the PI's standing rules: K562 and M3D are spent / tuning sets (a3r itself was chosen from about 20 readouts on K562 truths, so K562 is optimistic, and every stacked variant adds to that selection). No confirmatory claim, no p-values, no winner is selected, nothing here goes into the paper or the pre-read. After this commit the plan is not changed; every deviation is logged in `docs/explore_stack_results.md` with its timing relative to the first score. Question from the PI; definitions fixed by the lead and copied here.

Conventions inherited from `docs/plans/explore_idea3_plan.md` and `docs/plans/explore_ideas47_plan.md` (data, truths, scoring code, hub strata, paired delete-half jackknife, equivalence gate, jobs).

## 1. Question

Does a3r lend itself to multiple CLR passes? Taking the a3r output (a symmetric z-scale matrix) as the "CLR score" of a new attention step (softmax, walk, a3r readout) gives a stacked layer. Does a second or third layer improve on one layer, and does it improve on one layer run at the matched total depth (t = 18, 27), which is the control that says whether stacking is anything more than depth? Where do layers degrade (hub drift, collapsed rows, rows without candidates, reachability)?

## 2. Definitions (fixed by the lead; not changed after any score exists)

- **Layer-1 operator and a3r (L1, the reference).** MI (B-spline, bins "hg", transform none) -> CLR (normal calibration, Stouffer combination, truncated at 0, bilateral) -> softmax keff 10, cap 50 (`clr:::.softmax_rows`, row normalisation, self-loop for empty rows; `rz_operator`) -> lazy walk P = 0.5 I + 0.5 A -> M = P^9; readout = bilateral CLR on log P^9 with row statistics over reachable pairs only (reference code `docs/plans/corum_analysis.md`, "Reference code for the recommended readout (a3r ...)"; as `analysis/clr_logpt_test.R::a3r` and `analysis/explore_ideas47/ideas47_lib.R::a3r`). Output S1: symmetric z-scale matrix, zero diagonal, truncated at 0.
- **Stacked layer k + 1 (a3r2, a3r3, generically a3rN).** Take S_k in place of the CLR score in the operator step: candidates = up to 50 largest S_k > 0 per row, keff-10 temperature softmax, rows normalised, self-loop for rows with no positive candidate, lazy walk 0.5, power t_k, then the a3r readout, to produce S_{k+1}. Nothing else changes. **Primary: t_k = 9 at every layer** (no tuning; a transfer). N = 2 and 3.
- **Secondary descriptive grid, layer 2 only:** t_2 in {2, 4} (layer 1 stays at 9). Reported with the same columns, not rated pass/fail.
- **Control c1 (matched total reach).** a3r at ONE layer (the K10 operator) with deeper global t = 18 (matched to a3r2) and t = 27 (matched to a3r3).
- **Control c2 ("pure repeat"): a3r-style CLR applied again to S1 directly, no softmax, no walk.** Two fixed forms: **c2_z (primary)**: bilateral row-standardisation of S1 itself, i.e. `stouffer_sym(rowz_na(S1))` (row mean and SD over the off-diagonal entries of S1, zeros included because S1 is already on the z scale; negatives and non-finite z set to 0; the "a1" form applied to S1); **c2_log**: the literal a3r function applied to S1 (`a3r(S1)`: log of S1, zeros and non-finite entries treated as unreachable, row statistics over the positive entries). Both are scored; neither is selected.
- **Control c3 (hub control, as the robust_z and idea 3 runs), for every readout.** Hub-stratified AUPR with the within-stratum percentile in 10 strata of H_i + H_j (`rz_hub_stratify`): H = raw-MI row mean, (rowSums(M) - diag(M)) / (G - 1), **primary**; H = pi, the stationary distribution (teleport 1e-3) of the layer-1 lazy walk, **sensitivity** (the same H for every readout, as ideas 4 and 7).
- **Reference scores carried in every table:** attention K10 (t = 9), CLR, a3r.

Readouts scored (11 per dataset, each in three scopes): K10, CLR, a3r; a3r2 (primary), a3r3 (primary); a3r2_t2_2, a3r2_t2_4 (secondary grid); c1_t18, c1_t27; c2_z, c2_log.

## 3. Data and truths

- **K562** (first): essential Perturb-seq pseudobulk half A, `results/baselines_k562ess_truthgenes_genie3/expression.tsv`, 5,671 truth genes x 707 profiles (zero-range genes dropped). Truths: CORUM_all (= "CORUM raw"), CORUM_merged, CORUM_noRP, CORUM_merged_noRP, CORUM_RP, CORUM_small (3-6 genes), CORUM_large (>= 20), DoRothEA A-B, CollecTRI, KD response (near-duplicate caveat; reported, not in the pass rule). Constructed exactly as `analysis/explore_ideas47/ideas47_lib.R` (which follows `robust_z.R`), sourced unchanged.
- **M3D** (run regardless of the K562 outcome): 466 replicate-averaged experiments, 4,297 genes (constant rows dropped), as `analysis/explore_idea3/idea3_lib.R::i3_load` and `clr_logpt_test.R::load_expr`. Truths: RegulonDB strong evidence ("regulons_SC", base rate 0.0943), all evidence ("regulons_all", 0.1205), and the complex truth of the 26 Oct verdict ("complexes": merged Complex Portal complexes, `docs/clr_logpt_selections/ecoli_complexes_members.tsv` gated by its SHA-256, same-operon pairs removed unless the evaluability fallback applies, which it did not: 1,150 positive pairs, 54 complexes), built by the code of `clr_logpt_test.R::build_truths_ecoli` copied verbatim into the driver. No EcoMAC, DREAM5, RPE1 or other untouched truth is read.
- Truths are loaded only after the operator and layer matrices are built; no truth is used to choose anything (nothing is chosen).

## 4. Scoring

- Pooled co-membership AUPR (`ap_fast`, the repo tie convention), raw, hub MI, hub pi. Paired delete-half gene jackknife, 100 replicates, differences +/- 1.96 SD of the paired replicate differences (descriptive intervals, not a test). **K562**: `set.seed(20261016)`, `robust_z` stream order with the merged truths appended (`ideas47_lib.R`). **M3D**: `set.seed(20261022)`, halves drawn for all truths of `build_truths_ecoli` in table order, exactly as `clr_logpt_test.R::jack_dataset`, so that the stored replicates of `results/clr_logpt_test/truth/m3d/jack_replicates.csv` can be reproduced.
- Contrasts (AUPR units and enrichment units): each readout minus a3r, minus the matched-depth control (a3r2 vs c1_t18, a3r3 vs c1_t27), minus K10 and minus CLR; for every truth and scope.

## 5. Equivalence gate (first; the job aborts on failure and nothing else is scored or reported)

Tolerance as the source pipelines used: 1e-9 absolute wherever a stored per-pair or per-replicate file exists; the recorded jackknife intervals of `corum_analysis.md` are compared to 1e-9 against the full-precision record in `analysis/corum_exploratory/metrics_robust.csv` and to printed precision against the table of the document.

K562 (scored first: CLR, K10, a3r):
1. `e47_gate` of `ideas47_lib.R` unchanged: (a) CLR and K10, 7 `robust_z` truths x 3 scopes x (AUPR + 100 replicates) against `results/robust_z/k562/cache/V0__pair.rds` / `V0__att9.rds`; (b) CLR and attention AUPR against `results/baselines_k562ess_truthgenes_genie3/metrics.csv`; (c) CLR, K10, a3r AUPR on the seven `corum_expert2` truths against `results/corum_expert2/k562/t9_jackknife_EXPLORATORY.csv`, and the a3r CORUM_all enrichment printing as 7.95.
2. (d) a3r minus K10 on DoRothEA A-B, CollecTRI, CORUM (CORUM_all), KD response and CORUM noRP: point difference and the jackknife interval against `metrics_robust.csv` (variant `a_clr_bilateral_logP_t9_robust`: `d_vs_att9`, `att9_lo`, `att9_hi`), 1e-9; the CORUM line printing as +0.059 [+0.033, +0.084].
3. (e) generic stacking code with N = 1 (operator from CLR scores C through the stack function, walk to 9, a3r) equals the a3r matrix of the package path, max |diff| <= 1e-12.

M3D (scored first: CLR, K10, a3r):
1. (a) K10 AUPR for SC and all against `results/baselines_m3d/metrics.csv` (attention); hub-MI and hub-pi enrichments at t = 9 (K10) for SC and all against `results/robust_z/m3d/pooled_metrics.csv` (V0 att9), as `idea3_score.R`.
2. (b) CLR, K10 and a3r pooled AUPR on complexes, regulons_SC (and regulons_all, descriptive) against `results/clr_logpt_test/truth/m3d/metrics.csv` (`clr`, `attention_t9`, `a3r_t9`) and against `docs/clr_logpt_verdict/endpoints_all.csv` (a3r: 0.19209 complexes, 0.16222 SC), 1e-9; (c) all 100 stored jackknife replicates of those three methods x two primary truths against `jack_replicates.csv`, 1e-9. (d) the verdict's a3r minus attention differences and SEs (`endpoints_all.csv` P1, P3: d, se) reproduced to 1e-9.
3. (e) the N = 1 stack equals a3r, 1e-12.

## 6. Pass rule (literal; reported pass/fail for each layer count N = 2, 3 and for both hub stratifiers)

**K562.** Stacked readout X (a3r2 vs c1_t18, a3r3 vs c1_t27) **passes** only if ALL of:
1. X minus a3r on CORUM_all, raw: d >= +0.005 AUPR AND the paired-jackknife 95% interval excludes 0 (lower bound > 0);
2. X minus a3r on DoRothEA A-B and on CollecTRI, raw: d >= -0.005 each (point estimate);
3. X minus the matched-depth control c1 on CORUM_all, raw: d >= +0.005 AND the interval excludes 0;
4. the gain survives the hub stratifier with the same sign: X minus a3r on CORUM_all, hub-stratified, d > 0.

"Pass (MI)" evaluates condition 4 with the MI row-mean stratifier (primary); "Pass (pi)" evaluates condition 4 with the pi stratifier (sensitivity); conditions 1-3 use raw scores in both. A stricter reading is also tabulated ("strict MI", "strict pi"): conditions 1-3 evaluated on the hub-stratified scores of that stratifier (X, a3r and c1 all hub-stratified), condition 4 then being the sign of the condition-1 difference. The secondary t_2 grid and the c2 controls are not rated.

**M3D.** X passes only if RegulonDB SC and RegulonDB all are non-inferior to a3r at margin -0.005 (jackknife lower bound of X minus a3r >= -0.005 on each) AND the complex truth does not decrease (X minus a3r on "complexes" >= 0, point estimate). Evaluated on raw scores; the hub-stratified versions (MI primary, pi sensitivity) are tabulated the same way.

KD response is reported but is not part of either rule. Nothing is called a win in the results document: "passes the literal rule" is bookkeeping. If layer 2 or layer 3 fails, the document says so plainly.

## 7. Diagnostics (label-free, written per layer; layers 1-3 and, for descriptive purposes, layers 4-5 chained the same way on K562 and M3D without scoring them against any truth)

Per layer k (operator built from S_{k-1}; S_0 = the CLR score C): number and fraction of rows with no positive candidate (self-loop rows), number of rows with fewer than 10 positive candidates, median candidate count, median effective support 1/sum_j A_ij^2; for the walk at t = 9: share of unreachable off-diagonal pairs (P^9 = 0) and rows with < 2 reachable entries; hub drift: Spearman of pi_k (stationary distribution of the layer's lazy walk) with the MI row mean, Spearman of pi_k with pi_1, pi_k max / median, share of in-weight (column sums of A_k) held by the top 1% of genes; for S_k: fraction of nonzero upper-triangle entries, rows of all zero ("collapsed rows"), Spearman of the row mean of S_k with the MI row mean, and the overlap (Jaccard) of the top-1% pair sets of S_k with S_{k-1} and with the CLR score. These numbers say where things degrade; they carry no truth.

## 8. Outputs

`results/explore_stack/{k562,m3d}/`: `run.log`, `equivalence_gate.csv`, `diagnostics_layers.csv`, `pooled_levels.csv`, `contrasts.csv`, `pass_fail.csv`, `tables.md`; rds readout caches in `cache/` (not committed). Code in `analysis/explore_stack/` (`stack_lib.R`, `stack_k562.R`, `stack_m3d.R`, `stack_summarize.R`); sources the existing libraries read-only. Document: `docs/explore_stack_results.md`.

## 9. Compute (scheduler only: `clr-sched`, `jobs_pc/sched`; nice 19; per-job capped scopes at 1.3 x declared; `gc()` before forking)

Dense G x G walks are fine (P is sparse with about 51 non-zeros per row, so P^t by repeated sparse-dense multiplication is cheap); each layer builds a new softmax operator (about 1 min on 5,671 genes, dominated by the row-wise `order`) and a few dense G x G temporaries (257 MB each on K562, 148 MB on M3D). The memory is dominated by the scoring stage (10 truths, 2.5-5.3 million pairs each, 8 forked jackknife workers), as in ideas 4 and 7 whose peaks were 19.8-24.3 GB on K562.

| job | what | est. peak | declared (1.4 x) | cpus |
|---|---|---|---|---|
| `5500_stk_k562` | gate, 11 readouts x 3 scopes x 10 truths, layer diagnostics 1-5 | 24 GB | 34 GB | 8 |
| `5510_stk_m3d` | gate, 11 readouts x 3 scopes x 3 truths, layer diagnostics 1-5 | 12 GB | 17 GB | 8 |
| `5520_stk_summarize` | contrasts, pass/fail tables (after both) | 4 GB | 6 GB | 2 |

Time estimates: K562 about 3 min of set-up, about 130-170 s per readout for scoring (11 readouts, about 30 min), diagnostics layers 1-5 about 15 min, in total under 1.5 h; M3D about a third of that. An OOM is accepted (the scheduler retries once at 1.5 x) and is reported. Disk: caches a few MB per dataset (no matrices are kept), row in `docs/disk_budget.md`.

## 10. Caveats stated in advance

a3r was chosen from about 20 readouts on K562 truths; K562 is optimistic and its CORUM and DoRothEA were used in earlier tests; stacked variants add to that selection. M3D is the tuning set. The stacked a3r uses a z-scale input to a softmax whose temperature is set per row by perplexity, so the scale of S_k does not matter, but its top-50 ordering does. A positive result here would only motivate a preregistered test on untouched data, not a claim.
