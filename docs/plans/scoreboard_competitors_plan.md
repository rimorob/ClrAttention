# Competitor cells on the complex-truth rows of the raw scoreboard: plan

Written 16 Nov 2026 and committed BEFORE any scoring of the new cells. Purpose: every method of the short pre-read scoreboard (CLR, aCLR, tcaCLR, GENIE3, PIDC, GRNBoost2, |r|) scored on every complex-truth row, so that "not run" no longer decides those rows. This is descriptive scoring, not a test: no p-values, no intervals, no winner is declared. The PI authorised re-running competitors; any single job estimated at more than 2 hours of wall-clock is not started and is reported instead (section 5).

**Rule fixed now:** the scoring function, the symmetrisation, the truths and the gene universes below do not change after any new score has been seen. If a gate fails, the cell is not written and the failure is reported; the definition is not adjusted to make a gate pass.

## 1. Cells and sources

Metric of every cell: pooled AUPR / base rate over all gene pairs of the truth's universe, `ap_fast(S[cbind(i, j)], label)` divided by `mean(label)`, scope raw (the definition of the existing CLR cells in `docs/whitepapers/figures/scoreboard_data.csv`).

| row | truth | cells | score source | action |
|---|---|---|---|---|
| M3D (4,297 genes, 466 experiments) | Complex Portal complexes (`complexes`, 1,150 positive of 156,267 pairs) | GENIE3, PIDC, GRNBoost2 | saved symmetrised matrices `results/baselines_m3d/{genie3,pidc,grnboost2}_scores.rds` (GRNBoost2: legacy early stopping, seed 20261016, 27 Oct run; GENIE3 RF, 1,000 trees, K = sqrt; PIDC NetworkInference.jl Bayesian blocks) | reuse, score only |
| EcoMAC (4,189 genes, 2,198 arrays, ComBat) | Complex Portal complexes (1,148 positive of 146,902 pairs) | GENIE3, PIDC, GRNBoost2 | saved `results/baselines_ecomac/{genie3,pidc,grnboost2}_scores.rds` (GRNBoost2 legacy, 29 Oct) | reuse, score only |
| Yeast, DREAM5 network 4 expression (5,950 genes, 536 chips) | Complex Portal main truth (262 complexes, 13,438 positive of 1,282,352 pairs), frozen under `docs/complex_confirm_selections` and `results/complex_confirm/net4/truth` (read only) | GENIE3, PIDC, GRNBoost2, \|r\| | saved DREAM5 caches `results/dream5/net4/{genie3_all,grnboost2_all,pidc,abs_pearson}__replicate.scores` (all genes as regulators, replicate-level input, the same input as the CLR cell; GRNBoost2 legacy, seed 20261016; SHA-256 sidecars) | reuse, score only |
| RPE1 (8,747 genes, 1,023 pseudobulk profiles) | CORUM merged (678 complexes, 34,554 positive of 2,680,770 pairs) | GENIE3, PIDC, GRNBoost2 | none saved | **not started: each method estimated above 2 h (section 5)** |

\|r\| exists for M3D, EcoMAC and RPE1 (`results/clr_logpt_test/truth/*/metrics.csv`). The DREAM5 network-4 TF truth and any other candidate truth are not opened; the yeast Complex Portal truth is already consumed, so descriptive scoring there is allowed.

Symmetrisation of directed methods: `max(w_ij, w_ji)` with missing = 0 and zero diagonal, as `analysis/baselines.R::sym_max`; the saved E. coli matrices are already in that form, the yeast GENIE3 and GRNBoost2 matrices (5,950 x 5,950, directed, NA diagonal) are symmetrised the same way at scoring time; the yeast PIDC and |r| caches are packed symmetric upper triangles.

## 2. Code

`analysis/scoreboard/sb_cx_lib.R` (helpers; the truth builders and data loaders are the function definitions of `analysis/clr_logpt_test.R` evaluated verbatim, none of its modes or gates run; yeast truths and loaders are `analysis/complex_confirm/cc_lib.R::cc_build_truths(..., ym = FALSE)` and `analysis/dream5_lib.R::d5_read_net`, so the DREAM5 gold standard is never read) and `analysis/scoreboard/sb_cx_score.R` (gates and scoring; one job). Output `results/scoreboard/complex_competitors/` (`gates.csv`, `pooled.csv`, log). `analysis/scoreboard_matrix.R` gets a block that reads `pooled.csv` and adds the cells with status `computed_now`; the old rows must come out identical (checked by diff).

## 3. Gates (all must pass before any new cell is written; tolerance 1e-9 on AUPR unless stated)

(a) Saved competitor matrices re-scored on an existing truth reproduce the existing value: GENIE3, PIDC, GRNBoost2 on M3D regulons SC and all, EcoMAC SC and all (`results/baselines_{m3d,ecomac}/metrics.csv`), and on K562 CORUM (`results/baselines_k562ess_truthgenes_genie3/metrics.csv`, CORUM members from `regulon_members.csv`). Printed anchors: GENIE3 M3D SC 1.505 and K562 CORUM 6.32, PIDC 1.451 and 5.53.
(b) The scoring function applied to CLR's saved scores on each new row reproduces the existing CLR cell: M3D complexes 21.81, EcoMAC 22.95, RPE1 6.768 (`results/clr_logpt_test/truth/*/cache/fit.rds`), yeast 18.06 (`results/complex_confirm/net4/labelfree/cache/cc_clr.scores`); the yeast truth's `structure_md5` equals `truth_counts.csv` (5de536c4...).
(c) |r| through the same code path reproduces the existing |r| cells (M3D 29.18, EcoMAC 25.99, RPE1 2.844, computed fresh as `abs(cor(t(X)))` with NA = 0, the convention of `clr_logpt_test.R`), and on yeast the |r| recomputed by `d5_m_abs_pearson` (the function that made the cache the CLR/attention yeast pipeline uses) equals the cached `abs_pearson__replicate.scores` to 1e-12 before that cache is scored. The yeast competitor caches are also checked against their SHA-256 sidecars.

## 4. Job, safety, disk

One job `5700_sb_cx_score` through the scheduler (`# mem`, `# cpus`, nice 19, capped scope). Memory basis: the largest objects are the RPE1 CLR matrix (8,747^2 doubles, 0.61 GB, loaded for gate b) and the yeast 5,950^2 matrices (0.28 GB each, a handful live at once), so a few GB; a smoke run on M3D is measured first and the declaration is 1.4 x the measured peak. Wall time: minutes (scoring only; the EcoMAC ComBat and RPE1 loads dominate). No new score matrices are saved; outputs are CSVs and a log (< 1 MB). The scratch files made for the RPE1 timing (expression export 0.17 GB, two subsamples) sit in the session scratchpad and are deleted. `docs/disk_budget.md` gets a row.

## 5. RPE1 competitors: measured estimates, not started (2-hour rule)

Measured on the actual RPE1 matrix (8,745 genes after dropping the two `HSPA14` rows, which GENIE3 refuses; 1,023 profiles), all genes as regulators, from a small subsample of targets and extrapolated linearly in the number of targets (nice 19, 4 cores for GENIE3, 1 core for the others):

| method | measurement | estimated full run | at 40 cores |
|---|---|---|---|
| GENIE3 (RF, 1,000 trees) | 24 random targets in 525 s on 4 cores, peak 2.0 GB | about 210 core-hours | about 5.3 h |
| GRNBoost2 (legacy early stopping) | 6 random targets in 313 s on 1 core, peak 0.4 GB (trees per target median 332, noisy: M3D 494, EcoMAC 839) | about 130 core-hours (range 100-200) | about 3.2 h (2.5-5 h) |
| PIDC (Bayesian blocks, cubic in genes) | 600 genes in 106 s including start-up; the DREAM5 network-4 run (5,950 genes) took 43.6 h serial; parallel PIDC on that net ran out of memory at 72.7 GB | about 90-140 h serial; the parallel version would need well over 100 GB | not feasible as a single job |

Each exceeds the 2-hour limit, so none of the three RPE1 cells is started; they stay dashes and are reported to the lead for the PI. The RPE1 universe also contains the duplicated symbol `HSPA14` (two rows), which GENIE3 and PIDC cannot take; a run would drop both rows for the competitors (as `baselines.R` already does for K562) and the CLR cell would remain on 8,747 genes. A cheaper alternative that the lead may want to put to the PI (a different definition, not run here): all methods on the 2,316 CORUM genes only, the K562 precedent of `baselines.R --gene_subset truths`. Scaled from the measurements (GENIE3 roughly with G^1.5, GRNBoost2 with G^2, PIDC with G^3) that would take about 30 core-hours for GENIE3 (about 0.7 h at 40 cores), about 9 core-hours for GRNBoost2 (about 0.2 h), and about 2.6 h serial for PIDC (or about 20 minutes with 8 Distributed workers), but CLR and |r| would then also have to be rescored on that universe to stay like-for-like.

## 6. Outputs after scoring

New rows in `scoreboard_data.csv` (status `computed_now`, note naming the source matrices and the gates); not-run notes of the RPE1 competitor cells updated to say why; `scoreboard_crosscheck.csv` unchanged (it covers only K10 and CLR against `sota_table.md`; the pattern does not apply to these cells). Short-pre-read scoreboard figures regenerated by `analysis/fig_graphical_abstract_short_p3.R`; footers and the Figure 2 caption sentence corrected to what is now true; docx rebuilt and checked. Report: gates, wall times, disk, the new cells, the new tally over the 19 rows, rows whose winner changed.

## 7. Addendum after the run (16 Nov 2026; no definition changed)

Job `5701_sb_cx_score` (the first submission, `5700`, stopped on a coding error after all gates up to the yeast |r| gate had passed and before any cell was written; the error was a loop that dropped the method names, fixed, nothing else changed) passed all 54 gates (largest absolute difference 3.6e-14) and wrote the 13 cells; measured peak 3.5 GB (declared 10 GB), 3.4 min wall. A second PIDC timing was taken for the RPE1 estimate: 1,200 genes took 855 s (8.1 x the 600-gene time, as expected for a cubic cost), which extrapolates to about 92 h serial for 8,745 genes, in line with the 90-140 h of section 5. The three RPE1 competitor cells remain not run. `scoreboard_matrix.R` was re-run: the 670 rows are identical to before except the cells and notes named in section 6 (checked by comparing every row; `scoreboard_crosscheck.csv` is byte-identical).

## 8. Addendum 2 (16 Nov 2026, written and committed before any RPE1 competitor compute): RPE1 GENIE3, GRNBoost2 and PIDC, authorised by the PI

The PI authorised the three RPE1 competitor runs on all 8,745 genes ("multithreaded a little less than before, to avoid the memory crash"). Nothing about the scoring changes: same cell definition (pooled AUPR / base rate on the CORUM merged truth of `clr_logpt_test.R`, 34,554 positive of 2,680,770 pairs on the 8,747-gene universe of the existing CLR cell), same `score()` (pair extraction, non-finite = NA ranks last, `ap_fast`), same symmetrisation (max, missing = 0). The two `HSPA14` rows (which GENIE3 and PIDC cannot take, and which `baselines.R` already drops for K562) are removed from the competitors' input; pairs involving them get no score (NA, ranked last). Whether `HSPA14` is in any CORUM complex is checked and reported; if it is not, this is identical to scoring on the 8,745-gene truth. No definition, parameter or truth is changed after any new score exists; a failed gate means no cell is written.

**Runs, one after another (never concurrent), through the scheduler, 32 workers (57 % of 56 threads, below the 40 of the earlier measurement, per the PI), nice 19, capped scopes, declared memory 1.4 x measured on a subsample:**

1. Export of the RPE1 expression (8,745 genes x 1,023 profiles, the matrix `load_expr("rpe1")` gives, duplicated symbols dropped) to `results/scoreboard/complex_competitors/rpe1_work/` (scratch, deleted after scoring).
2. GENIE3: `analysis/baselines.R` settings exactly (R package GENIE3, RF, K = sqrt, 1,000 trees, all genes regulators, `set.seed(20261016)`; the doRNG path makes the result independent of the number of workers), estimated about 210 core-hours = about 6.6 h at 32 workers.
3. GRNBoost2: `tools/grnboost2.py` unchanged (legacy early stopping, seed 20261016, every gene a regulator), estimated about 130 core-hours (range 100-200) = about 4 h at 32 workers.
4. PIDC: see below.
5. One scoring job: all earlier gates repeated (a: saved M3D, EcoMAC and K562 matrices; b: CLR on RPE1 6.769; c: |r| 2.844) before any RPE1 competitor cell is written.

**PIDC, memory-bounded and fast, same definition.** Why the old parallel run needed 72.7 GB: `NetworkInference.jl` builds a G x G array of `NodePair` objects (a mutual information and a small specific-information vector per ordered pair, about 235 bytes each, 18 GB at 8,745 genes), and the Distributed version copies it to every worker. The serial cost is 4.5 microseconds per gene triple, mostly allocation. The new script `tools/pidc_chunked.jl` (Julia 1.6, one process, threads, no copies) keeps the package's own functions for everything that defines PIDC (`get_nodes` with Bayesian blocks, `get_joint_probabilities`, `apply_mutual_information_formula`, `apply_specific_information_formula`, the redundancy formula `sum(p_z .* min.(si1, si2))`, the proportional unique contribution, the Gamma-fit context `cdf(fit(Gamma, scores_i), s) + cdf(fit(Gamma, scores_j), s)` with the package's normal-CLR fallback when a fit fails) and changes only the bookkeeping:
- the specific information of all sources about target t is stored as one contiguous bins(t) x G matrix per target (about 4-9 GB in total instead of 18 GB and a per-pair heap object);
- the puc score of a pair {a, b} is accumulated by one thread over the third gene w in ascending order, adding the larger-index target's term first and then the smaller's, which is exactly the order in which the package's serial triple loop adds to that matrix cell, so the sums are bit-identical to the serial run (no races, no per-thread copies, no worker copies; each redundancy is computed twice instead of once, but allocation-free);
- the Gamma fit of each gene's score vector is done once per gene instead of once per pair (the same deterministic fit; the same fallback);
- the output is the packed upper triangle of weights (one double per pair; the package's sorted both-orientation edge list is not needed for scoring).
Peak memory is one process of about 15 GB (to be measured), threads 32.

PIDC gates before the RPE1 run (all must be bit-identical, tolerance 0): (i) 300- and 600-gene subsets of the RPE1 matrix: the new script against the unchanged `tools/pidc.jl`, every pair weight equal; (ii) the full M3D matrix (4,297 genes, the input of the saved `results/baselines_m3d/pidc_edges.tsv` that gave PIDC 1.451 on SC): the new script's weights against that saved edge file, every pair equal, and the re-scored AUPR equal to the stored one. If (ii) shows any difference, it is reported and the RPE1 PIDC cell is not written.

Estimate and go rule: the wall time is estimated from the timed subsets and the full M3D run (cubic in genes); if it is above about 6 h wall the PIDC run is not started and the estimate and memory plan are reported instead. The estimate is added to this section before the RPE1 PIDC run starts.

### 8a. PIDC estimate and gate results (16 Nov 2026, before the RPE1 PIDC run)

Gates (job `5712_cx_pidc_gate_final`, final `tools/pidc_chunked.jl`): the 300- and 600-gene RPE1 subsets against `tools/pidc.jl` and the full M3D matrix (4,297 genes, 9,229,956 pairs) against the saved serial package output `results/baselines_m3d/pidc_edges.tsv` are all bit-identical (0 differing pairs, max absolute difference 0.0). The first full-M3D run (before the garbage-collection blocks were added) was identical too, but peaked at 13.5 GB; with a full GC between 128-row blocks the peak is 3.0 GB at the same 138 s.
Estimate for the RPE1 run (8,745 genes, 32 threads): the M3D run took 138 s (puc stage 89 s, cubic in genes, 28 threads); scaled, about 12.5 min for the puc stage, about 1 min for reading the input, about 2 min for node pairs and about 1 min for the context: about 20 minutes wall. Memory: specific information 8,745^2 x 8.4 bins x 8 B = 5.1 GB, two G x G matrices 1.2 GB, the input and temporaries: about 10 GB estimated from the M3D peak, declared 16 GB. This is far below the 6-hour limit, so the run goes ahead.
Measured on the RPE1 matrix for GENIE3 (16 random targets, 8 workers: 180.5 s, peak 4.1 GB for the process tree, extrapolated to about 14 GB at 32 workers including the final matrix, declared 20 GB) and GRNBoost2 (600 genes, 8 workers: 212 s, peak 0.66 GB; edge-list peak scaled in section 8b).

### 8b. Outcome of the RPE1 runs (17 Nov 2026)

All 59 gates of job `5740_cx_rpe1_score` passed (the 54 of section 7 plus the five universe checks of the RPE1 competitor inputs; no CORUM positive pair involves the two dropped `HSPA14` rows, so the scoring is identical to scoring the 8,745-gene truth). Wall times (32 workers/threads): PIDC 19.3 min (peak 10.5 GB RSS, 6.7 GB by the scheduler; first gate attempt had peaked at 13.5 GB before the blocked garbage collection), GENIE3 435.6 min (the first attempt, declared 20 GB from an 8-worker extrapolation, was killed at 22.8 GB; the retry declared 35 GB peaked at 34.6 GB: forked workers grow as they finish targets, so a 16-target measurement underestimated memory about 2x), GRNBoost2 221.7 min (8,657,484 edges, peak 4.3 GB, trees per target median 237). RPE1 CORUM merged, pooled AUPR / base rate: PIDC 8.11, CLR 6.77, tcaCLR 6.14, GENIE3 5.47, GRNBoost2 3.85, aCLR 3.33, |r| 2.84.
