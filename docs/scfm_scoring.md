# Learned-attention comparator: preregistered scoring (2026-11-02)

Preregistration: `docs/preregistration_scfm.md` (primary fixed 2026-10-21; change log E1 to E15). Extraction: `docs/scfm_step3.md`. This document reports the scoring exactly as preregistered. Anything the preregistration does not define is labelled **exploratory**. Verdict wording follows the preregistration's decision table; it is not editorialised.

Scoring was released by the PI on 2026-11-02 ("GO"), together with permission to read CollecTRI only to count targets per TF for the K562 gene set (E14).

## 1. What was run

| Item | Value |
|---|---|
| Frozen inputs | The 239 step-3 files, re-verified before scoring: manifest SHA-256 `e27a8422...3f15` equals the E13 value and every file hash equals the manifest (0 mismatches; `analysis/scfm_score_export.py`). The 75 K562 files, SHA-256 of manifest `787cbf87...080b` (E15), likewise (0 mismatches). |
| Smoke test | Before any real truth was opened: scoring and statistics code run on synthetic truth planted from the networks (job 1011, `results/scfm_score_smoke/`): planted CLR-NS signals gave d_i > 0, planted R0 signals (orientation Aᵀ) gave d_i < 0 and Aᵀ beat A, random truths gave d_i near 0. The exact interval and p reproduce `wilcox.test(conf.int = TRUE, exact = TRUE)`. Code committed (ee9de39) before the real run. |
| BEELINE scoring | Run once on the real truths (job 1020, 11:23 to 11:46, commit 08fa10b). Scripts: `analysis/scfm_score_export.py`, `scfm_score.R` (metrics), `scfm_score_stats.R` (statistics). Outputs: `results/scfm_score/`. |
| Gate | The |r|, CLR, CLR-NS, GENIE3, GRNBoost2 and PIDC rows that `scfm_score.R` recomputes on the full gene sets reproduce `results/beeline_baselines/metrics.csv` exactly (432 rows, maximum absolute difference 0 in AUPRC ratio and EPR). The scoring function is the tie-aware `score_truth` of `analysis/beeline_baselines.R`, evaluated from that file. |
| K562 | Gene-set definition committed before extraction (E14, 916f0d1); extraction frozen before any K562 truth was opened (E15, 8f85558); scored once (job 1023, `analysis/scfm_score_k562.R`, `results/scfm_score_k562/`). |
| Repeat of the statistics only | After the real run a bug in `scfm_score_stats.R` was found (see section 9, item 1); the statistics stage was rerun on the unchanged metric table. No truth was reopened. |

## 2. Preregistered primary: CLR-NS vs scPRINT medium-v1.5 R0

d_i = for dataset i, the mean over gene sets and primary truths of log2(AUPRC ratio of CLR-NS / AUPRC ratio of scPRINT R0). Arm A, orientation S = Aᵀ (score of TF to gene g = A[g, TF]), each method on the intersection of the genes scPRINT covers (CLR-NS recomputed on it). Higher values favour CLR-NS. Margin delta = 0.15. 44 cells (dataset x gene set x truth) enter.

| Dataset | cells | AUPRC ratio CLR-NS (mean) | AUPRC ratio scPRINT R0 (mean) | d_i |
|---|---|---|---|---|
| hESC | 6 | 1.26 | 1.18 | +0.102 |
| hHep | 6 | 1.44 | 0.98 | +0.484 |
| mDC | 6 | 1.57 | 1.05 | +0.583 |
| mESC | 8 | 1.49 | 1.01 | +0.498 |
| mHSC-E | 6 | 4.57 | 1.11 | +1.524 |
| mHSC-GM | 6 | 7.44 | 1.23 | +1.897 |
| mHSC-L | 6 | 12.63 | 1.10 | +2.287 |

(Mean AUPRC ratios are plain means over the dataset's cells; the other two ratio-based methods are in `results/scfm_score/key_metrics_per_case.csv`: CLR 1.25, 1.40, 1.40, 1.43, 4.25, 6.92, 12.48; |r| 1.13, 1.35, 1.05, 1.19, 2.87, 6.56, 9.49 in the same dataset order.)

| Quantity | Value |
|---|---|
| Hodges-Lehmann estimate | **+1.032** log2 units |
| Exact Wilcoxon interval [L, U] (4th and 25th of 28 Walsh averages; coverage 92.2%) | **[+0.343, +1.897]** |
| Primary p (exact one-sided signed-rank, H0: median <= -0.15) | 0.00781 (the smallest attainable at n = 7, all seven d_i + 0.15 positive) |

**Decision-rule rows that apply (preregistered table).**

- L = +0.343 > -delta = -0.15: **CLR-NS is non-inferior to learned attention** (row 1).
- L = +0.343 > 0: **CLR-NS is superior** (row 2, tested next in the fixed sequence after NI held).
- Row 3 (tie, U < delta) does not apply: U = +1.897 > 0.15. Rows 4 and 5 (learned attention wins; U < 0, U < -delta) do not apply. The interval is not inconclusive.

**Robust NI (E11).** scPRINT takes, per dataset, the best of {R0, R0', R1, R2} x {S = Aᵀ, A, max(A, Aᵀ)} x {arm A, arm B where available}, chosen on the truths as the configuration with the smallest d_i. Best configurations: hESC arm B R2 (Aᵀ), hHep R1 (A), mDC R0' (A), mESC, mHSC-E, mHSC-GM and mHSC-L R1 (A); d_i = +0.092, +0.251, +0.375, +0.344, +1.219, +1.719, +2.093. HL +0.851, interval [+0.233, +1.656], p = 0.00781. **L > -delta: NI is robust.** (Arm-A-only version of the maximum: [+0.238, +1.656]; the original preregistered set {R0, R1, R2} x 3 orientations, arm A: [+0.251, +1.656].)

**Key secondaries (Holm at 0.05 over two tests; confirmatory because primary NI holds).**

| Test | d_i (hESC, hHep, mDC, mESC, mHSC-E, mHSC-GM, mHSC-L) | HL | interval | p (NI) | Holm |
|---|---|---|---|---|---|
| S1: CLR vs scPRINT R0, AUPRC ratio | +0.087, +0.463, +0.420, +0.471, +1.439, +1.799, +2.260 | +0.953 | [+0.279, +1.799] | 0.00781 | rejected (0.00781 <= 0.025) |
| S2: CLR-NS vs scPRINT R0, EPR | +1.495, +1.469, +2.263, +1.034, +3.160, +3.226, +13.093 | +2.337 | [+1.469, +7.678] | 0.00781 | rejected (0.00781 <= 0.05 after S1) |

Both key secondaries show non-inferiority (L > -delta) and the interval lies above 0. S2: in two cells (mHSC-L, STRING, both gene sets) scPRINT's EPR is exactly 0 and is floored at 1e-9 as in `beeline_baselines.R`, so mHSC-L's d_i (+13.09) and the upper end U = +7.68 depend on the floor; L = +1.469 does not (it is the second smallest d_i).

**Mixed model (descriptive sensitivity, preregistered; `nlme::lme` in place of `lme4`, see section 9).** log2 AUPRC ratio ~ method + N + (1 | dataset) + (1 | dataset:truth): CLR-NS minus scPRINT = +1.028, 95% CI [+0.678, +1.379].

## 3. Where the primary difference comes from (exploratory)

Exploratory, not preregistered (`analysis/scfm_score_exploratory.R`, `results/scfm_score/exploratory.txt`). Per primary truth, n = 7 datasets, HL and interval of d_i:

| Truth | HL | interval | cells where scPRINT R0 > CLR-NS |
|---|---|---|---|
| celltype_ChIP (14 cells) | -0.108 | [-0.357, +0.248] | 11 of 14 |
| nonspecific_ChIP (14 cells) | +1.066 | [+0.449, +1.732] | 0 of 14 |
| STRING (14 cells) | +2.192 | [+0.735, +3.963] | 0 of 14 |
| LOFGOF (mESC, 2 cells) | +0.484 | n = 1 | 0 of 2 |

So the primary difference is carried by the non-specific ChIP-seq and STRING truths; on the cell-type ChIP-seq truth scPRINT R0 scores higher than CLR-NS in 11 of 14 cells (largest in hESC: AUPRC ratio 1.57 against 1.03), and the interval for that truth includes 0. Other exploratory checks: taking mHSC-E from TFs+1000 only (the mHSC-E_500 organism-head flag) changes d_mHSC-E from +1.524 to +1.533 and leaves the interval unchanged; leaving out any one dataset (n = 6, coverage 90.6%) gives L between +0.300 and +0.498, all > -delta.

## 4. Descriptive readouts, orientations and arms (no claims)

Each row: CLR-NS minus the scPRINT readout in log2 AUPRC ratio (positive favours CLR-NS), HL and the same exact interval (n = 7, arm A). No p-values and no decision label are attached to these rows.

| Contrast | n | HL estimate and 92.2% interval |
|---|---|---|
| R0, S = Aᵀ | 7 | +1.03 [+0.34, +1.90] |
| R0, S = A | 7 | +1.10 [+0.35, +1.96] |
| R0, S = max(A, Aᵀ) | 7 | +1.07 [+0.32, +1.99] |
| R0', S = Aᵀ | 7 | +1.13 [+0.37, +1.95] |
| R0', S = A | 7 | +1.18 [+0.38, +2.06] |
| R0', S = max(A, Aᵀ) | 7 | +1.19 [+0.37, +2.04] |
| R2, S = Aᵀ | 7 | +1.20 [+0.41, +2.17] |
| R2, S = A | 7 | +1.25 [+0.46, +2.21] |
| R2, S = max(A, Aᵀ) | 7 | +1.24 [+0.43, +2.24] |
| R1, S = Aᵀ | 7 | +1.05 [+0.48, +1.82] |
| R1, S = A | 7 | +0.87 [+0.25, +1.66] |
| R1, S = max(A, Aᵀ) | 7 | +1.01 [+0.41, +1.82] |

**Arm B (GEO counts; hESC and mHSC-E/GM/L, n = 4; descriptive).** HL with the range of the four d_i; no interval or p-value (smallest one-sided p at n = 4 is 1/16).

| Contrast | n | HL estimate and 92.2% interval |
|---|---|---|
| R0, S = Aᵀ | 4 | +1.62 (range +0.10 to +2.28) |
| R0, S = A | 4 | +1.72 (range +0.24 to +2.37) |
| R0, S = max(A, Aᵀ) | 4 | +1.71 (range +0.16 to +2.39) |
| R0', S = Aᵀ | 4 | +1.76 (range +0.14 to +2.22) |
| R0', S = A | 4 | +1.82 (range +0.30 to +2.45) |
| R0', S = max(A, Aᵀ) | 4 | +1.84 (range +0.21 to +2.36) |
| R2, S = Aᵀ | 4 | +1.90 (range +0.09 to +2.47) |
| R2, S = A | 4 | +1.87 (range +0.36 to +2.65) |
| R2, S = max(A, Aᵀ) | 4 | +1.92 (range +0.24 to +2.61) |
| R1, S = Aᵀ | 4 | +1.54 (range +0.32 to +2.22) |
| R1, S = A | 4 | +1.35 (range +0.12 to +2.10) |
| R1, S = max(A, Aᵀ) | 4 | +1.51 (range +0.29 to +2.26) |

Input-scale sensitivity (E6): the paired difference Δ_i = log2 AUPRC ratio of scPRINT R0 under arm B minus under arm A is hESC -0.0003, mHSC-E -0.0079, mHSC-GM +0.0011, mHSC-L +0.0065; HL +0.0001, median +0.0004. The preregistered reading rule (median Δ_i > 0.15 means the primary result depends on the input scale on the count datasets) is **not met**. R0 under the two input scales is nearly identical (E13 QC: relative L1 0.04 to 0.09, Spearman 0.988 to 0.996).

## 5. Secondary checkpoints, controls and comparators (descriptive)

| Contrast | n | HL estimate and 92.2% interval |
|---|---|---|
| large-v1 R0 (reference CLR-NS) | 7 | +0.95 [+0.32, +1.74] |
| small-v1 R0 (reference CLR-NS) | 7 | +1.13 [+0.44, +1.75] |
| permutation control, mean of 3 seeds (reference CLR-NS) | 7 | +1.04 [+0.36, +1.88] |
| random weights, ESM2 embeddings kept, 3 seeds (reference CLR-NS) | 7 | +1.21 [+0.38, +2.16] |
| random weights, random embeddings, 3 seeds (reference CLR-NS) | 7 | +1.19 [+0.44, +2.03] |
| R0 minus permutation control (reference R0) | 7 | +0.01 [-0.00, +0.02] |
| R0 minus random weights, ESM2 kept (reference R0) | 7 | +0.13 [-0.02, +0.26] |
| R0 minus random weights, random embeddings (reference R0) | 7 | +0.10 [+0.01, +0.19] |
| R0 minus ESM2 gene-embedding cosine (reference R0) | 7 | -0.78 [-1.21, -0.35] |
| ESM2 gene-embedding cosine (reference CLR-NS) | 7 | +0.22 [-0.04, +0.62] |
| CLR applied to the symmetrised R0 (reference CLR-NS) | 7 | +1.05 [+0.34, +1.95] |

Reading, descriptive. The permutation control (each gene permuted independently across cells) scores the same as real-data R0 (R0 minus permutation +0.008, interval [-0.004, +0.018]). The random-weight networks score 0.10 to 0.13 log2 units below R0 on average (R0 minus random weights +0.130 [-0.020, +0.258] with ESM2 embeddings kept; +0.101 [+0.008, +0.193] with random embeddings). The data-free ESM2 gene-embedding cosine scores above R0 (R0 minus ESM2 cosine -0.782, [-1.208, -0.348]) while the interval for CLR-NS minus ESM2 cosine includes 0 (+0.224, [-0.042, +0.619]).

Comparators on the G-intersection and the other preregistered descriptive arms:

| Contrast | n | HL estimate and 92.2% interval |
|---|---|---|
| CLR-NS minus |r| | 7 | +0.28 [+0.13, +0.44] |
| CLR minus |r| | 7 | +0.23 [+0.09, +0.36] |
| CLR-NS minus CLR | 7 | +0.06 [+0.02, +0.10] |
| CLR-NS attention (soft10, held-out depth) minus scPRINT R0 (item 5) | 7 | +0.99 [+0.36, +1.83] |
| CLR-NS attention minus CLR-NS | 7 | -0.06 [-0.10, +0.01] |
| GENIE3 minus scPRINT R0 | 7 | +0.88 [+0.35, +1.69] |
| GRNBoost2 minus scPRINT R0 | 7 | +0.80 [+0.28, +1.54] |
| PIDC minus scPRINT R0 | 7 | +0.91 [+0.37, +1.64] |
| CLR-NS minus GENIE3 | 7 | +0.05 [-0.02, +0.21] |
| CLR-NS minus GRNBoost2 | 7 | +0.19 [+0.09, +0.32] |
| CLR-NS minus PIDC | 7 | +0.05 [-0.04, +0.16] |

GENIE3, GRNBoost2 and PIDC were fitted on the full gene sets (not refitted on the G-intersection) and scored on it.

Secondary truths (derived STRING variants) and gene sets, primary contrast:

| Contrast | n | HL estimate and 92.2% interval |
|---|---|---|
| primary + derived STRING truths | 7 | +1.15 [+0.51, +2.04] |
| derived STRING truths only | 7 | +1.28 [+0.59, +2.25] |
| TFs+500 only | 7 | +1.01 [+0.28, +1.91] |
| TFs+1000 only | 7 | +1.07 [+0.40, +1.82] |

## 6. scGPT whole-human (descriptive: hESC and hHep, n = 2; K562 in section 7)

n = 2 allows no interval or p-value (smallest one-sided p = 1/4) and no NI claim. Values are CLR-NS (or scPRINT R0) minus the scGPT readout in log2 AUPRC ratio, HL with range of the two d_i, on the scGPT G-intersection (CLR-NS) or the three-way intersection (scPRINT vs scGPT).

| Contrast | n | HL estimate and 92.2% interval |
|---|---|---|
| CLR-NS minus scGPT attention (last layer) | 2 | +0.37 (range +0.18 to +0.57) |
| CLR-NS minus scGPT token-embedding cosine | 2 | -0.78 (range -0.91 to -0.64) |
| CLR-NS minus scGPT hidden-state cosine | 2 | -0.29 (range -0.43 to -0.16) |
| scPRINT R0 minus scGPT attention | 2 | +0.08 (range +0.07 to +0.09) |
| scPRINT R0 minus scGPT token-embedding cosine | 2 | -1.07 (range -1.12 to -1.01) |
| scPRINT R0 minus scGPT hidden-state cosine | 2 | -0.58 (range -0.63 to -0.53) |

The LayerNorm'd token-embedding cosine (`tok_cos_ln`) was extracted in step 3 but is not a preregistered readout; its result (CLR-NS minus it: -0.754, range [-0.910, -0.598]) is exploratory.

## 7. K562 replication stratum (descriptive: no p-value, no NI claim)

Five subsamples of 1,024 non-targeting control cells, 1,244 genes (278 CollecTRI TFs with >= 10 targets plus the 1,000 most variable detected genes; scPRINT covers 1,218, scGPT 1,240, both 1,217). Truths: DoRothEA A-B (1,209 edges, 113 regulators in the gene set), CollecTRI (3,584 edges, 283 regulators), knockdown-response modules (3,235 edges, 166 regulators); on the scPRINT-covered genes 1,202, 3,582 and 2,874 edges remain. Directed TF-node scoring as for BEELINE (E14). d = log2(AUPRC ratio of the reference / that of the candidate), mean over the five subsamples; 95% percentile bootstrap over regulators (B = 1,000, one resample per truth shared by all subsamples and methods). The pooled interval is the percentile interval of the mean of the three truth-level draws. The bootstrap statistic on the full sample equals the `score_truth` value in all 21 contrast x truth cells.

Mean AUPRC ratio (EPR) over the five subsamples on the scPRINT-covered genes:

| Method | DoRothEA A-B | CollecTRI | Knockdown response |
|---|---|---|---|
| \|r\| | 1.70 (4.70) | 1.52 (3.47) | 1.12 (1.50) |
| CLR | 1.48 (2.59) | 1.46 (2.55) | 1.00 (1.06) |
| CLR-NS | 1.48 (3.43) | 1.37 (2.72) | 1.12 (1.35) |
| scPRINT medium-v1.5 R0 (Aᵀ) | 1.20 (0.10) | 1.06 (0.69) | 1.02 (1.08) |

Bootstrap contrasts (positive favours the first named method):

| Contrast | DoRothEA A-B | CollecTRI | Knockdown response | Mean of the three truths |
|---|---|---|---|---|
| CLR-NS minus scPRINT R0 | +0.310 [-0.163, +0.799] | +0.360 [+0.123, +0.654] | +0.134 [-0.013, +0.297] | **+0.268 [+0.090, +0.457]** |
| CLR-NS minus scGPT attention | +1.002 [+0.553, +1.627] | +0.712 [+0.403, +1.130] | +0.160 [+0.033, +0.332] | +0.625 [+0.425, +0.888] |
| CLR-NS minus scGPT token cosine | -0.714 [-1.257, -0.128] | -0.325 [-0.661, -0.016] | +0.096 [-0.039, +0.201] | -0.314 [-0.547, -0.092] |
| CLR-NS minus scGPT hidden cosine | +0.233 [-0.286, +0.861] | +0.214 [-0.080, +0.595] | +0.298 [+0.172, +0.440] | +0.248 [+0.044, +0.486] |
| scPRINT R0 minus scGPT attention (three-way set) | +0.722 [+0.154, +1.355] | +0.371 [+0.096, +0.672] | +0.012 [-0.092, +0.121] | +0.369 [+0.149, +0.616] |
| scPRINT R0 minus scGPT token cosine | -1.003 [-1.639, -0.397] | -0.664 [-1.015, -0.398] | -0.052 [-0.309, +0.111] | -0.573 [-0.815, -0.353] |
| scPRINT R0 minus scGPT hidden cosine | -0.045 [-0.643, +0.625] | -0.120 [-0.383, +0.134] | +0.201 [+0.069, +0.327] | +0.012 [-0.212, +0.252] |

Across the five subsamples the per-subsample d (CLR-NS minus scPRINT R0) ranges from +0.218 to +0.365 (DoRothEA), +0.228 to +0.437 (CollecTRI) and +0.112 to +0.149 (knockdown response). The other scPRINT readouts, the large-v1 and small-v1 checkpoints and the controls (point estimates, no interval) are in `results/scfm_score_k562/readouts_vs_clrns_mean_over_subsamples.csv`. In the knockdown-response truth the readouts R0 under orientation A, R1 and R2 have a higher AUPRC ratio than R0 under Aᵀ (1.89, 1.31, 1.17 against 1.02) and than CLR-NS (1.12); that is, on that truth the CLR-NS minus candidate difference is negative for R0 (S = A: -0.76), R1 (-0.23) and R2 (-0.07). These are descriptive, not preregistered K562 contrasts; the preregistered K562 contrast is R0 with S = Aᵀ. In the CSV, EPRs of exactly 0 (small-v1 on DoRothEA, 4 of 5 subsamples) are floored at 1e-9, so that row's EPR difference (+26.4) is a floor artefact.

## 8. Dropped genes and truth edges (BEELINE)

Per case: `results/scfm_score/dropped_genes_and_truth_edges.csv` (genes in the full case set, in G-intersection sets A, B, S and AS; truth edges before and after). On the primary truths, 276,960 truth edges lie in the full gene sets and 276,483 remain on the scPRINT-covered sets (477 dropped, 0.17%). The largest losses are mESC (98 cell-type ChIP, 155 non-specific ChIP, 73 LOF/GOF edges). TFs lost: 2 of 321 in mDC.

## 9. Deviations, interpretations and caveats

1. **Statistics bug found and fixed after the real run.** In the first statistics run, the arm-A scPRINT R0 (S = Aᵀ) rows were picked without restricting the gene set, so for hESC and hHep the descriptive R0 row, the robust-NI candidate list and the arm-B pairing also contained the same cells scored on the scGPT-restricted gene set (extra rows differing in the fifth decimal). The primary, S1 and S2 always used gene set A. After the fix the statistics stage (not the scoring) was rerun on the same metric table; the changes are at most 1.3e-4 in one descriptive row (HL, L) and 8.5e-5 in the arm-B row, and no decision changed. The committed outputs are the corrected ones.
2. **"CLR-NS" is the plain CLR-NS score matrix** (`clr_ns`: CLR of |r| of Blom scores), as in `results/beeline_baselines`. The CLR-NS soft10 attention at its held-out depth (preregistered item 5) is the descriptive row in section 5 (held-out depth recomputed on the G-intersection with seed 20261021 + case index).
3. **Robust-NI maximum.** "Best" is the configuration with the smallest d_i per dataset (the one most favourable to scPRINT in the contrast). For arm-B configurations the comparator is CLR-NS recomputed on the arm-B gene set (one gene fewer in mHSC-E_1000 and mHSC-GM_1000).
4. **CLR applied to A (hybrid)** is implemented as the CLR calibration of the symmetrised R0, (A + Aᵀ)/2, because the calibration expects a symmetric score.
5. **Mixed model:** `lme4` is not installed; the same model (random intercepts for dataset and for truth within dataset) was fitted with `nlme::lme` (REML).
6. **GENIE3, GRNBoost2, PIDC** come from `results/beeline_baselines` (fits on the full gene set) and are scored on the G-intersection; they are not refitted.
7. **Decision labels** are attached only to the preregistered tests (primary, S1, S2, robust NI); descriptive rows show intervals only. The tables in `contrasts_summary.csv` list p-values for some descriptive n = 7 rows because the code computes them; the preregistration attaches no claim to them.
8. **K562 definitions not fixed by the preregistration** were fixed in E14 before extraction (detection pooled over all controls, variance of log-normalised expression, E5 gene mapping, raw UMI counts as model input, TF-node scoring with edges from the three truths, bootstrap over regulators). No ESM2-cosine prior and no fine-tuning arm for K562; the fine-tuning arm was not run for any dataset.
9. **mHSC-E_500 organism-head flag** (E13: 0.67 of its cells are classified human by the organism head, every other case 0.99 to 1.00): no rule was applied to it. Exploratory check in section 3: the primary conclusion does not depend on it.
10. **Arm B has n = 4** and **scGPT has n = 2 on BEELINE** (hESC, hHep) plus K562; both are descriptive. mDC, mESC and hHep have no deposited counts (no arm B), and no SRA re-quantification was done (PI decision 2026-11-01).
11. **Leakage:** every scPRINT v1 and scGPT pair is classified "out, moderate confidence" (E9), so the leaked-dataset sensitivity analysis is empty. Residual limits of the audit are in E9.
12. **EPR floor:** two BEELINE cells (mHSC-L, STRING, both gene sets) have scPRINT R0 EPR = 0, and on K562 small-v1 has EPR = 0 on DoRothEA in four subsamples. Floored at 1e-9.
13. **Intervals are on dataset-level means.** With n = 7 the exact interval covers 92.2%, not 95%, by preregistration.
14. **Power:** the p of 0.00781 is the smallest value attainable at n = 7 and occurs whenever all seven differences exceed the margin, which is the case here; it does not measure the size of the effect.
15. **The truths are incomplete and the scores are small.** On the 44 primary cells the AUPRC ratio of scPRINT R0 ranges from 0.79 to 1.70 and that of CLR-NS from 0.94 to 32.3 (STRING, mHSC-L). The mean-attention map of scPRINT is close to the same prior with and without the data (permutation control, section 5).
16. **Scheduler and compute:** the K562 scPRINT runs took about 58 s per case (hESC_1000 in step 3: 11 s); cause not investigated.

## 10. Plain-language statement

1. We asked whether CLR attention, a formula with no training, finds known regulator-target links at least as well as the attention inside scPRINT, a model trained on tens of millions of cells.
2. On the seven BEELINE datasets the answer under the rule we fixed beforehand is yes: CLR-NS was at least as good by the 0.15 margin and the interval of its advantage lies entirely above zero (typical advantage: 1.03 log2 units, interval 0.34 to 1.90).
3. The advantage holds when scPRINT is given its best readout, head selection and input scale (robust non-inferiority), and the two key secondary tests agree.
4. It comes from the non-specific ChIP-seq and STRING truths; on cell-type ChIP-seq, scPRINT scored higher than CLR-NS in 11 of 14 cases, and the permutation control scored the same as the real data.
5. On the K562 replication (descriptive) the direction is the same but smaller (0.27 log2 units, interval 0.09 to 0.46); scGPT is mixed (CLR-NS is better than its attention, worse than its token-embedding cosine).
6. These are comparisons against known networks that are incomplete, with seven datasets and no claim that either method recovers causal regulation.

## Files

`results/scfm_score/` (metrics_long.csv, contrasts_summary.csv, stats.log, primary_cells.csv, robust_*.csv, armB_paired_cells.csv, key_metrics_per_case.csv, dropped_genes_and_truth_edges.csv, mixed_model_fixed.csv, exploratory.txt), `results/scfm_score_k562/`, `results/scfm_score_smoke/`, `results/scfm_k562/`. Scripts: `analysis/scfm_score_export.py`, `scfm_score.R`, `scfm_score_stats.R`, `scfm_score_k562.R`, `scfm_score_synth_truth.R`, `scfm_score_exploratory.R`, `scfm_k562_prep.py`, `scfm_k562_manifest.py`.
