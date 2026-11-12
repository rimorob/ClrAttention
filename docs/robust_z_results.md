# Robust row-standardisation in the CLR step (EXPLORATORY)

**Status.** Exploratory pathfinding. Not preregistered, not for the paper or any paper table. K562 CORUM is not held out (it was used in the keff test and in `analysis/corum_exploratory/`); M3D RegulonDB is the project's tuning set for the depth. No p-values are computed and none should be inferred from the intervals. No variant was chosen after looking at results and then re-tested on the same data; every variant, truth and read-out was fixed before the runs (task brief of 11 Nov 2026).

**Question (the PI's idea).** CLR z-scores each gene's MI row with the row mean and SD. A gene whose own module has many strong partners inflates its row SD, which lowers the z-scores of its true partners (this could hurt large complexes such as the ribosome). Does z_ij = (x_ij - median(x_i*)) / (IQR(x_i*)/1.349) change pairwise or attention performance?

## Variants

Everything after the row z is the package's default path: truncation at 0, symmetrisation (z_ij + z_ji)/sqrt(2), per-gene softmax (about 10 effective partners, cap 50), lazy walk alpha 0.5, symmetrised walk mass (P^t + P^t')/2. The variants are implemented in `analysis/robust_z/robust_z_lib.R` on the MI matrix and feed the package's own `clr:::.softmax_rows`; the package itself is unchanged. Names are distinct from the project's existing "rCLR".

| Name | Row centre / scale | Notes |
|---|---|---|
| V0 | mean / SD (default CLR) | package convention: the zeroed diagonal is part of the row |
| V0x | mean / SD, diagonal excluded | control for the diagonal convention only |
| V1 | median / (IQR/1.349) of the raw MI row | the PI's proposal; diagonal excluded |
| V2 | median / (IQR/1.349) of MI^(1/3) | Wilson-Hilferty scale, one-way; diagonal excluded |
| V3 | existing rCLR (`eb_additive`, two-way median polish of MI^(1/3), per-gene MAD) | reference only; as defined it is not truncated at 0 |
| V4 | median / (1.4826 MAD) of the raw MI row | diagnostic; diagonal excluded |

Quantiles are type 7 (R default). **Zero-scale guard.** A row whose scale is 0 gets z = 0 for the whole row (the package's own treatment of sd = 0). Rows that hit it: **0 for every variant in both datasets.**

## Equivalence check (before anything else)

V0 computed by `robust_z_lib.R` against the package (`ClrAttention$calibrate()`, `$build_operator(softmax_keff = 10, softmax_cap = 50)`, and the P^9 walk mass of `analysis/baselines.R`), full data, asserted to < 1e-8 by the script (it stops otherwise):

| Dataset | max abs diff, CLR scores | operator | walk mass t = 9 |
|---|---|---|---|
| M3D | 7.1e-14 | 7.9e-15 | 1.1e-15 |
| K562 | 3.4e-14 | 1.2e-14 | 1.1e-15 |

V0 pooled AUPR also equals the project's `results/baselines_m3d` and `results/baselines_k562ess_truthgenes_genie3` values for CLR and attention to < 3e-16 (see "V0 against the project's baselines run" in the appendix). The diagonal convention is immaterial: max |C(V0) - C(V0x)| = 0.010 (M3D), 0.0053 (K562); correlation > 0.99999998; V0x - V0 contrasts are 0.001 or smaller in every pooled table.

## Data, truths, metrics

- **M3D**: 466 experiments, 4,297 genes; RegulonDB strong evidence (SC, base rate 0.0943) and all evidence (all, 0.1205), the `analysis/baselines.R` constructions, same-operon pairs removed.
- **K562**: essential Perturb-seq pseudobulk, half A, 5,671 truth genes. CORUM all (1,185 complexes; base rate 0.0133), CORUM excluding the 8 RP-dominated complexes ("noRP"), the RP complexes alone ("RP"), small (3-6 genes), large (>= 20 genes, as `analysis/rclr_k562/rclr_k562.R`; the brief's ">= 25?" was not used), DoRothEA A-B, KD-response modules built from half B.
- **Pooled metric**: AUPR over base rate ("enrichment"), pairs within the genes of the truth. Contrasts are differences of enrichment with a paired delete-half gene jackknife (100 reps, seed 20261016, stream order of `baselines.R`); interval = difference +- 1.96 SD of the paired replicate differences. Not a confirmatory interval.
- **Macro**: per-regulon / per-complex AP (`regulon_ap_table`, regulons with >= 10 positive pairs), mean over regulons. The CORUM subsets (noRP, RP, small, large) are subsets of the CORUM per-complex table. Intervals: regulon bootstrap (percentile, 2,000 draws) and delete-20% gene jackknife (first 100 replicates of the `baselines.R` d20 stream; variance factor 4; regulon set fixed to those scored in every replicate).
- **Attention read-outs**: t = 9 (the project's depth) and the label-free held-out depth t* of `analysis/lab_splits.R` (AUPR of the attention mass built on one half for the other half's top-1% pairs of the *same variant's* pairwise score; 8 balanced lab splits for M3D, 6 random profile halves for K562, both directions, t = 1..20; t* = median over splits of the geometric mean of the two directions). That criterion is computed on each variant's own score, so it is not comparable across variants as an accuracy; it only fixes the depth. The V0 t* (9 in both datasets) equals the project's value.
- **Hub control**: within-stratum percentile in 10 strata of H_i + H_j. Primary stratifier: H = pi, the stationary distribution of the **V0 attention walk** (attention-derived; the stratifier of the sota table; the same for every variant). Second, method-neutral stratifier: H = the row mean of raw MI (does not depend on any variant). The project's CLR-degree stratifier is still only planned.
- **Hub drift**: among the 100,000 top-scoring pairs of a score matrix, the share with at least one endpoint among the 5% most connected genes, and the mean percentile rank of the two endpoints, with connectedness taken from pi (V0 attention) and from the MI row mean. Chance: share 0.0975, rank 0.5.

## Jobs

| Job | Runtime | Peak memory | Note |
|---|---|---|---|
| `4000_rz_m3d` | killed after 1.5 min | 16.7 GB (declared 16) | OOM kill in the depth stage (8 forked workers, about 2.3 GB each); scheduler auto-retry |
| `4000_rz_m3d_retry` | 57 min | 24.9 GB (declared 26) | the M3D run reported here |
| `4010_rz_k562` | 4 h 38 min | 40.3 GB (declared 36, cap 47) | the K562 run reported here; the declaration was raised from 30 to 36 GB before it started |

Nothing in the code changed between the first M3D attempt and its retry (the retry only raised the memory declaration, automatically). Outputs: `results/robust_z/{m3d,k562}/` (tables, CSV, logs; the rds caches of the replicate-level results are not committed). Code: `analysis/robust_z/robust_z_lib.R`, `robust_z.R`, `robust_z_summarize.R`.

## Reading, limited to what the tables show

All statements are descriptions of point estimates and intervals in the tables below (and the appendix); intervals resample genes only (networks held fixed), not experiments, cells or seeds.

**1. The PI's mechanism, as measured on the rows.** The row SD rises with a gene's number of globally strong partners (partners whose MI is in the top 1% of all cells): Spearman 0.86 (M3D) and 0.89 (K562). The robust scales rise too (IQR/1.349: 0.72 and 0.85; MAD: 0.71 and 0.83), so the robust scales are not independent of the number of strong partners. SD relative to the IQR scale rises with the number of strong partners but modestly: median SD/IQR-scale goes from 1.10 (lowest decile) to 1.19 (highest) in M3D and from 1.07 to 1.35 in K562 (Spearman of the ratio with the count: 0.35 and 0.77). In the top decile of K562 rows the SD is about 1.35 times the IQR scale and about 1.6 times the MAD scale. Against a gene's largest complex (K562) or regulon (M3D) the row scales are only weakly related (SD: 0.17 and -0.05).

**2. Pairwise score (CLR itself).** In M3D every robust variant is at or above V0 in pooled enrichment: V1 +0.013 [-0.004, +0.029] (SC) and +0.047 [+0.035, +0.059] (all); V2 +0.022 and +0.043; V3 +0.077 and +0.102; V4 +0.009 and +0.058. In K562: V1 raises CORUM all from 6.39 to 6.94 (+0.547 [+0.345, +0.749]) and KD response from 2.42 to 2.61 (+0.189 [+0.158, +0.220]); V2, V3 and V4 do not move CORUM all (intervals include 0, width about 0.5 to 0.9) but do raise KD response (+0.10 to +0.18). CORUM without ribosome: V1 +0.100 [+0.019, +0.182], V2 +0.158 [+0.050, +0.266], V3 +0.381 [+0.182, +0.581], V4 -0.019 [-0.155, +0.117]. Large complexes: V1 +0.023 [-0.000, +0.045], V3 +0.097 [+0.046, +0.149], V2 and V4 include 0. Small complexes: all four intervals include 0 (point estimates -0.23 to +0.07). The RP-only truth is saturated (base rate 0.75; enrichment 1.21 for every variant), so pooled enrichment cannot show a ribosome effect; the macro table (8 RP complexes, 7 scored) shows pairwise V1 +0.024 AP [gene interval +0.009, +0.038], the other variants smaller with intervals including 0. DoRothEA A-B: all within 0.03 of V0 (V3 +0.025 [+0.013, +0.037], V4 -0.008 [-0.013, -0.003]).

**3. Attention at t = 9.** M3D: no variant differs from V0 in pooled enrichment with an interval excluding 0 (V2 on all, -0.023 [-0.046, +0.001], is the nearest); at the variant's own held-out depth (6 to 7 instead of 9) V1 to V4 are lower on SC (-0.028 to -0.050, intervals exclude 0) and V2 is lower on all (-0.034 [-0.054, -0.014]). K562: V1 is lower on large complexes (-0.127 [-0.171, -0.082]), on CORUM all (-0.250 [-0.539, +0.039]) and DoRothEA (-0.017 [-0.027, -0.007]) and higher on KD response (+0.040 [+0.002, +0.078]); V2 is lower on KD response (-0.181 [-0.231, -0.130]) and includes 0 elsewhere; V4 is lower on CORUM all (-0.519 [-0.980, -0.059]), CORUM without ribosome (-0.196 [-0.390, -0.001]), large (-0.180 [-0.252, -0.108]) and DoRothEA (-0.025), higher on KD (+0.071 [+0.019, +0.123]); V3 is higher on large complexes (+0.158 [+0.051, +0.265]) and CORUM without ribosome (+0.216 [+0.031, +0.401]) and lower on KD (-0.302 [-0.386, -0.218]). In macro terms (per-complex mean AP) V1, V2 and V4 attention are below V0 on RP complexes (-0.030, -0.027, -0.020) and large complexes (-0.012, -0.010, -0.012), V3 above (+0.085, +0.016); M3D macro AP is lower than V0 for every robust variant in both truths for pairwise and attention (-0.005 to -0.014; for V3 attention the gene interval includes 0).

**4. Attention versus CLR.** K562 CORUM all and large complexes: for every variant (V0 to V4) attention minus the V0 pairwise score is negative with an interval excluding 0 (CORUM all -2.6 to -3.4 enrichment units, V0 -2.88 [-4.32, -1.44]; large -0.21 to -0.54). K562 KD response: V2 and V3 negative with intervals excluding 0 (-0.230, -0.351), V0, V1, V4 include 0. M3D: attention is above V0 CLR on both truths for every variant (SC +0.21 to +0.23; all +0.05 to +0.08, intervals exclude 0).

**5. Hub control.** With the attention-derived stratifier (pi of V0 attention), M3D pairwise gains shrink (V1 SC -0.004 [-0.018, +0.009]; all +0.020 [+0.010, +0.030]) and attention differences from V0 are -0.02 to -0.04 for V1, V2, V4 (several intervals exclude 0). With the MI-row-mean stratifier the pairwise gains are close to the unstratified ones (all: +0.044 to +0.073). In K562 the V1 pairwise gain on CORUM all vanishes under the pi stratifier (+0.014 [-0.178, +0.206]) and persists under the MI-row-mean stratifier (+0.483 [+0.234, +0.733]). The two stratifiers therefore disagree on this contrast; the pi stratifier is derived from V0's own attention walk.

**6. Hub drift.** By pi, the top pairs of V1, V2, V3 and V4 touch the most connected genes more often than V0's in M3D, and in K562 V1, V2 and V4 do while V3's touch them less often than V0's (share touching the top 5%; M3D pairwise: V0 0.211, V1 0.277, V2 0.279, V3 0.270, V4 0.276; attention t = 9: V0 0.257, V1 0.348, V2 0.328, V3 0.294, V4 0.348; K562 pairwise: 0.336, 0.373, 0.368, 0.328, 0.347; attention t = 9: 0.608, 0.711, 0.636, 0.543, 0.685; chance 0.0975). Ranked by the MI row mean V3 is below V0 in both datasets (K562 pairwise 0.259 against 0.419; M3D 0.034 against 0.169), while V1 and V4 are above V0 in K562 (0.507, 0.564 against 0.419) and close to V0 in M3D (0.175, 0.205 against 0.169). The share of non-zero pairwise scores is higher for V1, V2, V4 (0.62 to 0.64) than for V0 (0.53 to 0.54): the median-centred z are positive for more pairs.

**7. Held-out depth.** V0 gives t* = 9 in both datasets (M3D IQR 9-10; K562 8.25-9). V1, V2, V4 give 6 to 7 and V3 5 to 6 (K562: 6, 6, 6, 5). Attention at each variant's own t* is in the "attention t*" contrast tables (each variant at its own t* against V0 at 9); for several variant-truth cells the difference is more negative than at t = 9 (for example K562 CORUM all, V1: -0.398 [-0.704, -0.093] against -0.250 [-0.539, +0.039] at t = 9).

**What the tables do not show.** A consistent direction of the effect of the robust scale across datasets, read-outs or variants. In K562 attention at t = 9, V1, V2 and V4 are lower than V0 on large complexes and V3 is higher; in pairwise score V1 is higher on CORUM all and KD response. The pooled CORUM all estimate is dominated by the RP complexes, whose own pooled enrichment is saturated, so pooled CORUM all differences cannot be attributed to ribosome or non-ribosome complexes from the pooled number alone (the noRP and RP rows and the macro table are the split). Nothing here is a held-out test, and the order of the variants is not a recommendation.


## Appendix: all tables

Abbreviations: `truth` names as in the Data section (E. coli SC / all; K562 DoRothEA_AB, CORUM_all, KD_response, CORUM_noRP, CORUM_RP, CORUM_small, CORUM_large); pair = pairwise score, att9 = attention at t = 9, attT = attention at the variant's own held-out depth t*. Hub-stratified contrasts name the stratifier in the heading. 'Macro contrast' cells: mean difference in per-regulon AP; regulon-bootstrap interval; gene-jackknife interval. Generated by `analysis/robust_z/robust_z_summarize.R`; CSVs in `results/robust_z/<dataset>/` (`pooled_metrics.csv`, `contrasts_pooled.csv`, `contrasts_macro.csv`, `macro_levels.csv`, `hub_drift.csv`, `row_scale_*.csv`, `tstar_summary.csv`, `heldout_curves.csv`, `equivalence.csv`, `check_vs_baselines.csv`).

### Tables for m3d (generated by analysis/robust_z/robust_z_summarize.R; EXPLORATORY)

#### Pooled enrichment (AUPR / base rate), pairwise score

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC (base 0.0943) | 1.44 | 1.45 | 1.46 | 1.52 | 1.45 |
| all (base 0.1205) | 1.37 | 1.42 | 1.41 | 1.47 | 1.43 |

#### Pooled enrichment, attention t = 9

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC (base 0.0943) | 1.68 | 1.67 | 1.66 | 1.66 | 1.66 |
| all (base 0.1205) | 1.44 | 1.44 | 1.42 | 1.43 | 1.45 |

#### Pooled enrichment, attention at each variant's own held-out depth

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC (base 0.0943) | 1.68 | 1.65 | 1.63 | 1.63 | 1.63 |
| all (base 0.1205) | 1.44 | 1.43 | 1.41 | 1.41 | 1.43 |

#### Paired contrast, pairwise: variant - V0, raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V0x | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC | -0.000 [-0.000, -0.000] | +0.013 [-0.004, +0.029] | +0.022 [+0.003, +0.040] | +0.077 [+0.042, +0.112] | +0.009 [-0.008, +0.026] |
| all | -0.000 [-0.000, -0.000] | +0.047 [+0.035, +0.059] | +0.043 [+0.031, +0.055] | +0.102 [+0.072, +0.131] | +0.058 [+0.043, +0.073] |

#### Paired contrast, attention t=9: variant - V0, raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V0x | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC | +0.000 [-0.000, +0.000] | -0.009 [-0.037, +0.020] | -0.014 [-0.047, +0.020] | -0.020 [-0.070, +0.030] | -0.011 [-0.038, +0.016] |
| all | -0.000 [-0.000, +0.000] | +0.000 [-0.018, +0.019] | -0.023 [-0.046, +0.001] | -0.010 [-0.037, +0.018] | +0.011 [-0.009, +0.031] |

#### Paired contrast, attention t*: variant - V0, raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.028 [-0.052, -0.005] | -0.044 [-0.072, -0.017] | -0.050 [-0.095, -0.004] | -0.042 [-0.064, -0.020] |
| all | -0.011 [-0.027, +0.005] | -0.034 [-0.054, -0.014] | -0.024 [-0.050, +0.002] | -0.010 [-0.026, +0.007] |

#### Paired contrast, attention t=9 minus CLR(V0 pairwise), raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC | +0.234 [+0.169, +0.298] | +0.225 [+0.149, +0.301] | +0.220 [+0.141, +0.299] | +0.214 [+0.130, +0.298] | +0.223 [+0.150, +0.296] |
| all | +0.071 [+0.039, +0.103] | +0.071 [+0.035, +0.108] | +0.048 [+0.009, +0.088] | +0.061 [+0.020, +0.102] | +0.082 [+0.046, +0.118] |

#### Paired contrast, pairwise: variant - V0, hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.004 [-0.018, +0.009] | +0.003 [-0.011, +0.017] | +0.051 [+0.016, +0.086] | -0.007 [-0.023, +0.009] |
| all | +0.020 [+0.010, +0.030] | +0.016 [+0.006, +0.026] | +0.068 [+0.038, +0.097] | +0.031 [+0.016, +0.045] |

#### Paired contrast, attention t=9: variant - V0, hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.022 [-0.044, -0.001] | -0.024 [-0.053, +0.004] | -0.010 [-0.055, +0.034] | -0.022 [-0.043, -0.001] |
| all | -0.032 [-0.047, -0.017] | -0.044 [-0.064, -0.025] | -0.016 [-0.041, +0.009] | -0.026 [-0.042, -0.009] |

#### Paired contrast, attention t*: variant - V0, hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.029 [-0.048, -0.010] | -0.032 [-0.056, -0.008] | -0.016 [-0.055, +0.023] | -0.032 [-0.052, -0.013] |
| all | -0.028 [-0.042, -0.014] | -0.034 [-0.051, -0.017] | -0.007 [-0.030, +0.016] | -0.020 [-0.034, -0.005] |

#### Paired contrast, attention t=9 minus CLR(V0 pairwise), hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC | +0.157 [+0.100, +0.213] | +0.134 [+0.071, +0.197] | +0.132 [+0.064, +0.200] | +0.146 [+0.070, +0.222] | +0.135 [+0.074, +0.196] |
| all | +0.016 [-0.015, +0.048] | -0.016 [-0.053, +0.022] | -0.028 [-0.069, +0.012] | +0.000 [-0.043, +0.044] | -0.009 [-0.046, +0.027] |

#### Paired contrast, pairwise: variant - V0, hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | +0.007 [-0.010, +0.024] | +0.007 [-0.010, +0.024] | +0.013 [-0.012, +0.037] | +0.011 [-0.008, +0.030] |
| all | +0.050 [+0.037, +0.063] | +0.044 [+0.032, +0.056] | +0.073 [+0.055, +0.091] | +0.063 [+0.047, +0.079] |

#### Paired contrast, attention t=9: variant - V0, hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.013 [-0.041, +0.015] | -0.026 [-0.057, +0.006] | -0.037 [-0.087, +0.012] | -0.013 [-0.040, +0.014] |
| all | +0.004 [-0.015, +0.022] | -0.021 [-0.044, +0.002] | -0.012 [-0.038, +0.014] | +0.014 [-0.006, +0.034] |

#### Paired contrast, attention t*: variant - V0, hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.034 [-0.057, -0.010] | -0.057 [-0.083, -0.030] | -0.068 [-0.113, -0.023] | -0.045 [-0.068, -0.023] |
| all | -0.008 [-0.024, +0.008] | -0.033 [-0.052, -0.014] | -0.026 [-0.050, -0.003] | -0.007 [-0.024, +0.010] |

#### Paired contrast, attention t=9 minus CLR(V0 pairwise), hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC | +0.182 [+0.124, +0.240] | +0.169 [+0.101, +0.237] | +0.156 [+0.088, +0.225] | +0.145 [+0.071, +0.218] | +0.169 [+0.103, +0.235] |
| all | +0.031 [+0.002, +0.059] | +0.034 [+0.001, +0.067] | +0.010 [-0.025, +0.044] | +0.019 [-0.015, +0.052] | +0.044 [+0.011, +0.078] |

#### Macro: mean per-regulon AP, pairwise

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC (n=123) | 0.075 | 0.062 | 0.061 | 0.064 | 0.062 |
| all (n=157) | 0.068 | 0.060 | 0.058 | 0.059 | 0.059 |

#### Macro: mean per-regulon AP, attention t = 9

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC (n=123) | 0.083 | 0.075 | 0.074 | 0.076 | 0.074 |
| all (n=157) | 0.076 | 0.070 | 0.070 | 0.071 | 0.069 |

#### Macro contrast, pairwise: variant - V0 (mean difference in per-regulon AP; regulon bootstrap and delete-20% gene jackknife intervals)

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.012; reg [-0.021, -0.006]; gene [-0.019, -0.002] | -0.014; reg [-0.021, -0.008]; gene [-0.021, -0.004] | -0.011; reg [-0.017, -0.005]; gene [-0.019, -0.003] | -0.012; reg [-0.023, -0.006]; gene [-0.017, -0.002] |
| all | -0.007; reg [-0.012, -0.004]; gene [-0.010, -0.002] | -0.010; reg [-0.014, -0.006]; gene [-0.012, -0.003] | -0.008; reg [-0.013, -0.004]; gene [-0.013, -0.002] | -0.008; reg [-0.016, -0.004]; gene [-0.010, -0.002] |

#### Macro contrast, attention t=9: variant - V0 (mean difference in per-regulon AP; regulon bootstrap and delete-20% gene jackknife intervals)

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| SC | -0.008; reg [-0.012, -0.005]; gene [-0.012, -0.002] | -0.009; reg [-0.012, -0.005]; gene [-0.013, -0.003] | -0.007; reg [-0.012, -0.002]; gene [-0.014, +0.001] | -0.009; reg [-0.013, -0.005]; gene [-0.012, -0.002] |
| all | -0.006; reg [-0.009, -0.004]; gene [-0.008, -0.001] | -0.006; reg [-0.009, -0.004]; gene [-0.008, -0.000] | -0.005; reg [-0.009, -0.002]; gene [-0.009, +0.001] | -0.007; reg [-0.011, -0.005]; gene [-0.009, -0.001] |

#### Macro contrast, attention t=9 minus CLR(V0 pairwise) (mean difference in per-regulon AP; regulon bootstrap and delete-20% gene jackknife intervals)

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| SC | +0.008; reg [-0.000, +0.019]; gene [-0.007, +0.005] | +0.000; reg [-0.008, +0.011]; gene [-0.015, -0.000] | -0.000; reg [-0.009, +0.011]; gene [-0.016, -0.002] | +0.002; reg [-0.008, +0.013]; gene [-0.016, +0.001] | -0.000; reg [-0.008, +0.009]; gene [-0.016, -0.001] |
| all | +0.008; reg [+0.001, +0.017]; gene [-0.005, +0.008] | +0.002; reg [-0.004, +0.011]; gene [-0.010, +0.004] | +0.002; reg [-0.005, +0.010]; gene [-0.010, +0.004] | +0.003; reg [-0.005, +0.012]; gene [-0.010, +0.004] | +0.001; reg [-0.006, +0.009]; gene [-0.010, +0.003] |

#### Hub drift (top pairs of each score against two connectedness rankings; chance share = 0.0975, chance mean rank = 0.5)

| variant | kind | stratifier | share_touching_top5pct | mean_endpoint_rank | frac_positive_score |
|---|---|---|---|---|---|
| V0 | pair | att_pi | 0.211 | 0.675 | 0.529 |
| V0 | pair | mi_rowmean | 0.169 | 0.608 | 0.529 |
| V0 | att9 | att_pi | 0.257 | 0.656 | 1.000 |
| V0 | att9 | mi_rowmean | 0.116 | 0.552 | 1.000 |
| V0 | attT | att_pi | 0.257 | 0.656 | 1.000 |
| V0 | attT | mi_rowmean | 0.116 | 0.552 | 1.000 |
| V1 | pair | att_pi | 0.277 | 0.717 | 0.635 |
| V1 | pair | mi_rowmean | 0.175 | 0.590 | 0.635 |
| V1 | att9 | att_pi | 0.348 | 0.673 | 1.000 |
| V1 | att9 | mi_rowmean | 0.129 | 0.536 | 1.000 |
| V1 | attT | att_pi | 0.305 | 0.656 | 1.000 |
| V1 | attT | mi_rowmean | 0.125 | 0.536 | 1.000 |
| V2 | pair | att_pi | 0.279 | 0.720 | 0.635 |
| V2 | pair | mi_rowmean | 0.130 | 0.555 | 0.635 |
| V2 | att9 | att_pi | 0.328 | 0.669 | 1.000 |
| V2 | att9 | mi_rowmean | 0.112 | 0.516 | 1.000 |
| V2 | attT | att_pi | 0.270 | 0.641 | 1.000 |
| V2 | attT | mi_rowmean | 0.109 | 0.520 | 1.000 |
| V3 | pair | att_pi | 0.270 | 0.691 | 0.500 |
| V3 | pair | mi_rowmean | 0.034 | 0.469 | 0.500 |
| V3 | att9 | att_pi | 0.294 | 0.650 | 1.000 |
| V3 | att9 | mi_rowmean | 0.077 | 0.504 | 1.000 |
| V3 | attT | att_pi | 0.244 | 0.620 | 1.000 |
| V3 | attT | mi_rowmean | 0.077 | 0.506 | 1.000 |
| V4 | pair | att_pi | 0.276 | 0.714 | 0.635 |
| V4 | pair | mi_rowmean | 0.205 | 0.613 | 0.635 |
| V4 | att9 | att_pi | 0.348 | 0.674 | 1.000 |
| V4 | att9 | mi_rowmean | 0.133 | 0.543 | 1.000 |
| V4 | attT | att_pi | 0.283 | 0.645 | 1.000 |
| V4 | attT | mi_rowmean | 0.127 | 0.542 | 1.000 |
| V0x | pair | att_pi | 0.211 | 0.675 | 0.529 |
| V0x | pair | mi_rowmean | 0.169 | 0.608 | 0.529 |
| V0x | att9 | att_pi | 0.257 | 0.656 | 1.000 |
| V0x | att9 | mi_rowmean | 0.116 | 0.552 | 1.000 |

#### Zero-scale guard (rows with scale 0)

| variant | rows |
|---|---|
| V0 | 0 |
| V0x | 0 |
| V1 | 0 |
| V2 | 0 |
| V3 | 0 |
| V4 | 0 |

#### V0 against the project's baselines run

| truth | V0_pair_aupr | baselines_clr_aupr | V0_att_aupr | baselines_attention_aupr | max_abs_diff |
|---|---|---|---|---|---|
| SC | 0.1358703 | 0.1358703 | 0.1579082 | 0.1579082 | 2.775558e-16 |
| all | 0.1647927 | 0.1647927 | 0.1733404 | 0.1733404 | 1.665335e-16 |

#### Row scale versus number of strong partners (Spearman)

| row_quantity | versus | spearman | n_genes |
|---|---|---|---|
| mean_diag0 | n_strong_global | +0.774 | 4297 |
| mean_diag0 | n_z4_gt3 | +0.259 | 4297 |
| mean_diag0 | module_size | -0.064 | 1879 |
| sd_diag0 | n_strong_global | +0.863 | 4297 |
| sd_diag0 | n_z4_gt3 | +0.430 | 4297 |
| sd_diag0 | module_size | -0.047 | 1879 |
| median | n_strong_global | +0.713 | 4297 |
| median | n_z4_gt3 | +0.144 | 4297 |
| median | module_size | -0.077 | 1879 |
| iqr_scale | n_strong_global | +0.717 | 4297 |
| iqr_scale | n_z4_gt3 | +0.223 | 4297 |
| iqr_scale | module_size | -0.054 | 1879 |
| mad_scale | n_strong_global | +0.706 | 4297 |
| mad_scale | n_z4_gt3 | +0.174 | 4297 |
| mad_scale | module_size | -0.065 | 1879 |
| sd_over_iqr | n_strong_global | +0.346 | 4297 |
| sd_over_iqr | n_z4_gt3 | +0.785 | 4297 |
| sd_over_iqr | module_size | +0.096 | 1879 |
| sd_over_mad | n_strong_global | +0.494 | 4297 |
| sd_over_mad | n_z4_gt3 | +0.894 | 4297 |
| sd_over_mad | module_size | +0.099 | 1879 |
| mean_minus_median_in_iqr | n_strong_global | +0.616 | 4297 |
| mean_minus_median_in_iqr | n_z4_gt3 | +0.906 | 4297 |
| mean_minus_median_in_iqr | module_size | +0.041 | 1879 |

#### Row scale by decile of the number of globally strong partners

| decile | n_genes | n_strong_global_median | sd_median | iqr_scale_median | mad_scale_median | sd_over_iqr_median | sd_over_mad_median |
|---|---|---|---|---|---|---|---|
| 1 | 430 | 0 | 0.0223 | 0.0201 | 0.0189 | 1.1 | 1.16 |
| 2 | 430 | 0 | 0.022 | 0.0199 | 0.0188 | 1.08 | 1.14 |
| 3 | 429 | 0 | 0.0231 | 0.0209 | 0.0196 | 1.07 | 1.14 |
| 4 | 430 | 0 | 0.0229 | 0.0208 | 0.0194 | 1.09 | 1.15 |
| 5 | 430 | 1 | 0.0302 | 0.0267 | 0.0246 | 1.11 | 1.19 |
| 6 | 429 | 4 | 0.0345 | 0.0299 | 0.0274 | 1.13 | 1.23 |
| 7 | 430 | 9 | 0.0399 | 0.0349 | 0.0318 | 1.14 | 1.25 |
| 8 | 429 | 26 | 0.0463 | 0.0395 | 0.0362 | 1.16 | 1.28 |
| 9 | 430 | 74 | 0.0575 | 0.0464 | 0.0408 | 1.22 | 1.38 |
| 10 | 430 | 274 | 0.0747 | 0.064 | 0.0562 | 1.19 | 1.36 |

#### Label-free held-out depth

| variant | n_splits | t_median | t_q25 | t_q75 | t_min | t_max | tmax_scanned | within1pct_lo_median | within1pct_hi_median |
|---|---|---|---|---|---|---|---|---|---|
| V0 | 8 | 9 | 9.00 | 10 | 9 | 10 | 20 | 7.5 | 12.0 |
| V1 | 8 | 7 | 6.00 |  7 | 6 |  8 | 20 | 6.0 |  9.0 |
| V2 | 8 | 6 | 5.00 |  6 | 5 |  7 | 20 | 4.5 |  8.5 |
| V3 | 8 | 6 | 5.75 |  6 | 5 |  7 | 20 | 5.0 |  8.0 |
| V4 | 8 | 6 | 6.00 |  7 | 6 |  7 | 20 | 5.5 |  8.5 |


### Tables for k562 (generated by analysis/robust_z/robust_z_summarize.R; EXPLORATORY)

#### Pooled enrichment (AUPR / base rate), pairwise score

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.06 | 1.06 | 1.06 | 1.08 | 1.05 |
| CORUM_all (base 0.0133) | 6.39 | 6.94 | 6.37 | 6.07 | 6.44 |
| KD_response (base 0.0787) | 2.42 | 2.61 | 2.53 | 2.55 | 2.60 |
| CORUM_noRP (base 0.0119) | 2.18 | 2.28 | 2.34 | 2.57 | 2.17 |
| CORUM_RP (base 0.7502) | 1.21 | 1.21 | 1.21 | 1.21 | 1.21 |
| CORUM_small (base 0.0033) | 2.16 | 2.05 | 2.08 | 2.23 | 1.93 |
| CORUM_large (base 0.1263) | 2.29 | 2.31 | 2.30 | 2.39 | 2.27 |

#### Pooled enrichment, attention t = 9

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.09 | 1.08 | 1.08 | 1.08 | 1.07 |
| CORUM_all (base 0.0133) | 3.51 | 3.26 | 3.44 | 3.78 | 2.99 |
| KD_response (base 0.0787) | 2.38 | 2.42 | 2.19 | 2.07 | 2.45 |
| CORUM_noRP (base 0.0119) | 2.05 | 1.96 | 2.09 | 2.27 | 1.86 |
| CORUM_RP (base 0.7502) | 1.18 | 1.18 | 1.18 | 1.20 | 1.19 |
| CORUM_small (base 0.0033) | 2.01 | 1.86 | 2.05 | 2.22 | 1.81 |
| CORUM_large (base 0.1263) | 1.93 | 1.80 | 1.89 | 2.08 | 1.75 |

#### Pooled enrichment, attention at each variant's own held-out depth

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB (base 0.1830) | 1.09 | 1.08 | 1.09 | 1.09 | 1.08 |
| CORUM_all (base 0.0133) | 3.51 | 3.11 | 3.26 | 3.48 | 3.03 |
| KD_response (base 0.0787) | 2.38 | 2.31 | 2.11 | 1.99 | 2.34 |
| CORUM_noRP (base 0.0119) | 2.05 | 1.97 | 2.06 | 2.23 | 1.94 |
| CORUM_RP (base 0.7502) | 1.18 | 1.18 | 1.18 | 1.19 | 1.18 |
| CORUM_small (base 0.0033) | 2.01 | 1.84 | 2.00 | 2.06 | 1.81 |
| CORUM_large (base 0.1263) | 1.93 | 1.85 | 1.93 | 2.13 | 1.81 |

#### Paired contrast, pairwise: variant - V0, raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V0x | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB | -0.000 [-0.000, +0.000] | -0.002 [-0.007, +0.002] | +0.002 [-0.003, +0.006] | +0.025 [+0.013, +0.037] | -0.008 [-0.013, -0.003] |
| CORUM_all | -0.001 [-0.001, -0.000] | +0.547 [+0.345, +0.749] | -0.018 [-0.267, +0.232] | -0.318 [-0.784, +0.147] | +0.052 [-0.224, +0.327] |
| KD_response | -0.000 [-0.000, -0.000] | +0.189 [+0.158, +0.220] | +0.102 [+0.075, +0.128] | +0.127 [+0.057, +0.198] | +0.178 [+0.135, +0.221] |
| CORUM_noRP | -0.000 [-0.000, +0.000] | +0.100 [+0.019, +0.182] | +0.158 [+0.050, +0.266] | +0.381 [+0.182, +0.581] | -0.019 [-0.155, +0.117] |
| CORUM_RP | -0.000 [-0.000, +0.000] | +0.005 [+0.001, +0.009] | +0.003 [-0.003, +0.009] | +0.004 [-0.004, +0.012] | +0.004 [+0.001, +0.008] |
| CORUM_small | +0.000 [-0.000, +0.001] | -0.110 [-0.379, +0.158] | -0.083 [-0.399, +0.234] | +0.069 [-0.355, +0.493] | -0.230 [-0.581, +0.121] |
| CORUM_large | -0.000 [-0.000, +0.000] | +0.023 [-0.000, +0.045] | +0.010 [-0.010, +0.030] | +0.097 [+0.046, +0.149] | -0.021 [-0.056, +0.013] |

#### Paired contrast, attention t=9: variant - V0, raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V0x | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB | +0.000 [+0.000, +0.000] | -0.017 [-0.027, -0.007] | -0.011 [-0.023, +0.001] | -0.013 [-0.040, +0.014] | -0.025 [-0.039, -0.012] |
| CORUM_all | +0.000 [-0.000, +0.000] | -0.250 [-0.539, +0.039] | -0.067 [-0.447, +0.313] | +0.272 [-0.351, +0.895] | -0.519 [-0.980, -0.059] |
| KD_response | -0.000 [-0.000, -0.000] | +0.040 [+0.002, +0.078] | -0.181 [-0.231, -0.130] | -0.302 [-0.386, -0.218] | +0.071 [+0.019, +0.123] |
| CORUM_noRP | +0.000 [-0.000, +0.000] | -0.088 [-0.233, +0.057] | +0.040 [-0.112, +0.193] | +0.216 [+0.031, +0.401] | -0.196 [-0.390, -0.001] |
| CORUM_RP | +0.000 [-0.000, +0.000] | -0.002 [-0.013, +0.009] | -0.004 [-0.020, +0.013] | +0.013 [-0.006, +0.032] | +0.001 [-0.007, +0.010] |
| CORUM_small | +0.000 [-0.001, +0.001] | -0.154 [-0.504, +0.195] | +0.039 [-0.303, +0.382] | +0.206 [-0.141, +0.552] | -0.204 [-0.641, +0.232] |
| CORUM_large | +0.000 [+0.000, +0.000] | -0.127 [-0.171, -0.082] | -0.039 [-0.093, +0.015] | +0.158 [+0.051, +0.265] | -0.180 [-0.252, -0.108] |

#### Paired contrast, attention t*: variant - V0, raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.010 [-0.019, -0.002] | -0.005 [-0.016, +0.006] | -0.003 [-0.027, +0.021] | -0.018 [-0.030, -0.006] |
| CORUM_all | -0.398 [-0.704, -0.093] | -0.249 [-0.629, +0.131] | -0.034 [-0.547, +0.479] | -0.484 [-0.845, -0.123] |
| KD_response | -0.066 [-0.099, -0.032] | -0.262 [-0.307, -0.218] | -0.381 [-0.458, -0.304] | -0.034 [-0.076, +0.008] |
| CORUM_noRP | -0.078 [-0.199, +0.042] | +0.014 [-0.107, +0.135] | +0.177 [-0.005, +0.359] | -0.115 [-0.280, +0.050] |
| CORUM_RP | -0.003 [-0.012, +0.006] | -0.004 [-0.018, +0.010] | +0.007 [-0.011, +0.025] | -0.000 [-0.008, +0.007] |
| CORUM_small | -0.174 [-0.510, +0.162] | -0.014 [-0.244, +0.217] | +0.051 [-0.267, +0.369] | -0.205 [-0.598, +0.188] |
| CORUM_large | -0.079 [-0.120, -0.037] | +0.002 [-0.050, +0.054] | +0.203 [+0.098, +0.309] | -0.111 [-0.169, -0.054] |

#### Paired contrast, attention t=9 minus CLR(V0 pairwise), raw (difference in enrichment [95% delete-half gene jackknife])

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB | +0.035 [+0.007, +0.062] | +0.018 [-0.011, +0.046] | +0.024 [-0.004, +0.052] | +0.022 [-0.013, +0.057] | +0.010 [-0.019, +0.038] |
| CORUM_all | -2.880 [-4.316, -1.444] | -3.130 [-4.733, -1.526] | -2.947 [-4.565, -1.329] | -2.608 [-4.283, -0.932] | -3.399 [-5.060, -1.738] |
| KD_response | -0.049 [-0.153, +0.055] | -0.009 [-0.118, +0.100] | -0.230 [-0.345, -0.114] | -0.351 [-0.481, -0.221] | +0.022 [-0.091, +0.135] |
| CORUM_noRP | -0.134 [-0.350, +0.083] | -0.222 [-0.467, +0.024] | -0.093 [-0.345, +0.159] | +0.082 [-0.179, +0.343] | -0.329 [-0.572, -0.086] |
| CORUM_RP | -0.025 [-0.052, +0.003] | -0.027 [-0.058, +0.004] | -0.028 [-0.061, +0.005] | -0.012 [-0.039, +0.016] | -0.023 [-0.052, +0.006] |
| CORUM_small | -0.150 [-1.093, +0.794] | -0.304 [-1.300, +0.692] | -0.110 [-1.179, +0.958] | +0.056 [-0.989, +1.100] | -0.354 [-1.322, +0.615] |
| CORUM_large | -0.365 [-0.525, -0.205] | -0.491 [-0.664, -0.319] | -0.404 [-0.581, -0.226] | -0.207 [-0.383, -0.030] | -0.544 [-0.723, -0.366] |

#### Paired contrast, pairwise: variant - V0, hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.004 [-0.009, +0.000] | +0.000 [-0.004, +0.005] | +0.021 [+0.007, +0.034] | -0.010 [-0.015, -0.004] |
| CORUM_all | +0.014 [-0.178, +0.206] | -0.210 [-0.478, +0.058] | -0.220 [-0.624, +0.184] | -0.256 [-0.505, -0.006] |
| KD_response | +0.050 [+0.025, +0.075] | +0.013 [-0.009, +0.034] | +0.147 [+0.089, +0.206] | +0.044 [+0.011, +0.076] |
| CORUM_noRP | +0.023 [-0.039, +0.084] | +0.066 [+0.007, +0.124] | +0.244 [+0.132, +0.356] | -0.033 [-0.114, +0.048] |
| CORUM_RP | +0.004 [-0.000, +0.007] | +0.003 [-0.003, +0.009] | +0.004 [-0.004, +0.013] | +0.004 [-0.000, +0.008] |
| CORUM_small | -0.090 [-0.268, +0.087] | +0.059 [-0.114, +0.232] | +0.103 [-0.186, +0.392] | -0.178 [-0.476, +0.121] |
| CORUM_large | -0.034 [-0.074, +0.006] | -0.065 [-0.102, -0.028] | -0.042 [-0.128, +0.043] | -0.085 [-0.150, -0.019] |

#### Paired contrast, attention t=9: variant - V0, hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.013 [-0.023, -0.003] | -0.008 [-0.020, +0.005] | -0.014 [-0.040, +0.011] | -0.026 [-0.042, -0.011] |
| CORUM_all | -0.130 [-0.276, +0.016] | -0.028 [-0.231, +0.175] | +0.258 [-0.008, +0.524] | -0.228 [-0.392, -0.064] |
| KD_response | -0.052 [-0.079, -0.024] | -0.153 [-0.191, -0.116] | -0.112 [-0.174, -0.051] | -0.032 [-0.070, +0.007] |
| CORUM_noRP | -0.040 [-0.109, +0.030] | +0.038 [-0.045, +0.121] | +0.275 [+0.141, +0.409] | -0.089 [-0.197, +0.020] |
| CORUM_RP | -0.004 [-0.016, +0.007] | -0.005 [-0.022, +0.012] | +0.011 [-0.008, +0.030] | -0.001 [-0.010, +0.009] |
| CORUM_small | -0.051 [-0.215, +0.113] | +0.030 [-0.171, +0.231] | +0.165 [-0.193, +0.523] | -0.098 [-0.398, +0.201] |
| CORUM_large | -0.129 [-0.186, -0.072] | -0.059 [-0.129, +0.011] | +0.048 [-0.080, +0.175] | -0.186 [-0.262, -0.110] |

#### Paired contrast, attention t*: variant - V0, hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.012 [-0.020, -0.003] | -0.007 [-0.017, +0.004] | -0.005 [-0.027, +0.016] | -0.022 [-0.034, -0.009] |
| CORUM_all | -0.125 [-0.247, -0.002] | -0.012 [-0.172, +0.148] | +0.313 [+0.104, +0.523] | -0.161 [-0.305, -0.016] |
| KD_response | -0.030 [-0.054, -0.005] | -0.112 [-0.146, -0.079] | -0.039 [-0.093, +0.015] | -0.010 [-0.042, +0.022] |
| CORUM_noRP | -0.020 [-0.078, +0.038] | +0.052 [-0.020, +0.125] | +0.296 [+0.180, +0.412] | -0.043 [-0.130, +0.044] |
| CORUM_RP | -0.005 [-0.014, +0.004] | -0.005 [-0.019, +0.009] | +0.007 [-0.011, +0.026] | -0.002 [-0.011, +0.006] |
| CORUM_small | -0.034 [-0.149, +0.081] | +0.036 [-0.063, +0.134] | +0.133 [-0.026, +0.292] | -0.067 [-0.238, +0.104] |
| CORUM_large | -0.109 [-0.158, -0.060] | -0.045 [-0.107, +0.016] | +0.083 [-0.029, +0.195] | -0.138 [-0.201, -0.076] |

#### Paired contrast, attention t=9 minus CLR(V0 pairwise), hub-strat (pi of V0 attention) (difference in enrichment [95% delete-half gene jackknife])

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB | +0.031 [+0.009, +0.053] | +0.018 [-0.004, +0.040] | +0.023 [+0.001, +0.046] | +0.017 [-0.016, +0.049] | +0.005 [-0.020, +0.029] |
| CORUM_all | -1.325 [-2.153, -0.497] | -1.455 [-2.353, -0.557] | -1.353 [-2.268, -0.437] | -1.067 [-2.003, -0.131] | -1.553 [-2.472, -0.634] |
| KD_response | -0.311 [-0.395, -0.226] | -0.363 [-0.451, -0.274] | -0.464 [-0.558, -0.371] | -0.423 [-0.534, -0.312] | -0.343 [-0.433, -0.252] |
| CORUM_noRP | -0.129 [-0.269, +0.012] | -0.168 [-0.313, -0.024] | -0.091 [-0.246, +0.064] | +0.146 [-0.041, +0.334] | -0.218 [-0.374, -0.061] |
| CORUM_RP | -0.026 [-0.055, +0.003] | -0.030 [-0.063, +0.002] | -0.031 [-0.067, +0.004] | -0.015 [-0.044, +0.014] | -0.027 [-0.058, +0.004] |
| CORUM_small | -0.171 [-0.664, +0.322] | -0.222 [-0.652, +0.209] | -0.141 [-0.617, +0.335] | -0.006 [-0.572, +0.560] | -0.269 [-0.713, +0.175] |
| CORUM_large | -0.248 [-0.430, -0.066] | -0.377 [-0.578, -0.176] | -0.307 [-0.514, -0.100] | -0.200 [-0.421, +0.020] | -0.434 [-0.632, -0.236] |

#### Paired contrast, pairwise: variant - V0, hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | +0.003 [-0.002, +0.009] | +0.005 [-0.001, +0.010] | +0.001 [-0.007, +0.009] | +0.001 [-0.005, +0.007] |
| CORUM_all | +0.483 [+0.234, +0.733] | +0.139 [-0.115, +0.394] | +0.430 [+0.062, +0.799] | +0.195 [-0.033, +0.423] |
| KD_response | +0.177 [+0.151, +0.202] | +0.080 [+0.055, +0.106] | +0.173 [+0.129, +0.218] | +0.205 [+0.173, +0.236] |
| CORUM_noRP | +0.167 [+0.086, +0.249] | +0.131 [+0.048, +0.214] | +0.247 [+0.134, +0.360] | +0.166 [+0.061, +0.271] |
| CORUM_RP | +0.004 [-0.001, +0.008] | +0.002 [-0.004, +0.008] | +0.002 [-0.006, +0.009] | +0.003 [-0.001, +0.007] |
| CORUM_small | -0.037 [-0.297, +0.223] | -0.079 [-0.398, +0.240] | +0.040 [-0.257, +0.337] | -0.082 [-0.374, +0.210] |
| CORUM_large | +0.018 [-0.009, +0.046] | -0.006 [-0.035, +0.024] | +0.027 [-0.017, +0.071] | -0.022 [-0.060, +0.016] |

#### Paired contrast, attention t=9: variant - V0, hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.014 [-0.025, -0.003] | -0.013 [-0.025, -0.001] | -0.027 [-0.049, -0.004] | -0.021 [-0.034, -0.007] |
| CORUM_all | -0.159 [-0.393, +0.075] | -0.004 [-0.315, +0.307] | +0.543 [-0.036, +1.122] | -0.417 [-0.668, -0.166] |
| KD_response | +0.021 [-0.016, +0.058] | -0.154 [-0.199, -0.110] | -0.145 [-0.217, -0.074] | +0.053 [+0.006, +0.100] |
| CORUM_noRP | -0.121 [-0.239, -0.002] | -0.055 [-0.191, +0.080] | +0.098 [-0.102, +0.299] | -0.189 [-0.357, -0.022] |
| CORUM_RP | -0.004 [-0.017, +0.008] | -0.004 [-0.021, +0.013] | +0.015 [-0.003, +0.033] | -0.001 [-0.011, +0.009] |
| CORUM_small | -0.108 [-0.315, +0.099] | +0.005 [-0.329, +0.339] | +0.148 [-0.269, +0.565] | -0.142 [-0.391, +0.107] |
| CORUM_large | -0.128 [-0.177, -0.080] | -0.054 [-0.112, +0.004] | +0.142 [+0.033, +0.250] | -0.187 [-0.246, -0.129] |

#### Paired contrast, attention t*: variant - V0, hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.009 [-0.018, -0.000] | -0.009 [-0.020, +0.002] | -0.018 [-0.037, +0.001] | -0.015 [-0.027, -0.004] |
| CORUM_all | -0.292 [-0.507, -0.077] | -0.131 [-0.406, +0.145] | +0.219 [-0.164, +0.602] | -0.374 [-0.624, -0.124] |
| KD_response | -0.060 [-0.092, -0.027] | -0.216 [-0.256, -0.176] | -0.214 [-0.276, -0.152] | -0.027 [-0.067, +0.013] |
| CORUM_noRP | -0.129 [-0.257, -0.001] | -0.076 [-0.209, +0.057] | +0.073 [-0.130, +0.277] | -0.149 [-0.306, +0.007] |
| CORUM_RP | -0.004 [-0.014, +0.006] | -0.004 [-0.018, +0.011] | +0.011 [-0.007, +0.029] | -0.002 [-0.010, +0.007] |
| CORUM_small | -0.152 [-0.397, +0.094] | -0.056 [-0.286, +0.174] | +0.025 [-0.287, +0.338] | -0.162 [-0.449, +0.124] |
| CORUM_large | -0.083 [-0.126, -0.040] | -0.010 [-0.065, +0.045] | +0.178 [+0.077, +0.279] | -0.116 [-0.167, -0.065] |

#### Paired contrast, attention t=9 minus CLR(V0 pairwise), hub-strat (MI row mean) (difference in enrichment [95% delete-half gene jackknife])

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB | +0.016 [-0.010, +0.042] | +0.002 [-0.026, +0.030] | +0.003 [-0.025, +0.030] | -0.010 [-0.040, +0.020] | -0.005 [-0.033, +0.023] |
| CORUM_all | -2.066 [-3.288, -0.843] | -2.225 [-3.522, -0.928] | -2.070 [-3.372, -0.767] | -1.523 [-2.780, -0.265] | -2.483 [-3.848, -1.118] |
| KD_response | -0.006 [-0.087, +0.074] | +0.015 [-0.073, +0.102] | -0.161 [-0.247, -0.074] | -0.152 [-0.243, -0.060] | +0.047 [-0.044, +0.138] |
| CORUM_noRP | -0.119 [-0.329, +0.091] | -0.240 [-0.477, -0.003] | -0.174 [-0.412, +0.063] | -0.021 [-0.305, +0.263] | -0.308 [-0.551, -0.065] |
| CORUM_RP | -0.031 [-0.062, +0.000] | -0.035 [-0.069, -0.001] | -0.035 [-0.071, +0.001] | -0.016 [-0.045, +0.013] | -0.032 [-0.064, +0.000] |
| CORUM_small | -0.241 [-1.090, +0.608] | -0.349 [-1.251, +0.553] | -0.236 [-1.263, +0.791] | -0.093 [-1.014, +0.829] | -0.383 [-1.244, +0.478] |
| CORUM_large | -0.451 [-0.637, -0.266] | -0.580 [-0.788, -0.371] | -0.505 [-0.711, -0.300] | -0.310 [-0.491, -0.128] | -0.639 [-0.845, -0.433] |

#### Macro: mean per-regulon AP, pairwise

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB (n=105) | 0.011 | 0.011 | 0.011 | 0.012 | 0.010 |
| CORUM_all (n=499) | 0.012 | 0.012 | 0.012 | 0.012 | 0.011 |
| KD_response (n=566) | 0.162 | 0.173 | 0.159 | 0.162 | 0.172 |
| CORUM_noRP (n=492) | 0.006 | 0.006 | 0.006 | 0.006 | 0.005 |
| CORUM_RP (n=7) | 0.429 | 0.453 | 0.422 | 0.437 | 0.435 |
| CORUM_small (n=185) | 0.003 | 0.002 | 0.003 | 0.003 | 0.002 |
| CORUM_large (n=28) | 0.135 | 0.139 | 0.132 | 0.136 | 0.133 |

#### Macro: mean per-regulon AP, attention t = 9

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB (n=105) | 0.012 | 0.011 | 0.011 | 0.011 | 0.011 |
| CORUM_all (n=499) | 0.007 | 0.006 | 0.006 | 0.008 | 0.005 |
| KD_response (n=566) | 0.065 | 0.059 | 0.058 | 0.069 | 0.059 |
| CORUM_noRP (n=492) | 0.005 | 0.004 | 0.004 | 0.005 | 0.004 |
| CORUM_RP (n=7) | 0.142 | 0.112 | 0.115 | 0.227 | 0.122 |
| CORUM_small (n=185) | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| CORUM_large (n=28) | 0.064 | 0.051 | 0.053 | 0.079 | 0.051 |

#### Macro contrast, pairwise: variant - V0 (mean difference in per-regulon AP; regulon bootstrap and delete-20% gene jackknife intervals)

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.000; reg [-0.001, -0.000]; gene [-0.003, +0.002] | -0.000; reg [-0.001, +0.000]; gene [-0.002, +0.001] | +0.001; reg [+0.000, +0.001]; gene [-0.002, +0.004] | -0.001; reg [-0.001, -0.000]; gene [-0.004, +0.002] |
| CORUM_all | +0.000; reg [-0.000, +0.000]; gene [-0.000, +0.001] | -0.000; reg [-0.001, +0.000]; gene [-0.001, +0.001] | +0.000; reg [-0.000, +0.000]; gene [-0.002, +0.002] | -0.001; reg [-0.001, -0.000]; gene [-0.002, +0.000] |
| KD_response | +0.011; reg [+0.010, +0.012]; gene [+0.009, +0.014] | -0.003; reg [-0.004, -0.002]; gene [-0.006, -0.000] | +0.001; reg [-0.001, +0.002]; gene [-0.005, +0.006] | +0.011; reg [+0.009, +0.012]; gene [+0.008, +0.014] |
| CORUM_noRP | -0.000; reg [-0.000, -0.000]; gene [-0.001, +0.000] | -0.000; reg [-0.000, +0.000]; gene [-0.001, +0.001] | -0.000; reg [-0.000, +0.000]; gene [-0.001, +0.001] | -0.001; reg [-0.001, -0.000]; gene [-0.002, +0.000] |
| CORUM_RP | +0.024; reg [+0.019, +0.029]; gene [+0.009, +0.038] | -0.007; reg [-0.022, +0.006]; gene [-0.040, +0.026] | +0.008; reg [-0.008, +0.025]; gene [-0.031, +0.048] | +0.006; reg [+0.001, +0.011]; gene [-0.011, +0.023] |
| CORUM_small | -0.000; reg [-0.000, -0.000]; gene [-0.001, +0.000] | -0.000; reg [-0.000, +0.000]; gene [-0.000, +0.000] | +0.000; reg [-0.000, +0.000]; gene [-0.000, +0.000] | -0.000; reg [-0.001, -0.000]; gene [-0.001, +0.000] |
| CORUM_large | +0.004; reg [-0.000, +0.009]; gene [+0.000, +0.008] | -0.003; reg [-0.008, +0.000]; gene [-0.011, +0.005] | +0.001; reg [-0.004, +0.007]; gene [-0.009, +0.012] | -0.002; reg [-0.006, +0.001]; gene [-0.007, +0.003] |

#### Macro contrast, attention t=9: variant - V0 (mean difference in per-regulon AP; regulon bootstrap and delete-20% gene jackknife intervals)

| truth | V1 | V2 | V3 | V4 |
|---|---|---|---|---|
| DoRothEA_AB | -0.000; reg [-0.001, -0.000]; gene [-0.002, +0.000] | -0.000; reg [-0.001, +0.000]; gene [-0.002, +0.001] | -0.000; reg [-0.001, +0.000]; gene [-0.002, +0.001] | -0.001; reg [-0.001, -0.000]; gene [-0.002, +0.000] |
| CORUM_all | -0.001; reg [-0.002, -0.001]; gene [-0.002, -0.000] | -0.001; reg [-0.001, -0.000]; gene [-0.002, +0.000] | +0.001; reg [-0.000, +0.002]; gene [-0.001, +0.004] | -0.001; reg [-0.002, -0.001]; gene [-0.003, -0.000] |
| KD_response | -0.006; reg [-0.008, -0.004]; gene [-0.011, -0.001] | -0.007; reg [-0.009, -0.004]; gene [-0.015, +0.001] | +0.004; reg [+0.001, +0.008]; gene [-0.006, +0.014] | -0.006; reg [-0.008, -0.004]; gene [-0.011, +0.000] |
| CORUM_noRP | -0.001; reg [-0.001, -0.000]; gene [-0.002, +0.000] | -0.000; reg [-0.001, +0.000]; gene [-0.001, +0.001] | -0.000; reg [-0.001, +0.000]; gene [-0.002, +0.002] | -0.001; reg [-0.002, -0.001]; gene [-0.003, +0.000] |
| CORUM_RP | -0.030; reg [-0.042, -0.020]; gene [-0.049, -0.010] | -0.027; reg [-0.045, -0.005]; gene [-0.061, +0.008] | +0.085; reg [+0.044, +0.133]; gene [+0.042, +0.128] | -0.020; reg [-0.029, -0.012]; gene [-0.036, -0.005] |
| CORUM_small | -0.001; reg [-0.001, -0.000]; gene [-0.001, +0.001] | -0.000; reg [-0.001, +0.000]; gene [-0.001, +0.001] | -0.000; reg [-0.001, +0.000]; gene [-0.001, +0.001] | -0.001; reg [-0.001, -0.000]; gene [-0.001, +0.001] |
| CORUM_large | -0.012; reg [-0.020, -0.006]; gene [-0.020, -0.005] | -0.010; reg [-0.018, -0.003]; gene [-0.020, -0.000] | +0.016; reg [-0.002, +0.038]; gene [+0.002, +0.029] | -0.012; reg [-0.021, -0.005]; gene [-0.020, -0.004] |

#### Macro contrast, attention t=9 minus CLR(V0 pairwise) (mean difference in per-regulon AP; regulon bootstrap and delete-20% gene jackknife intervals)

| truth | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| DoRothEA_AB | +0.001; reg [-0.000, +0.001]; gene [-0.002, +0.004] | +0.000; reg [-0.001, +0.001]; gene [-0.003, +0.003] | +0.000; reg [-0.001, +0.001]; gene [-0.003, +0.003] | +0.000; reg [-0.001, +0.001]; gene [-0.003, +0.004] | -0.000; reg [-0.001, +0.001]; gene [-0.004, +0.003] |
| CORUM_all | -0.005; reg [-0.009, -0.002]; gene [-0.011, -0.006] | -0.006; reg [-0.010, -0.003]; gene [-0.013, -0.007] | -0.006; reg [-0.010, -0.003]; gene [-0.012, -0.006] | -0.004; reg [-0.007, -0.002]; gene [-0.009, -0.004] | -0.006; reg [-0.010, -0.003]; gene [-0.013, -0.007] |
| KD_response | -0.097; reg [-0.106, -0.088]; gene [-0.114, -0.081] | -0.103; reg [-0.113, -0.094]; gene [-0.121, -0.085] | -0.104; reg [-0.114, -0.094]; gene [-0.123, -0.085] | -0.093; reg [-0.100, -0.086]; gene [-0.110, -0.075] | -0.103; reg [-0.112, -0.093]; gene [-0.120, -0.086] |
| CORUM_noRP | -0.001; reg [-0.002, -0.001]; gene [-0.003, +0.000] | -0.002; reg [-0.002, -0.001]; gene [-0.004, -0.000] | -0.001; reg [-0.002, -0.001]; gene [-0.004, +0.000] | -0.001; reg [-0.002, -0.001]; gene [-0.004, +0.001] | -0.002; reg [-0.003, -0.001]; gene [-0.005, -0.001] |
| CORUM_RP | -0.287; reg [-0.379, -0.189]; gene [-0.375, -0.199] | -0.317; reg [-0.417, -0.217]; gene [-0.413, -0.220] | -0.314; reg [-0.407, -0.213]; gene [-0.408, -0.219] | -0.202; reg [-0.257, -0.135]; gene [-0.279, -0.124] | -0.307; reg [-0.399, -0.207]; gene [-0.403, -0.212] |
| CORUM_small | -0.000; reg [-0.001, +0.000]; gene [-0.001, +0.001] | -0.001; reg [-0.002, -0.000]; gene [-0.001, +0.000] | -0.001; reg [-0.001, -0.000]; gene [-0.001, +0.001] | -0.001; reg [-0.001, -0.000]; gene [-0.001, +0.001] | -0.001; reg [-0.002, -0.000]; gene [-0.001, +0.000] |
| CORUM_large | -0.071; reg [-0.128, -0.023]; gene [-0.093, -0.050] | -0.084; reg [-0.146, -0.033]; gene [-0.107, -0.060] | -0.082; reg [-0.141, -0.031]; gene [-0.105, -0.058] | -0.056; reg [-0.094, -0.023]; gene [-0.076, -0.036] | -0.084; reg [-0.141, -0.034]; gene [-0.107, -0.060] |

#### Hub drift (top pairs of each score against two connectedness rankings; chance share = 0.0975, chance mean rank = 0.5)

| variant | kind | stratifier | share_touching_top5pct | mean_endpoint_rank | frac_positive_score |
|---|---|---|---|---|---|
| V0 | pair | att_pi | 0.336 | 0.768 | 0.543 |
| V0 | pair | mi_rowmean | 0.419 | 0.742 | 0.543 |
| V0 | att9 | att_pi | 0.608 | 0.717 | 1.000 |
| V0 | att9 | mi_rowmean | 0.413 | 0.649 | 1.000 |
| V0 | attT | att_pi | 0.608 | 0.717 | 1.000 |
| V0 | attT | mi_rowmean | 0.413 | 0.649 | 1.000 |
| V1 | pair | att_pi | 0.373 | 0.776 | 0.617 |
| V1 | pair | mi_rowmean | 0.507 | 0.783 | 0.617 |
| V1 | att9 | att_pi | 0.711 | 0.724 | 1.000 |
| V1 | att9 | mi_rowmean | 0.433 | 0.645 | 1.000 |
| V1 | attT | att_pi | 0.571 | 0.701 | 1.000 |
| V1 | attT | mi_rowmean | 0.393 | 0.634 | 1.000 |
| V2 | pair | att_pi | 0.368 | 0.781 | 0.617 |
| V2 | pair | mi_rowmean | 0.419 | 0.750 | 0.617 |
| V2 | att9 | att_pi | 0.636 | 0.721 | 1.000 |
| V2 | att9 | mi_rowmean | 0.331 | 0.621 | 1.000 |
| V2 | attT | att_pi | 0.502 | 0.694 | 1.000 |
| V2 | attT | mi_rowmean | 0.305 | 0.613 | 1.000 |
| V3 | pair | att_pi | 0.328 | 0.761 | 0.500 |
| V3 | pair | mi_rowmean | 0.259 | 0.627 | 0.500 |
| V3 | att9 | att_pi | 0.543 | 0.707 | 1.000 |
| V3 | att9 | mi_rowmean | 0.171 | 0.526 | 1.000 |
| V3 | attT | att_pi | 0.382 | 0.661 | 1.000 |
| V3 | attT | mi_rowmean | 0.163 | 0.526 | 1.000 |
| V4 | pair | att_pi | 0.347 | 0.752 | 0.617 |
| V4 | pair | mi_rowmean | 0.564 | 0.814 | 0.617 |
| V4 | att9 | att_pi | 0.685 | 0.723 | 1.000 |
| V4 | att9 | mi_rowmean | 0.455 | 0.649 | 1.000 |
| V4 | attT | att_pi | 0.559 | 0.701 | 1.000 |
| V4 | attT | mi_rowmean | 0.417 | 0.639 | 1.000 |
| V0x | pair | att_pi | 0.336 | 0.767 | 0.543 |
| V0x | pair | mi_rowmean | 0.419 | 0.742 | 0.543 |
| V0x | att9 | att_pi | 0.607 | 0.717 | 1.000 |
| V0x | att9 | mi_rowmean | 0.413 | 0.649 | 1.000 |

#### Zero-scale guard (rows with scale 0)

| variant | rows |
|---|---|
| V0 | 0 |
| V0x | 0 |
| V1 | 0 |
| V2 | 0 |
| V3 | 0 |
| V4 | 0 |

#### V0 against the project's baselines run

| truth | V0_pair_aupr | baselines_clr_aupr | V0_att_aupr | baselines_attention_aupr | max_abs_diff |
|---|---|---|---|---|---|
| DoRothEA_AB | 0.19394986 | 0.19394986 | 0.20032497 | 0.20032497 | 3.885781e-16 |
| CORUM_all | 0.08471351 | 0.08471351 | 0.04653961 | 0.04653961 | 2.775558e-17 |
| KD_response | 0.19083871 | 0.19083871 | 0.18697540 | 0.18697540 | 3.608225e-16 |

#### Row scale versus number of strong partners (Spearman)

| row_quantity | versus | spearman | n_genes |
|---|---|---|---|
| mean_diag0 | n_strong_global | +0.797 | 5671 |
| mean_diag0 | n_z4_gt3 | +0.645 | 5671 |
| mean_diag0 | module_size | +0.106 | 2260 |
| sd_diag0 | n_strong_global | +0.890 | 5671 |
| sd_diag0 | n_z4_gt3 | +0.861 | 5671 |
| sd_diag0 | module_size | +0.168 | 2260 |
| median | n_strong_global | +0.650 | 5671 |
| median | n_z4_gt3 | +0.432 | 5671 |
| median | module_size | +0.063 | 2260 |
| iqr_scale | n_strong_global | +0.846 | 5671 |
| iqr_scale | n_z4_gt3 | +0.764 | 5671 |
| iqr_scale | module_size | +0.131 | 2260 |
| mad_scale | n_strong_global | +0.831 | 5671 |
| mad_scale | n_z4_gt3 | +0.718 | 5671 |
| mad_scale | module_size | +0.124 | 2260 |
| sd_over_iqr | n_strong_global | +0.772 | 5671 |
| sd_over_iqr | n_z4_gt3 | +0.941 | 5671 |
| sd_over_iqr | module_size | +0.244 | 2260 |
| sd_over_mad | n_strong_global | +0.820 | 5671 |
| sd_over_mad | n_z4_gt3 | +0.982 | 5671 |
| sd_over_mad | module_size | +0.221 | 2260 |
| mean_minus_median_in_iqr | n_strong_global | +0.826 | 5671 |
| mean_minus_median_in_iqr | n_z4_gt3 | +0.985 | 5671 |
| mean_minus_median_in_iqr | module_size | +0.191 | 2260 |

#### Row scale by decile of the number of globally strong partners

| decile | n_genes | n_strong_global_median | sd_median | iqr_scale_median | mad_scale_median | sd_over_iqr_median | sd_over_mad_median |
|---|---|---|---|---|---|---|---|
| 1 | 568 | 0 | 0.00792 | 0.00754 | 0.00742 | 1.07 | 1.09 |
| 2 | 567 | 0 | 0.00794 | 0.00753 | 0.00743 | 1.07 | 1.09 |
| 3 | 567 | 0 | 0.00776 | 0.00746 | 0.00739 | 1.05 | 1.06 |
| 4 | 567 | 0 | 0.00767 | 0.00737 | 0.00724 | 1.06 | 1.07 |
| 5 | 567 | 0 | 0.00755 | 0.00727 | 0.0072 | 1.05 | 1.06 |
| 6 | 567 | 1 | 0.0104 | 0.00879 | 0.00842 | 1.15 | 1.2 |
| 7 | 567 | 6 | 0.0133 | 0.011 | 0.0101 | 1.22 | 1.3 |
| 8 | 567 | 29 | 0.016 | 0.0127 | 0.0115 | 1.26 | 1.39 |
| 9 | 567 | 101 | 0.02 | 0.0156 | 0.0136 | 1.31 | 1.49 |
| 10 | 567 | 366 | 0.0293 | 0.0219 | 0.018 | 1.35 | 1.64 |

#### Label-free held-out depth

| variant | n_splits | t_median | t_q25 | t_q75 | t_min | t_max | tmax_scanned | within1pct_lo_median | within1pct_hi_median |
|---|---|---|---|---|---|---|---|---|---|
| V0 | 6 | 9 | 8.25 | 9 | 8 | 9 | 20 | 7.0 | 10 |
| V1 | 6 | 6 | 6.00 | 6 | 5 | 6 | 20 | 5.0 |  7 |
| V2 | 6 | 6 | 6.00 | 6 | 5 | 6 | 20 | 5.0 |  7 |
| V3 | 6 | 5 | 5.00 | 5 | 5 | 5 | 20 | 4.5 |  7 |
| V4 | 6 | 6 | 6.00 | 6 | 5 | 6 | 20 | 5.0 |  7 |

