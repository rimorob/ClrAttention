# EcoMAC re-score on arrays without an M3D counterpart (sensitivity analysis)

**Status: sensitivity analysis, explicitly NOT preregistered.** Approved by the PI after the overlap audit (`docs/plans/ecomac_m3d_overlap.md`). Date 2026-11-01. It reuses a consumed truth (the RegulonDB files that Table 2 already uses; no new truth was opened), so it cannot serve as a held-out test, and it must be disclosed as post hoc. The subsets were chosen by **similarity to M3D (expression profiles only), not by outcome**: the overlap calls in `results/ecomac_m3d_overlap/per_array_best_match.csv` were made before any of these results existed, with no truth file open; the headline call is "centred rule AND mutual nearest neighbour" (828 called, 1,370 not called).

## What was run

The headline EcoMAC pipeline, unchanged except for the array subset: `analysis/baselines.R --dataset ecomac` (the script behind the 0875/0876c/0882/0884 jobs and Table 2): `data/ecomac/EcoMAC.tsv`, ComBat by platform (`sva`), 4,189 genes, MI with hg bins, Stouffer-CLR, K10 softmax (cap 50, alpha 0.5) at t = 9 as transferred, pooled co-membership AUPR on RegulonDB (strong/confirmed "SC" and all evidence), hub stratification by pi_i + pi_j (10 strata), paired delete-half gene jackknife (100 halves, seed 20261016, the same halves for every subset: all subsets keep the same 4,189 genes). The only change to the code is an option `--array_subset <file>` (default empty = all arrays, behaviour unchanged); ComBat is re-fitted by platform within the subset and skipped when a single platform remains. Per-regulon (macro) tables were not produced (`--macro 0`). Interval method as in Table 2: point +/- 1.96 sd of the replicate differences. The ratio interval uses the same replicates (ratio of the replicate AUPRs). Intervals resample genes with networks held fixed; they do not resample arrays.

**Gate (passed).** Subset A (all 2,198 arrays) reproduces the existing `results/baselines_ecomac/metrics.csv` for every method and truth (CLR, |r|, attention, GENIE3 [reused scores], hub, hub-stratified variants, and all interval endpoints, TF-node numbers): max absolute difference 0 over all columns (tolerance was 1e-6).

Subsets (`analysis/ecomac_nonoverlap_subsets.R`, `results/ecomac_nonoverlap/subset_composition.csv`):

| Subset | Definition | Arrays | GPL199 (Antisense) | GPL3154 (Genome 2.0) | No GEO record | Post-build (GSE >= 20000) | Series |
|---|---|---|---|---|---|---|---|
| A | all | 2,198 | 844 | 735 | 619 | 518 | 127 |
| B | not called (no M3D counterpart) | 1,370 | 232 | 735 | 403 | 518 | 101 |
| C | called (828, contrast) | 828 | 612 | 0 | 216 | 0 | 35 |
| B-g3154 | B, GPL3154 only | 735 | 0 | 735 | 0 | 381 | 85 |
| B-mix | B, GPL199 uncalled + no-GEO | 635 | 232 | 0 | 403 | 137 | 16 |
| R1370 | (extra reference) random 1,370 of A, seed 20261101 | 1,370 | 503 | 462 | 405 | 313 | 122 |

B is 54% Genome 2.0 (GPL3154), 29% no-GEO, 17% GPL199; C has no Genome 2.0 array at all. B therefore mixes platform and overlap, which is why B-g3154 and B-mix are given. R1370 is not in the request: it shows how much AUPR moves when 38% of arrays are dropped at random, so that the B and C differences can be read against it.

## Results at t = 9 (as transferred)

Pooled AUPR; difference and ratio are attention minus / over CLR with delete-half gene jackknife 95% intervals. Base rates: SC 0.0939, all evidence 0.1210 (same truth pairs for every subset).

**Strong/confirmed (SC) truth**

| Subset | |r| | CLR | Attention | GENIE3 | Attention - CLR | Ratio |
|---|---|---|---|---|---|---|
| A (gate) | 0.110 | 0.130 | 0.158 | 0.141 | +0.0278 [+0.0213, +0.0343] | 1.21 [1.17, 1.26] |
| **B** | 0.106 | 0.129 | 0.157 | not run | **+0.0285 [+0.0218, +0.0353]** | **1.22 [1.17, 1.27]** |
| C | 0.117 | 0.134 | 0.157 | not run | +0.0228 [+0.0167, +0.0289] | 1.17 [1.13, 1.21] |
| B-g3154 | 0.116 | 0.136 | 0.159 | not run | +0.0230 [+0.0170, +0.0291] | 1.17 [1.13, 1.21] |
| B-mix | 0.095 | 0.114 | 0.133 | not run | +0.0197 [+0.0143, +0.0251] | 1.17 [1.13, 1.22] |
| R1370 (reference) | 0.110 | 0.130 | 0.157 | not run | +0.0273 [+0.0209, +0.0337] | 1.21 [1.16, 1.26] |

**Hub-controlled (attention vs CLR, each within-stratum percentile of pi_i + pi_j), SC**

| Subset | CLR | Attention | Attention - CLR | Ratio |
|---|---|---|---|---|
| A | 0.128 | 0.148 | +0.0195 [+0.0140, +0.0251] | 1.15 [1.11, 1.19] |
| **B** | 0.127 | 0.146 | **+0.0197 [+0.0141, +0.0254]** | 1.16 [1.11, 1.20] |
| C | 0.129 | 0.144 | +0.0150 [+0.0099, +0.0200] | 1.12 [1.08, 1.15] |
| B-g3154 | 0.131 | 0.148 | +0.0167 [+0.0114, +0.0221] | 1.13 [1.09, 1.17] |
| B-mix | 0.114 | 0.125 | +0.0108 [+0.0064, +0.0152] | 1.09 [1.06, 1.13] |
| R1370 | 0.129 | 0.147 | +0.0184 [+0.0128, +0.0240] | 1.14 [1.10, 1.19] |

**All-evidence truth**

| Subset | |r| | CLR | Attention | GENIE3 | Attention - CLR | Ratio |
|---|---|---|---|---|---|---|
| A (gate) | 0.144 | 0.161 | 0.178 | 0.168 | +0.0164 [+0.0118, +0.0211] | 1.10 [1.07, 1.13] |
| **B** | 0.138 | 0.158 | 0.178 | not run | **+0.0193 [+0.0145, +0.0242]** | 1.12 [1.09, 1.15] |
| C | 0.152 | 0.165 | 0.173 | not run | +0.0083 [+0.0038, +0.0128] | 1.05 [1.02, 1.08] |
| B-g3154 | 0.147 | 0.163 | 0.178 | not run | +0.0148 [+0.0108, +0.0189] | 1.09 [1.07, 1.12] |
| B-mix | 0.124 | 0.147 | 0.160 | not run | +0.0140 [+0.0090, +0.0190] | 1.10 [1.06, 1.13] |
| R1370 | 0.144 | 0.162 | 0.177 | not run | +0.0156 [+0.0111, +0.0201] | 1.10 [1.07, 1.12] |

All-evidence hub-controlled attention - CLR: A +0.0090 [+0.0047, +0.0133]; B +0.0104 [+0.0061, +0.0147]; C +0.0030 [-0.0011, +0.0072]; B-g3154 +0.0085 [+0.0047, +0.0123]; B-mix +0.0035 [-0.0004, +0.0074]; R1370 +0.0076 [+0.0034, +0.0118].

In every subset attention > CLR > |r| on both truths. GENIE3 was run only on A, where its scores were already cached (SC 0.141, all 0.168; below attention in both). On the subsets GENIE3 was skipped: it took 358 min on the full 2,198 arrays at 24 workers, which is not cheap, and GRNBoost2 was skipped as instructed. PIDC was not run.

## Label-free depth t* on each subset

Same criterion as the frozen EcoMAC selection (`analysis/keff_test.R`, variant K10, criterion A): 20 balanced author splits (same generator and seed; authors of the subset), training-half K10 attention mass scored by AP for the held-out half's top-1% CLR pairs, grid {1,2,3,4,5,6,8,10,12,16,20,24,32}, per-direction argmax, geometric mean over the two directions, median over splits. `analysis/ecomac_nonoverlap_tstar.R`; no truth read. **Validation:** on subset A the 520 curve points equal the frozen `results/keff_test/ecomac/heldout_curves.csv` (K10) to 5e-16 and give t* = 8, the frozen value.

| Subset | t* | split IQR | split range | median 1% plateau (t) |
|---|---|---|---|---|
| A | 8 | 7-8 | 6-9 | 5-10 |
| B | 8 | 8-8 | 7-9 | 6-10 |
| C | 9 | 8-10 | 8-11 | 8-10 |
| B-g3154 | 11 | 10-11 | 9-12 | 10-12 |
| B-mix | 8 | 6-8 | 5-9 | 8-10 |

Scores at t* (SC truth; hub-controlled in the last column; same code, `--t` set to the subset's t*; C is at t* = 9 = the main table):

| Subset (t*) | CLR | Attention | Attention - CLR (SC) | Ratio (SC) | Hub-controlled difference (SC) |
|---|---|---|---|---|---|
| A (8) | 0.130 | 0.157 | +0.0263 [+0.0201, +0.0325] | 1.20 [1.16, 1.25] | +0.0188 [+0.0134, +0.0242] |
| B (8) | 0.129 | 0.156 | +0.0271 [+0.0207, +0.0336] | 1.21 [1.16, 1.26] | +0.0192 [+0.0136, +0.0247] |
| B-g3154 (11) | 0.136 | 0.160 | +0.0247 [+0.0181, +0.0313] | 1.18 [1.14, 1.23] | +0.0173 [+0.0117, +0.0230] |
| B-mix (8) | 0.114 | 0.133 | +0.0190 [+0.0138, +0.0242] | 1.17 [1.12, 1.21] | +0.0105 [+0.0061, +0.0148] |

All-evidence differences at t*: A +0.0149, B +0.0177 [+0.0131, +0.0224], B-g3154 +0.0172 [+0.0128, +0.0216], B-mix +0.0129 [+0.0081, +0.0177] (hub-controlled B +0.0102 [+0.0060, +0.0145], B-g3154 +0.0085, B-mix +0.0033 [-0.0005, +0.0072]). Depth does not change any conclusion.

## Contrasts between subsets (descriptive)

Paired on the same gene halves (genes resampled only; arrays are not resampled, and B contains the arrays of B-g3154 and B-mix and is contained in A, so only disjoint pairs are interpretable). Difference of (attention - CLR):

| Contrast | SC | SC hub-controlled | All evidence | All evidence hub-controlled |
|---|---|---|---|---|
| B - C | +0.0057 [+0.0009, +0.0106] | +0.0048 [+0.0000, +0.0096] | +0.0110 [+0.0064, +0.0156] | +0.0074 [+0.0034, +0.0114] |
| B-g3154 - B-mix | +0.0034 [-0.0017, +0.0085] | +0.0059 [+0.0011, +0.0107] | +0.0009 [-0.0036, +0.0053] | +0.0050 [+0.0008, +0.0092] |
| B - R1370 | +0.0012 [-0.0012, +0.0036] | +0.0013 [-0.0011, +0.0038] | +0.0037 [+0.0014, +0.0061] | +0.0028 [+0.0005, +0.0052] |

## Reading

1. The attention gain over CLR on EcoMAC does not depend on the arrays that M3D also contains. On the 1,370 arrays with no M3D counterpart the SC difference is +0.0285 [+0.0218, +0.0353] (ratio 1.22), the same as on all arrays (+0.0278) and on a random 1,370-array subset (+0.0273); the hub-controlled difference is +0.0197 [+0.0141, +0.0254] versus +0.0195. On the called arrays alone the gain is smaller (+0.0228 SC; +0.0083 all evidence; hub-controlled all evidence +0.0030, interval including 0), so if anything the overlap arrays are the less favourable part for attention, and the headline is not carried by them.
2. The gain is not a Genome 2.0 platform effect: on GPL3154 only it is +0.0230 SC (hub-controlled +0.0167), and on GPL199-uncalled plus no-GEO arrays it is +0.0197 SC (hub-controlled +0.0108); the all-evidence hub-controlled gain on B-mix (+0.0035) is not distinguishable from 0, as on C.
3. This is a robustness statement about array dependence, not a held-out replication: the same RegulonDB truth, the same genes, the transferred depth and the same pipeline. It does not recover independence from M3D in the labs (the uncalled arrays still come from many of the same laboratories; 29 author surnames are shared) and it cannot exclude near-duplicates among uncalled arrays (the 95 pre-build GPL199 arrays that missed the rule are neither confirmed overlap nor independent).

## Caveats

- Not preregistered; the choice of the headline call (mutual nearest neighbour) was itself post hoc in the overlap audit. Subsets were defined by M3D similarity, not by outcome; no subset was dropped or added after seeing results, except that R1370 and the t* runs for all subsets were added at the outset of the analysis.
- Consumed truth; the AUPR levels differ across subsets in part because the subsets differ in platform and in genes' dynamic ranges (B-mix is lower for all methods).
- Intervals are gene jackknife only (network fixed), as in Table 2; no experiment-level or array-level resampling, and no multiplicity correction across the many cells above (the cells are descriptive). Contrast intervals likewise.
- ComBat re-fitted per subset (two platforms for B-mix, one for B-g3154: no correction applied).
- Compute: nice 19 scopes via the scheduler; baselines jobs peak 3.4-4.5 GB (declared 5), t* jobs 2.4-3.3 GB (declared 6).

## Proposed supplement note (not added to `supplement.md`)

> **Supplementary Note X | Re-scoring EcoMAC without M3D-matched arrays (post hoc sensitivity analysis).** An overlap audit found that 828 of the 2,198 EcoMAC arrays (37.7%) have a one-to-one counterpart in M3D (centred-correlation rule and mutual nearest neighbour), covering 91% of M3D chips, so EcoMAC is not an independent replication. As a sensitivity analysis, not preregistered and reusing the RegulonDB truth of Table 2, we re-ran the identical pipeline (ComBat by platform re-fitted within the subset, pooled AUPR, K10 attention at the transferred t = 9, gene-jackknife intervals) on the 1,370 EcoMAC arrays without an M3D counterpart, chosen by similarity to M3D and not by outcome. Attention remained ahead of CLR and |r|: AUPR 0.157 versus 0.129 and 0.106 on strong/confirmed regulons, a difference of +0.0285 (95% interval +0.0218 to +0.0353; ratio 1.22), against +0.0278 on all arrays and +0.0273 on a random subset of the same size; after hub stratification the difference was +0.0197 (+0.0141 to +0.0254). The gain held on Genome 2.0 arrays alone (+0.0230; hub-controlled +0.0167) and on the remaining GPL199 and unannotated arrays (+0.0197; +0.0108), and at the label-free depth of the subset (t* = 8; +0.0271). On the 828 matched arrays alone it was smaller (+0.0228; at all evidence levels +0.0083, hub-controlled +0.0030, interval including 0). Intervals resample genes only. This shows that the EcoMAC result is not carried by the arrays shared with M3D; it is not an independent test, because the truth, genes and pipeline are the same and many laboratories contribute to both compendia.

## Files

`analysis/baselines.R` (adds `--array_subset`), `analysis/ecomac_nonoverlap_subsets.R`, `analysis/ecomac_nonoverlap_tstar.R`, `analysis/ecomac_nonoverlap_summary.R`; results in `results/ecomac_nonoverlap/` (`subsets/`, `<subset>/metrics.csv`, `<subset>/jack_replicates.csv`, `<subset>_tstar/`, `tstar/`, `summary_metrics.csv`, `summary_contrasts.csv`, `subset_composition.csv`). Jobs 0960-0964 (`jobs_pc/done/`).
