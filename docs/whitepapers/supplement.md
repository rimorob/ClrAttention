# Supplementary Information

**Context-normalized attention without learning recovers gene regulatory programs**

*Draft, 2026-10-14.  Numbers await an independent QC pass; [pending] marks analyses now running.  Figures are in `figures/`.*

---

## Supplementary Note 1 | The dependence measure under batch effects

**Motivation.**  On the single-platform M3D compendium, CLR on mutual information (MI) and CLR on correlation perform alike (main text, Table 1).  Real compendia are rarely that homogeneous.  We therefore repeated the comparison on EcoMAC (Carrera et al. 2014), which pools 2,198 RMA-normalized *E. coli* arrays (4,189 genes) from GEO, ASAP and ArrayExpress, without cross-platform batch correction:

- 844 arrays on the Affymetrix E. coli Antisense array (GPL199);
- 735 arrays on the E. coli Genome 2.0 array (GPL3154);
- 619 arrays without a GEO record.

The platform of each array was taken from its GEO series or, for series that mix platforms, from its GEO sample record (`tools/get_ecomac.sh`).

**Batch structure bends pairwise relationships (Supplementary Fig. 1).**  In M3D, lexA and its SOS targets are related monotonically.  In EcoMAC, each platform group occupies its own region of the plane.  For lexA–sulA the result is a "Christmas tree":

- the platform groups spread horizontally at different heights;
- the SOS-induced arrays form a vertical trunk near high lexA.

Pearson r falls from 0.76 (M3D) to 0.52 (EcoMAC).  Within the Genome 2.0 arrays alone it is 0.35.  Pooled correlation in such data mixes a within-platform relationship with between-platform offsets.

**MI registers this structure; correlation largely does not (Supplementary Fig. 2).**  We compared B-spline MI with |Pearson r| for all gene pairs, against the MI of a bivariate Gaussian with the same correlation, −½ log₂(1 − r²) bits.

| Data | Pairs whose MI exceeds the Gaussian value by > 0.1 bits |
|---|---|
| M3D (one platform) | 0.2% |
| EcoMAC, uncorrected | **8.4%** |
| EcoMAC, per-platform centred | 0.3% |

In uncorrected EcoMAC, a band of pairs with moderate correlation (|r| 0.2–0.6) carries 0.3–0.6 bits.  These are the platform mixtures.  Per-platform centring removes nearly all of them.

The figure should not be read as MI "outperforming" Gaussian MI at low correlation, or "trailing" it at high correlation.  Those departures are biases of the estimator:

- **Near r = 0:** a positive finite-sample floor of about 0.05–0.15 bits, which is higher with fewer samples.
- **At high |r|:** saturation, because 7–9 joint bins cannot represent very strong dependence.

The informative feature is the excess band, which is absent in single-platform data.

**Regulon recovery (Supplementary Table 1).**  This is the main-text benchmark:  RegulonDB co-membership AUPR and median regulon coherence at strong/confirmed evidence.  Attention is soft10 at t = 8, the M3D held-out depth, not re-selected here.

**Supplementary Table 1 | EcoMAC: dependence measure, batch correction and platform.**  Each cell gives AUPR / coherence.

| Arrays | n | \|r\| | MI | CLR(\|r\|) | CLR(MI) | CLR-NS | Att. CLR(\|r\|) | Att. CLR(MI) | Att. CLR-NS |
|---|---|---|---|---|---|---|---|---|---|
| All arrays, uncorrected | 2198 | 0.097 / 0.526 | 0.099 / 0.534 | 0.119 / 0.537 | 0.120 / 0.556 | 0.120 / 0.535 | 0.148 / 0.595 | 0.152 / 0.604 | 0.151 / 0.587 |
| All arrays, per-platform centred | 2198 | 0.111 / 0.556 | 0.114 / 0.579 | 0.129 / 0.568 | 0.130 / 0.580 | 0.126 / 0.559 | 0.156 / 0.597 | 0.157 / 0.617 | 0.153 / 0.614 |
| Antisense (GPL199) only | 844 | 0.115 / 0.572 | 0.118 / 0.581 | 0.127 / 0.584 | 0.129 / 0.586 | 0.126 / 0.578 | 0.147 / 0.587 | 0.151 / 0.586 | 0.146 / 0.580 |
| Genome 2.0 (GPL3154) only | 735 | 0.116 / 0.569 | 0.118 / 0.593 | 0.133 / 0.568 | **0.136 / 0.599** | 0.132 / 0.573 | 0.154 / 0.603 | **0.158 / 0.618** | 0.156 / 0.604 |
| No GEO record only | 619 | 0.089 / 0.495 | 0.091 / 0.504 | 0.110 / 0.522 | 0.110 / 0.522 | 0.109 / 0.527 | 0.134 / 0.587 | 0.132 / 0.566 | 0.131 / 0.577 |

*CLR-NS: CLR on |r| of normal scores (Blom transform, average ranks for ties).*

**Reading.**

1. **Without batch correction, MI is the better dependence measure for CLR.**  CLR(MI) and CLR(|r|) tie on AUPR (0.120 against 0.119), but MI gives higher coherence (0.556 against 0.537).  The advantage persists with attention (0.604 against 0.595).
2. **Simple batch correction helps more than the choice of measure.**  Centring each gene within each platform group raises CLR(MI) from 0.120 to 0.130 AUPR and from 0.556 to 0.580 coherence.  Uncorrected MI is therefore below corrected correlation (0.129 / 0.568).  MI keeps a small edge after correction, mainly in coherence (0.580 against 0.568; with attention 0.617 against 0.597).
3. **More arrays are not automatically better.**
   - The 735 Genome 2.0 arrays alone give the best CLR(MI) of any configuration (0.136 / 0.599), above all 2,198 arrays even after centring.
   - Within that platform, MI clearly beats correlation on coherence (0.599 against 0.568).
   - The arrays without a GEO record are the weakest group.  Pooling heterogeneous sources dilutes the signal unless the heterogeneity is modelled.
4. **Attention helps in every configuration.**  The gain is largest on the uncorrected pool:  +0.032 AUPR and +0.048 coherence over CLR(MI), against +0.021 / +0.031 on M3D.

These are single fits without intervals.  [pending: held-out depth and experiment-level bootstrap intervals on EcoMAC, and ComBat-style correction as a stronger batch model.]

---

## Supplementary Note 2 | Zero inflation in single-cell data:  normal scores

In single-cell RNA-seq, many entries are exact zeros:  11–38% in BEELINE's selected gene sets, and more in wider gene contexts.  The B-spline MI estimator spends its resolution on that point mass.

We compared dependence measures for CLR on all seven BEELINE datasets and both gene sets, 72 dataset × gene set × truth combinations (`analysis/beeline_blom_test.R`).  Attention depths were taken from the raw-MI analysis and were not re-selected per variant.

**Supplementary Table 2 | CLR on BEELINE:  change in AUPRC ratio relative to CLR on raw-data MI.**

| CLR input | Median Δ AUPRC ratio | Combinations improved | Paired Wilcoxon P |
|---|---|---|---|
| \|r\| of normal scores, average ranks for ties (CLR-NS) | **+0.115** | 57/72 | 1 × 10⁻⁶ |
| \|r\| of normal scores, ties broken at random, mean of 10 draws | +0.071 | 50/72 | 0.002 |
| \|r\| of normal scores, ties broken at random, one draw | +0.053 | 49/72 | 0.005 |
| B-spline MI of normal scores, random ties | +0.021 | 46/72 | 0.42 |
| Gaussian MI, −½ log₂(1 − r²), of normal scores | +0.013 | 37/72 | 0.40 |
| B-spline MI of normal scores, average ties | −0.008 | 33/72 | 0.07 |

**Notes.**

- **Ties.**  With random tie-breaking, results for truths with few edges varied by up to 4.7 AUPRC-ratio units between seeds, so average ranks are preferred.  With average ranks, all zeros share one score.  Correlation handles that point mass well; the B-spline estimator does not.
- **Gaussian MI as a CLR input.**  It is worse than |r| itself.  The logarithm stretches each row's upper tail, so the per-gene z-scores become dominated by a few pairs.
- **Attention on CLR-NS.**  It improves on CLR for truths that reward co-regulated target sets (transcriptional-regulation STRING:  2.17 against 1.74).  It loses on truths dominated by co-expression evidence.

[pending: CLR-NS with its own held-out depth, in 5,000-gene and all-gene contexts, with nested cell-bootstrap intervals.]

Supplementary Notes 1 and 2 support a single conclusion:  the right dependence measure depends on the data.

- **Bulk, single platform:**  MI and correlation are interchangeable.
- **Pooled across platforms:**  MI registers batch mixtures that correlation misses, and it retains a coherence advantage even within a platform.
- **Zero-inflated counts:**  rank-based correlation is better.

The context normalization, and the attention built on it, help in every case.

---

## Supplementary Note 3 | Component ablation with held-out depth

Each operator applies the same softmax (10 effective neighbours, 50 candidates) to a different score matrix.  Its depth is chosen by its own held-out criterion over 20 random laboratory splits:  how well attention on one half of the laboratories reproduces the top 1% of the *same* score matrix computed on the other half (`analysis/lab_splits.R --ops ...`).  Intervals are paired delete-half gene jackknife (AUPR, ±1.96 SE) and regulon bootstrap (coherence) relative to CLR on MI without attention.

**Supplementary Table 5 | Component ablation at held-out depths** (466 experiments; SC = strong/confirmed evidence).

| Scores fed to the softmax | Context-normalized | t* | AUPR, SC | Δ AUPR, SC | Coherence, SC | Δ coherence, SC | Δ AUPR, all evidence |
|---|---|---|---|---|---|---|---|
| CLR scores ranked directly | yes | — | 0.136 | — | 0.597 | — | — |
| Raw \|r\|, no CLR | no | 8 | 0.143 | +0.008 [+0.003, +0.012] | 0.586 | +0.005 [−0.004, +0.017] | −0.005 [−0.010, −0.001] |
| Raw mutual information, no CLR | no | 6 | 0.143 | +0.007 [+0.001, +0.012] | 0.585 | +0.007 [−0.008, +0.018] | −0.002 [−0.007, +0.002] |
| CLR on \|r\| | yes | 7 | 0.153 | +0.017 [+0.011, +0.024] | 0.620 | +0.030 [+0.015, +0.044] | +0.002 [−0.003, +0.006] |
| CLR on \|r\| of normal scores | yes | 5 | 0.151 | +0.015 [+0.009, +0.021] | 0.632 | +0.032 [+0.018, +0.043] | +0.002 [−0.002, +0.006] |
| CLR on mutual information | yes | 9 | 0.158 | +0.022 [+0.016, +0.028] | 0.628 | +0.031 [+0.019, +0.044] | +0.009 [+0.005, +0.012] |

**Notes.**

- **The criteria differ by operator.**  Each t* is self-consistent, so the depths are internally valid but criterion values cannot be compared across operators.  These intervals do not include the uncertainty of choosing t*.
- **At all evidence levels,** only attention on CLR(MI) keeps a significant AUPR gain at its own t*.  At t = 8, CLR on normal-score |r| is also significant (+0.005; +0.001 to +0.010).  Raw |r| without CLR is significantly *below* CLR.
- **Coherence at all evidence levels** stays significant for all three context-normalized kernels.

## Supplementary Note 4 | Depth selection

Experiments of the 466-experiment M3D set were split 50 times at random into two halves by laboratory (39 labs; labs never straddle halves; half sizes 203–264 experiments).  In each split and direction, the held-out criterion was evaluated at every depth t = 1–64; t* per split is the geometric mean of the two directions' optima, and the reported depth is the median over splits (`analysis/lab_splits.R`; recomputed independently from the saved curves).

**Supplementary Table 3 | Held-out depth and softmax candidate cap.**  SC = strong/confirmed evidence; intervals are paired delete-half gene jackknife (AUPR) and regulon bootstrap (coherence) relative to CLR.

| Operator | Median t* (IQR; range) | Criterion within 1% of max | Δ AUPR vs CLR, SC | Δ coherence, SC | Δ AUPR, all evidence |
|---|---|---|---|---|---|
| soft10, cap 50 (default) | 9 (9–10; 8–11) | t = 7–11 | +0.022 [+0.017, +0.027] | +0.031 [+0.019, +0.043] | +0.009 [+0.005, +0.012] |
| soft10, cap 100 | 9 (9–10; 6–11) | t = 7–12 | +0.022 [+0.017, +0.026] | +0.033 [+0.023, +0.042] | +0.009 [+0.005, +0.013] |
| soft10, cap 200 | 9 (8–10; 5–11) | t = 6–12 | +0.020 [+0.016, +0.025] | +0.030 [+0.024, +0.044] | +0.009 [+0.005, +0.013] |
| soft10, no cap | 5 (4–6.75; 3–9) | t = 2.5–10 | +0.012 [+0.008, +0.016] | +0.030 [+0.021, +0.035] | +0.005 [+0.001, +0.008] |
| pearson_top10 (control) | 6 (5–6; 5–6) | t = 5–6 | — | — | — |

**Notes.**

- **Direction-level spread.**  The split-level range for soft10 (8–11) is over geometric means; individual directions range from 6 to 14, and in 17 of 50 splits the two directions differ by more than 2.
- **The uncapped softmax** has a bimodal depth (3–5 in 35 splits, 7–9 in 13).  With every positive score in the softmax, weight leaks to many weak neighbours and multi-step paths blur sooner.
- **The held-out depth is conservative.**  Pooled RegulonDB AUPR on the full data keeps rising past t* and peaks near t = 27 (SC), while coherence peaks at t = 4–8.  The criterion never sees RegulonDB, and it should not be read as benchmark-optimal.
- **The pearson_top10 control** is tuned against the same CLR-on-MI target, a criterion it does not share; its depth is therefore not optimized on its own terms.  The kernel ablation (Supplementary Note 3) uses each operator's own target.

## Supplementary Note 5 | Nested bootstrap

Each of 1,000 cluster draws (466 experiments resampled with replacement, one replicate array per occurrence; the same draws as the fixed-depth bootstrap, identical to within 10⁻¹² for shared methods) repeats the full procedure.  The draw's laboratories are split at random, t is selected on that split from t = 1–32, 40, 48, 56, 64, and attention is scored at the draw's own t* (`analysis/nested_bootstrap.R`).  No selected depth reached 32.

**Supplementary Table 6 | Nested bootstrap** (differences from CLR within draw; mean and 2.5–97.5% percentiles over 1,000 draws).

| Evidence | Method | Depth | Δ co-membership AUPR | Δ coherence | Δ TF-node AUPR (× random) | Δ median per-TF AUROC |
|---|---|---|---|---|---|---|
| SC | Attention, own t* | median 10 (8–12) | +0.022 [+0.019, +0.026] | +0.026 [+0.008, +0.045] | −1.71 [−2.16, −1.27] | +0.027 [−0.010, +0.061] |
| SC | Attention, t = 8 | 8 | +0.020 [+0.018, +0.023] | +0.028 [+0.010, +0.046] | −1.67 [−2.11, −1.25] | +0.031 [−0.004, +0.065] |
| SC | pearson_top10, own t* | 6 | −0.002 [−0.006, +0.001] | −0.026 [−0.046, −0.005] | −2.48 [−3.00, −1.97] | −0.039 [−0.079, +0.008] |
| All | Attention, own t* | median 10 | +0.009 [+0.006, +0.012] | +0.011 [−0.004, +0.025] | −1.10 [−1.39, −0.84] | +0.018 [−0.009, +0.046] |

Attention at its own depth exceeded CLR on SC co-membership AUPR in all 1,000 draws.  The per-TF AUROC column is a difference of medians across transcription factors, not a median of per-factor differences.

## Supplementary Note 6 | Perturbation targets

Influence scores are calibrated against null removals of the same number of experiments.  Three nulls answer different questions:

- **Random (500 draws per k):**  removals uniform over all experiments.  This can attribute laboratory effects to the perturbation.
- **Laboratory-matched (500 per perturbation):**  for each laboratory holding n of the perturbation's experiments, n of that laboratory's other experiments are removed.  Laboratories with fewer than 10 possible subsets fall back to all experiments.  Matching was complete for 10 of 15 perturbations and partial for RyhB (0%), HU (7%), Fis (20%), FNR (50%) and σ38 (56%).
- **Laboratory-matched excluding perturbation experiments** [pending, job 089]:  as above, with every experiment of any perturbation group excluded from the pools.

**Supplementary Table 7 | Variance-stratified AUROC of attention-based influence against differential expression** (15 perturbations).

| Null | Median AUROC | Beats DE | Paired Wilcoxon P |
|---|---|---|---|
| Random | 0.572 | 12/15 | 0.015 |
| Laboratory-matched | 0.550 | 9/15 | 0.095 |
| Laboratory-matched, excluding perturbation experiments | [pending] | | |

Differential expression scores 0.516 under every null.  The laboratory-matched null is conservative.  In the Covert laboratory, 10 of the 12 other experiments are knockouts of related regulators (ArcA, FNR, AppY) with overlapping regulons, so the null removes related biology as well as the laboratory effect.  AppY z-scores under the two nulls correlate at only r = 0.48, and the null rests on 19–74 distinct subsets for several groups.  Among the 10 fully matched groups the median falls from 0.571 to 0.544.

[pending: SSEM-Lasso (Cosgrove et al. 2008) re-implementation on the same perturbations, job 095.]

## Supplementary Note 7 | Removing global expression programs

Removing the top k principal components before co-expression network inference is a standard, label-free way to remove confounders and global programs such as growth rate (Parsana et al. 2019).  Because attention propagates along multi-step paths, a reviewer may ask whether its gain merely duplicates such a correction.  We removed k = 0–20 components, and the number chosen by the permutation test of `sva::num.sv` (Buja–Eyuboglu), from M3D and from EcoMAC after ComBat correction for platform.  We then scored CLR and attention (t = 9) at each k (`analysis/pc_removal.R`).  Attention's gain over CLR has a paired delete-half gene jackknife interval.

**Supplementary Table 4 | Attention after removal of k principal components** (strong/confirmed evidence; AUPR / median coherence).

| Data | k | CLR(MI) AUPR / coh. | Att. CLR(MI) AUPR / coh. | Δ AUPR [95% CI] | CLR(|r|) AUPR / coh. | Att. CLR(|r|) AUPR / coh. | Δ AUPR [95% CI] |
|---|---|---|---|---|---|---|---|
| M3D | 0 | 0.136 / 0.597 | 0.158 / 0.628 | +0.022 [+0.017, +0.027] | 0.138 / 0.590 | 0.156 / 0.621 | +0.018 [+0.013, +0.024] |
| M3D | 1 | 0.142 / 0.612 | 0.163 / 0.635 | +0.021 [+0.016, +0.026] | 0.143 / 0.605 | 0.160 / 0.606 | +0.018 [+0.012, +0.023] |
| M3D | 2 | 0.145 / 0.622 | 0.163 / 0.634 | +0.018 [+0.013, +0.023] | 0.145 / 0.613 | 0.160 / 0.629 | +0.014 [+0.010, +0.019] |
| M3D | 3 | 0.146 / 0.624 | 0.162 / 0.625 | +0.017 [+0.012, +0.021] | 0.146 / 0.613 | 0.160 / 0.611 | +0.014 [+0.009, +0.019] |
| M3D | 5 | 0.137 / 0.618 | 0.154 / 0.648 | +0.018 [+0.014, +0.022] | 0.135 / 0.592 | 0.151 / 0.653 | +0.016 [+0.012, +0.021] |
| M3D | 8 | 0.130 / 0.605 | 0.151 / 0.639 | +0.021 [+0.017, +0.025] | 0.127 / 0.591 | 0.143 / 0.626 | +0.016 [+0.012, +0.020] |
| M3D | 12 | 0.130 / 0.593 | 0.149 / 0.632 | +0.018 [+0.014, +0.023] | 0.128 / 0.579 | 0.145 / 0.627 | +0.017 [+0.012, +0.021] |
| M3D | 20 | 0.127 / 0.611 | 0.147 / 0.658 | +0.020 [+0.016, +0.024] | 0.125 / 0.587 | 0.147 / 0.657 | +0.022 [+0.018, +0.026] |
| M3D | 31 (num.sv) | 0.124 / 0.593 | 0.145 / 0.633 | +0.021 [+0.018, +0.025] | 0.121 / 0.577 | 0.144 / 0.637 | +0.023 [+0.019, +0.027] |
| EcoMAC (ComBat by platform) | 0 | 0.130 / 0.580 | 0.158 / 0.617 | +0.028 [+0.022, +0.034] | 0.129 / 0.567 | 0.157 / 0.598 | +0.028 [+0.022, +0.034] |
| EcoMAC (ComBat by platform) | 1 | 0.143 / 0.637 | 0.166 / 0.652 | +0.023 [+0.017, +0.029] | 0.140 / 0.605 | 0.163 / 0.625 | +0.023 [+0.017, +0.028] |
| EcoMAC (ComBat by platform) | 2 | 0.146 / 0.619 | 0.169 / 0.652 | +0.022 [+0.018, +0.027] | 0.143 / 0.604 | 0.164 / 0.620 | +0.022 [+0.016, +0.027] |
| EcoMAC (ComBat by platform) | 3 | 0.150 / 0.645 | 0.172 / 0.639 | +0.022 [+0.017, +0.026] | 0.146 / 0.607 | 0.167 / 0.621 | +0.021 [+0.016, +0.026] |
| EcoMAC (ComBat by platform) | 5 | 0.142 / 0.643 | 0.162 / 0.637 | +0.020 [+0.016, +0.025] | 0.138 / 0.611 | 0.156 / 0.609 | +0.018 [+0.013, +0.023] |
| EcoMAC (ComBat by platform) | 8 | 0.130 / 0.613 | 0.152 / 0.634 | +0.022 [+0.017, +0.027] | 0.127 / 0.598 | 0.148 / 0.616 | +0.021 [+0.017, +0.026] |
| EcoMAC (ComBat by platform) | 12 | 0.129 / 0.618 | 0.153 / 0.641 | +0.024 [+0.019, +0.029] | 0.126 / 0.591 | 0.147 / 0.638 | +0.021 [+0.017, +0.026] |
| EcoMAC (ComBat by platform) | 20 | 0.134 / 0.632 | 0.163 / 0.686 | +0.029 [+0.024, +0.035] | 0.128 / 0.590 | 0.154 / 0.643 | +0.026 [+0.021, +0.030] |
| EcoMAC (ComBat by platform) | 68 (num.sv) | 0.124 / 0.639 | 0.153 / 0.709 | +0.029 [+0.024, +0.034] | 0.117 / 0.604 | 0.143 / 0.645 | +0.026 [+0.021, +0.031] |

**Reading.**

1. **The attention gain survives component removal.**  At every k, on both compendia, attention on CLR(MI) beats CLR(MI), with intervals excluding zero (M3D +0.017 to +0.022; EcoMAC +0.020 to +0.029).  Attention does not merely duplicate the removal of global programs.
2. **Removing a few components helps both methods.**  With k = 2–3, CLR improves (M3D 0.136 → 0.146; EcoMAC 0.130 → 0.150) and attention on top of it gives the highest AUPR of any configuration.  Beyond about five components AUPR declines, as regulatory programs themselves start to be removed.
3. **The permutation rule over-removes for this purpose.**  It selects k = 31 (M3D) and 68 (EcoMAC), where CLR falls below the uncorrected AUPR, although coherence with attention keeps rising.  We therefore keep k = 0 in the main analyses.  Choosing k = 2–3 because it scores best against RegulonDB would tune the preprocessing on the gold standard.
4. **Mutual information is at least as good as |r| at every k,** and its advantage on EcoMAC widens as components are removed.

A label-free criterion based on held-out reproducibility could not choose k either:  it favours keeping reproducible global programs, which are real but not regulon-specific (`analysis/ecomac_batch.R`).

## Supplementary Note 8 | The top principal components are sigma-factor and growth programs

Removing the first two or three principal components improves both CLR and attention on M3D (Supplementary Note 7), so it matters what those components contain.  For each of the first five components of the 466-experiment compendium, we asked how well the absolute gene loading separates the members of a regulon from all other genes (AUROC; two-sided Mann–Whitney test, Benjamini–Hochberg adjusted; `analysis/pc_regulons.R`).  An AUROC of 0.5 means no association.

**Supplementary Table 8 | Regulons on the top principal components** (AUROC of |loading|; * BH-adjusted P < 0.05; all evidence levels; percentage of variance explained in the header).

| Gene set | n | PC1 (24.4%) | PC2 (9.6%) | PC3 (7.9%) | PC4 (6.6%) | PC5 (4.3%) |
|---|---|---|---|---|---|---|
| Ribosomal proteins (growth rate) | 54 | 0.98* | 0.22* | 0.64* | 0.46 | 0.42 |
| σ28 (flagella, motility) | 73 | 0.55 | 0.80* | 0.82* | 0.66* | 0.45 |
| FlhDC (flagellar master regulator) | 92 | 0.47 | 0.72* | 0.67* | 0.69* | 0.39* |
| σ38 / σS (stationary phase, stress) | 331 | 0.64* | 0.56* | 0.63* | 0.58* | 0.59* |
| Fis (growth phase) | 188 | 0.64* | 0.55* | 0.52 | 0.69* | 0.62* |
| ArcA (respiration) | 232 | 0.60* | 0.59* | 0.54 | 0.64* | 0.67* |
| σ70 (housekeeping) | 1596 | 0.56* | 0.55* | 0.54* | 0.62* | 0.60* |
| σ32 (heat shock) | 146 | 0.53 | 0.45 | 0.50 | 0.50 | 0.46 |
| σ24 (envelope stress) | 128 | 0.54 | 0.44* | 0.56* | 0.46 | 0.50 |
| σ54 (nitrogen) | 136 | 0.48 | 0.44* | 0.37* | 0.48 | 0.53 |
| CRP (carbon catabolite) | 551 | 0.51 | 0.50 | 0.50 | 0.63* | 0.54* |

**PC1 (24% of variance) is growth rate.**  Ribosomal-protein genes separate almost perfectly (AUROC 0.98).  The experiments at the two extremes are 15–24 h biofilms on one side and early time points of rapidly growing cultures on the other.  The σS regulon (0.64) and Fis (0.64), whose levels track growth phase, lie along the same axis.

**PC2 and PC3 are motility and stationary phase.**  The flagellar σ28 regulon and its master regulator FlhDC load strongly on both (σ28 0.80 and 0.82; FlhDC 0.72 and 0.67).  PC3 is anchored by dilution time courses into stationary phase and by σS (rpoS) experiments, and σS loads on it as well (0.63).  The end of PC2 opposite the motility genes is formed by experiments from a single overexpression series grown in 96-well plates, so this component mixes biology with a laboratory-specific culture condition.  The heat-shock (σ32), envelope-stress (σ24) and nitrogen (σ54) sigmulons show no enrichment, and neither do most specific transcription-factor regulons.

**Implication.**  Bacterial expression variance is dominated by a few condition-specific programs, chiefly growth rate and the sigma-factor responses that accompany it (σS for stationary phase and stress, σ28 for motility).  Removing the top components therefore deletes genuine regulatory biology, not technical noise.  It helps the co-membership benchmark because these large, diffuse programs create many co-expressed pairs outside any specific regulon, and those pairs dilute the precision of regulon recovery.  Their removal is therefore a choice about which biology to emphasize, not a correction.  We keep all components in the main analyses.  Attention's gain over CLR persists after removal (Supplementary Note 7), so it is not explained by these global programs.  Whether the same holds in organisms without sigma-factor-dominated transcription is a question for human data [pending].

---

## Supplementary Figures

**Supplementary Fig. 1 | Batch structure in pooled compendia.**  Expression of lexA against two SOS genes in:

1. M3D build 6 (907 arrays, one platform);
2. EcoMAC, all 2,198 arrays, coloured by platform;
3. EcoMAC Genome 2.0 arrays only.

Each panel reports n, Pearson r and B-spline MI.  (a) lexA–sulA, `figures/scatter_lexA_sulA.png`.  (b) lexA–recA, `figures/scatter_lexA_recA.png`.

**Supplementary Fig. 2 | Mutual information against correlation for all gene pairs.**  Two-dimensional histograms (log counts) of B-spline MI against |Pearson r| for M3D, EcoMAC uncorrected and EcoMAC with per-platform centring.  The curve is the MI of a bivariate Gaussian, −½ log₂(1 − r²).  `figures/mi_vs_rho.png`.
