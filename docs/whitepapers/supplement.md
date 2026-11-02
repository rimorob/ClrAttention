# Supplementary Information

**Context-normalized attention without learning recovers gene regulatory programs**

*Draft, updated 2026-11-02.  Results have been checked by an independent statistical review unless marked otherwise (record in `docs/paper_plan.md`); [PENDING] marks analyses that are preregistered or running.  Figures are in `figures/`.  Notes 1–8 keep their earlier numbering; Notes 9–18 are new (Note 16, the tie-convention sensitivity analysis, added 2026-10-31; Note 17, the EcoMAC/M3D overlap check, and Note 18, the preregistered yeast attempt that could not be run, added 2026-11-02).*

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

These are single fits without intervals.  ComBat correction by platform, a stronger batch model than centring, is used in Supplementary Note 7 and in the method comparison (main text, Table 2); like centring, it helps both dependence measures more than switching between them, and mutual information keeps a coherence edge after correction.  [PENDING:  experiment-level bootstrap intervals on EcoMAC.  The audit of arrays shared with M3D is done (Supplementary Note 17): 828 of the 2,198 EcoMAC arrays (37.7%) have a one-to-one M3D counterpart, so EcoMAC is not an independent replication, and its +0.032 is a point estimate.]

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

**Inference.**  The P values in Supplementary Table 2 treat the 72 combinations as independent, but they come from 7 datasets.  The inferential basis is the per-dataset blocked comparison (Supplementary Note 13):  CLR-NS improves on CLR in 7 of 7 datasets (mean log2 AUPRC ratio +0.06 on BEELINE's primary truths; exact P = 0.016, the smallest attainable with 7 blocks).  CLR-NS was chosen on these truths.

**Attention and wider contexts.**  Overall, attention adds nothing over CLR-NS.  Building CLR-NS and attention on wider gene contexts (the 5,000 most variable genes, all expressed genes, or genes selected as driven by a multimodality or principal-component filter) gives no gain; one dataset drives a nominal loss.  The driven-gene filter is better than variance ranking, which drops well-expressed genes, and Pearson residuals are no better than log-normalization for single-cell co-expression.  [PENDING:  nested cell-bootstrap intervals for every single-cell statement.]

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

## Supplementary Note 6 | Perturbation targets (preliminary; the subject of separate work)

*This analysis extends the operator to a different task and is reported here for completeness.  Its results are mixed, and a full treatment, with differential CLR and single-cell Perturb-seq contrasts, is planned as separate work.*

**Task.**  A stringent test of a network is whether it says something about an intervention that the measurements of the intervention do not already say.  We posed it on M3D, where the identification of perturbation targets was first benchmarked (Cosgrove et al. 2008).  For 15 groups of experiments in which a regulator was deleted, overexpressed or mutated, the targets are the regulator's RegulonDB regulon, mapped by a table fixed in advance.  Each gene is scored by the influence of the perturbation's experiments on its row of the network, that is, the root-mean-square change in the row when those k experiments are removed, z-scored against null removals of k other experiments.

**A variance confound.**  Gene variance, which ignores the perturbation entirely, reaches a median AUROC of 0.646 and beats differential expression (DE) in 14 of 15 perturbations, because regulon members are more variable than other genes.  We therefore compare targets only with non-targets of similar variance:  the variance-stratified AUROC pools Mann–Whitney statistics within deciles of gene variance, computed without the perturbation's own experiments.  Under this control gene variance falls to 0.521 and DE scores 0.516, whereas attention-based influence keeps a median AUROC of 0.572 and beats DE in 12 of 15 perturbations (paired Wilcoxon P = 0.015); influence on CLR and on correlation reaches 0.543 (P = 0.15) and 0.530 (P = 0.80).  The effect is concentrated in specific regulators (HU, AppY, LexA/SOS and RyhB) and absent for the global regulators ppGpp and CRP, whose regulons are diffuse.

**Nulls.**  Influence scores are calibrated against null removals of the same number of experiments.  Three nulls answer different questions:

- **Random (500 draws per k):**  removals uniform over all experiments.  This can attribute laboratory effects to the perturbation.
- **Laboratory-matched (500 per perturbation):**  for each laboratory holding n of the perturbation's experiments, n of that laboratory's other experiments are removed.  Laboratories with fewer than 10 possible subsets fall back to all experiments.  Matching was complete for 10 of 15 perturbations and partial for RyhB (0%), HU (7%), Fis (20%), FNR (50%) and σ38 (56%).
- **Laboratory-matched excluding perturbation experiments** [PENDING, job 089]:  as above, with every experiment of any perturbation group excluded from the pools.

**Supplementary Table 7 | Variance-stratified AUROC of attention-based influence against differential expression** (15 perturbations).

| Null | Median AUROC | Beats DE | Paired Wilcoxon P |
|---|---|---|---|
| Random | 0.572 | 12/15 | 0.015 |
| Laboratory-matched | 0.550 | 9/15 | 0.095 |
| Laboratory-matched, excluding perturbation experiments | [PENDING] | | |

Differential expression scores 0.516 under every null.  The laboratory-matched null is conservative.  In the Covert laboratory, 10 of the 12 other experiments are knockouts of related regulators (ArcA, FNR, AppY) with overlapping regulons, so the null removes related biology as well as the laboratory effect.  AppY z-scores under the two nulls correlate at only r = 0.48, and the null rests on 19–74 distinct subsets for several groups.  Among the 10 fully matched groups the median falls from 0.571 to 0.544.  We therefore regard the advantage over DE as suggestive rather than established.

**SSEM-Lasso.**  We re-implemented the SSEM-Lasso network filter of Cosgrove et al. (2008) on the same perturbations and truths, with each perturbation's experiments held out of training, so there is no leakage.  Variance-stratified median AUROC is 0.572 for attention-based influence (beats DE in 12/15, P = 0.015), 0.549 for SSEM (12/15, P = 0.041) and 0.516 for DE.  Head to head, attention and SSEM are not distinguishable (SSEM better in 6/15, P = 0.45), and they carry different signal (per-gene Spearman correlation 0.18).  SSEM ranks the perturbed regulator itself first in 11 of 15 perturbations (DE 2/15, attention 3/15), whereas attention leads on small regulons (HU, RyhB, AppY).  The two are complementary:  SSEM finds the regulator and attention-based influence finds its regulon.  We make no claim that attention beats SSEM.  This is not a fair regulon test for SSEM, which was designed for direct targets and keeps the regulator as a predictor.  Caveats:  Holm correction over the six methods gives P ≈ 0.09; the effective number of perturbations is below 15 (two Covert-laboratory perturbations cluster, and ArcA and FNR overlap); without HU and AppY, P = 0.057; and laboratory matching covers only 5 of 15 perturbations.  [PENDING:  SSEM with the regulator dropped from the predictors and trained on the full compendium as published; the same laboratory hold-out as the attention null; a preregistered combination of the two; more perturbations.]

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

1. **The attention gain survives component removal.**  At every k, on both (partly overlapping; Note 17) compendia, attention on CLR(MI) beats CLR(MI), with intervals excluding zero (M3D +0.017 to +0.022; EcoMAC +0.020 to +0.029).  Attention does not merely duplicate the removal of global programs.
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

**Implication.**  Bacterial expression variance is dominated by a few condition-specific programs, chiefly growth rate and the sigma-factor responses that accompany it (σS for stationary phase and stress, σ28 for motility).  Removing the top components therefore deletes genuine regulatory biology, not technical noise.  It helps the co-membership benchmark because these large, diffuse programs create many co-expressed pairs outside any specific regulon, and those pairs dilute the precision of regulon recovery.  Their removal is therefore a choice about which biology to emphasize, not a correction.  We keep all components in the main analyses.  Attention's gain over CLR persists after removal (Supplementary Note 7), so it is not explained by these global programs.  Whether the same holds in organisms without sigma-factor-dominated transcription is a question for human data [PENDING].

## Supplementary Note 9 | Hub control

**Why it is needed.**  The attention walk converges toward its stationary distribution π, so at large depth attention mass between two genes approaches a function of π_i and π_j.  On M3D, π is very skewed (maximum 310 times the median), and hubness carries real signal:  ranking pairs by H_ij = π_i + π_j alone gives an AUROC of 0.546 on strong/confirmed (SC) regulons.  The pooled AUPR rise toward t ≈ 27 could therefore partly be convergence to hubs.

**Division over-corrects.**  Dividing each method's scores by H removes the hub signal from every method, not only from attention:  it strips CLR to 1.05 times the base rate.  An earlier finding that large-t attention is "mostly hubs" used this division and must be re-tested with stratification before it is used to constrain t [PENDING].

**Within-stratum scoring.**  Pairs are assigned to ten quantile strata of H, and each method's scores are converted to percentile ranks within each stratum (`analysis/baselines.R`, `<method>_hubstrat`), so that a method gains only by ranking pairs better than other pairs of similar hubness.  Every method in the comparison is scored this way (main text, Table 2b).

**Results (M3D, t = 9).**

- *Against CLR, SC regulons:*  attention stays ahead under every adjustment (division by H +0.026; π-residual +0.020; within-stratum percentile +0.016, 95% interval 0.012 to 0.020), and ahead in all ten hub deciles in every jackknife replicate.  In the later comparison run the within-stratum lead was 0.015.
- *All-evidence regulons:*  the hub-controlled lead over CLR is small (+0.004 to +0.008).  Hub control shrinks the raw lead over CLR on M3D by one third (SC) to three quarters (all evidence).
- *Against competitors, hub-controlled pooled AUPR:*  attention leads CLR, GENIE3 and PIDC on M3D SC (by 0.015, 0.009 and 0.015) and CLR and GENIE3 on EcoMAC SC (0.020 and 0.009), and it is first in the EcoMAC SC row of Table 2b.  At all evidence levels it leads only CLR, and only on EcoMAC (0.009); it ties GENIE3 on both compendia and PIDC on EcoMAC.
- *Per regulon, AP against CLR (SC):*  EcoMAC is consistent (attention better in 66% of regulons; P = 0.002), M3D is not significant (54%; P = 0.11).  We do not claim that attention improves most regulons on M3D.
- *K562:*  the DoRothEA A–B lead over CLR survives (+0.006; 0.002 to 0.010); CollecTRI has no signal under hub control (the hub score alone is at least as good as attention); the CORUM and knockdown-response losses persist.

**Limitation.**  The stratifier is derived from the attention walk.  A method-neutral stratifier (for example CLR degree) and the per-regulon (macro) view under hub control are planned sensitivity checks [PENDING].

## Supplementary Note 10 | Per-regulon (macro) average precision

**Definition (fixed with the statistics reviewer before scoring; `analysis/macro_aupr.R`).**  For every regulon with at least 10 positive non-operon pairs, AP is computed over the same candidates as coherence, with ties broken in expectation.  A gene in several regulons counts in each; regulons are not merged, because multi-membership is meaningful.  The paired difference ΔAP = AP(comparator) − AP(attention) is computed within one run, so negative values favour attention.  Uncertainty has two sources and both must exclude 0:  regulons (BCa bootstrap for the mean, Wilcoxon signed rank for the median) and genes (delete-20% gene jackknife, 200 replicates, regulon set fixed to regulons keeping at least three positive pairs in every replicate; delete-50% as a stability check).  Holm correction is applied to the larger P within families fixed in advance (*E. coli* 40 tests; K562 primary 50; K562 secondary 40; in the current version all 40 secondary tests are present, because the secondary run was repeated with GENIE3 rows and has four comparators, |r|, CLR, GENIE3 and GRNBoost2, on five truths).  The GRNBoost2 inputs were replaced after the first run (Supplementary Note 13), so the adjusted P values of the other comparators moved through the step-down (for example |r| on M3D SC, mean 0.056 to 0.062, and on M3D all, mean 0.037 to 0.040, median 0.028 to 0.030; CLR on EcoMAC SC, median 0.86 to 0.94; PIDC on EcoMAC all, mean 0.0066 to 0.0072; |r| on EcoMAC all, median 0.63 to 0.69); their point estimates and intervals are identical to the earlier version, and no decision changed.  The gene jackknife keeps only regulons scored in every replicate, which drops small regulons, so the two intervals of the decision rule cover different regulon sets.

**Supplementary Table 9 | Per-regulon ΔAP against attention, *E. coli*** (`results/macro_aupr_v3`, GRNBoost2 regenerated 29–31 Oct 2026; t = 9; negative favours attention).  Intervals are 95% over regulons and over genes; P values are Holm-adjusted within the family of 40; values of 0.004 are at the resampling floor and are bounds.

| Data | Evidence | Comparator | Regulons | Mean ΔAP | Regulon interval | Gene interval | Holm P, mean | Median ΔAP | Holm P, median | Attention wins |
|---|---|---|---|---|---|---|---|---|---|---|
| M3D | SC | \|r\| | 123 | −0.0205 | [−0.0387, −0.0109] | [−0.0213, −0.0049] | 0.062 | −0.0011 | 0.16 | 75% |
| M3D | SC | CLR | 123 | −0.0084 | [−0.0229, −0.0013] | [−0.0047, +0.0069] | 1 | −0.0001 | 1 | 54% |
| M3D | SC | GENIE3 | 123 | −0.0041 | [−0.0173, +0.0013] | [−0.0062, +0.0082] | 1 | +0.0001 | 1 | 46% |
| M3D | SC | GRNBoost2 | 123 | −0.0153 | [−0.0284, −0.0084] | [−0.0202, −0.0042] | 0.092 | −0.0002 | 1 | 56% |
| M3D | SC | PIDC | 123 | −0.0221 | [−0.0433, −0.0115] | [−0.0216, −0.0034] | 0.20 | −0.0003 | 1 | 62% |
| M3D | all | \|r\| | 157 | −0.0204 | [−0.0350, −0.0117] | [−0.0243, −0.0060] | 0.040 | −0.0017 | 0.030 | 76% |
| M3D | all | CLR | 157 | −0.0084 | [−0.0197, −0.0021] | [−0.0080, +0.0054] | 1 | −0.0001 | 1 | 53% |
| M3D | all | GENIE3 | 157 | −0.0016 | [−0.0115, +0.0029] | [−0.0066, +0.0089] | 1 | +0.0001 | 1 | 43% |
| M3D | all | GRNBoost2 | 157 | −0.0117 | [−0.0213, −0.0066] | [−0.0175, −0.0029] | 0.18 | −0.0002 | 1 | 56% |
| M3D | all | PIDC | 157 | −0.0199 | [−0.0363, −0.0098] | [−0.0180, −0.0026] | 0.23 | −0.0002 | 1 | 55% |
| EcoMAC | SC | \|r\| | 119 | −0.0281 | [−0.0435, −0.0188] | [−0.0307, −0.0121] | 0.004 | −0.0024 | 0.21 | 81% |
| EcoMAC | SC | CLR | 119 | −0.0131 | [−0.0210, −0.0071] | [−0.0119, +0.0001] | 1 | −0.0008 | 0.94 | 66% |
| EcoMAC | SC | GENIE3 | 119 | −0.0012 | [−0.0071, +0.0038] | [−0.0081, +0.0055] | 1 | +0.0003 | 1 | 43% |
| EcoMAC | SC | GRNBoost2 | 119 | −0.0181 | [−0.0290, −0.0102] | [−0.0252, −0.0097] | 0.004 | −0.0004 | 1 | 56% |
| EcoMAC | SC | PIDC | 119 | −0.0182 | [−0.0311, −0.0100] | [−0.0196, −0.0035] | 0.15 | −0.0006 | 1 | 66% |
| EcoMAC | all | \|r\| | 155 | −0.0230 | [−0.0356, −0.0153] | [−0.0338, −0.0154] | 0.004 | −0.0021 | 0.69 | 74% |
| EcoMAC | all | CLR | 155 | −0.0104 | [−0.0176, −0.0052] | [−0.0155, −0.0020] | 0.28 | −0.0005 | 1 | 61% |
| EcoMAC | all | GENIE3 | 155 | +0.0016 | [−0.0029, +0.0084] | [−0.0073, +0.0049] | 1 | +0.0002 | 1 | 44% |
| EcoMAC | all | GRNBoost2 | 155 | −0.0114 | [−0.0196, −0.0040] | [−0.0224, −0.0083] | 0.092 | −0.0001 | 1 | 54% |
| EcoMAC | all | PIDC | 155 | −0.0173 | [−0.0277, −0.0100] | [−0.0235, −0.0073] | 0.0072 | −0.0003 | 1 | 59% |

**Reading.**

1. **Parity with CLR and GENIE3.**  No difference from CLR or GENIE3 survives the decision rule in any cell.  Against CLR the regulon interval of the mean often excludes 0 but the gene interval does not.
2. **Attention ahead of |r| (3 of 4 cells), PIDC (1 of 4, EcoMAC all evidence) and GRNBoost2 (1 of 4, EcoMAC SC evidence, Holm P = 0.004) on the mean.**  In the remaining cells the mean differences from |r|, PIDC and GRNBoost2 also favour attention and both unadjusted intervals exclude 0, but the Holm-adjusted P is above 0.05 (GRNBoost2: 0.092, 0.18 and 0.092), so under the decision rule these cells show no difference.  Medians are about zero throughout (Holm P = 1 for GRNBoost2 in 4 of 4).
3. **Where the differences sit** (`macro_strata.csv`).  Against |r|, CLR and PIDC the mean gaps come mostly from regulons of 5–24 genes (mean ΔAP there −0.011 to −0.031, against +0.004 to −0.021 in regulons of 25 genes or more).  Against GRNBoost2 they do not:  the mean ΔAP is −0.009 to −0.017 in regulons of 5–24 genes and −0.013 to −0.021 in regulons of 25 or more, and +0.0011 in regulons of 5–9 genes on EcoMAC at all evidence levels.  Attention's win rate is highest in regulons of 25 genes or more for every comparator (in the 25–99 stratum for 6 of the 20 comparator–cell combinations, otherwise in the ≥ 100 stratum); in regulons of 100 genes or more it is 60–100% (GRNBoost2 100%, 88%, 90% and 88%; n = 10, 17, 10, 16), and in regulons of 5–9 genes it is 28–77% (GRNBoost2 46%, 38%, 47% and 30%).  The pooled gain therefore comes from ranking across regulons, not within them.
4. **Stability.**  Across all families, three cells are flagged by the delete-50% check, all in the K562 secondary (all-genes) family on DoRothEA A–B (CLR, GENIE3 and GRNBoost2; the GRNBoost2 cell has mean ΔAP −0.0001, Holm P = 1).  The previous version of Table 9 (`macro_aupr_v2`) flagged one cell (CLR), and about 13 were flagged in the first run.  In the current run 2 BCa intervals were undefined and are counted as not excluding 0.

**K562.**  On transcription-factor truths the per-regulon differences are tiny (|ΔAP| ≤ 0.0025) and none is established.  On CORUM, PIDC (+0.007), CLR (+0.005) and GENIE3 (+0.004) beat attention; the CLR and GENIE3 differences were fragile in the first run (the delete-50% gene interval included 0), and no CORUM comparison is flagged in the rerun.  On knockdown-response modules, CLR (+0.097), GENIE3 (+0.084), |r| (+0.091), PIDC (+0.099) and GRNBoost2 (+0.035, Holm P = 0.0045) all beat attention per module, which wins only 2% of modules against CLR and 9% against GRNBoost2 (Supplementary Note 11); in pooled AUPR GRNBoost2 is nevertheless below attention on this truth (Table 2).  For GRNBoost2 minus attention the CORUM mean difference is +0.0013 (Holm P = 1) and the transcription-factor truths give −0.0005, −0.0002 and −0.0004 (DoRothEA A–B, A–C, CollecTRI; Holm P = 1).  The secondary all-genes K562 run (8,559 genes; |r|, CLR, GENIE3 and GRNBoost2, no PIDC, 40 tests) was rerun with GENIE3 rows; there GRNBoost2 minus attention is +0.0290 on knockdown-response modules (Holm P = 0.0035; attention wins 12%), +0.0011 on CORUM and −0.0001 to −0.0005 on the transcription-factor truths (Holm P = 1).  These K562 GRNBoost2 values come from a single seed (`results/macro_aupr_v3`).

## Supplementary Note 11 | Resolution against reach:  where attention gains and loses

**Global regulons carry the pooled gain.**  CRP, σS, Fur, IHF and FNR account for 62% of positive pair incidences in the *E. coli* benchmark, and attention's pooled lead is concentrated in them.  Per regulon (Supplementary Note 10) its win rate is highest in regulons of 25 or more genes.

**Direct edges.**  On the TF-node benchmark (199 transcription factors, Y's mRNA as the proxy for its activity) the pooled AUPR ratio on M3D (SC truths) is 6.75 for CLR, 6.54 for GENIE3, 5.53 for PIDC, 5.17 for GRNBoost2 (5.14–5.39 with seeds 1–3, `docs/grnboost2_diagnostic/m3d_main_eval.csv`), 5.17 for attention and 3.74 for |r|; on EcoMAC it is 5.82 for GRNBoost2, 5.77 for GENIE3, 5.17 for CLR, 4.79 for attention, 4.66 for PIDC and 2.43 for |r|.  GENIE3 has the best median per-TF AUROC (EcoMAC 0.690, then GRNBoost2 0.679, against 0.596 for attention) and the best median per-regulon coherence.  In the nested bootstrap (Supplementary Table 6) attention trails CLR by 1.71 times random (−2.16 to −1.27), and the difference of median per-TF AUROC is not significant.

**Complexes.**  In K562, 86% of attention's CORUM loss comes from ten ribosome and pre-rRNA complexes.  For the 60S ribosomal subunit, AP is 0.19 for attention against 0.66 for CLR.  These complexes are among the most tightly co-expressed gene sets in the cell, so the loss is consistent with the walk trading their sharp pairwise dependences for aggregation over wider neighbourhoods.  It is not specific to aggregation:  GRNBoost2, a single-step method, is also poor on these complexes (60S AP 0.29; pooled CORUM 3.66 times the base rate against attention's 3.51, 2.51 against 2.08 after hub control), and no GRNBoost2–attention difference on CORUM survives the corrections used in the paper (pooled +0.0020, −0.0057 to +0.0096; hub-controlled +0.0057, 0.0007 to 0.0108 uncorrected and −0.0009 to +0.0124 Bonferroni-corrected for the five comparisons in the Table 2 row; per complex +0.0013, Holm P = 1).

**Knockdown-response modules.**  The raw tie in Table 2a is explained by hubness:  the hub score alone reaches 2.42 times the base rate, and under hub control attention falls to 1.60 against 1.88–2.00 for CLR, GENIE3, PIDC and |r| (GRNBoost2 1.55; difference from attention −0.0039, −0.0091 to +0.0013).  Per module, attention loses to all five methods.  Against CLR, GENIE3, PIDC and |r| the loss grows with module size, so it is not a small-module effect; it is driven by tight co-expression cliques.  Against GRNBoost2 the loss is about constant in size (mean ΔAP +0.033 in modules of 10–24 genes and +0.035 in modules of 25–99, `macro_strata.csv`).  Two properties of this truth favour direct similarity:  62% of the modules are near-duplicates of one another, and they are built from a held-out half of the same knockdown signal.  [PENDING:  a cluster bootstrap over groups of near-duplicate modules and a sensitivity analysis without ribosomal genes.]

**Candidate fix.**  The preregistered multi-head operator, which averages walks at keff 3, 10 and 30, could keep the sharp neighbourhoods of tight complexes while aggregating over large regulons (Supplementary Note 15) [PENDING].

## Supplementary Note 12 | Granularity, multiscale structure and the depth–size relation

These analyses tested whether depth acts as a resolution parameter.  They did not support that claim, and we report them as negatives.

**Granularity sweep** (M3D and K562; truths from fine to coarse:  operons, transcription-factor regulons, sigmulons without σ70 and σS; CORUM, DoRothEA, CollecTRI, knockdown-response modules).

- Only operons show a sharp interior optimum:  AUPR rises from 0.30 to 0.50 at t = 4–5 and falls beyond t ≈ 8.
- Coarser truths are flat (AUPR changes of at most 0.02 over t = 8–32), with optima censored at the edge of the grid.
- There is no fine-to-coarse ordering of optimal depths (*E. coli* argmax 4, 32 and 13 for operons, regulons and sigmulons; K562 0, 26, 4, 7 and 32 for CORUM, knockdown-response modules, DoRothEA A–B, DoRothEA A–C and CollecTRI; ordering test P = 0.12).
- Single sets drive the optima:  dropping σS moves the sigmulon argmax from 24 to 13.
- No depth helps CORUM.

**Multiscale pilot.**  If regulation were organized in discrete levels, the number of effective modules N(t) of the walk would show plateaus as t grows.  It does not, for any real operator:  the smallest |d log N / d log t| is 0.43–0.64.  N(t) decays as a near power law (spectral dimension about 1–1.5), far more slowly than for a permuted operator (N = 3 at t = 409, against 30 for the real operator), so the slow structure is real but continuous.  Changing keff reshapes the curve but reveals no level, and the per-gene module size k_i(t) is not a module size at small t, because π is hub-skewed.  Plateau detection is therefore dropped; any multiscale method needs a scale-free summary (slopes, spectral gaps against the null) and π-normalized scores.

**Depth against set size** (per set, M3D).

| Sets | Spearman correlation of best depth with size |
|---|---|
| RegulonDB regulons | 0.34 (t ∝ n^0.43) |
| Label-free correlation modules | 0.44 (0.73 excluding sets at the depth ceiling) |
| Random gene sets | 0.12 |
| Operons | −0.03 once ceiling ties are removed |

Size, not regulatory type, predicts the optimum.  This is a generic property of diffusion (covering time grows with module size) and the mechanism that any multiscale claim must respect, not a biological finding.  The extra preference of the largest regulons (Fis, H-NS, IHF, FNR and Fur, at the t = 64 edge) for depth appears to be convergence toward hubs; it must be re-tested with hub-stratified scores (Supplementary Note 9).  [PENDING:  the full size–depth analysis on a logarithmic depth grid to t = 128, with 20 size-matched control draws per set, de-hubbed per-set scores, per-set bootstrap intervals on the best depth, a K562 repeat and keff 3, 10 and 30.]

## Supplementary Note 13 | BEELINE:  competitor methods

GENIE3, GRNBoost2 and PIDC were run on all 14 BEELINE cases (7 datasets × 2 gene sets) with the same genes and truths as CLR, CLR-NS and |r| (`analysis/beeline_baselines.R`).  Per dataset, log2 AUPRC ratios are averaged over gene sets and truths, and methods are compared over the seven datasets with exact Wilcoxon signed-rank and sign tests (smallest attainable P = 0.0156).  CLR was fixed in advance; CLR-NS was chosen on these truths, so comparisons against it are not confirmatory.  GRNBoost2 was run with arboreto 0.1.6 and legacy early stopping (`tools/grnboost2.py --early_stop legacy`), production seed 20261021, and rerun with seeds 1–3; the rerun was gated on the other methods' rows reproducing exactly (see the order of events below).
**Supplementary Table 10 | BEELINE primary truths (44 dataset × gene set × truth cells, 7 datasets).**

| Method | Median AUPRC ratio | Mean log2 vs CLR | Better than CLR | P (Wilcoxon) | Mean log2 vs CLR-NS | Better than CLR-NS | P (Wilcoxon) |
|---|---|---|---|---|---|---|---|
| CLR | 1.48 | — | — | — | −0.06 | 0/7 | 0.016 |
| CLR-NS | 1.60 | +0.06 | 7/7 | 0.016 | — | — | — |
| GENIE3 | 1.69 | −0.04 | 3/7 | 0.69 | −0.10 | 1/7 | 0.11 |
| PIDC | 1.56 | −0.02 | 4/7 | 0.81 | −0.08 | 1/7 | 0.16 |
| GRNBoost2 | 1.47 | −0.15 | 1/7 | 0.031 | −0.21 | 0/7 | 0.016 |
| \|r\| | 1.20 | −0.23 | 0/7 | 0.016 | −0.29 | 0/7 | 0.016 |

**Reading.**

- CLR is comparable to GENIE3 and PIDC and beats |r| in 7 of 7 datasets.  GRNBoost2 scores below CLR in 6 of 7 datasets (mean log2 −0.15; hESC is the exception, +0.015; Wilcoxon P = 0.031, sign P = 0.125, uncorrected).  We report this as a descriptive count:  the sign test P is 0.125, and with seven datasets and five methods compared with CLR no Holm- or Bonferroni-corrected P could fall below 0.05.  With seeds 1, 2 and 3 the mean log2 relative to CLR is −0.157, −0.157 and −0.159, hESC is again the only dataset above CLR in every seed, and the median AUPRC ratio is 1.456, 1.485 and 1.483 (`results/beeline_baselines/summary_gb2s{1,2,3}.csv`; no gate on these runs).  Relative to CLR-NS, GRNBoost2 is better in 0 of 7 datasets for the production seed and seed 1 and in 1 of 7 for seeds 2 and 3.  GENIE3 has the highest median over the 44 cells, but its per-dataset mean is below CLR's.
- CLR-NS has the highest per-dataset mean, but it was selected on these truths.
- The human datasets split:  on hESC, GENIE3 and PIDC beat CLR-NS, driven by STRING; on hHep, CLR-NS beats all three.  We do not claim that CLR-NS beats learned methods on human data.
- Cell-type-specific ChIP-seq truths are near random for every method (median AUPRC ratio 0.98–1.02; |r| equals CLR-NS), as BEELINE reported.

**GRNBoost2: order of events, fault and rerun.**  GRNBoost2 was first run and reported with arboreto 0.1.6 under scikit-learn 1.9.1, and those results appeared in our working tables between 16 and 21 October 2026 (Table 2, Supplementary Tables 9 and 10, and a progress checkpoint).  In those runs GRNBoost2 ranked last everywhere, on bulk and single-cell data and on early precision.  On 22 October 2026 we found that arboreto's early-stopping rule is broken under scikit-learn 1.3 and later:  the out-of-bag improvement it monitors changed meaning, so most targets stopped after 25–36 trees.  A diagnostic designed after that discovery (job 0894, scored 26–27 October; arboreto's own `grnboost2()`, the current wrapper, a legacy monitor and no early stopping; three seeds) showed that neither our wrapper nor the seed explained the low scores.  On 27 October 2026, before any rerun, we fixed the rule in writing (commit 050695d, `docs/plans/grnboost2_diagnostic.md`):  all GRNBoost2 results are rerun with legacy early stopping (arboreto's published stopping rule with the pre-1.3 out-of-bag improvement); no early stopping is a sensitivity analysis only; the diagnostic's variant results had been seen when this was chosen; no other comparator, truth, split, seed or scoring script is changed.  Legacy early stopping reproduced stock arboreto 0.1.6 under scikit-learn 1.2.2 on BEELINE mDC_500 (3 seeds; trees per target and AUPRC ratios within one seed SD; commit 45a183b, `docs/grnboost2_equivalence/`).  The rerun (27–31 October 2026) replaced every GRNBoost2 result.  Gate:  in every rerun that recomputed other methods (M3D, EcoMAC, K562 on 5,671 genes and on 8,559 genes, and the 14 BEELINE cases) every non-GRNBoost2 row reproduced the earlier result to within 1e-6 (observed difference 0, `analysis/rerun_gate.R`).  The seed runs on M3D and BEELINE (seeds 1–3) had no gate.  The K562 8,559-gene run was added on 28 October.  The per-regulon comparison was recomputed on 31 October (`results/macro_aupr_v3`).  The superseded outputs are kept in `results/grnboost2_superseded/`.  The arboreto fault concerns the early-stopping monitor under scikit-learn 1.3 and later (observed here with 1.9.1).

Share of off-diagonal pairs tied at zero after symmetrization, first run and rerun:  M3D 94.05% and 45.67%, EcoMAC 93.72% and 29.39%, K562 (5,671 genes) 94.93% and 70.56% (`grnboost2_scores.rds`).  Trees per target in the rerun (median, interquartile range):  M3D 494 (450–535), EcoMAC 839 (747–936), K562 (5,671 genes) 240 (173–310), K562 (8,559 genes) 233 (170–304) (`baselines.log`).

**Supplementary Table 10a | GRNBoost2, first run and rerun: pooled AUPR over the base rate (Table 2 metric).**  The EcoMAC, M3D and K562 5,671-gene rows are those of Table 2; the 8,559-gene K562 run (secondary) is not in Table 2.

| Run | Truth | Raw, first run | Raw, rerun | Hub-controlled, first run | Hub-controlled, rerun |
|---|---|---|---|---|---|
| M3D | SC | 1.18 | 1.36 | 1.28 | 1.37 |
| M3D | all | 1.13 | 1.26 | 1.18 | 1.28 |
| EcoMAC | SC | 1.21 | 1.36 | 1.19 | 1.37 |
| EcoMAC | all | 1.14 | 1.26 | 1.09 | 1.27 |
| K562, 5,671 genes | DoRothEA A–B | 1.02 | 1.05 | 1.02 | 1.05 |
| K562, 5,671 genes | DoRothEA A–C | 1.01 | 1.03 | 1.00 | 1.03 |
| K562, 5,671 genes | CollecTRI | 1.01 | 1.04 | 0.94 | 1.01 |
| K562, 5,671 genes | CORUM | 2.18 | 3.66 | 1.70 | 2.51 |
| K562, 5,671 genes | KD response | 1.37 | 1.84 | 1.08 | 1.55 |
| K562, 8,559 genes | DoRothEA A–B | 1.01 | 1.04 | 1.00 | 1.04 |
| K562, 8,559 genes | DoRothEA A–C | 1.01 | 1.03 | 0.99 | 1.03 |
| K562, 8,559 genes | CollecTRI | 1.01 | 1.04 | 0.94 | 1.01 |
| K562, 8,559 genes | CORUM | 2.21 | 3.76 | 1.72 | 2.61 |
| K562, 8,559 genes | KD response | 1.33 | 1.80 | 1.03 | 1.49 |

Source: `metrics.csv` of the four baseline runs, first run in `results/grnboost2_superseded/`.  Attention (Table 2) is 1.68, 1.44, 1.68 and 1.47 on the four *E. coli* cells, and 3.51 on K562 CORUM against GRNBoost2 3.66 (paired difference GRNBoost2 minus attention +0.0020, jackknife 95% interval −0.0057 to +0.0096).

**Supplementary Table 10b | GRNBoost2 minus attention per regulon, first run and rerun** (`results/macro_aupr_v2` and `results/macro_aupr_v3`; positive values favour GRNBoost2; Holm P for the mean).

| Family | Truth | Regulons | Mean ΔAP, first run | Mean ΔAP, rerun | Holm P, first run | Holm P, rerun | Attention wins, first run | Attention wins, rerun |
|---|---|---|---|---|---|---|---|---|
| *E. coli* M3D | SC | 123 | −0.0269 | −0.0153 | 0.004 | 0.0916 | 71% | 56% |
| *E. coli* M3D | all | 157 | −0.0240 | −0.0117 | 0.004 | 0.184 | 71% | 56% |
| *E. coli* EcoMAC | SC | 119 | −0.0292 | −0.0181 | 0.004 | 0.004 | 68% | 56% |
| *E. coli* EcoMAC | all | 155 | −0.0224 | −0.0114 | 0.004 | 0.0916 | 65% | 54% |
| K562 primary | DoRothEA A–B | 105 | −0.0025 | −0.0005 | 0.664 | 1 | 79% | 66% |
| K562 primary | DoRothEA A–C | 271 | −0.0008 | −0.0002 | 1 | 1 | 77% | 63% |
| K562 primary | CollecTRI | 491 | −0.0016 | −0.0004 | 0.132 | 1 | 69% | 60% |
| K562 primary | CORUM | 499 | −0.0022 | +0.0013 | 0.164 | 1 | 46% | 36% |
| K562 primary | KD response | 566 | −0.0178 | +0.0348 | 0.0127 | 0.0045 | 88% | 9% |
| K562 secondary | DoRothEA A–B | 105 | −0.0021 | −0.0001 | 0.729 | 1 | 79% | 66% |
| K562 secondary | DoRothEA A–C | 269 | −0.0008 | −0.0001 | 0.306 | 1 | 82% | 72% |
| K562 secondary | CollecTRI | 490 | −0.0010 | −0.0005 | 0.141 | 1 | 72% | 64% |
| K562 secondary | CORUM | 499 | −0.0018 | +0.0011 | 0.254 | 1 | 44% | 38% |
| K562 secondary | KD response | 585 | −0.0173 | +0.0290 | 0.00903 | 0.0035 | 90% | 12% |

**Supplementary Table 10c | GRNBoost2 on BEELINE (primary truths, 44 cells, 7 datasets), first run and rerun** (`results/grnboost2_superseded/beeline_baselines/summary.csv`, `results/beeline_baselines/summary.csv`).

| Statistic | First run | Rerun, seed 20261021 | Rerun, seeds 1 / 2 / 3 |
|---|---|---|---|
| Median AUPRC ratio | 1.26 | 1.47 | 1.456 / 1.485 / 1.483 |
| Mean log2 vs CLR | −0.353 | −0.152 | −0.157 / −0.157 / −0.159 |
| Datasets better than CLR | 0/7 | 1/7 (hESC) | 1 / 1 / 1 |
| Wilcoxon P, sign P (vs CLR) | 0.0156, 0.0156 | 0.0312, 0.125 | 0.0312, 0.125 (all three) |
| Mean log2 vs CLR-NS | −0.415 | −0.213 | −0.219 / −0.218 / −0.221 |
| Datasets better than CLR-NS | 0/7 | 0/7 | 0 / 1 / 1 |
| Wilcoxon P (vs CLR-NS) | 0.0156 | 0.0156 | 0.0156 / 0.0312 / 0.0312 |

**PIDC.**  On human single-cell data, PIDC falls back from a gamma to a normal approximation for under 1% of pairs, including pairs involving SOX3 and PAX6.

## Supplementary Note 14 | Hard FDR operators inside the cluster bootstrap

CLR scores carry a permutation null, so hard selection at a false-discovery rate is the natural alternative to a softmax.  Hard selection alone isolates 62–77% of genes, so the operators tested add each gene's ten strongest neighbours (FDR 0.05 + top 10, at its held-out t* = 7) or read out the absolute correlation of profiles diffused over an FDR 0.20 graph (|cor(Pᵗ Z)| at t = 8 and 32; no held-out t*).  Every draw of the cluster bootstrap (200 draws, experiments resampled with replicate arrays, same seeds as the main run so the draws pair) recomputes its own permutation null (`analysis/chips_bootstrap.R --fdr_B`).

**Supplementary Table 11 | Hard FDR operators against CLR** (cluster bootstrap, 200 draws; mean and 2.5–97.5% percentiles; marginal intervals).

| Evidence | Operator | Δ AUPR vs CLR | Δ coherence vs CLR |
|---|---|---|---|
| SC | FDR 0.05 + top 10, t = 7 | +0.012 [+0.010, +0.015] | +0.002 [−0.019, +0.024] |
| SC | FDR 0.20, diffused-profile \|r\|, t = 8 | +0.001 [−0.002, +0.004] | −0.023 [−0.041, −0.004] |
| SC | FDR 0.20, diffused-profile \|r\|, t = 32 | +0.004 [+0.001, +0.008] | −0.022 [−0.042, −0.006] |
| All | FDR 0.05 + top 10, t = 7 | −0.001 [−0.003, +0.002] | −0.009 [−0.026, +0.006] |
| All | FDR 0.20, diffused-profile \|r\|, t = 8 | −0.003 [−0.006, +0.000] | −0.023 [−0.039, −0.007] |
| All | FDR 0.20, diffused-profile \|r\|, t = 32 | −0.003 [−0.006, +0.001] | −0.025 [−0.043, −0.007] |

**Reading.**  The best hard operator gains over CLR on SC regulons (+0.012) without changing coherence or helping all-evidence regulons.  Paired with the same draws, soft attention (keff 10, t = 8) beats it by 0.0075 AUPR (0.0053 to 0.0095), in all 200 draws, with 0.026 higher coherence.  Diffused-profile readouts lower coherence below CLR.

**Caveats.**  The intervals are marginal over about 12 comparisons, and the small FDR 0.20 effects do not survive Bonferroni correction.  The FDR 0.20 depths were not held out.  The fraction of isolated genes was not logged in this run.  Point estimates on the full data are reported separately for the 907-array (+0.018) and 466-experiment (+0.013) sets rather than averaged.

## Supplementary Note 15 | Preregistered analyses [PENDING]

**Attention width (`docs/preregistration_keff.md`; fixed 2026-10-16, implementation choices logged 2026-10-21 before any full run).**  Three data-driven width rules were reviewed before running and each was expected to fail:  one global temperature on CLR z-scores gives a bimodal width split (29% of rows with fewer than 3 partners, 16% near-uniform) that tracks the top score rather than regulon breadth; per-gene widths from local false-discovery rates violate the method's assumptions (a point mass at 0, a heavy tail, dependent rows); and width chosen on held-out data is biased toward wide attention (keff 25 > 10 > 3 on one split).  They are being tested empirically anyway, next to a keff sweep (3, 5, 25) and a multi-head operator (walks at keff 3, 10 and 30 averaged), each at its own held-out depth.

- *Primary endpoint:*  difference in pooled AUPR from keff = 10 on M3D and EcoMAC transcription-factor regulons and on K562 CollecTRI and CORUM; paired delete-half gene jackknife; Holm over 8 variants × 4 endpoints = 32 tests; smallest effect of interest ±0.005.
- *Decision rule:*  the headline changes only if a variant wins by a Holm-significant margin of at least 0.005 on at least three of the four endpoints, loses on none and replicates on both M3D and EcoMAC (which are not independent; Note 17).  The multi-head operator becomes the main operator if it gains on CORUM without losing on regulons.
- *Protocol:*  held-out selections are computed and committed before any truth is read; the truth phase refuses to run unless they match the committed copies.  Selections for M3D and EcoMAC are committed (`docs/keff_selections/`).  Smoke tests on 300-gene subsets scored truths before the change log was fixed and are disclosed there.
- *A degeneracy found before scoring:*  on the CLR scale the two-group local-FDR model assigns local FDR ≈ 1 to most pairs (83% of rows at the floor of 2 partners), and on raw mutual information 100% of rows reach the ceiling of 30.  The local-FDR variant runs as preregistered; if it meets the decision rule it will be reported as a sharp-attention (about keff 2) result, not as a calibrated width.

**Learned attention (`docs/preregistration_scfm.md`; fixed 2026-10-21, before any foundation-model network was computed).**  The primary contrast is CLR-NS against zero-shot scPRINT v1 attention (mean over layers, heads and cells; special tokens removed; rows renormalized; target as query and transcription factor as key) on the seven BEELINE datasets, all cells, scored on genes that every method covers.  The primary hypothesis is non-inferiority with a margin of 0.15 log2 units of AUPRC ratio, read from the exact Wilcoxon (Hodges–Lehmann) interval of the per-dataset differences (coverage 92.2%), with superiority tested next in a fixed sequence.  Key secondary contrasts (CLR against scPRINT; CLR-NS against scPRINT on early precision) are confirmatory only if the primary holds, with Holm over two tests.  scGPT, readouts other than the primary one, the controls (random weights, permuted expression, data-free gene embeddings), CLR applied to learned attention, fine-tuning and K562 control cells are descriptive.  A leakage audit of each model's pretraining corpus against the evaluation accessions precedes scoring; leakage would favour the foundation models and so make a non-inferiority result conservative.  CLR-NS was chosen on these truths, which the preregistration discloses.

## Supplementary Note 16 | Sensitivity of Table 2 to the tie convention

**Why.**  Several scores contain many exact ties, and a reviewer asked whether the tie convention of the pooled AUPR changes Table 2.  We ran the check on 31 October 2026.  It is a sensitivity analysis:  the primary convention is unchanged, was fixed before this check, and every number in Table 2 and the main text is as reported.

**Conventions (`analysis/tie_check/`; results in `results/tie_check/`).**  *Existing* (the convention of `ap_fast` in `analysis/regulondb.R`, used throughout):  every pair in a block of tied scores receives the precision at the end of the block, which is pessimistic for methods with ties.  *Expected ties*:  the exact expectation of average precision when the order within each tied block is uniformly random (a hypergeometric number of earlier positives within the block, evaluated analytically, no Monte Carlo).  *Optimistic* (positives first in every block) and *strict* (negatives first) are the extreme orders and are reported as bounds only.  The marks (bold, †) were recomputed under each convention by the rule of `analysis/sota_table.R` (paired delete-half jackknife, the same 100 splits).

**The existing convention is reproduced exactly.**  Rescoring from the cached scores gives every cell of the pooled AUPR and of the 100 jackknife replicates of the existing convention to within 5e-16 of `metrics.csv` (172 method-truth-dataset checks).

**Where ties occur.**  Over the benchmark pairs (pairs between annotated genes of a truth), attention, GENIE3, PIDC and |r| have no ties, except that attention has 2 to 30 tied positives on the K562 truths (4 to 33 on the 8,559-gene run); the hub-only score has at most 2.2% tied positives.  CLR has 43.5 to 49.7% of benchmark pairs tied (43.0 to 49.7% including the 8,559-gene K562 run), one large block of near-zero scores containing about 45% of the positives on the *E. coli* and K562 transcription-factor truths (27 to 31% on CORUM and knockdown-response modules).  GRNBoost2 has 28 to 71% of benchmark pairs tied (up to 79% on the 8,559-gene K562 run).

**Result under expected ties.**  No cell of Table 2 (raw or hub-controlled) changes by more than 0.11 in ratio over the base rate (largest:  GRNBoost2 on K562 CORUM, 3.76 to 3.87 on the 8,559-gene run and 3.66 to 3.74 on the 5,671-gene run); the ratios of GENIE3, PIDC, |r| and attention are unchanged.  The changes that alter a mark or a rank are:

- *Marks (5,671-gene K562):*  GRNBoost2 gains a † on DoRothEA A–B and A–C; PIDC on knockdown-response modules loses its †.  On the 8,559-gene run CLR gains a † on DoRothEA A–B and |r| on knockdown-response modules loses its †.
- *Rank swaps without a change of mark:*  CLR and PIDC on M3D SC (CLR 4th to 3rd); CLR and GRNBoost2 on EcoMAC SC after hub control; GRNBoost2 on K562 CollecTRI from 5th to 3rd (raw); and a number of swaps among near-tied K562 columns (all changes are listed in `results/tie_check/table2_changes_vs_old.csv`).
- *Knockdown-response modules:*  CLR rises from 3rd to 1st on the 5,671-gene run (2.42 to 2.46), with GENIE3 and |r| shifting down by one place; all three are marked best in both conventions.
- *K562 CORUM:*  the order is unchanged (CLR 6.39 to 6.42; attention 3.51).
- *Attention against CLR (paired jackknife, pooled AUPR):*  on M3D strong-evidence regulons the difference goes from +0.0220 (95% interval +0.0159 to +0.0281; ratio 1.162) to +0.0205 (+0.0144 to +0.0265; ratio 1.149), so the gain remains clearly significant.  Every other attention-over-CLR difference that excluded 0 still does (for example EcoMAC SC +0.0278 to +0.0263; K562 DoRothEA A–B +0.0064 to +0.0056, whose lower limit moves from +0.0013 to +0.0005), the K562 CORUM loss is unchanged (−0.0382 to −0.0385), and the knockdown-response difference includes 0 under both conventions (−0.0039 to −0.0065).

**The extreme conventions are not credible estimates.**  Because the largest block of CLR's scores contains about 45% of its positives and 45 to 50% of all benchmark pairs, the order within it is decisive, and the two extremes differ greatly:  with positives first CLR on M3D SC reaches 1.65 (attention 1.68; on M3D all evidence levels 1.56 against attention's 1.44), with negatives first 1.35.  Neither order is realistic, since the positives are not systematically first or last among pairs with an identical score.  We report them as bounds on what the tie-breaking rule could do, not as alternatives to the reported numbers.

**Raw Table 2 under both conventions** (pooled AUPR over the base rate; bold and † as in Table 2, recomputed per convention; the hub-controlled version and the two extreme conventions are in `results/tie_check/summary.md`).

*M3D*

| Truth | Convention | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|
| SC | existing | **1.68** | 1.44 | 1.51 | 1.45 | 1.36 | 1.23 | 1.08 |
| SC | expected ties | **1.68** | 1.46 | 1.51 | 1.45 | 1.38 | 1.23 | 1.08 |
| all | existing | **1.44** | 1.37 | **1.40**† | **1.40**† | 1.26 | 1.23 | 1.22 |
| all | expected ties | **1.44** | 1.38 | **1.40**† | **1.40**† | 1.27 | 1.23 | 1.22 |

*EcoMAC*

| Truth | Convention | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|
| SC | existing | **1.68** | 1.39 | 1.50 | 1.43 | 1.36 | 1.18 | 1.04 |
| SC | expected ties | **1.68** | 1.40 | 1.50 | 1.43 | 1.37 | 1.18 | 1.04 |
| all | existing | **1.47** | 1.33 | 1.38 | 1.39 | 1.26 | 1.19 | 1.20 |
| all | expected ties | **1.47** | 1.35 | 1.38 | 1.39 | 1.26 | 1.19 | 1.20 |

*K562 (5,671 genes)*

| Truth | Convention | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|
| DoRothEA AB | existing | **1.09** | **1.06**† | **1.08** | **1.07** | 1.05 | **1.05**† | 1.00 |
| DoRothEA AB | expected ties | **1.09** | **1.06**† | **1.08** | **1.07** | **1.06**† | **1.05**† | 1.00 |
| DoRothEA ABC | existing | **1.07** | 1.03 | **1.04**† | **1.04**† | 1.03 | 1.03 | 1.02 |
| DoRothEA ABC | expected ties | **1.07** | 1.03 | **1.04**† | **1.04**† | **1.04**† | 1.03 | 1.02 |
| CollecTRI | existing | **1.10** | 1.03 | 1.04 | 1.05 | 1.04 | 1.01 | **1.14** |
| CollecTRI | expected ties | **1.10** | 1.03 | 1.04 | 1.05 | 1.05 | 1.01 | **1.14** |
| CORUM | existing | 3.51 | **6.39** | 6.32 | **6.33** | 3.66 | **7.00** | 1.55 |
| CORUM | expected ties | 3.51 | **6.42** | 6.32 | **6.33** | 3.74 | **7.00** | 1.55 |
| KD response | existing | **2.38** | **2.42** | **2.45** | **2.35**† | 1.84 | **2.45** | **2.42** |
| KD response | expected ties | **2.38** | **2.46** | **2.45** | 2.35 | 1.91 | **2.45** | **2.42** |

*K562 (8,559 genes; no PIDC)*

| Truth | Convention | Attention | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | Hub only |
|---|---|---|---|---|---|---|---|---|
| DoRothEA AB | existing | **1.10** | 1.06 | **1.08** | - | 1.04 | 1.05 | 1.01 |
| DoRothEA AB | expected ties | **1.10** | **1.07**† | **1.08** | - | 1.05 | 1.05 | 1.01 |
| DoRothEA ABC | existing | **1.08** | 1.04 | **1.05** | - | 1.03 | 1.03 | 1.03 |
| DoRothEA ABC | expected ties | **1.08** | 1.04 | **1.05** | - | 1.04 | 1.03 | 1.03 |
| CollecTRI | existing | **1.13** | 1.04 | **1.06**† | - | 1.04 | 1.02 | **1.15** |
| CollecTRI | expected ties | **1.13** | 1.04 | **1.06**† | - | 1.06 | 1.02 | **1.15** |
| CORUM | existing | 3.72 | **6.75** | **6.63** | - | 3.76 | **7.00** | 1.54 |
| CORUM | expected ties | 3.72 | **6.78** | **6.63** | - | 3.87 | **7.00** | 1.54 |
| KD response | existing | **2.54** | **2.55** | **2.55** | - | 1.80 | **2.47**† | **2.47** |
| KD response | expected ties | **2.54** | **2.58** | **2.55** | - | 1.90 | 2.47 | **2.47** |

---

## Supplementary Note 17 | EcoMAC and M3D share arrays: overlap check and a non-overlap sensitivity analysis (post hoc)

**Why this was checked.**  EcoMAC (Carrera et al. 2014; 2,198 RMA-normalized arrays, 4,189 genes, two Affymetrix platforms) and M3D (907 chips from 466 experiments, 4,297 genes) are both compiled largely from public *E. coli* Affymetrix data, so EcoMAC cannot be called an independent replication of M3D without an audit.  M3D carries no GEO accessions (0 matches in chip names, experiment names and descriptions), so shared CEL files cannot be proved from identifiers; the audit rests on expression similarity and cannot distinguish identical CEL files from independent replicates of identical conditions by the same laboratories (29 author surnames are shared between the two compendia).  No truth file was opened (`analysis/ecomac_m3d_overlap.R`, `results/ecomac_m3d_overlap/`, `docs/plans/ecomac_m3d_overlap.md`).

**Overlap audit.**  On the 4,056 genes shared by the two compendia, we computed the Pearson correlation of per-array standardized log2 profiles between every EcoMAC array and every M3D chip.  The primary rule (r > 0.99, with 0.98 and 0.995 as sensitivity checks) flags 406 EcoMAC arrays (18.5%), and 0.98 flags 815 (37.1%); 0.995 flags none because the two files are normalized differently (the maximum raw correlation is 0.9934), so it is uninformative.  Raw correlation alone is a weak discriminator, because similar conditions within one compendium also reach 0.98.  A second rule used profiles centred on each compendium's gene means and flagged arrays above the maximum (0.7036) seen in a time-matched null, the 518 EcoMAC arrays from series GSE20000 or later, released after the September 2009 build of M3D and therefore not in it; it flags 904 arrays (41.1%) and none in the null.  **Post hoc additions.**  After inspecting these distributions we added the requirement that the EcoMAC array be the mutual nearest neighbour of its best M3D chip (the array is also that chip's best EcoMAC match).  This was not in the rules fixed before computing, and it is the headline call: 828 of 2,198 EcoMAC arrays (37.7%) are mutual nearest neighbours, all of them inside the centred rule, and none of the 518 post-2009 arrays is.  The one-to-one matched pairs cover 91% of M3D chips (828 of 907) and 92% of its experiments (431 of 466).  The calls are whole series: of 127 EcoMAC series or sources, 35 have a called array and 29 are called in at least 90% of their arrays, each array matched to a different M3D chip (for example GSE6836 266 of 284; GSE10158 43 of 43; GSE1121 43 of 43).  By group, 612 of 844 GPL199 (Antisense) arrays are called (86.6% of the 707 with GSE below 20000), 216 of 619 arrays without a GEO record, and none of the 735 GPL3154 (Genome 2.0) arrays or of the post-2009 arrays.  The 95 pre-build GPL199 arrays that miss the rule (median raw r 0.976) are neither confirmed overlap nor confirmed independent, and the 79 M3D chips with no EcoMAC counterpart are unexplained.  EcoMAC is therefore a larger, partly overlapping compendium (about 38% of its arrays and about 91% of M3D's chips are the same or matched hybridizations; it adds Genome 2.0 and post-2009 arrays and batch handling by ComBat per platform), not an independent replication.

**Re-scoring without M3D-matched arrays (post hoc sensitivity analysis, not preregistered).**  We re-ran the identical pipeline on the 1,370 EcoMAC arrays without an M3D counterpart (ComBat by platform re-fitted within the subset, pooled AUPR, K10 attention at the transferred t = 9, gene-jackknife intervals with 100 halves).  The subsets were defined by similarity to M3D, not by outcome, and the calls were made before any of these results existed; no subset was dropped or added after seeing results, except that a random subset of the same size and the label-free depth runs were added at the outset.  The analysis reuses the RegulonDB truth of Table 2, so it cannot serve as a held-out test.  A gate on the full set reproduced the existing EcoMAC table exactly (maximum absolute difference 0).

| EcoMAC subset (arrays) | CLR | Attention | Attention − CLR (SC), 95% gene-jackknife interval | Hub-controlled difference (SC) |
|---|---|---|---|---|
| all (2,198) | 0.130 | 0.158 | +0.0278 [+0.0213, +0.0343] | +0.0195 [+0.0140, +0.0251] |
| no M3D counterpart (1,370) | 0.129 | 0.157 | **+0.0285 [+0.0218, +0.0353]** | +0.0197 [+0.0141, +0.0254] |
| random 1,370 of all (reference) | 0.130 | 0.157 | +0.0273 [+0.0209, +0.0337] | +0.0184 [+0.0128, +0.0240] |
| M3D-matched only (828) | 0.134 | 0.157 | +0.0228 [+0.0167, +0.0289] | +0.0150 [+0.0099, +0.0200] |
| no counterpart, Genome 2.0 only (735) | 0.136 | 0.159 | +0.0230 [+0.0170, +0.0291] | +0.0167 [+0.0114, +0.0221] |
| no counterpart, GPL199 plus no GEO record (635) | 0.114 | 0.133 | +0.0197 [+0.0143, +0.0251] | +0.0108 [+0.0064, +0.0152] |

Strong/confirmed (SC) regulons, base rate 0.0939; |r| scored 0.106 on the 1,370 non-matched arrays (0.110 on all).  In every subset attention exceeded CLR, and CLR exceeded |r|, on both truths.  At all evidence levels the differences were +0.0164 (all), +0.0193 (no counterpart) and +0.0083 (matched only; hub-controlled +0.0030, interval including 0).  The label-free depth, chosen on each subset as for the frozen EcoMAC selection (20 author splits, criterion A, no truth read), was t* = 8 for all arrays and for the non-matched arrays (the value of the frozen selection), 9 for the matched arrays and 11 for Genome 2.0 only; at each subset's own t* the SC difference on the non-matched arrays was +0.0271 [+0.0207, +0.0336].  No conclusion changed with depth.  Intervals resample genes with the networks held fixed; arrays are not resampled and there is no correction for the many descriptive cells.

**Reading.**  The EcoMAC gain of attention over CLR does not depend on the arrays that M3D also contains: on the 1,370 arrays without a counterpart it is the same size as on all arrays and on a random subset of the same size, and if anything the matched arrays are the less favourable part.  It is not a Genome 2.0 platform effect either.  This is a robustness statement about array dependence, not a held-out replication: the truth, genes and pipeline are the same, many laboratories contribute to both compendia, and near-duplicates among the unmatched arrays cannot be excluded.  The main text therefore treats M3D and EcoMAC as two partly overlapping compendia, not as two independent confirmatory sets.

---

## Supplementary Note 18 | A preregistered held-out test of reciprocity in yeast that could not be run

**Purpose.**  Reciprocal attention (keeping attention from gene i to gene j only when each is among the other's softmax candidates) helped small complexes in synthetic data but has never been scored on real data.  It was therefore preregistered (`docs/preregistration_reciprocity_yeast.md`, fixed in commit 02cbcfd, 1 November 2026) as a held-out test on *Saccharomyces cerevisiae*, an organism and truths not used elsewhere in this project, using the Kemmeren et al. (2014) deleteome (GSE42527 mutants and GSE45115 wild-type cultures; 6,169 genes by 652 mutants), with a second arm testing whether attention beats CLR on regulon co-membership (H1).  The protocol had a sanity gate (gate F) that had to pass before any truth was downloaded or opened; if it failed, the preregistered outcome was verdict 1, "not run: no adequately powered held-out test of reciprocity was possible", with nothing else run.

**What happened.**  The deleteome's deposited orientation metadata proved inconsistent.  By each sample's own value definition every GSE42527 value is already log2(mutant/reference), yet the published matrix agrees with the values as deposited only for mutants of the older sub-series and, for most mutants of the `[hs1991]` sub-series, only after negating hybridizations with the reference in channel 1; some mutants fit neither.  No single global orientation rule reproduced the published matrix (best median r 0.793, below the preregistered floor of 0.8).  Two deviations were logged, each before any truth was downloaded or opened, and both are in the preregistration's change log.  Deviation 1 (743d24e) applied the structural dye-balance check to balanced groups only (15 of 701 mutants lack one hybridization of each orientation).  Deviation 2 (25a5a8d) replaced the global rule by a per-sub-series rule fitted on a fixed 50% fitting set of balanced mutants (seed 20261106) and evaluated on the held-out half, with the same floor.  The held-out median r was then 0.820 (0.928 in sub-series `none`, 0.422 in `[hs1991]`), so profiles were rebuilt (384 mutants) and a wild-type matrix of 169 cultures (50.8% of balanced wild-type pairs passed the orientation pair check) was built.

**Gate F failed.**  The gate passes if the median chance-corrected reproducibility R′ of the mutant network (halves of mutants) exceeds that of the wild-type network (halves of wild-type cultures) and the share of mutant top-1% pairs also in the wild-type top 1% is at most 0.5.  The median R′ was 0.122 for mutants (0.093 to 0.161 across the 10 draws) against 0.296 for wild types (0.280 to 0.316), so the first condition failed in all 10 draws; the second held (0.205).  By the fixed preregistration the result is verdict 1.  Gate F ran on rebuilt mutants because the rebuild did not reproduce the published matrix exactly, and misoriented hybridizations cancel signal when averaged, so the failure may partly reflect the orientation problem and not the deleteome network itself; this cannot be tested without changing the gate after seeing its result.  Any rerun, for example on the published matrix, would be a post hoc deviation and exploratory.

**Status.**  No truth was downloaded or opened, no human arm was started, and the yeast truths remain unconsumed (although the design cannot be rerun on them without disclosing this failure).  Reciprocity therefore remains a synthetic-only candidate, and it is not part of the method or of any claim in the main text.  Records: `docs/reciprocity_yeast_selections/` (`step3_stop_report.md`, `step3_stop_report_2.md`, `step4_gateF_fail_report.md`, `feasibility_gate.csv`, `feasibility_draws.csv`, `orientation_rules.csv`, `download_log.md`).

---

## Supplementary Figures

**Supplementary Fig. 1 | Mixing platforms creates the "Christmas tree" and dilutes a known LexA–target relationship.**  Expression (log2, RMA) of lexA against two genes it represses, recA (top row) and sulA (bottom row), in:

1. M3D build 6:  907 arrays from one platform (the 466 experiments used elsewhere are these arrays averaged over replicates);
2. EcoMAC, all 2,198 arrays, coloured and shaped by platform;
3. EcoMAC, the 735 Genome 2.0 arrays only.

Axis ranges are shared along each row.  "No GEO record" marks 619 EcoMAC arrays (M3D and ASAP chip identifiers) without a GEO series, whose platform therefore could not be looked up.  The dashed 90% normal ellipse outlines arrays annotated with a DNA-damaging treatment (quinolone, mitomycin C, UV or nalidixic acid):  the SOS-induced cluster at high lexA.  Each panel gives n, Pearson r and B-spline MI in bits.  For recA, r is 0.85, 0.74 and 0.55 and MI 0.71, 0.46 and 0.10 bits across the three columns; for sulA, r is 0.76, 0.52 and 0.35 and MI 0.57, 0.28 and 0.08 bits.  The Genome 2.0 arrays include only one DNA-damage array, so their weak relation reflects the missing SOS induction as well as the single platform.  `figures/scatter_lexA.png` (`analysis/ecomac_scatter.R`).

**Supplementary Fig. 2 | Excess MI over correlation comes from mixing platforms (8.4% of pairs vs 0.2–0.3%).**  Two-dimensional histograms of B-spline MI (shared histogram bins; bits) against |Pearson r| for every gene pair in M3D (466 experiments, one platform; 9,229,956 pairs), EcoMAC uncorrected and EcoMAC with per-platform centring (2,198 arrays; 8,771,766 pairs each), on the full MI range.  Colour gives gene pairs per bin on a log scale.  The dashed curve is the MI expected for a linear (Gaussian) relationship with the same correlation, −½ log₂(1 − r²).  The large number in each panel is the share of pairs whose MI exceeds that curve by more than 0.1 bits:  0.2% (M3D), 8.4% (EcoMAC uncorrected) and 0.3% (EcoMAC centred); Spearman correlations between MI and |r| are 0.76, 0.77 and 0.77.  Mass below the curve (4.7% of M3D pairs by more than 0.1 bits) and the floor near r = 0 are estimator biases (bin saturation and finite-sample bias), not weaker or extra dependence.  `figures/mi_vs_rho.png` (`analysis/mi_vs_rho.R`).

**Supplementary Fig. 3 | Component contributions at all evidence levels.**  Shapley decomposition as in main-text Fig. 3, for RegulonDB regulons at all evidence levels.  The two evidence levels disagree in sign (with an interval excluding 0 at one level at least) only on co-membership AUPR:  the Shapley share of mutual information (strong/confirmed −0.003, all +0.002) and of shared bins (+0.003, 0.000), and the drop-one losses of shared bins (+0.003, −0.003) and CLR normalization (+0.006, −0.001).  `figures/fig_shapley_all.png`.
