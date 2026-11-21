---
title: "Focus on what matters: summarizing attention in context likelihood of relatedness"
---

**Abstract.** Attention normalizes each query's compatibilities within its own context, an idea older than the transformer (Vaswani et al. 2017), with roots in kernel regression (Nadaraya 1964; Watson 1964) and fast weights (Schmidhuber 1992; Schlag et al. 2021). Context likelihood of relatedness (CLR; Faith et al. 2007), its first widely adopted form in computational biology, remains standard for unsupervised network inference but stops at the scores. We add the rest of the attention operation: a temperature-controlled softmax that converts scores to weights, and a weighted sum of neighbours' values, repeated as a short walk whose depth is chosen on held-out data without labels ({{ACLR_NAME}}). On strong-evidence *Escherichia coli* regulons, {{ACLR_NAME}} raises recovery 16% over CLR and, under hub control, leads GENIE3 and PIDC, but it loses to CLR on tight complexes such as ribosomes. A second CLR on the walk's log ({{LEAD_NAME}}) recovers complexes in K562 (selected among about 20 read-outs). In a preregistered test on HepG2 Perturb-seq it gained 0.027 AUPR over {{ACLR_NAME}} (95% CI 0.011 to 0.042), carried by a few large complexes; the *Escherichia coli* test missed its threshold.

{{LEAD_NAME}} is {{LEAD_DESC}}; like {{ACLR_NAME}}, it has no fitted parameters, and it is our lead candidate.[^names] The appendix has more.

**Evidence for {{LEAD_NAME}}.** In a preregistered test on *E. coli*, {{LEAD_NAME}}'s gain over {{ACLR_NAME}} on protein complexes was +0.026 AUPR on M3D (p = 0.079) and +0.029 on EcoMAC (p = 0.051), which did not reach the preregistered threshold (one-sided p < 0.025; Appendix C). A second preregistered test, on HepG2 Perturb-seq pseudobulk with Complex Portal human complexes, was met. This is a separate fixed-truth result on a new truth, not a replication of the *E. coli* test, which was not met. HepG2 has the same laboratory and screen family as the K562 data on which {{LEAD_NAME}} was selected, so this is a cross-cell-type test and not a cross-technology one. The gain over {{ACLR_NAME}} is +0.027 AUPR (0.0669 against 0.0401; total 95% interval +0.011 to +0.042; one-sided p = 0.0004), above the smallest effect of interest (+0.005). It is carried by complexes of 20 or more genes (89% of it), and three complexes (the cytoplasmic 60S and 40S ribosomes and a proteasome complex) carry 81%. 113 of 211 complexes improve (sign test p = 0.17; the median complex does not), the complex-level interval (−0.014 to +0.067) includes 0, and without ribosomal complexes the gain is +0.004 (−0.002 to +0.010). {{LEAD_NAME}} has the highest pooled AUPR of the seven methods but is not shown better than CLR, |r|, GENIE3 or PIDC (Holm-adjusted one-sided p = 0.15, 0.07, 0.052 and 0.41; PIDC is tied, +0.001); it is ahead of {{ACLR_NAME}} and GRNBoost2. Error bars for GENIE3, GRNBoost2 and PIDC comparisons include gene resampling only and are narrower than they should be. On the yeast Kemmeren deletion compendium with CYC2008 complexes (estimate only, about 10 effectively independent complexes), {{LEAD_NAME}} is +0.057 (0.039 to 0.074) above {{ACLR_NAME}}, carried by the same kind of complex (ribosomes). The size pattern predicted in advance held on HepG2 (gain at 20 or more genes, about zero below 20) and partly on Kemmeren (gain at 20 or more genes; the 7-to-19 class is +0.011, with an interval that includes 0). An earlier exploratory estimate on DREAM5 yeast was +0.096 (Appendix D). Hub-controlled (Figure 3), {{LEAD_NAME}} is above CLR in all 12 rows.[^hub]

**Losses and caveats.** {{LEAD_NAME}} was selected among about 20 read-outs on K562 truths, so its K562 numbers are a selection and not a fit (it has no fitted parameters, but the best of many noisy scores is biased upward). Its depth (t = 9) was set for aCLR. A depth and width search over a grid, run on truths already used, found a flat response (changing either moves {{LEAD_NAME}}'s AUPR by at most about 10% relative, even when chosen with the answer in hand), so the setting is unchanged. A variant that sets depth from the noise level of the data ({{LEAD_NAME}}2) matched {{LEAD_NAME}} within 5% on four of the five used truths but was about 8% lower on HepG2, and has not been tested on fresh data; on the scoreboard it takes the highest printed value from {{LEAD_NAME}} in 8 rows and is below {{LEAD_NAME}} on several non-pooled metrics (Figure 2). CLR and GENIE3 are ahead of both methods on direct TF-to-target edges in *E. coli* (AUPR ratio: CLR 6.75, GENIE3 6.54, {{LEAD_NAME}} 5.84, {{ACLR_NAME}} 5.17). {{ACLR_NAME}} is below CLR on K562 protein complexes (3.51 against 6.39; 86% of the loss comes from ten ribosomal and pre-rRNA complexes), and per regulon it ties CLR and GENIE3. Every positive result for {{ACLR_NAME}} comes from the *E. coli* development lineage (M3D with RegulonDB, and the overlapping EcoMAC), where its settings were chosen; in a preregistered DREAM5 yeast test its co-membership gain over CLR (5.8%) stayed below the 10% margin. All methods are within 20% of random on K562 transcription-factor truths. On the two new complex truths {{LEAD_NAME}} is first in pooled AUPR on both, but it is tied with PIDC on HepG2. HepG2 shares its laboratory and screen family with the K562 tuning data. 418 of the 2,040 HepG2 truth genes are themselves perturbed in the panel (the gain is also present in pairs where neither gene is perturbed: +0.014, +0.003 to +0.024). The gain is not broad across complexes. Appendix A gives the remaining caveats.

{{PAGEBREAK}}

![**Figure 1. Graphical abstract.** The 2007 paper describes only the Euclidean combination, sqrt(z~ij~² + z~ji~²); the released 2007 software offered it and the Stouffer combination, (z~ij~ + z~ji~)/sqrt(2), which is the one used here. Dotted green = new in this work. {{LEAD_NAME}}2 is exploratory. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p1.png){width=6.36in}

{{PAGEBREAK}}

![**Figure 2. Scoreboard, raw scores** (pooled AUPR over base rate unless a row says otherwise). A bold dot marks a row's highest printed value, not a statistical winner; no intervals or tests are shown. {{LEAD_NAME}}2 is exploratory and outside the preregistration, with depth set by the data, and is lower than {{LEAD_NAME}} on several non-pooled metrics. Details and per-dataset depths: figure footnote, Appendix. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p3.png){width=6.69in}

{{PAGEBREAK}}

![**Figure 3. Hub-controlled scoreboard.** The rows of Figure 2 that have a hub-controlled version. Hub-controlled: each method's scores are converted to percentiles within 10 strata of the attention walk's stationary distribution (the hub score of a gene pair, pi~i~ + pi~j~), so a method is not rewarded merely for favouring hubs. Marks as in Figure 2; {{LEAD_NAME}}2 is exploratory. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p3b.png){width=7.0in}

# Potential follow-ups and decision points

- **Perturbation targets.** Attention-based influence ranks a perturbed transcription factor's targets above differential expression in 12 of 15 perturbations (p = 0.015); under a laboratory-matched null the figure is 9 of 15 (p = 0.095). The result is suggestive and not established.
- **Differential CLR.** If perturbations are pursued, differential CLR is a possible extension. It is planned and not yet ported, and there is no result.
- **Decision points.** (1) Which method is presented as the headline: {{ACLR_NAME}}, {{LEAD_NAME}}, or both with their trade-offs. (2) How much weight the breadth qualifications carry in the complex result: it is worded as a separate fixed-truth result on a new truth (the HepG2 test was met), not as a replication of the *E. coli* test, which was not met, and the gain is carried by a few large complexes. (3) The naming of the methods.

{{PAGEBREAK}}

# Appendix

**Notes to Figures 1, 2 and 3.**

*Figure 1 (continued).* Middle, the additions outlined in dotted green: each gene attends over its own top 50 positive CLR scores, a per-gene temperature sets about 10 effective partners, and each walk step follows the rows of the attention matrix; {{LEAD_NAME}} then adds the log of the walk and a second CLR. Right, pooled precision-recall area as a multiple of random; the dashed bars are {{LEAD_NAME}} (lead candidate) and {{LEAD_NAME}}2, which is exploratory and outside the preregistration: {{LEAD_NAME}} with the walk depth set by the data instead of t = 9.

*Figure 2 (continued).* Every metric we have for CLR, {{ACLR_NAME}}, {{LEAD_NAME}}, {{LEAD_NAME}}2 and the competitors. Bars start at random and are scaled within each row. A marked (bold) value is the highest printed value in its row, not a statistical winner, and no intervals or tests are shown; equal printed values share the mark. The bottom row counts, for each method, the rows (of 21) where it has the highest printed value or ties it; every method was run on every row of this figure. Dotted green outlines mark the human-data groups. {{LEAD_NAME}} was picked among about 20 read-outs by its K562 scores, and K562 and M3D are tuning sets for it. The last two rows (HepG2 and Kemmeren) come from later preregistered runs, on which {{ACLR_NAME}} and {{LEAD_NAME}} were not tuned; on HepG2, {{LEAD_NAME}} is not shown better than CLR, |r|, GENIE3 or PIDC (Appendix C), and the Kemmeren row is an estimate only. {{LEAD_NAME}}2 is exploratory and outside the preregistration: {{LEAD_NAME}} with the depth set by the data (t2 per data set in the figure), not t = 9. It carries no statistical-tie marks in the figure (Table A1 has them). With its column present it holds the highest printed value in 8 of the 21 rows, all of which {{LEAD_NAME}} would hold without it ({{LEAD_NAME}}: 6 rows with the column, 14 without), and it is lower than {{LEAD_NAME}} on several non-pooled metrics (TF-edge AUPR, per-TF AUROC, per-regulon AP, RPE1, BEELINE medians). The hub-controlled version of the earlier rows is Figure 3; the two new rows are not in it.

*Figure 3 (continued).* Hub control: the walk converges toward well-connected genes (hubs), so each method's scores are converted to percentiles within 10 strata of the attention walk's stationary distribution (the hub score of a gene pair, pi_i + pi_j), and methods are compared within strata. A method is therefore not rewarded merely for favouring hubs. Bars start at random and are scaled within each row, so compare the printed values across a row, not bar lengths across rows. A bold value with a dot is the highest printed value in its row (equal printed values share the mark), not a statistical winner; no intervals or tests are shown. Hub control was not run for direct edges, per-regulon AP and BEELINE; the DREAM5 yeast row is left out because its hub-controlled scores were stored with a different stratifier, and the HepG2 and Kemmeren rows were not hub-controlled. Dotted green outlines mark the human-data groups. {{LEAD_NAME}}2 is exploratory and outside the preregistration ({{LEAD_NAME}} with the walk depth set by the data, not t = 9); it carries no statistical-tie marks in the figure (Table A1 has them), and it takes the highest printed value in 9 of the 12 rows, tied with {{LEAD_NAME}} in three. It is above CLR in 11 of the 12 rows (RPE1 is the exception).

# A. Methods compared, and further results for {{ACLR_NAME}}

**Methods compared.** Components of each method, with its *E. coli* AUPR (the component ablation: strong-evidence regulons, 466 experiments, RegulonDB):

- CLR: mutual information, CLR z; scores ranked directly. AUPR 0.136.
- Attention on |r|: |r|, softmax, walk (t\* = 8); no CLR normalization. AUPR 0.143.
- Attention on mutual information: mutual information, softmax, walk (t\* = 6); no CLR normalization. AUPR 0.143.
- {{ACLR_NAME}}: mutual information, CLR z, softmax, walk (t\* = 9). AUPR 0.158.
- {{LEAD_NAME}}: {{ACLR_NAME}}, then log and a second CLR (t = 9). AUPR 0.162.
- {{LEAD_NAME}}2: {{LEAD_NAME}} with the walk depth set by the data (exploratory). AUPR 0.168.

The two attention rows without CLR normalization show that the softmax alone gains less. Here t\* is the held-out depth, chosen per method; the method comparison itself fixes t = 9, and each ablation value uses its method's own t\*, except {{LEAD_NAME}} (t = 9). {{LEAD_NAME}}2 is exploratory and outside the preregistration: {{LEAD_NAME}} with the walk depth set by the data instead of t = 9 (t2 = 25 on M3D, 29 on K562, 17 per BEELINE case); its *E. coli* value is the pooled AUPR on M3D strong-evidence regulons at t2. {{LEAD_NAME}} was picked among about 20 read-outs by its K562 scores (a selection, not a fit: it has no fitted parameters); its preregistered *E. coli* test did not reach the preregistered threshold, and its preregistered HepG2 test did (Appendix C). On BEELINE (median AUPRC ratio over the 44 case-truth cells of the primary truths, 14 dataset-and-gene-set cases), CLR is plain CLR, and {{ACLR_NAME}} and {{LEAD_NAME}} use the CLR-NS operator (CLR on the absolute correlation of normal scores; {{LEAD_NAME}} at t = 9, not tuned, exploratory). The K562 CORUM and BEELINE values of these methods are in Figure 2. Competitors, BEELINE median ratio: GENIE3 1.69, PIDC 1.56, GRNBoost2 1.47. Against zero-shot attention from the single-cell foundation model scPRINT, CLR-NS is ahead in 7 of 7 datasets (+1.03 log2 AUPR; A.3). Method names are temporary.

**A.1 Method comparison against established methods.**

| Data | Truth | {{ACLR_NAME}} | {{LEAD_NAME}} | {{LEAD_NAME}}2 | CLR | GENIE3 | PIDC | GRNBoost2 | \|r\| | {{LEAD_NAME}} minus CLR (95% CI) |
|:--|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Raw scores** | | | | | | | | | | |
| M3D | Regulons, strong evidence | 1.68 | 1.72 | **1.78** | 1.44 | 1.51 | 1.45 | 1.36 | 1.23 | +0.28 (+0.21, +0.35) |
| M3D | Regulons, all evidence | 1.44 | 1.56 | **1.61** | 1.37 | 1.40 | 1.40 | 1.26 | 1.23 | +0.19 (+0.15, +0.23) |
| EcoMAC | Regulons, strong evidence | 1.68 | 1.72 | **1.77** | 1.39 | 1.50 | 1.43 | 1.36 | 1.18 | +0.33 (+0.26, +0.41) |
| EcoMAC | Regulons, all evidence | 1.47 | 1.58 | **1.64** | 1.33 | 1.38 | 1.39 | 1.26 | 1.19 | +0.25 (+0.21, +0.29) |
| K562 | DoRothEA AB | **1.09**† | **1.11** | 1.08 | 1.06 | **1.08**† | 1.07 | 1.05 | 1.05 | +0.05 (+0.02, +0.08) |
| K562 | DoRothEA ABC | **1.07** | **1.07** | **1.06**† | 1.03 | 1.04 | 1.04 | 1.03 | 1.03 | +0.05 (+0.03, +0.07) |
| K562 | CollecTRI | 1.10 | 1.13 | **1.18** | 1.03 | 1.04 | 1.05 | 1.04 | 1.01 | +0.11 (+0.07, +0.14) |
| K562 | CORUM complexes | 3.51 | **7.95** | **7.83** | 6.39 | 6.32 | 6.33 | 3.66 | **7.00**† | +1.56 (+0.85, +2.27) |
| K562 | KD response modules | 2.38 | 3.13 | **3.55** | 2.42 | 2.45 | 2.35 | 1.84 | 2.45 | +0.70 (+0.61, +0.79) |
| **Hub-controlled (percentiles within 10 strata of the {{ACLR_NAME}} walk's stationary distribution)** | | | | | | | | | | |
| M3D | Regulons, strong evidence | 1.56 | 1.63 | **1.67** | 1.40 | 1.47 | 1.40 | 1.37 | 1.22 | +0.23 (+0.17, +0.28) |
| M3D | Regulons, all evidence | 1.35 | **1.46** | **1.47** | 1.33 | 1.37 | 1.34 | 1.28 | 1.21 | +0.12 (+0.09, +0.15) |
| EcoMAC | Regulons, strong evidence | 1.57 | 1.65 | **1.69** | 1.37 | 1.47 | 1.40 | 1.37 | 1.16 | +0.28 (+0.22, +0.34) |
| EcoMAC | Regulons, all evidence | 1.36 | 1.46 | **1.48** | 1.29 | 1.35 | 1.33 | 1.27 | 1.15 | +0.18 (+0.14, +0.21) |
| K562 | DoRothEA AB | **1.09**† | **1.11** | **1.11** | 1.06 | 1.07 | 1.07 | 1.05 | 1.05 | +0.05 (+0.04, +0.07) |
| K562 | DoRothEA ABC | **1.06** | **1.06** | **1.06** | 1.02 | 1.04 | 1.03 | 1.03 | 1.02 | +0.04 (+0.03, +0.05) |
| K562 | CollecTRI | 1.02 | **1.05** | **1.05** | 1.01 | **1.03**† | **1.03**† | 1.01 | 0.99 | +0.04 (+0.02, +0.06) |
| K562 | CORUM complexes | 2.08 | **4.15** | **4.01** | 3.40 | **3.47**† | **3.82** | 2.51 | **3.76** | +0.75 (+0.26, +1.23) |
| K562 | KD response modules | 1.60 | **2.09**† | **2.14** | 1.91 | 1.93 | 2.00 | 1.55 | 1.88 | +0.18 (+0.12, +0.24) |

: Table A1. Pooled AUPR divided by the base rate (1.00 = random), nine raw and nine hub-controlled rows. Hub-controlled: percentiles within 10 strata of the {{ACLR_NAME}} walk's stationary distribution. Bold: best in the row or statistically tied with it among the eight methods ({{ACLR_NAME}}, {{LEAD_NAME}}, {{LEAD_NAME}}2, CLR, GENIE3, PIDC, GRNBoost2 and |r|; paired gene jackknife, Bonferroni-corrected 95%); † tied only after the correction. Adding {{LEAD_NAME}} and {{LEAD_NAME}}2 removes some marks that the original six methods had: 36 of their 47 marks (18 of 22 in the raw rows, 18 of 25 in the hub-controlled rows) are lost, and none is gained. {{LEAD_NAME}} is the lead candidate; {{LEAD_NAME}}2 is exploratory: its rule was fixed in advance, but the depth 24 and 32 curves of the earlier grid had been seen on these truths. It is {{LEAD_NAME}} with the walk depth set by the data instead of t = 9 (t2 = 25 on M3D, 18 on EcoMAC, 29 on K562, the same in the hub-controlled rows). The last column is the pooled AUPR difference {{LEAD_NAME}} minus CLR divided by the base rate (enrichment units), with its uncorrected 95% jackknife interval; Table A2 has the other paired differences. Intervals resample genes only, so they understate the uncertainty (in the DREAM5 yeast test the experiment-resampling part of the standard error was as large as the gene part or larger), and some bold marks may change. M3D and K562 were tuning sets for the earlier methods. The K562 transcription-factor rows (DoRothEA and CollecTRI) are within 20% of random for every method (the largest value is {{LEAD_NAME}}2 on CollecTRI, 1.18). EcoMAC is not an independent replication of M3D: 828 of its 2,198 arrays (37.7%) have a one-to-one M3D counterpart, covering 91% of M3D's chips. On the 1,370 EcoMAC arrays without an M3D counterpart, {{ACLR_NAME}}'s lead over CLR is unchanged in size (+0.0285 AUPR, against +0.0278 on all arrays; a post hoc sensitivity analysis). GRNBoost2 uses legacy early stopping. Method names are temporary.

Counted from the bold marks of Table A1, {{LEAD_NAME}} is in the bold set in 3 of the 9 raw rows and 6 of the 9 hub-controlled rows, {{LEAD_NAME}}2 in 8 of 9 raw and 9 of 9 hub-controlled, and {{ACLR_NAME}} in 2 of 9 raw and 2 of 9 hub-controlled (K562 DoRothEA AB and ABC in both); at least one of {{LEAD_NAME}} and {{LEAD_NAME}}2 is in the bold set in 9 of 9 raw and 9 of 9 hub-controlled rows. In the paired differences (last column of Table A1 and Table A2), {{LEAD_NAME}} exceeds CLR in all 18 rows with the 95% interval excluding zero, and {{LEAD_NAME}}2 exceeds CLR in 17 of 18 (not separable in K562 DoRothEA AB, raw). {{LEAD_NAME}} exceeds {{ACLR_NAME}} in 15 of 18 rows and is not separable from it in 3 (EcoMAC strong-evidence regulons, raw; K562 DoRothEA ABC, raw and hub-controlled), never significantly lower; {{LEAD_NAME}}2 exceeds {{ACLR_NAME}} in 14 of 18 and is not separable in 4 (K562 DoRothEA AB and ABC, raw and hub-controlled), never significantly lower. These intervals are uncorrected for multiplicity and resample genes only, and they cover pooled co-membership AUPR rows only. On other metrics of Figure 2 (direct TF-target edges, per-regulon AP, per-TF AUROC) {{LEAD_NAME}}2 is lower than {{ACLR_NAME}} in 11 of the 12 rows (the exception is the EcoMAC all-evidence per-TF AUROC, 0.594 against 0.591).

| Data | Truth | {{LEAD_NAME}} minus CLR | {{LEAD_NAME}} minus {{ACLR_NAME}} | {{LEAD_NAME}}2 minus CLR | {{LEAD_NAME}}2 minus {{ACLR_NAME}} |
|:--|:--|:--:|:--:|:--:|:--:|
| **Raw scores** | | | | | |
| M3D | Regulons, strong evidence | +0.28 (+0.21, +0.35) | +0.05 (+0.01, +0.09) | +0.34 (+0.25, +0.44) | +0.11 (+0.04, +0.17) |
| M3D | Regulons, all evidence | +0.19 (+0.15, +0.23) | +0.12 (+0.09, +0.15) | +0.24 (+0.19, +0.29) | +0.17 (+0.12, +0.21) |
| EcoMAC | Regulons, strong evidence | +0.33 (+0.26, +0.41) | +0.04 (−0.00, +0.07)† | +0.38 (+0.28, +0.48) | +0.08 (+0.03, +0.14) |
| EcoMAC | Regulons, all evidence | +0.25 (+0.21, +0.29) | +0.11 (+0.08, +0.14) | +0.31 (+0.25, +0.36) | +0.17 (+0.13, +0.21) |
| K562 | DoRothEA AB | +0.05 (+0.02, +0.08) | +0.02 (+0.00, +0.03) | +0.02 (−0.02, +0.06)† | −0.01 (−0.04, +0.01)† |
| K562 | DoRothEA ABC | +0.05 (+0.03, +0.07) | +0.01 (−0.01, +0.02)† | +0.03 (+0.00, +0.06) | −0.01 (−0.03, +0.01)† |
| K562 | CollecTRI | +0.11 (+0.07, +0.14) | +0.04 (+0.02, +0.06) | +0.15 (+0.09, +0.21) | +0.08 (+0.04, +0.12) |
| K562 | CORUM complexes | +1.56 (+0.85, +2.27) | +4.44 (+2.52, +6.37) | +1.44 (+0.61, +2.27) | +4.32 (+2.33, +6.31) |
| K562 | KD response modules | +0.70 (+0.61, +0.79) | +0.75 (+0.65, +0.85) | +1.12 (+0.98, +1.27) | +1.17 (+1.03, +1.31) |
| **Hub-controlled (percentiles within 10 strata of the {{ACLR_NAME}} walk's stationary distribution)** | | | | | |
| M3D | Regulons, strong evidence | +0.23 (+0.17, +0.28) | +0.07 (+0.04, +0.10) | +0.27 (+0.20, +0.34) | +0.11 (+0.07, +0.15) |
| M3D | Regulons, all evidence | +0.12 (+0.09, +0.15) | +0.11 (+0.08, +0.13) | +0.13 (+0.09, +0.17) | +0.12 (+0.09, +0.14) |
| EcoMAC | Regulons, strong evidence | +0.28 (+0.22, +0.34) | +0.07 (+0.04, +0.10) | +0.33 (+0.25, +0.40) | +0.12 (+0.08, +0.16) |
| EcoMAC | Regulons, all evidence | +0.18 (+0.14, +0.21) | +0.10 (+0.08, +0.12) | +0.19 (+0.15, +0.23) | +0.12 (+0.09, +0.15) |
| K562 | DoRothEA AB | +0.05 (+0.04, +0.07) | +0.02 (+0.01, +0.04) | +0.05 (+0.03, +0.07) | +0.02 (−0.00, +0.04)† |
| K562 | DoRothEA ABC | +0.04 (+0.03, +0.05) | +0.01 (−0.01, +0.02)† | +0.03 (+0.02, +0.05) | +0.00 (−0.02, +0.02)† |
| K562 | CollecTRI | +0.04 (+0.02, +0.06) | +0.03 (+0.02, +0.05) | +0.04 (+0.02, +0.07) | +0.03 (+0.01, +0.05) |
| K562 | CORUM complexes | +0.75 (+0.26, +1.23) | +2.07 (+0.91, +3.24) | +0.61 (+0.12, +1.10) | +1.93 (+0.81, +3.06) |
| K562 | KD response modules | +0.18 (+0.12, +0.24) | +0.49 (+0.42, +0.56) | +0.23 (+0.15, +0.31) | +0.54 (+0.45, +0.64) |

: Table A2. Paired differences in pooled AUPR divided by the base rate (enrichment units), each with an uncorrected 95% gene-jackknife interval, in the 18 rows of Table A1 (the first difference column repeats the last column of Table A1). † The interval includes zero (the methods are not separable in that row); the test uses the unrounded interval, so a printed bound of 0.00 or −0.00 can fall on either side of zero. Intervals resample genes only, so they understate the uncertainty, and they are not corrected for multiplicity (18 rows, four contrasts). {{LEAD_NAME}}2 is exploratory. Method names are temporary.

**A.2 *E. coli* M3D, 466 experiments, RegulonDB regulon co-membership (the headline).**

- {{ACLR_NAME}} beats CLR by +0.022 AUPR [+0.019, +0.026]. Depth is re-selected without labels in each of 1,000 nested bootstrap draws, and {{ACLR_NAME}} is ahead in all 1,000.
- The gain is consistent across 50 laboratory splits (+0.022 [+0.017, +0.027]) and a replicate-aware bootstrap on the 907 chips (+0.020).
- Per regulon (each regulon weighted equally), {{ACLR_NAME}} ties CLR and GENIE3; its win rate is highest in regulons of 25 or more genes.
- Context normalization is essential: a softmax on uncalibrated |r| or MI gets about a third of the gain and none of the coherence gain. A Shapley decomposition over six components attributes the AUPR gain mainly to depth (+0.018) and CLR normalization (+0.013).
- Removing 0 to 20 global principal components leaves the gain intact; the top components are biology (growth rate, motility, the stationary-phase program).
- Hub control: the walk converges toward well-connected genes, so every method is also scored within strata of hubness (Figure 3). {{ACLR_NAME}} stays first on strong-evidence truths (ahead of CLR in every hub decile on M3D) and ties GENIE3 on all-evidence truths. Under a stratifier that does not use the walk (row mean of raw MI; exploratory), it stays ahead of CLR on M3D strong evidence (+0.18 enrichment [+0.12, +0.24]); that covers {{ACLR_NAME}} against CLR only, not GENIE3 or PIDC.

**A.3 Further caveats.**

- **Human data.** On K562 transcription-factor truths every method is within 20% of random. On single cells (BEELINE), CLR on normal scores (CLR-NS) has the highest per-dataset mean of the methods tested, but it was chosen on these truths and is not significantly ahead of GENIE3 or PIDC, and GENIE3 has the higher median; {{ACLR_NAME}} adds nothing there. Against zero-shot attention from the single-cell foundation model scPRINT (preregistered, run once), CLR-NS is non-inferior and superior on all seven BEELINE datasets (Hodges-Lehmann estimate of the log2 AUPRC-ratio difference +1.03, exact interval +0.34 to +1.90), but scPRINT's attention scores the same as the same model run on permuted expression, so the comparator is weak here.
- **Complexes.** GRNBoost2 does about as poorly as {{ACLR_NAME}} on K562 complexes, so the loss is not specific to the walk. {{LEAD_NAME}} recovers it on K562 (7.95) but is below PIDC and CLR on RPE1 merged CORUM (6.14 against 8.11 and 6.77) and below |r| on *E. coli* complexes (26.1 against 29.2). Per regulon, {{LEAD_NAME}} (25.2) is level with {{ACLR_NAME}} (25.6) and above CLR (23.0).
- **Direct edges.** {{ACLR_NAME}} also trails CLR on DREAM5 network 3, which shares the *E. coli* lineage.
- **Yeast, DREAM5 network 4 (preregistered, run once).** For {{ACLR_NAME}} against CLR, the co-membership gain was 5.8% (interval +1.8% to +10.0% with group-level resampling; it includes zero with series-level resampling), below the pre-specified 10% margin: the verdict "equivalent within 10%" sits on the boundary with "inconclusive", and the result is not "superior". The edge read-out was inconclusive. Every method is near chance on yeast (AUROC 0.51 to 0.52 for all seven methods).
- **Where the evidence comes from.** The settings (10 effective partners, 50 candidates, lazy weight 0.5, the MI estimator, the truth construction) were settled on the *E. coli* data after many looks; the nested bootstrap protects only the depth choice. Only the yeast tests lie outside that lineage.

{{PAGEBREAK}}

# B. Walk depth

![**Figure A1. AUPR gain over CLR at each depth, one panel per truth.** Walk depth *t* on a log2 axis from 1 to 256, with the best depth marked, the within-1% interval shaded and the grid edge 256 dashed. On K562 CORUM complexes, CLR is best and every attention depth is below it.](/home/boris/Projects/ClrAttention/results/granularity_log/depth_log_curves.png){width=7.0in}

The held-out depth is stable: median *t* = 9, with the criterion flat within 1% for *t* = 7 to 11.

| Organism | Truth, fine to coarse | Median size (sets) | Best t | AUPR gain over CLR |
|:-----------|:-------------------|:----------:|:------------:|:------------:|
| *E. coli* | 1 operons | 2 (786) | 4 | +0.208 |
| *E. coli* | 2 TF regulons | 11 (158) | 96 | +0.027 |
| *E. coli* | 3 sigmulons | 128 (5) | 12 | +0.043 |
| K562 | 1 CORUM complexes | 4 (1,186) | CLR (0) | 0.000 |
| K562 | 2 KD-response modules | 100 (585) | 24 | +0.017 |
| K562 | 3 DoRothEA AB | 13 (105) | 4 | +0.008 |
| K562 | 4 DoRothEA A-C | 28 (269) | 8 | +0.005 |
| K562 | 5 CollecTRI | 16 (490) | 48 | +0.014 |

: Table B1. Best t is the highest-AUPR point over CLR and the geometric grid t = 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, plus CLR as t = 0. AUPR gain: absolute AUPR at the best depth minus CLR's. K562 CORUM: CLR beats all 16 attention depths; the best attention depth is t = 16, with an AUPR 0.037 below CLR's (0.0894). The pre-specified test of whether best depth rises from fine to coarse truths (Spearman rho at most 0.5 in at least one organism) came out as "not resolution-like": *E. coli* rho = 0.500 (n = 3), K562 rho = 0.700 (n = 5). With 3 and 5 truths and flat depth curves, the test could neither establish nor exclude the ordering.

# C. The preregistered tests of {{LEAD_NAME}}

**C.1 The *E. coli* test.**

- Complex superiority over {{ACLR_NAME}} at *t* = 9 was required at one-sided p < 0.025 on both *E. coli* datasets. It was +0.026 AUPR on M3D (p = 0.079) and +0.029 on EcoMAC (p = 0.051), so the threshold was not reached. The jackknife standard error on complexes was 0.018, and the smallest effect the design could detect was about 0.05.
- Regulon non-inferiority (margin 0.005 AUPR) passed on both: +0.004 and +0.003.
- {{LEAD_NAME}} has no fitted parameters, so it cannot memorise the truth, but it was picked as the best of about 20 read-outs by its K562 scores, and the best of many noisy scores is biased upward. The candidates were highly correlated and the gain is large, so the bias is probably a fraction of the gap, but it has not been measured. On RPE1, another cell line (descriptive only), {{LEAD_NAME}} is +0.036 [+0.013, +0.059] AUPR above {{ACLR_NAME}} on merged CORUM.

**C.2 The HepG2 and Kemmeren runs.**

Both runs followed a plan fixed before any score was computed, and each was run once. Arm A (HepG2) is the confirmatory test: {{LEAD_NAME}} minus {{ACLR_NAME}} in pooled co-membership AUPR (t = 9, raw scores), one-sided p < 0.025, and an interval above the smallest effect of interest (+0.005). Arm B (Kemmeren) is an estimate only, with no p-value and no claim. The intervals combine a delete-half jackknife over genes and one over experiments (100 experiment replicates). {{ACLR_NAME}} and {{LEAD_NAME}} were not tuned on either dataset.

- **Arm A data.** HepG2 essential-gene Perturb-seq of Nadig et al. (2025, *Nature Genetics*; GEO GSE264667), pseudobulk built by the K562 recipe (9,622 genes). Truth: EBI Complex Portal human complexes (Meldal et al. 2022, *Nucleic Acids Research* 50:D578; taxon 9606), merged to 359 complexes of 3 to 200 genes, 17,417 positive pairs among 2,040 genes. It has the same laboratory and screen family as the K562 data on which {{LEAD_NAME}} was selected, so the test is cross-cell-type and not cross-technology.
- **Arm B data.** The yeast deletion compendium of Kemmeren et al. (2014, *Cell* 157:740; GEO GSE42528), 6,169 genes and 1,283 mutant profiles. Truth: CYC2008 (Pu et al. 2009, *Nucleic Acids Research* 37:825), taken from an archived capture of the authors' file, merged to 213 complexes, 10,308 positive pairs among 1,324 genes. The truth overlaps in content with the Complex Portal yeast truth used earlier.
- **Built-in link between data and truth.** In both screens the perturbed genes include subunits of the truth complexes (418 of 2,040 truth genes on HepG2; 340 of 1,324 on Kemmeren), so co-membership is partly built into the design. The gain is also present in pairs where neither gene is perturbed (Table C1).
- **Reachability.** No confinement of the walk's reachable set was seen on either arm: the median gene reaches 9,542 of 9,621 other genes on HepG2 and 6,150 of 6,168 on Kemmeren, and no complex was closed off from the rest.

| Quantity | Arm A: HepG2, Complex Portal (confirmatory) | Arm B: Kemmeren, CYC2008 (estimate only) |
|:-----------------------------|:-------------------|:-------------------|
| {{LEAD_NAME}} minus {{ACLR_NAME}}, pooled AUPR | +0.027 (0.0669 against 0.0401) | +0.057 (0.2076 against 0.1510) |
| Total 95% interval | +0.011 to +0.042 | +0.039 to +0.074 |
| One-sided p (threshold 0.025) | 0.0004 | none (estimate only) |
| Complex-level interval (sensitivity) | −0.014 to +0.067 | −0.009 to +0.122 |
| Effective independent complexes | 34.0 | 10.4 |
| Complexes (10 or more positive pairs) that improve; sign test p | 113 of 211; p = 0.17 | 51 of 96; p = 0.31 |
| Share of the gain from complexes of 20 or more genes | 89% | 95% |
| Share of the gain from the three largest contributors | 81% | 82% |
| Complexes of 3 to 6 genes | +0.0003 (−0.0007 to +0.0014) | −0.004 (−0.013 to +0.005) |
| Complexes of 7 to 19 genes | +0.0014 (−0.0024 to +0.0052) | +0.011 (−0.006 to +0.028) |
| Complexes of 20 or more genes | +0.035 (+0.012 to +0.059) | +0.073 (+0.048 to +0.099) |
| Without ribosomal complexes | +0.004 (−0.002 to +0.010) | +0.035 (+0.015 to +0.055) |
| Pairs in which neither gene is perturbed in the panel | +0.014 (+0.003 to +0.024) | +0.049 (+0.030 to +0.067) |
| Exploratory: {{LEAD_NAME}}2 minus {{ACLR_NAME}}, pooled AUPR | +0.021 (0.0615 against 0.0401) | +0.058 (0.2087 against 0.1510) |
| Exploratory: {{LEAD_NAME}}2 minus {{ACLR_NAME}}, interval and p-value | not computed | not computed |

: Table C1. {{LEAD_NAME}} minus {{ACLR_NAME}}, pooled AUPR difference, with total 95% intervals (gene and experiment resampling). The three largest contributors are the cytoplasmic 60S and 40S ribosomal subunits and a merged proteasome complex (the 19S-20S-PA28 hybrid proteasomes) on HepG2, and the cytoplasmic large and small ribosomal subunits and the mitochondrial large ribosomal subunit on Kemmeren. The size classes are those predicted in advance: a gain at 20 or more genes and about zero below 20. On HepG2 this held. On Kemmeren it held for the 3-to-6 and 20-or-more classes; the 7-to-19 class is +0.011 with an interval that includes 0, so it is not marked as near zero. The last two rows are {{LEAD_NAME}}2, {{LEAD_NAME}} with the walk depth set by the data (t2 = 32 on HepG2, 23 on Kemmeren), exploratory and outside the preregistration: point estimates only. Method names are temporary.

| Method compared | Arm A: pooled AUPR | Arm A: {{LEAD_NAME}} minus method (total 95% interval) | Arm A: Holm-adjusted one-sided p | Arm B: pooled AUPR | Arm B: {{LEAD_NAME}} minus method (total 95% interval) |
|:----------|:------:|:-------------------|:------:|:------:|:-------------------|
| {{LEAD_NAME}} | 0.0669 | | | 0.2076 | |
| {{LEAD_NAME}}2 (exploratory) | 0.0615 | not computed | not computed | 0.2087 | not computed |
| {{ACLR_NAME}} | 0.0401 | +0.027 (+0.011 to +0.042) | 0.002 | 0.1510 | +0.057 (+0.039 to +0.074) |
| CLR | 0.0582 | +0.009 (−0.003 to +0.021) | 0.15 | 0.1211 | +0.087 (+0.058 to +0.115) |
| \|r\| | 0.0516 | +0.015 (+0.0002 to +0.030) | 0.070 | 0.0754 | +0.132 (+0.091 to +0.174) |
| GENIE3 | 0.0592 | +0.008 (+0.001 to +0.015)\* | 0.052 | 0.1092 | +0.098 (+0.071 to +0.126)\* |
| PIDC | 0.0660 | +0.001 (−0.007 to +0.009)\* | 0.41 | 0.1422 | +0.065 (+0.041 to +0.090)\* |
| GRNBoost2 | 0.0392 | +0.028 (+0.015 to +0.040)\* | < 0.0001 | 0.0780 | +0.130 (+0.097 to +0.162)\* |

: Table C2. Pooled AUPR and the difference from {{LEAD_NAME}} on each truth. Holm adjustment over six comparisons; Arm B has no p-values. {{LEAD_NAME}}2 (the walk depth set by the data) is outside the seven methods and the registered comparisons: its AUPR is shown, below {{LEAD_NAME}}'s on Arm A and slightly above it on Arm B, and its differences, intervals and p-values were not computed. \* Error bars for GENIE3, GRNBoost2 and PIDC comparisons include gene resampling only and are narrower than they should be. On HepG2, {{LEAD_NAME}} is not shown better than CLR, |r|, GENIE3 or PIDC after the Holm adjustment (unadjusted one-sided p = 0.075, 0.023, 0.013 and 0.41), and the "best of the seven" statement fails because CLR and PIDC do not reach 0.025 even unadjusted; the intervals for the differences from CLR and PIDC include 0. Method names are temporary.

# D. The exploratory yeast estimate for {{LEAD_NAME}}

How large is the pooled AUPR difference between {{LEAD_NAME}} and {{ACLR_NAME}} on protein complexes, in a dataset and truth on which {{LEAD_NAME}} was neither selected nor tuned? We ran an estimation study once, with the plan fixed before compute: no p-value and no confirmed or not-confirmed label. Expression is the DREAM5 network 4 yeast compendium (536 arrays, 5,950 genes); the truth is EBI Complex Portal yeast complexes (262 merged complexes, 13,438 positive gene pairs, about 14.5 effectively independent complexes). The reported interval combines a gene jackknife and an experiment-series jackknife (59 series).

| Contrast (AUPR difference) | Estimate | Total interval (gene + series) | Complex-level interval (sensitivity) |
|:----------------------------|:--------:|:------------------:|:------------------:|
| {{LEAD_NAME}} minus {{ACLR_NAME}} | +0.096 | +0.067 to +0.125 | −0.023 to +0.215 |
| {{LEAD_NAME}} minus CLR | +0.039 | +0.010 to +0.067 | −0.026 to +0.103 |
| {{LEAD_NAME}}2 minus {{ACLR_NAME}} (exploratory) | +0.093 | not computed | not computed |
| {{ACLR_NAME}} minus CLR | −0.058 | −0.098 to −0.018 | −0.122 to +0.007 |

: Table D1. Pooled enrichment over random: CLR 18.1, {{ACLR_NAME}} 12.6, {{LEAD_NAME}} 21.7, {{LEAD_NAME}}2 21.4 (base rate 0.0105). {{LEAD_NAME}}2 is {{LEAD_NAME}} with the walk depth set by the data (t2 = 21 here), exploratory and outside the preregistration; its row is a point estimate only. The complex-level interval treats complexes as the sampling unit (delete-half jackknife); it is a sensitivity line, and the total interval is the reported one. Method names are temporary.

- The total interval for {{LEAD_NAME}} minus {{ACLR_NAME}} lies above 0, above the smallest effect of interest (+0.005) and above the *E. coli* estimate (+0.027). With complexes as the sampling unit the interval includes 0.
- The gain is carried by large complexes: +0.137 [+0.106, +0.168] for complexes of 20 or more genes, against −0.004 [−0.011, +0.004] for 3 to 6 genes and −0.003 [−0.009, +0.002] for 7 to 19 genes. Without the ribosomal complexes it is +0.046 [+0.025, +0.066]. Per complex (125 complexes with at least 10 positive pairs), 60 improve and 65 worsen. So the gain is about how large complexes compare with the rest of the genome, not about ranking within a complex.
- This is an exploratory estimate in a second organism, not a confirmation, replication or validation. The yeast Complex Portal truth is used up for any later complex candidate. A code review was done after the run, and the reviewer saw the headline estimate by accident, so it is not fully blind; it found no bug that could change a reported number.

# E. The method in pseudocode

{{ACLR_NAME}} has no learned parameters. Given an expression compendium X (G genes by E experiments), the operator is:

```
# 1. dependence between every pair of genes
D[i,j] = MI(X[i], X[j])            # B-spline MI; single cells: |r| of Blom normal scores (CLR-NS)

# 2. CLR: z-score each gene's row against its own background
for each gene i:
    z[i,j] = max(0, (D[i,j] - mean(D[i,:])) / sd(D[i,:]))
S[i,j] = (z[i,j] + z[j,i]) / sqrt(2)           # pair score; S[i,i] = 0
# (Stouffer combination, as used here; the 2007 paper describes sqrt(z[i,j]^2 + z[j,i]^2))

# 3. attention: per-gene-temperature softmax over each gene's strongest partners
for each gene i:
    C = the (at most 50) j with largest S[i,j] and S[i,j] > 0
    if |C| <= 10:  A[i,j] = 1/|C| for j in C
    else:  find temperature T_i by bisection so that
              exp( entropy( softmax(S[i,C] / T_i) ) ) = 10      # ~10 effective partners
           A[i,j] = softmax(S[i,C] / T_i)[j]
    (a gene with no candidates keeps its own mass)

# 4. lazy multi-step walk
P = 0.5 * I + 0.5 * A               # rows sum to 1: half the mass stays, half attends
M = P^t, with t = 9                  # depth chosen without labels (below)

# 5. read-out
score[i,j] = 0.5 * (M[i,j] + M[j,i])   # attention mass from i to j, symmetrised

# depth t* (label-free)
split experiments into two halves by laboratory; for each depth t:
    agreement(t) = AUPR of half-1's M_t for the top 1% of half-2's S scores
t* = geometric mean of the two directions' best t; median over 50 splits = 9
```

Defaults are fixed in advance and never tuned on a truth set: 10 effective partners per gene, 50 candidates at most, lazy weight 0.5, depth 9. Step 2 with the Euclidean combination, sqrt(z[i,j]^2 + z[j,i]^2), is CLR as published in 2007; the Stouffer combination used here (and offered by the released 2007 software) ranks nearly identically; steps 3 to 5 are what this work adds. Transformer analogy: MI is the query-key score, the CLR z-score is the scaled logit, the softmax is the softmax, and the walk is repeated attention with a residual connection.

**{{LEAD_NAME}}** keeps steps 1 to 4 unchanged. In step 5, it replaces the raw walk mass M = P^t by its logarithm, which is the log-probability of reaching j from i in t steps (a multi-hop attention logit). It then z-scores each gene's row of the log matrix against its own background, using only pairs the walk reached, and combines the two directions as in step 2, clipped at 0. The pipeline is CLR, softmax, walk, log, CLR: one normalization layer on each side of the attention operator. No parameter is added and t = 9 is unchanged.

```
Pt = P^9;  L = log(Pt)  # unreachable pairs (Pt = 0) are set to missing
for each gene i:  Z[i,j] = max(0, (L[i,j] - mean_j L[i,:]) / sd_j L[i,:])   # over reachable j only
score[i,j] = (Z[i,j] + Z[j,i]) / sqrt(2)
```

# F. AI disclosure

The analysis code, runs and drafts were produced with Claude (Anthropic), with independent AI statistical reviews, under the author's direction. A contributions statement is drafted to match Nature Portfolio policy.

# G. References for the prior art

Faith, J. J. *et al.* Large-scale mapping and validation of *Escherichia coli* transcriptional regulation from a compendium of expression profiles. *PLoS Biol.* **5**, e8 (2007).

Nadaraya, E. A. On estimating regression. *Theory Probab. Appl.* **9**, 141–142 (1964).

Schlag, I., Irie, K. & Schmidhuber, J. Linear transformers are secretly fast weight programmers. In *Proc. 38th Int. Conf. Mach. Learn.*, PMLR **139**, 9355–9366 (2021).

Schmidhuber, J. Learning to control fast-weight memories: an alternative to dynamic recurrent networks. *Neural Comput.* **4**, 131–139 (1992).

Vaswani, A. *et al.* Attention is all you need. In *Adv. Neural Inf. Process. Syst.* **30** (2017).

Watson, G. S. Smooth regression analysis. *Sankhyā A* **26**, 359–372 (1964).

[^names]: {{ACLR_NAME}} is CLR with the attention summary added; method names are temporary.

[^hub]: Printed values, no tests; M3D and K562 were tuning sets.
