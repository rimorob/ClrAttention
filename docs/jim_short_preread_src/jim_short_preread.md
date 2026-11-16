---
title: "CLR with attention: a short pre-read for Jim Collins"
---

CLR (Faith et al. 2007) standardizes each gene's mutual information against that gene's own background. **{{ACLR_NAME}}** is CLR with attention: CLR, then a softmax over each gene's top 50 positive CLR scores, then a short lazy walk, with no learned parameters. **{{LEAD_NAME}}** is {{LEAD_DESC}}; it also has no fitted parameters and is our lead candidate.[^names] On *E. coli* regulons from RegulonDB, {{ACLR_NAME}} scores 1.68 times random against 1.44 for CLR and 1.51 for GENIE3, and {{LEAD_NAME}} scores 1.72; on K562 protein complexes {{LEAD_NAME}} scores 7.95 against 6.39 for CLR and 3.51 for {{ACLR_NAME}}. Methods, controls and further results are in the appendix.

![**Figure 1. Graphical abstract.** Left, CLR as published (Faith et al. 2007); the Euclidean combination is the one in the 2007 packages we found, the Stouffer combination is the one used here. Middle, the additions outlined in dotted green: each gene attends over its own top 50 positive CLR scores, a per-gene temperature sets about 10 effective partners, and each walk step follows the rows of the (non-symmetric) attention matrix; {{LEAD_NAME}} then adds the log of the walk and a second CLR. Right, pooled precision-recall area as a multiple of random. Dotted green = new in this work, including all human data. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p1.png){width=8.7in}

{{PAGEBREAK}}

![**Figure 2. Scoreboard, raw scores.** Every metric we have for CLR, {{ACLR_NAME}}, {{LEAD_NAME}} and the competitors. Bars start at random and are scaled within each row; the best value in each row is marked; no intervals or tests are shown; the bottom row counts, for each method, the rows (of 19) where it has the highest value or ties it at the printed precision, and a method that was not run on a row cannot win it; a dash is a cell that was not run. Dotted green outlines mark the human-data groups. {{LEAD_NAME}} was picked among about 20 read-outs by its K562 scores, and K562 and M3D are tuning sets for it. The hub-controlled version of the same rows is in Appendix F. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p3.png){width=8.8in}

**Evidence for {{LEAD_NAME}}.** In a preregistered test on *E. coli*, {{LEAD_NAME}}'s gain over {{ACLR_NAME}} on protein complexes was +0.026 AUPR on M3D (p = 0.079) and +0.029 on EcoMAC (p = 0.051), which did not reach the preregistered threshold (one-sided p < 0.025; Appendix C). An exploratory estimate on yeast protein complexes gives +0.096 AUPR over {{ACLR_NAME}} (total 95% interval +0.067 to +0.125); the complex-level interval, −0.023 to +0.215, includes 0, and the gain sits in complexes of 20 or more genes (Appendix D).

{{PAGEBREAK}}

**Losses and caveats.** {{LEAD_NAME}} was selected among about 20 read-outs on K562 truths, so its K562 numbers are a selection and not a fit (it has no fitted parameters, but the best of many noisy scores is biased upward). CLR and GENIE3 are ahead of both methods on direct TF-to-target edges in *E. coli* (AUPR ratio: CLR 6.75, GENIE3 6.54, {{LEAD_NAME}} 5.84, {{ACLR_NAME}} 5.17). {{ACLR_NAME}} is below CLR on K562 protein complexes (3.51 against 6.39; 86% of the loss comes from ten ribosomal and pre-rRNA complexes), and per regulon it ties CLR and GENIE3. Every positive result for {{ACLR_NAME}} comes from the *E. coli* development lineage (M3D with RegulonDB, and the overlapping EcoMAC), where its settings were chosen; in a preregistered DREAM5 yeast test its co-membership gain over CLR (5.8%) stayed below the 10% margin. All methods are within 10% of random on K562 transcription-factor truths. Appendix A gives the remaining caveats.

# Potential follow-ups and decision points

- **Perturbation targets.** Attention-based influence ranks a perturbed transcription factor's targets above differential expression in 12 of 15 perturbations (p = 0.015); under a laboratory-matched null the figure is 9 of 15 (p = 0.095). The result is suggestive and not established.
- **Differential CLR.** If perturbations are pursued, differential CLR is a possible extension. It is planned and not yet ported, and there is no result.
- **Decision points.** (1) Which method is presented as the headline: {{ACLR_NAME}}, {{LEAD_NAME}}, or both with their trade-offs. (2) Whether a further evaluation of {{LEAD_NAME}} on a new complex truth is warranted. (3) The naming of the methods.

{{PAGEBREAK}}

# Appendix

# A. Methods compared, and further results for {{ACLR_NAME}}

![**Figure A1. The methods compared.** *E. coli* AUPR is the component ablation (strong-evidence regulons, 466 experiments, RegulonDB; each method at its own held-out depth t\*, except {{LEAD_NAME}} at t = 9). K562 CORUM is pooled AUPR as a multiple of random. BEELINE values are median AUPRC ratios over the 44 case-truth cells of the primary truths (14 dataset-and-gene-set cases); for {{ACLR_NAME}} and {{LEAD_NAME}} they use the CLR-NS operator (CLR on the absolute correlation of normal scores; {{LEAD_NAME}} at t = 9, not tuned, exploratory). The two ablation rows show that the softmax without CLR normalization gains less. {{LEAD_NAME}} was picked among about 20 read-outs by its K562 scores; its preregistered *E. coli* test did not reach the preregistered threshold. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p2.png){width=9.2in}

**A.1 {{ACLR_NAME}} against established methods.**

| Data | Truth | {{ACLR_NAME}} | CLR | GENIE3 | PIDC | GRNBoost2 |
|:-----------------------|:------------------------|:------:|:------:|:------:|:------:|:------:|
| M3D | RegulonDB, strong evidence | **1.68** | 1.44 | 1.51 | 1.45 | 1.36 |
| EcoMAC (2,198 arrays; overlaps M3D) | strong evidence | **1.68** | 1.39 | 1.50 | 1.43 | 1.36 |
| M3D, hub-controlled | strong evidence | **1.56** | 1.40 | 1.47 | 1.40 | 1.37 |
| Human K562 Perturb-seq | DoRothEA A-B | **1.09** | **1.06**† | **1.08** | **1.07** | 1.05 |
| Human K562 Perturb-seq | CORUM complexes | 3.51 | **6.39** | 6.32 | **6.33** | 3.66 |

: Table A1. Pooled AUPR divided by the base rate (1.00 = random). Bold: best in the full row or statistically tied with it (paired gene jackknife, Bonferroni-corrected); † tied only after the correction. On CORUM, |r| is best (7.00). EcoMAC is not an independent replication of M3D: 828 of its 2,198 arrays (37.7%) have a one-to-one M3D counterpart, covering 91% of M3D's chips. On the 1,370 EcoMAC arrays without an M3D counterpart, {{ACLR_NAME}}'s lead over CLR is unchanged in size (+0.0285 AUPR, against +0.0278 on all arrays; a post hoc sensitivity analysis). Intervals resample genes only; in the DREAM5 yeast test the experiment-resampling part of the standard error was as large as the gene part or larger, so these intervals understate the total uncertainty and some differences marked as wins may be ties. GRNBoost2 uses legacy early stopping. Method names are temporary.

**A.2 *E. coli* M3D, 466 experiments, RegulonDB regulon co-membership (the headline).**

- {{ACLR_NAME}} beats CLR by +0.022 AUPR [+0.019, +0.026]. Depth is re-selected without labels in each of 1,000 nested bootstrap draws, and {{ACLR_NAME}} is ahead in all 1,000.
- The gain is consistent across 50 laboratory splits (+0.022 [+0.017, +0.027]) and a replicate-aware bootstrap on the 907 chips (+0.020).
- Per regulon (each regulon weighted equally), {{ACLR_NAME}} ties CLR and GENIE3; its win rate is highest in regulons of 25 or more genes.
- Context normalization is essential: a softmax on uncalibrated |r| or MI gets about a third of the gain and none of the coherence gain. A Shapley decomposition over six components attributes the AUPR gain mainly to depth (+0.018) and CLR normalization (+0.013).
- Removing 0 to 20 global principal components leaves the gain intact; the top components are biology (growth rate, motility, the stationary-phase program).
- Hub control: the walk converges toward well-connected genes, so every method is also scored within strata of hubness. {{ACLR_NAME}} stays first on strong-evidence truths (ahead of CLR in every hub decile on M3D) and ties GENIE3 on all-evidence truths. Under a stratifier that does not use the walk (row mean of raw MI; exploratory), it stays ahead of CLR on M3D strong evidence (+0.18 enrichment [+0.12, +0.24]); that covers {{ACLR_NAME}} against CLR only, not GENIE3 or PIDC.

**A.3 Further caveats.**

- **Human data.** On K562 transcription-factor truths every method is within 10% of random. On single cells (BEELINE), CLR on normal scores (CLR-NS) has the highest per-dataset mean of the methods tested, but it was chosen on these truths and is not significantly ahead of GENIE3 or PIDC, and GENIE3 has the higher median; {{ACLR_NAME}} adds nothing there. Against zero-shot attention from the single-cell foundation model scPRINT (preregistered, run once), CLR-NS is non-inferior and superior on all seven BEELINE datasets (Hodges-Lehmann estimate of the log2 AUPRC-ratio difference +1.03, exact interval +0.34 to +1.90), but scPRINT's attention scores the same as the same model run on permuted expression, so the comparator is weak here.
- **Complexes.** GRNBoost2 does about as poorly as {{ACLR_NAME}} on K562 complexes, so the loss is not specific to the walk. {{LEAD_NAME}} recovers it on K562 (7.95) but is below CLR on RPE1 merged CORUM (6.14 against 6.77) and below |r| on *E. coli* complexes (26.1 against 29.2). Per regulon, {{LEAD_NAME}} (25.2) is level with {{ACLR_NAME}} (25.6) and above CLR (23.0).
- **Direct edges.** {{ACLR_NAME}} also trails CLR on DREAM5 network 3, which shares the *E. coli* lineage.
- **Yeast, DREAM5 network 4 (preregistered, run once).** For {{ACLR_NAME}} against CLR, the co-membership gain was 5.8% (interval +1.8% to +10.0% with group-level resampling; it includes zero with series-level resampling), below the pre-specified 10% margin: the verdict "equivalent within 10%" sits on the boundary with "inconclusive", and the result is not "superior". The edge read-out was inconclusive. Every method is near chance on yeast (AUROC 0.51 to 0.52 for all seven methods).
- **Where the evidence comes from.** The settings (10 effective partners, 50 candidates, lazy weight 0.5, the MI estimator, the truth construction) were settled on the *E. coli* data after many looks; the nested bootstrap protects only the depth choice. Only the yeast tests lie outside that lineage.

{{PAGEBREAK}}

# B. Walk depth

![**Figure A2. AUPR gain over CLR at each depth, one panel per truth.** Walk depth *t* on a log2 axis from 1 to 256, with the best depth marked, the within-1% interval shaded and the grid edge 256 dashed. On K562 CORUM complexes, CLR is best and every attention depth is below it.](/home/boris/Projects/ClrAttention/results/granularity_log/depth_log_curves.png){width=8.4in}

The held-out depth is stable: median *t* = 9, with the criterion flat within 1% for *t* = 7 to 11.

| Organism | Truth, fine to coarse | Median size (sets) | Best t | AUPR gain over CLR |
|:-----------|:-------------------|:----------:|:------------:|:------------:|
| *E. coli* | 1 operons | 2 (786) | 4 | +0.208 |
| *E. coli* | 2 TF regulons | 11 (158) | 96 | +0.027 |
| *E. coli* | 3 sigmulons | 128 (5) | 12 | +0.043 |
| K562 | 1 CORUM complexes | 4 (1,186) | CLR (0) | 0.000 |
| K562 | 2 KD-response modules | 100 (585) | 24 | +0.017 |
| K562 | 3 DoRothEA A-B | 13 (105) | 4 | +0.008 |
| K562 | 4 DoRothEA A-C | 28 (269) | 8 | +0.005 |
| K562 | 5 CollecTRI | 16 (490) | 48 | +0.014 |

: Table B1. Best t is the highest-AUPR point over CLR and the geometric grid t = 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, plus CLR as t = 0. AUPR gain: absolute AUPR at the best depth minus CLR's. K562 CORUM: CLR beats all 16 attention depths; the best attention depth is t = 16, with an AUPR 0.037 below CLR's (0.0894). The pre-specified test of whether best depth rises from fine to coarse truths (Spearman rho at most 0.5 in at least one organism) came out as "not resolution-like": *E. coli* rho = 0.500 (n = 3), K562 rho = 0.700 (n = 5). With 3 and 5 truths and flat depth curves, the test could neither establish nor exclude the ordering.

# C. The preregistered *E. coli* test of {{LEAD_NAME}}

- Complex superiority over {{ACLR_NAME}} at *t* = 9 was required at one-sided p < 0.025 on both *E. coli* datasets. It was +0.026 AUPR on M3D (p = 0.079) and +0.029 on EcoMAC (p = 0.051), so the threshold was not reached. The jackknife standard error on complexes was 0.018, and the smallest effect the design could detect was about 0.05.
- Regulon non-inferiority (margin 0.005 AUPR) passed on both: +0.004 and +0.003.
- {{LEAD_NAME}} has no fitted parameters, so it cannot memorise the truth, but it was picked as the best of about 20 read-outs by its K562 scores, and the best of many noisy scores is biased upward. The candidates were highly correlated and the gain is large, so the bias is probably a fraction of the gap, but it has not been measured. On RPE1, another cell line (descriptive only), {{LEAD_NAME}} is +0.036 [+0.013, +0.059] AUPR above {{ACLR_NAME}} on merged CORUM.

# D. The exploratory yeast estimate for {{LEAD_NAME}}

How large is the pooled AUPR difference between {{LEAD_NAME}} and {{ACLR_NAME}} on protein complexes, in a dataset and truth on which {{LEAD_NAME}} was neither selected nor tuned? We ran an estimation study once, with the plan fixed before compute: no p-value and no confirmed or not-confirmed label. Expression is the DREAM5 network 4 yeast compendium (536 arrays, 5,950 genes); the truth is EBI Complex Portal yeast complexes (262 merged complexes, 13,438 positive gene pairs, about 14.5 effectively independent complexes). The reported interval combines a gene jackknife and an experiment-series jackknife (59 series).

| Contrast (AUPR difference) | Estimate | Total interval (gene + series) | Complex-level interval (sensitivity) |
|:----------------------------|:--------:|:------------------:|:------------------:|
| {{LEAD_NAME}} minus {{ACLR_NAME}} | +0.096 | +0.067 to +0.125 | −0.023 to +0.215 |
| {{LEAD_NAME}} minus CLR | +0.039 | +0.010 to +0.067 | −0.026 to +0.103 |
| {{ACLR_NAME}} minus CLR | −0.058 | −0.098 to −0.018 | −0.122 to +0.007 |

: Table D1. Pooled enrichment over random: CLR 18.1, {{ACLR_NAME}} 12.6, {{LEAD_NAME}} 21.7 (base rate 0.0105). The complex-level interval treats complexes as the sampling unit (delete-half jackknife); it is a sensitivity line, and the total interval is the reported one. Method names are temporary.

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

Defaults are fixed in advance and never tuned on a truth set: 10 effective partners per gene, 50 candidates at most, lazy weight 0.5, depth 9. Step 2 alone is CLR; steps 3 to 5 are what this work adds. Transformer analogy: MI is the query-key score, the CLR z-score is the scaled logit, the softmax is the softmax, and the walk is repeated attention with a residual connection.

**{{LEAD_NAME}}** keeps steps 1 to 4 unchanged. In step 5, it replaces the raw walk mass M = P^t by its logarithm, which is the log-probability of reaching j from i in t steps (a multi-hop attention logit). It then z-scores each gene's row of the log matrix against its own background, using only pairs the walk reached, and combines the two directions as in step 2, clipped at 0. The pipeline is CLR, softmax, walk, log, CLR: one normalization layer on each side of the attention operator. No parameter is added and t = 9 is unchanged.

```
Pt = P^9;  L = log(Pt)  # unreachable pairs (Pt = 0) are set to missing
for each gene i:  Z[i,j] = max(0, (L[i,j] - mean_j L[i,:]) / sd_j L[i,:])   # over reachable j only
score[i,j] = (Z[i,j] + Z[j,i]) / sqrt(2)
```

# F. Hub-controlled scoreboard

Hub control: the walk converges toward well-connected genes (hubs), so each method's scores are converted to percentiles within 10 strata of the attention walk's stationary distribution (the hub score of a gene pair, pi_i + pi_j), and methods are compared within strata. A method is therefore not rewarded merely for favouring hubs. Bars start at random and are scaled within each row, so compare the printed values across a row, not bar lengths across rows. The best value in each row is marked with a dot; ties and intervals are not shown.

![**Figure A3. Hub-controlled scoreboard.** The same rows as Figure 2 that have a hub-controlled version (not run for direct edges, per-regulon AP and BEELINE; the yeast row is left out because its hub-controlled scores were stored with a different stratifier). Dotted green outlines mark the human-data groups. Method names are temporary.](/home/boris/Projects/ClrAttention/docs/whitepapers/figures/graphical_abstract_short_p3b.png){width=8.6in}

# G. AI disclosure

The analysis code, runs and drafts were produced with Claude (Anthropic), with independent AI statistical reviews, under the author's direction. A contributions statement is drafted to match Nature Portfolio policy.

[^names]: Method names are temporary.
