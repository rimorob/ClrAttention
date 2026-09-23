# CLR as attention in *E. coli*:  M3D × RegulonDB

*Working draft and findings log, 2026-09-23.*

## 1. Question

Does completing CLR as an attention operator improve recovery of *E. coli* regulons from a large expression compendium?  "Completing" means softmax weighting plus depth chosen on held-out data.  Does the gain depend on the calibrated MI kernel, rather than on propagation alone?

## 2. Data

- **Expression:**  the M3D *E. coli* compendium, v4 build 6 (Faith et al. 2008).  It has 4,297 genes, measured on 907 Affymetrix chips from 466 experiments.  The primary analysis uses the 466 replicate-averaged experiments, which avoids pseudo-replication.  The 907-chip set is analysed with a replicate-aware bootstrap.
- **Gold standard:**  RegulonDB release 14.5 (RISet, TUSet, OperonSet, NetworkSigmaGene).  Every regulator class is included:
  - 155 regulons at strong/confirmed evidence:  126 TF, 18 sRNA, 6 sigma, 5 other protein.
  - 196 regulons at any evidence level, which adds one ppGpp regulon.
  - 1,879 annotated genes and 1.76 M scored pairs at strong/confirmed evidence, with a base rate of 9.4%.
  - σ70, CRP and Nac exceed the 500-gene size cap and are excluded.
- **Name mapping:**  2008-era names are mapped to b-numbers through NCBI synonyms.  2,946 of 3,095 regulated genes are on the array.

## 3. Topline findings

All numbers are on the 466 experiments with strong/confirmed regulons, unless stated otherwise.  Intervals are 95%.

1. **CLR beats correlation and raw MI.**  Co-membership AUPR:
   - CLR (HG bins, Stouffer):  0.136.
   - Raw MI:  0.115, a difference of −0.021 [−0.026, −0.016].
   - |Pearson|:  0.116, a difference of −0.020 [−0.024, −0.016].
   - The 2007 parity settings:  0.134, significantly but slightly below the current default (−0.002 [−0.003, −0.001]).
2. **Bias-controlled binning matters a little, in the expected direction.**  Coarse common bins (6–10) beat per-gene Freedman–Diaconis bins (median 33), and CLR AUPR falls monotonically as bins increase from 8 to 16.
3. **CLR attention at the held-out depth improves on single-step CLR.**
   - **Depth selection:**  the two lab-split directions agreed, choosing t* = 10 and 7, a geometric mean of 8.
   - **Pooled co-membership AUPR:**  `soft10` scored 0.157, **+0.021 [+0.015, +0.027]** over CLR, a 15% relative gain.
   - **Per-regulon coherence:**  median AUROC rose from 0.597 to 0.632, **+0.031 [+0.018, +0.046]**.  69% of 144 regulons improved (Wilcoxon p = 4×10⁻⁴).
   - **All-evidence benchmark:**  +0.007 AUPR and +0.026 coherence, both significant.
4. **The gain requires the CLR kernel.**  The identical construction on |Pearson| (`pearson_top10` at its held-out t* = 6) does not beat CLR on AUPR (−0.001 [−0.005, +0.003]).  It is significantly *worse* on coherence (−0.014 [−0.025, −0.003]).
5. **The softmax matters, not just density.**  The hard-thresholded CLR operator (`fdr05_top10`) gains less:  +0.013 AUPR, and no coherence gain.
6. **Depth has an optimum.**  Pooled AUPR keeps rising until about t = 32, but regulon coherence peaks at t = 4–8 and collapses by t = 128.  Deep propagation blurs large, diffuse programs at the expense of specific regulons.  The held-out criterion, which never sees RegulonDB, selects the depth where both metrics improve.
7. **Regulator classes.**  Most of the gain comes from TF regulons:  AUPR 0.114 → 0.145 at t = 20, and coherence 0.59 → 0.63.  Sigma regulons gain a little.  sRNA regulons are nearly unrecoverable from mRNA co-expression under any method.  The diffuse ppGpp regulon *loses* under propagation.
8. **Negative results worth reporting:**
   - Diffusing profiles and then correlating them adds little.  The best case is +0.007 at t = 32 with the `fdr20` operator.
   - Signed and conditional-expectation value transforms do not help, because 99% of selected edges are positively correlated.
   - Stouffer and Euclidean combination are equivalent.
   - The permutation-FDR operator alone leaves 62–77% of genes without neighbours, which makes it unusable for propagation.

## 4. 907-chip compendium and replicate-aware bootstrap

Every bootstrap draw recomputes MI, CLR and attention from scratch.  There are two schemes:

- **Cluster:**  40 draws.  The 466 experiments are resampled with replacement, with one replicate chip per occurrence.  This captures condition-level and technical variability, and it gives the headline intervals.
- **Replicate:**  20 draws.  One randomly chosen chip is used per experiment, which captures technical variability only.

Co-membership AUPR and per-regulon coherence at strong/confirmed evidence.  Differences are paired against CLR within each draw, with 2.5–97.5% bootstrap percentiles:

| Method | Point AUPR (907 chips) | ΔAUPR vs CLR, cluster | ΔAUPR vs CLR, replicate | ΔCoherence vs CLR, cluster |
|---|---|---|---|---|
| CLR (HG, Stouffer) | 0.135 | 0 | 0 | 0 |
| `soft10` at t = 8 | 0.158 | **+0.020 [+0.018, +0.022]** | +0.021 [+0.020, +0.022] | **+0.030 [+0.014, +0.051]** |
| `pearson_top10` at t = 6 | 0.136 | −0.002 [−0.004, +0.002] | −0.002 [−0.003, −0.001] | −0.026 [−0.049, −0.011] |

With all evidence levels, `soft10` gains +0.007 [+0.005, +0.010] AUPR and +0.014 [+0.004, +0.026] coherence in the cluster bootstrap.  `pearson_top10` loses on both (−0.013 AUPR and −0.034 coherence).

**Reading.**  The 907-chip analysis reproduces the 466-experiment result almost exactly, and resampling whole experiments does not shrink the gain.  The attention gain is therefore not an artefact of pseudo-replication or of a few influential experiments.  Propagation on |Pearson| gives no gain on AUPR and costs coherence, so the gain again depends on the CLR kernel.

*An extension to 200 cluster and 100 replicate draws is queued, to firm up the tail percentiles.*

## 5. Secondary:  TF-node edges (Faith, Hayete et al. 2007 convention)

*To be summarized from `tfnode_edge_pr_by_config.csv`.*  This convention uses the regulator's own mRNA as a stand-in for its activity.  It is reported for continuity with the 2007 paper.

## 6. Perturbation-target identification (running)

**Question:**  can the network identify a perturbed regulator's targets better than differential expression?  Kendiukhov (2026) found that attention maps from single-cell foundation models do not beat such trivial features.

**Design:**

- **Perturbations:**  14 M3D perturbation groups with a RegulonDB-mapped regulon.  The gene-to-regulator mapping was fixed in advance:
  - rpoS → σ38
  - relA → ppGpp
  - recA → LexA/SOS
  - hupA/B → HU
  - fnr, arcA, fis, crp, soxS, oxyR, appY, cpxR, hns and ryhB each map to their own regulon.
- **Score:**  the delete-k jackknife influence of P's experiments on |Pearson|, CLR and CLR-attention networks.  It is z-scored per gene against random delete-k subsets.
- **Baselines:**  differential expression (mean |z| of P's experiments) and gene variance.

**Reference method.**  The then state of the art on this exact task and compendium is SSEM-Lasso network filtering (Cosgrove, Zhou, Gardner & Kolaczyk, *Bioinformatics* 2008, [doi:10.1093/bioinformatics/btn476](https://doi.org/10.1093/bioinformatics/btn476)).  It fits a sparse regression network and scores genes by how much the perturbation departs from what the network predicts.  It was evaluated on M3D genetic perturbations by sensitivity among the top 100 ranked genes, a metric we now report as well.  Boris's original KL-divergence approach matched it but did not beat it.  Re-running SSEM-Lasso, and porting the KLD method once the MATLAB source is located, are the natural comparisons.

**First full run (all 4,297 genes, 30 null draws per k, 15 perturbation groups).**  Median AUROC for recovering the perturbed regulator's regulon:

| Score | Median AUROC | Wins vs DE | Paired Wilcoxon p vs DE |
|---|---|---|---|
| Gene variance (ignores the perturbation) | 0.646 | 14/15 | 2×10⁻⁴ |
| CLR-attention influence | 0.590 | 13/15 | 0.013 |
| CLR influence | 0.573 | 13/15 | 0.008 |
| DE + CLR-attention (Stouffer) | 0.564 | 12/15 | 0.015 |
| \|Pearson\| influence | 0.526 | 9/15 | 0.85 |
| Differential expression | 0.519 | — | — |

**This run exposes a confound, and it is not yet a positive result.**  Gene variance knows nothing about the perturbation, yet it beats every score.  Regulon members are simply more variable than other genes, so any score correlated with variance gets credit.  The network influence scores beat differential expression, but they may partly be tracking variance too.  The perturbation-specific cases are encouraging:  CLR-attention influence exceeds gene variance for LexA/SOS (0.73 against 0.55), HU (0.83 against 0.75) and AppY (0.90 against 0.82).  Across the board, however, the variance prior wins.  Sensitivity among the top 100 genes is low for every method (at most 0.07).

**Variance-controlled run (100 null draws per k).**  The stratified AUROC compares targets only with non-targets in the same gene-variance decile, so a score that merely tracks variance gets 0.5.

| Score | Median AUROC, raw | Median AUROC, variance-stratified | Wins vs DE (stratified) | Paired p vs DE (stratified) |
|---|---|---|---|---|
| CLR-attention influence | 0.572 | **0.573** | 13/15 | 0.010 |
| CLR influence | 0.569 | 0.553 | 12/15 | 0.073 |
| DE + CLR-attention (Stouffer) | 0.552 | 0.539 | 12/15 | 0.026 |
| \|Pearson\| influence | 0.543 | 0.537 | 8/15 | 0.56 |
| Gene variance | 0.646 | 0.523 | 8/15 | 1.0 |
| Differential expression | 0.519 | 0.516 | — | — |

**Verdict.**  The variance prior's advantage was entirely the confound:  once variance is held fixed it drops to 0.52.  CLR-attention influence keeps its advantage (0.573), and it beats differential expression in 13 of 15 perturbations (p = 0.010).  The network influence therefore carries perturbation-specific information beyond gene variance.  The effect is modest, and it is concentrated in a few regulators:

- HU:  0.82 against 0.56 for DE.
- AppY:  0.77 against 0.19.
- LexA/SOS:  0.68 against 0.59.
- RyhB:  0.76 against 0.73.
- ArcA:  0.63 against 0.53.

It fails for ppGpp (0.42), CRP (0.43) and SoxS (0.46).  ppGpp and CRP are diffuse, global regulators, and the diffuse ppGpp regulon also resists attention in the co-membership benchmark.  CLR-attention influence beats CLR influence in 10 of 15 perturbations (median difference +0.015), but that difference is not significant (p = 0.13).  |Pearson| influence adds nothing over DE, so, as in the regulon benchmark, the signal needs the CLR kernel.  Sensitivity among the top 100 genes stays low for every method (at most 0.07).  That is well below the SSEM-Lasso sensitivities reported by Cosgrove et al. (2008), although their gold standard and gene universe differ.  A direct re-run of SSEM-Lasso is needed before comparing.

## 7. Limitations

- **Circularity checks.**  RegulonDB evidence partly comes from expression studies, possibly including M3D experiments.  The strong/confirmed tier limits this, and it applies equally to every method compared.
- **Depth criterion.**  The held-out criterion predicts CLR structure in the other half of the data, so it favours depths that stabilize CLR structure.  That is a reasonable proxy, but not identical to biological accuracy.
- **Sample structure.**  Experiments within a lab are correlated.  The lab split handles this for depth selection, but the gene jackknife does not model it.  The cluster bootstrap on the 907 chips addresses it.
