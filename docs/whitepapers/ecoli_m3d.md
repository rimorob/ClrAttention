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

## 4. 907-chip compendium and replicate-aware bootstrap (running)

Point estimates on all 907 chips reproduce the 466-experiment result.  Co-membership AUPR at strong/confirmed evidence:

| Method | AUPR |
|---|---|
| CLR | 0.135 |
| `soft10` at t = 8 | 0.159 |
| `pearson_top10` at t = 6 | 0.137 |

Each bootstrap draw recomputes MI from scratch.  There are two schemes:

- **Cluster:**  40 draws, resampling experiments with one replicate chip per occurrence.
- **Replicate:**  20 draws, one chip per experiment.

*Results to be added:*  `results/chips_bootstrap/summary.csv`.

## 5. Secondary:  TF-node edges (Faith 2007 convention)

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

**Early signal, 800-gene cloud test with 5 null draws (not a result):**  CLR-attention influence had a median AUROC of 0.63, against 0.55 for differential expression.  It won in 10 of 14 perturbations, with the largest gains for LexA/SOS (0.88 against 0.51), RyhB, SoxS and AppY.  It failed for ppGpp, which is again the diffuse regulon.  *The full run is queued:  `results/perturbation/`.*

## 7. Limitations

- **Circularity checks.**  RegulonDB evidence partly comes from expression studies, possibly including M3D experiments.  The strong/confirmed tier limits this, and it applies equally to every method compared.
- **Depth criterion.**  The held-out criterion predicts CLR structure in the other half of the data, so it favours depths that stabilize CLR structure.  That is a reasonable proxy, but not identical to biological accuracy.
- **Sample structure.**  Experiments within a lab are correlated.  The lab split handles this for depth selection, but the gene jackknife does not model it.  The cluster bootstrap on the 907 chips addresses it.
