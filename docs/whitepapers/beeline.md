# CLR as attention on single-cell data:  BEELINE case study

*Working draft and findings log, 2026-09-23.  Full run completed on the local workstation (`results/beeline/`).*

## 1. Question

Does the gain from CLR attention seen on bulk *E. coli* expression carry over to sparse, noisy single-cell data?  The test uses the standard single-cell GRN benchmark, where most published methods perform close to random (Pratapa et al. 2020).  A second question:  does calibrated, *estimated* attention do better than simple correlation and linear baselines?  Learned attention from single-cell foundation models reportedly does not (Kendiukhov 2026).

## 2. Data and protocol (BEELINE)

**Datasets.**  The seven experimental scRNA-seq datasets from BEELINE (Zenodo 3701939), in log-normalized counts:

| Dataset | Species | Cells |
|---|---|---|
| hESC | human | 758 |
| hHep | human | 425 |
| mDC | mouse | 383 |
| mESC | mouse | 421 |
| mHSC-E | mouse | 1,071 |
| mHSC-GM | mouse | 889 |
| mHSC-L | mouse | 847 |

**Gene sets.**  BEELINE's "TFs + 500" and "TFs + 1000" selections:  genes with a Bonferroni-corrected VGAM p below 0.01, keeping every TF among them plus the N most variable.

**Ground truths.**  Cell-type-specific ChIP-seq, non-specific ChIP-seq and STRING for every dataset.  mESC also has a LOF/GOF network.

**Metrics.**  These follow BEELINE:

- Candidate edges are (TF, gene) pairs.
- AUPRC ratio is AUPRC divided by the edge density.
- EPR is precision among the top k edges, where k is the number of true edges, divided by density.

A secondary metric is co-membership AUPR ratio, using TF target sets from the cell-type ChIP network as regulons.

**Methods.**  These are the same as in the E. coli study, fixed in advance:

- |Pearson| and |Spearman|.
- Ridge-regularized partial correlation, as the linear-model baseline.
- Raw B-spline MI.
- CLR (HG bins, Stouffer).
- CLR attention (`soft10`) and |Pearson| attention (`pearson_top10`), each at depths chosen on held-out random halves of cells.

## 3. Early signal (mESC, TFs + 500, one cloud test; not a result)

- 1,071 genes (620 TFs) × 421 cells.
- Held-out depth:  t* = 7 for CLR attention (the two halves gave 7 and 8) and t* = 2 for Pearson attention.
- Against the cell-type ChIP network, every method is close to random, which matches BEELINE's own findings for this truth.

| Method | AUPRC ratio | EPR |
|---|---|---|
| CLR attention | 1.107 | 1.133 |
| \|Pearson\| | 1.088 | 1.114 |
| CLR | 1.075 | 1.094 |
| Partial correlation | 0.973 | 0.965 |

- Across all four mESC ground truths, plain CLR had the highest median AUPRC ratio (1.46; EPR 2.38).  CLR attention was next (1.41), then |Pearson| (1.33).
- Partial correlation, the linear-model baseline, was the weakest.

## 4. Results (7 datasets × 2 gene sets × 3–4 ground truths = 44 settings)

**Held-out depth.**  The held-out cell split chose t* = 6–11 for CLR attention in every dataset.  That matches the *E. coli* depth of 8, although the depth was never tuned on single-cell data.  Pearson attention always chose t* = 2.

**Median AUPRC ratio over all 44 settings:**

| Method | Median AUPRC ratio | Median EPR | Beats CLR (AUPRC) |
|---|---|---|---|
| CLR | 1.48 | 2.84 | — |
| CLR attention (`soft10`) | 1.44 | 2.29 | 27% |
| \|Pearson\| attention | 1.23 | 2.08 | 20% |
| \|Spearman\| | 1.21 | 1.67 | 25% |
| \|Pearson\| | 1.20 | 1.60 | 36% |
| Raw MI | 1.14 | 1.69 | 20% |
| Partial correlation (ridge) | 1.06 | 1.22 | 23% |

**By ground truth (median AUPRC ratio):**

| Truth | CLR | CLR attention | \|Pearson\| | Partial correlation |
|---|---|---|---|---|
| Cell-type ChIP-seq (14) | 1.00 | 0.93 | 1.02 | 1.01 |
| Non-specific ChIP-seq (14) | 1.54 | 1.59 | 1.32 | 1.08 |
| STRING (14) | 2.11 | 1.76 | 1.90 | 1.19 |
| mESC LOF/GOF (2) | 1.45 | 1.54 | 1.34 | 0.93 |

**Co-membership** (targets sharing a TF in the cell-type ChIP network):  every method scores 0.99–1.00, which is random.

**Topline.**

1. CLR is the best method overall on BEELINE.  It beats correlation, MI and the linear-model baseline by about 20–40% in AUPRC ratio and more in EPR.  The MI-plus-context-normalization idea carries over to single-cell data.
2. CLR attention does *not* improve on CLR here.  It is about equal or slightly better on non-specific ChIP and LOF/GOF, and worse on STRING.  The bulk *E. coli* gain does not transfer.
3. Against cell-type-specific ChIP-seq, and on regulon co-membership, every method is at chance.  This agrees with BEELINE's own finding that these truths are barely recoverable from scRNA-seq co-expression.  Nothing in this benchmark can separate the methods on those truths.
4. The ridge partial correlation (the linear-model baseline) is the weakest method, which fits the pattern that linear conditional-independence models struggle on sparse single-cell data.

**Possible reasons attention helps in bulk data but not here.**  These are hypotheses and have not been tested.

- Single-cell co-expression is dominated by broad cell-state programs and by dropout.  Multi-hop propagation spreads program-level signal, which may blur the specific edges that STRING rewards.
- BEELINE scores directed TF-to-target edges.  Attention mass is a symmetric multi-hop quantity, better matched to co-membership, and that metric is at chance for every method here.

## 5. Notes for interpretation

- **The BEELINE truths are TF-to-target networks, which is the "regulator as node" convention.**  In single-cell data a TF's mRNA is often a poor proxy for its activity.  The co-membership metric (targets sharing a TF) is the analogue of the regulator-agnostic *E. coli* benchmark, and may be the fairer test of attention.
- **BEELINE's inputs are not perturbation experiments.**  The mESC LOF/GOF network is a ground truth only.  Perturbation-target identification is therefore evaluated on M3D instead (*E. coli* white paper, section 6).
- **Baselines for comparison.**  Published BEELINE results for PIDC (MI-based), GENIE3, GRNBoost2 and PPCOR should be tabulated from the paper's supplementary tables rather than rerun.  That still needs doing, and needs the supplement.
