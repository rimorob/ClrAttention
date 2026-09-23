# CLR as attention on single-cell data:  BEELINE case study

*Working draft and findings log, 2026-09-23.  Results pending; the full run is queued on the local workstation (`results/beeline/`).*

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

## 4. Results

*To be filled from `results/beeline/summary.csv` and `tfnode_metrics.csv`:*

- median AUPRC ratio and EPR per method across 7 datasets × 2 gene sets × 3–4 truths;
- ratio to |Pearson|;
- the fraction of settings in which each method beats CLR;
- held-out depths per dataset.

## 5. Notes for interpretation

- **The BEELINE truths are TF-to-target networks, which is the "regulator as node" convention.**  In single-cell data a TF's mRNA is often a poor proxy for its activity.  The co-membership metric (targets sharing a TF) is the analogue of the regulator-agnostic *E. coli* benchmark, and may be the fairer test of attention.
- **BEELINE's inputs are not perturbation experiments.**  The mESC LOF/GOF network is a ground truth only.  Perturbation-target identification is therefore evaluated on M3D instead (*E. coli* white paper, section 6).
- **Baselines for comparison.**  Published BEELINE results for PIDC (MI-based), GENIE3, GRNBoost2 and PPCOR should be tabulated from the paper's supplementary tables rather than rerun.  That still needs doing, and needs the supplement.
