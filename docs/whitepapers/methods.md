# CLR as attention:  methods

*Working draft and findings log, 2026-09-23.  See `CITATION_LOG.md` (D17–D28) for the full decision record.*

## 1. Summary

CLR (context likelihood of relatedness; Faith, Hayete et al. 2007) turns pairwise mutual information between genes into context-normalized scores.  For each gene, every other gene's MI is expressed as a z-score against that gene's own MI background, and the two directions are combined.  Read as a machine-learning operator, this is an attention score matrix:  each gene (query) receives a calibrated, data-dependent weight over every other gene (key), normalized within its own context.  CLR stopped at the scores.  It never used the rest of the attention machinery:  weight normalization with a temperature, aggregation, and depth.

This project completes the operator.  The completed method has three stages:

1. Estimate MI with a bias-controlled B-spline estimator.
2. Calibrate it into CLR scores with a permutation null.
3. Turn the scores into softmax attention and apply it for a depth chosen on held-out data.

The historical claim we can defend is narrower than "CLR was the first attention mechanism".  Kernel smoothing (Nadaraya–Watson 1964) and fast-weight programmers (Schmidhuber 1991) are earlier instances of the pattern.  The defensible claim is that CLR independently computed attention-style scores over continuous molecular measurements, seven years before neural attention.  Adding the components that attention research made standard gives a significant, controlled improvement.  See `docs/paper_framing.md`.

## 2. Mutual-information estimation

MI is estimated with the B-spline estimator of Daub et al. (2004), spline order 3, ported from the original C sources and verified to machine precision (≤ 3×10⁻¹⁵) against an independent implementation.  The production kernel exploits the fact that each sample has at most three non-zero spline weights.  It accumulates the joint histogram sparsely and is bit-identical to the dense kernel at default compiler flags, but 5–20× faster.

**Bin count.**  Every gene uses one common bin count, chosen for the joint (2-D) histogram rather than for either marginal.  The default is the Hacine-Gharbi et al. (2012) low-bias rule for histogram MI, evaluated at ρ = 0, which gives 9 bins at N = 907 and 7 at N = 466.  Two findings motivate this:

- Per-gene adaptive bins bias MI upward in proportion to bins_i × bins_j.  Freedman–Diaconis gives heavy-tailed genes up to 50 bins, so they became false hubs, and CLR calibration did not remove the bias.
- Bins sized for the marginal leave the joint table under-populated.  At 33 bins per axis there are 1,089 cells for 907 points.

On M3D, the common coarse-bin rules (6–10 bins) are indistinguishable from one another, and all beat both per-gene FD and fine bins.

**Transform.**  MI is computed on raw values.  A rank (copula) transform was tested.  It helped on synthetic monotone distortions, but it was slightly worse on real data (co-membership AUPR 0.127 against 0.134), so raw values are the default.

## 3. CLR calibration

Each row of the MI matrix is z-scored against its own background (sample standard deviation, negative z clipped to 0).  The two directions are combined with Stouffer's rule, (z_ij + z_ji)/√2, which has a standard-normal null.  The 2007 Euclidean rule, √(z_ij² + z_ji²), gives almost identical rankings:  Spearman correlation 0.993, and the top-5,000 lists share 4,864 pairs.

**Permutation null for edge selection.**  Each gene's samples are shuffled independently, B = 100 times, and MI and CLR are recomputed with the observed bin counts.  The null scores are pooled in a streaming histogram, and edges are selected by Benjamini–Hochberg at q = 0.05.

The null is built on CLR scores rather than raw MI.  This makes selection ask whether an edge is exceptional *relative to each gene's own background*, which is CLR's contribution.  It also keeps selection robust to global programs such as growth rate or batch.  With a global confounder, the MI null kept 45% of pairs at precision 0.03, while the CLR null kept 1.2% at precision 0.99.  A permutation-calibrated higher-criticism threshold is also available.  Efron's empirical null is logged as a future alternative.

## 4. Attention operators

The calibrated scores S define directed attention matrices A, where row i holds gene i's attention over the other genes:

| Operator | Construction |
|---|---|
| `soft10` (primary) | Softmax over each row's top 50 CLR scores.  Each row's temperature is set by bisection so that its effective neighbour count, exp(entropy), equals 10. |
| `fdr05_top10` | BH-selected CLR edges, plus each gene's own top-10 CLR neighbours (directed).  Weights are the CLR scores. |
| `pearson_top10` (control) | The same top-10 construction on \|Pearson r\|.  It separates the value of the CLR kernel from generic propagation. |

Rows are normalized to sum to 1, and genes with no neighbours keep a self-loop.

## 5. Readouts and depth

Depth t is iterated attention through the lazy operator P = (1 − α)I + αA, with α = 0.5.  Two readouts were compared:

- **Attention mass**, (Pᵗ + Pᵗᵀ)/2:  the total multi-hop attention flowing between two genes.  This is the primary readout.
- **Diffused-profile similarity**, |cor(Pᵗ Z)|, where Z is the standardized expression.  It adds little over single-step CLR.

Value transforms were also tested:

- A signed value, sign(r)·E_j.
- A conditional-expectation value (a continuous W_V), E[z_i | z_j], read off the MI estimator's own B-spline basis.

Both work on synthetic data:  they restore a repressed target's correlation from +0.78 to −0.97.  Neither helps on the real benchmark, because 99% of selected edges are positively correlated.

**Held-out depth selection.**  Experiments are split into two halves by lab (experimenter), or by random halves of cells for single-cell data.  The operator is built on one half, and t is chosen to maximize how well attention mass predicts the other half's top 1% of CLR pairs.  The search is a geometric scan (t = 1, 2, 4, …, 128) followed by golden-section refinement on the bracketing interval.  The gold standard is never used to choose t.

## 6. Evaluation

**Regulator-agnostic regulon benchmark (E. coli).**  RegulonDB regulons for every regulator class are expanded to genes:  transcription factors, sRNAs, ppGpp, other proteins and sigma factors (sigmulons).  Regulons of 5–500 genes are kept.  A gene pair is positive if the two genes share a regulon, and same-operon pairs are excluded.  The regulator is never used as a node, because regulators are generally latent.  Two metrics are reported:

- co-membership AUPR, pooled over annotated gene pairs;
- per-regulon coherence, the AUROC of within-regulon pairs against member-to-unrelated pairs.

Both are reported at strong/confirmed evidence and at all evidence.

**TF-node benchmark (BEELINE, and E. coli as a secondary result).**  The Faith, Hayete et al. (2007) and BEELINE convention scores candidate (TF, gene) edges.  The metrics are AUPRC ratio and early precision ratio (EPR).

**Uncertainty.**  Three resampling schemes are used:

- **Delete-half jackknife over genes.**  Paired across methods, for pooled AUPR.
- **Bootstrap over regulons.**  Paired across methods, for coherence.
- **Replicate-aware bootstrap on the 907-chip compendium.**  Every draw recomputes MI from scratch.  In the cluster version, experiments are resampled with replacement and one replicate chip is used per occurrence.  In the replicate version, one chip is drawn per experiment, which removes pseudo-replication.

## 7. Perturbation targets (in progress)

The influence of a perturbation P on gene g is how much g's network row changes when P's k experiments are removed.  It is calibrated per gene against random delete-k subsets, which makes it a delete-k jackknife influence.  No Gaussian (Hotelling) assumption is needed.  The comparison baseline is differential expression.  See the E. coli white paper, section 6.

## 8. Software

This is the R package `clr` on the `accuracy-review` branch:  an Rcpp core with OpenMP, and a macOS build helper.  Resampling fits run in parallel with `foreach` over the local cores minus two, with the worker count capped by RAM and the remaining cores given to OpenMP inside each worker (`analysis/parallel.R`; CITATION_LOG D28).  The analysis scripts are in `analysis/`:

- `run_m3d_regulondb.R`
- `followup.R`
- `chips_bootstrap.R`
- `perturbation.R`
- `beeline.R`
