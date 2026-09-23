# Citable decisions log — CLR-as-attention project

Running record of every citable methodological decision, with a one-paragraph
justification each. BibTeX entries live in `references.bib`; each section below
names its bib key(s). Decisions proposed by the user are marked **[user]**;
decisions made by the assistant are marked **[assistant]**; jointly developed
ones **[joint]**. Considered-and-rejected options are logged with their
rejection rationale, so the reasoning is recoverable. Internal empirical or
implementation decisions with no external citation are logged without a bib key.
Bibliographic details verified against publisher/PubMed records 2026-09-22;
where the most relevant citation was genuinely ambiguous it is flagged for the
user rather than guessed.

---

## D1. B-spline mutual-information estimator (historical parity) [assistant]

The MI engine is the B-spline-smoothed estimator recovered from the historical
CLR sources: continuous expression values are softly assigned to bins via
cubic B-splines (order 3) rather than hard-binned, which suppresses the
discretization bias and variance that plague naive histogram MI on the
small-to-moderate sample sizes typical of expression compendia. Keeping the
historical estimator as the parity default preserves a continuous line of
evidence back to the original CLR study: any performance difference between
the new pipeline and the historical one can then be attributed to the new
calibration, thresholding, and diffusion machinery rather than to a changed MI
engine. *Refs: `daub2004`, `faith2007`.*

## D2. Row-wise z-score calibration, clip negatives, Euclidean bilateral combination [assistant]

Historical CLR converts each row of the MI matrix to z-scores against that
gene's own MI background, zeroes negative z-scores *before* combining the two
directed scores, and fuses the pair with the Euclidean norm
sqrt(z_ij^2 + z_ji^2). Row-wise calibration is what makes CLR a
"context likelihood": a raw MI value is meaningless without knowing whether it
is surprising relative to the genes that gene i otherwise correlates with, and
the pre-combination clipping ensures that a strong signal in one direction
cannot be cancelled by an uninformative reverse direction. This exact pipeline
was verified byte-for-byte against the recovered 2007/2008 sources, so it
stays the parity default. *Ref: `faith2007`.*

## D3. Stouffer combination as a first-class alternative [user]

As an alternative to the historical Euclidean bilateral combination, the
pipeline offers Stouffer's sum-of-z-scores divided by sqrt(2). Under the null
that both directed z-scores are standard normal and independent, the Stouffer
statistic is itself exactly standard normal, which gives a clean calibration
story the Euclidean norm lacks (the Euclidean norm's null is Rayleigh, not
normal). The cost is symmetry of cancellation: a strongly positive z_ij can
be diluted by a weak z_ji, whereas the historical clipping-then-Euclidean
form never lets that happen. Offering both, with Euclidean as parity default,
turns a 2007 implementation detail into a testable modeling choice.
*Ref: `stouffer1949`.*

## D4. Adaptive per-gene bin counts, Freedman–Diaconis default [user]

The number of B-spline bins is not a fixed constant: it is computed per gene
from that gene's own empirical distribution, defaulting to the
Freedman–Diaconis rule (bin width 2·IQR·n^(−1/3)) with Scott's rule and
Sturges' rule as options, capped to a sane range. A fixed bin count is
indefensible once the method is framed as operating on *continuous*
distributions: genes with tight, near-degenerate profiles and genes with
heavy-tailed, multimodal profiles need different resolutions, and a
one-size-fits-all bin count either washes out structure or manufactures
spurious MI. FD is the default because its L2-optimality argument for
histogram density estimation is the closest classical theory to what the
B-spline MI estimator needs — a discretization that tracks the data's own
scale. *Refs: `freedmandiaconis1981`, `scott1979`, `sturges1926`.*

## D5. OpenMP parallelization of the MI triangle, threads default cores−2 [assistant]

Pairwise MI is embarrassingly parallel over the upper triangle of the gene
matrix, so the C++ core parallelizes it with OpenMP dynamic scheduling (pairs
have non-uniform cost once bin counts are adaptive, hence dynamic rather
than static chunks). The default thread count is max(1, hardware cores − 2),
reserving headroom on shared machines — a pragmatic choice for a VM and for
the user's workstation alike, since MI at 7,000 genes is a ~30-minute
single-threaded job that becomes a ~2-minute job at 16+ threads. This is an
engineering decision rather than a statistical one, cited to the OpenMP
specification for reproducibility of the parallel semantics.
*Ref: `openmp2021`.*

## D6. Attention operator: row-stochastic matrix from thresholded CLR scores [joint]

The conceptual reframe of the whole project: CLR's calibrated pairwise MI
matrix is an attention matrix over *continuous* distributions, where
conventional attention operates over nominal/token or dictionary
distributions. The attention operator is built by thresholding the CLR score
matrix and row-normalizing, so each gene's row is a probability distribution
over the genes it attends to — exactly the query-key-value pattern with MI in
place of the scaled dot product and, in the current design, the gene
expression profiles themselves as the values. The reframe is what justifies
importing the transformer-era toolkit (sparsity, diffusion, depth) into a
2007 network-inference algorithm. *Ref: `vaswani2017`.*

## D7. Permutation null: independent per-gene shuffles [user]

Significance of a CLR edge is judged against a null built by independently
permuting each gene's sample vector, recomputing MI and the full CLR
calibration, and pooling scores across permutations. Shuffling each gene
separately destroys all genuine co-expression while preserving every gene's
marginal distribution, which is precisely the right null for the question
"is this pair's MI surprising given what these two genes look like
individually." This follows the precedent set in information-theoretic network
inference by ARACNE, which established permutation-based MI significance as
the standard way to separate direct statistical dependencies from background.
Confirmed with the user 2026-09-22: ARACNE stands as the citation.
*Ref: `margolin2006`.*

## D8. Pooled empirical null via streaming histogram; beta and KDE rejected for tail inference [assistant]

The null distribution is the pooled empirical survival function over all
permutation CLR scores, accumulated in a fixed 20,000-bin streaming
histogram — no parametric fit, no cached null matrices. The user proposed
fitting a beta distribution and, later, a KDE to the null; both were rejected
for the specific use we need them for: extreme-tail probabilities. A beta fit
is dominated by the bulk of the null and can be arbitrarily wrong orders of
magnitude out in the tail where the p-values actually live, and a KDE's tail
is an artifact of its bandwidth choice, not of the data. The honest
alternative for tail inference is peaks-over-threshold extreme-value theory,
noted here as the principled fallback if the empirical null's finite
resolution (B permutations × pairs) ever becomes the binding constraint.
**Considered and rejected:** parametric beta null, KDE-based tail p-values.
*Ref: `coles2001`.*

## D9. Higher-criticism threshold with sqrt(2 log log M) significance guard [assistant]

The attention threshold is selected by Tukey's higher criticism: order the M
empirical p-values, maximize the standardized exceedance
sqrt(M)·(i/M − p_(i))/sqrt(p_(i)(1−p_(i))) over the first half of the order
statistics, and keep edges up to the maximizing index — but only if the
maximized HC* clears the sqrt(2 log log M) asymptotic null scale, otherwise
warn. HC is the right tool because attention sparsification is a sparse
heterogeneous mixture problem: a small unknown fraction of pairs carry real
signal amid a vast null background, which is exactly the regime where HC is
known to achieve the optimal detection boundary adaptively, without knowing
the sparsity level or signal strength in advance. The sqrt(2 log log M)
cutoff is not a tuning knob: it is the almost-sure growth rate of the HC
objective under the global null, so requiring HC* to exceed it is what keeps
the procedure from hallucinating structure in pure noise.
*Ref: `donohojin2004`.*

## D10. Benjamini–Hochberg FDR as the alternative threshold [user]

Alongside higher criticism, the pipeline offers the Benjamini–Hochberg
procedure at q = 0.05. Where HC is tuned for *detection* — is there any
signal at all, and where does it concentrate — BH is tuned for *selection
with an error budget*: it returns the largest edge set whose expected false
discovery proportion is controlled. The two answer different questions a
practitioner actually asks ("which edges can I trust at 5% FDR" vs. "where is
the evidence of structure strongest"), and on the synthetic benchmark they
behave as theory predicts: FDR is more liberal at the margin (5/5 true edges,
1 false positive) while HC is conservative (4/5, 0 false positives). Offering
both keeps the method honest about the detection/selection distinction.
*Ref: `benjaminihochberg1995`.*

## D11. Bonferroni rejected for attention sparsification [assistant]

Family-wise error control was considered and rejected: at M ~ 25M pairs for
7,000 genes, Bonferroni's threshold is so conservative that only the most
extreme edges survive, which defeats the purpose of building a dense-enough
attention operator for diffusion. Attention needs a *population* of edges,
not a handful of certainties; FDR and HC both serve that need with stated
error semantics. **Considered and rejected:** Bonferroni correction. (No new
citation; follows from D9/D10.)

## D12. Diffusion iteration E^{(t+1)} = ((1−α)I + αÂ)E^{(t)} [joint]

The value/concept aggregation that historical CLR lacked is supplied by lazy
diffusion of the expression matrix through the row-stochastic attention
operator: each step mixes a gene's profile with the attention-weighted average
of its neighbors, with α controlling how much mass moves per step. The cited
formulation is Zhou and Schölkopf's regularization framework for learning
from graph data, which derives exactly this iteration —
f ← αSf + (1−α)y, with the lazy random walk transition
P = (1−α)I + αD⁻¹W — from a smoothness-plus-fidelity objective on the graph,
gives its closed form, and analyzes the stationary distribution it converges
to. That stationary distribution is precisely the oversmoothing failure mode:
iterated diffusion collapses toward a representation constant on each
connected component, which is why a depth-diagnosis statistic (still open) is
needed to stop before the signal is washed out. Confirmed with the user
2026-09-22: this is the intended diffusion citation. *Ref: `zhouscholkopf2004`.*

## D13. Synthetic benchmark: 12 genes, planted nonlinear modules [assistant]

Validation uses a 12-gene × 300-sample synthetic with five planted undirected
edges mixing quadratic, sinusoidal, saturating, and signed-linear
dependencies against seven independent null genes. The design is deliberately
adversarial to linear methods: the quadratic edge exists specifically so that
Pearson correlation fails where MI succeeds, which is the load-bearing claim
of an information-theoretic attention. Twelve genes keeps the full
permutation-threshold pipeline runnable in seconds so the test suite stays
deterministic and fast; scaling behavior is measured separately (D14), not in
the unit tests. (Internal design decision; no external citation.)

## D14. Log-time scaling model for MI runtime extrapolation [assistant]

Single-threaded MI build time was measured at 40/80/160/320 genes (900
samples, fixed bins, four replicates) and fit in log-log space: slope
1.91 ± 0.01, residual SE 0.02, extrapolating to ~30 min at 7,000 genes with a
95% prediction interval of roughly 28–33 min. The near-quadratic exponent is
expected — the work is dominated by the ~n²/2 pair evaluations, with a slight
sub-quadratic discount from per-pair overhead amortization. The extrapolation
spans ~22× beyond the largest measured matrix and must be re-checked at
1,000+ genes before it is used for production planning; it exists to size the
move to the user's 28-core workstation, where the same computation should
take minutes. (Internal empirical result; no external citation.)

## D15. M3D E. coli v4 build 6 as the lead dataset [joint]

The lead real-data setting is the Many Microbe Microarrays Database E. coli
compendium (4,297 genes × 907 arrays, uniformly RMA-normalized), chosen
because it is the largest uniformly processed bulk bacterial expression
resource with structured experimental metadata and a curated gold-standard
regulatory network for out-of-cohort validation. Uniform normalization across
laboratories is what makes the cross-condition/cohort/lab generalization
claims testable without new experiments. Provenance caveat, recorded for the
paper's methods section: the probe-to-gene summarization step is undocumented
in the mirror and needs publication-grade validation before results rest on
it. *Ref: `faith2008m3d`.*

## D16. Compressed sensing and significant-components theory as the theoretical backbone [user]

The user expects compressed-sensing theory to do heavy lifting, and chose
the phase-transition line over the L1-recovery/RIP line: Donoho and Tanner's
universality of phase transitions, which locates the sharp boundary in the
undersampling/sparsity tradeoff where sparse structure goes from recoverable
to unrecoverable — the natural language for asking when a thresholded
attention operator preserves regulatory signal and when it collapses.
Alongside it, two results on the number of significant components: the
Gavish–Donoho optimal hard threshold (4/√3, scaled by noise level and matrix
aspect ratio) for deciding how many singular values of a noisy matrix are
real, and Horn's parallel analysis, the empirical procedure of retaining
only eigenvalues that exceed those of matched random data. Both bear
directly on the open diffusion-depth question — they are principled answers
to "how many modes are signal before the representation is noise," which is
what the depth-diagnosis statistic must operationalize to stop diffusion
before oversmoothing. Confirmed with the user 2026-09-22.
*Refs: `donohotanner2009`, `gavishdonoho2014`, `horn1965`.*

---

# Accuracy review, 2026-09-23

The entries below record the decisions taken in the accuracy review that
preceded the M3D x RegulonDB run. Each was motivated by a defect that was
confirmed empirically. The scripts are in `validation/` and the regression
tests are in `tests/testthat/`. Where an entry changes an earlier decision,
the earlier entry is left as written and the change is recorded here.

## D17. Rank (empirical-copula) transform before MI; FD bins computed on the ranks [joint]

Each gene is replaced by its within-gene ranks before B-spline MI
estimation. Mutual information is invariant to strictly monotone
transformations of either variable, so the transform leaves the quantity
being estimated unchanged and changes only the estimator. It also gives
every gene the same (uniform) marginal, which removes a confirmed bias
pathway in the D4 design. Freedman–Diaconis on raw values gives
heavy-tailed genes many bins: the median was 23 bins for Gaussian genes and
50, the clamp, for t3 and lognormal genes. The finite-sample upward bias of
MI grows with bins_i × bins_j: at N = 907 the null MI was 0.012 bits with
10 bins, 0.09 with 25 and 0.32 with 50. That bias survived CLR calibration
(Spearman correlation between null CLR score and bins_i × bins_j was 0.32–0.51
across gene mixes), so
heavy-tailed genes became false hubs. Freedman–Diaconis applied to the ranks
gives ⌈N^(1/3)⌉-type counts that are identical across genes (10 at N = 907).
On a synthetic with M3D-like N = 907, with monotone heavy-tail distortions,
the edge AUPR was 0.24 for raw values with FD bins, 0.57 for raw values with
10 bins (2007 parity) and 0.955 for ranks with FD bins. An independent
QC replication got 0.31, 0.58 and 0.94. Caveat: the synthetic uses
monotone distortions, which favour the rank transform by construction. The
real-data comparison of the four configurations is part of the M3D run. The
Gaussian-copula MI literature makes the same argument: separate the
marginals from the dependence. Historical parity remains available as
`transform = "none", bins = 10`. **Default changed with user approval.**
*Ref: `ince2017`.*

## D18. Benjamini–Hochberg is the default threshold; HC is permutation-calibrated [assistant, user-delegated]

The earlier HC implementation (D9) selected noise. On 60 independent genes
it kept 11–32% of pairs, and the sqrt(2 log log M) guard fired in only 1 of
6 runs, because that bound is the typical size of the null HC*, not a
significance cutoff. The rewrite makes four changes:

- HC is evaluated at the edges of the histogram bins, so tied scores are
  never split.
- The search is restricted to i ≤ max(α0·M, 10), with α0 = 0.1.
- HC* is gated against the (1 − level) quantile of leave-one-out HC* values
  computed from the B permutation replicates. Each replicate is scored
  against the pool of the other B − 1 replicates.
- The HC+ floor p ≥ 1/M is dropped: with permutation p-values it removes
  exactly the strongest edges, and the calibration already absorbs the
  extreme-p blow-up the floor guarded against.

On 20 pure-noise runs, calibrated HC selected any edge 0–10% of the time,
consistent with the nominal 5%. BH at q = 0.05 found slightly more true
edges at a similar false-positive count, so BH is the default. The user
delegated this choice. Donoho & Jin (2008) is the right citation for
using HC as a *threshold*; Donoho & Jin (2004) covers detection.
*Refs: `donohojin2008`, `benjaminihochberg1995`.*

## D19. The permutation null stays on CLR scores; the MI null is an option only [user]

For one iteration the null was moved to raw MI, because CLR scores are not
exactly comparable between observed and permuted data. In a 24-gene toy
with 6–8-gene modules, a gene's own module inflated its row background, so
true-edge CLR scores (median 2.3) fell below the null per-replicate maximum
(4.4) and nothing was selected. The user pointed out that the selection has
to answer "exceptional relative to each gene's background", which is
CLR's contribution. With a global confounder that answer is the only usable
one. In a synthetic with a global factor (loading 0.35), BH on the MI null
kept 45% of pairs at precision 0.03, while BH on the CLR null kept 1.2% at
precision 0.99. That is seed 3, which is enshrined as a test. An
independent QC replication with other seeds and designs found the same
direction at a smaller magnitude: the MI null kept 7–18% of pairs at
precision 0.07–0.19, and the CLR null kept 0.8–0.9% at precision
0.94–0.98. The CLR null is the default. The small-G compression is
documented as a limitation of the regime and is not treated as a reason to
switch statistics. If the permutation null misbehaves at compendium scale,
the principled refinement is an empirical null fitted to the central bulk of
the observed CLR scores, which is not implemented. *Ref: `efron2004`.*

## D20. The permutation replicates reuse the observed per-gene bin counts [user]

`estimate_mi()` records the bin vector it actually used, and every
permutation replicate reuses it verbatim, together with the spline order
and the transform. Because FD depends only on the IQR and the range, which a
permutation preserves, the replicates already reproduced the same counts.
The change makes that guarantee explicit and independent of the binning
rule. It does not fix D17's bias on its own: pooling nulls across pairs
with different bins_i × bins_j still mixes different null distributions.

## D21. Correctness fixes with no methodological choice [assistant]

- **Constant genes.** A constant gene divided by zero in `x_to_z`. MI then
  became H(Y) against every other gene, so the gene became the top hub.
  Constant genes are now rejected in both R and the C++ core.
- **Isolated genes.** Genes with no selected edge now get a self-loop.
  Previously their row of P was (1 − α)e_i, so they decayed as (1 − α)^t,
  and tau = Inf zeroed the whole embedding.
- **Standardization.** `diffuse()` row-standardizes E^(0) by default, so
  diffusion mixes expression shapes rather than baseline levels.
- **Histogram p-values.** These now count the null values that share the
  observed score's bin. The earlier code was anti-conservative by up to one
  bin.
- **BH ties.** The cutoff returns the smallest score in the whole rejected
  set, so tied histogram p-values no longer drop members of the last
  rejected block. The QC agent found this.
- **Calibration guard.** `select_threshold(statistic = "clr")` refuses KDE
  and Rayleigh calibrations. KDE scores are ≤ 0 and put every null value in
  one bin.
- **Sparse kernel.** Each sample has at most `spline_order` nonzero
  B-spline weights. The joint histogram now accumulates only those, in the
  same sample order, so the result is bit-identical to the dense kernel
  (max |Δ| = 0 on 200 genes) and 5–20× faster. Memory drops from
  O(G·N·bins) to O(G·N·order).

## Sign blindness of diffusion (the missing W_V): implemented 2026-09-23 as `diffuse(values = ...)`, see D24

MI is sign-agnostic, but diffusion averages raw profiles. A repressor and
its target therefore partially cancel: after five steps in the synthetic
test, g7 (−1.5·tfb) correlated +0.78 with tfb. Transformer attention avoids
this through a learned value projection W_V. The continuous analogue
proposed here is value_{j→i} = E[x_i | x_j], read off the B-spline joint
histogram that MI already computes. It would handle negative and
non-monotone links (such as the quadratic edge) without needing "meanings".
Not implemented; flagged for a design decision.

## D22. The default MI is back on raw values, with one common bin count taken as the median of per-gene Scott counts [user]

D17 is **reversed for the default**. The rank (copula) transform stays
available as `transform = "rank"`.

**The genome-wide result that decided it.** On the full M3D compendium
(4,297 genes × 907 arrays), scored by co-membership AUPR on the
strong/confirmed RegulonDB regulator-agnostic benchmark (base rate 0.094),
CLR on raw values beat CLR on ranks:

| Setting | CLR AUPR |
|---|---|
| Raw values, 10 common bins | 0.134 |
| Raw values, per-gene FD bins | 0.129 |
| Ranks, 10 bins | 0.127 |

This matches the user's long-standing experience that MI on raw values
outperforms the copula transform on average. D17's synthetic evidence was
built from monotone distortions, which favour the rank transform by
construction, and should not have decided a default.

**What D17 got right, and what is kept.** The real defect was that
per-gene bin counts differed between genes. The fix kept from D17 is a
common bin count for every gene. The historical rule, as the user
reports it, is to compute each gene's optimal count and use the median
for all genes, with spline order 3. `bins = "median_scott"` is the
default. Scott's constant of 3.49 is the same as Wand's zero-stage rule,
which the recovered sources implemented as `calcNumBins`.

**Bin counts on M3D.** The median rule gives 21 bins at N = 907 with Scott
and 33 with FD. The historical "6–10 bins" corresponds to the smaller N of
2007-era datasets. Whether 21 bins beats a fixed 10 at N = 907 is measured
by the `none_median_scott` configuration of the M3D run.

## D23. The bin count is chosen for the 2-D (joint) histogram, and `hg` is the default [joint]

The user objected to 21–33 bins per axis at N = 907. At 33 bins the joint
table has 1,089 cells for 907 points, fewer than one point per cell. The
bin budget belongs to the joint histogram, because MI's variance and bias
come from H(X,Y). The user recalled that Daub et al. used about 6 bins.

Two rules for the joint histogram are implemented:

- **Scott's multivariate normal-reference rule** (`median_scott2d`).
  For d = 2, h = 3.5·sd·N^(−1/(d+2)) = 3.5·sd·N^(−1/4). The count comes out
  at 12 at N = 907 and 9 at N = 466, as the median over M3D genes.
- **Hacine-Gharbi et al. (2012)'s low-bias rule** (`hg`). This rule is
  derived specifically for histogram MI. It is evaluated at ρ = 0 so that
  every pair shares one discretization, and it is exact for the null every
  pair is tested against: B = round(√((1 + √(1 + 24N))/2)). That gives 9
  bins at N = 907 and 7 at N = 466.

The formula is quoted from memory of the paper and of López de Prado's
treatment of it. It must be checked against the paper before publication.

**Why this does not change the bias argument of D17 and D22.** B-spline
order 3 spreads each point over 3 × 3 cells, which lowers variance compared
with a hard histogram. The rules are therefore, if anything, conservative
for this estimator.

**Default.** `hg` is now the default. It is also the pre-registered
primary configuration of the M3D run, so the choice is fixed before seeing
that run's results. Scott-2d, the historical median-Scott rule and fixed 6,
8, 10, 12 and 16 bins are reported as sensitivity analyses. They are not
used to pick the rule after the fact.
*Refs: `hacinegharbi2012`, `scott1979`, `daub2004`.*

## D24. Value transforms for diffusion: signed and conditional-expectation [joint]

This implements the W_V proposal. `diffuse()` now takes one of three
values:

- `values = "raw"` passes E_j, the previous behaviour.
- `values = "signed"` passes sign(cor_ij)·E_j.
- `values = "conditional"` passes f_ij(E_j), where
  f_ij(u) = E[z_i | z_j = u].

f_ij is the nonparametric regression of gene i on gene j. It is read off
the same B-spline basis the MI estimator uses: gene j's bin count and
spline order, with `weights_at()` exported from the C++ core. f_ij is
fitted once on the data and re-applied to the diffused profiles at every
step. Values outside z_j's range are clamped to it. The messages are
pair-specific, which makes this edge-conditioned message passing. No
parameters are learned.

**Synthetic check (12 genes, top-3 operator, 10 steps).** The table gives
the correlation of each diffused target with its true driver.

| Target | raw | signed | conditional | t = 0 |
|---|---|---|---|---|
| Repressed target (g7 vs tfb) | +0.78 (flipped) | −0.95 | −0.97 | −0.83 |
| Quadratic target (g2 vs tfa²) | 0.26 | 0.25 | 0.51 | 0.90 |

**M3D, 1,500-gene quick run.** Only 2 of 3,032 selected edges are
negatively correlated. Signed diffusion is therefore identical to raw
diffusion to 4 decimals, and conditional diffusion is slightly below raw
(co-membership AUPR on strong/confirmed evidence 0.192 vs 0.196 at t = 8).
Sign cancellation is not what limits diffusion on this benchmark with this
operator. Results on the full compendium are pending.

## D25. Stouffer is the default bilateral combination, and the default run uses the 466-experiment compendium [user]

The default changes from `combine = "euclidean"` to `"stouffer"` on
principle. With independent standard-normal directed z-scores, Stouffer's
null is standard normal; the Euclidean norm's null is Rayleigh-type. The
2007 Euclidean combination remains available (`combine = "euclidean"`, and
the `parity2007` and `none_hg_eu` run configurations).

The change has almost no empirical effect on M3D (466 experiments, default
MI settings):

- Across pairs with any positive score, the Spearman correlation between
  the two combinations' scores is 0.993.
- The two top-5,000 lists share 4,864 pairs.
- The share of negatively correlated pairs among the top 100,000 is 5.9%
  under Stouffer and 6.3% under Euclidean.
- On the full 907-array data, co-membership AUPR differs by 0.0003.

The scarcity of anti-correlated edges comes from the data, not from the
combination. MI is sign-blind, and both combinations use clipped,
non-negative z-scores.

`tools/run_m3d.sh` now runs only the 466 replicate-averaged experiments by
default (`SET=chips` runs the 907-array set). The 907-array set counts
technical replicates as independent samples, which inflates N: it pushes
the bin rule from 7 to 9 bins and overstates evidence for condition-specific
genes. Both sets cover the same experiments and all 4,297 genes. The
optional MI-null comparison (`--mi-null`) is no longer part of the default
run.
