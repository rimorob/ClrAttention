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
