# CLR R Port — Port Notes

Package `clr`, at the root of this repository.
R6 port of the CLR (Context Likelihood of Relatedness) algorithm, reframed as
iterative attention over continuous data distributions.

## 1. Layering

```
R/                          R6 orchestration + functional core (pure R)
  clr_attention.R             ClrAttention class (pipeline, validation, caching)
  mi.R / calibrate.R / bins.R functional core: bspline_mi(), clr_calibrate(),
                              fd_bins(), bins_for_genes()
src/rcpp_bindings.cpp       THIN Rcpp wrappers only. Convert R column-major
                            matrices <-> core row-major buffers at the boundary,
                            call the core, return. No math lives here.
src/RcppExports.cpp         Hand-maintained .Call registration
                            (mirrors Rcpp::compileAttributes() output, so no
                            build-time codegen step is needed).
src/core/                   PURE C++17. No Rcpp, no R headers, no R assumptions.
  clr_core.{hpp,cpp}          B-spline MI + CLR calibration, ported from
                            InfoKit2.c and clr.m.
  tests/test_core.cpp         Standalone g++ harness (no R needed).
```

**Backend cleanliness rules** (hard constraint — a Python fork is planned):
- `src/core` signatures use only `double*`, `int*`, `std::vector`, enums.
  No `Rcpp::` types, no `SEXP`, no R API anywhere in core.
- Layout is explicit at every boundary and documented in `clr_core.hpp`:
  **row-major everywhere**, `buf[i * ncols + j]`, one gene = one contiguous row.
  The R bindings transpose R's column-major matrices at the boundary
  (exactly as the original `mi.c` MEX did for MATLAB).
- MI is computed in `double` throughout; the original cast the final matrix
  to `float` for MATLAB. We keep `double` — strictly more accurate, and the
  calibration step is scale-sensitive.
- New code must not assume R shenanigans: no 1-based indexing, no column-major
  shortcuts, no R RNG, no R error handling in core (core throws
  `std::invalid_argument`; each binding layer translates).

## 2. The parallel seam

The unit of parallel work is **`clr_core::mi_pair()`**:

```cpp
double mi_pair(const double* wx, const double* wy, double hx, double hy,
               std::size_t n_samples, int nbx, int nby);
```

Given two genes' precomputed B-spline marginal weights (`wx`, `wy`) and their
marginal entropies, it computes the joint entropy and returns MI. It is
embarrassingly parallel over gene pairs: no shared state, no writes except the
caller's output cell. `mi_matrix()` is the serial driver (precompute weights
once per gene, then loop pairs) — it exists so there is exactly one place to
swap in a parallel implementation.

**CUDA plan** (future, not implemented): map gene pairs `(i, j)` to threads;
each thread calls the `mi_pair` logic on the precomputed weight buffers
(uploaded once). The joint-entropy inner loop is a reduction over samples —
the natural CUDA shape is one block per pair or a strided grid over pairs.
Preferred route is R's **torch** (`rtorch`) for the R package; fallback is a
Python interface via **reticulate** driving a CUDA kernel or torch there.
Either way, `mi_pair` is the function being ported to the device.

**Python fork plan**: `src/core/` compiles unchanged under pybind11 —
`mi_matrix`, `clr_calibrate`, `mi_pair` take raw buffers, so the pybind11 layer
is the same kind of thin boundary code as `rcpp_bindings.cpp` (NumPy arrays
are row-major by default, so the transpose step the R side needs mostly
disappears). Keep core free of Rcpp-isms so this stays true.

## 3. Parity decisions (vs. verified `clr.m` / `InfoKit2.c`, byte-identical in both backups)

- `estimate_mi(bins = 10)` + `calibrate(method = "normal", combine = "euclidean")`
  reproduces `clr.m` exactly: diagonal-zeroed MI → row-wise z-scores with the
  **sample** sd (n − 1, matching MATLAB `std`) → **negatives clipped to 0
  before combining** → bilateral Euclidean `sqrt(z_ij^2 + z_ji^2)` →
  diagonal zeroed. Verified by independent re-derivation in both the C++
  harness and the testthat suite.
- B-spline kernel (`SplineBlend`, `SplineKnots`, `xToZ`, `findWeights`,
  `hist1d`/`hist2d`, `entropy1d`/`entropy2d`, `miSubMarix`) ported faithfully,
  including the right-edge inclusive case (`|v − u[k+1]| < 1e-10`) and the
  negative-clamp on blend values. `spline_order = 3`, the historical default.
- **Deviations** (deliberate, documented in code):
  1. `double` instead of `float` for the MI matrix.
  2. Constant MI row (sd = 0): MATLAB produced NaN z-scores (0/0); we emit 0.
     NaNs poison every downstream stage; a constant row carries no signal.
  3. `method = "kde"` is a structural port of `clr.m`'s KDE branch, not
     bit-exact: R `density()` and MATLAB `ksdensity` use different default
     bandwidth selectors. Same grid (1000 pts, row min→max), same normalized
     empirical CDF, same `log(A) + log(A')` combination with `-Inf → 0`,
     positives → 0.
- `method = "rayleigh"` is a faithful pure-R port: `raylfit` MLE
  `σ = sqrt(mean(x^2)/2)`, `raylcdf = 1 − exp(−x²/2σ²)`,
  combination `A + A' − A·A'`, diagonal zeroed.

## 4. Stouffer addition (new, not in the historical code)

`combine = "stouffer"` computes `(z_ij + z_ji) / sqrt(2)` **after the same
negative clip** as the Euclidean path. Euclidean stays the default for
historical parity. Rationale for the clip-then-combine order: it preserves the
historical semantics where only above-background likelihoods contribute; the
combination choice then only changes how the two directional evidences merge.

## 5. Adaptive binning (Freedman–Diaconis default)

The historical code fixed 10 bins for every gene (the old `calcNumBins` /
Wand zero-stage-rule machinery existed but `mi()` defaulted to fixed 10).
The R port estimates the count **per gene** from the gene's own distribution:

- `"fd"` (default): `h = 2·IQR·n^(−1/3)`, `bins = ceil(range/h)`.
- `"scott"`: `h = 3.5·sd·n^(−1/3)`.
- `"sturges"`: `ceil(log2(n)) + 1` (also the fallback for degenerate genes:
  zero IQR/sd, < 2 finite values).
- Integer: fixed count for all genes (parity mode; `10` = historical default).
- Adaptive counts are clamped to **[5, 50]** (user-approved "sensible" default).

Per-gene counts are sound for the B-spline estimator because the joint
histogram is the outer product of the two marginal B-spline bases: pair
`(i, j)` uses an `n_bins[i] × n_bins[j]` grid, so each gene keeps a consistent
marginal discretization across all its pairs. When all counts are equal, the
computation is exactly the historical fixed-bin one.

## 6. R6 API (per approved proposal; Hadley good-form rules)

```r
fit <- ClrAttention$new(data)   # genes x samples; validated in $initialize()
fit$estimate_mi(bins = "fd", spline_order = 3)$
    calibrate(method = "normal", combine = "euclidean")$
    build_operator(k = 50, tau = NULL, alpha = 0.5)$
    diffuse(steps = 10)
fit$mi; fit$clr_scores; fit$operator; fit$embedding; fit$trajectory; fit$params
```

- One class for the pipeline; UpperCamelCase class, snake_case methods.
- Strict read-only active bindings: accessing a stage that hasn't run errors
  with `not available yet: run $estimate_mi() first` (etc.).
- Side-effect methods return `invisible(self)` → chaining; re-running an
  early stage invalidates downstream caches.
- `build_operator`: top-`k` per row (positive scores only) or threshold
  `tau`, then row-normalized to a stochastic matrix; all-zero rows stay zero.
- `diffuse`: `E <- ((1−α)I + αÂ)E`, full trajectory cached privately;
  `$embedding` is the final iterate. `reset_diffusion()` drops the trajectory.
- `diagnose_depth()` (test/validate log-ratio stopping statistic) is
  **not yet implemented** — the statistic design is still open; it will land
  here once finalized.

## 7. Verification status

- **C++ core**: verified standalone with g++ 13.3 (`-std=c++17 -Wall -Wextra`),
  `src/core/tests/test_core.cpp`: MI symmetry / non-negativity / self-MI,
  `mi_pair` seam == driver, adaptive bins symmetric, Euclidean/Stouffer match
  independently re-derived z-score formulas, correlated pairs outrank noise.
  **ALL CORE TESTS PASSED.**
- **R package**: testthat suite written (`tests/testthat/`); `R CMD check`-style
  load + full suite run pending R availability (R was installing in the
  background at build time).
- Not yet done: real-data run (e.g. M3D E. coli), R-level parity against the
  MATLAB outputs, performance profiling, CUDA/torch work, Python fork.

## 8. Open items / before release

- Fill in `DESCRIPTION` authorship (currently a placeholder).
- Decide: sparse `Matrix` for the operator at scale (currently dense).
- `diagnose_depth()` once the stopping statistic is finalized.
- R-level parity check against historical MATLAB `clr` outputs on a small
  fixture (needs MATLAB or Octave run of the recovered code).

## 9. Permutation threshold selection (new, 2026-09-22)

`ClrAttention$select_threshold()` replaces the hand-tuned k/tau with a
data-driven attention threshold, per user-approved design:

- B independent shuffles (each gene's expression vector permuted
  independently), MI + CLR rebuilt with the fitted settings each time.
- Pooled empirical null accumulated in a *streaming* fixed-bin histogram
  (20k bins over [0, hi), hi from the first bootstrap; overflow bin;
  double counts since B*M can exceed 2^31). Null matrices are never kept.
- P-values from the empirical survival function with a
  (1 + #{null > s})/(1 + n_null) continuity correction. No parametric fit:
  beta was rejected (wrong support/tail); KDE rejected for tail inference
  (kernel-driven, anti-conservative at depth); Rayleigh/normal are the
  natural parametric nulls if one is ever wanted.
- Cutoff by Tukey's higher criticism ("hc", default -- Donoho & Jin,
  argmax of sqrt(M)(i/M - p_(i))/sqrt(p_(i)(1-p_(i))) over i <= M/2, with a
  sqrt(2 log log M) no-signal warning) or Benjamini-Hochberg FDR ("fdr",
  level q). Bonferroni was rejected as too conservative for attention.
- Each gene keeps however many connections survive: per-gene degree
  adapts, no k. `build_operator()` uses the selected threshold when tau is
  NULL; tau = Inf (FDR with no discoveries) yields an empty operator.
  **Corrected 2026-09-23:** the (1-alpha)I residual did *not* handle this
  gracefully -- isolated genes decayed as (1-alpha)^t. Genes with no
  selected edge now carry a self-loop (section 10).
- Cost is ~B MI builds; embarrassingly parallel across bootstraps only
  via threads within each build. Uses R's RNG (set.seed for repro).
- Enshrined in tests/testthat/test-synthetic.R: on the 12-gene planted
  network (B=20), HC keeps 4/5 true edges with 0 false positives of 66
  pairs; FDR keeps 5/5 true + 1 FP.

## 10. Accuracy review (2026-09-23)

The review was independent: a separate QC agent reproduced every defect.
Scripts are in `validation/` and the reasoning is in `CITATION_LOG.md`,
entries D17 to D21.

**The estimator is correct.**  The B-spline MI estimator matches an
independent `splines::splineDesign` implementation of Daub et al. 2004 to
within 3e-15.

**Defaults that changed**

| Setting | Old default | New default |
|---|---|---|
| MI transform | none (raw values) | `transform = "rank"` (empirical copula) |
| Bins | FD on raw values | FD computed on the ranks, which gives equal bins for every gene |
| Threshold method | `hc` | `method = "fdr"` (BH, q = 0.05) |

Two further changes accompany these.  HC is now calibrated against
leave-one-out permutation HC* (arguments `hc_alpha0`, `hc_level`).
`diffuse(standardize = TRUE)` is the new default.  To reproduce the 2007
results, use `transform = "none", bins = 10`.

**The permutation null.**  The default is `statistic = "clr"`.
`statistic = "mi"` is an option intended for small or confounder-free data.
Section 9's shuffles are unchanged apart from the following:

- They reuse the observed per-gene bin vector (`params$mi$bins_used`) and
  transform.
- P-values count the null values in the observed score's own bin.
- One histogram is kept per replicate, which the HC gate needs.
- `hi` is set from the first replicate's maximum, with no floor of 10;
  observed scores above it are compared against the overflow count.
- The new `$edges` binding gives the selected pair set.

**Operator and diffusion**

- A row with no selected edges gets a self-loop, so every row of P is
  stochastic.
- E^(0) is the row-standardized data.
- Sign blindness (anticorrelated neighbours cancel) is still open; see the
  W_V note at the end of `CITATION_LOG.md`.

**Robustness**

- Constant genes are rejected in R, and the C++ core throws on them.
- `clr_openmp_info()` reports whether the build is parallel. Apple clang
  without libomp builds serially; `tools/setup_mac.sh` handles this.

**Performance.**  `mi_pair_sparse()` and `SparseWeights` are the production
kernel.  The result is bit-identical to `mi_pair()`, which is checked in
`test_core.cpp`, and 5-20x faster.  Memory drops to O(G·N·order).  The dense
`mi_pair()` stays as the reference seam.

**Real-data run.**  `analysis/run_m3d_regulondb.R`, driven by
`tools/run_m3d.sh`, does the following:

- Computes edge-level precision-recall in the Faith, Hayete et al. 2007 style (every known
  TF × gene pair) for four configurations.
- Selects edges by permutation.
- Runs a Design-A depth sweep that scores per-TF regulon average precision
  from attention mass P^t and from the correlation of diffused profiles.
- Saves artifacts following section 9 of the design review.
