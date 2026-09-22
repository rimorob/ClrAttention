# ClrAttention

R port of the **context likelihood of relatedness (CLR)** network-inference
algorithm, reframed as an early form of attention over continuous data
distributions. Provides:

- a B-spline mutual-information estimator (ported from the original CLR C
  sources; `bspline_mi`),
- CLR calibration with the historical Euclidean bilateral combination, the
  Stouffer `(z_ij + z_ji)/√2` variant, plus Rayleigh and KDE variants
  (`clr_calibrate`),
- data-adaptive per-gene bin counts via Freedman–Diaconis with caps
  (`fd_bins`, `bins_for_genes`),
- a permutation-bootstrap attention threshold (independent per-gene shuffles,
  streaming-histogram empirical null, Tukey higher criticism or BH-FDR),
  sparse row-stochastic attention operators, and iterative diffusion,
  orchestrated by the R6 class **`ClrAttention`**,
- a module finder (`find_modules`) and a module-evaluation harness
  (`ModuleTester`: Hungarian and independent-best-match F1, soft Sinkhorn
  assignment uncertainty).

Large expression compendia (e.g. the 65 MB M3D matrix used in development)
are **not** vendored; the demos below are fully synthetic and self-contained.

## Install

Requirements: R ≥ 4.1, a C++17 compiler, and R headers.

```sh
# Debian/Ubuntu
sudo apt install r-base r-base-dev

# macOS (with Homebrew)
brew install r
```

Clone and install R dependencies:

```sh
git clone https://github.com/rimorob/ClrAttention.git
cd ClrAttention
Rscript -e 'install.packages(c("R6", "Rcpp", "testthat", "pkgload"))'
```

Install the package from the checkout (run from the repo root):

```sh
R CMD INSTALL .
```

(Alternative: skip installing and load the package in place with
`pkgload::load_all(".")` from the repo root — the demos do this
automatically if `clr` is not installed.)

## Run the demos

Both demos are synthetic and write their plots/CSVs next to the script.

```sh
cd demo
Rscript demo_small_case.R        # t = 0 pipeline: MI -> CLR -> modules -> ModuleTester
Rscript demo_diffusion_depth.R   # Design-A sweep: fixed MI/operator, depth t = 0..30
```

`demo_small_case.R` reproduces the baseline: 24 genes × 100 samples, three
planted modules, Hungarian weighted F1 ≈ 0.96.

`demo_diffusion_depth.R` reproduces the diffusion-depth experiment: 26 genes ×
120 samples with two combinatorial bridge genes. MI and the attention operator
are estimated once; only the diffusion depth `t` varies. Outputs
(`depth_sweep/`) show hard-module recovery degrading modestly then plateauing
while bridge-gene affiliation entropy rises monotonically — diffusion
gradually reveals mixed membership.

## Run the tests

```r
testthat::test_dir("tests/testthat")
# or
R CMD build . && R CMD check clr_*.tar.gz --no-manual
```

## Layout

| Path | Contents |
|---|---|
| `R/`, `src/`, `man/`, `tests/` | the `clr` R package |
| `demo/` | runnable synthetic demos + their reference outputs |
| `docs/diffusion_depth_design_review.md` | design review for the fixed-MI diffusion experiment |
| `references.bib`, `CITATION_LOG.md` | citation index: one BibTeX entry + one justification paragraph per citable technical decision |
| `PORT_NOTES.md` | notes on porting the original MATLAB/C CLR sources |

## Status notes

- The per-depth module finder in `demo_diffusion_depth.R` (complete-linkage
  on 1−|cor|, fixed height) is **provisional scaffolding**; the exact
  per-depth similarity/module procedure is not settled. See
  `docs/diffusion_depth_design_review.md`.
- Module citations are on hold until the module/affiliation evaluation design
  is settled.
