// clr_core.hpp — Pure C++17 CLR computational core.
//
// MULTI-LANGUAGE CONTRACT (see PORT_NOTES.md):
//   * No Rcpp, no R headers, no Python headers, no R/Python assumptions.
//   * All matrices are ROW-MAJOR: element (i, j) of an (r x c) matrix lives
//     at buf[i * c + j]. A gene's samples are one contiguous row.
//   * The unit of parallel work is mi_pair(): given two genes' precomputed
//     B-spline marginal weights, compute their mutual information. A future
//     CUDA kernel (or R torch / Python equivalent) maps gene pairs to threads
//     calling exactly this function; mi_matrix() is the serial driver.
//   * Public API uses only std::vector / raw pointers / enums — nothing a
//     pybind11 binding or an Rcpp wrapper cannot call directly.
//
// Ported from InfoKit2.c (CLR backup, 2008) and clr.m 'normal' calibration.
// Where the port deviates from the original, it is marked DEVIATION.

#pragma once

#include <cstddef>
#include <vector>

namespace clr_core {

// Bilateral combination of row z-scores after negative clipping.
enum class Combine : int {
  EUCLIDEAN = 0,  // historical clr.m: sqrt(z_ij^2 + z_ji^2)
  STOUFFER  = 1   // new: (z_ij + z_ji) / sqrt(2)
};

// ---------------------------------------------------------------------------
// B-spline mutual-information machinery (ported from InfoKit2.c)
// ---------------------------------------------------------------------------

double log2d(double x);

// Uniform clamped knot vector of length num_bins + spline_order.
std::vector<int> spline_knots(int num_bins, int spline_order);

// Recursive B-spline blending value. Faithful port of SplineBlend();
// `n` is num_bins (used only for the right-edge inclusive case).
double spline_blend(int k, int t, const int *knots, double v, int n);

// Map one gene's samples into spline parameter space [0, num_bins - spline_order + 1].
void x_to_z(const double *x, double *z, std::size_t n_samples,
            int spline_order, int num_bins);

// Marginal B-spline weights for one gene: weights_out is [num_bins][n_samples]
// row-major, i.e. weights_out[bin * n_samples + s].
void gene_weights(const double *x, std::size_t n_samples, int spline_order,
                  int num_bins, double *weights_out);

// Shannon entropy (bits) of a marginal from its weights.
double marginal_entropy(const double *weights, std::size_t n_samples, int num_bins);

// Mutual information of one gene pair from precomputed marginal weights.
//   wx: [nbx][n_samples] weights of gene x,  wy: [nby][n_samples] of gene y.
//   hx, hy: their marginal entropies.
// THIS IS THE PARALLEL SEAM: embarrassingly parallel over pairs.
double mi_pair(const double *wx, const double *wy, double hx, double hy,
               std::size_t n_samples, int nbx, int nby);

// ---------------------------------------------------------------------------
// Sparse B-spline weights (the production kernel)
// ---------------------------------------------------------------------------
//
// A sample's B-spline weights are nonzero in at most `spline_order`
// CONSECUTIVE bins. SparseWeights stores, per sample s, the first bin
// first[s] and the spline_order values vals[s * order + a] for bins
// first[s] + a. The joint histogram of a pair then costs
// O(n_samples * order^2 + nbx * nby) instead of O(n_samples * nbx * nby).
// Each joint cell receives its contributions in ascending-sample order, the
// same order as the dense loop, so mi_pair_sparse() is BIT-IDENTICAL to
// mi_pair() (verified in tests/test_core.cpp).
struct SparseWeights {
  int num_bins = 0;
  int order = 0;
  std::vector<int> first;     // [n_samples]
  std::vector<double> vals;   // [n_samples * order]
};

// Compress dense [num_bins][n_samples] weights. Throws if any sample has
// nonzero weight outside a window of `order` consecutive bins.
SparseWeights sparsify_weights(const double *weights, std::size_t n_samples,
                               int num_bins, int order);

// Same contract as mi_pair(). `joint` is caller-owned scratch of at least
// nbx * nby doubles (one per thread); it is overwritten.
double mi_pair_sparse(const SparseWeights &wx, const SparseWeights &wy,
                      double hx, double hy, std::size_t n_samples,
                      double *joint);

// Serial driver over all pairs. data is [n_vars][n_samples] row-major;
// bins_per_var[v] is gene v's own bin count (adaptive binning); mi_out is
// [n_vars][n_vars] row-major, symmetric. The diagonal holds the genes'
// self-MI (as in the original mi MEX); clr_calibrate() zeroes it.
//
// Parallelism (OpenMP): the outer gene loop is parallelized; iterations are
// independent (pair (i,j), j>=i, is written exactly once, by iteration i).
// n_threads <= 0 selects default_num_threads(). Workers never call back
// into R/Python; only pure C++ runs inside the parallel region.
void mi_matrix(const double *data, std::size_t n_vars, std::size_t n_samples,
               const int *bins_per_var, int spline_order, double *mi_out,
               int n_threads = 0);

// Default thread budget: machine cores minus 2, floored at 1. Kept in the
// core (not the bindings) so R, Python, and CUDA-adjacent callers share it.
int default_num_threads();

// ---------------------------------------------------------------------------
// CLR calibration (ported from clr.m, method == 'normal')
// ---------------------------------------------------------------------------
//
//   1. Zero the MI diagonal.
//   2. Row-wise z-scores with the sample standard deviation (n - 1),
//      matching MATLAB's std() default.
//   3. Clip negative z-scores to 0 BEFORE combining (historical behavior).
//   4. Bilateral combination (Euclidean or Stouffer).
//   5. Zero the output diagonal.
//
// DEVIATION: if a row's standard deviation is 0 (constant row), the original
// MATLAB code produced NaN z-scores (0/0); we emit 0 instead. Rationale:
// NaNs poison every downstream stage; a constant MI row carries no signal.
void clr_calibrate(const double *mi, std::size_t n_vars, Combine combine,
                   double *scores_out);

}  // namespace clr_core
