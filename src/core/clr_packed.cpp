// clr_packed.cpp — tiled (band-parallel) MI into a packed upper triangle and
// CLR calibration in place. Pure C++17; see clr_core.hpp for the layout and
// contracts. Built on the band primitive of the original miSubMarix()
// (InfoKit2.c): a band of consecutive rows is a contiguous slice.

#include "clr_core.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace clr_core {

std::vector<std::size_t> equal_area_bands(std::size_t n, std::size_t n_bands) {
  if (n < 2) throw std::invalid_argument("n must be >= 2");
  if (n_bands < 1) n_bands = 1;
  if (n_bands > n - 1) n_bands = n - 1;
  const double total = static_cast<double>(packed_size(n));
  std::vector<std::size_t> start;
  start.push_back(0);
  double cum = 0.0;
  std::size_t b = 1;
  for (std::size_t i = 0; i + 1 < n && b < n_bands; ++i) {
    cum += static_cast<double>(n - 1 - i);          // pairs in row i
    if (cum >= total * static_cast<double>(b) / static_cast<double>(n_bands)) {
      start.push_back(i + 1);
      ++b;
    }
  }
  start.push_back(n);
  return start;
}

template <typename T>
void mi_packed(const double *data, std::size_t n_vars, std::size_t n_samples,
               const int *bins_per_var, int spline_order, T *out,
               int n_threads) {
  if (n_vars < 2) throw std::invalid_argument("n_vars must be >= 2");
  if (n_samples < 2) throw std::invalid_argument("n_samples must be >= 2");
  // Marginal weights and entropies once per gene (as mi_matrix / miSubMarix).
  std::vector<SparseWeights> W(n_vars);
  std::vector<double> H(n_vars);
  int max_bins = 0;
  for (std::size_t v = 0; v < n_vars; ++v) {
    const int nb = bins_per_var[v];
    if (nb < 2) throw std::invalid_argument("each bin count must be >= 2");
    if (spline_order > nb)
      throw std::invalid_argument("spline_order must be <= every bin count");
    std::vector<double> dense(static_cast<std::size_t>(nb) * n_samples);
    gene_weights(data + v * n_samples, n_samples, spline_order, nb, dense.data());
    H[v] = marginal_entropy(dense.data(), n_samples, nb);
    W[v] = sparsify_weights(dense.data(), n_samples, nb, spline_order);
    if (nb > max_bins) max_bins = nb;
  }
  const std::size_t scratch = static_cast<std::size_t>(max_bins) * max_bins;
#ifdef _OPENMP
  const int nt = (n_threads > 0) ? n_threads : default_num_threads();
#else
  const int nt = 1;
#endif
  // Several bands per thread so the dynamic schedule can absorb the
  // remaining imbalance (the kernel cost varies with the bin counts).
  const std::vector<std::size_t> band = equal_area_bands(n_vars, 8 * static_cast<std::size_t>(nt));
  const long n_bands = static_cast<long>(band.size()) - 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
  {
    std::vector<double> joint(scratch);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (long b = 0; b < n_bands; ++b) {
      for (std::size_t i = band[b]; i < band[b + 1]; ++i) {
        T *row = out + packed_offset(i, n_vars);              // pair (i, j) at row[j - i - 1]
        for (std::size_t j = i + 1; j < n_vars; ++j)
          row[j - i - 1] = static_cast<T>(mi_pair_sparse(W[i], W[j], H[i], H[j], n_samples,
                                                 joint.data()));
      }
    }
  }
  (void)n_threads;
}

template <typename T>
void clr_calibrate_packed(T *packed, std::size_t n_vars, Combine combine) {
  if (n_vars < 2) throw std::invalid_argument("n_vars must be >= 2");
  const std::size_t n = n_vars;
  const double dn = static_cast<double>(n);
  // Pass 1: row means. Each gene's partners are visited in ascending index
  // order (rows i < j arrive first, then the gene's own row), which is the
  // dense row order; the diagonal adds 0.
  std::vector<double> mean(n, 0.0), var(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const T *row = packed + packed_offset(i, n);
    for (std::size_t j = i + 1; j < n; ++j) {
      const double x = static_cast<double>(row[j - i - 1]);
      mean[i] += x;
      mean[j] += x;
    }
  }
  for (std::size_t i = 0; i < n; ++i) mean[i] /= dn;
  // Pass 2: sample variances, same order; the diagonal term (0 - mean)^2 is
  // added at the gene's own position, after its lower partners.
  for (std::size_t i = 0; i < n; ++i) {
    var[i] += (0.0 - mean[i]) * (0.0 - mean[i]);
    const T *row = packed + packed_offset(i, n);
    for (std::size_t j = i + 1; j < n; ++j) {
      const double x = static_cast<double>(row[j - i - 1]);
      const double di = x - mean[i], dj = x - mean[j];
      var[i] += di * di;
      var[j] += dj * dj;
    }
  }
  std::vector<double> sd(n);
  for (std::size_t i = 0; i < n; ++i) sd[i] = std::sqrt(var[i] / (dn - 1.0));
  // Pass 3: z-scores both ways, clip, combine, overwrite in place.
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
  for (long ii = 0; ii < static_cast<long>(n); ++ii) {
    const std::size_t i = static_cast<std::size_t>(ii);
    T *row = packed + packed_offset(i, n);
    for (std::size_t j = i + 1; j < n; ++j) {
      const double x = static_cast<double>(row[j - i - 1]);
      double a = (sd[i] > 0.0) ? (x - mean[i]) / sd[i] : 0.0;   // DEVIATION as clr_calibrate
      double b = (sd[j] > 0.0) ? (x - mean[j]) / sd[j] : 0.0;
      if (a < 0.0) a = 0.0;
      if (b < 0.0) b = 0.0;
      row[j - i - 1] = static_cast<T>((combine == Combine::STOUFFER) ? (a + b) * inv_sqrt2
                                                             : std::sqrt(a * a + b * b));
    }
  }
}

template <typename T>
void packed_to_dense(const T *packed, std::size_t n, double *out) {
  for (std::size_t i = 0; i < n; ++i) {
    out[i * n + i] = 0.0;
    const T *row = packed + packed_offset(i, n);
    for (std::size_t j = i + 1; j < n; ++j) {
      const double x = static_cast<double>(row[j - i - 1]);
      out[i * n + j] = x;
      out[j * n + i] = x;
    }
  }
}

template <typename T>
void packed_block(const T *packed, std::size_t n, const std::size_t *idx,
                  std::size_t m, double *out) {
  for (std::size_t a = 0; a < m; ++a) {
    if (idx[a] >= n) throw std::invalid_argument("block index out of range");
    out[a * m + a] = 0.0;
    for (std::size_t b = a + 1; b < m; ++b) {
      const double x = (idx[a] == idx[b]) ? 0.0
                       : static_cast<double>(packed[packed_index(idx[a], idx[b], n)]);
      out[a * m + b] = x;
      out[b * m + a] = x;
    }
  }
}

template <typename T>
void packed_topk(const T *packed, std::size_t n, std::size_t k, int *idx_out,
                 double *val_out) {
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    std::vector<std::pair<double, int>> cand;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
    for (long ii = 0; ii < static_cast<long>(n); ++ii) {
      const std::size_t i = static_cast<std::size_t>(ii);
      cand.clear();
      for (std::size_t j = 0; j < n; ++j) {
        if (j == i) continue;
        const double x = static_cast<double>(packed[packed_index(i, j, n)]);
        if (x > 0.0) cand.emplace_back(x, static_cast<int>(j));
      }
      const std::size_t kk = std::min(k, cand.size());
      // larger value first; ties: lower index first
      auto cmp = [](const std::pair<double, int> &p, const std::pair<double, int> &q) {
        return (p.first > q.first) || (p.first == q.first && p.second < q.second);
      };
      std::partial_sort(cand.begin(), cand.begin() + static_cast<long>(kk), cand.end(), cmp);
      for (std::size_t r = 0; r < k; ++r) {
        idx_out[i * k + r] = (r < kk) ? cand[r].second : -1;
        val_out[i * k + r] = (r < kk) ? cand[r].first : 0.0;
      }
    }
  }
}

template void mi_packed<float>(const double *, std::size_t, std::size_t, const int *, int, float *, int);
template void mi_packed<double>(const double *, std::size_t, std::size_t, const int *, int, double *, int);
template void clr_calibrate_packed<float>(float *, std::size_t, Combine);
template void clr_calibrate_packed<double>(double *, std::size_t, Combine);
template void packed_to_dense<float>(const float *, std::size_t, double *);
template void packed_to_dense<double>(const double *, std::size_t, double *);
template void packed_block<float>(const float *, std::size_t, const std::size_t *, std::size_t, double *);
template void packed_block<double>(const double *, std::size_t, const std::size_t *, std::size_t, double *);
template void packed_topk<float>(const float *, std::size_t, std::size_t, int *, double *);
template void packed_topk<double>(const double *, std::size_t, std::size_t, int *, double *);

}  // namespace clr_core
