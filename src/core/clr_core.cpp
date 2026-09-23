// clr_core.cpp — Pure C++17 implementation. See clr_core.hpp for the
// multi-language contract. Ported from InfoKit2.c (functions SplineBlend,
// SplineKnots, xToZ, findWeights, hist1d/hist2d, entropy1d/entropy2d,
// miSubMarix) and from clr.m's 'normal' calibration path.

#include "clr_core.hpp"

#include <cmath>
#include <stdexcept>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace clr_core {

double log2d(double x) { return std::log(x) / std::log(2.0); }

std::vector<int> spline_knots(int num_bins, int spline_order) {
  if (num_bins < 2) throw std::invalid_argument("num_bins must be >= 2");
  if (spline_order < 1 || spline_order > num_bins)
    throw std::invalid_argument("spline_order must be in [1, num_bins]");
  const int d = num_bins - 1;
  std::vector<int> u(static_cast<std::size_t>(num_bins + spline_order));
  for (int j = 0; j <= d + spline_order; ++j) {
    if (j < spline_order)
      u[static_cast<std::size_t>(j)] = 0;
    else if (j <= d)
      u[static_cast<std::size_t>(j)] = u[static_cast<std::size_t>(j - 1)] + 1;
    else
      u[static_cast<std::size_t>(j)] = u[static_cast<std::size_t>(d)] + 1;
  }
  return u;
}

double spline_blend(int k, int t, const int *u, double v, int n) {
  double value = 0.0;
  if (t == 1) {
    if ((u[k] <= v && v < u[k + 1]) ||
        (std::fabs(v - u[k + 1]) < 1e-10 && (k + 1 == n)))
      value = 1.0;
    else
      value = 0.0;
  } else {
    const double d1 = static_cast<double>(u[k + t - 1] - u[k]);
    const double d2 = static_cast<double>(u[k + t] - u[k + 1]);
    if (d1 == 0.0 && d2 == 0.0) {
      value = 0.0;
    } else if (d1 == 0.0) {
      value = (u[k + t] - v) / d2 * spline_blend(k + 1, t - 1, u, v, n);
    } else if (d2 == 0.0) {
      value = (v - u[k]) / d1 * spline_blend(k, t - 1, u, v, n);
    } else {
      value = (v - u[k]) / d1 * spline_blend(k, t - 1, u, v, n) +
              (u[k + t] - v) / d2 * spline_blend(k + 1, t - 1, u, v, n);
    }
  }
  if (value < 0.0) value = 0.0;  // rounding sometimes makes this < 0
  return value;
}

void x_to_z(const double *x, double *z, std::size_t n_samples, int spline_order,
            int num_bins) {
  double xmin = x[0], xmax = x[0];
  for (std::size_t s = 1; s < n_samples; ++s) {
    if (x[s] < xmin) xmin = x[s];
    if (x[s] > xmax) xmax = x[s];
  }
  // DEVIATION: the original divided by (xmax - xmin) unguarded; a constant
  // gene then produced NaN z, all-zero weights, H(X) = 0 and MI(X, Y) = H(Y)
  // -- i.e. the constant gene became the strongest hub. Refuse instead.
  if (!(xmax > xmin))
    throw std::invalid_argument("constant variable: max == min (no information)");
  const double scale =
      static_cast<double>(num_bins - spline_order + 1) / (xmax - xmin);
  for (std::size_t s = 0; s < n_samples; ++s) z[s] = (x[s] - xmin) * scale;
}

void weights_at(const double *v, std::size_t n, double xmin, double xmax,
                int spline_order, int num_bins, double *weights_out) {
  if (!(xmax > xmin))
    throw std::invalid_argument("weights_at: need xmax > xmin");
  const std::vector<int> knots = spline_knots(num_bins, spline_order);
  const double top = static_cast<double>(num_bins - spline_order + 1);
  const double scale = top / (xmax - xmin);
  for (std::size_t s = 0; s < n; ++s) {
    double x = v[s];
    if (x < xmin) x = xmin;          // clamp: values outside the fitted
    if (x > xmax) x = xmax;          // range take the boundary basis
    const double z = (x - xmin) * scale;
    for (int b = 0; b < num_bins; ++b)
      weights_out[static_cast<std::size_t>(b) * n + s] =
          spline_blend(b, spline_order, knots.data(), z, num_bins);
  }
}

void gene_weights(const double *x, std::size_t n_samples, int spline_order,
                  int num_bins, double *weights_out) {
  const std::vector<int> knots = spline_knots(num_bins, spline_order);
  std::vector<double> z(n_samples);
  x_to_z(x, z.data(), n_samples, spline_order, num_bins);
  for (std::size_t s = 0; s < n_samples; ++s) {
    for (int b = 0; b < num_bins; ++b) {
      weights_out[static_cast<std::size_t>(b) * n_samples + s] =
          spline_blend(b, spline_order, knots.data(), z[s], num_bins);
    }
  }
}

double marginal_entropy(const double *weights, std::size_t n_samples,
                        int num_bins) {
  double H = 0.0;
  const double n = static_cast<double>(n_samples);
  for (int b = 0; b < num_bins; ++b) {
    double h = 0.0;
    const double *row = weights + static_cast<std::size_t>(b) * n_samples;
    for (std::size_t s = 0; s < n_samples; ++s) h += row[s];
    h /= n;
    if (h > 0.0) H -= h * log2d(h);
  }
  return H;
}

double mi_pair(const double *wx, const double *wy, double hx, double hy,
               std::size_t n_samples, int nbx, int nby) {
  double H = 0.0;  // joint entropy H(X, Y) in bits
  const double n = static_cast<double>(n_samples);
  for (int bx = 0; bx < nbx; ++bx) {
    const double *rowx = wx + static_cast<std::size_t>(bx) * n_samples;
    for (int by = 0; by < nby; ++by) {
      const double *rowy = wy + static_cast<std::size_t>(by) * n_samples;
      double h = 0.0;
      for (std::size_t s = 0; s < n_samples; ++s) h += rowx[s] * rowy[s];
      h /= n;
      if (h > 0.0) H -= h * log2d(h);
    }
  }
  return hx + hy - H;
}

SparseWeights sparsify_weights(const double *weights, std::size_t n_samples,
                               int num_bins, int order) {
  if (order < 1 || order > num_bins)
    throw std::invalid_argument("order must be in [1, num_bins]");
  SparseWeights sw;
  sw.num_bins = num_bins;
  sw.order = order;
  sw.first.assign(n_samples, 0);
  sw.vals.assign(n_samples * static_cast<std::size_t>(order), 0.0);
  for (std::size_t s = 0; s < n_samples; ++s) {
    int lo = -1, hi = -1;
    for (int b = 0; b < num_bins; ++b) {
      if (weights[static_cast<std::size_t>(b) * n_samples + s] != 0.0) {
        if (lo < 0) lo = b;
        hi = b;
      }
    }
    int f = (lo < 0) ? 0 : lo;
    if (f + order > num_bins) f = num_bins - order;
    if (lo >= 0 && hi >= f + order)
      throw std::logic_error("B-spline support wider than spline order");
    sw.first[s] = f;
    for (int a = 0; a < order; ++a)
      sw.vals[s * static_cast<std::size_t>(order) + static_cast<std::size_t>(a)] =
          weights[static_cast<std::size_t>(f + a) * n_samples + s];
  }
  return sw;
}

double mi_pair_sparse(const SparseWeights &wx, const SparseWeights &wy,
                      double hx, double hy, std::size_t n_samples,
                      double *joint) {
  const int nbx = wx.num_bins, nby = wy.num_bins;
  const int kx = wx.order, ky = wy.order;
  const std::size_t cells = static_cast<std::size_t>(nbx) * nby;
  for (std::size_t c = 0; c < cells; ++c) joint[c] = 0.0;
  for (std::size_t s = 0; s < n_samples; ++s) {
    const double *vx = wx.vals.data() + s * static_cast<std::size_t>(kx);
    const double *vy = wy.vals.data() + s * static_cast<std::size_t>(ky);
    const int fx = wx.first[s], fy = wy.first[s];
    for (int a = 0; a < kx; ++a) {
      double *row = joint + static_cast<std::size_t>(fx + a) * nby + fy;
      for (int b = 0; b < ky; ++b) row[b] += vx[a] * vy[b];
    }
  }
  double H = 0.0;
  const double n = static_cast<double>(n_samples);
  for (std::size_t c = 0; c < cells; ++c) {
    const double h = joint[c] / n;
    if (h > 0.0) H -= h * log2d(h);
  }
  return hx + hy - H;
}

int openmp_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 0;
#endif
}

int default_num_threads() {
  const unsigned hc = std::thread::hardware_concurrency();
  if (hc <= 3) return 1;  // tiny machines (or unknown): stay serial
  return static_cast<int>(hc) - 2;
}

void mi_matrix(const double *data, std::size_t n_vars, std::size_t n_samples,
               const int *bins_per_var, int spline_order, double *mi_out,
               int n_threads) {
  if (n_vars < 2) throw std::invalid_argument("n_vars must be >= 2");
  if (n_samples < 2) throw std::invalid_argument("n_samples must be >= 2");

  // Precompute marginal weights + entropies once per gene (as miSubMarix
  // did), then keep only the sparse form: memory is O(G * N * order)
  // instead of O(G * N * bins).
  std::vector<SparseWeights> W(n_vars);
  std::vector<double> H(n_vars);
  int max_bins = 0;
  for (std::size_t v = 0; v < n_vars; ++v) {
    const int nb = bins_per_var[v];
    if (nb < 2) throw std::invalid_argument("each bin count must be >= 2");
    if (spline_order > nb)
      throw std::invalid_argument("spline_order must be <= every bin count");
    std::vector<double> dense(static_cast<std::size_t>(nb) * n_samples);
    gene_weights(data + v * n_samples, n_samples, spline_order, nb,
                 dense.data());
    H[v] = marginal_entropy(dense.data(), n_samples, nb);
    W[v] = sparsify_weights(dense.data(), n_samples, nb, spline_order);
    if (nb > max_bins) max_bins = nb;
  }
  const std::size_t scratch = static_cast<std::size_t>(max_bins) * max_bins;

  // Triangular pair loop: iteration i owns pairs (i, j>=i) and writes each
  // entry exactly once, so iterations are independent. Dynamic scheduling
  // absorbs the triangular load imbalance. Pure C++ in here: no R/Python
  // API calls from worker threads. Each thread owns its joint scratch.
#ifdef _OPENMP
  const int nt = (n_threads > 0) ? n_threads : default_num_threads();
#pragma omp parallel num_threads(nt)
#endif
  {
    std::vector<double> joint(scratch);
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (std::size_t i = 0; i < n_vars; ++i) {
      for (std::size_t j = i; j < n_vars; ++j) {
        const double m = mi_pair_sparse(W[i], W[j], H[i], H[j], n_samples,
                                        joint.data());
        mi_out[i * n_vars + j] = m;
        mi_out[j * n_vars + i] = m;
      }
    }
  }
  (void)n_threads;
}

void clr_calibrate(const double *mi, std::size_t n_vars, Combine combine,
                   double *scores_out) {
  if (n_vars < 2) throw std::invalid_argument("n_vars must be >= 2");

  // Work on a copy with the diagonal zeroed (clr.m: MI - diag(diag(MI))).
  std::vector<double> m(mi, mi + n_vars * n_vars);
  for (std::size_t i = 0; i < n_vars; ++i) m[i * n_vars + i] = 0.0;

  // Row-wise z-scores, sample std (n - 1), negatives clipped to 0.
  std::vector<double> z(n_vars * n_vars);
  const double dn = static_cast<double>(n_vars);
  for (std::size_t i = 0; i < n_vars; ++i) {
    const double *row = m.data() + i * n_vars;
    double mean = 0.0;
    for (std::size_t j = 0; j < n_vars; ++j) mean += row[j];
    mean /= dn;
    double var = 0.0;
    for (std::size_t j = 0; j < n_vars; ++j) {
      const double d = row[j] - mean;
      var += d * d;
    }
    const double sd = std::sqrt(var / (dn - 1.0));
    for (std::size_t j = 0; j < n_vars; ++j) {
      double zij = (sd > 0.0) ? (row[j] - mean) / sd : 0.0;  // DEVIATION: see header
      if (zij < 0.0) zij = 0.0;
      z[i * n_vars + j] = zij;
    }
  }

  // Bilateral combination, then zero the diagonal.
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  for (std::size_t i = 0; i < n_vars; ++i) {
    scores_out[i * n_vars + i] = 0.0;
    for (std::size_t j = i + 1; j < n_vars; ++j) {
      const double a = z[i * n_vars + j];
      const double b = z[j * n_vars + i];
      const double s = (combine == Combine::STOUFFER)
                           ? (a + b) * inv_sqrt2
                           : std::sqrt(a * a + b * b);
      scores_out[i * n_vars + j] = s;
      scores_out[j * n_vars + i] = s;
    }
  }
}

}  // namespace clr_core
