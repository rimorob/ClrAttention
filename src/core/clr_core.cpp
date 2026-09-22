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
  const double scale =
      static_cast<double>(num_bins - spline_order + 1) / (xmax - xmin);
  for (std::size_t s = 0; s < n_samples; ++s) z[s] = (x[s] - xmin) * scale;
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

  // Precompute marginal weights + entropies once per gene (as miSubMarix did).
  std::vector<std::vector<double>> W(n_vars);
  std::vector<double> H(n_vars);
  for (std::size_t v = 0; v < n_vars; ++v) {
    const int nb = bins_per_var[v];
    if (nb < 2) throw std::invalid_argument("each bin count must be >= 2");
    W[v].resize(static_cast<std::size_t>(nb) * n_samples);
    gene_weights(data + v * n_samples, n_samples, spline_order, nb,
                 W[v].data());
    H[v] = marginal_entropy(W[v].data(), n_samples, nb);
  }

  // Triangular pair loop: iteration i owns pairs (i, j>=i) and writes each
  // entry exactly once, so iterations are independent. Dynamic scheduling
  // absorbs the triangular load imbalance. Pure C++ in here: no R/Python
  // API calls from worker threads.
#ifdef _OPENMP
  const int nt = (n_threads > 0) ? n_threads : default_num_threads();
#pragma omp parallel for schedule(dynamic) num_threads(nt)
#endif
  for (std::size_t i = 0; i < n_vars; ++i) {
    for (std::size_t j = i; j < n_vars; ++j) {
      const double m = mi_pair(W[i].data(), W[j].data(), H[i], H[j], n_samples,
                               bins_per_var[i], bins_per_var[j]);
      mi_out[i * n_vars + j] = m;
      mi_out[j * n_vars + i] = m;
    }
  }
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
