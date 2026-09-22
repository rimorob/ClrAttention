// rcpp_bindings.cpp — THIN Rcpp wrappers. No math lives here: convert R's
// column-major matrices to the core's row-major buffers at the boundary,
// call clr_core, convert back. All numerics are in src/core/.

#include <Rcpp.h>

#include "core/clr_core.hpp"

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_mi_matrix(Rcpp::NumericMatrix data,
                                  Rcpp::IntegerVector bins_per_var,
                                  int spline_order,
                                  int n_threads) {
  const int G = data.nrow();
  const int E = data.ncol();
  if (bins_per_var.size() != G)
    Rcpp::stop("bins_per_var length (%d) != number of genes (%d)",
               bins_per_var.size(), G);

  // R column-major -> core row-major (genes are rows in both).
  std::vector<double> buf(static_cast<std::size_t>(G) * E);
  for (int i = 0; i < G; ++i)
    for (int e = 0; e < E; ++e)
      buf[static_cast<std::size_t>(i) * E + e] = data(i, e);

  std::vector<int> bins(G);
  for (int i = 0; i < G; ++i) bins[static_cast<std::size_t>(i)] = bins_per_var[i];

  std::vector<double> out(static_cast<std::size_t>(G) * G);
  clr_core::mi_matrix(buf.data(), static_cast<std::size_t>(G),
                      static_cast<std::size_t>(E), bins.data(), spline_order,
                      out.data(), n_threads);

  Rcpp::NumericMatrix res(G, G);
  for (int i = 0; i < G; ++i)
    for (int j = 0; j < G; ++j)
      res(i, j) = out[static_cast<std::size_t>(i) * G + j];
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_clr_calibrate(Rcpp::NumericMatrix mi, int combine) {
  const int G = mi.nrow();
  if (mi.ncol() != G) Rcpp::stop("mi must be square");
  if (combine != 0 && combine != 1)
    Rcpp::stop("combine must be 0 (euclidean) or 1 (stouffer)");

  std::vector<double> buf(static_cast<std::size_t>(G) * G);
  for (int i = 0; i < G; ++i)
    for (int j = 0; j < G; ++j)
      buf[static_cast<std::size_t>(i) * G + j] = mi(i, j);

  std::vector<double> out(static_cast<std::size_t>(G) * G);
  clr_core::clr_calibrate(buf.data(), static_cast<std::size_t>(G),
                         static_cast<clr_core::Combine>(combine), out.data());

  Rcpp::NumericMatrix res(G, G);
  for (int i = 0; i < G; ++i)
    for (int j = 0; j < G; ++j)
      res(i, j) = out[static_cast<std::size_t>(i) * G + j];
  return res;
}
