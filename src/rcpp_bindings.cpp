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

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_openmp_info() {
  return Rcpp::IntegerVector::create(
      Rcpp::Named("openmp_max_threads") = clr_core::openmp_max_threads(),
      Rcpp::Named("default_threads") = clr_core::default_num_threads());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_weights_at(Rcpp::NumericVector v, double xmin,
                                   double xmax, int spline_order,
                                   int num_bins) {
  const std::size_t n = static_cast<std::size_t>(v.size());
  // core layout [bin][n] row-major == R column-major n x num_bins
  Rcpp::NumericMatrix out(static_cast<int>(n), num_bins);
  clr_core::weights_at(v.begin(), n, xmin, xmax, spline_order, num_bins,
                       out.begin());
  return out;
}

// ---------------------------------------------------------------------------
// Packed (tiled) MI -> CLR. The packed float32 scores stay in C++ behind an
// external pointer, so R never holds a G x G matrix unless it asks for one.
// ---------------------------------------------------------------------------
struct PackedClr {
  std::size_t n;
  std::vector<float> v;
};

static std::vector<double> to_row_major(const Rcpp::NumericMatrix &data) {
  const int G = data.nrow(), E = data.ncol();
  std::vector<double> buf(static_cast<std::size_t>(G) * E);
  for (int i = 0; i < G; ++i)
    for (int e = 0; e < E; ++e) buf[static_cast<std::size_t>(i) * E + e] = data(i, e);
  return buf;
}

// [[Rcpp::export]]
SEXP cpp_clr_packed(Rcpp::NumericMatrix data, Rcpp::IntegerVector bins_per_var,
                    int spline_order, int combine, int n_threads) {
  const int G = data.nrow(), E = data.ncol();
  if (bins_per_var.size() != G) Rcpp::stop("bins_per_var length != number of genes");
  if (combine != 0 && combine != 1) Rcpp::stop("combine must be 0 (euclidean) or 1 (stouffer)");
  std::vector<double> buf = to_row_major(data);
  std::vector<int> bins(bins_per_var.begin(), bins_per_var.end());
  PackedClr *p = new PackedClr{static_cast<std::size_t>(G),
                               std::vector<float>(clr_core::packed_size(static_cast<std::size_t>(G)))};
  try {
    clr_core::mi_packed<float>(buf.data(), p->n, static_cast<std::size_t>(E), bins.data(),
                               spline_order, p->v.data(), n_threads);
    clr_core::clr_calibrate_packed<float>(p->v.data(), p->n,
                                          static_cast<clr_core::Combine>(combine));
  } catch (...) { delete p; throw; }
  return Rcpp::XPtr<PackedClr>(p, true);
}

static PackedClr *packed_ptr(SEXP ptr) {
  Rcpp::XPtr<PackedClr> xp(ptr);
  if (xp.get() == nullptr) Rcpp::stop("invalid (released or unserialized) clr_packed pointer");
  return xp.get();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_packed_dense(SEXP ptr) {
  PackedClr *p = packed_ptr(ptr);
  const int G = static_cast<int>(p->n);
  std::vector<double> out(p->n * p->n);
  clr_core::packed_to_dense<float>(p->v.data(), p->n, out.data());
  Rcpp::NumericMatrix res(G, G);
  for (int i = 0; i < G; ++i)
    for (int j = 0; j < G; ++j) res(i, j) = out[static_cast<std::size_t>(i) * G + j];
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_packed_block(SEXP ptr, Rcpp::IntegerVector idx) {
  PackedClr *p = packed_ptr(ptr);
  const std::size_t m = static_cast<std::size_t>(idx.size());
  std::vector<std::size_t> ix(m);
  for (std::size_t a = 0; a < m; ++a) {
    if (idx[a] == NA_INTEGER || idx[a] < 1 || static_cast<std::size_t>(idx[a]) > p->n)
      Rcpp::stop("block index out of range");
    ix[a] = static_cast<std::size_t>(idx[a] - 1);
  }
  std::vector<double> out(m * m);
  clr_core::packed_block<float>(p->v.data(), p->n, ix.data(), m, out.data());
  Rcpp::NumericMatrix res(static_cast<int>(m), static_cast<int>(m));
  for (std::size_t a = 0; a < m; ++a)
    for (std::size_t b = 0; b < m; ++b) res(static_cast<int>(a), static_cast<int>(b)) = out[a * m + b];
  return res;
}

// [[Rcpp::export]]
Rcpp::List cpp_packed_topk(SEXP ptr, int k) {
  PackedClr *p = packed_ptr(ptr);
  if (k < 1) Rcpp::stop("k must be >= 1");
  const std::size_t K = static_cast<std::size_t>(k);
  std::vector<int> ix(p->n * K);
  std::vector<double> vl(p->n * K);
  clr_core::packed_topk<float>(p->v.data(), p->n, K, ix.data(), vl.data());
  const int G = static_cast<int>(p->n);
  Rcpp::IntegerMatrix I(G, k);
  Rcpp::NumericMatrix V(G, k);
  for (int i = 0; i < G; ++i)
    for (int r = 0; r < k; ++r) {
      const int j = ix[static_cast<std::size_t>(i) * K + r];
      I(i, r) = (j < 0) ? NA_INTEGER : j + 1;                  // 1-based for R
      V(i, r) = vl[static_cast<std::size_t>(i) * K + r];
    }
  return Rcpp::List::create(Rcpp::Named("idx") = I, Rcpp::Named("val") = V);
}

// Exactness check: the packed path in double precision, expanded to dense.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_clr_packed_double_dense(Rcpp::NumericMatrix data,
                                                Rcpp::IntegerVector bins_per_var,
                                                int spline_order, int combine, int n_threads) {
  const int G = data.nrow(), E = data.ncol();
  std::vector<double> buf = to_row_major(data);
  std::vector<int> bins(bins_per_var.begin(), bins_per_var.end());
  std::vector<double> pk(clr_core::packed_size(static_cast<std::size_t>(G)));
  clr_core::mi_packed<double>(buf.data(), static_cast<std::size_t>(G), static_cast<std::size_t>(E),
                              bins.data(), spline_order, pk.data(), n_threads);
  clr_core::clr_calibrate_packed<double>(pk.data(), static_cast<std::size_t>(G),
                                         static_cast<clr_core::Combine>(combine));
  std::vector<double> out(static_cast<std::size_t>(G) * G);
  clr_core::packed_to_dense<double>(pk.data(), static_cast<std::size_t>(G), out.data());
  Rcpp::NumericMatrix res(G, G);
  for (int i = 0; i < G; ++i)
    for (int j = 0; j < G; ++j) res(i, j) = out[static_cast<std::size_t>(i) * G + j];
  return res;
}
