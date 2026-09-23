// test_core.cpp — standalone verification of the pure C++ core (no R needed).
// Compile: g++ -std=c++17 -O2 -I.. test_core.cpp ../clr_core.cpp -o test_core
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "../clr_core.hpp"

static int failures = 0;
#define CHECK(cond, msg)                                              \
  do {                                                                \
    if (!(cond)) {                                                    \
      std::printf("FAIL: %s\n", msg);                                 \
      ++failures;                                                     \
    }                                                                 \
  } while (0)

int main() {
  const std::size_t G = 6, E = 80;
  std::mt19937 rng(42);
  std::normal_distribution<double> gauss(0.0, 1.0);

  // Synthetic genes x samples, row-major: genes 0-1 correlated, 2-3
  // correlated, 4-5 independent noise.
  std::vector<double> data(G * E);
  for (std::size_t e = 0; e < E; ++e) {
    const double a = gauss(rng), b = gauss(rng);
    data[0 * E + e] = a + 0.2 * gauss(rng);
    data[1 * E + e] = a + 0.2 * gauss(rng);
    data[2 * E + e] = b + 0.2 * gauss(rng);
    data[3 * E + e] = b + 0.2 * gauss(rng);
    data[4 * E + e] = gauss(rng);
    data[5 * E + e] = gauss(rng);
  }

  // --- MI matrix, historical parity settings (10 bins, order 3) ---
  std::vector<int> bins(G, 10);
  std::vector<double> mi(G * G, -1.0);
  clr_core::mi_matrix(data.data(), G, E, bins.data(), 3, mi.data());

  for (std::size_t i = 0; i < G; ++i) {
    for (std::size_t j = 0; j < G; ++j) {
      CHECK(std::fabs(mi[i * G + j] - mi[j * G + i]) < 1e-12,
            "MI matrix not symmetric");
      CHECK(mi[i * G + j] > -1e-6, "MI unexpectedly negative");
    }
    CHECK(mi[i * G + i] > 0.0, "self-MI should be positive");
  }
  // Correlated pairs must carry more MI than noise pairs.
  CHECK(mi[0 * G + 1] > mi[0 * G + 4] + 0.1, "correlated MI < noise MI (0,1)");
  CHECK(mi[2 * G + 3] > mi[2 * G + 5] + 0.1, "correlated MI < noise MI (2,3)");

  // --- mi_pair seam consistency: driver == direct per-pair call ---
  {
    const int nb = 10;
    std::vector<double> w0(nb * E), w1(nb * E);
    clr_core::gene_weights(data.data() + 0 * E, E, 3, nb, w0.data());
    clr_core::gene_weights(data.data() + 1 * E, E, 3, nb, w1.data());
    const double h0 = clr_core::marginal_entropy(w0.data(), E, nb);
    const double h1 = clr_core::marginal_entropy(w1.data(), E, nb);
    const double direct = clr_core::mi_pair(w0.data(), w1.data(), h0, h1, E, nb, nb);
    CHECK(std::fabs(direct - mi[0 * G + 1]) < 1e-12,
          "mi_pair seam disagrees with mi_matrix driver");
  }

  // --- sparse kernel is BIT-IDENTICAL to the dense reference kernel ---
  {
    std::uniform_int_distribution<int> ub(3, 40);
    std::vector<double> scratch(40 * 40);
    int mismatches = 0;
    for (int order = 2; order <= 3; ++order) {
      for (int rep = 0; rep < 40; ++rep) {
        const int nbx = ub(rng), nby = ub(rng);
        const std::size_t gx = static_cast<std::size_t>(rep) % G;
        const std::size_t gy = (gx + 1 + static_cast<std::size_t>(rep) % (G - 1)) % G;
        std::vector<double> wx(nbx * E), wy(nby * E);
        clr_core::gene_weights(data.data() + gx * E, E, order, nbx, wx.data());
        clr_core::gene_weights(data.data() + gy * E, E, order, nby, wy.data());
        const double hx = clr_core::marginal_entropy(wx.data(), E, nbx);
        const double hy = clr_core::marginal_entropy(wy.data(), E, nby);
        const double dense = clr_core::mi_pair(wx.data(), wy.data(), hx, hy, E, nbx, nby);
        const auto sx = clr_core::sparsify_weights(wx.data(), E, nbx, order);
        const auto sy = clr_core::sparsify_weights(wy.data(), E, nby, order);
        const double sparse = clr_core::mi_pair_sparse(sx, sy, hx, hy, E, scratch.data());
        if (dense != sparse) ++mismatches;
      }
    }
    CHECK(mismatches == 0, "sparse MI kernel not bit-identical to dense kernel");
  }

  // --- adaptive bins: different per-gene counts still work ---
  {
    std::vector<int> abins = {6, 8, 10, 12, 15, 20};
    std::vector<double> mi2(G * G);
    clr_core::mi_matrix(data.data(), G, E, abins.data(), 3, mi2.data());
    for (std::size_t i = 0; i < G; ++i)
      for (std::size_t j = 0; j < G; ++j)
        CHECK(std::fabs(mi2[i * G + j] - mi2[j * G + i]) < 1e-12,
              "adaptive-bin MI not symmetric");
    CHECK(mi2[0 * G + 1] > mi2[0 * G + 4], "adaptive bins lost signal");
  }

  // --- CLR calibration: euclidean (historical) ---
  std::vector<double> se(G * G), ss(G * G);
  clr_core::clr_calibrate(mi.data(), G, clr_core::Combine::EUCLIDEAN, se.data());
  clr_core::clr_calibrate(mi.data(), G, clr_core::Combine::STOUFFER, ss.data());
  for (std::size_t i = 0; i < G; ++i) {
    CHECK(se[i * G + i] == 0.0, "euclidean diag not zero");
    CHECK(ss[i * G + i] == 0.0, "stouffer diag not zero");
    for (std::size_t j = 0; j < G; ++j) {
      CHECK(se[i * G + j] >= 0.0, "euclidean score negative");
      CHECK(ss[i * G + j] >= 0.0, "stouffer score negative");
      CHECK(std::fabs(se[i * G + j] - se[j * G + i]) < 1e-12,
            "euclidean not symmetric");
      CHECK(std::fabs(ss[i * G + j] - ss[j * G + i]) < 1e-12,
            "stouffer not symmetric");
    }
  }

  // --- Independent re-derivation of z-scores: verify combination formulas ---
  {
    std::vector<double> m(mi.begin(), mi.end());
    for (std::size_t i = 0; i < G; ++i) m[i * G + i] = 0.0;
    std::vector<std::vector<double>> z(G, std::vector<double>(G));
    for (std::size_t i = 0; i < G; ++i) {
      double mean = 0.0;
      for (std::size_t j = 0; j < G; ++j) mean += m[i * G + j];
      mean /= G;
      double var = 0.0;
      for (std::size_t j = 0; j < G; ++j) {
        const double d = m[i * G + j] - mean;
        var += d * d;
      }
      const double sd = std::sqrt(var / (G - 1.0));
      for (std::size_t j = 0; j < G; ++j) {
        double zij = (m[i * G + j] - mean) / sd;
        if (zij < 0.0) zij = 0.0;
        z[i][j] = zij;
      }
    }
    const double tol = 1e-9;
    for (std::size_t i = 0; i < G; ++i)
      for (std::size_t j = i + 1; j < G; ++j) {
        const double eu = std::sqrt(z[i][j] * z[i][j] + z[j][i] * z[j][i]);
        const double st = (z[i][j] + z[j][i]) / std::sqrt(2.0);
        char msg[128];
        std::snprintf(msg, sizeof(msg), "euclidean != sqrt(z^2+z'^2) at (%zu,%zu)", i, j);
        CHECK(std::fabs(se[i * G + j] - eu) < tol, msg);
        std::snprintf(msg, sizeof(msg), "stouffer != (z+z')/sqrt(2) at (%zu,%zu)", i, j);
        CHECK(std::fabs(ss[i * G + j] - st) < tol, msg);
      }
  }

  // --- Correlated pairs should top the CLR ranking ---
  {
    double top_corr = se[0 * G + 1] > se[2 * G + 3] ? se[2 * G + 3] : se[0 * G + 1];
    double top_noise = 0.0;
    for (std::size_t i = 0; i < G; ++i)
      for (std::size_t j = i + 1; j < G; ++j)
        if (!((i < 2 && j < 2) || (i >= 2 && i < 4 && j >= 2 && j < 4)))
          if (se[i * G + j] > top_noise) top_noise = se[i * G + j];
    CHECK(top_corr > top_noise, "CLR did not rank correlated pairs above noise");
  }

  if (failures == 0) {
    std::printf("ALL CORE TESTS PASSED\n");
    return 0;
  }
  std::printf("%d FAILURES\n", failures);
  return 1;
}
