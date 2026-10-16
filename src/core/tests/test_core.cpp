// test_core.cpp — standalone verification of the pure C++ core (no R needed).
// Compile: g++ -std=c++17 -O2 -I.. test_core.cpp ../clr_core.cpp -o test_core
#include <algorithm>
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
        // Bit-identical at default flags (same per-cell summation order);
        // aggressive FP contraction (-ffp-contract=fast/-march=native) can
        // introduce ~1e-15 differences, hence the tolerance.
        if (std::fabs(dense - sparse) > 1e-12) ++mismatches;
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

  // --- packed (tiled) MI -> CLR: equal to the dense path -------------------
  {
    const std::size_t PG = 300, PE = 120;
    std::vector<double> pd(PG * PE);
    for (std::size_t e = 0; e < PE; ++e) {
      double f[5];
      for (int q = 0; q < 5; ++q) f[q] = gauss(rng);
      for (std::size_t g = 0; g < PG; ++g)
        pd[g * PE + e] = f[g % 5] * ((g % 7) < 3 ? 1.0 : 0.3) + gauss(rng);
    }
    std::vector<int> pb(PG, 7);
    std::vector<double> dmi(PG * PG), dclr(PG * PG);
    clr_core::mi_matrix(pd.data(), PG, PE, pb.data(), 3, dmi.data(), 4);
    clr_core::clr_calibrate(dmi.data(), PG, clr_core::Combine::STOUFFER, dclr.data());

    // bands: cover 0..n, increasing, roughly equal area
    const auto band = clr_core::equal_area_bands(PG, 16);
    CHECK(band.front() == 0 && band.back() == PG, "bands do not cover all rows");
    std::size_t maxa = 0, mina = clr_core::packed_size(PG);
    for (std::size_t b = 0; b + 1 < band.size(); ++b) {
      CHECK(band[b] < band[b + 1], "bands not strictly increasing");
      std::size_t area = clr_core::packed_offset(band[b + 1], PG) - clr_core::packed_offset(band[b], PG);
      maxa = std::max(maxa, area); mina = std::min(mina, area);
    }
    CHECK(maxa < 2 * mina + PG, "bands badly unbalanced");

    // double: bit-identical MI and CLR
    std::vector<double> pk(clr_core::packed_size(PG));
    clr_core::mi_packed<double>(pd.data(), PG, PE, pb.data(), 3, pk.data(), 4);
    bool same_mi = true;
    for (std::size_t i = 0; i < PG; ++i)
      for (std::size_t j = i + 1; j < PG; ++j)
        if (pk[clr_core::packed_index(i, j, PG)] != dmi[i * PG + j]) same_mi = false;
    CHECK(same_mi, "packed MI (double) differs from mi_matrix");
    clr_core::clr_calibrate_packed<double>(pk.data(), PG, clr_core::Combine::STOUFFER);
    std::vector<double> back(PG * PG);
    clr_core::packed_to_dense<double>(pk.data(), PG, back.data());
    bool same_clr = true;
    for (std::size_t q = 0; q < PG * PG; ++q) if (back[q] != dclr[q]) same_clr = false;
    CHECK(same_clr, "packed CLR (double) is not bit-identical to clr_calibrate");

    // float: close
    std::vector<float> pf(clr_core::packed_size(PG));
    clr_core::mi_packed<float>(pd.data(), PG, PE, pb.data(), 3, pf.data(), 4);
    clr_core::clr_calibrate_packed<float>(pf.data(), PG, clr_core::Combine::STOUFFER);
    double maxd = 0.0;
    for (std::size_t i = 0; i < PG; ++i)
      for (std::size_t j = i + 1; j < PG; ++j)
        maxd = std::max(maxd, std::fabs(static_cast<double>(pf[clr_core::packed_index(i, j, PG)]) - dclr[i * PG + j]));
    CHECK(maxd < 1e-4, "packed CLR (float) too far from double");

    // top-k against brute force on the dense scores (ties: lower index first)
    const std::size_t K = 20;
    std::vector<int> ti(PG * K); std::vector<double> tv(PG * K);
    clr_core::packed_topk<double>(pk.data(), PG, K, ti.data(), tv.data());
    bool topk_ok = true;
    for (std::size_t i = 0; i < PG; ++i) {
      std::vector<std::pair<double, int>> c;
      for (std::size_t j = 0; j < PG; ++j)
        if (j != i && dclr[i * PG + j] > 0) c.emplace_back(dclr[i * PG + j], static_cast<int>(j));
      std::stable_sort(c.begin(), c.end(), [](const std::pair<double, int> &p, const std::pair<double, int> &q) { return p.first > q.first; });
      for (std::size_t r = 0; r < K; ++r) {
        const int want = r < c.size() ? c[r].second : -1;
        if (ti[i * K + r] != want) topk_ok = false;
      }
    }
    CHECK(topk_ok, "packed_topk disagrees with brute force");

    // block
    const std::size_t idx[4] = {7, 3, 150, 299};
    std::vector<double> blk(16);
    clr_core::packed_block<double>(pk.data(), PG, idx, 4, blk.data());
    bool blk_ok = true;
    for (int a = 0; a < 4; ++a)
      for (int b = 0; b < 4; ++b)
        if (blk[a * 4 + b] != (a == b ? 0.0 : dclr[idx[a] * PG + idx[b]])) blk_ok = false;
    CHECK(blk_ok, "packed_block disagrees with dense");
    std::printf("packed path: MI identical %d, CLR identical %d, float max |diff| %.2e\n",
                static_cast<int>(same_mi), static_cast<int>(same_clr), maxd);
  }

  if (failures == 0) {
    std::printf("ALL CORE TESTS PASSED\n");
    return 0;
  }
  std::printf("%d FAILURES\n", failures);
  return 1;
}
