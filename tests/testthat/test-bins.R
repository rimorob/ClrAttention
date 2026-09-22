test_that("fd_bins behaves sensibly on known distributions", {
  set.seed(1)
  x <- rnorm(1000)
  nb <- fd_bins(x, rule = "fd")
  expect_true(nb >= 5 && nb <= 50)
  # Freedman-Diaconis on N(0,1), n=1000: h ~ 2*1.35*0.1 = 0.27, range ~ 6 -> ~22
  expect_true(nb >= 10 && nb <= 40)

  # constant vector -> degenerate -> Sturges fallback, then clamped
  nb_const <- fd_bins(rep(3, 200), rule = "fd")
  expect_equal(nb_const, as.integer(min(max(ceiling(log2(200)) + 1, 5), 50)))

  # scott and sturges run and respect caps
  expect_true(fd_bins(x, rule = "scott") >= 5)
  expect_equal(fd_bins(rnorm(50), rule = "sturges"), as.integer(ceiling(log2(50)) + 1))

  # caps are honored
  expect_equal(fd_bins(rnorm(1e5), rule = "fd", max_bins = 12), 12L)
  expect_error(fd_bins(1, rule = "fd"), "at least 2 finite values")
  expect_error(fd_bins(c(NA_real_, NaN), rule = "fd"), "at least 2 finite values")
})

test_that("bins_for_genes handles rules and integers", {
  set.seed(2)
  d <- matrix(rnorm(4 * 60), nrow = 4)
  b <- bins_for_genes(d, bins = "fd")
  expect_length(b, 4)
  expect_true(all(b >= 5 & b <= 50))
  expect_equal(bins_for_genes(d, bins = 10), rep(10L, 4))
  expect_error(bins_for_genes(d, bins = 1), ">= 2")
  expect_error(bins_for_genes(d, bins = "nope"), "match.arg|should be one of")
})
