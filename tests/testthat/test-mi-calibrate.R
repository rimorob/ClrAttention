make_toy <- function() {
  set.seed(42)
  G <- 8; E <- 60
  d <- matrix(rnorm(G * E), nrow = G)
  # inject two correlated pairs
  sig <- rnorm(E)
  d[1, ] <- sig + 0.2 * rnorm(E)
  d[2, ] <- sig + 0.2 * rnorm(E)
  sig2 <- rnorm(E)
  d[3, ] <- sig2 + 0.2 * rnorm(E)
  d[4, ] <- sig2 + 0.2 * rnorm(E)
  d
}

test_that("bspline_mi returns a valid symmetric MI matrix", {
  d <- make_toy()
  mi <- bspline_mi(d, bins = 10, spline_order = 3, transform = "none")
  expect_true(is.matrix(mi) && nrow(mi) == 8 && ncol(mi) == 8)
  expect_equal(mi, t(mi), tolerance = 1e-12)
  expect_true(all(mi >= -1e-6))
  expect_true(all(diag(mi) > 0))  # self-MI, zeroed later by calibration
  expect_gt(mi[1, 2], mi[1, 5] + 0.1)
  expect_gt(mi[3, 4], mi[3, 6] + 0.1)
})

test_that("bspline_mi works with adaptive bins and propagates names", {
  d <- make_toy()
  rownames(d) <- paste0("g", 1:8)
  mi <- bspline_mi(d, bins = "fd")
  expect_equal(rownames(mi), paste0("g", 1:8))
  expect_equal(mi, t(mi), tolerance = 1e-12)
  mi2 <- bspline_mi(d, bins = "sturges", spline_order = 2)
  expect_equal(dim(mi2), c(8L, 8L))
})

test_that("clr_calibrate normal/euclidean matches clr.m semantics", {
  d <- make_toy()
  mi <- bspline_mi(d, bins = 10)
  A <- clr_calibrate(mi, method = "normal", combine = "euclidean")
  expect_equal(A, t(A), tolerance = 1e-12)
  expect_true(all(diag(A) == 0))
  expect_true(all(A >= 0))

  # independent re-derivation: zero diag, row z (n-1), clip, euclidean
  m <- mi; diag(m) <- 0
  z <- t(apply(m, 1, function(r) {
    s <- sd(r)
    zz <- if (s > 0) (r - mean(r)) / s else rep(0, length(r))
    pmax(zz, 0)
  }))
  ref <- sqrt(z^2 + t(z)^2); diag(ref) <- 0
  expect_equal(A, ref, tolerance = 1e-9)
  expect_equal(rownames(A), rownames(mi))
})

test_that("clr_calibrate stouffer equals (z+z')/sqrt(2) after the same clip", {
  d <- make_toy()
  mi <- bspline_mi(d, bins = 10)
  S <- clr_calibrate(mi, method = "normal", combine = "stouffer")
  m <- mi; diag(m) <- 0
  z <- t(apply(m, 1, function(r) {
    s <- sd(r)
    zz <- if (s > 0) (r - mean(r)) / s else rep(0, length(r))
    pmax(zz, 0)
  }))
  ref <- (z + t(z)) / sqrt(2); diag(ref) <- 0
  expect_equal(S, ref, tolerance = 1e-9)
  expect_true(all(S >= 0))
  expect_equal(S, t(S), tolerance = 1e-12)
})

test_that("rayleigh and kde variants run with correct shape and combination", {
  d <- make_toy()
  mi <- bspline_mi(d, bins = 10)
  R <- clr_calibrate(mi, method = "rayleigh")
  expect_equal(dim(R), c(8L, 8L))
  expect_true(all(diag(R) == 0))
  expect_true(all(R >= 0 & R <= 1 + 1e-9))
  # combination identity A + A' - A*A' recomputed from marginal CDFs
  G <- 8
  C <- matrix(0, G, G)
  for (i in seq_len(G)) {
    sg <- sqrt(mean(mi[i, ]^2) / 2)
    C[i, ] <- 1 - exp(-mi[i, ]^2 / (2 * sg^2))
  }
  expect_equal(R, C + t(C) - C * t(C) - diag(diag(C + t(C) - C * t(C))),
               tolerance = 1e-9, ignore_attr = TRUE)
  # (the trailing -diag(diag(...)) mirrors clr.m's final A - diag(diag(A)))

  K <- clr_calibrate(mi, method = "kde")
  expect_equal(dim(K), c(8L, 8L))
  expect_true(all(diag(K) == 0))
  expect_true(all(K <= 0))  # log-CDF combination is non-positive

  expect_error(clr_calibrate(mi, method = "kde", combine = "stouffer"),
               "only meaningful")
})

test_that("calibration validates inputs", {
  expect_error(clr_calibrate(matrix(1, 3, 4)), "square")
  expect_error(bspline_mi(matrix(rnorm(6), nrow = 1)), "at least 2 genes")
  bad <- matrix(rnorm(16), 4); bad[1, 1] <- NA
  expect_error(bspline_mi(bad), "NA/NaN/Inf")
})

test_that("threaded MI matches serial MI exactly", {
  d <- make_toy()
  mi1 <- bspline_mi(d, bins = 10, threads = 1)
  mi2 <- bspline_mi(d, bins = 10, threads = 2)
  mi_def <- bspline_mi(d, bins = 10)  # default: cores - 2
  expect_equal(mi2, mi1, tolerance = 1e-12)
  expect_equal(mi_def, mi1, tolerance = 1e-12)
  expect_error(bspline_mi(d, bins = 10, threads = 0), "positive integer")
  expect_error(bspline_mi(d, bins = 10, threads = -3), "positive integer")
})
