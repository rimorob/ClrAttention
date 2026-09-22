make_toy <- function() {
  set.seed(7)
  G <- 10; E <- 50
  d <- matrix(rnorm(G * E), nrow = G)
  rownames(d) <- paste0("g", seq_len(G))
  d
}

test_that("ClrAttention validates input in initialize()", {
  expect_error(ClrAttention$new("nope"), "numeric matrix")
  expect_error(ClrAttention$new(matrix(rnorm(6), nrow = 1)), "at least 2 genes")
  expect_error(ClrAttention$new(matrix(rnorm(6), ncol = 1)), "at least 2 samples")
  bad <- matrix(rnorm(12), 3, 4); bad[2, 2] <- Inf
  expect_error(ClrAttention$new(bad), "NA/NaN/Inf")
  fit <- ClrAttention$new(make_toy())
  expect_s3_class(fit, "R6")
  expect_equal(class(fit)[1], "ClrAttention")
})

test_that("active bindings are strict: error before the stage runs", {
  fit <- ClrAttention$new(make_toy())
  expect_error(fit$mi, "estimate_mi")
  expect_error(fit$clr_scores, "calibrate")
  expect_error(fit$operator, "build_operator")
  expect_error(fit$embedding, "diffuse")
  expect_error(fit$calibrate(), "estimate_mi")
  expect_error(fit$build_operator(), "calibrate")
  expect_error(fit$diffuse(), "build_operator")
})

test_that("pipeline chains and each stage populates its binding", {
  fit <- ClrAttention$new(make_toy())
  out <- fit$estimate_mi(bins = 10)$calibrate()$build_operator(k = 3)$diffuse(steps = 5)
  expect_true(identical(out, fit))  # invisible(self) chaining

  G <- 10
  expect_equal(dim(fit$mi), c(G, G))
  expect_equal(dim(fit$clr_scores), c(G, G))
  expect_true(all(diag(fit$clr_scores) == 0))
  A <- fit$operator
  expect_equal(dim(A), c(G, G))
  expect_equal(unname(rowSums(A)), rep(1, G), tolerance = 1e-12)
  # top-k sparsity: at most k nonzeros per row
  expect_true(all(rowSums(A > 0) <= 3))
  expect_equal(dim(fit$embedding), c(G, 50))
  expect_length(fit$trajectory, 6)  # E^(0..5)

  p <- fit$params
  expect_equal(p$n_genes, G)
  expect_equal(p$diffuse$steps, 5L)
  expect_equal(p$operator$alpha, 0.5)
})

test_that("trajectory starts at the raw data", {
  d <- make_toy()
  fit <- ClrAttention$new(d)$estimate_mi(bins = 10)$calibrate()$
    build_operator(k = 4)$diffuse(steps = 3)
  expect_equal(unname(fit$trajectory[[1]]), unname(d), ignore_attr = TRUE)
})

test_that("diffusion converges toward a stationary embedding on a toy", {
  set.seed(11)
  d <- matrix(rnorm(6 * 40), nrow = 6)
  fit <- ClrAttention$new(d)$estimate_mi(bins = 10)$calibrate()$
    build_operator(k = 5, alpha = 0.9)$diffuse(steps = 200)
  tr <- fit$trajectory
  gap_late <- max(abs(tr[[201]] - tr[[200]]))
  gap_early <- max(abs(tr[[2]] - tr[[1]]))
  expect_lt(gap_late, gap_early)  # increments shrink: converging, not diverging
  expect_true(all(is.finite(fit$embedding)))
})

test_that("build_operator validates and supports tau", {
  fit <- ClrAttention$new(make_toy())$estimate_mi(bins = 10)$calibrate()
  expect_error(fit$build_operator(k = 0), "positive integer")
  expect_error(fit$build_operator(alpha = 0), "in \\(0, 1\\]")
  expect_error(fit$build_operator(alpha = 1.5), "in \\(0, 1\\]")
  expect_error(fit$build_operator(tau = -1), "non-negative")
  fit$build_operator(tau = 1.0)
  A <- fit$operator
  expect_true(all(A[A > 0] <= 1.0 + 1e-12))  # row-stochastic: no entry exceeds 1
  rs <- rowSums(A)
  expect_equal(unname(rs[rs > 0]), rep(1, sum(rs > 0)), tolerance = 1e-12)
})

test_that("re-running a stage invalidates downstream stages", {
  fit <- ClrAttention$new(make_toy())$estimate_mi(bins = 10)$calibrate()$
    build_operator()$diffuse(steps = 2)
  expect_silent(fit$embedding)
  fit$estimate_mi(bins = 12)
  expect_error(fit$clr_scores, "calibrate")
  expect_error(fit$embedding, "diffuse")
})

test_that("reset_diffusion clears only the trajectory", {
  fit <- ClrAttention$new(make_toy())$estimate_mi(bins = 10)$calibrate()$
    build_operator()$diffuse(steps = 2)
  fit$reset_diffusion()
  expect_error(fit$embedding, "diffuse")
  expect_silent(fit$operator)
  expect_null(fit$params$diffuse)
})

test_that("print() summarizes state", {
  fit <- ClrAttention$new(make_toy())
  expect_output(fit$print(), "ClrAttention.*10 genes x 50 samples")
  fit$estimate_mi(bins = 10)
  expect_output(fit$print(), "MI:\\s+estimated")
})
