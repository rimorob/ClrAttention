# Synthetic benchmark: 12 genes with a planted regulatory network.
#
# Enshrined regression test for the CLR-as-attention work. Compares
# historical-parity CLR (fixed 10 bins, Euclidean combine -- the 2008
# algorithm) against the new pipeline (adaptive FD bins, Stouffer combine,
# attention operator + diffusion), on data where ground truth is known.
#
# Design (genes x samples, 12 x 300):
#   g1 = TF_A drives g2 (quadratic), g3 (sinusoidal), g4 (saturating)
#   g5 = TF_B drives g6 (linear +), g7 (linear -)
#   g8..g12 = independent null genes
# True undirected edges: (1,2) (1,3) (1,4) (5,6) (5,7) -- 5 of 66 pairs.
# The quadratic edge has ~zero Pearson correlation, so this also guards
# the MI estimator's nonlinearity advantage.
#
# threads = 1 is deliberate: OpenMP reduction order is not bit-identical
# across thread counts, and a unit test must be deterministic.

synthetic_clr_data <- function(seed = 7, n = 300) {
  set.seed(seed)
  tfa <- rnorm(n)
  tfb <- rnorm(n)
  X <- rbind(
    tfa,
    tfa^2 + rnorm(n, sd = 0.7),
    sin(2 * tfa) + rnorm(n, sd = 0.7),
    tanh(3 * tfa) + rnorm(n, sd = 0.7),
    tfb,
    2 * tfb + rnorm(n, sd = 1.0),
    -1.5 * tfb + rnorm(n, sd = 1.0),
    matrix(rnorm(n * 5), nrow = 5)
  )
  rownames(X) <- paste0("g", seq_len(12))
  truth <- matrix(0, 12, 12)
  edges <- rbind(c(1, 2), c(1, 3), c(1, 4), c(5, 6), c(5, 7))
  for (e in seq_len(nrow(edges))) {
    truth[edges[e, 1], edges[e, 2]] <- truth[edges[e, 2], edges[e, 1]] <- 1
  }
  list(X = X, truth = truth, modules = list(1:4, 5:7, 8:12))
}

# AUROC of a symmetric score matrix against a 0/1 adjacency, via the
# Mann-Whitney U statistic. Upper triangle only (undirected pairs).
pair_auroc <- function(scores, labels) {
  s <- scores[upper.tri(scores)]
  l <- labels[upper.tri(labels)]
  n1 <- sum(l == 1)
  n0 <- sum(l == 0)
  r <- rank(s)
  (sum(r[l == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

test_that("historical-parity CLR recovers the planted nonlinear network", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, threads = 1)$
    calibrate(method = "normal", combine = "euclidean")
  # Observed 1.000 across seeds; 0.95 leaves wide margin for platform noise.
  expect_gt(pair_auroc(fit$clr_scores, d$truth), 0.95)
})

test_that("new pipeline (FD bins + Stouffer) does not regress vs parity", {
  d <- synthetic_clr_data()
  old <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, threads = 1)$
    calibrate(method = "normal", combine = "euclidean")
  new <- ClrAttention$new(d$X)$
    estimate_mi(bins = "fd", threads = 1)$
    calibrate(method = "normal", combine = "stouffer")
  a_old <- pair_auroc(old$clr_scores, d$truth)
  a_new <- pair_auroc(new$clr_scores, d$truth)
  expect_gt(a_new, 0.95)
  expect_gte(a_new, a_old - 0.02)
})

test_that("permutation HC threshold keeps true edges and sparsifies", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, threads = 1)$
    calibrate(method = "normal", combine = "euclidean")$
    select_threshold(B = 20, method = "hc", threads = 1)
  tau <- fit$threshold
  expect_true(is.finite(tau) && tau > 0)
  expect_equal(fit$params$threshold$method, "hc")
  expect_equal(fit$params$threshold$B, 20)
  S <- fit$clr_scores
  ut <- upper.tri(S)
  kept <- S[ut] >= tau
  expect_lte(sum(kept), 10)                      # sparse: 4 of 66 observed
  expect_gte(sum(kept[d$truth[ut] == 1]), 4)     # 4 of 5 true edges observed
})

test_that("FDR threshold is a valid alternative selector", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, threads = 1)$
    calibrate(method = "normal", combine = "euclidean")$
    select_threshold(B = 20, method = "fdr", q = 0.05, threads = 1)
  tau <- fit$threshold
  expect_true(is.finite(tau) && tau > 0)
  S <- fit$clr_scores
  ut <- upper.tri(S)
  expect_gte(sum((S[ut] >= tau)[d$truth[ut] == 1]), 4)
})

test_that("build_operator() picks up the selected threshold automatically", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, threads = 1)$
    calibrate(method = "normal", combine = "euclidean")$
    select_threshold(B = 20, method = "hc", threads = 1)$
    build_operator(alpha = 0.5)
  op <- fit$operator
  S <- fit$clr_scores
  off <- row(op) != col(op)
  expect_equal(sum(op[off] > 0), sum(S[off] >= fit$threshold))
  rs <- rowSums(op)
  expect_true(all(abs(rs - 1) < 1e-12))
  # genes with no surviving edge carry exactly a self-loop
  iso <- rowSums(S >= fit$threshold & off) == 0
  expect_true(all(diag(op)[iso] == 1))
})

test_that("threshold needs select_threshold() first", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, threads = 1)$
    calibrate(method = "normal", combine = "euclidean")
  expect_error(fit$threshold, "select_threshold")
})

test_that("attention diffusion aggregates genes by planted module", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = "fd", threads = 1)$
    calibrate(method = "normal", combine = "stouffer")$
    build_operator(k = 4, alpha = 0.5)$
    diffuse(steps = 5)
  D <- as.matrix(dist(fit$embedding))  # gene-gene distances, genes are rows
  within <- mean(unlist(lapply(d$modules, function(m) {
    dm <- D[m, m]
    dm[upper.tri(dm)]
  })))
  block <- matrix(0, 12, 12)
  for (m in d$modules) block[m, m] <- 1
  up <- upper.tri(D)
  between <- mean(D[up][block[up] == 0])
  # Observed between/within ratio ~4-5x across seeds; 2x is a safe floor.
  expect_gt(between, 2 * within)
})
