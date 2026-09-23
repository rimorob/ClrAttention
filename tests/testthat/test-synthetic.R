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

# 1 for pairs inside a planted regulatory module (TF + its targets).
same_module <- function(d) {
  lab <- integer(12)
  lab[1:4] <- 1L; lab[5:7] <- 2L; lab[8:12] <- 3:7
  outer(lab, lab, "==") * 1
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
    estimate_mi(bins = 10, transform = "none", threads = 1)$
    calibrate(method = "normal", combine = "euclidean")
  # Observed 1.000 across seeds; 0.95 leaves wide margin for platform noise.
  expect_gt(pair_auroc(fit$clr_scores, d$truth), 0.95)
})

test_that("new pipeline (FD bins + Stouffer) does not regress vs parity", {
  d <- synthetic_clr_data()
  old <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, transform = "none", threads = 1)$
    calibrate(method = "normal", combine = "euclidean")
  new <- ClrAttention$new(d$X)$
    estimate_mi(bins = "fd", transform = "none", threads = 1)$
    calibrate(method = "normal", combine = "stouffer")
  a_old <- pair_auroc(old$clr_scores, d$truth)
  a_new <- pair_auroc(new$clr_scores, d$truth)
  expect_gt(a_new, 0.95)
  expect_gte(a_new, a_old - 0.02)
})

test_that("permutation HC threshold keeps true edges and sparsifies", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(threads = 1)$calibrate()$
    select_threshold(B = 20, method = "hc", threads = 1)
  expect_true(is.finite(fit$threshold) && fit$threshold > 0)
  expect_equal(fit$params$threshold$method, "hc")
  expect_equal(fit$params$threshold$statistic, "clr")
  expect_equal(fit$params$threshold$B, 20)
  ut <- upper.tri(d$truth)
  kept <- fit$edges[ut]
  # MI selects dependence, so co-regulated sibling pairs (3-4, 6-7) are
  # legitimate selections; false positives are pairs across modules.
  expect_lte(sum(kept[same_module(d)[ut] == 0]), 1)
  expect_gte(sum(kept[d$truth[ut] == 1]), 4)     # >= 4 of 5 true edges
})

test_that("FDR threshold on the MI null selects the planted edges (small G)", {
  # 12 genes is far below the regime CLR's context null is built for, so
  # this small-G unit test uses the MI null; the CLR-null default is tested
  # at G = 200 below.
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(threads = 1)$calibrate()$
    select_threshold(B = 20, statistic = "mi", threads = 1)
  expect_equal(fit$params$threshold$method, "fdr")
  ut <- upper.tri(d$truth)
  expect_gte(sum(fit$edges[ut][d$truth[ut] == 1]), 4)
  expect_lte(sum(fit$edges[ut][same_module(d)[ut] == 0]), 2)
})

test_that("statistic = 'clr' still works on the parity pipeline", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, transform = "none", threads = 1)$
    calibrate(method = "normal", combine = "euclidean")$
    select_threshold(B = 20, statistic = "clr", threads = 1)
  S <- fit$clr_scores
  expect_identical(fit$edges, S >= fit$threshold & row(S) != col(S))
  ut <- upper.tri(S)
  expect_gte(sum(fit$edges[ut][d$truth[ut] == 1]), 4)
})

test_that("MI-statistic null recovers large modules the CLR null misses", {
  # 24 genes, three planted modules of 8/7/6: modules are large relative to
  # G, so CLR row backgrounds are dominated by the module itself and CLR
  # scores are compressed below the (pure-noise) permuted CLR null.
  set.seed(7)
  n <- 100
  X <- matrix(rnorm(24 * n), 24)
  mods <- list(1:8, 9:15, 16:21)
  for (m in mods) {
    f <- rnorm(n)
    for (g in m) X[g, ] <- 0.85 * f + sqrt(1 - 0.85^2) * X[g, ]
  }
  lab <- c(rep(1, 8), rep(2, 7), rep(3, 6), 4:6)
  truth <- outer(lab, lab, "==") & !diag(24)
  a <- ClrAttention$new(X)$estimate_mi(threads = 1)$calibrate()
  a$select_threshold(B = 20, statistic = "clr", threads = 1)
  n_clr <- sum(a$edges)
  a$select_threshold(B = 20, statistic = "mi", threads = 1)
  E <- a$edges
  expect_gt(sum(E & truth), 0.8 * sum(truth))
  expect_lte(sum(E & !truth) / sum(E), 0.10)  # BH at q = 0.05: FDP ~ q
  expect_lt(n_clr, sum(E))
})

test_that("build_operator() picks up the selected edge set automatically", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(threads = 1)$calibrate()$
    select_threshold(B = 20, statistic = "mi", threads = 1)$
    build_operator(alpha = 0.5)
  op <- fit$operator
  S <- fit$clr_scores
  off <- row(op) != col(op)
  expect_identical(op[off] > 0, (fit$edges & S > 0)[off])
  expect_true(all(abs(rowSums(op) - 1) < 1e-12))
  iso <- rowSums(fit$edges & S > 0) == 0   # no weighted edge -> self-loop
  expect_true(all(diag(op)[iso] == 1))
  expect_equal(fit$params$operator$selection, "mi")
})

test_that("threshold needs select_threshold() first", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = 10, transform = "none", threads = 1)$
    calibrate(method = "normal", combine = "euclidean")
  expect_error(fit$threshold, "select_threshold")
})

test_that("attention diffusion aggregates genes by planted module", {
  d <- synthetic_clr_data()
  fit <- ClrAttention$new(d$X)$
    estimate_mi(bins = "fd", transform = "none", threads = 1)$
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

test_that("calibrated HC selects nothing on pure noise (seeded)", {
  set.seed(101)
  X <- rbind(matrix(rnorm(20 * 300), 20), matrix(rt(20 * 300, 3), 20))
  fit <- ClrAttention$new(X)$
    estimate_mi(bins = "fd", transform = "rank", threads = 1)$
    calibrate()
  expect_message(fit$select_threshold(B = 20, method = "hc", threads = 1),
                 "no signal")
  expect_identical(fit$threshold, Inf)
  hc <- fit$params$threshold$hc
  expect_length(hc$hc_null, 20)
  expect_lte(hc$hc_star, hc$hc_crit)
})

test_that("rank transform makes MI invariant to monotone distortions", {
  d <- synthetic_clr_data()
  X2 <- d$X
  X2[2, ] <- exp(X2[2, ])            # strictly monotone per-gene maps
  X2[6, ] <- X2[6, ]^3
  a <- bspline_mi(d$X, bins = "fd", transform = "rank", threads = 1)
  b <- bspline_mi(X2, bins = "fd", transform = "rank", threads = 1)
  expect_equal(a, b)
})

test_that("default CLR null stays sparse under a global confounder; MI null does not", {
  # 200 genes x 907 samples, ten planted 8-gene modules, plus a global factor
  # loading every gene (growth-rate-like program). Selection must answer
  # "exceptional relative to each gene's background", not "any dependence".
  set.seed(3)
  n <- 907; G <- 200
  X <- matrix(rnorm(G * n), G)
  lab <- c(rep(1:10, each = 8), 11:(10 + G - 80))
  for (m in 1:10) {
    f <- rnorm(n)
    for (i in which(lab == m)) X[i, ] <- 0.5 * f + sqrt(0.75) * X[i, ]
  }
  X <- X + 0.35 * matrix(rnorm(n), G, n, byrow = TRUE)
  truth <- outer(lab, lab, "==") & !diag(G)
  a <- ClrAttention$new(X)$estimate_mi(threads = 2)$calibrate()
  set.seed(1)
  a$select_threshold(B = 10, threads = 2)            # default: FDR, CLR null
  E <- a$edges
  frac_clr <- sum(E) / (G * (G - 1))
  expect_equal(a$params$threshold$statistic, "clr")
  expect_gt(sum(E & truth) / sum(E), 0.9)            # precise
  expect_gt(sum(E & truth) / sum(truth), 0.6)        # and substantive recall
  set.seed(1)
  a$select_threshold(B = 10, statistic = "mi", threads = 2)
  frac_mi <- sum(a$edges) / (G * (G - 1))
  expect_gt(frac_mi, 0.1)                            # MI null: "everything"
  expect_gt(frac_mi, 5 * frac_clr)
})

test_that("BH cutoff keeps every member of a tied rejected block", {
  # six scores share the smallest (histogram) p-value; order() may put any
  # of them last in p-order, but tau must keep all six
  p <- c(rep(1e-4, 6), 0.3, 0.6, 0.9)
  s <- c(6.9, 6.7, 7.4, 6.8, 7.1, 7.0, 1, 0.5, 0.1)
  o <- order(p)
  tau <- clr:::.fdr_cutoff(p[o], s[o], M = 9, q = 0.05)
  expect_equal(tau, 6.7)
  expect_equal(sum(s >= tau), 6)
})
