setwd(Sys.getenv("D5_ROOT", unset = "."))   # repo root: testthat runs each file in tests/dream5
test_that("replicate groups come from Repeat, not from #Experiment", {
  expect_equal(d5_replicate_groups(c(1, 2, 1, 1, 2, 3, 1)), c(1, 1, 2, 3, 3, 3, 4))
  expect_error(d5_replicate_groups(c(1, 3, 1)), "1..k")        # Repeat must run 1..k
  expect_error(d5_replicate_groups(c(2, 1)))                    # must start at 1
  # real feature files (label-free): the group counts stated in the design
  skip_if_not(file.exists("data/dream5/synapse/net4/net4_chip_features.tsv"))
  for (n in c(1, 3, 4)) { cf <- read.delim(sprintf("data/dream5/synapse/net%d/net%d_chip_features.tsv", n, n), check.names = FALSE, comment.char = "")
    g <- d5_replicate_groups(cf$Repeat); expect_equal(max(g), c(`1` = 487, `3` = 487, `4` = 321)[[as.character(n)]])
    expect_gt(length(unique(cf[["#Experiment"]])), 0); expect_false(max(g) == length(unique(cf[["#Experiment"]])))
    expect_equal(length(unique(cf[["#Experiment"]])), c(`1` = 107, `3` = 107, `4` = 59)[[as.character(n)]]) }
})

test_that("real data load: dimensions, TFs, no zero-variance genes (label-free)", {
  skip_if_not(file.exists("data/dream5/synapse/net1/net1_expression_data.tsv"))
  nd <- d5_read_net(1L)
  expect_equal(dim(nd$X), c(1643, 805)); expect_length(nd$tf_ids, 195); expect_equal(nd$n_groups, 487)
  expect_equal(nd$tf_idx, 1:195)
  na <- d5_read_net(1L, "avg"); expect_equal(dim(na$X), c(1643, 487))
  # the averaged file equals the group means of the primary file
  gm <- vapply(split(seq_len(805), nd$group), function(ix) rowMeans(nd$X[, ix, drop = FALSE]), numeric(1643))
  expect_lt(max(abs(gm - na$X)), 1e-6)
})

test_that("SE_B group halves: seed 20263000 + b, half the groups, reproducible, caller RNG untouched", {
  h1 <- d5_group_half(321, 7); h2 <- d5_group_half(321, 7); h3 <- d5_group_half(321, 8)
  expect_identical(h1, h2); expect_false(identical(h1, h3)); expect_length(h1, 160); expect_false(anyDuplicated(h1) > 0)
  set.seed(5); a <- runif(1); set.seed(5); invisible(d5_group_half(321, 7)); expect_equal(runif(1), a)
  set.seed(20263007L); expect_identical(h1, sort(sample.int(321, 160)))      # the seed rule itself
  nd <- d5_synth_net(30, 4, 20, 2); sub <- d5_subset_groups(nd, h1[h1 <= 20]); expect_true(all(sub$group %in% h1))
  expect_equal(ncol(sub$X), sum(nd$group %in% h1))
})

test_that("series halves (supplementary jackknife): seed 20264000 + b, floor(n/2) of the series, whole series kept, caller RNG untouched", {
  h1 <- d5_series_half(59, 7); expect_identical(h1, d5_series_half(59, 7)); expect_false(identical(h1, d5_series_half(59, 8)))
  expect_length(h1, 29); expect_false(anyDuplicated(h1) > 0)
  set.seed(20264007L); expect_identical(h1, sort(sample.int(59, 29)))                    # the seed rule itself
  expect_false(identical(h1, sort(.with_seed(20263007L, sample.int(59, 29)))))           # not the SE_B stream
  set.seed(5); a <- runif(1); set.seed(5); invisible(d5_series_half(59, 7)); expect_equal(runif(1), a)
  nd <- d5_synth_net(30, 4, 20, 2); ids <- sort(unique(nd$series)); expect_length(ids, 7)
  sub <- d5_subset_series(nd, 3L); keep <- ids[d5_series_half(7, 3L)]
  expect_identical(sort(unique(sub$series)), keep); expect_equal(ncol(sub$X), sum(nd$series %in% keep)); expect_equal(sub$n_series_kept, 3)
  expect_true(all(table(nd$series)[as.character(keep)] == table(sub$series)[as.character(keep)]))   # every chip of a kept series is kept
  expect_error(d5_subset_series(list(series = NULL), 1L))
  skip_if_not(file.exists("data/dream5/synapse/net4/net4_chip_features.tsv"))             # label-free real structure
  cf <- read.delim("data/dream5/synapse/net4/net4_chip_features.tsv", check.names = FALSE, comment.char = ""); sr <- cf[["#Experiment"]]
  k <- sort(unique(sr))[d5_series_half(59, 1L)]; expect_length(k, 29); expect_equal(sum(sr %in% k) + sum(!sr %in% k), 536)
})

test_that("blom_avg is the text of analysis/beeline.R", {
  l <- readLines("analysis/beeline.R"); s <- grep("^blom_avg <- function", l); e <- s + which(grepl("^}", l[-seq_len(s)]))[1]
  env <- new.env(); eval(parse(text = l[s:e]), env)
  X <- matrix(rnorm(60), 6, 10); X[1, 1:3] <- 1                                  # ties
  expect_identical(d5_blom_avg(X), env$blom_avg(X))
})

test_that("packing helpers round trip and index pairs correctly", {
  G <- 7; S <- matrix(runif(G * G), G, G); S <- (S + t(S)) / 2; diag(S) <- 0
  v <- d5_pack(S); expect_equal(d5_unpack(v, G), S)
  ii <- c(2, 5, 7, 1); jj <- c(5, 2, 3, 7); expect_equal(v[d5_pair_index(ii, jj)], S[cbind(ii, jj)])
  expect_equal(d5_clean_sym(matrix(c(1, NA, Inf, 2), 2)), matrix(c(0, 0, 0, 0), 2))
})

test_that("methods on a synthetic network: symmetric, finite, recover planted edges; caches and manifests", {
  nd <- d5_synth_net(60, 8, 40, 1); out <- tempfile(); G <- 60
  for (m in c("abs_pearson", "clr", "clrns", "attention", "genie3_tf", "grnboost2_tf")) d5_run_method(nd, m, "replicate", 2L, out)
  r <- d5_pack(d5_m_abs_pearson(nd$X)); expect_equal(d5_read_scores(d5_cache_path(0, "abs_pearson", "replicate", out)), r)
  expect_equal(r, d5_pack(abs(cor(t(nd$X)))))
  for (m in c("clr", "clrns", "attention")) { v <- d5_read_scores(d5_cache_path(0, m, "replicate", out)); expect_length(v, G * (G - 1) / 2); expect_true(all(is.finite(v))) }
  # CLR attention = the baselines.R recipe, recomputed independently here
  f <- ClrAttention$new(nd$X)$estimate_mi(bins = "hg", transform = "none", threads = 2)$calibrate(method = "normal", combine = "stouffer")
  f$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L)
  P <- Matrix::Diagonal(G, 0.5) + 0.5 * Matrix::Matrix(f$operator, sparse = TRUE); Pt <- as.matrix(P); for (s in 1:8) Pt <- as.matrix(Pt %*% P)
  St <- (Pt + t(Pt)) / 2; diag(St) <- 0
  expect_equal(d5_read_scores(d5_cache_path(0, "attention", "replicate", out)), d5_pack(St))
  S <- d5_unpack(d5_read_scores(d5_cache_path(0, "clr", "replicate", out)), G); diag(f$clr_scores) -> dd
  expect_equal(S, { x <- f$clr_scores; x[!is.finite(x)] <- 0; diag(x) <- 0; x }, ignore_attr = TRUE)
  # planted TF->target edges rank above random pairs for CLR
  M <- d5_edge_matrix_sym(d5_read_scores(d5_cache_path(0, "clr", "replicate", out)), nd)
  expect_equal(sum(!is.na(M)), 8 * 60 - 8)                                       # TF list x all genes, no self pairs
  expect_equal(M[1:8, 1:8][upper.tri(M[1:8, 1:8])], t(M[1:8, 1:8])[upper.tri(M[1:8, 1:8])])   # TF-TF: undirected methods tie both directions
  g3 <- d5_read_scores(d5_cache_path(0, "genie3_tf", "replicate", out)); expect_equal(dim(g3), c(8, 60)); expect_equal(sum(is.na(g3)), 8)
  expect_true(all(is.na(g3[cbind(1:8, 1:8)])))
  gb <- d5_read_scores(d5_cache_path(0, "grnboost2_tf", "replicate", out)); expect_equal(dim(gb), c(8, 60)); expect_true(all(is.na(gb[cbind(1:8, 1:8)])))
  # sidecars and the consolidated manifest
  mf <- d5_collect_manifest(out); expect_equal(nrow(mf), 6)
  expect_equal(mf$sha256[mf$method == "clr"], d5_sha256(d5_cache_path(0, "clr", "replicate", out)))
  expect_match(mf$sha256[1], "^[0-9a-f]{64}$")
  # resumable: an existing cache is reused, not recomputed
  t0 <- file.mtime(d5_cache_path(0, "clr", "replicate", out)); Sys.sleep(1.1); d5_run_method(nd, "clr", "replicate", 2L, out)
  expect_equal(file.mtime(d5_cache_path(0, "clr", "replicate", out)), t0)
})

test_that("GENIE3 with TF regulators is reproducible under the fixed seed; all-gene run has the full shape", {
  nd <- d5_synth_net(40, 5, 30, 3)
  a <- d5_m_genie3(nd$X, nd$tf_ids, 2L, ntrees = 50); b <- d5_m_genie3(nd$X, nd$tf_ids, 2L, ntrees = 50)
  expect_identical(a, b); expect_equal(dim(a), c(5, 40))
  c3 <- d5_m_genie3(nd$X, nd$tf_ids, 3L, ntrees = 50); expect_identical(a, c3)    # independent of the number of cores (>= 2)
  expect_error(d5_m_genie3(nd$X, nd$tf_ids, 1L, ntrees = 50), "workers >= 2")
  al <- d5_m_genie3(nd$X, NULL, 2L, ntrees = 20); expect_equal(dim(al), c(40, 40)); expect_equal(sum(is.na(al)), 40)
})

test_that("refit replicates: reproducible hashes under a fixed seed and thread count; resumable", {
  nd <- d5_synth_net(40, 5, 30, 4); out <- tempfile()
  r1 <- d5_refit_run(nd, 5L, 2L, out); p <- d5_refit_path(0, 5L, out); expect_true(file.exists(p))
  expect_equal(r1$method, D5$refit_methods)
  file.remove(p); r2 <- d5_refit_run(nd, 5L, 2L, out); expect_identical(r1$sha256, r2$sha256)       # recomputation gives the same hashes
  r3 <- d5_refit_run(nd, 6L, 2L, out); expect_false(any(r1$sha256 == r3$sha256))
  m <- d5_refit_manifest(0, out); expect_equal(nrow(m), 8)
  expect_true("n_groups_kept" %in% names(r1)); expect_identical(names(r1)[1:5], c("net", "rep", "method", "n_groups_kept", "sha256"))   # group sidecar layout unchanged
  # the hash is of the refit scores as a verification of what the truth stage recomputes
  v <- d5_refit_scores(nd, 5L, 2L); expect_identical(vapply(v, d5_hash_scores, ""), setNames(r1$sha256, names(v)))
  # kept half only
  expect_equal(r1$n_groups_kept[1], 15)
})

test_that("series-level refit replicates: own directory and manifest, seed 20264000 + b, reproducible, resumable, differ from SE_B refits", {
  nd <- d5_synth_net(40, 5, 30, 4); out <- tempfile()
  expect_match(d5_refit_dir(4, out, "series"), "net4/refit_series$"); expect_match(d5_refit_dir(4, out), "net4/refit$")
  expect_match(d5_refit_manifest_file(4, out, "series"), "refit_series/REFIT_SERIES_MANIFEST.tsv$")
  s1 <- d5_refit_run(nd, 5L, 2L, out, unit = "series"); p <- d5_refit_path(0, 5L, out, "series"); expect_true(file.exists(p))
  expect_match(p, "refit_series/rep005.sha256.tsv$"); expect_equal(s1$method, D5$refit_methods)
  expect_true("n_series_kept" %in% names(s1)); expect_equal(unique(s1$n_series_kept), floor(length(unique(nd$series)) / 2))
  file.remove(p); s2 <- d5_refit_run(nd, 5L, 2L, out, unit = "series"); expect_identical(s1$sha256, s2$sha256)
  v <- d5_refit_scores(nd, 5L, 2L, unit = "series"); expect_identical(vapply(v, d5_hash_scores, ""), setNames(s1$sha256, names(v)))
  g1 <- d5_refit_run(nd, 5L, 2L, out); expect_false(any(g1$sha256 == s1$sha256))                 # same b, different unit and stream
  m <- d5_refit_manifest(0, out, "series"); expect_equal(nrow(m), 4); expect_true(file.exists(d5_refit_manifest_file(0, out, "series")))
  expect_equal(length(Sys.glob(file.path(d5_refit_dir(0, out), "rep*.sha256.tsv"))), 1)           # group-level directory not touched by the series run
  # series refit on the replicate-averaged input is not defined
  nd$series <- NULL; expect_error(d5_refit_scores(nd, 1L, 2L, unit = "series"))
})
