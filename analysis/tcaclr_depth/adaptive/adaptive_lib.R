# Shared helpers for the adaptive (label-free) setting studies of docs/plans/tcaclr_depth_grid_plan.md, Revision 2 items 3, 8, 9, 11, 12.
# EXPLORATORY. Source from the repo root:  source("analysis/tcaclr_depth/adaptive/adaptive_lib.R")
# The operator is built by the grid library (gl_make_P: clr:::.softmax_rows + the registered finish step + lazy weight 0.5); nothing here builds a second one.
suppressPackageStartupMessages({ source("analysis/tcaclr_depth/grid_lib.R"); library(RSpectra) })

AD <- list(seed = 20261121L, alpha = GL$alpha, eps_edge = c(0.05, 0.01), unit_tol = 1e-6)

# NOTE: forking (mclapply) after the calling process has run clr's estimate_mi (OpenMP threads) deadlocks the children (observed). Keep every clr call inside the
# tasks passed here when workers > 1; the calling process may do BLAS work (gl_neff, svds) but never fits.
ad_lapply <- function(x, f, workers = 1L)           # errors are returned as strings so that one failure does not lose the rest
  if (workers > 1L) parallel::mclapply(x, function(i) tryCatch(f(i), error = function(e) paste("ERROR", conditionMessage(e))), mc.cores = workers, mc.preschedule = FALSE) else lapply(x, f)
ad_mem <- function() gl_vmhwm_gb()

## ---- data -----------------------------------------------------------------------------------------------------------------------
ad_load <- function(path, genes = 0L, seed = AD$seed) {       # as stage1_labelfree.R: constant genes dropped; --genes = a fixed random slice
  X <- gl_read_matrix(path); X <- X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]
  if (genes > 0L && genes < nrow(X)) { set.seed(seed); X <- X[sort(sample(nrow(X), genes)), , drop = FALSE] }
  X
}
ad_permute_genes <- function(X, seed) {                       # each gene's values permuted across experiments independently
  set.seed(seed); P <- X; for (i in seq_len(nrow(X))) P[i, ] <- X[i, sample.int(ncol(X))]; P
}
ad_rank_normal <- function(X) {                               # per gene across experiments, qnorm((rank - 0.5)/N), re-standardised (ties)
  Z <- t(apply(X, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 0.5) / length(x))))
  Z <- Z - rowMeans(Z); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1)); s[s == 0] <- 1; Z / s
}
ad_fit <- function(X, threads = 8L) {                          # as lab_splits.R build_fit: constant rows get 1e-9 noise (can occur in a half)
  rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0)) X[rng <= 0, ] <- X[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9)
  gl_fit(X, threads)
}

## ---- blocks and halves (as stage1_labelfree.R) -----------------------------------------------------------------------------------
ad_bal_split <- function(unit, seed) {                         # shuffle units, give each to the half with fewer columns so far
  set.seed(seed); sz <- table(unit); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (u in names(sz)) { h <- which.min(tot); half[u] <- h; tot[h] <- tot[h] + sz[[u]] }
  unname(half[as.character(unit)])
}
# kind "lab": the stored laboratory splits (results/lab_splits/split_assignments.csv; unit = laboratory, halves = the stored split columns)
# kind "kmeans": resampling units = k-means clusters of the knockdown profiles (columns): top 30 principal components of the gene-centred,
#   experiment-standardised matrix, k = round(N_eff) clusters (--clusters overrides), seed fixed; each split shuffles clusters (seed + split).
# kind "random": experiments are the units in 8 equal groups (no block structure available; used by the planted simulation)
ad_units <- function(X, kind, splits, seed = AD$seed, lab_splits = "results/lab_splits/split_assignments.csv", clusters = 0L, nblocks = 8L) {
  if (kind == "lab") {
    sa <- utils::read.csv(lab_splits, check.names = FALSE, stringsAsFactors = FALSE); stopifnot(all(colnames(X) %in% sa$experiment))
    sa <- sa[match(colnames(X), sa$experiment), ]; sc <- grep("^split_", names(sa), value = TRUE)[seq_len(splits)]
    return(list(unit = sa$lab, splits = lapply(sc, function(s) sa[[s]]), k = length(unique(sa$lab))))
  }
  if (kind == "kmeans") {
    k <- if (clusters > 0L) clusters else max(2L, as.integer(round(gl_neff(X)[["N_eff"]])))
    Zc <- scale(X - rowMeans(X)); Zc[!is.finite(Zc)] <- 0
    pcs <- RSpectra::svds(Zc, k = 30)$v                                   # experiments x 30 loadings (right singular vectors)
    set.seed(seed); km <- stats::kmeans(pcs, centers = k, nstart = 3L, iter.max = 100L)
    unit <- paste0("km", km$cluster)
  } else { set.seed(seed); unit <- paste0("blk", sample(rep_len(seq_len(nblocks), ncol(X)))); k <- nblocks }
  list(unit = unit, splits = lapply(seq_len(splits), function(s) ad_bal_split(unit, seed + s)), k = k)
}

## ---- the operator in symmetrised-by-degree form ------------------------------------------------------------------------------------
# P = 0.5 I + 0.5 A (gl_make_P). A is a row-stochastic attention matrix and is NOT reversible, so D^{1/2} P D^{-1/2} is not symmetric with
# D = stationary mass. Symmetrised form used here: W = (A + A') / 2 (the undirected attention graph), d = rowSums(W),
#   P_sym = 0.5 I + 0.5 D^{-1/2} W D^{-1/2}  (symmetric, real spectrum in [0, 1], trivial eigenvalue 1 with eigenvector ~ sqrt(d)).
# The eigenvalues of the non-symmetric P itself are returned as a cross-check (ad_edge_nonsym).
ad_P_from_fit <- function(f, keff, cap = gl_cap(keff), clipS = TRUE) gl_make_P(gl_pair_score(f, clipS), keff, cap)
ad_Psym <- function(P, alpha = AD$alpha) {   # returns a symmetric dgCMatrix
  G <- nrow(P); A <- (P - Matrix::Diagonal(G, 1 - alpha)) / alpha
  A <- methods::as(A, "generalMatrix"); W <- (A + Matrix::t(A)) / 2; d <- Matrix::rowSums(W)
  Dm <- Matrix::Diagonal(G, 1 / sqrt(d))
  Ps <- Matrix::Diagonal(G, 1 - alpha) + alpha * (Dm %*% W %*% Dm)
  methods::as(Matrix::drop0(Matrix::forceSymmetric(Ps)), "generalMatrix")   # dgCMatrix, exactly symmetric (eigs_sym reads the lower triangle)
}
# top-k eigenpairs (largest first); the trivial 1 is value 1. vectors = FALSE returns values only.
ad_eigs <- function(Ps, k, vectors = FALSE) {
  k <- min(k, nrow(Ps) - 2L)
  r <- RSpectra::eigs_sym(Ps, k = k, which = "LA", opts = list(ncv = min(nrow(Ps) - 1L, max(2L * k + 1L, 40L)), maxitr = 5000L, tol = 1e-9, retvec = vectors))
  o <- order(r$values, decreasing = TRUE); r$values <- r$values[o]; if (vectors) r$vectors <- r$vectors[, o, drop = FALSE]; r
}
# second-largest modulus of the non-symmetric P (diagnostic): eigenvalues of P and of D^1/2 P D^-1/2 coincide for reversible chains only
ad_edge_nonsym <- function(P, k = 6L) {
  v <- tryCatch(RSpectra::eigs(methods::as(P, "generalMatrix"), k = k, which = "LM", opts = list(ncv = max(2L * k + 1L, 30L), maxitr = 5000L, tol = 1e-9))$values, error = function(e) NA_real_)
  m <- sort(Mod(v), decreasing = TRUE); if (length(m) >= 2L) m[2L] else NA_real_
}
ad_tstar <- function(edge, eps) ifelse(edge >= 1 | edge <= 0, NA_real_, log(eps) / log(edge))
