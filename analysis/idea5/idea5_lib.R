# EXPLORATORY (docs/plans/explore_idea5_plan.md). Shared helpers for the learned-attention test on M3D.
# Not preregistered, not for the paper. Everything after the learned logits is the package / robust_z code path.

# ---- data (as analysis/robust_z/robust_z.R, M3D branch) ---------------------------------------------
idea5_load_m3d <- function(quick = 0L, quick_seed = 1L) {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
  exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description); exper[!grepl("experimenter:", desc$description)] <- NA
  lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])[colnames(X)]
  lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab)))
  storage.mode(X) <- "double"
  if (quick > 0L) { set.seed(quick_seed); X <- X[sort(sample(nrow(X), quick)), ] }   # random genes, no truth used
  list(X = X, lab = unname(lab))
}

# balanced lab split of analysis/lab_splits.R (shuffle labs, give each to the half with fewer experiments)
idea5_lab_split <- function(lab) {
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab])
}

# CLR score matrix exactly as the package pipeline (analysis/lab_splits.R build_fit)
idea5_fit <- function(Xm, threads = 8L) {
  rng <- apply(Xm, 1, function(x) diff(range(x)))
  if (any(rng <= 0)) { set.seed(1L); Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9) }
  f <- clr::ClrAttention$new(Xm)$estimate_mi(bins = "hg", transform = "none", threads = threads)
  f$calibrate(method = "normal", combine = "stouffer")
  f
}

# ---- softmax temperatures: a port of the bisection in clr:::.softmax_rows, returning T_i --------------
idea5_row_temps <- function(S, keff = 10, cap = 50L) {
  G <- nrow(S); cap <- as.integer(cap); Tm <- rep(NA_real_, G)
  for (i in seq_len(G)) {
    row <- S[i, ]; row[i] <- -Inf
    top <- order(row, decreasing = TRUE)[seq_len(min(cap, G - 1L))]
    top <- top[row[top] > 0]
    if (!length(top) || length(top) <= keff) next
    s <- row[top]
    neff <- function(lt) { w <- exp((s - max(s)) / exp(lt)); w <- w / sum(w); exp(-sum(w[w > 0] * log(w[w > 0]))) }
    lo <- log(1e-4); hi <- log(1e4)
    for (it in 1:60) { mid <- (lo + hi) / 2; if (neff(mid) < keff) lo <- mid else hi <- mid }
    Tm[i] <- exp((lo + hi) / 2)
  }
  Tm
}
# weights implied by the temperatures (for the equivalence check against .softmax_rows)
idea5_softmax_from_temps <- function(S, Tm, keff = 10, cap = 50L) {
  G <- nrow(S); A <- matrix(0, G, G)
  for (i in seq_len(G)) {
    row <- S[i, ]; row[i] <- -Inf
    top <- order(row, decreasing = TRUE)[seq_len(min(cap, G - 1L))]; top <- top[row[top] > 0]
    if (!length(top)) next
    s <- row[top]
    if (length(top) <= keff) { A[i, top] <- 1 / length(top); next }
    w <- exp((s - max(s)) / Tm[i]); A[i, top] <- w / sum(w)
  }
  A
}

# ---- binary matrices shared with python ---------------------------------------------------------------
idea5_write_f64 <- function(M, path) { con <- file(path, "wb"); writeBin(as.vector(M), con, size = 8L); close(con) }
idea5_read_mat <- function(path, G, size = 4L) {     # python writes D.T (C order) = column-major D
  con <- file(path, "rb"); on.exit(close(con)); x <- readBin(con, "numeric", n = G * G, size = size)
  if (length(x) != G * G) stop("short read: ", path); matrix(x, G, G)
}
idea5_read_vec <- function(path, n) { con <- file(path, "rb"); on.exit(close(con)); readBin(con, "numeric", n = n, size = 8L) }

# ---- the learned-logit readout (plan section 2) -------------------------------------------------------
# L_ij = a_i Z_ij / tau0 + D_ij (Z NULL for scratch models); C = pmax(tau0 (L + L')/2, 0), zero diagonal.
idea5_score_matrix <- function(D, Z = NULL, a = NULL, tau0) {
  L <- D; diag(L) <- 0
  if (!is.null(Z)) L <- L + (a * Z) / tau0              # a recycles down columns: element [i, j] times a[i]
  C <- tau0 * (L + t(L)) / 2; rm(L)
  C[C < 0] <- 0; diag(C) <- 0
  C
}
