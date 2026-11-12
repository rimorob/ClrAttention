# EXPLORATORY. Row-standardisation variants for the CLR step of CLR attention.
# Not preregistered, not for the paper or any paper table. See analysis/robust_z/robust_z.R.
#
# The package CLR is: z_ij = (x_ij - mean_i) / sd_i  (row i of the MI matrix, diagonal zeroed
# BEFORE the row mean and the sample sd, n - 1, are taken: src/core/clr_core.cpp), then
# negatives are truncated at 0, then the bilateral combination (z_ij + z_ji) / sqrt(2), zero
# diagonal. Every variant here changes ONLY the (centre, scale) of the row; truncation,
# symmetrisation, the softmax operator, the lazy walk and the symmetrised walk mass are the
# package's own code path (clr:::.softmax_rows, the .finish_operator row normalisation and
# self-loop, the lazy walk of analysis/baselines.R).
#
# Variant names (distinct from the project's existing 'rCLR'):
#   V0   mean / SD, the default CLR (package convention: the zeroed diagonal is part of the row)
#   V0x  mean / SD with the diagonal EXCLUDED from the row (control for the diagonal convention)
#   V1   median / (IQR / 1.349) of the raw MI row (the PI's proposal), diagonal excluded
#   V2   median / (IQR / 1.349) of MI^(1/3) (Wilson-Hilferty scale), one-way, diagonal excluded
#   V3   the existing 'rCLR' (eb_additive in analysis/synthetic_width/sim_width.R: two-way median
#        polish of MI^(1/3), per-gene MAD, bilateral (r_ij/s_i + r_ij/s_j)/sqrt(2), NOT truncated);
#        reference only
#   V4   median / (1.4826 MAD) of the raw MI row (diagnostic), diagonal excluded
# Zero-scale guard: a row whose scale is 0 (more than half of its values tied for IQR, MAD) gets
# z = 0 for the whole row, exactly as the package treats sd = 0 (src/core/clr_core.cpp); the
# number of such rows is returned per variant. No fallback to another scale (that would mix
# definitions inside one score matrix).

RZ_VARIANTS <- c("V0", "V1", "V2", "V3", "V4")

# G x G matrix -> (G - 1) x G matrix whose column i is row i of M without its diagonal element.
rz_offdiag <- function(M) {
  G <- nrow(M); Mt <- t(M)
  matrix(Mt[-((seq_len(G) - 1) * (G + 1) + 1)], G - 1L, G)
}

# column-wise quantiles (type 7, as stats::quantile) of a (G-1) x G matrix: rows = probs
rz_colq <- function(Y, probs) {
  if (requireNamespace("matrixStats", quietly = TRUE))
    t(matrixStats::colQuantiles(Y, probs = probs, type = 7, useNames = FALSE))
  else apply(Y, 2, stats::quantile, probs = probs, names = FALSE, type = 7)
}

# Row centre and scale of variant `v` on the MI matrix M (diagonal ignored/excluded as
# described above). Returns list(center, scale, transform) with transform in "raw" / "cube".
rz_center_scale <- function(M, v) {
  G <- nrow(M)
  if (v == "V0") {                     # package convention: diagonal zeroed and kept in the row
    M0 <- M; diag(M0) <- 0
    mu <- rowMeans(M0); sdv <- sqrt(rowSums((M0 - mu)^2) / (G - 1))
    return(list(center = mu, scale = sdv, transform = "raw"))
  }
  Y <- rz_offdiag(M)
  if (v == "V0x") {
    mu <- colMeans(Y); return(list(center = mu, scale = sqrt(colSums((Y - rep(mu, each = G - 1L))^2) / (G - 2)), transform = "raw"))
  }
  if (v == "V2") Y <- pmax(Y, 0)^(1 / 3)
  if (v %in% c("V1", "V2")) {
    q <- rz_colq(Y, c(0.25, 0.5, 0.75))
    return(list(center = q[2, ], scale = (q[3, ] - q[1, ]) / 1.349, transform = if (v == "V2") "cube" else "raw"))
  }
  if (v == "V4") {
    md <- rz_colq(Y, 0.5)[1, ]
    ad <- abs(Y - rep(md, each = G - 1L))
    return(list(center = md, scale = 1.4826 * rz_colq(ad, 0.5)[1, ], transform = "raw"))
  }
  stop("variant has no (centre, scale) form: ", v)
}

# existing rCLR (eb_additive, U only), copied from analysis/rclr_k562/rclr_k562.R rclr_of();
# returns the bilateral score matrix (untruncated, zero diagonal).
rz_rclr_existing <- function(M) {
  G <- nrow(M); Y <- M^(1 / 3); diag(Y) <- NA
  mu <- stats::median(Y, na.rm = TRUE); a <- rep(0, G)
  for (it in 1:8) { R <- Y - mu - outer(a, a, "+"); da <- apply(R, 1, stats::median, na.rm = TRUE)
    a <- a + da; mu <- mu + stats::median(a); a <- a - stats::median(a) }
  R <- Y - mu - outer(a, a, "+"); s <- apply(R, 1, stats::mad, na.rm = TRUE)
  U <- (sweep(R, 1, s, "/") + sweep(R, 2, s, "/")) / sqrt(2); diag(U) <- NA; U[is.na(U)] <- 0
  list(U = U, s = s)
}

# The row-standardised matrix Z (before truncation and symmetrisation; diagonal 0) of a variant,
# for V0, V0x, V1, V2, V4. Rows with scale <= 0 are zero (guard).
rz_zmatrix <- function(M, v) {
  cs <- rz_center_scale(M, v)
  Y <- if (cs$transform == "cube") pmax(M, 0)^(1 / 3) else M
  ok <- is.finite(cs$scale) & cs$scale > 0
  Z <- (Y - cs$center) / ifelse(ok, cs$scale, 1)
  Z[!ok, ] <- 0; diag(Z) <- 0
  attr(Z, "n_zero_scale") <- sum(!ok)
  Z
}

# Pairwise score matrix of a variant: truncate at 0, symmetrise (z_ij + z_ji)/sqrt(2), zero
# diagonal (V3: the existing rCLR bilateral matrix, as is). Returns the matrix with attribute
# n_zero_scale (rows hit by the zero-scale guard).
rz_scores <- function(M, v) {
  if (v == "V3") { r <- rz_rclr_existing(M); U <- r$U; attr(U, "n_zero_scale") <- sum(!(r$s > 0)); return(U) }
  Z <- rz_zmatrix(M, v); nz <- attr(Z, "n_zero_scale")
  Z[Z < 0] <- 0
  Z <- (Z + t(Z)) / sqrt(2); diag(Z) <- 0
  attr(Z, "n_zero_scale") <- nz
  Z
}

# Attention operator from a score matrix: the package's softmax (keff 10, cap 50), then the
# row normalisation and self-loop of ClrAttention$.finish_operator.
rz_operator <- function(C, keff = 10, cap = 50L) {
  A <- clr:::.softmax_rows(C, keff = keff, cap = cap)
  rs <- rowSums(A); nz <- rs > 0
  A[nz, ] <- A[nz, ] / rs[nz]
  if (any(!nz)) A[cbind(which(!nz), which(!nz))] <- 1
  A
}

rz_lazy_P <- function(A, alpha = 0.5)
  Matrix::Diagonal(nrow(A), 1 - alpha) + alpha * Matrix::Matrix(A, sparse = TRUE)

# Symmetrised walk mass S_t = (P^t + P^t')/2 (diagonal kept, as analysis/baselines.R; scoring
# zeroes it) for every t in `ts`. Returns a list named by t.
rz_walk <- function(P, ts) {
  out <- list(); Pt <- as.matrix(P)
  if (1L %in% ts) out[["1"]] <- (Pt + t(Pt)) / 2
  for (t in 2:max(ts)) { Pt <- as.matrix(Pt %*% P)
    if (t %in% ts) out[[as.character(t)]] <- (Pt + t(Pt)) / 2 }
  out
}

# stationary distribution of the lazy walk with teleport 1e-3 (analysis/baselines.R hub diagnostic)
rz_stationary <- function(P, eps = 1e-3) {
  G <- nrow(P); pi <- rep(1 / G, G)
  for (it in 1:20000) { q <- (1 - eps) * as.numeric(pi %*% P) + eps / G; if (sum(abs(q - pi)) < 1e-13) break; pi <- q }
  q
}

# within-stratum percentile in k strata of H over the upper triangle (analysis/baselines.R)
rz_hub_stratify <- function(S, H, k = 10L) {
  S[!is.finite(S)] <- 0; ut <- which(upper.tri(S)); h <- H[ut]
  g <- cut(h, unique(stats::quantile(h, seq(0, 1, length.out = k + 1L))), include.lowest = TRUE)
  r <- stats::ave(S[ut], g, FUN = function(x) rank(x) / length(x))
  R <- matrix(0, nrow(S), ncol(S)); R[ut] <- r; R + t(R)
}
