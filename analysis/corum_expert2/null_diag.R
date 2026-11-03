# Helpers for docs/plans/corum_expert2.md, Task B (EXPLORATORY, label-free diagnostics).
# Row-null diagnostics of a dependence matrix D (G x G, symmetric, zero diagonal) against
# a within-gene permutation null Dp (same construction on data whose cells/profiles were
# permuted independently for each gene), plus label-free alternative row calibrations.
# No truth is read here.

# CLR's own row statistics: mean and sample sd over all G entries of the row, the zero
# diagonal included (src/core/clr_packed.cpp; clr.m). Returns the unclipped directed z
# matrix Z[i, j] = (D[i, j] - mu_i) / sd_i (row i judges pair (i, j)).
clr_rowz <- function(D) {
  G <- nrow(D); mu <- rowSums(D) / G
  sdv <- sqrt(rowSums((D - mu)^2) / (G - 1)); sdv[sdv <= 0] <- 1
  (D - mu) / sdv
}
stouffer <- function(Z) { Z[!is.finite(Z) | Z < 0] <- 0; S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S }
offd <- function(M) { diag(M) <- NA; M }
row_moments <- function(D) {         # off-diagonal skewness and excess kurtosis per row
  X <- offd(D); mu <- rowMeans(X, na.rm = TRUE); C <- X - mu
  m2 <- rowMeans(C^2, na.rm = TRUE); m3 <- rowMeans(C^3, na.rm = TRUE); m4 <- rowMeans(C^4, na.rm = TRUE)
  cbind(skew = m3 / m2^1.5, exkurt = m4 / m2^2 - 3)
}
q3 <- function(v) stats::quantile(v, c(.1, .5, .9), na.rm = TRUE, names = FALSE)

# One-row data frame of diagnostics. sparsity: per-gene fraction of zero values (or NULL).
row_null_diag <- function(D, Dp, label, sparsity = NULL) {
  G <- nrow(D); ut <- upper.tri(D)
  Z <- clr_rowz(D); Zp <- clr_rowz(Dp)
  mr <- row_moments(D); mp <- row_moments(Dp)
  mu <- rowSums(D) / G; mup <- rowSums(Dp) / G
  sdr <- sqrt(rowSums((D - mu)^2) / (G - 1)); sdp <- sqrt(rowSums((Dp - mup)^2) / (G - 1))
  zo <- offd(Z); zpo <- offd(Zp)
  clipped <- rowMeans(zo < 0, na.rm = TRUE)
  Sp <- stouffer(Zp); S <- stouffer(Z)
  ex1 <- function(z, c) mean(z > c, na.rm = TRUE)
  out <- data.frame(label = label, G = G,
    skew_med = stats::median(mr[, 1]), skew_q90 = q3(mr[, 1])[3], exkurt_med = stats::median(mr[, 2]), exkurt_q90 = q3(mr[, 2])[3],
    skew_null_med = stats::median(mp[, 1]), exkurt_null_med = stats::median(mp[, 2]),
    frac_clipped_med = stats::median(clipped), frac_clipped_q10 = q3(clipped)[1], frac_clipped_q90 = q3(clipped)[3],
    mean_ratio_real_null_med = stats::median(mu / mup), mean_ratio_q90 = q3(mu / mup)[3],
    sd_ratio_real_null_med = stats::median(sdr / sdp), sd_ratio_q90 = q3(sdr / sdp)[3])
  for (c in c(2, 3, 4)) {          # one-sided row z under the permutation null vs N(0, 1)
    out[[sprintf("null_z_gt%d", c)]] <- ex1(zpo, c)
    out[[sprintf("gauss_gt%d", c)]] <- stats::pnorm(c, lower.tail = FALSE)
    out[[sprintf("null_S_gt%d", c)]] <- mean(Sp[ut] > c)
    out[[sprintf("real_S_gt%d", c)]] <- mean(S[ut] > c)
  }
  out$null_S_sd <- stats::sd(Sp[ut]); out$null_S_mean <- mean(Sp[ut])
  out$null_cor_zij_zji <- stats::cor(Zp[ut], t(Zp)[ut])
  pg <- rowMeans(zpo > 3, na.rm = TRUE)            # per-gene null exceedance at z > 3
  out$null_gene_gt3_q10 <- q3(pg)[1]; out$null_gene_gt3_med <- q3(pg)[2]; out$null_gene_gt3_q90 <- q3(pg)[3]
  out$null_gene_gt3_max <- max(pg)
  if (!is.null(sparsity)) {
    out$sparsity_med <- stats::median(sparsity); out$frac_genes_ge50pct_zero <- mean(sparsity >= 0.5)
    out$spearman_null_gt3_vs_sparsity <- suppressWarnings(stats::cor(pg, sparsity, method = "spearman"))
    out$spearman_sdratio_vs_sparsity <- suppressWarnings(stats::cor(sdr / sdp, sparsity, method = "spearman"))
  }
  out$real_gene_gt3_med <- stats::median(rowMeans(zo > 3, na.rm = TRUE))
  out
}

# Label-free alternative calibrations of a SIGNED correlation matrix R (zero diagonal ignored).
# Every alternative ends with clip-at-0 and the Stouffer pair combination, as CLR-NS does.
offrow <- function(M, f) { X <- offd(M); apply(X, 1, function(v) f(v[!is.na(v)])) }
alt_scores <- function(R, mi = NULL) {
  G <- nrow(R); A <- abs(R); diag(A) <- 0
  out <- list()
  out$abs_r <- A
  out$clr_ns <- stouffer(clr_rowz(A))                                     # reference (equals clr_calibrate)
  Fz <- abs(atanh(pmin(pmax(R, -1 + 1e-12), 1 - 1e-12))); diag(Fz) <- 0
  out$fisher <- stouffer(clr_rowz(Fz))                                    # Fisher z, then |.|, then CLR
  med <- offrow(A, stats::median); mad <- offrow(A, stats::mad)
  mad[mad <= 0] <- 1
  out$robust <- stouffer((A - med) / mad)                                 # row median / MAD
  dl <- offrow(R, stats::median); sg <- offrow(R, function(v) stats::IQR(v) / 1.349); sg[sg <= 0] <- 1
  P2 <- 2 * stats::pnorm(-abs(R - dl) / sg)                               # Efron-style central matching on signed r
  out$efron <- stouffer(stats::qnorm(P2, lower.tail = FALSE))
  s2 <- offrow(R, function(v) mean(v^2))                                  # half-normal (chi_1) row fit, MLE
  Ph <- 2 * stats::pnorm(-A / sqrt(s2))
  out$halfnormal <- stouffer(stats::qnorm(Ph, lower.tail = FALSE))
  out$rayleigh_absr <- clr::clr_calibrate(A, method = "rayleigh")         # clr.m 'rayleigh' on |r|
  Rk <- offd(A); Rk <- t(apply(Rk, 1, function(v) { o <- v; k <- !is.na(v)
    o[k] <- stats::qnorm((rank(v[k]) - 3 / 8) / (sum(k) + 1 / 4)); o }))
  out$rowrank <- stouffer(Rk)                                             # per-row rank to normal
  if (!is.null(mi)) {
    out$clr_mi <- clr::clr_calibrate(mi, method = "normal", combine = "stouffer")
    out$rayleigh_mi <- clr::clr_calibrate(mi, method = "rayleigh")       # clr.m 'rayleigh' on MI
  }
  for (k in names(out)) { s <- out[[k]]; s[!is.finite(s)] <- 0; diag(s) <- 0; out[[k]] <- s }
  out
}

# How much the pair ranking changes relative to a reference score matrix.
rank_change <- function(S, Sref, label, ref = "clr_ns") {
  ut <- upper.tri(S); a <- S[ut]; b <- Sref[ut]; n1 <- ceiling(0.01 * length(a))
  ta <- order(-a)[seq_len(n1)]; tb <- order(-b)[seq_len(n1)]
  G <- nrow(S)
  top10 <- function(M) lapply(seq_len(G), function(i) { r <- M[i, ]; r[i] <- -Inf; order(-r)[1:10] })
  A10 <- top10(S); B10 <- top10(Sref)
  data.frame(label = label, method = NA, reference = ref, spearman = suppressWarnings(stats::cor(a, b, method = "spearman")),
             jaccard_top1pct = length(intersect(ta, tb)) / length(union(ta, tb)),
             row_top10_overlap = mean(vapply(seq_len(G), function(i) length(intersect(A10[[i]], B10[[i]])) / 10, 0)),
             ties_at_top_frac = mean(a[ta] == max(a)))
}
