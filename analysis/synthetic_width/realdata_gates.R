#!/usr/bin/env Rscript
# LABEL-FREE, real expression only (no truth set is read, nothing is scored).
# Round 2 of docs/plans/background_width_exploration.md (section 8.3): agreement
# between the two partner-less-gene calls on M3D and K562.
#   GE  empirical null: per-row Efron central matching on the unclipped Stouffer
#       CLR row; E_i = #{u_ij > delta0 + c sigma0}, c = qnorm(1 - 0.5/(G - 1));
#       open if E_i >= 3.
#   GR  reproducibility over a k-way non-hierarchical cluster cut of the samples
#       (PI decision): k-means (k = 5, nstart 25) on the top 10 PCs of the
#       gene-standardized samples; for each held-out cluster c, repro_c =
#       |top-50 positive CLR of the other clusters & top-50 of cluster c|;
#       S_i = sum_c repro_c; open if S_i >= smallest s with
#       P(Poisson(5 * 50^2 / (G - 1)) >= s) < 0.05 (approximate: folds dependent).
# Clusters separate conditions, not batches: a stricter test than lab halves.
# Usage: Rscript realdata_gates.R <m3d|k562> [threads=8]
suppressPackageStartupMessages(library(clr))
a <- commandArgs(TRUE); DS <- a[1]; THR <- if (length(a) >= 2) as.integer(a[2]) else 8L
stopifnot(DS %in% c("m3d", "k562"))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
if (DS == "m3d") {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(normalizePath(m3d))
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE, stringsAsFactors = FALSE)
} else {
  raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
}
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"; rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
G <- nrow(X); N <- ncol(X); say(sprintf("%s: %d genes x %d samples", DS, G, N))
mi_of <- function(Xm) bspline_mi(Xm, bins = "hg", spline_order = 3, threads = THR)
clr_of <- function(M) clr_calibrate(M, method = "normal", combine = "stouffer")
top50 <- function(C) { T <- matrix(FALSE, G, G)
  for (i in seq_len(G)) { r <- C[i, ]; r[i] <- -Inf; o <- order(r, decreasing = TRUE)[1:50]; T[i, o[r[o] > 0]] <- TRUE }
  T }
efron_E <- function(M) {
  diag(M) <- NA; mu <- rowMeans(M, na.rm = TRUE); s <- sqrt(rowSums((M - mu)^2, na.rm = TRUE) / (G - 2))
  Z <- (M - mu) / s; rm(M); U <- (Z + t(Z)) / sqrt(2); rm(Z); cc <- stats::qnorm(1 - 0.5 / (G - 1))
  vapply(seq_len(G), function(i) {
    u <- U[i, -i]; d <- stats::density(u, n = 512); q <- stats::quantile(u, c(0.25, 0.75), names = FALSE)
    sel <- d$x >= q[1] & d$x <= q[2] & d$y > 0
    b <- if (sum(sel) >= 3) stats::.lm.fit(cbind(1, d$x[sel], d$x[sel]^2), log(d$y[sel]))$coefficients else c(NA, NA, NA)
    if (all(is.finite(b)) && b[3] < 0) { s2 <- -1 / (2 * b[3]); d0 <- b[2] * s2; s0 <- sqrt(s2) } else {
      d0 <- stats::median(u); s0 <- diff(q) / 1.349 }
    sum(u > d0 + cc * s0) }, 0)
}
M <- mi_of(X); E <- efron_E(M); rm(M); invisible(gc(verbose = FALSE)); say("full MI and GE done")
Xs <- (X - rowMeans(X)) / apply(X, 1, stats::sd)
pc <- stats::prcomp(t(Xs), rank. = 10)$x; set.seed(20261031L)
cl <- stats::kmeans(pc, centers = 5, nstart = 25, iter.max = 100)$cluster; rm(Xs, pc)
say("cluster sizes: ", paste(table(cl), collapse = " / "))
S <- numeric(G)
for (k in 1:5) {
  Ta <- top50(clr_of(mi_of(X[, cl != k]))); Tb <- top50(clr_of(mi_of(X[, cl == k])))
  S <- S + rowSums(Ta & Tb); rm(Ta, Tb); invisible(gc(verbose = FALSE)); say("fold ", k, " done") }
lam <- 5 * 50^2 / (G - 1); s0 <- min(which(stats::ppois(0:200 - 1, lam, lower.tail = FALSE) < 0.05)) - 1L
oE <- E >= 3; oR <- S >= s0
po <- mean(oE == oR); pe <- mean(oE) * mean(oR) + (1 - mean(oE)) * (1 - mean(oR))
qs <- function(v) paste(sprintf("%.1f", stats::quantile(v, c(0.1, 0.5, 0.9))), collapse = " / ")
res <- data.frame(dataset = DS, G = G, N = N, clusters = paste(table(cl), collapse = "/"),
  GE_E_q10_med_q90 = qs(E), GE_open = mean(oE), GE_E0 = mean(E == 0),
  GR_S_q10_med_q90 = qs(S), GR_threshold = s0, GR_open = mean(oR), GR_S0 = mean(S == 0),
  both_open = mean(oE & oR), GE_only = mean(oE & !oR), GR_only = mean(!oE & oR), both_closed = mean(!oE & !oR),
  kappa = (po - pe) / (1 - pe), spearman_E_S = stats::cor(E, S, method = "spearman"))
dir.create("analysis/synthetic_width/out/round2", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, sprintf("analysis/synthetic_width/out/round2/realdata_gates_%s.csv", DS), row.names = FALSE)
print(t(res)); say("done")
