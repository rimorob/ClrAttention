#!/usr/bin/env Rscript
# Reachability-confinement diagnostic (docs/plans/fable_sim_reading.md section 2) on one N = 400 Gaussian/identity compendium
# (seed 1, regenerated: matrices are not saved). Synthetic only. Usage (repo root): Rscript analysis/sim_planted/sim_reach_diag.R [seed]
args <- commandArgs(TRUE); s <- if (length(args)) as.integer(args[1]) else 1L
source("analysis/sim_planted/sim_lib.R"); U <- sim_ut(); G <- SIM_G
st <- sim_structure(s); tr <- sim_truth(st, U); X <- sim_data(st, 400L, "gaussian", "identity", 2L)
fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = 4L)$calibrate(method = "normal", combine = "stouffer")
fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A <- fit$operator; fit$release(); P <- sim_P(A)
Pt <- as.matrix(P); for (t in 2:9) Pt <- step_walk(Pt, P)
# robust (registered) tcaCLR
pr <- tca_parts(Pt); Sr <- stouffer0(pr$Z)
# non-robust a3: floor log(1e-300), row statistics over all off-diagonal entries
L <- log(pmax(Pt, 1e-300)); diag(L) <- NA; mu <- rowMeans(L, na.rm = TRUE); sg <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (G - 2))
Za <- (L - mu) / sg; Za[is.na(Za) | Za < 0] <- 0; Sa <- stouffer0(Za)
Lr <- log(Pt); Lr[!is.finite(Lr)] <- NA; diag(Lr) <- NA
rows <- list()
for (k in which(st$sizes >= 70L)) { m <- st$members[[k]]; n <- length(m)
  reach <- !is.na(Lr[m, ]); inm <- reach[, m, drop = FALSE]; inm_n <- rowSums(inm)           # reachable entries (off-diagonal) and in-module ones
  zm <- function(Z) mean(Z[m, m][upper.tri(diag(n))] * 1)                                       # mean one-sided Z over member entries (both directions averaged)
  w <- tr$modpos[[k]]; o <- tr$outpos[[k]]; ap <- function(S) { v <- S[U$UT]; ap_fast(c(v[w], v[o]), c(rep(TRUE, length(w)), rep(FALSE, length(o)))) }
  win <- function(S) mean(S[U$UT][w])
  rows[[length(rows) + 1L]] <- data.frame(size = n, lambda = st$lam[k], reachable_offdiag = mean(rowSums(reach)), in_module_reachable = mean(inm_n), module_size_minus1 = n - 1,
    inmodule_share = mean(inm_n / rowSums(reach)), sigma_log_P9_robust = mean(pr$sigma[m]), sigma_all_entries = mean(sg[m]),
    mean_logP_member_entries = mean(Lr[m, m], na.rm = TRUE), row_mean_logP = mean(mu[m]), member_oneZ_robust = zm(pr$Z), member_oneZ_all = zm(Za),
    member_S_robust = win(Sr), member_S_all = win(Sa), AP_robust = ap(Sr), AP_all = ap(Sa)) }
d <- do.call(rbind, rows); d <- d[order(d$size, d$lambda), ]
cat(sprintf("seed %d, N = 400; unreachable share of all pairs %.4f; rows with fewer than 2 reachable: %d\n", s, 1 - sum(pr$nreach) / (G * (G - 1)), sum(pr$nreach < 2)))
print(format(d, digits = 3), row.names = FALSE)
cat(sprintf("all-modules reference: median reachable off-diagonal entries per row %d of %d\n", as.integer(median(pr$nreach)), G - 1))
write.csv(d, sprintf("results/sim_planted/addendum/reach_diag_N400_s%02d.csv", s), row.names = FALSE)
