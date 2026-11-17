# Planted-module simulation (docs/plans/fable_theory_memo.md section 5.4; choices and scoring rules fixed in
# docs/plans/sim_planted_choices.md before any run). Synthetic data only. No method parameter is tuned: keff 10, cap 50,
# alpha 0.5, t = 9 (wide arm keff 40, cap 200). Sourced by sim_run.R and sim_summarize.R from the repo root.
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R")               # ap_fast (the project's AUPR)
source("analysis/robust_z/robust_z_lib.R")   # rz_lazy_P, rz_stationary, rz_hub_stratify

SIM_G <- 3000L
SIM_SIZES <- c(5L, 10L, 20L, 40L, 70L, 100L)
SIM_COUNTS <- c(40L, 20L, 10L, 5L, 3L, 2L)

## ---- generator ---------------------------------------------------------------------------------------------------
sim_structure <- function(s, G = SIM_G) {
  set.seed(20261117L + s)
  sizes <- rep(SIM_SIZES, SIM_COUNTS); nmod <- length(sizes)               # 80 modules
  perm <- sample(G); pos <- 0L
  members <- vector("list", nmod); lam <- numeric(nmod)
  kin <- ave(seq_len(nmod), sizes, FUN = seq_along)                         # index of the module within its size
  for (k in seq_len(nmod)) { members[[k]] <- perm[pos + seq_len(sizes[k])]; pos <- pos + sizes[k]
    lam[k] <- c(0.55, 0.70, 0.85)[((kin[k] - 1L + s - 1L) %% 3L) + 1L] }
  big <- which(sizes >= 20L)
  sats <- lapply(big, function(k) { v <- perm[pos + 1:3]; pos <<- pos + 3L; v }); names(sats) <- big
  hubs <- perm[pos + 1:50]; pos <- pos + 50L
  regs <- perm[pos + 1:30]; pos <- pos + 30L
  pool <- perm[(pos + 1L):G]
  nreg <- 30L; mm <- c(5L, 15L, 40L)[((seq_len(nreg) - 1L) %% 3L) + 1L]
  bb <- c(0.35, 0.50)[(((seq_len(nreg) - 1L) %/% 3L) %% 2L) + 1L]
  cnt <- integer(G); targets <- vector("list", nreg)
  for (r in seq_len(nreg)) { av <- pool[cnt[pool] < 2L]; tg <- av[sample.int(length(av), mm[r])]
    cnt[tg] <- cnt[tg] + 1L; targets[[r]] <- tg }
  nF <- nmod + nreg + 3L
  L <- matrix(0, G, nF)
  for (k in seq_len(nmod)) L[members[[k]], k] <- lam[k]
  for (i in seq_along(big)) L[sats[[i]], big[i]] <- 0.8 * lam[big[i]]
  for (r in seq_len(nreg)) { L[regs[r], nmod + r] <- 0.6; L[targets[[r]], nmod + r] <- bb[r] }
  for (h in hubs) L[h, sample.int(nmod + nreg, 33L)] <- 0.25
  L[, (nmod + nreg + 1L):nF] <- matrix(stats::rnorm(G * 3L, 0, 0.15), G, 3L)
  sigma_g <- stats::runif(G, 0.3, 1.0)
  ss <- rowSums(L^2); res <- pmax(0.05, 1 - ss); sdg <- sqrt(ss + res)
  list(s = s, G = G, nF = nF, sizes = sizes, lam = lam, members = members, sats = sats, big = big, hubs = hubs,
       regs = regs, targets = targets, reg_m = mm, reg_b = bb, L = L, res = res, sdg = sdg, sigma_g = sigma_g,
       hubsum = rowSums(abs(L / sdg)))
}

sim_data <- function(st, N, noise = c("gaussian", "t4"), marginal = c("identity", "copula"), Nidx) {
  noise <- match.arg(noise); marginal <- match.arg(marginal); G <- st$G
  set.seed(20262000L + 100L * st$s + 10L * Nidx + as.integer(noise == "t4") + 1L)
  F <- matrix(stats::rnorm(st$nF * N), st$nF, N)
  E <- if (noise == "gaussian") matrix(stats::rnorm(G * N), G, N) else matrix(stats::rt(G * N, 4) / sqrt(2), G, N)
  X <- (st$L %*% F + sqrt(st$res) * E) / st$sdg
  if (marginal == "copula") X <- exp(st$sigma_g * X)
  X
}

## ---- pair indexing (upper triangle, column-major as which(upper.tri)) --------------------------------------------
sim_ut <- function(G = SIM_G) {
  UT <- which(upper.tri(matrix(FALSE, G, G))); ij <- arrayInd(UT, c(G, G))
  list(UT = UT, ri = ij[, 1], ci = ij[, 2], np = length(UT))
}
pair_pos <- function(i, j) { lo <- pmin(i, j); hi <- pmax(i, j); as.integer((hi - 1) * (hi - 2) / 2 + lo) }
within_pos <- function(m) { m <- sort(m); ut <- which(upper.tri(diag(length(m))), arr.ind = TRUE); pair_pos(m[ut[, 1]], m[ut[, 2]]) }

# truths and the position sets every metric needs (all in terms of upper-triangle positions)
sim_truth <- function(st, U) {
  np <- U$np; nmod <- length(st$sizes); G <- st$G
  modpos <- lapply(st$members, within_pos)
  bin <- integer(np); binid <- match(st$sizes, SIM_SIZES)
  for (k in seq_len(nmod)) bin[modpos[[k]]] <- binid[k]
  nreg <- length(st$regs)
  regpos <- lapply(seq_len(nreg), function(r) within_pos(c(st$regs[r], st$targets[[r]])))
  strat_id <- (match(st$reg_m, c(5L, 15L, 40L)) - 1L) * 2L + match(st$reg_b, c(0.35, 0.50))      # 1..6 (m, b) strata
  stpos <- lapply(1:6, function(k) sort(unique(unlist(regpos[strat_id == k]))))
  t3pos <- sort(unique(unlist(regpos)))
  t4pos <- sort(unique(unlist(lapply(seq_len(nreg), function(r) pair_pos(rep(st$regs[r], length(st$targets[[r]])), st$targets[[r]])))))
  satpos <- lapply(seq_along(st$big), function(i) { k <- st$big[i]
    pair_pos(rep(st$sats[[i]], each = length(st$members[[k]])), rep(st$members[[k]], times = 3L)) })
  outpos <- lapply(st$members, function(m) { out <- setdiff(seq_len(G), m)
    pair_pos(rep(m, each = length(out)), rep(out, times = length(m))) })
  list(bin = bin, modpos = modpos, binid = binid, t3pos = t3pos, stpos = stpos, strat_id = strat_id, t4pos = t4pos,
       satpos = satpos, outpos = outpos, all_sat = sort(unique(unlist(satpos))))
}

## ---- walk, readouts ----------------------------------------------------------------------------------------------
sim_P <- function(A, alpha = 0.5) rz_lazy_P(A, alpha)
step_walk <- function(Pt, P) as.matrix(Pt %*% P)

# tcaCLR (a3r) reference code, docs/plans/corum_analysis.md; also returns the row SD (sigma_bg) and reachable counts
tca_parts <- function(Pt) {
  L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  sg <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1))
  Z <- (L - mu) / sg; Z[is.na(Z) | Z < 0] <- 0
  list(Z = Z, sigma = sg, nreach = n)
}
stouffer0 <- function(Z) { S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S }
tca_S <- function(Pt) { p <- tca_parts(Pt); stouffer0(p$Z) }
sym0 <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }
# the project's clr_R (corum_followup.R): row z, diagonal excluded, SD denominator ncol - 2, negatives to 0, Stouffer
clr_R <- function(S) { diag(S) <- NA; mu <- rowMeans(S, na.rm = TRUE)
  sdv <- sqrt(rowSums((S - mu)^2, na.rm = TRUE) / (ncol(S) - 2)); Z <- (S - mu) / sdv; Z[is.na(Z) | Z < 0] <- 0
  out <- (Z + t(Z)) / sqrt(2); diag(out) <- 0; out }

sinkhorn_P <- function(P, tol = 1e-9, maxit = 2000L) {
  P <- Matrix::Matrix(P, sparse = TRUE); G <- nrow(P); r <- rep(1, G); cc <- rep(1, G); res <- Inf; it <- 0L
  Pt <- Matrix::t(P)
  while (it < maxit && res > tol) { it <- it + 1L
    cc <- 1 / as.numeric(Pt %*% r); r <- 1 / as.numeric(P %*% cc)
    res <- max(abs(as.numeric(Matrix::t(P) %*% r) * cc - 1)) }
  B <- Matrix::Diagonal(G, r) %*% P %*% Matrix::Diagonal(G, cc)
  list(P = B, iters = it, resid = res)
}

hub_groups <- function(h, k = 10L) {      # the strata of rz_hub_stratify, as integer codes
  g <- cut(h, unique(stats::quantile(h, seq(0, 1, length.out = k + 1L))), include.lowest = TRUE)
  split(seq_along(h), g)
}
stratify_vec <- function(v, groups) {
  v[!is.finite(v)] <- 0; r <- numeric(length(v))
  for (ix in groups) r[ix] <- rank(v[ix]) / length(ix)
  r
}

## ---- metrics from one upper-triangle score vector ----------------------------------------------------------------
auc_rank <- function(x, y) {              # P(x > y) + 0.5 P(x = y)
  r <- rank(c(x, y)); (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y)) }

sim_metrics <- function(v, tr, U, st, hubsets, groups = NULL, full = TRUE) {
  np <- U$np; v[!is.finite(v)] <- 0
  out <- list(); add <- function(metric, key, value) out[[length(out) + 1L]] <<- data.frame(metric = metric, key = as.character(key), value = value, stringsAsFactors = FALSE)
  lab1 <- tr$bin > 0L; nm <- sum(lab1)
  add("aupr_T1", "all", ap_fast(v, lab1)); add("base_T1", "all", nm / np)
  q <- ceiling(0.01 * np); top <- order(v, decreasing = TRUE, method = "radix")[seq_len(q)]
  ish <- hubsets$pi; add("hubshare_pi", "top1", mean(ish[U$ri[top]] | ish[U$ci[top]]))
  ish <- hubsets$planted; add("hubshare_planted", "top1", mean(ish[U$ri[top]] | ish[U$ci[top]]))
  add("top1_module_share", "top1", mean(lab1[top]))
  if (!full) return(do.call(rbind, out))
  for (b in seq_along(SIM_SIZES)) { sel <- which(tr$bin == 0L | tr$bin == b)
    add("aupr_T2", SIM_SIZES[b], ap_fast(v[sel], tr$bin[sel] == b)); add("base_T2", SIM_SIZES[b], sum(tr$bin == b) / length(sel)) }
  lab3 <- logical(np); lab3[tr$t3pos] <- TRUE
  add("aupr_T3", "all", ap_fast(v, lab3)); add("base_T3", "all", mean(lab3))
  for (k in 1:6) { excl <- setdiff(tr$t3pos, tr$stpos[[k]]); keep <- rep(TRUE, np); keep[excl] <- FALSE
    l <- logical(np); l[tr$stpos[[k]]] <- TRUE
    add("aupr_T3s", paste0("m", c(5, 15, 40)[(k - 1) %/% 2 + 1], "_b", c(0.35, 0.50)[(k - 1) %% 2 + 1]), ap_fast(v[keep], l[keep])) }
  lab4 <- logical(np); lab4[tr$t4pos] <- TRUE
  add("aupr_T4", "all", ap_fast(v, lab4)); add("base_T4", "all", mean(lab4))
  # macro per-module AP: within-module pairs against member-to-outside pairs
  for (k in seq_along(tr$modpos)) { w <- tr$modpos[[k]]; o <- tr$outpos[[k]]
    add("map_module", k, ap_fast(c(v[w], v[o]), c(rep(TRUE, length(w)), rep(FALSE, length(o))))) }
  # mechanism: mean within-module score per module; satellite vs member-member probability
  for (k in seq_along(tr$modpos)) add("mod_mean_score", k, mean(v[tr$modpos[[k]]]))
  for (i in seq_along(tr$satpos)) { k <- st$big[i]; add("sat_auc", k, auc_rank(v[tr$satpos[[i]]], v[tr$modpos[[k]]])) }
  add("sat_share_top1", "top1", sum(top %in% tr$all_sat) / q)
  # hub-stratified AUPR (T1 and T3) under both stratifiers
  if (!is.null(groups)) for (nmg in names(groups)) { r <- stratify_vec(v, groups[[nmg]])
    add(paste0("aupr_T1_hub_", nmg), "all", ap_fast(r, lab1)); add(paste0("aupr_T3_hub_", nmg), "all", ap_fast(r, lab3)) }
  do.call(rbind, out)
}

# per-row mechanism quantities of one walk matrix Pt (dense P^t): sigma_bg, unreachable share, retention
walk_mech <- function(Pt, st, sigma = NULL, nreach = NULL) {
  G <- nrow(Pt)
  if (is.null(sigma)) { p <- tca_parts(Pt); sigma <- p$sigma; nreach <- p$nreach }
  ret <- vapply(seq_along(st$members), function(k) { m <- st$members[[k]]; sum(Pt[m, m]) / length(m) }, 0)
  sbg <- vapply(seq_along(st$members), function(k) stats::median(sigma[st$members[[k]]], na.rm = TRUE), 0)
  list(sigma_median = stats::median(sigma, na.rm = TRUE), sigma_cv = stats::sd(sigma, na.rm = TRUE) / mean(sigma, na.rm = TRUE),
       unreach_share = 1 - sum(nreach) / (G * (G - 1)), rows_lt2_reach = sum(nreach < 2L), retention = ret, sigma_members = sbg)
}
