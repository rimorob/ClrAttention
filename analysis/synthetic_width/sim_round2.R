#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Round 2 of the background-width exploration
# (docs/plans/background_width_exploration.md, section 8; rule fixed before
# running). One run = one (generator, partner-less fraction f, batch strength b,
# seed). Writes metrics and gate diagnostics for that run.
# Usage: Rscript sim_round2.R <gen: std|regheavy> <f: 0.05..0.5> <b: 0|1|2> <seed>
args <- commandArgs(TRUE)
GEN <- args[1]; FPL <- as.numeric(args[2]); BAT <- as.numeric(args[3]); SEED <- as.integer(args[4])
stopifnot(GEN %in% c("std", "regheavy"))
src <- readLines("analysis/synthetic_width/sim_width.R")
src <- src[seq_len(grep("^## ---- run", src) - 1L)]
src <- src[!grepl("^args <-|^SEEDS <-|^SCEN <-|^NSCALE <-|^OUT <-|^dir.create", src)]
SCEN <- "messy"; eval(parse(text = src))     # G = 2000, N = 500 and all helpers
OUT <- "analysis/synthetic_width/out/round2/runs"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TAG <- sprintf("%s_f%02d_b%d_s%d", GEN, round(100 * FPL), BAT, SEED)

simulate2 <- function(seed) {
  set.seed(seed)
  facs <- list(); load <- vector("list", G); role <- rep("background", G); csize <- rep(NA_integer_, G)
  add_load <- function(g, f, b) for (k in seq_along(g)) load[[g[k]]] <<- c(load[[g[k]]], stats::setNames(b[k], f))
  nxt <- 1L; take <- function(n) { idx <- nxt:(nxt + n - 1L); nxt <<- nxt + n; idx }
  spec <- c(`3` = 10, `4` = 8, `5` = 8, `6` = 5, `8` = 5, `10` = 4, `12` = 3, `15` = 2, `20` = 3, `30` = 2, `40` = 1)
  target_mod <- round((1 - FPL) * G) - 105L
  if (GEN == "std") {
    unit <- sum(spec * as.integer(names(spec))) + 12 + 500 + 30; m <- target_mod / unit
    cnt <- spec; cnt[] <- pmax(1, round(spec * m)); n_reg <- max(1, round(12 * m)); n_hub <- round(30 * m)
    cplx_genes <- sum(cnt * as.integer(names(spec))); pool_n <- target_mod - cplx_genes - n_reg - n_hub
    rsz <- rep_len(c(10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150), n_reg)
    tl <- c(0.2, 0.5); rl <- c(0.5, 0.8)
  } else {
    cnt <- spec; cnt[] <- pmax(1, round(spec * 0.5)); n_reg <- 30L; n_hub <- 30L
    cplx_genes <- sum(cnt * as.integer(names(spec))); pool_n <- target_mod - cplx_genes - n_reg - n_hub
    rsz <- round(exp(stats::runif(n_reg, log(50), log(400))))
    tl <- c(0.15, 0.45); rl <- c(0.4, 0.8)
  }
  cliques <- list(); cinfo <- data.frame()
  for (s in names(spec)) for (r in seq_len(cnt[[s]])) {
    n <- as.integer(s); lam <- stats::runif(1, 0.35, 0.85); nm <- sprintf("C%d_%d", n, r)
    g <- take(n); facs[[nm]] <- stats::rnorm(N); add_load(g, nm, rep(lam, n))
    role[g] <- "complex"; csize[g] <- n; cliques[[nm]] <- g
    cinfo <- rbind(cinfo, data.frame(clique = nm, size = n, loading = lam)) }
  g <- take(80); facs$RIBO <- stats::rnorm(N); facs$R60 <- stats::rnorm(N); facs$R40 <- stats::rnorm(N)
  add_load(g, "RIBO", rep(0.75, 80)); add_load(g[1:47], "R60", rep(0.35, 47)); add_load(g[48:80], "R40", rep(0.35, 33))
  role[g] <- "ribosome"; csize[g] <- 80L; cliques$RIBO <- g
  cinfo <- rbind(cinfo, data.frame(clique = "RIBO", size = 80, loading = 0.75))
  sat <- take(25); add_load(sat, "RIBO", rep(0.55, 25)); role[sat] <- "satellite"
  tfs <- take(n_reg); pool <- take(pool_n)
  part <- split(sample(pool), rep_len(seq_len(n_reg), pool_n))   # every pool gene in >= 1 regulon
  regulons <- list(); edges <- data.frame()
  for (r in seq_len(n_reg)) { nm <- sprintf("TF%d", r); facs[[nm]] <- stats::rnorm(N)
    add_load(tfs[r], nm, stats::runif(1, rl[1], rl[2])); role[tfs[r]] <- "regulator"
    extra <- max(0, rsz[r] - length(part[[r]]))
    tg <- sort(unique(c(part[[r]], sample(setdiff(pool, part[[r]]), min(extra, pool_n - length(part[[r]]))))))
    add_load(tg, nm, stats::runif(length(tg), tl[1], tl[2]))
    regulons[[nm]] <- c(tfs[r], tg); edges <- rbind(edges, data.frame(tf = tfs[r], target = tg)) }
  role[pool] <- "target"
  hubs <- if (n_hub > 0) take(n_hub) else integer(0); fn <- names(facs)
  for (h in hubs) for (ff in sample(fn, 5)) add_load(h, ff, 0.30)
  role[hubs] <- "hub"
  facs$GLOBAL <- stats::rnorm(N); gl <- sample(G, round(0.3 * G)); add_load(gl, "GLOBAL", rep(0.20, length(gl)))
  for (k in 1:2) { nm <- paste0("PROG", k); facs[[nm]] <- stats::rnorm(N); add_load(seq_len(G), nm, stats::rnorm(G, 0, 0.25)) }
  Fm <- do.call(rbind, facs); X <- matrix(0, G, N)
  for (i in seq_len(G)) { b <- load[[i]]; e <- stats::rt(N, 4) / sqrt(2)
    if (sum(b^2) > 0.9) b <- b * sqrt(0.9 / sum(b^2))
    X[i, ] <- colSums(b * Fm[names(b), , drop = FALSE]) + sqrt(1 - sum(b^2)) * e }
  lab <- rep_len(1:10, N)
  if (BAT > 0) { sb <- stats::runif(G, 0, 0.8 * BAT); D <- matrix(stats::rnorm(G * 10), G, 10) * sb; X <- X + D[, lab] }
  for (i in sample(G, round(0.2 * G))) { q <- stats::quantile(X[i, ], stats::runif(1, 0.3, 0.6)); X[i, X[i, ] < q] <- q }
  X <- (X - rowMeans(X)) / apply(X, 1, stats::sd)
  # "background" = loads only on GLOBAL/PROG (no module factor): partner-less
  list(X = X, role = role, csize = csize, cliques = cliques, cinfo = cinfo, regulons = regulons,
       edges = edges, tfs = tfs, sat = sat, lab = lab)
}

# Efron central matching per row of the unclipped Stouffer CLR scores.
efron_gate <- function(M) {
  Mz <- M; diag(Mz) <- NA
  mu <- rowMeans(Mz, na.rm = TRUE); s <- sqrt(rowSums((Mz - mu)^2, na.rm = TRUE) / (G - 2))
  Z <- (Mz - mu) / s; U <- (Z + t(Z)) / sqrt(2)
  cc <- stats::qnorm(1 - 0.5 / (G - 1))
  out <- t(vapply(seq_len(G), function(i) {
    u <- U[i, -i]; d <- stats::density(u, n = 512); q <- stats::quantile(u, c(0.25, 0.75), names = FALSE)
    sel <- d$x >= q[1] & d$x <= q[2] & d$y > 0
    b <- if (sum(sel) >= 3) stats::.lm.fit(cbind(1, d$x[sel], d$x[sel]^2), log(d$y[sel]))$coefficients else c(NA, NA, NA)
    if (all(is.finite(b)) && b[3] < 0) { s2 <- -1 / (2 * b[3]); d0 <- b[2] * s2
      pi0 <- min(1, exp(b[1] + d0^2 / (2 * s2) + log(sqrt(2 * pi * s2)))); s0 <- sqrt(s2)
    } else { d0 <- stats::median(u); s0 <- diff(q) / 1.349; pi0 <- 1 }
    E <- sum(u > d0 + cc * s0)
    r <- u > d0; f <- stats::approx(d$x, d$y, u[r], rule = 2)$y
    lsum <- sum(1 - ifelse(f > 0, pmin(1, pi0 * stats::dnorm(u[r], d0, s0) / f), 1))
    c(E = E, lfdr_sum = lsum) }, c(E = 0, lfdr_sum = 0)))
  out
}
auc <- function(pos, neg) { if (!length(pos) || !length(neg)) return(NA_real_); r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg)) }

t0 <- Sys.time()
sim <- simulate2(20261031L + SEED); X <- sim$X; tr <- make_truth(sim)
istf <- seq_len(G) %in% sim$tfs; Ed <- matrix(FALSE, G, G); Ed[cbind(sim$edges$tf, sim$edges$target)] <- TRUE
Ed <- Ed | t(Ed); ii <- (tr$idx - 1L) %% G + 1L; jj <- (tr$idx - 1L) %/% G + 1L
tr$tfpair <- xor(istf[ii], istf[jj]); tr$edge <- Ed[tr$idx]; rm(Ed)
M <- mi_of(X); C <- clr_of(M); eg <- efron_gate(M); rm(M)
hA <- sim$lab <= 5L
TA <- top50(clr_of(mi_of(X[, hA]))); TB <- top50(clr_of(mi_of(X[, !hA]))); rep_n <- rowSums(TA & TB); rm(TA, TB)
lambda0 <- 50 * 50 / (G - 1); k0 <- min(which(stats::ppois(0:50 - 1, lambda0, lower.tail = FALSE) < 0.05)) - 1L
rho_R <- as.numeric(rep_n >= k0); rho_E <- as.numeric(eg[, "E"] >= 3); rho_O <- as.numeric(sim$role != "background")
K10 <- selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L)); RC <- recip(K10)
ops <- list(K10 = K10, RECIP = RC, GR = sink2(K10, rho_R), GR_RECIP = sink2(RC, rho_R),
            GE = sink2(K10, rho_E), GE_RECIP = sink2(RC, rho_E), ORACLE = sink2(K10, rho_O), ORACLE_RECIP = sink2(RC, rho_O))
extra <- c("K10", "GR_RECIP", "GE_RECIP")
score2 <- function(S) { b <- score(S, tr); s <- S[tr$idx]; s[!is.finite(s)] <- -Inf; k <- tr$tfpair
  o <- order(s[k], decreasing = TRUE, method = "radix")
  b$edge_x <- ap_sorted(s[k][o], tr$edge[k][o]) / mean(tr$edge[k]); b }
rows <- list(); add <- function(name, S) rows[[length(rows) + 1L]] <<- cbind(op = name, score2(S))
add("CLR", C)
for (nm in names(ops)) { P9 <- walk9(lazy(ops[[nm]])); add(paste0(nm, "|att_t9"), att(P9))
  if (nm %in% extra) { add(paste0(nm, "|a3r_t9"), a3r(P9)); add(paste0(nm, "|cos_t9"), cosr(P9)) } }
res <- cbind(gen = GEN, f = FPL, b = BAT, seed = SEED, do.call(rbind, rows))
utils::write.csv(res, file.path(OUT, paste0(TAG, "_metrics.csv")), row.names = FALSE)
pl <- sim$role == "background"; sm <- sim$role == "complex" & sim$csize %in% 3:6
roles <- ifelse(sim$role == "complex", ifelse(sim$csize <= 6, "complex_3-6", ifelse(sim$csize <= 15, "complex_8-15", "complex_20-40")), sim$role)
gd <- data.frame(gen = GEN, f = FPL, b = BAT, seed = SEED, f_real = mean(pl), k0 = k0,
  auc_repro = auc(rep_n[sm], rep_n[pl]), auc_efronE = auc(eg[sm, "E"], eg[pl, "E"]),
  auc_efron_lfdr = auc(eg[sm, "lfdr_sum"], eg[pl, "lfdr_sum"]),
  spearman_E_repro = stats::cor(eg[, "E"], rep_n, method = "spearman"),
  agree_open = mean(rho_E == rho_R), kappa = { po <- mean(rho_E == rho_R); pe <- mean(rho_E) * mean(rho_R) + (1 - mean(rho_E)) * (1 - mean(rho_R)); (po - pe) / (1 - pe) })
op <- do.call(rbind, lapply(split(seq_len(G), roles), function(g) data.frame(gen = GEN, f = FPL, b = BAT, seed = SEED,
  role = roles[g[1]], n = length(g), open_GR = mean(rho_R[g]), open_GE = mean(rho_E[g]),
  E_median = stats::median(eg[g, "E"]), lfdr_sum_median = stats::median(eg[g, "lfdr_sum"]), repro_median = stats::median(rep_n[g]))))
utils::write.csv(gd, file.path(OUT, paste0(TAG, "_gates.csv")), row.names = FALSE)
utils::write.csv(op, file.path(OUT, paste0(TAG, "_open.csv")), row.names = FALSE)
cat(sprintf("%s done in %.0f s; f_real %.3f; AUC repro %.3f efronE %.3f\n", TAG,
            as.numeric(difftime(Sys.time(), t0, units = "secs")), gd$f_real, gd$auc_repro, gd$auc_efronE))
