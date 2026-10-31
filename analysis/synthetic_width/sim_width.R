#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Exploration of the PI's "background-adaptive attention width"
# idea (docs/plans/background_width_exploration.md). No real expression data,
# truth set or result file is read.
#
# Question: if each gene's attention width followed how many partners rise
# above its (non-Gaussian, gene- and pair-specific) MI background, would
# complexes of size 5, 20 and 80 improve without costing regulons? And is
# "width" the right lever at all, as opposed to the budget (where non-signal
# mass goes) or the readout scale?
#
# Generator (per seed): G = 2,000 genes, N = 500 samples, factor model.
#   - complexes with a CORUM-like size spectrum: sizes 3 (x10), 4 (x8),
#     5 (x8), 6 (x5), 8 (x5), 10 (x4), 12 (x3), 15 (x2), 20 (x3), 30 (x2),
#     40 (x1), loading per complex U(0.35, 0.85) (many weak small complexes,
#     as in K562, where the median within-complex CLR of 5-24-gene complexes
#     is 0.4);
#   - a ribosome of 80 genes: common factor 0.75, plus a 60S (47) / 40S (33)
#     subunit factor at 0.35; it holds about half of all complex pairs, as the
#     RP complexes do in K562 CORUM;
#   - 25 satellites loading 0.55 on the ribosome factor (not members);
#   - 12 regulators (mRNA loading U(0.5, 0.8)) with 10-150 targets each from a
#     pool of 500 (overlapping), target loadings U(0.2, 0.5);
#   - 30 hubs loading 0.30 on 5 random module factors; a weak global factor
#     (0.20 on 30% of genes); the rest is background.
# Scenarios:
#   clean: Gaussian noise, nothing else;
#   messy: the background that broke O2 (permutation null): 10 labs with
#     gene-specific lab-effect SD U(0, 0.8) (a rank-10, gene-specific,
#     non-Gaussian pairwise background), two global programs with loadings
#     N(0, 0.25^2), t4 idiosyncratic noise, and zero-inflation in 20% of genes
#     (values below a gene-specific 30-60% quantile set to the floor).
# MI: the package's B-spline estimator (bins "hg", order 3, one thread);
# CLR: the package's Stouffer CLR. Softmax K10/K25: the package's
# .softmax_rows; per-gene keff vectors: the same bisection (checked equal).
# Lazy walk alpha = 0.5, depth t = 9 (the method's transferred depth; never
# selected on truth here), t = 1 also reported.
#
# Usage: Rscript sim_width.R [seeds=1:3] [scenario=messy] [nscale=0|1]
# Run with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19.
suppressPackageStartupMessages(library(clr))
args <- commandArgs(TRUE)
SEEDS <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:3
SCEN <- if (length(args) >= 2) args[2] else "messy"
NSCALE <- if (length(args) >= 3) as.integer(args[3]) == 1L else FALSE
stopifnot(SCEN %in% c("clean", "messy"))
OUT <- file.path("analysis/synthetic_width/out", SCEN)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
G <- 2000L; N <- 500L
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

## ---- simulation ---------------------------------------------------------------
simulate <- function(seed, n_samp = N) {
  set.seed(seed)
  facs <- list(); load <- vector("list", G)
  role <- rep("background", G); csize <- rep(NA_integer_, G)
  add_load <- function(g, f, b) for (k in seq_along(g)) load[[g[k]]] <<- c(load[[g[k]]], stats::setNames(b[k], f))
  nxt <- 1L; take <- function(n) { idx <- nxt:(nxt + n - 1L); nxt <<- nxt + n; idx }
  spec <- c(`3` = 10, `4` = 8, `5` = 8, `6` = 5, `8` = 5, `10` = 4, `12` = 3, `15` = 2, `20` = 3, `30` = 2, `40` = 1)
  cliques <- list(); cinfo <- data.frame()
  for (s in names(spec)) for (r in seq_len(spec[[s]])) {
    n <- as.integer(s); lam <- stats::runif(1, 0.35, 0.85); nm <- sprintf("C%d_%d", n, r)
    g <- take(n); facs[[nm]] <- stats::rnorm(n_samp); add_load(g, nm, rep(lam, n))
    role[g] <- "complex"; csize[g] <- n; cliques[[nm]] <- g
    cinfo <- rbind(cinfo, data.frame(clique = nm, size = n, loading = lam)) }
  g <- take(80); facs$RIBO <- stats::rnorm(n_samp); facs$R60 <- stats::rnorm(n_samp); facs$R40 <- stats::rnorm(n_samp)
  add_load(g, "RIBO", rep(0.75, 80)); add_load(g[1:47], "R60", rep(0.35, 47)); add_load(g[48:80], "R40", rep(0.35, 33))
  role[g] <- "ribosome"; csize[g] <- 80L; cliques$RIBO <- g
  cinfo <- rbind(cinfo, data.frame(clique = "RIBO", size = 80, loading = 0.75))
  sat <- take(25); add_load(sat, "RIBO", rep(0.55, 25)); role[sat] <- "satellite"
  tfs <- take(12); pool <- take(500); rsz <- c(10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150)
  regulons <- list()
  for (r in seq_along(rsz)) { nm <- sprintf("TF%d", r); facs[[nm]] <- stats::rnorm(n_samp)
    add_load(tfs[r], nm, stats::runif(1, 0.5, 0.8)); role[tfs[r]] <- "regulator"
    tg <- sort(sample(pool, rsz[r])); add_load(tg, nm, stats::runif(length(tg), 0.2, 0.5))
    role[tg] <- "target"; regulons[[nm]] <- c(tfs[r], tg) }
  hubs <- take(30); fn <- names(facs)
  for (h in hubs) for (ff in sample(fn, 5)) add_load(h, ff, 0.30)
  role[hubs] <- "hub"
  facs$GLOBAL <- stats::rnorm(n_samp); gl <- sample(G, round(0.3 * G)); add_load(gl, "GLOBAL", rep(0.20, length(gl)))
  if (SCEN == "messy") for (k in 1:2) { nm <- paste0("PROG", k); facs[[nm]] <- stats::rnorm(n_samp)
    add_load(seq_len(G), nm, stats::rnorm(G, 0, 0.25)) }
  Fm <- do.call(rbind, facs)
  X <- matrix(0, G, n_samp)
  noise <- function() if (SCEN == "messy") stats::rt(n_samp, 4) / sqrt(2) else stats::rnorm(n_samp)
  for (i in seq_len(G)) { b <- load[[i]]
    if (length(b)) { if (sum(b^2) > 0.9) b <- b * sqrt(0.9 / sum(b^2))
      X[i, ] <- colSums(b * Fm[names(b), , drop = FALSE]) + sqrt(1 - sum(b^2)) * noise()
    } else X[i, ] <- noise() }
  sb <- rep(0, G)
  if (SCEN == "messy") {
    lab <- rep_len(1:10, n_samp); sb <- stats::runif(G, 0, 0.8)
    D <- matrix(stats::rnorm(G * 10), G, 10) * sb          # gene x lab shifts
    X <- X + D[, lab]
    drop <- sample(G, round(0.2 * G))
    for (i in drop) { q <- stats::quantile(X[i, ], stats::runif(1, 0.3, 0.6)); X[i, X[i, ] < q] <- q }
  }
  X <- (X - rowMeans(X)) / apply(X, 1, stats::sd)
  list(X = X, role = role, csize = csize, cliques = cliques, cinfo = cinfo, regulons = regulons,
       tfs = tfs, sat = sat, sb = sb, lab = rep_len(1:10, n_samp))
}

## ---- background models (the PI's idea and its relatives) -------------------------
mi_of <- function(X) bspline_mi(X, bins = "hg", spline_order = 3, threads = 1)
clr_of <- function(M) clr_calibrate(M, method = "normal", combine = "stouffer")
ut_idx <- function() which(upper.tri(matrix(0, G, G)))
sym_from_ut <- function(v, idx) { Q <- matrix(0, G, G); Q[idx] <- v; Q + t(Q) }

# Mirror local fdr on a pooled standardized score u (background assumed
# symmetric about its mode on this scale; signal only on the right):
# lfdr(x) = min(1, f(2 m0 - x) / f(x)) for x > m0, else 1. Returns q = 1 - lfdr.
mirror_q <- function(u) {
  d <- stats::density(u, n = 4096, bw = "nrd0")
  m0 <- d$x[which.max(d$y)]
  f <- function(x) stats::approx(d$x, d$y, x, rule = 2, yleft = 0, yright = 0)$y
  fx <- f(u); f0 <- f(2 * m0 - u)
  q <- ifelse(u > m0 & fx > 0, 1 - pmin(1, f0 / fx), 0)
  list(q = q, m0 = m0, pi0 = min(1, 2 * mean(u <= m0)))
}

# EB1: two-way robust background on the cube-root MI scale (Wilson-Hilferty:
# 2N MI is roughly chi-square under independence, so MI^(1/3) is roughly
# Gaussian), y_ij = mu + a_i + a_j + r_ij (symmetric median polish), per-gene
# robust scale s_i = MAD of row i's residuals, bilateral Stouffer of the two
# row-standardized residuals, mirror lfdr pooled over pairs.
eb_additive <- function(M) {
  Y <- M^(1 / 3); diag(Y) <- NA
  mu <- stats::median(Y, na.rm = TRUE); a <- rep(0, G)
  for (it in 1:8) { R <- Y - mu - outer(a, a, "+"); da <- apply(R, 1, stats::median, na.rm = TRUE)
    a <- a + da; mu <- mu + stats::median(a); a <- a - stats::median(a) }
  R <- Y - mu - outer(a, a, "+"); s <- apply(R, 1, stats::mad, na.rm = TRUE)
  U <- (sweep(R, 1, s, "/") + sweep(R, 2, s, "/")) / sqrt(2); diag(U) <- NA
  idx <- ut_idx(); mq <- mirror_q(U[idx])
  Q <- sym_from_ut(mq$q, idx); U[is.na(U)] <- 0
  list(U = U, Q = Q, pi0 = mq$pi0, a = a, s = s)
}
# EBG: "graphon" background, the literal X-by-Y form of the PI's idea. Genes are
# binned (20 quantile bins) by their background level theta_i = row median of
# MI^(1/3); the background of pair (i, j) is the empirical distribution of all
# pairs in cell (bin_i, bin_j) (about 10,000 pairs per cell), standardized by
# the cell's median and MAD; one pooled shape (mirror lfdr) across cells.
# Captures interactions (two batch-heavy genes share more background MI than
# the sum of their row effects predicts).
eb_graphon <- function(M, nb = 20L) {
  Y <- M^(1 / 3); diag(Y) <- NA
  th <- apply(Y, 1, stats::median, na.rm = TRUE)
  b <- cut(rank(th, ties.method = "first"), nb, labels = FALSE)
  idx <- ut_idx(); i <- (idx - 1L) %% G + 1L; j <- (idx - 1L) %/% G + 1L
  y <- Y[idx]; cell <- pmin(b[i], b[j]) * 100L + pmax(b[i], b[j])
  med <- tapply(y, cell, stats::median); md <- tapply(y, cell, stats::mad)
  k <- match(cell, as.integer(names(med))); u <- (y - med[k]) / md[k]
  mq <- mirror_q(u)
  list(U = sym_from_ut(u, idx), Q = sym_from_ut(mq$q, idx), pi0 = mq$pi0, theta = th)
}
# O2 as preregistered (permutation null on the CLR scale, pooled, z > 0,
# Storey pi0, sum of (1 - lfdr) over the cap-50 candidates, clamp 2-30).
o2_keff <- function(C, C0) {
  ut <- upper.tri(C); z <- C[ut]; z0 <- sort(C0[ut]); N0 <- length(z0)
  p <- (1 + N0 - findInterval(z, z0, left.open = TRUE)) / (1 + N0)
  pi0 <- min(1, mean(p > 0.5) / 0.5)
  zp <- z[z > 0]; z0p <- z0[z0 > 0]; hi <- max(zp, z0p)
  d1 <- stats::density(zp, from = 0, to = hi, n = 4096); d0 <- stats::density(z0p, from = 0, to = hi, n = 4096)
  lf <- function(x) { f <- stats::approx(d1$x, d1$y, x, rule = 2)$y; f0 <- stats::approx(d0$x, d0$y, x, rule = 2)$y
    out <- ifelse(f > 0, pmin(1, pi0 * f0 / f), 1); out[x <= 0] <- 1; out }
  keff <- vapply(seq_len(G), function(i) { row <- C[i, ]; row[i] <- -Inf
    top <- order(row, decreasing = TRUE)[1:50]; top <- top[row[top] > 0]
    if (!length(top)) 0 else sum(1 - lf(row[top])) }, 0)
  list(keff = pmin(pmax(keff, 2), 30), pi0 = pi0)
}
# RMT width: eigenvectors of the normal-score correlation matrix with
# eigenvalues above the Marchenko-Pastur edge (1 + sqrt(G/N))^2 and
# localized (inverse participation ratio <= 150 genes). A gene that
# participates (v_ik^2 * IPR_k >= 0.25) in such a spike gets keff = IPR_k - 1;
# others keep 10.
ipr_keff <- function(X) {
  Z <- t(apply(X, 1, function(x) stats::qnorm(rank(x) / (length(x) + 1))))
  e <- eigen(stats::cor(t(Z)), symmetric = TRUE)
  edge <- (1 + sqrt(G / ncol(X)))^2
  sp <- which(e$values > edge); V <- e$vectors[, sp, drop = FALSE]
  ipr <- 1 / colSums(V^4); loc <- ipr <= 150
  keff <- rep(10, G)
  if (any(loc)) { W <- sweep(V[, loc, drop = FALSE]^2, 2, ipr[loc], "*")
    k <- max.col(W); part <- W[cbind(seq_len(G), k)] >= 0.25
    keff[part] <- pmin(pmax(ipr[loc][k[part]] - 1, 2), 100) }
  list(keff = keff, n_spikes = length(sp), n_local = sum(loc))
}

# EVT: per-row non-Gaussian background with an extreme-value tail. For row i,
# y = MI^(1/3) over j != i; the background tail above the row's 90th percentile
# is modelled as exponential (GPD with shape 0), with scale fitted from the
# 90th-95th percentile spacing (a module of n genes is the top n/G of the row:
# 4% for the 80-gene ribosome here, so the fitting window is clean). Expected
# background count above x: 0.1 (G - 1) exp(-(x - q90) / beta). Tail FDR at
# the k-th largest value = expected count / k, made monotone (q-value style);
# q_ij = 1 - FDR. Returns the G x G matrix of q (row-directed).
evt_q <- function(M) {
  Y <- M^(1 / 3); diag(Y) <- NA; Q <- matrix(0, G, G)
  for (i in seq_len(G)) {
    y <- Y[i, ]; y[i] <- -Inf; qs <- stats::quantile(y[-i], c(0.90, 0.95), names = FALSE)
    beta <- (qs[2] - qs[1]) / log(2); o <- order(y, decreasing = TRUE)[1:200]
    en <- 0.1 * (G - 1) * exp(-(y[o] - qs[1]) / beta); fdr <- pmin(1, en / seq_along(o))
    fdr <- rev(cummin(rev(fdr))); Q[i, o] <- 1 - fdr }
  Q
}

## ---- operators ---------------------------------------------------------------------
# Softmax over each row's top-`cap` positive scores with per-row perplexity keff[i]
# (same bisection as clr:::.softmax_rows; rows with <= keff candidates uniform).
softmax_keff <- function(S, keff, cap = 50L) {
  A <- matrix(0, G, G)
  for (i in seq_len(G)) {
    row <- S[i, ]; row[i] <- -Inf
    top <- order(row, decreasing = TRUE)[seq_len(min(cap, G - 1L))]
    top <- top[row[top] > 0]; if (!length(top)) next
    s <- row[top]; k <- keff[i]
    if (length(top) <= k) { A[i, top] <- 1 / length(top); next }
    neff <- function(lt) { w <- exp((s - max(s)) / exp(lt)); w <- w / sum(w)
      exp(-sum(w[w > 0] * log(w[w > 0]))) }
    lo <- log(1e-4); hi <- log(1e4)
    for (it in 1:60) { mid <- (lo + hi) / 2; if (neff(mid) < k) lo <- mid else hi <- mid }
    w <- exp((s - max(s)) / exp((lo + hi) / 2)); A[i, top] <- w / sum(w)
  }
  A
}
selfloop <- function(A) { e <- rowSums(A) == 0; A[cbind(which(e), which(e))] <- 1; A }
# Background sink: row i sends rho_i = min(1, Q_i) of its budget to partners
# (shaped by A) and keeps 1 - rho_i on itself, where Q_i = sum_j q_ij is the
# expected number of above-background partners.
sink <- function(A, Q) { rho <- pmin(1, rowSums(Q)); A <- A * rho; diag(A) <- diag(A) + (1 - rho); A }
# q-attention: weight = posterior probability of being above background,
# divided by max(1, Q_i); the remainder stays on the gene (no forced budget).
qatt <- function(Q) { Qs <- rowSums(Q); A <- Q / pmax(1, Qs); diag(A) <- diag(A) + (1 - rowSums(A)); A }
# Symmetric sink: attention i -> j is scaled by rho_i * rho_j (neither a
# partner-less gene's row nor its column carries mass); the remainder stays on i.
sink2 <- function(A, rho) { B <- A * outer(rho, rho); diag(B) <- 0; diag(B) <- 1 - rowSums(B); B }
# Reproducibility background (IDR-like): split the samples by lab into two
# halves (labs 1-5 vs 6-10), compute CLR in each; candidate j of row i is kept
# only if it is in row i's top-50 positive CLR in BOTH halves. Softmax K over
# the kept candidates of the full-data CLR (rows with <= K kept candidates are
# uniform over them, rows with none keep their mass), so the width adapts
# (min(K, #reproducible partners)) and partner-less rows fall silent.
top50 <- function(C) { T <- matrix(FALSE, G, G)
  for (i in seq_len(G)) { r <- C[i, ]; r[i] <- -Inf; o <- order(r, decreasing = TRUE)[1:50]; T[i, o[r[o] > 0]] <- TRUE }
  T }
lazy <- function(A) { P <- 0.5 * A; diag(P) <- diag(P) + 0.5; P }
walk9 <- function(P) { P2 <- P %*% P; P4 <- P2 %*% P2; P8 <- P4 %*% P4; P8 %*% P }
# readouts
att <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }
# Profile cosine: cosine similarity of the two genes' t-step walk distributions
# (rows of P^t). Bounded in [0, 1]; two members of a mixed clique of any size n
# have nearly the same distribution (uniform on the clique), so cos -> 1
# independently of n ("inherit, don't divide"). (A first version divided the
# symmetrized P^t by its diagonal; that is not a valid cosine for a
# non-reversible walk and was replaced before the multi-seed runs.)
cosr <- function(Pt) { Rn <- Pt / sqrt(rowSums(Pt^2)); R <- tcrossprod(Rn); diag(R) <- 0; R }
# Reciprocity sink: keep attention i -> j only if j also attends to i at t = 1;
# the removed (unreturned) mass stays on i. Label-free "background" = attention
# that is not returned (one-sided satellites, rows of genes without partners).
recip <- function(A) { B <- A * (t(A) > 0); diag(B) <- 0; B <- selfloop(B)
  diag(B) <- diag(B) + pmax(0, 1 - rowSums(B)); B }
zna <- function(L) { mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  s <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1)); Z <- (L - mu) / s; Z[!is.finite(Z) | Z < 0] <- 0; Z }
a3r <- function(Pt) { L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA; Z <- zna(L)
  S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S }

## ---- metrics -------------------------------------------------------------------------
ap_sorted <- function(s, l) { n <- length(s); if (!n || !sum(l)) return(NA_real_)
  last <- c(s[-1L] != s[-n], TRUE); tp <- cumsum(l)[last]; k <- which(last)
  sum(diff(c(0, tp)) * tp / k) / sum(l) }
auroc <- function(pos, neg) { r <- rank(c(pos, neg)); np <- length(pos); nn <- length(neg)
  (sum(r[seq_len(np)]) - np * (np + 1) / 2) / (np * nn) }
make_truth <- function(sim) {
  idx <- ut_idx(); i <- (idx - 1L) %% G + 1L; j <- (idx - 1L) %/% G + 1L
  cl <- integer(G); for (k in seq_along(sim$cliques)) cl[sim$cliques[[k]]] <- k
  cpx <- cl[i] > 0 & cl[i] == cl[j]
  sz <- ifelse(cpx, sim$csize[i], 0L)
  inribo <- sim$role == "ribosome"
  RM <- matrix(FALSE, G, length(sim$regulons)); for (k in seq_along(sim$regulons)) RM[sim$regulons[[k]], k] <- TRUE
  reg <- rowSums(RM[i, , drop = FALSE] & RM[j, , drop = FALSE]) > 0
  issat <- seq_len(G) %in% sim$sat
  list(idx = idx, cpx = cpx, sz = sz, noribo = !(inribo[i] | inribo[j]), reg = reg,
       ribo_pair = inribo[i] & inribo[j], sat_ribo = (issat[i] & inribo[j]) | (issat[j] & inribo[i]))
}
score <- function(S, tr) {
  s <- S[tr$idx]; s[!is.finite(s)] <- -Inf
  o <- order(s, decreasing = TRUE, method = "radix"); so <- s[o]
  apm <- function(lab, keep = NULL) { if (is.null(keep)) keep <- rep(TRUE, length(lab))
    k <- keep[o]; ap_sorted(so[k], lab[o][k]) / mean(lab[keep]) }
  neg <- !tr$cpx
  data.frame(cpx_x = apm(tr$cpx), cpx_noribo_x = apm(tr$cpx, tr$noribo),
    s5_x = apm(tr$sz == 5L, neg | tr$sz == 5L), s20_x = apm(tr$sz == 20L, neg | tr$sz == 20L),
    s80_x = apm(tr$sz == 80L, neg | tr$sz == 80L),
    small_x = apm(tr$sz %in% 3:6, neg | tr$sz %in% 3:6), mid_x = apm(tr$sz %in% 8:15, neg | tr$sz %in% 8:15),
    reg_x = apm(tr$reg), sat_auc = auroc(s[tr$ribo_pair], s[tr$sat_ribo]))
}

## ---- run -------------------------------------------------------------------------------
res <- list(); wid <- list(); diagn <- list()
for (seed in SEEDS) {
  t0 <- Sys.time(); sd0 <- 20261030L + seed
  sim <- simulate(sd0); X <- sim$X; tr <- make_truth(sim)
  M <- mi_of(X); C <- clr_of(M)
  set.seed(sd0 + 1000L); C0 <- clr_of(mi_of(t(apply(X, 1, sample))))
  o2 <- o2_keff(C, C0); rm(C0)
  e1 <- eb_additive(M); eg <- eb_graphon(M); ip <- ipr_keff(X)
  Qg <- eg$Q; Qs <- rowSums(Qg); Qe <- evt_q(M); Qes <- rowSums(Qe)
  hA <- sim$lab <= 5L
  TA <- top50(clr_of(mi_of(X[, hA]))); TB <- top50(clr_of(mi_of(X[, !hA])))
  Crep <- C * (TA & TB); rep_n <- rowSums(TA & TB); rm(TA, TB)
  oracle <- ifelse(is.na(sim$csize), 10, sim$csize - 1)
  oracle_small <- ifelse(is.na(sim$csize), 10, pmin(10, sim$csize - 1))
  isbg <- sim$role == "background"
  keff_eb <- pmin(pmax(Qs, 2), 100); keff_evt <- pmin(pmax(Qes, 2), 100)
  K10 <- selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L))
  stopifnot(max(abs(K10 - selfloop(softmax_keff(C, rep(10, G), 50L)))) < 1e-8)
  ops <- list(
    K10 = K10,
    K25 = selfloop(clr:::.softmax_rows(C, keff = 25, cap = 50L)),
    O2_perm_width = selfloop(softmax_keff(C, o2$keff, 50L)),
    # oracle decompositions (not realizable; they bound what width can do)
    ORACLE_width = selfloop(softmax_keff(C, oracle, 100L)),          # own-row width = n - 1
    ORACLE_small = selfloop(softmax_keff(C, oracle_small, 50L)),     # narrow only small complexes
    ORACLE_noise_wide = selfloop(softmax_keff(C, ifelse(isbg, 50, 10), 50L)),  # spread no-partner rows thin
    ORACLE_sink = sink(K10, matrix(as.numeric(!isbg), G, 1)),        # no-partner rows keep their mass
    ORACLE_sink2 = sink2(K10, as.numeric(!isbg)),                     # ... and receive none
    REPRO10 = selfloop(clr:::.softmax_rows(Crep, keff = 10, cap = 50L)),
    REPRO25 = selfloop(clr:::.softmax_rows(Crep, keff = 25, cap = 50L)),
    # the PI's idea, two background models
    EB_width = selfloop(softmax_keff(C, keff_eb, 100L)),             # graphon (X-by-Y) background
    EVT_width = selfloop(softmax_keff(C, keff_evt, 100L)),           # per-row EVT tail background
    EVT_sink = sink(K10, Qe),
    RMT_width = selfloop(softmax_keff(C, ip$keff, 100L)),
    RECIP = recip(K10),
    q_attention = qatt(Qg))
  rows <- list()
  add <- function(name, S) rows[[length(rows) + 1L]] <<- cbind(op = name, score(S, tr))
  add("absr", abs(stats::cor(t(X))) - diag(G)); add("CLR", C)
  add("EB_graphon_score", eg$U); add("EB_additive_score", e1$U)
  for (nm in names(ops)) {
    P <- lazy(ops[[nm]]); P9 <- walk9(P)
    add(paste0(nm, "|att_t1"), att(P)); add(paste0(nm, "|att_t9"), att(P9))
    add(paste0(nm, "|cos_t9"), cosr(P9)); add(paste0(nm, "|a3r_t9"), a3r(P9))
    # mechanism: per-pair mass and in-complex mass at t = 9 for sizes 5, 20, 80
    for (s in c(5L, 20L, 80L)) for (g in sim$cliques[sim$cinfo$size == s]) {
      B <- P9[g, g]; offd <- B; diag(offd) <- NA
      diagn[[length(diagn) + 1L]] <- data.frame(seed = seed, op = nm, size = s,
        pair_mass_t9 = mean(offd, na.rm = TRUE),
        in_complex_t9 = mean((rowSums(B) - diag(B)) / (1 - diag(P9)[g])),
        A_pair_t1 = { b <- ops[[nm]][g, g]; diag(b) <- NA; mean(b, na.rm = TRUE) }) }
    rm(P, P9)
  }
  r <- do.call(rbind, rows); r$seed <- seed; res[[length(res) + 1L]] <- r
  roleg <- ifelse(sim$role %in% c("complex", "ribosome"), paste0("complex_", sim$csize), sim$role)
  roleg[roleg %in% paste0("complex_", 3:6)] <- "complex_3-6"
  roleg[roleg %in% paste0("complex_", 8:15)] <- "complex_8-15"
  roleg[roleg %in% paste0("complex_", c(20, 30, 40))] <- "complex_20-40"
  wid[[length(wid) + 1L]] <- data.frame(seed = seed, gene = seq_len(G), role = roleg, sb = sim$sb,
    keff_O2 = o2$keff, keff_EB = keff_eb, Q_EB = Qs, Q_EB_additive = rowSums(e1$Q),
    Q_EVT = Qes, keff_RMT = ip$keff, keff_oracle = oracle,
    recip_kept = rowSums(ops$RECIP) - diag(ops$RECIP), repro_n = rep_n)
  say(sprintf("%s seed %d: %.0f s; O2 pi0 %.3f; EB graphon pi0 %.3f, additive pi0 %.3f; RMT spikes %d (localized %d)",
              SCEN, seed, as.numeric(difftime(Sys.time(), t0, units = "secs")), o2$pi0, eg$pi0, e1$pi0,
              ip$n_spikes, ip$n_local))
  rm(M, C, ops, e1, eg, Qg, K10); invisible(gc(verbose = FALSE))
}
res <- do.call(rbind, res); wid <- do.call(rbind, wid); diagn <- do.call(rbind, diagn)
utils::write.csv(res, file.path(OUT, "metrics_by_seed.csv"), row.names = FALSE)
wsum <- do.call(rbind, lapply(split(wid, wid$role), function(d) data.frame(role = d$role[1], n_genes_per_seed = nrow(d) / length(SEEDS),
  O2 = stats::median(d$keff_O2), EB = stats::median(d$keff_EB), EB_q90 = stats::quantile(d$keff_EB, 0.9),
  Q_EB = stats::median(d$Q_EB), Q_EB_additive = stats::median(d$Q_EB_additive),
  Q_EVT = stats::median(d$Q_EVT), Q_EVT_q90 = stats::quantile(d$Q_EVT, 0.9), RMT = stats::median(d$keff_RMT),
  recip_kept = stats::median(d$recip_kept), repro_n = stats::median(d$repro_n),
  repro_zero = mean(d$repro_n == 0), oracle = stats::median(d$keff_oracle))))
utils::write.csv(wsum, file.path(OUT, "widths_by_role.csv"), row.names = FALSE)
if (SCEN == "messy") {
  bg <- wid[wid$role == "background", ]
  sp <- function(v) stats::cor(bg[[v]], bg$sb, method = "spearman")
  cat(sprintf("background genes: Spearman with lab-effect SD: Q_EB %.2f, Q_EB_additive %.2f, Q_EVT %.2f, keff_O2 %.2f\n",
              sp("Q_EB"), sp("Q_EB_additive"), sp("Q_EVT"), sp("keff_O2")))
}
dsum <- stats::aggregate(cbind(pair_mass_t9, in_complex_t9, A_pair_t1) ~ op + size, data = diagn, FUN = mean)
utils::write.csv(dsum, file.path(OUT, "mechanism_by_size.csv"), row.names = FALSE)
print(wsum, row.names = FALSE, digits = 3)

## ---- optional: does the background-defined width track biology or power? ------------
if (NSCALE) {
  sim <- simulate(20261030L + SEEDS[1]); roleg <- ifelse(sim$role %in% c("complex", "ribosome"),
    ifelse(sim$csize >= 20, "complex_>=20", "complex_<20"), sim$role)
  ns <- list()
  for (n_sub in c(125L, 250L, 500L)) {
    set.seed(7L); cols <- sort(sample(N, n_sub)); Ms <- mi_of(sim$X[, cols]); eg <- eb_graphon(Ms); qe <- rowSums(evt_q(Ms))
    Cs <- clr_of(Ms); set.seed(8L); C0 <- clr_of(mi_of(t(apply(sim$X[, cols], 1, sample)))); o2 <- o2_keff(Cs, C0)
    ns[[length(ns) + 1L]] <- do.call(rbind, lapply(split(seq_len(G), roleg), function(g)
      data.frame(N = n_sub, role = roleg[g[1]], Q_EB_median = stats::median(rowSums(eg$Q)[g]),
                 Q_EVT_median = stats::median(qe[g]),
                 keff_O2_median = stats::median(o2$keff[g]), pi0_EB = eg$pi0)))
    rm(Ms, Cs, C0, eg); invisible(gc(verbose = FALSE))
  }
  ns <- do.call(rbind, ns); utils::write.csv(ns, file.path(OUT, "width_vs_N.csv"), row.names = FALSE)
  print(ns, row.names = FALSE, digits = 3)
}
say("done")
