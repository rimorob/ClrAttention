#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Evaluates the PI's asymmetry proposal (columns = regulators /
# latent factors, rows = observed genes, as in the gene x TF form of CLR, Faith
# et al. 2007) against CLR, attention and the a3r readout, on simulated data
# with planted latents and known truth. No real data or truth set is read.
#
# Generative model (per seed): G genes x N samples, Gaussian factor model.
#   - cliques ("complexes"): sizes 4, 6, 10, 20, 40 at loadings 0.60, 0.75, 0.90
#     (pairwise r = 0.36, 0.56, 0.81), plus one ribosome-like clique of 80 at 0.85;
#   - satellites: 25 genes loading 0.6 on the ribosome factor only (not members;
#     they mimic INTS8/EIF4B, which spend most of their budget on the 60S);
#   - regulons: 10 regulator activities; regulator mRNA loads 0.5-0.8 on its own
#     activity; 10-150 targets each drawn from a pool of 380 (overlapping), target
#     loadings U(0.25, 0.60) (loose);
#   - hubs: 20 genes loading 0.30 on 6 random factors each;
#   - global factor 0.20 on a random 30% of genes; the rest is background noise.
# Every gene has unit variance. Dependence = Gaussian MI, -0.5 log(1 - r^2)
# (monotone in |r|; the method's B-spline MI is not needed for Gaussian data).
#
# Operators (all from the same CLR matrix and the same softmax, keff 10, cap 50,
# lazy alpha 0.5, as the method):
#   absr, CLR, att (P^t + P^t')/2, a3r (bilateral row-CLR on log P^t, reachable);
#   C1 a3c: Faith-2007 asymmetric CLR on log P^t: each entry i->j is z-scored in
#        its row (gene i's background over keys) and in its column (key j's
#        background over all querying genes), Stouffer; symmetric = (D + D')/2;
#   C2 sinkhorn: doubly stochastic attention (Sinkhorn scaling of the lazy
#        operator: every key receives the same total attention), att and a3r readouts;
#   C3 lcclr: latent-column CLR: each column j defines a latent f_j = attention-
#        weighted mean of the genes that attend to j (column-normalized lazy P,
#        leave-i-out), D_ij = MI(x_i, f_j), asymmetric CLR (row z over latents,
#        column z over genes), Stouffer; symmetric = (D + D')/2;
#   C4 residual: tight units (average linkage on 1 - r, r >= 0.6), unit factors
#        regressed out of every gene, MI/CLR/attention on residuals (C4a); C4b puts
#        within-unit pairs first (by original CLR), then C4a.
# Scenarios: base (as above); global (adds 2 strong global programs with gene-
# specific loadings N(0, 0.25^2), the structure CLR's per-gene background removes);
# tf_matched (regulator mRNA loading drawn like the targets', U(0.25, 0.60): the
# control for the direction metric).
# Usage: Rscript sim_asymmetry.R [seeds=1:4] [scenario=base] [outdir]
suppressPackageStartupMessages(library(clr))
args <- commandArgs(TRUE)
SEEDS <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:4
SCEN <- if (length(args) >= 2) args[2] else "base"          # base | global | tf_matched
OUT <- if (length(args) >= 3) args[3] else file.path("analysis/synthetic_asymmetry/out", SCEN)
stopifnot(SCEN %in% c("base", "global", "tf_matched"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
G <- 1000L; N <- 400L; TS <- c(4L, 9L)
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

## ---- simulation ---------------------------------------------------------------
simulate <- function(seed) {
  set.seed(seed)
  role <- rep("background", G); grp <- rep(NA_character_, G)
  facs <- list(); load <- list()          # factor name -> N-vector; gene -> named loadings
  add_load <- function(g, f, b) for (k in seq_along(g)) load[[g[k]]] <<- c(load[[g[k]]], stats::setNames(b[k], f))
  load <- vector("list", G)
  nxt <- 1L; take <- function(n) { idx <- nxt:(nxt + n - 1L); nxt <<- nxt + n; idx }
  cliques <- list(); cinfo <- data.frame()
  for (lam in c(0.60, 0.75, 0.90)) for (n in c(4, 6, 10, 20, 40)) {
    nm <- sprintf("C%d_%.2f", n, lam); g <- take(n); facs[[nm]] <- stats::rnorm(N)
    add_load(g, nm, rep(lam, n)); role[g] <- "clique"; grp[g] <- nm; cliques[[nm]] <- g
    cinfo <- rbind(cinfo, data.frame(clique = nm, size = n, loading = lam)) }
  g <- take(80); facs$RIBO <- stats::rnorm(N); add_load(g, "RIBO", rep(0.85, 80))
  role[g] <- "ribosome"; grp[g] <- "RIBO"; cliques$RIBO <- g
  cinfo <- rbind(cinfo, data.frame(clique = "RIBO", size = 80, loading = 0.85))
  g <- take(25); add_load(g, "RIBO", rep(0.6, 25)); role[g] <- "satellite"
  tfs <- take(10); pool <- take(380); sizes <- c(10, 15, 20, 30, 40, 60, 80, 100, 120, 150)
  regulons <- list(); edges <- data.frame()
  for (r in 1:10) { nm <- sprintf("R%d", r); facs[[nm]] <- stats::rnorm(N)
    add_load(tfs[r], nm, if (SCEN == "tf_matched") stats::runif(1, 0.25, 0.60) else stats::runif(1, 0.5, 0.8))
    role[tfs[r]] <- "regulator"
    tg <- sort(sample(pool, sizes[r])); add_load(tg, nm, stats::runif(length(tg), 0.25, 0.60))
    regulons[[nm]] <- c(tfs[r], tg); edges <- rbind(edges, data.frame(tf = tfs[r], target = tg)) }
  role[pool] <- "target"
  hubs <- take(20); fn <- names(facs)
  hubfac <- list(); for (h in hubs) { f <- sample(fn, 6); hubfac[[as.character(h)]] <- f; for (ff in f) add_load(h, ff, 0.30) }
  role[hubs] <- "hub"
  facs$GLOBAL <- stats::rnorm(N); gl <- sample(G, round(0.3 * G)); add_load(gl, "GLOBAL", rep(0.20, length(gl)))
  if (SCEN == "global") for (k in 1:2) { nm <- paste0("PROG", k); facs[[nm]] <- stats::rnorm(N)
    add_load(seq_len(G), nm, stats::rnorm(G, 0, 0.25)) }   # strong global programs, gene-specific loadings
  Fm <- do.call(rbind, facs)
  X <- matrix(0, G, N)
  for (i in seq_len(G)) { b <- load[[i]]
    if (length(b)) { if (sum(b^2) > 0.9) b <- b * sqrt(0.9 / sum(b^2))
      X[i, ] <- colSums(b * Fm[names(b), , drop = FALSE]) + sqrt(1 - sum(b^2)) * stats::rnorm(N)
    } else X[i, ] <- stats::rnorm(N) }
  X <- (X - rowMeans(X)) / apply(X, 1, stats::sd)
  list(X = X, role = role, grp = grp, cliques = cliques, cinfo = cinfo, regulons = regulons,
       edges = edges, tfs = tfs, hubs = hubs, hubfac = hubfac)
}

## ---- operators -----------------------------------------------------------------
gmi <- function(R) { R <- pmin(abs(R), 1 - 1e-12); M <- -0.5 * log(1 - R^2); diag(M) <- 0; M }
clr_of <- function(X) clr_calibrate(gmi(stats::cor(t(X))), method = "normal", combine = "stouffer")
softmax_A <- function(S) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L)
  e <- rowSums(A) == 0; A[cbind(which(e), which(e))] <- 1; A }
lazy <- function(A) { P <- 0.5 * A; diag(P) <- diag(P) + 0.5; P }
# z-scores of non-missing entries along rows (margin 1) or columns (margin 2); negatives and NA -> 0
zna <- function(L, margin = 1) {
  if (margin == 2) return(t(zna(t(L), 1)))
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  s <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1))
  Z <- (L - mu) / s; Z[!is.finite(Z) | Z < 0] <- 0; Z }
logL <- function(Pt) { L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA; L }
a3r <- function(Pt) { Z <- zna(logL(Pt), 1); S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; list(S = S, D = Z) }
a3c <- function(Pt) { L <- logL(Pt); D <- (zna(L, 1) + zna(L, 2)) / sqrt(2); diag(D) <- 0
  list(S = (D + t(D)) / 2, D = D) }
att <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; list(S = S, D = Pt) }
sinkhorn <- function(K, tol = 1e-9, maxit = 5000) {
  v <- rep(1, ncol(K))
  for (it in seq_len(maxit)) { u <- 1 / drop(K %*% v); v <- 1 / drop(crossprod(K, u))
    if (it %% 10 == 0 && max(abs(u * drop(K %*% v) - 1)) < tol) break }
  K * outer(u, v) }
# latent-column CLR: D_ij = MI(x_i, f_j^(-i)), f_j = sum_k W_kj x_k with W column-normalized
lcclr <- function(X, P) {
  C <- stats::cor(t(X)); W <- sweep(P, 2, colSums(P), "/")
  CW <- C %*% W; V <- colSums(W * CW)
  num <- CW - W; vv <- sweep(-2 * W * CW + W^2, 2, V, "+")
  R <- num / sqrt(pmax(vv, 1e-12)); R[vv < 1e-9] <- 0
  M <- gmi(R); D <- (zna(M, 1) + zna(M, 2)) / sqrt(2); diag(D) <- 0
  list(S = (D + t(D)) / 2, D = D, R = R) }
walk <- function(P, ts) { out <- list(); Pt <- P
  for (s in seq_len(max(ts))) { if (s > 1) Pt <- Pt %*% P; if (s %in% ts) out[[as.character(s)]] <- Pt }
  out }

## ---- metrics ---------------------------------------------------------------------
ap_sorted <- function(s, l) { n <- length(s); if (!n || !sum(l)) return(NA_real_)
  last <- c(s[-1L] != s[-n], TRUE); tp <- cumsum(l)[last]; k <- which(last)
  sum(diff(c(0, tp)) * tp / k) / sum(l) }
auroc <- function(pos, neg) { r <- rank(c(pos, neg)); np <- length(pos); nn <- length(neg)
  (sum(r[seq_len(np)]) - np * (np + 1) / 2) / (np * nn) }
ap <- function(s, l) { s[!is.finite(s)] <- -Inf; o <- order(s, decreasing = TRUE, method = "radix"); ap_sorted(s[o], l[o]) }

make_truth <- function(sim) {
  ut <- which(upper.tri(matrix(0, G, G)), arr.ind = TRUE); i <- ut[, 1]; j <- ut[, 2]
  cl <- integer(G); for (k in seq_along(sim$cliques)) cl[sim$cliques[[k]]] <- k
  cpx <- cl[i] > 0 & cl[i] == cl[j]
  ribo <- sim$cliques$RIBO; inribo <- seq_len(G) %in% ribo
  RM <- matrix(FALSE, G, length(sim$regulons)); for (k in seq_along(sim$regulons)) RM[sim$regulons[[k]], k] <- TRUE
  reg <- rowSums(RM[i, , drop = FALSE] & RM[j, , drop = FALSE]) > 0
  E <- matrix(FALSE, G, G); E[cbind(sim$edges$tf, sim$edges$target)] <- TRUE; E <- E | t(E)
  istf <- seq_len(G) %in% sim$tfs
  mods <- c(sim$cliques, sim$regulons); HM <- matrix(FALSE, G, G); WM <- matrix(FALSE, G, G)
  for (h in sim$hubs) for (f in sim$hubfac[[as.character(h)]]) { g <- mods[[f]]; HM[h, g] <- TRUE; WM[g, g] <- TRUE }
  diag(WM) <- FALSE
  list(i = i, j = j, cpx = cpx, noribo = !(inribo[i] | inribo[j]), reg = reg, edge = E[cbind(i, j)],
       tfpair = xor(istf[i], istf[j]), cl = cl,
       sat_ribo = (sim$role[i] == "satellite" & inribo[j]) | (sim$role[j] == "satellite" & inribo[i]),
       hubpair = sim$role[i] == "hub" | sim$role[j] == "hub",
       ribo_pair = inribo[i] & inribo[j],
       # hub controls: (hub, member of a module the hub loads on) vs within-pairs of those modules
       hub_mem = HM[cbind(i, j)] | HM[cbind(j, i)], hubmod_pos = WM[cbind(i, j)])
}
score_all <- function(name, S, D, tr, sim) {
  s <- S[cbind(tr$i, tr$j)]; s[!is.finite(s)] <- -Inf
  o <- order(s, decreasing = TRUE, method = "radix")
  K <- sum(tr$cpx); topc <- o[seq_len(K)]; Kr <- sum(tr$reg); topr <- o[seq_len(Kr)]
  m <- tr$noribo
  pc <- sapply(seq_along(sim$cliques), function(k) {     # per-clique AP: pairs touching a member
    touch <- (tr$cl[tr$i] == k) | (tr$cl[tr$j] == k); lab <- tr$cpx & tr$cl[tr$i] == k
    ap(s[touch], lab[touch]) })
  dir <- NA_real_
  if (!is.null(D)) { e <- sim$edges
    dir <- mean(D[cbind(e$target, e$tf)] > D[cbind(e$tf, e$target)]) }
  row <- data.frame(op = name,
    cpx_x = ap(s, tr$cpx) / mean(tr$cpx),
    cpx_noribo_x = ap(s[m], tr$cpx[m]) / mean(tr$cpx[m]),
    reg_x = ap(s, tr$reg) / mean(tr$reg),
    edge_x = ap(s[tr$tfpair], tr$edge[tr$tfpair]) / mean(tr$edge[tr$tfpair]),
    sat_in_topK = mean(tr$sat_ribo[topc]), hub_in_topKreg = mean(tr$hubpair[topr]),
    sat_auc = auroc(s[tr$ribo_pair], s[tr$sat_ribo]), hub_auc = auroc(s[tr$hubmod_pos], s[tr$hub_mem]),
    dir_target_to_tf = dir)
  pcd <- cbind(sim$cinfo, op = name, ap_clique = pc)
  list(row = row, pc = pcd)
}

## ---- run -------------------------------------------------------------------------
res <- list(); pcs <- list(); diag_rows <- list()
for (seed in SEEDS) {
  t0 <- Sys.time(); sim <- simulate(20261026L + seed); X <- sim$X; tr <- make_truth(sim)
  ops <- list()
  C <- clr_of(X); R <- abs(stats::cor(t(X))); diag(R) <- 0
  Zc <- zna({ M <- gmi(R); diag(M) <- NA; M }, 1)              # CLR's directed half (row z of MI)
  ops$absr <- list(S = R, D = NULL); ops$CLR <- list(S = C, D = Zc)
  A <- softmax_A(C); P <- lazy(A); Pw <- walk(P, TS)
  Pds <- sinkhorn(P); Pdw <- walk(Pds, TS)
  for (t in TS) { k <- as.character(t)
    ops[[paste0("att_t", t)]] <- att(Pw[[k]]); ops[[paste0("a3r_t", t)]] <- a3r(Pw[[k]])
    ops[[paste0("C1_a3c_t", t)]] <- a3c(Pw[[k]])
    ops[[paste0("C2_sink_att_t", t)]] <- att(Pdw[[k]]); ops[[paste0("C2_sink_a3r_t", t)]] <- a3r(Pdw[[k]]) }
  lc <- lcclr(X, P); ops$C3_lcclr <- lc[c("S", "D")]
  lc2 <- lcclr(X, Pw[["4"]]); ops$C3_lcclr_t4 <- lc2[c("S", "D")]
  # C4: residual attention after removing tight-unit factors
  hc <- stats::hclust(stats::as.dist(1 - stats::cor(t(X))), "average"); u <- stats::cutree(hc, h = 0.4)
  big <- as.integer(names(which(table(u) >= 3)))
  Fu <- t(sapply(big, function(k) colMeans(X[u == k, , drop = FALSE])))
  Xr <- X - (X %*% t(Fu)) %*% solve(Fu %*% t(Fu), Fu)
  Xr <- (Xr - rowMeans(Xr)) / pmax(apply(Xr, 1, stats::sd), 1e-12)
  Pr <- walk(lazy(softmax_A(clr_of(Xr))), 9L)[["9"]]
  ops$C4a_resid_att_t9 <- att(Pr)
  inunit <- outer(u, u, "==") & outer(u %in% big, u %in% big, "&"); diag(inunit) <- FALSE
  S4 <- ops$C4a_resid_att_t9$S / max(ops$C4a_resid_att_t9$S); S4[inunit] <- 1e6 + C[inunit]
  ops$C4b_units_first <- list(S = S4, D = NULL)
  for (nm in names(ops)) { sc <- score_all(nm, ops[[nm]]$S, ops[[nm]]$D, tr, sim)
    res[[length(res) + 1]] <- cbind(seed = seed, sc$row); pcs[[length(pcs) + 1]] <- cbind(seed = seed, sc$pc) }
  # mechanism diagnostics: per-clique mean member-member weight, mass retained in clique
  for (k in seq_along(sim$cliques)) { g <- sim$cliques[[k]]; n <- length(g)
    off <- function(M) { m <- M[g, g]; diag(m) <- NA; mean(m, na.rm = TRUE) }
    diag_rows[[length(diag_rows) + 1]] <- data.frame(seed = seed, sim$cinfo[k, ],
      A_pair = off(A), A_inclique = mean(rowSums(A[g, g])),
      P9_pair = off(Pw[["9"]]), P9_inclique = mean(rowSums(Pw[["9"]][g, g]) - diag(Pw[["9"]][g, g])) /
        mean(1 - diag(Pw[["9"]][g, g])),
      Pds9_pair = off(Pdw[["9"]]), in_degree_A = mean(colSums(A)[g]),
      clr_pair = off(C), a3r9_pair = off(ops$a3r_t9$S), lc_R_pair = off(lc$R), r_pair = off(R),
      n_units_hit = length(unique(u[g])), units_ge3 = sum(u[g] %in% big) / n) }
  pi_ <- colSums(Pw[["9"]]) / G
  say(sprintf("seed %d done in %.0f s (units >= 3: %d covering %d genes; pi max/median %.1f; sinkhorn colsum err %.1e)",
              seed, as.numeric(difftime(Sys.time(), t0, units = "secs")), length(big), sum(u %in% big),
              max(pi_) / stats::median(pi_), max(abs(colSums(Pds) - 1))))
}
res <- do.call(rbind, res); pcs <- do.call(rbind, pcs); dg <- do.call(rbind, diag_rows)
utils::write.csv(res, file.path(OUT, "metrics_by_seed.csv"), row.names = FALSE)
num <- setdiff(names(res), c("seed", "op"))
agg <- do.call(rbind, lapply(split(res, res$op), function(d)
  data.frame(op = d$op[1], t(sapply(num, function(v) sprintf("%.3f ± %.3f", mean(d[[v]]), stats::sd(d[[v]])))))))
agg <- agg[match(unique(res$op), agg$op), ]
utils::write.csv(agg, file.path(OUT, "metrics_summary.csv"), row.names = FALSE)
pcs$size_class <- ifelse(pcs$clique == "RIBO", "80 (ribosome)", as.character(pcs$size))
pca <- stats::aggregate(ap_clique ~ op + loading + size_class, data = pcs, FUN = mean)
utils::write.csv(pca, file.path(OUT, "per_clique_ap.csv"), row.names = FALSE)
dga <- stats::aggregate(. ~ clique + size + loading, data = dg[, setdiff(names(dg), "seed")], FUN = mean)
utils::write.csv(dga[order(dga$loading, dga$size), ], file.path(OUT, "mechanism_by_clique.csv"), row.names = FALSE)
print(agg, row.names = FALSE)
say("done")
