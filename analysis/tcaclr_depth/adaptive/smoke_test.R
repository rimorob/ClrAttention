#!/usr/bin/env Rscript
# Smoke tests of the adaptive scripts (gene slice, seconds):  OMP_NUM_THREADS=8 nice -n 19 Rscript analysis/tcaclr_depth/adaptive/smoke_test.R [--genes 500]
suppressPackageStartupMessages({ source("analysis/tcaclr_depth/adaptive/adaptive_lib.R"); source("analysis/tcaclr_depth/adaptive/planted_sim.R") })
opt <- gl_options(commandArgs(TRUE), list(genes = 500L, matrix = "results/baselines_m3d/expression.tsv", threads = 8L)); say <- gl_say(); ok_all <- TRUE
chk <- function(what, ok, val = "") { say(sprintf("%-4s %s %s", if (ok) "PASS" else "FAIL", what, val)); if (!ok) ok_all <<- FALSE }
X <- ad_load(opt$matrix, opt$genes); G <- nrow(X); say(sprintf("slice %d genes x %d experiments", G, ncol(X)))
f <- ad_fit(X, opt$threads); P <- ad_P_from_fit(f, 10, 40L)
## 1. operator is the grid library's (not a second builder) and its symmetrised form is a valid symmetric lazy operator
chk("ad_P_from_fit == gl_make_P", identical(P, gl_make_P(f$clr_scores, 10, 40L)))
Ps <- ad_Psym(P); chk("P_sym symmetric", Matrix::isSymmetric(Ps)); D <- as.matrix(Ps); ev <- eigen(D, symmetric = TRUE)
chk("dense eigenvalues in [0, 1] (lazy => nonnegative)", min(ev$values) > -1e-10 && max(ev$values) < 1 + 1e-10, sprintf("min %.2g max %.6f", min(ev$values), max(ev$values)))
e <- ad_eigs(Ps, 25L, vectors = TRUE); chk("RSpectra top 25 == dense eigen", max(abs(e$values - ev$values[1:25])) < 1e-8, format(max(abs(e$values - ev$values[1:25]))))
A <- (P - Matrix::Diagonal(G, 0.5)) / 0.5; W <- (A + Matrix::t(A)) / 2; d <- Matrix::rowSums(W); v1 <- e$vectors[, 1]; chk("trivial eigenvector proportional to sqrt(degree)", abs(abs(sum(v1 * sqrt(d))) / sqrt(sum(d))) > 1 - 1e-8)
chk("non-symmetric second modulus finite", is.finite(ad_edge_nonsym(P)), sprintf("%.4f vs symmetrised %.4f", ad_edge_nonsym(P), e$values[2]))
## 2. planting: marginals exact, pair correlation = dose (linear), no Pearson correlation (nonlin), AP from positives == ap_fast
X0 <- X[1:300, ]; pl <- ps_plant(X0, 10, 3, 0.5, "linear", 11L); m <- pl$modules[[1]]
chk("planted genes keep exactly their marginals", all(vapply(m, function(g) identical(unname(sort(pl$X[g, ])), unname(sort(ad_permute_genes(X0, 11L)[g, ]))), NA)))
R <- stats::cor(t(ad_rank_normal(pl$X[m, ]))); chk("linear plant: mean within-module correlation near dose 0.5", abs(mean(R[upper.tri(R)]) - 0.5) < 0.08, sprintf("%.3f", mean(R[upper.tri(R)])))
pn <- ps_plant(X0, 10, 3, 0.5, "nonlin", 11L); Zn <- ad_rank_normal(pn$X[pn$modules[[1]], ]); Rn <- stats::cor(t(Zn)); Rq <- stats::cor(t(Zn), t(Zn^2)); odd <- seq(1, 10, 2); ev <- seq(2, 10, 2)
chk("nonlin plant: linear-vs-U-shaped pairs have Pearson near 0 but dependence on the squared score", abs(mean(Rn[odd, ev])) < 0.1 && mean(Rq[ev, odd]) > 0.1, sprintf("%.3f / %.3f", mean(Rn[odd, ev]), mean(Rq[ev, odd])))
U <- gl_ut(300L); lab <- logical(length(U$ut)); lab[pl$truth] <- TRUE; chk("truth positions are within-module pairs", all(paste(U$ri[pl$truth], U$cj[pl$truth]) %in% unlist(lapply(pl$modules, function(g) { g <- sort(g); ij <- utils::combn(g, 2); paste(ij[1, ], ij[2, ]) }))))
set.seed(3); v <- round(rnorm(length(lab)), 1); chk("gl_ap_T == ap_fast on the planted truth (tied scores)", identical(gl_ap_T(v, gl_sorted(v), pl$truth), ap_fast(v, as.integer(lab))))
## 3. PMI read-out finite, and aCLR/tcaCLR variants of one P^t
fs <- ad_fit(pl$X, opt$threads); Pp <- ad_P_from_fit(fs, 10, 40L); pi <- gl_stationary(Pp); Pt <- as.matrix(Pp %*% Pp); Z <- ps_pmi_S(Pt, pi)
chk("diffPMI Z finite and non-negative", all(is.finite(Z)) && min(Z) >= 0, sprintf("share positive %.3f", mean(Z > 0)))
## 4. permuted control through the edge rule (3 permutations of a permuted slice) and no structure in the likelihood
se <- se_run(ad_permute_genes(X, 5L), c(10), 3L, 30L, 30L, opt$threads, 1L, 5L, "ctrl")$summary
chk("permuted control: edge computed; counts above edge reported", is.finite(se$edge_mean), sprintf("above mean-edge %d, above max-edge %d, permutation p %.2f", se$n_above_mean, se$n_above_max, se$obs_lambda2_perm_p))
say(if (ok_all) "ALL PASSED" else "SOME FAILED"); say(sprintf("peak RSS %.2f GB", ad_mem()))
