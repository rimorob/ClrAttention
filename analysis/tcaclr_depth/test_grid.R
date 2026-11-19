#!/usr/bin/env Rscript
# Smoke tests of grid_lib.R on a gene slice (default: first 800 genes of the M3D baselines matrix). Repo root:
#   OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 nice -n 19 Rscript analysis/tcaclr_depth/test_grid.R [--genes 800] [--matrix results/baselines_m3d/expression.tsv]
suppressPackageStartupMessages({ source("analysis/tcaclr_confirm/tc_lib.R"); source("analysis/tcaclr_depth/grid_lib.R") })
opt <- gl_options(commandArgs(TRUE), list(genes = 800L, matrix = "results/baselines_m3d/expression.tsv", threads = 8L))
say <- gl_say(); ok_all <- TRUE
chk <- function(what, ok, val = "") { say(sprintf("%-4s %s %s", if (ok) "PASS" else "FAIL", what, val)); if (!ok) ok_all <<- FALSE }
X <- gl_read_matrix(opt$matrix); X <- X[seq_len(opt$genes), ]; X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; G <- nrow(X)
say(sprintf("slice: %d genes x %d profiles", G, ncol(X)))
tm <- function(label, expr) { t0 <- Sys.time(); r <- force(expr); say(sprintf("time %-34s %.2f s", label, as.numeric(difftime(Sys.time(), t0, units = "secs")))); r }
f <- tm("fit (MI + calibrate)", gl_fit(X, opt$threads)); U <- gl_ut(G)

## 1. ut / utT = upper.tri order
ut0 <- which(upper.tri(matrix(0, G, G))); chk("gl_ut equals which(upper.tri)", identical(U$ut, ut0))
M <- matrix(seq_len(G * G), G, G); chk("gl_ut transpose index", identical(M[U$utT], t(M)[U$ut]))

## 2. clipped operator identical to ClrAttention$build_operator (registered path)
f$build_operator(alpha = GL$alpha, softmax_keff = 10, softmax_cap = 50L); A_reg <- f$operator; dimnames(A_reg) <- NULL
A_me <- tm("operator (clipped, keff 10 cap 50)", gl_operator_A(f$clr_scores, 10, 50L)); chk("operator A bit-identical (keff 10, cap 50)", identical(A_me, A_reg))
f$build_operator(alpha = GL$alpha, softmax_keff = 28, softmax_cap = 150L); A_reg <- f$operator; dimnames(A_reg) <- NULL
chk("operator A bit-identical (keff 28, cap 150)", identical(gl_operator_A(f$clr_scores, 28, 150L), A_reg))
f$build_operator(alpha = GL$alpha, softmax_keff = 10, softmax_cap = 50L)
P <- Matrix::Diagonal(G, 1 - GL$alpha) + GL$alpha * Matrix::Matrix(f$operator, sparse = TRUE)

## 3. registered chain and read-outs at t = 9 vs the grid chain
Pt_reg <- tm("registered chain t = 9 (8 sparse products)", gl_power_seq(P, 9L))
a_reg <- tc_pack(tc_clean_sym((Pt_reg + t(Pt_reg)) / 2)); tca_reg <- tc_pack(tc_clean_sym(sk_a3r(Pt_reg)))
v <- gl_variants(Pt_reg, U)
chk("aCLR read-out bit-identical on the registered chain", identical(v$aCLR, a_reg)); chk("tcaCLR (clipO) read-out bit-identical on the registered chain", identical(v$tca_clip, tca_reg))
got <- list(); tm("grid walk t = 9", gl_walk(P, 9L, function(t, Pt) got[[1]] <<- gl_variants(Pt, U)))
d_a <- max(abs(got[[1]]$aCLR - a_reg)); d_t <- max(abs(got[[1]]$tca_clip - tca_reg))
chk("grid chain (squarings) vs registered chain: aCLR max abs diff <= 1e-12", d_a <= 1e-12, format(d_a)); chk("grid chain vs registered chain: tcaCLR max abs diff <= 1e-12", d_t <= 1e-12, format(d_t))
say("dense products for the 16 grid depths: ", gl_n_dense(GL$depths), "; for t = 9: ", gl_n_dense(9L))

## 4. all 16 depths: the walk against the sequential chain at a few depths
sel <- c(1L, 3L, 6L, 12L, 24L); res <- list(); Pseq <- NULL
gl_walk(P, GL$depths, function(t, Pt) if (t %in% sel) res[[as.character(t)]] <<- gl_variants(Pt, U)$tca_clip)
for (t in sel) { Pq <- gl_power_seq(P, t); d <- max(abs(res[[as.character(t)]] - tc_pack(tc_clean_sym(sk_a3r(Pq)))))
  chk(sprintf("tcaCLR at t = %d: grid vs sequential chain, max abs diff <= 1e-9", t), d <= 1e-9, format(d)) }

## 5. AP from positives only equals ap_fast bit for bit (with heavy ties)
set.seed(1); v <- round(rnorm(50000), 1); v[sample(50000, 5000)] <- 0; Tx <- sample(50000, 700); lab <- integer(50000); lab[Tx] <- 1L
chk("gl_ap_T == ap_fast on tied scores", identical(gl_ap_T(v, gl_sorted(v), Tx), ap_fast(v, lab)))
v2 <- v; v2[1:5] <- NA; v2[Tx][1:3] <- NA   # NA are ranked last by ap_fast (-Inf); not used in the grid, skipped
vv <- got[[1]]$tca_clip; Tx <- gl_top(vv, gl_sorted(vv)); lab <- integer(length(vv)); lab[Tx] <- 1L
chk("gl_ap_T == ap_fast on a real tcaCLR vector", identical(gl_ap_T(vv, gl_sorted(vv), Tx), ap_fast(vv, lab)), format(ap_fast(vv, lab)))

## 6. unclipped vs clipped pair score
Sc <- f$clr_scores; Su <- gl_unclipped_clr(f$mi)
Z <- { M0 <- f$mi; diag(M0) <- 0; mu <- rowSums(M0) / G; M0 <- M0 - mu; M0 / sqrt(rowSums(M0^2) / (G - 1)) }
neg <- (Z < 0) | t(Z < 0); diag(neg) <- FALSE; same <- abs(Sc - Su) <= 1e-12; diag(same) <- TRUE
chk("unclipped == clipped wherever neither z is negative (tol 1e-12)", all(same[!neg]), sprintf("max |diff| %.2g", max(abs(Sc - Su)[!neg])))
chk("unclipped != clipped only where a z is negative", all(neg[!same]))
ut <- upper.tri(Sc); say(sprintf("pairs with a negative z: %.1f%%; pairs where the two scores differ: %.1f%%", 100 * mean(neg[ut]), 100 * mean(!same[ut])))
rowtop <- function(S, cap = 150L) t(vapply(seq_len(G), function(i) { r <- S[i, ]; r[i] <- -Inf; o <- order(r, decreasing = TRUE)[seq_len(min(cap, G - 1L))]; o[r[o] > 0][seq_len(min(cap, G - 1L))] }, integer(min(cap, G - 1L))))
Ac <- gl_operator_A(Sc, 10, 150L); Au <- gl_operator_A(Su, 10, 150L)
say(sprintf("operator keff 10 cap 150: rows with a different attended set %.1f%%; rows with different weights %.1f%%; mean |A_u - A_c| row L1 %.4f",
  100 * mean(vapply(seq_len(G), function(i) !identical(which(Ac[i, ] > 0), which(Au[i, ] > 0)), NA)), 100 * mean(rowSums(abs(Ac - Au)) > 1e-12), mean(rowSums(abs(Ac - Au)))))
cand <- vapply(seq_len(G), function(i) { r <- Sc[i, ]; r[i] <- -Inf; o <- order(r, decreasing = TRUE)[seq_len(150L)]; o <- o[r[o] > 0]; mean(Z[i, o] < 0 | Z[o, i] < 0) }, 0)
say(sprintf("share of the clipped path's top-150 candidates with a negative z in either direction: mean %.2f%%", 100 * mean(cand)))

## 7. stage 1 split function on the slice; neff identity; criteria sanity
ne <- gl_neff(X); say(sprintf("N = %d, N_eff = %.1f, G = %d, G_eff = %.1f", ne[["N"]], ne[["N_eff"]], ne[["G"]], ne[["G_eff"]]))
Zr <- t(scale(t(X))); Cg <- tcrossprod(Zr) / (ncol(X) - 1); chk("G_eff equals (tr R)^2 / tr(R^2) of the explicit gene correlation matrix", abs(ne[["G_eff"]] - sum(diag(Cg))^2 / sum(Cg^2)) < 1e-6 * ne[["G_eff"]])
Cc <- cor(X - rowMeans(X)); chk("N_eff equals (tr C)^2 / tr(C^2) of the explicit experiment correlation matrix", abs(ne[["N_eff"]] - sum(diag(Cc))^2 / sum(Cc^2)) < 1e-8 * ne[["N_eff"]])
chk("grid axes", length(GL$depths) == 17 && 9L %in% GL$depths && identical(gl_cap(GL$keffs), c(28L, 40L, 56L, 80L, 112L, 160L, 224L, 320L)))
chk("registered depths closed under the chain (dense products for the grid)", gl_n_dense(GL$depths) == 12, gl_n_dense(GL$depths))
set.seed(2); h <- sample(c(1L, 2L), ncol(X), TRUE); X1 <- X[, h == 1]; X2 <- X[, h == 2]
f1 <- gl_fit(X1, opt$threads); f2 <- gl_fit(X2, opt$threads); nes <- c(gl_neff(X1)[["N_eff"]], gl_neff(X2)[["N_eff"]])
ops <- gl_grid_ops(reference = TRUE, keffs = c(10, 40), depths = c(1L, 4L, 9L, 18L)); ops[[length(ops)]]$depths <- 9L
r <- tm("stage 1 split, 5 operators x 4 depths", gl_stage1_split(f1, f2, ops, nes, say))
chk("stage 1 rows: finite C1, C2, C2_null, C3, C3_bic", all(is.finite(as.matrix(r[, c("C1", "C2", "C2_null", "C3", "C3_bic")]))))
chk("C3 <= log-prob bound and C3_e10 finite", all(r$C3 < 0) && all(is.finite(r$C3_e10)))
chk("tcaCLR clipO TRUE/FALSE share C3 and k_reach; aCLR rows have clipO NA", isTRUE(all.equal(r$C3[r$method == "tcaCLR" & r$clipO %in% TRUE], r$C3[r$method == "tcaCLR" & r$clipO %in% FALSE])) && all(is.na(r$clipO[r$method == "aCLR"])))
print(r[r$direction == 1 & r$method == "tcaCLR" & r$clipO %in% TRUE, c("keff","cap","clipS","t","C1","C2","C2_null","C2_adj","C3","C3_e10","C3_bic","k_reach","c3_unreach")], row.names = FALSE, digits = 4)
say(if (ok_all) "ALL PASSED" else "SOME FAILED"); say(sprintf("peak RSS %.2f GB", gl_vmhwm_gb()))
