#!/usr/bin/env Rscript
# Reference check: the grid code at the registered cell (keff 10, cap 50, t 9, clipped both) must reproduce the registered aCLR and tcaCLR scores.
#   Rscript analysis/tcaclr_depth/check_reference.R [--expr results/tcaclr_confirm/armB/expression.tsv] [--scores results/tcaclr_confirm/armB/scores] [--genes 0] [--threads 8]
# Full mode (--genes 0): refit MI -> CLR on the expression file, rebuild the operator, and compare, for aCLR and tcaCLR,
#   (a) the grid path (gl_make_P + gl_walk power chain by squarings + gl_variants) against the STORED packed score files,
#   (b) a registered-style recomputation (make_scores.R fit logic: build_operator, 8 sequential sparse products, sk_a3r) against the stored files,
#   (c) the grid path against the registered-style recomputation, and the operator matrix A (identical?).
# Slice mode (--genes N, no stored files): the first N genes of the expression file; only (b)-(c) and the operator identity are checked.
# Passes when every maximum absolute difference is <= 1e-12 (the squaring chain differs from the sequential chain only by floating-point rounding, about 1e-15).
suppressPackageStartupMessages({ source("analysis/tcaclr_confirm/tc_lib.R"); source("analysis/tcaclr_depth/grid_lib.R") })
opt <- gl_options(commandArgs(TRUE), list(expr = "results/tcaclr_confirm/armB/expression.tsv", scores = "results/tcaclr_confirm/armB/scores", genes = 0L, threads = 8L))
say <- gl_say(); ok <- TRUE; chk <- function(what, pass, val = "") { say(sprintf("%-4s %s %s", if (pass) "PASS" else "FAIL", what, val)); if (!pass) ok <<- FALSE }
X <- gl_read_matrix(opt$expr); if (opt$genes > 0L) X <- X[seq_len(opt$genes), ]; genes <- rownames(X); G <- nrow(X); full <- opt$genes == 0L
say(sprintf("%s: %d genes x %d profiles; mode %s", opt$expr, G, ncol(X), if (full) "full (stored scores)" else "slice"))
if (full) { st <- tc_scores_open(opt$scores); chk("gene order identical to genes.txt of the score directory", identical(st$genes, genes)); stored <- list(aCLR = tc_scores_packed(st, "aCLR"), tcaCLR = tc_scores_packed(st, "tcaCLR")) }
t0 <- Sys.time(); f <- gl_fit(X, opt$threads); say(sprintf("fit %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
## registered path (make_scores.R fit logic)
f$build_operator(alpha = TC$alpha, softmax_keff = TC$keff, softmax_cap = as.integer(TC$cap)); A_reg <- f$operator; dimnames(A_reg) <- NULL
P <- Matrix::Diagonal(G, 1 - TC$alpha) + TC$alpha * Matrix::Matrix(f$operator, sparse = TRUE)
Pt <- gl_power_seq(P, TC$t); reg <- list(aCLR = tc_pack(tc_clean_sym((Pt + t(Pt)) / 2)), tcaCLR = tc_pack(tc_clean_sym(sk_a3r(Pt)))); rm(Pt); invisible(gc(verbose = FALSE))
## grid path
S <- f$clr_scores; chk("clipped pair score is the registered CLR score", identical(S, f$clr_scores))
A_me <- gl_operator_A(S, GL$ref$keff, GL$ref$cap); chk("operator A bit-identical to ClrAttention$build_operator (keff 10, cap 50)", identical(A_me, A_reg)); rm(A_me, A_reg)
Pg <- gl_make_P(S, GL$ref$keff, GL$ref$cap); chk("lazy operator P identical to make_scores.R", max(abs(as.matrix(Pg - P))) == 0)
U <- gl_ut(G); me <- NULL; t1 <- Sys.time(); gl_walk(Pg, GL$ref$t, function(t, Pt) me <<- gl_variants(Pt, U)); say(sprintf("grid chain t = 9: %.1f s", as.numeric(difftime(Sys.time(), t1, units = "secs"))))
res <- list(); md <- function(a, b) max(abs(a - b))
for (m in c("aCLR", "tcaCLR")) { g <- if (m == "aCLR") me$aCLR else me$tca_clip
  if (full) { res[[length(res) + 1L]] <- data.frame(method = m, comparison = "grid path vs STORED", max_abs_diff = md(g, stored[[m]]), n = length(g))
              res[[length(res) + 1L]] <- data.frame(method = m, comparison = "registered-style recomputation vs STORED", max_abs_diff = md(reg[[m]], stored[[m]]), n = length(g)) }
  res[[length(res) + 1L]] <- data.frame(method = m, comparison = "grid path vs registered-style recomputation", max_abs_diff = md(g, reg[[m]]), n = length(g)) }
res <- do.call(rbind, res); print(res, row.names = FALSE, digits = 4)
for (k in seq_len(nrow(res))) chk(sprintf("%s, %s: max abs diff <= 1e-12", res$method[k], res$comparison[k]), res$max_abs_diff[k] <= 1e-12, format(res$max_abs_diff[k]))
say(if (ok) "REFERENCE CHECK PASSED" else "REFERENCE CHECK FAILED"); say(sprintf("peak RSS %.2f GB", gl_vmhwm_gb())); if (!ok) quit(status = 1)
