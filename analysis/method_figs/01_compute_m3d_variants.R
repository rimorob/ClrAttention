# Second worked-example script: the same nine-gene M3D neighbourhood (from 00_compute_m3d.R) read through the other rungs of the ladder (tcaCLR2 depth, global temperature),
# and the profiles of one pair for the |r| panel. Output: analysis/method_figs/cache/m3d_variants.rds. Run after 00_compute_m3d.R:
#   nice -n 19 Rscript analysis/method_figs/01_compute_m3d_variants.R
suppressPackageStartupMessages(source("analysis/global_temp_ladder/gt_lib.R"))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds")
X <- gl_truth_matrix("m3d_tf_regulons"); G <- nrow(X); genes <- rownames(X); Nb <- ex$Nb
f <- ad_fit(X, 8L); S <- f$clr_scores; Tm <- gt_topmat(S, 50L)
sec <- function(Pt) { Z <- gl_zmat(Pt); S2 <- (pmax(Z, 0) + pmax(t(Z), 0)) / sqrt(2); diag(S2) <- 0; S2 }
P <- gl_make_P(S, 10, 50)
S2_t2 <- sec(gl_power_seq(P, ex$t2)); say("tcaCLR2 (t2 = ", ex$t2, ") done")
S2_ref9 <- sec(gl_power_seq(P, 9L)); stopifnot(max(abs(S2_ref9[Nb, Nb] - ex$S2_sub)) < 1e-9)
Pv <- as.matrix(gt_make_P(Tm, ex$tau_v2)); S2_v2 <- sec(gl_power_seq(Pv, 9L)); say("V2 done")
gtab <- utils::read.delim("analysis/ecoli_k12_genes.tsv", stringsAsFactors = FALSE); sy <- gtab$symbol[match(genes, gtab$bnumber)]
stopifnot(identical(sy[Nb], ex$names))
a <- Nb[1]; b <- Nb[4]   # lexA, recA
r_ab <- cor(X[a, ], X[b, ]); say(sprintf("Pearson r(%s, %s) = %.4f over %d profiles", sy[a], sy[b], r_ab, ncol(X)))
out <- list(S2_t2_sub = S2_t2[Nb, Nb], S2_v2_sub = S2_v2[Nb, Nb], x = X[a, ], y = X[b, ], r = r_ab, pair = c(sy[a], sy[b]), n = ncol(X),
            absr_sub = abs(cor(t(X[Nb, ]))), t2 = ex$t2)
saveRDS(out, "analysis/method_figs/cache/m3d_variants.rds"); say("saved")
