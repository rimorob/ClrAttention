#!/usr/bin/env Rscript
# EXPLORATORY. Conditional stage (idea 1 passed): M3D equivalence gate vs robust_z m3d V0 caches (1e-9), reference scoring, 6 lab splits,
# half-dataset K10 operators and attention graphs for leiden_sweep.py. STOPS on failed equivalence.
args <- commandArgs(TRUE); threads <- if (length(args)) as.integer(args[1]) else 8L
source("analysis/explore_ideas12/ideas12_lib.R"); source("analysis/explore_ideas12/ideas12_m3d_lib.R")
out <- file.path(EI_OUT, "m3d"); gd <- file.path(out, "graphs"); dir.create(file.path(out, "cache"), recursive = TRUE, showWarnings = FALSE); dir.create(gd, showWarnings = FALSE)
say <- ei_say(file.path(out, "stage_m3d_prep.log"))
env <- ei_setup_m3d(say); G <- env$G; X <- env$X; lab <- env$lab
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads); f$calibrate(method = "normal", combine = "stouffer")
M <- f$mi; C <- f$clr_scores; f$build_operator(alpha = EI_ALPHA, softmax_keff = 10, softmax_cap = 50L); A <- f$operator
P <- Matrix::Diagonal(G, 1 - EI_ALPHA) + EI_ALPHA * Matrix::Matrix(A, sparse = TRUE); S9 <- ei_attention(A)
q_att <- rz_stationary(P); H_att <- outer(q_att, q_att, "+"); mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+")
saveRDS(list(q_att = q_att, mi_rm = mi_rm, genes = env$genes), file.path(out, "cache", "m3d_core.rds")); f$release(); rm(f, M, P); invisible(gc(verbose = FALSE))
ev <- ei_make_eval(env, H_att, H_mi, say); rm(H_att, H_mi)
r_clr <- ev(C); r_k10 <- ev(S9); saveRDS(r_clr, file.path(out, "cache", "ref_CLR.rds")); saveRDS(r_k10, file.path(out, "cache", "ref_K10.rds"))
rz <- "results/robust_z/m3d/cache"; rows <- list()
for (nm in c("CLR", "K10")) { ref <- readRDS(file.path(rz, c(CLR = "V0__pair.rds", K10 = "V0__att9.rds")[[nm]])); mine <- if (nm == "CLR") r_clr else r_k10
  for (sc in c("pooled", "hub_att", "hub_mi")) for (e in names(ref[[sc]])) rows[[length(rows) + 1L]] <- data.frame(readout = nm, scope = sc, truth = e,
    max_abs_diff_aupr = abs(mine[[sc]][[e]]$aupr - ref[[sc]][[e]]$aupr), max_abs_diff_jack = max(abs(mine[[sc]][[e]]$jack - ref[[sc]][[e]]$jack))) }
eq <- do.call(rbind, rows); eq$pass <- eq$max_abs_diff_aupr < 1e-9 & eq$max_abs_diff_jack < 1e-9
utils::write.csv(eq, file.path(out, "equivalence.csv"), row.names = FALSE)
say(sprintf("equivalence: %d checks, max |dAUPR| %.3g, max |djack| %.3g", nrow(eq), max(eq$max_abs_diff_aupr), max(eq$max_abs_diff_jack)))
if (!all(eq$pass)) { print(eq[!eq$pass, ]); say("EQUIVALENCE FAILED: STOP"); quit(status = 2) }
say("EQUIVALENCE PASSED (1e-9)")
write_graph <- function(A, name) { W <- (A + t(A)) / 2; diag(W) <- 0; W[lower.tri(W, diag = TRUE)] <- 0; ix <- which(W > 0, arr.ind = TRUE)
  utils::write.table(data.frame(i = ix[, 1] - 1L, j = ix[, 2] - 1L, w = format(W[ix], digits = 17)), file.path(gd, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); nrow(ix) }
say(sprintf("full graph: %d edges", write_graph(as.matrix(A), "full"))); rm(A, C, S9)
lab_split <- function() { sz <- table(lab); sz <- sz[sample(length(sz))]; half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }; unname(half[lab]) }
set.seed(EI_DEPTH_SEED); splits <- lapply(1:6, function(s) lab_split())
for (si in 1:6) { h <- splits[[si]]; set.seed(EI_DEPTH_SEED + si)
  for (d in 1:2) { Xm <- X[, h == d]; rng <- apply(Xm, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9)
    Mh <- bspline_mi(Xm, bins = "hg", spline_order = 3, threads = threads, transform = "none"); Ah <- rz_operator(rz_scores(Mh, "V0")); rm(Mh)
    say(sprintf("split %d half %d: %d experiments, %d edges", si, d, ncol(Xm), write_graph(Ah, sprintf("s%dh%d", si, d)))); rm(Ah) } }
say("done")
