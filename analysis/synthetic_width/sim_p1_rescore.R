#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Rescore of the round-2 runs (sim_round2.R, seeds 20261031 + s) with the endpoint
# construction of the DRAFT yeast preregistration (docs/preregistration_reciprocity_yeast.md, P1/S1),
# requested by the review (docs/plans/yeast_prereg_review.md, S1). No real data, no truth file.
#
# Per run (regenerated bit-for-bit from the round-2 generator; sanity: K10/RECIP small_x at t = 9 must
# equal the stored round-2 metrics):
#   P1  universe = pairs of genes that are both in some complex (cliques incl. the 80-gene ribosome);
#       positives = co-members of a 3-6-gene complex; pairs co-members only of a larger complex excluded;
#       negatives = the rest (pairs of complex genes in different complexes).
#   S1  universe = pairs of complex genes; positives = co-members of any complex.
#   REG regulon co-membership as round 2 (all pairs).
# Operators: CLR (pairwise), K10, R (= recip(K10)); attention readout (P^t + P^t')/2.
# Depths: t = 9 and each operator's own t* by criterion A on the lab halves (labs 1-5 vs 6-10):
#   train on one half, AP of S_t for the other half's top-1% CLR pairs, grid {1,...,32}; per direction
#   argmax; geometric mean of the two directions; snapped to the grid (one split, as in keff_test.R).
# Usage: Rscript sim_p1_rescore.R <runs: e.g. std:0.1:1:1,std:0.1:1:2 | all> <workers>
suppressPackageStartupMessages(library(clr))
args <- commandArgs(TRUE)
RUNS <- if (args[1] == "all") {
  g <- expand.grid(s = 1:4, b = 0:2, f = c(0.05, 0.10, 0.25, 0.50)); r <- sprintf("std:%s:%d:%d", g$f, g$b, g$s)
  c(r, sprintf("regheavy:0.1:%d:%d", rep(0:1, each = 4), rep(1:4, 2)))
} else strsplit(args[1], ",")[[1]]
WORKERS <- as.integer(args[2])
src <- readLines("analysis/synthetic_width/sim_width.R")
src <- src[seq_len(grep("^## ---- run", src) - 1L)]
src <- src[!grepl("^args <-|^SEEDS <-|^SCEN <-|^NSCALE <-|^OUT <-|^dir.create", src)]
SCEN <- "messy"; eval(parse(text = src))
r2 <- readLines("analysis/synthetic_width/sim_round2.R")
eval(parse(text = r2[grep("^simulate2 <- function", r2):(grep("^# Efron central matching", r2) - 1L)]))
OUT <- "analysis/synthetic_width/out/p1_rescore"; dir.create(file.path(OUT, "runs"), showWarnings = FALSE, recursive = TRUE)
TGRID <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32)
ap_abs <- function(s, lab) { s[!is.finite(s)] <- -Inf; o <- order(s, decreasing = TRUE, method = "radix"); ap_sorted(s[o], lab[o]) }

one_run <- function(run) {
  p <- strsplit(run, ":")[[1]]; GEN <<- p[1]; FPL <<- as.numeric(p[2]); BAT <<- as.numeric(p[3]); SEED <- as.integer(p[4])
  TAG <- sprintf("%s_f%02d_b%d_s%d", GEN, round(100 * FPL), BAT, SEED); fout <- file.path(OUT, "runs", paste0(TAG, ".csv"))
  if (file.exists(fout)) return(TAG)
  t0 <- Sys.time(); sim <- simulate2(20261031L + SEED); X <- sim$X; tr <- make_truth(sim)
  cl <- integer(G); for (k in seq_along(sim$cliques)) cl[sim$cliques[[k]]] <- k
  ii <- (tr$idx - 1L) %% G + 1L; jj <- (tr$idx - 1L) %/% G + 1L
  inc <- cl[ii] > 0 & cl[jj] > 0; same <- inc & cl[ii] == cl[jj]; szp <- ifelse(same, sim$csize[ii], 0L)
  p1_keep <- inc & (!same | szp %in% 3:6); p1_lab <- same & szp %in% 3:6
  s1_keep <- inc; s1_lab <- same
  ends <- list(P1 = list(k = p1_keep, l = p1_lab), S1 = list(k = s1_keep, l = s1_lab), REG = list(k = rep(TRUE, length(tr$idx)), l = tr$reg))
  scr <- function(S) { s <- S[tr$idx]; vapply(ends, function(e) ap_abs(s[e$k], e$l[e$k]), 0) }
  ops_of <- function(C) { K <- selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L)); list(K10 = K, R = recip(K)) }
  # criterion A on lab halves
  hA <- sim$lab <= 5L; CA <- clr_of(mi_of(X[, hA])); CB <- clr_of(mi_of(X[, !hA]))
  ut <- tr$idx; top1 <- function(C) { v <- C[ut]; v >= sort(v, decreasing = TRUE)[ceiling(0.01 * length(v))] }
  crit <- function(Ctr, target) { o <- ops_of(Ctr); sapply(o, function(A) { P <- lazy(A); Pt <- P; out <- numeric(0)
    for (t in seq_len(max(TGRID))) { if (t > 1) Pt <- Pt %*% P
      if (t %in% TGRID) out <- c(out, ap_abs(att(Pt)[ut], target)) }; out }) }
  cAB <- crit(CA, top1(CB)); cBA <- crit(CB, top1(CA)); rm(CA, CB)
  tstar <- sapply(colnames(cAB), function(op) { g <- sqrt(TGRID[which.max(cAB[, op])] * TGRID[which.max(cBA[, op])])
    TGRID[which.min(abs(log(TGRID) - log(g)))] })
  # full data
  C <- clr_of(mi_of(X)); ops <- ops_of(C)
  rows <- list(data.frame(op = "CLR", depth = "na", t = NA, t(scr(C))))
  sanity <- c()
  for (nm in names(ops)) { P <- lazy(ops[[nm]]); Pt <- P; need <- sort(unique(c(9, tstar[[nm]])))
    for (t in seq_len(max(need))) { if (t > 1) Pt <- Pt %*% P
      if (t %in% need) { S <- att(Pt); v <- scr(S)
        if (t == 9) { rows[[length(rows) + 1]] <- data.frame(op = nm, depth = "t9", t = 9, t(v)); sanity[nm] <- score(S, tr)$small_x }
        if (t == tstar[[nm]]) rows[[length(rows) + 1]] <- data.frame(op = nm, depth = "tstar", t = t, t(v)) } } }
  stored <- utils::read.csv(sprintf("analysis/synthetic_width/out/round2/runs/%s_metrics.csv", TAG))
  ref <- c(K10 = stored$small_x[stored$op == "K10|att_t9"], R = stored$small_x[stored$op == "RECIP|att_t9"])
  res <- cbind(gen = GEN, f = FPL, b = BAT, seed = SEED, do.call(rbind, rows),
               base_P1 = mean(p1_lab[p1_keep]), base_S1 = mean(s1_lab[s1_keep]), base_REG = mean(tr$reg),
               n_pos_P1 = sum(p1_lab), sanity_maxreldiff = max(abs(sanity[names(ref)] / ref - 1)))
  utils::write.csv(res, fout, row.names = FALSE)
  cat(sprintf("%s done in %.0f s; t* K10 %d R %d; sanity %.2e\n", TAG, as.numeric(difftime(Sys.time(), t0, units = "secs")),
              tstar[["K10"]], tstar[["R"]], res$sanity_maxreldiff[1]))
  TAG
}
invisible(parallel::mclapply(RUNS, function(r) tryCatch(one_run(r), error = function(e) { cat("ERROR", r, conditionMessage(e), "\n"); NA }),
                             mc.cores = WORKERS, mc.preschedule = FALSE))
