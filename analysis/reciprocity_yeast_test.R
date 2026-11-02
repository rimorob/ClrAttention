#!/usr/bin/env Rscript
# Test script for docs/preregistration_reciprocity_yeast.md (FIXED 2026-11-01, commit 02cbcfd).
# THIS VERSION CONTAINS THE LABEL-FREE PARTS ONLY (protocol steps 5 and 6). The truth phase (step 10) is
# added before step 7; until then --truth refuses. No mode of this version reads any truth file.
#
# Modes:
#   --check_R                 protocol step 5: (i) unit test of the operator B on a hand-built 6-gene example with an
#                             underflow case; (ii) synthetic seed std_f10_b1_s1: the RECIP|att_t9 metrics of
#                             analysis/synthetic_width/out/round2/runs/std_f10_b1_s1_metrics.csv reproduced to relative
#                             1e-8 with the preregistered candidate-set operator; positivity vs candidate-set support
#                             compared and logged. Writes docs/reciprocity_yeast_selections/R_check.csv.
#   --select yeast|k562       protocol step 6: label-free t* for K10 and R (criterion A, 10 random halves, both
#                             directions, grid 1..32; per split-direction argmax, per split geometric mean, median over
#                             splits snapped to the grid), label-free diagnostics, and for yeast the SE_B proxy (first 50
#                             experiment-jackknife replicates). Writes <ds>_selections.csv (+ curves, per-row diagnostics).
#   --diagnostics m3d         label-free diagnostics only (no t*), as listed in the preregistration.
# Options: --workers 5 --threads 3
invisible(Sys.setlocale("LC_COLLATE", "C")); RNGkind("Mersenne-Twister", "Inversion", "Rejection")
suppressPackageStartupMessages({ library(clr); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
opt <- list(workers = 5L, threads = 3L)
i <- 1L; mode <- NULL; ds <- NULL
while (i <= length(args)) { a <- args[i]
  if (a == "--check_R") { mode <- "check_R"; i <- i + 1L; next }
  if (a %in% c("--select", "--diagnostics", "--truth")) { mode <- sub("^--", "", a); ds <- args[i + 1L]; i <- i + 2L; next }
  opt[[sub("^--", "", a)]] <- as.integer(args[i + 1L]); i <- i + 2L }
if (is.null(mode)) stop("usage: --check_R | --select yeast|k562 | --diagnostics m3d")
if (mode == "truth") stop("truth phase not yet written (protocol step 7); this version is label-free only")
SEL <- "docs/reciprocity_yeast_selections"; D <- "data/yeast"; OUTD <- "results/reciprocity_yeast"
dir.create(SEL, showWarnings = FALSE, recursive = TRUE); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
env_info <- function() data.frame(R = R.version.string, RNG = paste(RNGkind(), collapse = "/"), collate = Sys.getlocale("LC_COLLATE"),
                                  script_sha256 = digest::digest(file = "analysis/reciprocity_yeast_test.R", algo = "sha256"))
TGRID <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32)

## ---- operators (preregistration "Operators") ------------------------------------------------
selfloop <- function(A) { e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1; A }
cand_mask <- function(C, cap = 50L) {   # psi_ij = 1[j in K_i]: up to `cap` largest positive CLR scores of row i
  n <- nrow(C); idx <- lapply(seq_len(n), function(i) { r <- C[i, ]; r[i] <- -Inf
    o <- order(r, decreasing = TRUE)[seq_len(min(cap, n - 1L))]; o[r[o] > 0] })
  Matrix::sparseMatrix(i = rep(seq_len(n), lengths(idx)), j = unlist(idx), x = TRUE, dims = c(n, n)) }
K10_of <- function(C) selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L))
R_of <- function(A, psi) { M <- as.matrix(psi & Matrix::t(psi)); B <- A * M; diag(B) <- 0; B <- selfloop(B)
  diag(B) <- diag(B) + pmax(0, 1 - rowSums(B)); B }
lazyP <- function(B) { P <- 0.5 * B; diag(P) <- diag(P) + 0.5; P }
att <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }
ap_end <- function(s, l) { o <- order(s, decreasing = TRUE, method = "radix"); s <- s[o]; l <- l[o]; n <- length(s)
  if (!n || !sum(l)) return(NA_real_); last <- c(s[-1L] != s[-n], TRUE); tp <- cumsum(l)[last]; k <- which(last)
  sum(diff(c(0, tp)) * tp / k) / sum(l) }   # tie blocks counted at their end, as ap_fast
fitC <- function(X, threads) { rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0)) X[rng <= 0, ] <- X[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9)
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads)
  f$calibrate(method = "normal", combine = "stouffer"); C <- f$clr_scores; dimnames(C) <- NULL; C }

## ---- mode: check_R (step 5) -------------------------------------------------------------------
if (mode == "check_R") {
  rows <- list(); add <- function(check, value, tol, passed) rows[[length(rows) + 1L]] <<- data.frame(check, value, tol, passed)
  # (i) hand-built example: 6 genes; gene 1's attention to gene 2 is returned only through a candidate whose
  # softmax weight underflows (gene 2 -> gene 1). Positivity would drop 1 -> 2; the candidate-set rule keeps it.
  A <- matrix(0, 6, 6); A[1, 2] <- 0.6; A[1, 3] <- 0.4; A[2, 1] <- 0; A[2, 4] <- 1; A[3, 1] <- 0.5; A[3, 5] <- 0.5
  A[4, 2] <- 1; A[5, 6] <- 1; A[6, 6] <- 1
  psi <- Matrix::Matrix(FALSE, 6, 6, sparse = TRUE); psi[1, 2] <- psi[1, 3] <- psi[2, 1] <- psi[2, 4] <- psi[3, 1] <- psi[3, 5] <- TRUE
  psi[4, 2] <- psi[5, 6] <- TRUE
  B <- R_of(A, psi)
  Bexp <- matrix(0, 6, 6); Bexp[1, 2] <- 0.6; Bexp[1, 3] <- 0.4; Bexp[2, 4] <- 1; Bexp[3, 1] <- 0.5; Bexp[3, 3] <- 0.5
  Bexp[4, 2] <- 1; Bexp[5, 5] <- 1; Bexp[6, 6] <- 1
  add("unit_B_matches_expected", max(abs(B - Bexp)), 1e-15, max(abs(B - Bexp)) <= 1e-15)
  add("unit_rows_sum_to_1", max(abs(rowSums(B) - 1)), 1e-15, max(abs(rowSums(B) - 1)) <= 1e-15)
  Bpos <- A * (t(A) > 0); diag(Bpos) <- 0; Bpos <- selfloop(Bpos); diag(Bpos) <- diag(Bpos) + pmax(0, 1 - rowSums(Bpos))
  add("unit_underflow_case_differs_from_positivity", abs(B[1, 2] - Bpos[1, 2]), 0, B[1, 2] == 0.6 && Bpos[1, 2] == 0)
  # (ii) synthetic reproduction
  src <- readLines("analysis/synthetic_width/sim_width.R"); src <- src[seq_len(grep("^## ---- run", src) - 1L)]
  src <- src[!grepl("^args <-|^SEEDS <-|^SCEN <-|^NSCALE <-|^OUT <-|^dir.create", src)]
  SCEN <- "messy"; eval(parse(text = src))
  r2 <- readLines("analysis/synthetic_width/sim_round2.R")
  eval(parse(text = r2[grep("^simulate2 <- function", r2):(grep("^# Efron central matching", r2) - 1L)]))
  GEN <- "std"; FPL <- 0.1; BAT <- 1
  sim <- simulate2(20261031L + 1L); X <- sim$X; tr <- make_truth(sim)
  istf <- seq_len(G) %in% sim$tfs; Ed <- matrix(FALSE, G, G); Ed[cbind(sim$edges$tf, sim$edges$target)] <- TRUE
  Ed <- Ed | t(Ed); ii <- (tr$idx - 1L) %% G + 1L; jj <- (tr$idx - 1L) %/% G + 1L
  tr$tfpair <- xor(istf[ii], istf[jj]); tr$edge <- Ed[tr$idx]; rm(Ed)
  C <- clr_of(mi_of(X)); Ak <- selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L))
  psi <- cand_mask(C); Bp <- R_of(Ak, psi); Bs <- recip(Ak)
  add("synthetic_candidate_vs_positivity_maxdiff", max(abs(Bp - Bs)), 0, max(abs(Bp - Bs)) == 0)
  add("synthetic_underflow_entries", sum(as.matrix(psi) & Ak == 0), NA, TRUE)
  S <- att(walk9(lazy(Bp))); b <- score(S, tr); s <- S[tr$idx]; s[!is.finite(s)] <- -Inf; k <- tr$tfpair
  o <- order(s[k], decreasing = TRUE, method = "radix"); b$edge_x <- ap_sorted(s[k][o], tr$edge[k][o]) / mean(tr$edge[k])
  ref <- utils::read.csv("analysis/synthetic_width/out/round2/runs/std_f10_b1_s1_metrics.csv")
  ref <- ref[ref$op == "RECIP|att_t9", ]
  for (m in c("cpx_x", "cpx_noribo_x", "s5_x", "s20_x", "s80_x", "small_x", "mid_x", "reg_x", "sat_auc", "edge_x")) {
    rel <- abs(b[[m]] / ref[[m]] - 1); add(paste0("synthetic_", m), rel, 1e-8, is.finite(rel) && rel <= 1e-8) }
  out <- cbind(do.call(rbind, rows), env_info())
  utils::write.csv(out, file.path(SEL, "R_check.csv"), row.names = FALSE)
  print(out[, 1:4]); say("R check: ", if (all(out$passed)) "PASSED" else "FAILED (no verdict; fix and log)")
  quit(save = "no", status = if (all(out$passed)) 0L else 1L)
}

## ---- data ---------------------------------------------------------------------------------------
load_ds <- function(ds) {
  if (ds == "yeast") { x <- utils::read.delim(file.path(D, "deleteome_responsive_700.tsv"), check.names = FALSE)
    want <- strsplit(readLines(file.path(SEL, "deleteome_matrix.sha256"), 1L), "[[:space:]]+")[[1]][1]
    if (!identical(want, digest::digest(file = file.path(D, "deleteome_responsive_700.tsv"), algo = "sha256"))) stop("matrix SHA-256 mismatch")
  } else if (ds == "k562") x <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE)
  else if (ds == "m3d") x <- utils::read.delim("data/E_coli_v4_Build_6/avg_E_coli_v4_Build_6_exps466probes4297.tab", check.names = FALSE)
  X <- as.matrix(x[, -1]); storage.mode(X) <- "double"; rownames(X) <- x[[1]]
  X[apply(X, 1, function(v) diff(range(v))) > 0, , drop = FALSE] }

diagnostics <- function(X, C, threads) {     # preregistration "Label-free diagnostics"
  G <- nrow(X); A <- K10_of(C); psi <- cand_mask(C); B <- R_of(A, psi)
  offA <- rowSums(A) - diag(A); offB <- rowSums(B) - diag(B)
  removed_row <- ifelse(offA > 0, 1 - offB / offA, NA)
  Xs <- (X - rowMeans(X)) / apply(X, 1, stats::sd); pc <- stats::prcomp(t(Xs), rank. = 10)$x
  set.seed(20261031L); cl <- stats::kmeans(pc, centers = 5, nstart = 25, iter.max = 100)$cluster
  top50 <- function(C) { T <- matrix(FALSE, G, G)
    for (i in seq_len(G)) { r <- C[i, ]; r[i] <- -Inf; o <- order(r, decreasing = TRUE)[1:50]; T[i, o[r[o] > 0]] <- TRUE }; T }
  S <- numeric(G)
  for (k in 1:5) { S <- S + rowSums(top50(fitC(X[, cl != k], threads)) & top50(fitC(X[, cl == k], threads))); invisible(gc()) }
  list(summary = data.frame(G = G, N = ncol(X), kmeans_sizes = paste(table(cl), collapse = "/"),
         mass_removed_share = 1 - sum(offB) / sum(offA),
         removed_row_q10 = stats::quantile(removed_row, 0.1, na.rm = TRUE), removed_row_median = stats::median(removed_row, na.rm = TRUE),
         removed_row_q90 = stats::quantile(removed_row, 0.9, na.rm = TRUE), pure_selfloop_rows_R = mean(offB == 0),
         pure_selfloop_rows_K10 = mean(offA == 0), underflow_entries = sum(as.matrix(psi) & A == 0),
         repro_S_median = stats::median(S), spearman_removed_vs_repro = stats::cor(removed_row, S, method = "spearman", use = "complete.obs")),
       per_row = data.frame(gene = rownames(X), mass_removed = removed_row, repro_S = S)) }

## ---- mode: diagnostics (m3d) -----------------------------------------------------------------
if (mode == "diagnostics") {
  X <- load_ds(ds); C <- fitC(X, opt$threads * opt$workers); dg <- diagnostics(X, C, opt$threads * opt$workers)
  utils::write.csv(cbind(dataset = ds, dg$summary, env_info()), file.path(SEL, paste0(ds, "_label_free_diagnostics.csv")), row.names = FALSE)
  utils::write.csv(dg$per_row, file.path(OUTD, paste0(ds, "_diagnostics_per_row.csv")), row.names = FALSE)
  print(t(dg$summary)); quit(save = "no")
}

## ---- mode: select (step 6) --------------------------------------------------------------------
if (mode == "select") {
  X <- load_ds(ds); G <- nrow(X); n <- ncol(X); say(sprintf("%s: %d genes x %d samples", ds, G, n))
  if (ds == "yeast") { sp <- utils::read.csv(file.path(SEL, "splits.csv")); stopifnot(identical(sp$mutant, colnames(X)))
    halves <- lapply(1:10, function(s) sp[[sprintf("split_%02d", s)]])
  } else { set.seed(20261101); halves <- lapply(1:10, function(s) { h <- integer(n); h[sample.int(n, floor(n / 2))] <- 1L; h + 1L })
    utils::write.csv(data.frame(sample = colnames(X), do.call(cbind, stats::setNames(halves, sprintf("split_%02d", 1:10)))),
                     file.path(SEL, paste0(ds, "_splits.csv")), row.names = FALSE) }
  ut <- which(upper.tri(matrix(0, G, G)))
  curve_of <- function(C, target) { A <- K10_of(C); B <- R_of(A, cand_mask(C))
    vapply(list(K10 = A, R = B), function(M) { P <- lazyP(M); Pt <- P; out <- numeric(0)
      for (t in seq_len(max(TGRID))) { if (t > 1) Pt <- Pt %*% P; if (t %in% TGRID) out <- c(out, ap_end(att(Pt)[ut], target)) }
      out }, numeric(length(TGRID))) }
  top1 <- function(C) { v <- C[ut]; as.integer(v >= stats::quantile(v, 0.99)) }
  run_split <- function(s) { h <- halves[[s]]; C1 <- fitC(X[, h == 1], opt$threads); C2 <- fitC(X[, h == 2], opt$threads)
    c12 <- curve_of(C1, top1(C2)); c21 <- curve_of(C2, top1(C1)); rm(C1, C2); invisible(gc())
    rbind(data.frame(split = s, direction = 1, t = TGRID, K10 = c12[, "K10"], R = c12[, "R"]),
          data.frame(split = s, direction = 2, t = TGRID, K10 = c21[, "K10"], R = c21[, "R"])) }
  cf <- file.path(OUTD, paste0(ds, "_criterionA_curves.csv"))
  if (file.exists(cf)) cur <- utils::read.csv(cf) else {
    cur <- do.call(rbind, parallel::mclapply(1:10, run_split, mc.cores = opt$workers, mc.preschedule = FALSE))
    utils::write.csv(cur, cf, row.names = FALSE) }
  snap <- function(x) TGRID[which.min(abs(log(TGRID) - log(x)))]
  sel_of <- function(op) {
    best <- do.call(rbind, lapply(split(cur, list(cur$split, cur$direction), drop = TRUE), function(y)
      data.frame(split = y$split[1], direction = y$direction[1], t_best = y$t[which.max(y[[op]])])))
    per_s <- aggregate(t_best ~ split, best, function(v) exp(mean(log(v))))
    data.frame(operator = op, t_star = snap(stats::median(per_s$t_best)), t_split_q25 = stats::quantile(per_s$t_best, .25),
               t_split_median = stats::median(per_s$t_best), t_split_q75 = stats::quantile(per_s$t_best, .75),
               t_splits = paste(round(per_s$t_best, 2), collapse = "/")) }
  sel <- rbind(sel_of("K10"), sel_of("R"))
  sel$flag_grid_edge <- sel$t_star %in% range(TGRID)
  sel$flag_tR_gt_2tK10 <- sel$t_star[sel$operator == "R"] > 2 * sel$t_star[sel$operator == "K10"]
  say("t*: K10 ", sel$t_star[1], ", R ", sel$t_star[2])
  C <- fitC(X, opt$threads * opt$workers); dg <- diagnostics(X, C, opt$threads * opt$workers); rm(C); invisible(gc())
  utils::write.csv(dg$per_row, file.path(OUTD, paste0(ds, "_diagnostics_per_row.csv")), row.names = FALSE)
  seb <- NA_real_; seb_n <- 0L
  if (ds == "yeast") {   # SE_B proxy: first 50 experiment-jackknife replicates (preregistration "Power")
    rp <- utils::read.csv(file.path(SEL, "experiment_jackknife_replicates.csv")); stopifnot(identical(rp$mutant, colnames(X)))
    tK <- sel$t_star[1]; tR <- sel$t_star[2]
    one <- function(b) { keep <- rp[[sprintf("rep_%03d", b)]] == 1L
      Ck <- fitC(X[, keep], opt$threads); Cd <- fitC(X[, !keep], opt$threads); tg <- top1(Cd); rm(Cd)
      A <- K10_of(Ck); B <- R_of(A, cand_mask(Ck)); rm(Ck)
      wt <- function(M, t) { P <- lazyP(M); Pt <- P; if (t > 1) for (k in 2:t) Pt <- Pt %*% P; att(Pt)[ut] }
      d <- ap_end(wt(B, tR), tg) - ap_end(wt(A, tK), tg); invisible(gc()); d }
    dd <- unlist(parallel::mclapply(1:50, one, mc.cores = opt$workers, mc.preschedule = FALSE))
    utils::write.csv(data.frame(rep = 1:50, d_critA_R_minus_K10 = dd), file.path(OUTD, "yeast_seb_proxy_reps.csv"), row.names = FALSE)
    seb <- stats::sd(dd); seb_n <- sum(is.finite(dd)) }
  out <- cbind(dataset = ds, sel, dg$summary, seb_proxy_sd_critA = seb, seb_proxy_reps = seb_n,
               split_seed = if (ds == "yeast") "20261101 (yeast_prep splits)" else "20261101",
               jackknife_seeds = "P1 20261102; P2/H1 20261112; S1 20261122; S6 20261132; experiment 20262000+b",
               env_info(), row.names = NULL)
  utils::write.csv(out, file.path(SEL, paste0(ds, "_selections.csv")), row.names = FALSE)
  print(t(out[, 1:8])); say("selections written: ", file.path(SEL, paste0(ds, "_selections.csv")), " (commit and push)")
}
