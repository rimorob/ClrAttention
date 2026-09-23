#!/usr/bin/env Rscript
# 907-chip M3D compendium: point estimates plus replicate-aware bootstrap.
#
# The 907 chips are technical/biological replicates of 466 experiments (207
# experiments have 1 chip, 107 have 2, 135 have 3, 17 have 4-5). Treating chips
# as independent samples is pseudo-replication. Two resampling schemes, both
# recomputing MI/CLR/attention from scratch on every draw:
#   replicate  one randomly chosen chip per experiment (all 466 experiments):
#              technical-noise variability only; no pseudo-replication and
#              no duplicated columns.
#   cluster    466 experiments drawn WITH replacement (conditions resampled),
#              one chip per occurrence, taking a different replicate chip on
#              repeated draws while any remain (duplicate columns only when
#              an experiment is drawn more often than it has chips): condition
#              + technical variability -- the headline error bars.
# Point estimates are reported on (i) all 907 chips (as in the 2007 paper) and
# (ii) the 466 replicate-averaged experiments (analysis/followup.R).
#
# Methods (fixed in advance; depths are the held-out-lab t* from
# results/avg_followup): |Pearson|, raw MI, CLR parity 2007 (10 bins,
# Euclidean), CLR HG-Stouffer (default), soft10 attention at t = 4, 8*, 16,
# pearson_top10 attention at t = 6*.
#
# Usage (repo root):
#   Rscript analysis/chips_bootstrap.R [--Bcluster 40] [--Brep 20]
#       [--out results/chips_bootstrap] [--threads N] [--quick 0]
# Results are appended per draw to draws.csv, so a partial run is usable.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            Bcluster = 40L, Brep = 20L, out = "results/chips_bootstrap",
            threads = NULL, quick = 0L, seed = 20260924L, alpha = 0.5)
ints <- c("Bcluster", "Brep", "threads", "quick", "seed")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key %in% ints) as.integer(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "bootstrap.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))
source(file.path(script_dir, "regulons.R"))
set.seed(opt$seed)

## ---- data ------------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
rd <- function(f) {
  raw <- utils::read.delim(file.path(opt$m3d, f), check.names = FALSE,
                           stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
  rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
  X
}
Xc <- rd(paste0(base, "_chips907probes4297.tab"))
Xa <- rd(paste0("avg_", base, "_exps466probes4297.tab"))
stopifnot(identical(rownames(Xc), rownames(Xa)))
keep <- apply(Xc, 1, function(x) diff(range(x))) > 0 &
  apply(Xa, 1, function(x) diff(range(x))) > 0
if (opt$quick > 0L) {
  v <- apply(Xa, 1, stats::var); v[!keep] <- -Inf
  keep <- seq_len(nrow(Xa)) %in% order(-v)[seq_len(opt$quick)]
}
Xc <- Xc[keep, ]; Xa <- Xa[keep, ]
genes <- rownames(Xc); G <- length(genes)
desc <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_descriptions")),
                          stringsAsFactors = FALSE)
exp_of <- stats::setNames(desc$experiment_name, desc$chip_name)[colnames(Xc)]
stopifnot(!anyNA(exp_of))
chips_by_exp <- split(colnames(Xc), exp_of)
say(sprintf("%d genes; %d chips of %d experiments (replicates per experiment: %s)",
            G, ncol(Xc), length(chips_by_exp),
            paste(names(table(lengths(chips_by_exp))), table(lengths(chips_by_exp)),
                  sep = "x", collapse = " ")))

to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)
bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
  ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
  reg <- make_regulons(rtab, genes, conf_keep = ck)
  p <- comembership_pairs(reg, operon_of, genes)
  list(reg = reg, pairs = p, co = comember_matrix(reg, G),
       idx = cbind(p$i, p$j))
})

## ---- one full evaluation of a sample matrix -----------------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
attention <- function(A, ts) {
  P <- Matrix::Matrix((1 - opt$alpha) * diag(G) + opt$alpha * A, sparse = TRUE)
  out <- list(); Pt <- diag(G)
  for (t in seq_len(max(ts))) {
    Pt <- as.matrix(Pt %*% P)
    if (t %in% ts) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; out[[as.character(t)]] <- S }
  }
  out
}
evaluate <- function(X, label) {
  t0 <- Sys.time()
  rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0)) {  # a gene constant in this draw: tiny jitter, logged
    X[rng <= 0, ] <- X[rng <= 0, ] + matrix(stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9),
                                            sum(rng <= 0))
  }
  Z <- Zs(X); R <- tcrossprod(Z) / (ncol(Z) - 1)
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none",
                                       threads = opt$threads)
  f$calibrate(method = "normal", combine = "stouffer")
  M <- f$mi; diag(M) <- 0
  par <- ClrAttention$new(X)$estimate_mi(bins = 10, transform = "none",
                                         threads = opt$threads)$
    calibrate(method = "normal", combine = "euclidean")$clr_scores
  f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
  sa <- attention(f$operator, c(4, 8, 16))
  Ar <- abs(R); diag(Ar) <- 0
  f$set_operator(clr:::.topk_rows(Ar, 10L), alpha = opt$alpha, label = "pearson_top10")
  pa <- attention(f$operator, 6)
  S <- list(abs_pearson = abs(R), mi = M, clr_parity2007 = par,
            clr_hg_stouffer = f$clr_scores,
            soft10_t4 = sa[["4"]], soft10_t8 = sa[["8"]], soft10_t16 = sa[["16"]],
            pearson_top10_t6 = pa[["6"]])
  S <- lapply(S, function(s) { diag(s) <- 0; s })
  rows <- list()
  for (ev in names(bench)) {
    b <- bench[[ev]]
    for (m in names(S)) {
      ap <- pr_summary(S[[m]][b$idx], b$pairs$label)[["aupr"]]
      coh <- stats::median(regulon_coherence(S[[m]], b$reg, b$co, operon_of, genes),
                           na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(draw = label, evidence = ev,
                                              method = m, aupr = ap,
                                              coherence = coh)
    }
  }
  out <- do.call(rbind, rows)
  say(sprintf("%-14s N=%d  %.0f s  | SC AUPR: CLR %.4f  soft10_t8 %.4f  pearson_t6 %.4f",
              label, ncol(X), as.numeric(difftime(Sys.time(), t0, units = "secs")),
              out$aupr[out$evidence == "SC" & out$method == "clr_hg_stouffer"],
              out$aupr[out$evidence == "SC" & out$method == "soft10_t8"],
              out$aupr[out$evidence == "SC" & out$method == "pearson_top10_t6"]))
  out
}
draws_f <- file.path(opt$out, "draws.csv")
append_rows <- function(d) utils::write.table(d, draws_f, sep = ",", row.names = FALSE,
                                             col.names = !file.exists(draws_f),
                                             append = file.exists(draws_f))
done <- if (file.exists(draws_f)) unique(utils::read.csv(draws_f)$draw) else character()

## ---- point estimates ---------------------------------------------------------
for (pe in list(list("point_chips907", Xc), list("point_avg466", Xa)))
  if (!pe[[1]] %in% done) append_rows(evaluate(pe[[2]], pe[[1]]))

## ---- resampling ----------------------------------------------------------------
draw_cols <- function(mode) {
  ex <- names(chips_by_exp)
  pick <- if (mode == "replicate") ex else sample(ex, length(ex), replace = TRUE)
  used <- list()
  vapply(pick, function(e) {
    ch <- chips_by_exp[[e]]
    avail <- setdiff(ch, used[[e]])
    if (!length(avail)) avail <- ch
    c1 <- if (length(avail) == 1L) avail else sample(avail, 1L)
    used[[e]] <<- c(used[[e]], c1)
    c1
  }, "")
}
for (mode in c("cluster", "replicate")) {
  B <- if (mode == "cluster") opt$Bcluster else opt$Brep
  for (b in seq_len(B)) {
    lab <- sprintf("%s_%03d", mode, b)
    set.seed(opt$seed + 1000L * (mode == "replicate") + b)
    cols <- draw_cols(mode)          # drawn before any skip: seeds stay aligned
    if (lab %in% done) next
    append_rows(evaluate(Xc[, cols, drop = FALSE], lab))
  }
}

## ---- summary ---------------------------------------------------------------
d <- utils::read.csv(draws_f)
d$mode <- sub("_.*$", "", d$draw)
ref <- d[d$method == "clr_hg_stouffer", c("draw", "evidence", "aupr", "coherence")]
names(ref)[3:4] <- c("aupr_ref", "coh_ref")
d <- merge(d, ref, by = c("draw", "evidence"))
d$d_aupr <- d$aupr - d$aupr_ref; d$d_coh <- d$coherence - d$coh_ref
q <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)
summ <- do.call(rbind, lapply(split(d, list(d$mode, d$evidence, d$method), drop = TRUE),
  function(x) {
    boot <- x$mode[1] %in% c("cluster", "replicate")
    data.frame(mode = x$mode[1], evidence = x$evidence[1], method = x$method[1],
               n = nrow(x), aupr = mean(x$aupr),
               aupr_lo = if (boot) q(x$aupr, .025) else NA,
               aupr_hi = if (boot) q(x$aupr, .975) else NA,
               d_aupr = mean(x$d_aupr),
               d_aupr_lo = if (boot) q(x$d_aupr, .025) else NA,
               d_aupr_hi = if (boot) q(x$d_aupr, .975) else NA,
               coherence = mean(x$coherence),
               d_coh = mean(x$d_coh),
               d_coh_lo = if (boot) q(x$d_coh, .025) else NA,
               d_coh_hi = if (boot) q(x$d_coh, .975) else NA)
  }))
summ <- summ[order(summ$mode, summ$evidence, -summ$aupr), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
for (m in unique(summ$mode)) for (ev in c("SC", "all")) {
  s <- summ[summ$mode == m & summ$evidence == ev, ]
  say(sprintf("== %s [%s], n = %d", m, ev, s$n[1]))
  for (k in seq_len(nrow(s)))
    say(sprintf("   %-18s AUPR %.4f  dAUPR vs CLR %+.4f [%+.4f, %+.4f]  coherence %.3f  dCoh %+.3f [%+.3f, %+.3f]",
                s$method[k], s$aupr[k], s$d_aupr[k], s$d_aupr_lo[k], s$d_aupr_hi[k],
                s$coherence[k], s$d_coh[k], s$d_coh_lo[k], s$d_coh_hi[k]))
}
say("done")
