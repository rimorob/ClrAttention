#!/usr/bin/env Rscript
# Sensitivity analysis (NOT preregistered): label-free depth t* for K10 attention on an EcoMAC array
# subset, by the same criterion as the frozen EcoMAC selection (analysis/keff_test.R, variant K10,
# criterion A): 20 balanced author splits (same generator and seed, authors of the subset), each
# direction's training-half K10 attention mass S_t scored by AP for the held-out half's top-1% CLR
# pairs, t grid {1,2,3,4,5,6,8,10,12,16,20,24,32}; per direction argmax, per split geometric mean of
# the two directions (rounded), t* = median over splits (floor(x + 0.5)). No truth is read.
# ComBat by platform is re-fitted within the subset (skipped for a single platform). Functions
# lab_split-style helpers are imported from the analysis scripts as keff_test.R does.
# Validation: on --subset A the curves must equal results/keff_test/ecomac/heldout_curves.csv (K10).
#   Rscript analysis/ecomac_nonoverlap_tstar.R --subset B [--threads 8] [--nsplit 20]
args <- commandArgs(trailingOnly = TRUE)
opt <- list(subset = "B", threads = 8L, nsplit = 20L, seed = 20261014L, alpha = 0.5, cap = 50L,
            target_frac = 0.01, tgrid = "1,2,3,4,5,6,8,10,12,16,20,24,32", out = "results/ecomac_nonoverlap/tstar")
i <- 1L
while (i <= length(args)) { key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]; opt[[key]] <- if (key %in% c("threads", "nsplit", "seed", "cap")) as.integer(v) else
    if (key %in% c("alpha", "target_frac")) as.numeric(v) else v; i <- i + 2L }
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, paste0("tstar_", opt$subset, ".log"))
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R")   # ap_fast
tgrid <- as.integer(strsplit(opt$tgrid, ",")[[1]])
import_defs <- function(file, nms) {
  got <- character(0)
  for (e in parse(file.path("analysis", file), keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% nms && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      eval(e, globalenv()); got <- c(got, as.character(e[[2]])) }
  if (length(setdiff(nms, got))) stop("not found in ", file, ": ", paste(setdiff(nms, got), collapse = ", "))
}
import_defs("lab_splits.R", c("build_fit", "lazy_P")); import_defs("ecomac_batch.R", "author_split")
nthr <- function() opt$threads
ai <- sort(unique(as.integer(readLines(file.path("results/ecomac_nonoverlap/subsets", paste0(opt$subset, ".txt"))))))
X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
rownames(X) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
author_all <- utils::read.delim("data/ecomac/samples.tsv", stringsAsFactors = FALSE, check.names = FALSE)[["Author"]]
X <- X[, ai, drop = FALSE]; plat <- plat[ai]; author <- author_all[ai]
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (length(unique(plat)) > 1L) X <- suppressMessages(sva::ComBat(X, batch = factor(plat)))
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; storage.mode(X) <- "double"
genes <- rownames(X); G <- length(genes)
say(sprintf("subset %s: %d genes x %d arrays; platforms %s", opt$subset, G, ncol(X),
            paste(sprintf("%s=%d", names(table(plat)), as.integer(table(plat))), collapse = "; ")))
set.seed(opt$seed)
splits <- lapply(seq_len(opt$nsplit), function(s) author_split())
selfloop <- function(A) { e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1; A }
soft_k <- function(C, k) selfloop(clr:::.softmax_rows(C, keff = k, cap = opt$cap))
walk_mass <- function(As, ts, fn) {
  Ps <- lapply(As, function(A) lazy_P(as.matrix(A))); rm(As); H <- length(Ps); Pts <- vector("list", H)
  for (t in seq_len(max(ts))) {
    for (h in seq_len(H)) Pts[[h]] <- if (t == 1L) as.matrix(Ps[[h]]) else as.matrix(Pts[[h]] %*% Ps[[h]])
    if (t %in% ts) { S <- if (H == 1L) Pts[[1]] else Reduce(`+`, Pts) / H
                     S <- (S + t(S)) / 2; diag(S) <- 0; fn(S, t); rm(S) }
    if (t %% 8L == 0L) invisible(gc(verbose = FALSE))
  }
  invisible(NULL)
}
ut <- which(upper.tri(matrix(0, G, G)))
rows <- list()
for (si in seq_len(opt$nsplit)) {
  t0 <- Sys.time(); h <- splits[[si]]
  Cs <- lapply(1:2, function(d) { set.seed(opt$seed + 100L * si + d); f <- build_fit(X[, h == d, drop = FALSE])
    C <- f$clr_scores; f$release(); C })
  for (dir in 1:2) {
    te <- 3L - dir; sv <- Cs[[te]][ut]; tA <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)); rm(sv)
    As <- list(Matrix::Matrix(soft_k(Cs[[dir]], 10), sparse = TRUE)); cr <- numeric(length(tgrid))
    walk_mass(As, tgrid, function(S, t) cr[match(t, tgrid)] <<- ap_fast(S[ut], tA))
    rows[[length(rows) + 1L]] <- data.frame(split = si, direction = dir, t = tgrid, crit_A = cr)
  }
  say(sprintf("split %02d (halves %d/%d arrays) %.1f min", si, sum(h == 1L), sum(h == 2L),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cur <- do.call(rbind, rows)
utils::write.csv(cur, file.path(opt$out, paste0("curves_", opt$subset, ".csv")), row.names = FALSE)
per_sd <- do.call(rbind, lapply(split(cur, list(cur$split, cur$direction), drop = TRUE), function(y) {
  y <- y[order(y$t), ]; k <- which.max(y$crit_A)
  data.frame(split = y$split[1], direction = y$direction[1], t_best = y$t[k], crit_best = y$crit_A[k],
             flat_lo = min(y$t[y$crit_A >= 0.99 * y$crit_A[k]]), flat_hi = max(y$t[y$crit_A >= 0.99 * y$crit_A[k]])) }))
per_s <- do.call(rbind, lapply(split(per_sd, per_sd$split), function(y)
  data.frame(split = y$split[1], t_split = as.integer(round(exp(mean(log(y$t_best))))))))
t_star <- as.integer(floor(stats::median(per_s$t_split) + 0.5))
utils::write.csv(per_sd, file.path(opt$out, paste0("per_split_direction_", opt$subset, ".csv")), row.names = FALSE)
mean_curve <- stats::aggregate(crit_A ~ t, cur, mean)
sel <- data.frame(subset = opt$subset, n_arrays = ncol(X), n_genes = G, nsplit = opt$nsplit, t_star = t_star,
                  t_q25 = stats::quantile(per_s$t_split, 0.25, names = FALSE), t_q75 = stats::quantile(per_s$t_split, 0.75, names = FALSE),
                  t_min = min(per_s$t_split), t_max = max(per_s$t_split), t_mean_curve_argmax = mean_curve$t[which.max(mean_curve$crit_A)],
                  flat_lo_median = stats::median(per_sd$flat_lo), flat_hi_median = stats::median(per_sd$flat_hi))
utils::write.csv(sel, file.path(opt$out, paste0("selection_", opt$subset, ".csv")), row.names = FALSE)
say("selection: ", paste(names(sel), unlist(sel), sep = "=", collapse = " "))
if (opt$subset == "A" && opt$nsplit == 20L) {
  ref <- utils::read.csv("results/keff_test/ecomac/heldout_curves.csv"); ref <- ref[ref$variant == "keff" & ref$param == 10, ]
  m <- merge(cur, ref[, c("split", "direction", "t", "crit_A")], by = c("split", "direction", "t"), suffixes = c("", "_ref"))
  say(sprintf("validation vs frozen K10 curves: %d rows, max |diff| = %.3g", nrow(m), max(abs(m$crit_A - m$crit_A_ref))))
}
say("done")
