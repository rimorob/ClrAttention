#!/usr/bin/env Rscript
# LABEL-FREE diagnostic on real M3D expression (466 replicate-averaged
# experiments, 4,297 genes, 39+ labs). No truth set (RegulonDB, complexes,
# operons) is read, and no network is scored. Question: on real data, what
# "width" do the background models of docs/plans/background_width_exploration.md
# produce, how many genes have no reproducible partner, and does the
# background-defined width track sample size (power) rather than biology?
#
# Quantities per gene:
#   Q_EB    sum of (1 - lfdr) under the graphon (X-by-Y) background
#   Q_EVT   sum of (1 - tail FDR) under the per-row EVT tail background
#   repro_n number of top-50 positive-CLR candidates shared by both halves of
#           one balanced lab split (seed 20261030; labs as in lab_splits.R)
# plus Q_EVT on random experiment subsets of 25%, 50% and 100% (seed 7).
# Usage: Rscript analysis/synthetic_width/m3d_label_free.R [m3d|k562]
# Run with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19, <= 4 GB.
suppressPackageStartupMessages(library(clr))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
DS <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "m3d"
stopifnot(DS %in% c("m3d", "k562"))
if (DS == "m3d") {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(normalizePath(m3d))
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                           check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"; rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
  exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description)
  exper[!grepl("experimenter:", desc$description)] <- NA
  lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])[colnames(X)]
  lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab))); rm(desc, exper)
} else {
  # K562 essential Perturb-seq pseudobulk (5,671 genes x 707 knockdown profiles,
  # half A, as in docs/plans/corum_analysis.md). Expression only; the truth files
  # in the same folder are not opened. No labs: each profile is its own "lab",
  # so the split below is a random half of the knockdown profiles.
  raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv",
                           check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"; rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  lab <- colnames(X)
}
G <- nrow(X); N <- ncol(X)
say(sprintf("%s: %d genes x %d experiments, %d labs", DS, G, N, length(unique(lab))))

# the same background functions as sim_width.R (G taken from this data)
mi_of <- function(X) bspline_mi(X, bins = "hg", spline_order = 3, threads = 1)
clr_of <- function(M) clr_calibrate(M, method = "normal", combine = "stouffer")
mirror_q <- function(u) {
  d <- stats::density(u, n = 4096, bw = "nrd0"); m0 <- d$x[which.max(d$y)]
  f <- function(x) stats::approx(d$x, d$y, x, rule = 2, yleft = 0, yright = 0)$y
  fx <- f(u); f0 <- f(2 * m0 - u)
  list(q = ifelse(u > m0 & fx > 0, 1 - pmin(1, f0 / fx), 0), pi0 = min(1, 2 * mean(u <= m0)))
}
q_graphon_rowsum <- function(M, nb = 20L) {
  Y <- M^(1 / 3); diag(Y) <- NA; th <- apply(Y, 1, stats::median, na.rm = TRUE)
  b <- cut(rank(th, ties.method = "first"), nb, labels = FALSE)
  idx <- which(upper.tri(Y)); i <- (idx - 1L) %% G + 1L; j <- (idx - 1L) %/% G + 1L
  y <- Y[idx]; rm(Y); cell <- pmin(b[i], b[j]) * 100L + pmax(b[i], b[j])
  med <- tapply(y, cell, stats::median); md <- tapply(y, cell, stats::mad)
  k <- match(cell, as.integer(names(med))); u <- (y - med[k]) / md[k]; mq <- mirror_q(u)
  list(Q = as.numeric(rowsum(c(mq$q, mq$q), c(i, j))), pi0 = mq$pi0)
}
q_evt_rowsum <- function(M) {
  Y <- M^(1 / 3); diag(Y) <- NA
  vapply(seq_len(G), function(i) { y <- Y[i, ]; y[i] <- -Inf
    qs <- stats::quantile(y[-i], c(0.90, 0.95), names = FALSE); beta <- (qs[2] - qs[1]) / log(2)
    o <- order(y, decreasing = TRUE)[1:200]
    fdr <- pmin(1, 0.1 * (G - 1) * exp(-(y[o] - qs[1]) / beta) / seq_along(o))
    sum(1 - rev(cummin(rev(fdr)))) }, 0)
}
top50 <- function(C) { T <- matrix(FALSE, G, G)
  for (i in seq_len(G)) { r <- C[i, ]; r[i] <- -Inf; o <- order(r, decreasing = TRUE)[1:50]; T[i, o[r[o] > 0]] <- TRUE }
  T }
qs <- function(v) paste(sprintf("%.1f", stats::quantile(v, c(0.1, 0.25, 0.5, 0.75, 0.9))), collapse = " / ")

out <- list()
M <- mi_of(X); say("full MI done")
gq <- q_graphon_rowsum(M); qe <- q_evt_rowsum(M)
out$full <- data.frame(quantity = c("Q_EB (graphon)", "Q_EVT"), N = N,
  q10_q25_med_q75_q90 = c(qs(gq$Q), qs(qe)), frac_below_1 = c(mean(gq$Q < 1), mean(qe < 1)),
  frac_above_30 = c(mean(gq$Q > 30), mean(qe > 30)), pi0 = c(gq$pi0, NA))
rm(M); invisible(gc(verbose = FALSE))

# one balanced lab split (as lab_splits.R's lab_split(), seed 20261030)
set.seed(20261030L)
sz <- table(lab); sz <- sz[sample(length(sz))]
half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
h <- unname(half[lab])
TA <- top50(clr_of(mi_of(X[, h == 1]))); TB <- top50(clr_of(mi_of(X[, h == 2])))
rep_n <- rowSums(TA & TB); rm(TA, TB); invisible(gc(verbose = FALSE))
say(sprintf("lab split: %d / %d experiments", sum(h == 1), sum(h == 2)))
out$repro <- data.frame(quantity = "repro_n (top-50 candidates shared across lab halves)", N = N,
  q10_q25_med_q75_q90 = qs(rep_n), frac_below_1 = mean(rep_n == 0), frac_above_30 = mean(rep_n > 30), pi0 = NA)
rb <- data.frame(repro_n_bin = c("0", "1-2", "3-5", "6-10", "11-20", ">20"),
  frac_genes = as.numeric(table(cut(rep_n, c(-1, 0, 2, 5, 10, 20, 50)))) / G)

# Q_EVT against sample size (random experiment subsets)
for (f in c(0.25, 0.5)) { set.seed(7L); cols <- sort(sample(N, round(f * N)))
  qv <- q_evt_rowsum(mi_of(X[, cols])); invisible(gc(verbose = FALSE))
  out[[paste0("evt", f)]] <- data.frame(quantity = "Q_EVT", N = length(cols), q10_q25_med_q75_q90 = qs(qv),
    frac_below_1 = mean(qv < 1), frac_above_30 = mean(qv > 30), pi0 = NA) }
res <- do.call(rbind, out)
dir.create("analysis/synthetic_width/out", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(res, sprintf("analysis/synthetic_width/out/%s_label_free_widths.csv", DS), row.names = FALSE)
utils::write.csv(rb, sprintf("analysis/synthetic_width/out/%s_label_free_repro_bins.csv", DS), row.names = FALSE)
print(res, row.names = FALSE); print(rb, row.names = FALSE)
say("done")
