#!/usr/bin/env Rscript
# EXPLORATORY (docs/plans/explore_idea5_plan.md). Prep for the learned-attention test: data, lab split, CLR score
# matrices Z for the full data and for each lab half, tau0 per Z. No truth is read here.
# Usage (repo root): Rscript analysis/idea5/idea5_prep.R --out results/idea5/prep [--quick N] [--threads 8]
args <- commandArgs(TRUE)
opt <- list(out = "results/idea5/prep", quick = 0L, threads = 8L, split_seed = 20261112L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); if (!key %in% names(opt)) stop("unknown argument: ", args[k])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[k + 1]) else args[k + 1] }
for (k in c("quick", "threads", "split_seed")) opt[[k]] <- as.integer(opt[[k]])
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = file.path(opt$out, "prep.log"), append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/idea5/idea5_lib.R")
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
d <- idea5_load_m3d(opt$quick); X <- d$X; lab <- d$lab; G <- nrow(X); N <- ncol(X)
say(sprintf("%d genes x %d experiments, %d labs", G, N, length(unique(lab))))
set.seed(opt$split_seed); half <- idea5_lab_split(lab)
say(sprintf("lab halves: %d and %d experiments", sum(half == 1), sum(half == 2)))
utils::write.csv(data.frame(experiment = colnames(X), lab = lab, half = half), file.path(opt$out, "split.csv"), row.names = FALSE)
writeLines(rownames(X), file.path(opt$out, "genes.txt"))
idea5_write_f64(X, file.path(opt$out, "X.f64"))
sets <- list(full = rep(TRUE, N), h1 = half == 1, h2 = half == 2)
tau0 <- c()
for (nm in names(sets)) {
  t0 <- Sys.time(); f <- idea5_fit(X[, sets[[nm]]], opt$threads); Z <- f$clr_scores
  Tm <- idea5_row_temps(Z); tau0[nm] <- stats::median(Tm, na.rm = TRUE)
  idea5_write_f64(Z, file.path(opt$out, sprintf("Z_%s.f64", nm)))
  if (nm == "full") { M <- f$mi; utils::write.csv(data.frame(gene = rownames(X), mi_rowmean = (rowSums(M) - diag(M)) / (G - 1)), file.path(opt$out, "mi_rowmean_full.csv"), row.names = FALSE) }
  say(sprintf("%s: %d experiments, CLR in %.0f s; tau0 = %.5f (rows with a temperature: %d of %d; positive scores per row median %d)", nm, sum(sets[[nm]]),
              as.numeric(difftime(Sys.time(), t0, units = "secs")), tau0[nm], sum(is.finite(Tm)), G, as.integer(stats::median(rowSums(Z > 0)))))
  f$release(); rm(f, Z, Tm); invisible(gc(verbose = FALSE))
}
writeLines(jsonlite::toJSON(list(G = G, N = N, tau0 = as.list(tau0), split_seed = opt$split_seed, quick = opt$quick,
                                 n_half = c(sum(half == 1), sum(half == 2))), auto_unbox = TRUE, digits = 17), file.path(opt$out, "meta.json"))
say("prep done")
