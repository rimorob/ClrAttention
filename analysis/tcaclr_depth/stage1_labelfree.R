#!/usr/bin/env Rscript
# Stage 1 of docs/plans/tcaclr_depth_grid_plan.md (Revision 2): label-free criteria for every cell of the depth x width x clip grid on one
# development matrix. NO truth is read. EXPLORATORY.
#
#   Rscript analysis/tcaclr_depth/stage1_labelfree.R --matrix results/baselines_m3d/expression.tsv --kind lab --splits 8 --out results/tcaclr_depth/stage1_m3d
#   Rscript analysis/tcaclr_depth/stage1_labelfree.R --matrix results/baselines_k562ess_truthgenes_genie3/expression.tsv --kind random --splits 8 --out results/tcaclr_depth/stage1_k562
#
# --kind lab     the stored laboratory splits of results/lab_splits/split_assignments.csv (first --splits of them), halves = unions of laboratories
# --kind random  resampling units are knockdown CLUSTERS (Ward (ward.D2) hclust on 1 - r of the knockdown profiles, number of clusters =
#                round(N_eff of the full matrix) unless --clusters is given); each split shuffles the clusters (seed --seed + split) and gives each
#                to the half with fewer knockdowns so far (as lab_splits.R does with laboratories)
# Output: <out>/split_NN.csv (one per split, restartable), <out>/stage1_<dataset>.csv (merged), <out>/neff.csv, <out>/split_assignments.csv.
# Columns: dataset, split, direction, method (aCLR|tcaCLR), keff, cap, clipS, clipO, t, C1, C2, C2_null, C2_adj, C3, C3_e10, C3_bic, k_reach,
#   c3_unreach, n_eff. Direction 1 = half 1's score against half 2's targets, direction 2 the reverse. The reference operator
#   (keff 10, cap 50, clipped, t = 9) is included with cap = 50.
suppressPackageStartupMessages(source("analysis/tcaclr_depth/grid_lib.R"))
opt <- gl_options(commandArgs(TRUE), list(matrix = "", kind = "lab", splits = 8L, out = "", dataset = "", workers = 1L, threads = 4L, seed = 20261120L,
  lab_splits = "results/lab_splits/split_assignments.csv", clusters = 0L, genes = 0L, keffs = "", depths = "", reference = 1L))
if (!nzchar(opt$matrix) || !nzchar(opt$out)) stop("--matrix and --out are required"); stopifnot(opt$kind %in% c("lab", "random"))
if (!nzchar(opt$dataset)) opt$dataset <- if (opt$kind == "lab") "M3D" else "K562"
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, "stage1.log"))
say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))

X <- gl_read_matrix(opt$matrix); if (opt$genes > 0L) X <- X[seq_len(opt$genes), ]
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; G <- nrow(X); say(sprintf("%s: %d genes x %d experiments", opt$dataset, G, ncol(X)))
ne_all <- gl_neff(X); say(sprintf("full matrix: N_eff %.1f of %d experiments; G_eff %.1f of %d genes", ne_all[["N_eff"]], ne_all[["N"]], ne_all[["G_eff"]], G))

## ---- splits: halves made of whole resampling units ----------------------------------------------------------------------------
bal_split <- function(unit, seed) {           # shuffle units, give each to the half with fewer columns so far
  set.seed(seed); sz <- table(unit); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (u in names(sz)) { h <- which.min(tot); half[u] <- h; tot[h] <- tot[h] + sz[[u]] }
  unname(half[as.character(unit)])
}
if (opt$kind == "lab") {
  sa <- utils::read.csv(opt$lab_splits, check.names = FALSE, stringsAsFactors = FALSE); stopifnot(all(colnames(X) %in% sa$experiment), !anyDuplicated(sa$experiment))
  sa <- sa[match(colnames(X), sa$experiment), ]; sc <- grep("^split_", names(sa), value = TRUE)[seq_len(opt$splits)]
  splits <- lapply(sc, function(s) sa[[s]]); unit <- sa$lab
  say(sprintf("laboratory splits %s from %s; %d laboratories", paste(sc, collapse = " "), opt$lab_splits, length(unique(unit))))
} else {
  k <- if (opt$clusters > 0L) opt$clusters else max(2L, as.integer(round(ne_all[["N_eff"]])))
  hc <- stats::hclust(stats::as.dist(1 - stats::cor(X)), method = "ward.D2"); unit <- paste0("cl", stats::cutree(hc, k = k))
  say(sprintf("knockdown clusters: %d (Ward ward.D2 on 1 - r; N_eff-based count); sizes: median %g, max %d", k, stats::median(table(unit)), max(table(unit))))
  splits <- lapply(seq_len(opt$splits), function(s) bal_split(unit, opt$seed + s))
}
gl_wcsv(data.frame(experiment = colnames(X), unit = unit, do.call(cbind, stats::setNames(splits, sprintf("split_%02d", seq_along(splits))))), file.path(opt$out, "split_assignments.csv"))

## ---- operators ---------------------------------------------------------------------------------------------------------------
keffs <- if (nzchar(opt$keffs)) as.numeric(strsplit(opt$keffs, ",")[[1]]) else GL$keffs
depths <- if (nzchar(opt$depths)) as.integer(strsplit(opt$depths, ",")[[1]]) else GL$depths
ops <- gl_grid_ops(reference = opt$reference == 1L, keffs = keffs, depths = depths)
say(sprintf("%d operators (%d keff x 2 clipS%s), %d depths", length(ops), length(keffs), if (opt$reference == 1L) " + reference" else "", length(depths)))

run_split <- function(si) {
  f_out <- file.path(opt$out, sprintf("split_%02d.csv", si)); if (file.exists(f_out)) return(paste("split", si, "already done"))
  t0 <- Sys.time(); h <- splits[[si]]; set.seed(opt$seed + 1000L + si)
  halves <- lapply(1:2, function(d) { Xh <- X[, h == d]; rng <- apply(Xh, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) Xh[rng <= 0, ] <- Xh[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xh), sd = 1e-9)   # as lab_splits.R build_fit
    Xh })
  nes <- vapply(halves, function(Xh) gl_neff(Xh)[["N_eff"]], 0)
  fits <- lapply(halves, gl_fit, threads = opt$threads); rm(halves)
  r <- gl_stage1_split(fits[[1]], fits[[2]], ops, nes, say = function(...) say(sprintf("split %02d: ", si), ...))
  for (f in fits) f$release()
  r <- cbind(dataset = opt$dataset, split = si, r)
  gl_wcsv(r, f_out); sprintf("split %02d: %.1f min, peak RSS %.1f GB", si, as.numeric(difftime(Sys.time(), t0, units = "mins")), gl_vmhwm_gb())
}
todo <- seq_along(splits)
res <- if (opt$workers > 1L) parallel::mclapply(todo, function(s) tryCatch(run_split(s), error = function(e) paste("ERROR split", s, conditionMessage(e))), mc.cores = opt$workers, mc.preschedule = FALSE) else lapply(todo, run_split)
for (m in res) say(m)
fs <- list.files(opt$out, "^split_[0-9]+\\.csv$", full.names = TRUE)
if (length(fs)) { all <- do.call(rbind, lapply(fs, utils::read.csv)); gl_wcsv(all, file.path(opt$out, sprintf("stage1_%s.csv", opt$dataset))); say("merged ", length(fs), " splits, ", nrow(all), " rows") }
gl_wcsv(data.frame(what = "full matrix", t(ne_all)), file.path(opt$out, "neff.csv"))
if (any(grepl("^ERROR", unlist(res)))) quit(status = 1)
