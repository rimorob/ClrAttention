#!/usr/bin/env Rscript
# Runner for the label-free adaptive-setting studies on one development matrix (plan Revision 2, items 3, 8, 11). EXPLORATORY; no truth is read.
#   Rscript analysis/tcaclr_depth/adaptive/run_adaptive.R --dataset m3d  --out results/tcaclr_depth/adaptive_m3d
#   Rscript analysis/tcaclr_depth/adaptive/run_adaptive.R --dataset k562 --out results/tcaclr_depth/adaptive_k562
# Steps (each restartable: a step whose output file exists is skipped): (1) spectral_edge.R on the full matrix for every keff; (2) graph_signal_lik.R over the --splits
# blocked splits (M3D: stored laboratory splits; K562: k-means clusters of the knockdown profiles); (3) negative_control.R on --ncperms per-gene permuted copies.
# Negative controls use --nckeffs (default 10, 28, 80) and the leave-one-out edge (--seperms 0) on the first --ncsplits splits. Defaults: the 8 grid widths (7, 10, 14, 20, 28, 40, 56, 80), the 17 grid depths. Output: the files of the three scripts in --out, plus run_adaptive.log (time and peak RSS per step).
if (!exists("ad_Psym")) source("analysis/tcaclr_depth/adaptive/adaptive_lib.R")
for (f in c("spectral_edge", "graph_signal_lik", "negative_control")) source(sprintf("analysis/tcaclr_depth/adaptive/%s.R", f))
opt <- gl_options(commandArgs(TRUE), list(dataset = "m3d", matrix = "", kind = "", genes = 0L, splits = 8L, keffs = "", depths = "", r = 400L, perms = 10L, ncperms = 20L, seperms = 0L, nckeffs = "10,28,80", ncsplits = 2L,
  kobs = 300L, clusters = 0L, threads = 4L, workers = 4L, seed = AD$seed, out = "", lab_splits = "results/lab_splits/split_assignments.csv"))
if (!nzchar(opt$out)) stop("--out is required"); stopifnot(opt$dataset %in% c("m3d", "k562"))
if (!nzchar(opt$matrix)) opt$matrix <- if (opt$dataset == "m3d") "results/baselines_m3d/expression.tsv" else "results/baselines_k562ess_truthgenes_genie3/expression.tsv"
if (!nzchar(opt$kind)) opt$kind <- if (opt$dataset == "m3d") "lab" else "kmeans"
tag <- toupper(opt$dataset); dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, "run_adaptive.log"))
say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
keffs <- if (nzchar(opt$keffs)) as.numeric(strsplit(opt$keffs, ",")[[1]]) else GL$keffs
nckeffs <- as.numeric(strsplit(opt$nckeffs, ",")[[1]])
depths <- if (nzchar(opt$depths)) as.integer(strsplit(opt$depths, ",")[[1]]) else GL$depths
X <- ad_load(opt$matrix, opt$genes, opt$seed); say(sprintf("%s: %d genes x %d experiments", tag, nrow(X), ncol(X)))
units <- ad_units(X, opt$kind, opt$splits, opt$seed, opt$lab_splits, opt$clusters); say(sprintf("%d resampling units (%s), %d splits", units$k, opt$kind, length(units$splits)))
step <- function(name, file, expr) { if (file.exists(file.path(opt$out, file))) { say("skip ", name, " (output exists)"); return(invisible(NULL)) }
  t0 <- Sys.time(); force(expr); say(sprintf("%s: %.1f min, peak RSS of the main process %.2f GB", name, as.numeric(difftime(Sys.time(), t0, units = "mins")), ad_mem())) }
step("spectral edge", sprintf("spectral_edge_%s.tsv", tag), {
  res <- se_run(X, keffs, opt$perms, opt$kobs, 30L, opt$threads, opt$workers, opt$seed, tag, TRUE, say); se_write(res, opt$out, tag) })
step("graph-signal likelihood", sprintf("gs_selected_%s.tsv", tag), {
  res <- gs_run(X, units, keffs, depths, opt$r, opt$threads, opt$workers, say); gs_write(res, opt$out, tag) })
step("negative control", sprintf("nc_summary_%s.tsv", tag), {
  units_nc <- units; units_nc$splits <- units$splits[seq_len(min(opt$ncsplits, length(units$splits)))]
  res <- nc_run(X, units_nc, nckeffs, depths, opt$ncperms, opt$r, opt$seperms, 40L, opt$threads, opt$workers, opt$seed, say)
  w <- function(x, nm) utils::write.table(x, file.path(opt$out, sprintf("nc_%s_%s.tsv", nm, tag)), sep = "\t", quote = FALSE, row.names = FALSE); w(res$perm, "perm"); w(res$summary, "summary") })
say("done")
