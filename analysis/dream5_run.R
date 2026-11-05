#!/usr/bin/env Rscript
# DREAM5 pipeline driver (docs/plans/dream5_design.md, steps 4-5). Label-free: reads only the expression, chip
# features, gene ids and TF lists. One pipeline for nets 1, 3, 4 and all methods of the design's table.
#
#   Rscript analysis/dream5_run.R --net 4 --methods abs_pearson,clr,clrns,attention [--workers 16] [--input replicate]
#   Rscript analysis/dream5_run.R --net 4 --refit_from 1 --refit_to 100 --workers 16     # SE_B refits (hash sidecars; resumable)
#   Rscript analysis/dream5_run.R --net 4 --refit_unit series --refit_from 1 --refit_to 100 --workers 8   # supplementary series-level refits (net4 only)
#   Rscript analysis/dream5_run.R --collect 1                                              # consolidate SCORE_MANIFEST.tsv
#   Rscript analysis/dream5_run.R --net 1 --methods pidc --export_only 1                   # write expression file for the PIDC lane
#
# methods: abs_pearson, clr, clrns, attention, genie3_tf, grnboost2_tf, pidc, genie3_all, grnboost2_all
#   (the *_all runs use every gene as regulator; descriptive module readout; symmetrised by max at scoring time)
# Caches: results/dream5/net<n>/<method>__<input>.scores plus a .sha256 sidecar; existing caches are reused (resumable).
# Refit replicates keep hashes only unless --keep scores (see d5_refit_run in analysis/dream5_lib.R).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(net = "", methods = "", input = "replicate", workers = 16L, out = "results/dream5", collect = 0L,
            refit_from = 0L, refit_to = 0L, refit_unit = "group", keep = "none", export_only = 0L, force = 0L)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]; opt[[key]] <- if (is.integer(opt[[key]])) as.integer(v) else v; i <- i + 2L
}
source("analysis/dream5_lib.R"); D5$out <- opt$out
say <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

if (opt$collect == 1L) {
  r <- d5_collect_manifest(opt$out); say("manifest: ", nrow(r), " score files"); print(r[, c("net", "method", "input", "bytes", "seconds")])
  for (u in c("group", "series")) for (n in D5$nets) if (length(Sys.glob(file.path(d5_refit_dir(n, opt$out, u), "rep*.sha256.tsv")))) {
    x <- d5_refit_manifest(n, opt$out, u); say(u, "-level refit manifest net", n, ": ", length(unique(x$rep)), " replicates") }
  quit(save = "no")
}
net <- as.integer(opt$net); stopifnot(net %in% D5$nets)
say(sprintf("net%d input=%s workers=%d  R %s  RNG %s", net, opt$input, opt$workers, getRversion(), paste(RNGkind(), collapse = "/")))
nd <- d5_read_net(net, opt$input)
say(sprintf("%d genes x %d samples, %d TFs, %d replicate groups", nrow(nd$X), ncol(nd$X), length(nd$tf_ids), nd$n_groups))
dir.create(d5_dir(net, opt$out), recursive = TRUE, showWarnings = FALSE)
if (opt$export_only == 1L) { p <- file.path(d5_dir(net, opt$out), sprintf("expression__%s.tsv", opt$input)); d5_expr_tsv(nd$X, p); say("exported ", p); quit(save = "no") }

if (nzchar(opt$methods)) for (m in strsplit(opt$methods, ",")[[1]]) {
  d5_run_method(nd, m, opt$input, opt$workers, opt$out, force = opt$force == 1L); invisible(gc(verbose = FALSE))
}
if (opt$refit_to > 0L) {
  stopifnot(opt$input == "replicate", opt$refit_unit %in% c("group", "series"), opt$refit_unit == "group" || net == 4L)
  for (b in seq(max(1L, opt$refit_from), opt$refit_to)) {
    done <- file.exists(d5_refit_path(net, b, opt$out, opt$refit_unit)); t0 <- Sys.time()
    d5_refit_run(nd, b, opt$workers, opt$out, keep = opt$keep, unit = opt$refit_unit); invisible(gc(verbose = FALSE))
    say(sprintf("%s-level refit %d/%d %s (%.1f min)", opt$refit_unit, b, opt$refit_to, if (done) "reused" else "done", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  d5_refit_manifest(net, opt$out, opt$refit_unit)
}
say("done")
