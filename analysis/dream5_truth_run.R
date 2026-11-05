#!/usr/bin/env Rscript
# DREAM5 TRUTH-PHASE driver (docs/plans/dream5_design.md, protocol step 7). NOT TO BE RUN BEFORE THE TRUTH STAGE.
# It refuses to run unless (1) the design is FIXED, tracked, unchanged and pushed, (2) the score manifest and the refit
# manifest are committed, pushed and unchanged, and (3) every score file still has its manifest SHA-256. The guards run
# before the truth file is opened.
#
# The truth-stage details (cap ties, TF-TF pairs, pair universe, jackknife cap) have no defaults: they are passed on the
# command line after the lead's dated change-log entry resolves them from the official scoring scripts.
#
#   Rscript analysis/dream5_truth_run.R --net 4 --truth_edges <file> --out_dir results/dream5/truth_net4 \
#       --ties_at_cap keep_tied|drop_tied|list_order --ties_in_list end_of_block|expected --tftf keep_both|collapse_max|drop \
#       --universe all|gs [--universe_tf <ids file> --universe_gene <ids file> --cap_order cap_then_universe|universe_then_cap] \
#       --jackknife_cap reapply_cap|scaled_cap|restrict_after_cap [--structure_only 1] [--workers 16] [--nrep 200]
#   --truth_edges: TSV, no header; two columns (regulator gene ID, target gene ID; every row a positive edge) or three with a
#   0/1 label last (only label-1 rows are kept, the dropped ones are reported); self-edges are dropped and reported. The
#   extraction from the official gold-standard files is part of the truth stage and is not done here.
#   --structure_only 1: only the structure counts and the Y-M evaluability gate (no score is read). All switches are still
#   required and validated (a switch or file that would be ignored is an error, e.g. --cap_order or --universe_tf with
#   --universe all).
#   --series_n 0: no supplementary series-level jackknife (its manifest is then not demanded by the guard). With --series_n > 0
#   the series part runs only AFTER verdict.csv and the primary outputs are written, inside tryCatch; a failure there is written
#   to SUPPLEMENTARY_SERIES_FAILED.txt and cannot block or change the primary verdict.

source("analysis/dream5_score_lib.R")

## 1. every option and switch is parsed and validated before any guard, net or truth file is touched
o <- d5_truth_options(commandArgs(trailingOnly = TRUE)); opt <- o$opt; sw <- o$sw; net <- opt$net

## 2. guards: nothing below runs, and the truth file is not opened, unless they pass
d5_design_guard()
do.call(d5_truth_guard, d5_truth_guard_args(opt$refit_n, opt$series_n))

nd <- d5_read_net(net)
universe <- d5_universe_from_files(opt, nd)
if (!file.exists(opt$truth_edges)) stop("truth edge file not found: ", opt$truth_edges)
tab <- utils::read.delim(opt$truth_edges, header = FALSE, stringsAsFactors = FALSE)
edges <- d5_edges_from_ids(tab, nd)
message(sprintf("truth: %d edges kept; dropped: %d with IDs not in the gene or TF lists, %d label-0 rows, %d self-edges, %d repeated rows",
                nrow(edges), attr(edges, "n_dropped"), attr(edges, "n_label0"), attr(edges, "n_self"), attr(edges, "n_dup")))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(d5_structure_counts(edges, nd), file.path(opt$out_dir, "structure_counts.csv"), row.names = FALSE)
if (opt$structure_only == 1L) { print(utils::read.csv(file.path(opt$out_dir, "structure_counts.csv"))); quit(save = "no") }

refit_manifest <- file.path(d5_refit_dir(4L), "REFIT_MANIFEST.tsv")
series_manifest <- file.path(d5_refit_dir(4L, unit = "series"), "REFIT_SERIES_MANIFEST.tsv")
scorer <- if (net == 4L) function(b) d5_refit_scores(nd, b, opt$workers) else NULL
verify <- if (net == 4L) function(v, b) d5_refit_verify(v, b, 4L, refit_manifest) else NULL
## 3. primary pipeline (no series part) and primary outputs first
res <- d5_truth_pipeline(nd, edges, sw, nrep = opt$nrep, refit_reps = if (net == 4L) seq_len(opt$refit_n) else integer(0),
                         refit_scorer = scorer, refit_verify = verify, universe = universe, cap = opt$cap)
d5_write_results(res, opt$out_dir)
print(utils::read.csv(file.path(opt$out_dir, "verdict.csv")))

## 4. supplementary series-level jackknife: after the primary outputs, inside tryCatch
if (net == 4L && opt$series_n > 0L) {
  sr <- d5_run_series_supplementary(res, opt$out_dir, nd, edges, sw, seq_len(opt$series_n),
          series_scorer = function(b) d5_refit_scores(nd, b, opt$workers, unit = "series"),
          series_verify = function(v, b) d5_refit_verify(v, b, 4L, series_manifest), universe = universe, cap = opt$cap,
          say = function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))))
  if (sr$ok) { cat("\nSUPPLEMENTARY (series-level jackknife; descriptive; not used in verdict.csv):\n")
    print(sr$sup$series_line[, c("net", "endpoint", "contrast", "n_rep", "SD_series", "SE_B_group")], row.names = FALSE) }
}
