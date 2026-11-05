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
#   --truth_edges: TSV, two columns (regulator gene ID, target gene ID; G1...), no header; the extraction from the
#   official gold-standard files is part of the truth stage and is not done here.
#   --structure_only 1: only the structure counts and the Y-M evaluability gate (no score is read).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(net = "", truth_edges = "", out_dir = "", ties_at_cap = NA_character_, ties_in_list = NA_character_, tftf = NA_character_,
            universe = NA_character_, universe_tf = "", universe_gene = "", cap_order = NA_character_, jackknife_cap = NA_character_,
            structure_only = 0L, workers = 16L, nrep = 200L, cap = 100000L, refit_n = 100L)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]; opt[[key]] <- if (is.integer(opt[[key]])) as.integer(v) else v; i <- i + 2L
}
source("analysis/dream5_score_lib.R")
net <- as.integer(opt$net); stopifnot(net %in% D5$nets, nzchar(opt$truth_edges), nzchar(opt$out_dir))

## guards first: nothing below runs, and the truth file is not opened, unless they pass
d5_design_guard()
refit_manifest <- file.path(d5_refit_dir(4L), "REFIT_MANIFEST.tsv")
d5_truth_guard(refit_manifest = refit_manifest,
               refit_required = expand.grid(net = 4L, rep = seq_len(opt$refit_n), method = D5$refit_methods, stringsAsFactors = FALSE))

nd <- d5_read_net(net)
tab <- utils::read.delim(opt$truth_edges, header = FALSE, stringsAsFactors = FALSE)
edges <- d5_edges_from_ids(tab, nd); message(sprintf("truth: %d edges kept, %d dropped (IDs not in the gene or TF lists)", nrow(edges), attr(edges, "n_dropped")))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(d5_structure_counts(edges, nd), file.path(opt$out_dir, "structure_counts.csv"), row.names = FALSE)
if (opt$structure_only == 1L) { print(utils::read.csv(file.path(opt$out_dir, "structure_counts.csv"))); quit(save = "no") }

sw <- d5_edge_switches(ties_at_cap = opt$ties_at_cap, ties_in_list = opt$ties_in_list, tftf = opt$tftf, universe = opt$universe,
                       cap_order = opt$cap_order, jackknife_cap = opt$jackknife_cap)
universe <- if (nzchar(opt$universe_tf)) list(tf = match(readLines(opt$universe_tf), nd$tf_ids), gene = match(readLines(opt$universe_gene), nd$genes)) else NULL
scorer <- if (net == 4L) function(b) d5_refit_scores(nd, b, opt$workers) else NULL
verify <- if (net == 4L) function(v, b) d5_refit_verify(v, b, 4L, refit_manifest) else NULL
res <- d5_truth_pipeline(nd, edges, sw, nrep = opt$nrep, refit_reps = if (net == 4L) seq_len(opt$refit_n) else integer(0),
                         refit_scorer = scorer, refit_verify = verify, universe = universe, cap = opt$cap)
d5_write_results(res, opt$out_dir)
print(utils::read.csv(file.path(opt$out_dir, "verdict.csv")))
