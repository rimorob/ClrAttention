#!/usr/bin/env Rscript
# RPE1 GENIE3 for the scoreboard (docs/plans/scoreboard_competitors_plan.md, section 8): analysis/baselines.R settings exactly
# (GENIE3 R package, RF, K = sqrt, 1,000 trees, every gene a regulator and a target, set.seed(20261016); the doRNG path (nCores >= 2) makes the
# result independent of the number of workers). Saves the directed G x G weight matrix (uncompressed) for the scoring job.
#   Rscript analysis/scoreboard/sb_cx_rpe1_genie3.R [workers]
args <- commandArgs(trailingOnly = TRUE); workers <- if (length(args)) as.integer(args[1]) else 32L
suppressPackageStartupMessages(library(GENIE3))
w <- "results/scoreboard/complex_competitors/rpe1_work"
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
ex <- utils::read.delim(file.path(w, "expression.tsv"), check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; rm(ex)
genes <- rownames(X); stopifnot(!anyDuplicated(genes)); say("GENIE3: ", nrow(X), " genes x ", ncol(X), " profiles, ", workers, " workers")
invisible(gc()); set.seed(20261016L); t0 <- Sys.time()
W <- GENIE3(X, regulators = NULL, targets = NULL, treeMethod = "RF", K = "sqrt", nTrees = 1000L, nCores = workers, verbose = FALSE)
W <- W[genes, genes]; say(sprintf("GENIE3 computed in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
tmp <- file.path(w, "genie3_W.rds.part"); saveRDS(W, tmp, compress = FALSE); file.rename(tmp, file.path(w, "genie3_W.rds")); say("saved")
