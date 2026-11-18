#!/usr/bin/env Rscript
# GENIE3 for one arm of the frozen tcaCLR plan (section 3): analysis/baselines.R settings exactly (GENIE3 R package, RF, K = sqrt, 1,000 trees,
# every gene a regulator and a target, set.seed(20261016); the doRNG path (nCores >= 2) makes the result independent of the number of workers).
# Same as analysis/scoreboard/sb_cx_rpe1_genie3.R with the work directory as an argument. Saves the directed G x G weight matrix, uncompressed.
#   Rscript analysis/tcaclr_confirm/run_genie3.R <work dir holding expression.tsv> [workers]
args <- commandArgs(trailingOnly = TRUE); w <- args[1]; workers <- if (length(args) > 1L) as.integer(args[2]) else 32L
suppressPackageStartupMessages(library(GENIE3))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
ex <- utils::read.delim(file.path(w, "expression.tsv"), check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; rm(ex)
genes <- rownames(X); stopifnot(!anyDuplicated(genes)); say("GENIE3: ", nrow(X), " genes x ", ncol(X), " profiles, ", workers, " workers")
invisible(gc()); set.seed(20261016L); t0 <- Sys.time()
W <- GENIE3(X, regulators = NULL, targets = NULL, treeMethod = "RF", K = "sqrt", nTrees = 1000L, nCores = workers, verbose = FALSE)
W <- W[genes, genes]; say(sprintf("GENIE3 computed in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
tmp <- file.path(w, "genie3_W.rds.part"); saveRDS(W, tmp, compress = FALSE); file.rename(tmp, file.path(w, "genie3_W.rds")); say("saved")
