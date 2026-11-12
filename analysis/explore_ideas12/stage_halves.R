#!/usr/bin/env Rscript
# EXPLORATORY. Step 3a of docs/plans/explore_ideas12_plan.md: 6 profile splits (the robust_z depth-stage splits), MI -> CLR -> K10 operator
# on each half, undirected attention graphs W = (A + A')/2 written as edge lists for leiden_sweep.py; top-20 lists for the half-split detector.
args <- commandArgs(TRUE); threads <- if (length(args)) as.integer(args[1]) else 8L
source("analysis/explore_ideas12/ideas12_lib.R")
out <- file.path(EI_OUT, "k562"); gd <- file.path(out, "graphs"); dir.create(gd, recursive = TRUE, showWarnings = FALSE)
say <- ei_say(file.path(out, "stage_halves.log"))
core <- readRDS(file.path(out, "cache", "k562_core.rds")); genes <- core$genes; G <- length(genes)
raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; stopifnot(identical(rownames(X), genes)); lab <- colnames(X)
write_graph <- function(A, name) { W <- (A + t(A)) / 2; diag(W) <- 0; W[lower.tri(W, diag = TRUE)] <- 0; ix <- which(W > 0, arr.ind = TRUE)
  utils::write.table(data.frame(i = ix[, 1] - 1L, j = ix[, 2] - 1L, w = format(W[ix], digits = 17)), file.path(gd, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  nrow(ix) }
ne <- write_graph(as.matrix(core$A), "full"); say(sprintf("full graph: %d edges", ne))
lab_split <- function() {          # copied from robust_z.R
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab]) }
set.seed(EI_DEPTH_SEED); splits <- lapply(1:6, function(s) lab_split())
top20 <- list(); info <- list()
for (si in 1:6) { h <- splits[[si]]; set.seed(EI_DEPTH_SEED + si); top20[[si]] <- list()
  for (d in 1:2) { t0 <- Sys.time(); Xm <- X[, h == d]; rng <- apply(Xm, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9)
    M <- bspline_mi(Xm, bins = "hg", spline_order = 3, threads = threads, transform = "none")
    C <- rz_scores(M, "V0"); rm(M); A <- rz_operator(C); top20[[si]][[d]] <- ei_topk(C, seq_len(G), 20L); rm(C)
    ne <- write_graph(A, sprintf("s%dh%d", si, d)); rm(A); invisible(gc(verbose = FALSE))
    info[[length(info) + 1L]] <- data.frame(split = si, half = d, n_profiles = ncol(Xm), n_zero_range_noised = sum(rng <= 0), edges = ne)
    say(sprintf("split %d half %d: %d profiles, %d edges (%.0f s)", si, d, ncol(Xm), ne, as.numeric(difftime(Sys.time(), t0, units = "secs")))) } }
saveRDS(list(top20 = top20, splits = splits), file.path(out, "cache", "halves_top20.rds"))
utils::write.csv(do.call(rbind, info), file.path(out, "halves_info.csv"), row.names = FALSE)
say("done")
