#!/usr/bin/env Rscript
# SYNTHETIC ground truth for the smoke test of analysis/scfm_score.R and analysis/scfm_score_stats.R. It never reads a real
# truth file: the "truths" are planted from the networks themselves (and are therefore unrelated to any real regulation).
# Written in the BEELINE edge-list format (Gene1 = TF, Gene2 = target) with the real file names of truths(ds), under --out/<ds>/.
#
# Scenarios (--scenario):
#   mixed    celltype_ChIP, STRING_regulatory: planted from CLR-NS-like |r| of normal scores (CLR-NS should win);
#            nonspecific_ChIP, STRING_curated: planted from scPRINT R0 in the orientation S = t(A) (R0 should win and S = t(A)
#            should beat S = A); STRING: half and half; LOFGOF: random
#   scprint  all planted from scPRINT R0 (S = t(A))
#   clrns    all planted from the normal-score correlation
#   null     all random edges (every AUPRC ratio is about 1)
# Usage:  Rscript analysis/scfm_score_synth_truth.R --out DIR [--scenario mixed] [--cases all|...] [--seed 7]
args0 <- commandArgs(trailingOnly = TRUE)
o <- list(out = NA_character_, scenario = "mixed", cases = "all", seed = 7L, exports = "~/scfm_models/step4_scoring")
i <- 1L; while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(o))
  o[[k]] <- if (k == "seed") as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
stopifnot(!is.na(o$out)); o$exports <- path.expand(o$exports)
src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]      # beeline.R setup and helpers (no truths read)
args <- c("--out", file.path(o$out, "_tmp"), "--threads", "2")
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  if (identical(deparse(e)[1], "args <- commandArgs(trailingOnly = TRUE)")) next
  eval(e, globalenv()) }
blom <- function(X) t(apply(X, 1, function(x) qnorm((rank(x, ties.method = "average") - 3 / 8) / (length(x) + 1 / 4))))
set.seed(o$seed)
# one truth set per dataset (as the real truths(ds), shared by both gene sets): planted from the TFs+1000 case
CASES <- paste0(c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L"), "_1000")
if (o$cases != "all") CASES <- strsplit(o$cases, ",")[[1]]
for (case in CASES) {
  ds <- sub("_[0-9]+$", "", case); N <- as.integer(sub("^.*_", "", case))
  ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1, check.names = FALSE))
  g <- intersect(select_genes(ds, N)$genes, rownames(ex)); X <- ex[g, , drop = FALSE]
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]; rm(ex)
  genes <- readLines(file.path(o$exports, "scprint_medium_A", paste0(case, ".genes.txt"))); G <- length(genes)
  A <- matrix(readBin(file.path(o$exports, "scprint_medium_A", paste0(case, ".R0.f32")), "numeric", n = G * G, size = 4L), G, G)
  sR0 <- t(A); sNS <- abs(stats::cor(t(blom(X[genes, ])))); rm(A)
  tfs <- intersect(toupper(utils::read.csv(tf_file(species(ds)), stringsAsFactors = FALSE)$TF), toupper(genes))
  tfi <- sort(sample(match(tfs, toupper(genes)), min(60L, length(tfs))))
  plant <- function(kind, m = 25L, noise = 0.3) {
    e <- lapply(tfi, function(r) {
      s <- switch(kind, R0 = rank(sR0[r, ]) / G, NS = rank(sNS[r, ]) / G, half = 0.5 * rank(sR0[r, ]) / G + 0.5 * rank(sNS[r, ]) / G,
                  random = stats::runif(G))
      if (kind != "random") s <- s + stats::rnorm(G, sd = noise)
      s[r] <- -Inf; t <- order(s, decreasing = TRUE)[seq_len(m)]; data.frame(Gene1 = genes[r], Gene2 = genes[t]) })
    do.call(rbind, e) }
  kinds <- switch(o$scenario,
    mixed = c(celltype_ChIP = "NS", nonspecific_ChIP = "R0", STRING = "half", LOFGOF = "random", STRING_regulatory = "NS", STRING_curated = "R0"),
    scprint = c(celltype_ChIP = "R0", nonspecific_ChIP = "R0", STRING = "R0", LOFGOF = "R0", STRING_regulatory = "R0", STRING_curated = "R0"),
    clrns = c(celltype_ChIP = "NS", nonspecific_ChIP = "NS", STRING = "NS", LOFGOF = "NS", STRING_regulatory = "NS", STRING_curated = "NS"),
    null = c(celltype_ChIP = "random", nonspecific_ChIP = "random", STRING = "random", LOFGOF = "random", STRING_regulatory = "random", STRING_curated = "random"))
  d <- file.path(o$out, ds); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(kinds)) utils::write.csv(plant(kinds[[nm]]), file.path(d, paste0(nm, ".csv")), row.names = FALSE, quote = FALSE)
  cat(sprintf("%s: synthetic truths written (%s) for %d genes, %d planted TFs\n", case, o$scenario, G, length(tfi)))
}
