#!/usr/bin/env Rscript
# GRNBoost2 diagnostic (review 2026-10-22): our GRNBoost2 ranks last on bulk
# E. coli and on BEELINE, while Pratapa et al. 2020 report it near GENIE3.
# Audit finding: tools/grnboost2.py reproduces arboreto 0.1.6 faithfully, but
# arboreto's early-stopping monitor (core.EarlyStopMonitor: stop when the mean
# of oob_improvement_ over the last 25 stages is < 0) assumes scikit-learn < 1.3,
# where oob_improvement_[i] is the OOB-loss drop of stage i on ONE OOB sample.
# Since 1.3 it is oob_scores_[i-1] - oob_scores_[i] on two different random 10%
# OOB subsets, so the window mean is sampling noise and most targets stop at the
# earliest allowed stage (25-30 trees instead of a few hundred). Variants
# (tools/grnboost2_variants.py), each with 3 seeds:
#   current     tools/grnboost2.py as used for the reported results
#   legacy_es   arboreto profile with the pre-1.3 (same-sample) OOB improvement
#   no_es       arboreto profile without early stopping (5000 stages)
#   dask        arboreto.algo.grnboost2 itself, local Client(processes=False)
# References in the same table: |Pearson r|; the GRNBoost2 and GENIE3 scores
# already reported (production seed). Scoring is beeline_baselines.R's
# tie-aware score_truth (its definitions evaluated here, not copied): directed
# TF -> target weights, and also sym-max (M3D's GENIE3 is stored sym-max).
# Truths: BEELINE primary (cell-type ChIP, non-specific ChIP, STRING, mESC
# LOF/GOF); M3D: RegulonDB TF -> gene, strong + confirmed evidence (as baselines.R).
# Usage (repo root):
#   Rscript analysis/grnboost2_check.R [--cases hESC_1000,mHSC-E_500,m3d]
#       [--variants current,legacy_es,no_es,dask] [--seeds 1,2,3] [--workers 8]
#       [--out results/grnboost2_check] [--smoke 0] [--m3d_no_es_seeds 1]
# no_es on M3D costs about 200 CPU-hours per seed (4,297 targets x 5,000 stages),
# so it runs only for the first --m3d_no_es_seeds seeds.
# --smoke n: each case reduced to the n genes with the most primary-truth edges
#   (for testing; writes its own expression.tsv under --out).
# Edge files are cached under <out>/<case>/<variant>_s<seed>.tsv.

args0 <- commandArgs(trailingOnly = TRUE)
gopt <- list(cases = "hESC_1000,mHSC-E_500,m3d", variants = "current,legacy_es,no_es,dask", seeds = "1,2,3",
             workers = 8L, out = "results/grnboost2_check", smoke = 0L, m3d_no_es_seeds = 1L, python = "~/venvs/clr-baselines/bin/python")
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i]); if (!k %in% names(gopt)) stop("unknown argument: ", args0[i])
  gopt[[k]] <- if (k %in% c("workers", "smoke", "m3d_no_es_seeds")) as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
dir.create(gopt$out, recursive = TRUE, showWarnings = FALSE)

# beeline_baselines.R's setup (which evaluates beeline.R's) up to and including score_truth
bsrc <- parse("analysis/beeline_baselines.R", keep.source = TRUE)
stop_at <- which(vapply(bsrc, function(e) grepl("^say\\(\"beeline_baselines: stage", deparse(e)[1]), TRUE))[1]
args0 <- c("--out", gopt$out, "--threads", as.character(gopt$workers))
for (e in as.list(bsrc)[seq_len(stop_at - 1L)]) {
  if (identical(deparse(e)[1], "args0 <- commandArgs(trailingOnly = TRUE)")) next
  eval(e, globalenv())
}
logf <- file.path(gopt$out, "grnboost2_check.log")
say("grnboost2_check: ", paste(names(gopt), unlist(lapply(gopt, format)), sep = "=", collapse = " "))
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")

# M3D truth: RegulonDB TF -> target, evidence S/C (as analysis/baselines.R)
m3d_truth <- function() {
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
  tfr <- rtab[rtab$reg_class == "TF" & !is.na(rtab$bnumber), ]
  tp <- strsplit(.classic_tf_to_genes(tfr$regulator), ";", fixed = TRUE)
  net <- data.frame(Gene1 = to_bn(unlist(tp)), Gene2 = rep(tfr$bnumber, lengths(tp)),
                    confidence = rep(tfr$confidence, lengths(tp)), stringsAsFactors = FALSE)
  net <- net[!is.na(net$Gene1) & net$confidence %in% c("S", "C"), ]
  list(RegulonDB_SC = net)
}
case_info <- function(cs) {
  if (cs == "m3d") return(list(src = "results/baselines_m3d", gts = m3d_truth(),
                               existing = c(grnboost2 = "results/baselines_m3d/grnboost2_edges.tsv",
                                            genie3_symmax = "results/baselines_m3d/genie3_scores.rds")))
  ds <- sub("_[0-9]+$", "", cs); d <- file.path("results/beeline_baselines", cs)
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  list(src = d, gts = gts[intersect(names(gts), primary)], existing = c(grnboost2 = file.path(d, "grnboost2_edges.tsv")))
}
edges_mat <- function(f, genes) { e <- utils::read.delim(f, stringsAsFactors = FALSE, na.strings = character(0))
  M <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  ia <- match(e$TF, genes); ib <- match(e$target, genes); ok <- !is.na(ia) & !is.na(ib) & ia != ib
  M[cbind(ia[ok], ib[ok])] <- e$importance[ok]; M }
score_all <- function(S, genes, gts, lab) {
  S[!is.finite(S)] <- 0; diag(S) <- 0; out <- list()
  for (orient in c("directed", "symmax")) { SS <- if (orient == "symmax") pmax(S, t(S)) else S
    for (tn in names(gts)) { r <- score_truth(SS, genes, NULL, gts[[tn]]); if (is.null(r)) next
      out[[length(out) + 1L]] <- data.frame(lab, truth = tn, orientation = orient, t(r), check.names = FALSE) } }
  do.call(rbind, out)
}

rows <- list(); seeds <- as.integer(strsplit(gopt$seeds, ",")[[1]]); variants <- strsplit(gopt$variants, ",")[[1]]
for (cs in strsplit(gopt$cases, ",")[[1]]) {
  ci <- case_info(cs); cd <- file.path(gopt$out, cs); dir.create(cd, showWarnings = FALSE)
  inp <- file.path(ci$src, "expression.tsv")
  X <- as.matrix(utils::read.delim(inp, row.names = 1, check.names = FALSE, na.strings = character(0)))
  if (gopt$smoke > 0L) {        # genes with the most truth edges (both ends counted)
    deg <- table(toupper(unlist(lapply(ci$gts, function(g) c(g$Gene1, g$Gene2)))))
    d <- as.numeric(deg[toupper(rownames(X))]); d[is.na(d)] <- 0
    X <- X[sort(order(-d)[seq_len(min(gopt$smoke, nrow(X)))]), , drop = FALSE]
    inp <- file.path(cd, "expression.tsv")
    utils::write.table(data.frame(gene = rownames(X), X, check.names = FALSE), inp, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  genes <- rownames(X); say(sprintf("%s: %d genes x %d samples", cs, nrow(X), ncol(X)))
  r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0
  rows[[length(rows) + 1L]] <- score_all(r, genes, ci$gts, data.frame(case = cs, variant = "abs_pearson", seed = NA,
                                         n_nonzero = NA, minutes = NA)); rm(r)
  if (gopt$smoke == 0L) for (m in names(ci$existing)) { f <- ci$existing[[m]]; if (!file.exists(f)) next
    S <- if (grepl("\\.rds$", f)) readRDS(f)[genes, genes] else edges_mat(f, genes)
    rows[[length(rows) + 1L]] <- score_all(S, genes, ci$gts, data.frame(case = cs, variant = paste0(m, "_reported"),
                                           seed = NA, n_nonzero = sum(S > 0), minutes = NA)) }
  for (v in variants) for (sd in seeds) {
    if (cs == "m3d" && v == "no_es" && match(sd, seeds) > gopt$m3d_no_es_seeds) next
    of <- file.path(cd, sprintf("%s_s%d.tsv", v, sd)); lf <- sub("\\.tsv$", ".log", of); mins <- NA
    if (!file.exists(of)) {
      a <- if (v == "current") c("tools/grnboost2.py", inp, paste0(of, ".part"), "--workers", gopt$workers, "--seed", sd) else
        c("tools/grnboost2_variants.py", inp, paste0(of, ".part"), "--mode", v, "--workers", gopt$workers, "--seed", sd)
      t0 <- Sys.time(); st <- system2(path.expand(gopt$python), a, stdout = lf, stderr = lf)
      mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      if (st != 0 || !file.exists(paste0(of, ".part"))) { say(sprintf("%s %s seed %d FAILED (see %s)", cs, v, sd, lf)); next }
      file.rename(paste0(of, ".part"), of)
    }
    lg <- readLines(lf); say(sprintf("%s %s seed %d: %s", cs, v, sd, utils::tail(lg[grepl("grnboost2", lg)], 1)))
    S <- edges_mat(of, genes)
    rows[[length(rows) + 1L]] <- score_all(S, genes, ci$gts, data.frame(case = cs, variant = v, seed = sd,
                                           n_nonzero = sum(S > 0), minutes = mins))
  }
  res <- do.call(rbind, rows); utils::write.csv(res, file.path(gopt$out, "metrics.csv"), row.names = FALSE)
  x <- res[res$case == cs & res$orientation == "directed", ]
  for (v in unique(x$variant)) { y <- x[x$variant == v, ]
    say(sprintf("  %-10s %-20s mean over truths x seeds: AUPRC ratio %.3f  EPR %.3f  (n = %d)", cs, v,
                mean(y$auprc_ratio), mean(y$epr), nrow(y))) }
}
say("done")
