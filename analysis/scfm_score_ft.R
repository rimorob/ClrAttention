#!/usr/bin/env Rscript
# Descriptive fine-tuning arm (E10/E17/E18). Stage 1 of 2: score the FINE-TUNED attention networks on the real truths with the same
# scoring function, gene sets, orientations and truths as analysis/scfm_score.R. Stage 2: analysis/scfm_score_ft_stats.R.
# The CLR-NS and zero-shot rows are NOT recomputed: they are read from results/scfm_score/metrics_long.csv (E16) in stage 2.
# The scoring definitions are EVALUATED from the repository, as in scfm_score.R: analysis/beeline.R (setup, helpers, truths()), and
# score_truth / edges_to / case_X of analysis/beeline_baselines.R.
#
# Inputs (no truth file is read before the scoring loop):
#   <exports>/<src>/<case>.<readout>.f32  from analysis/scfm_score_export_ft.py (after verifying the fine-tuning manifest, E17 freeze)
#   results/scfm_score/cases_meta.csv     gene-set keys of the E16 scoring (key_S for scGPT)
# Scored (all S = t(A) unless stated): scPRINT fine-tuned R0 (orientations T, A, max), R0p (T), R2 (T) on gene set A;
#   scGPT fine-tuned attn_last_mean (T), tok_cos / tok_cos_ln / hid_cos (sym) on the scGPT gene set (hESC, hHep).
# Usage (repo root): Rscript analysis/scfm_score_ft.R --out results/scfm_score_ft [--truth real|synthetic --synth_dir DIR]
#      [--exports ~/scfm_models/step5_ft/scoring] [--src_scprint scprint_ft_A] [--src_scgpt scgpt_ft_A]
# (--src_* only relocate the source directory for the synthetic-truth smoke test; the "source" label in the table is always scprint_ft_A / scgpt_ft_A)
args0 <- commandArgs(trailingOnly = TRUE)
sopt <- list(out = "results/scfm_score_ft", truth = "real", synth_dir = NA_character_, cases = "all", threads = 4L,
             exports = "~/scfm_models/step5_ft/scoring", src_scprint = "scprint_ft_A", src_scgpt = "scgpt_ft_A", meta = "results/scfm_score/cases_meta.csv")
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(sopt))
  sopt[[k]] <- if (k == "threads") as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
sopt$exports <- path.expand(sopt$exports)
src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]
args <- c("--out", sopt$out, "--threads", as.character(sopt$threads))
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  if (identical(deparse(e)[1], "args <- commandArgs(trailingOnly = TRUE)")) next
  eval(e, globalenv())
}
bb <- parse("analysis/beeline_baselines.R")
for (e in bb) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] %in% c("score_truth", "edges_to", "case_X")) eval(e, globalenv())
stopifnot(exists("score_truth"), exists("edges_to"), exists("case_X"))
if (sopt$truth == "synthetic") {
  stopifnot(!is.na(sopt$synth_dir))
  truths <- function(ds) { d <- file.path(sopt$synth_dir, ds)
    nm <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", if (ds == "mESC") "LOFGOF", "STRING_regulatory", "STRING_curated")
    stats::setNames(file.path(d, paste0(nm, ".csv")), nm) }
}
dir.create(sopt$out, recursive = TRUE, showWarnings = FALSE)
say(sprintf("scfm_score_ft: truth=%s out=%s exports=%s", sopt$truth, sopt$out, sopt$exports))
CM <- utils::read.csv(sopt$meta, stringsAsFactors = FALSE)
CASES <- c(t(outer(c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L"), c(500, 1000), paste, sep = "_")))
if (sopt$cases != "all") CASES <- strsplit(sopt$cases, ",")[[1]]
rd_mat <- function(srcdir, case, readout) {
  f <- file.path(sopt$exports, srcdir, sprintf("%s.%s.f32", case, readout)); if (!file.exists(f)) return(NULL)
  g <- readLines(file.path(sopt$exports, srcdir, sprintf("%s.genes.txt", case))); G <- length(g)
  m <- matrix(readBin(f, "numeric", n = G * G, size = 4L), G, G); dimnames(m) <- list(g, g); stopifnot(all(is.finite(m))); m
}
rd_genes <- function(srcdir, case) { f <- file.path(sopt$exports, srcdir, sprintf("%s.genes.txt", case)); if (file.exists(f)) readLines(f) else NULL }
orient <- function(A, o) switch(o, T = t(A), A = A, max = pmax(A, t(A)), sym = A)
clean <- function(s) { s[!is.finite(s)] <- 0; diag(s) <- 0; s }
rows <- list()
score_all <- function(S, genes, gts, meta) {
  S <- clean(S)
  for (tn in names(gts)) { r <- score_truth(S, genes, NULL, gts[[tn]]); if (is.null(r)) next
    rows[[length(rows) + 1L]] <<- data.frame(meta, truth = tn, t(r), row.names = NULL) }
}
for (case in CASES) {
  ds <- sub("_[0-9]+$", "", case); N <- as.integer(sub("^.*_", "", case)); t0 <- Sys.time()
  X <- case_X(ds, N); gfull <- rownames(X)
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  cm <- CM[CM$dataset == ds & CM$N == N, ]; stopifnot(nrow(cm) == 1L)
  gA <- rd_genes(sopt$src_scprint, case); stopifnot(!is.null(gA), all(gA %in% gfull), length(gA) == cm$n_A)
  mk <- function(key, g, srcname, readout, o) data.frame(dataset = ds, N = N, gset = key, G = length(g), source = srcname, readout = readout, orient = o)
  for (rr in c("R0", "R0p", "R2")) { A <- rd_mat(sopt$src_scprint, case, rr); stopifnot(!is.null(A))
    for (o in if (rr == "R0") c("T", "A", "max") else "T") score_all(orient(A, o), gA, gts, mk("A", gA, "scprint_ft_A", rr, o)) }
  gS <- rd_genes(sopt$src_scgpt, case)
  if (!is.null(gS)) { stopifnot(length(gS) == cm$n_S, all(gS %in% gfull))
    for (rr in c("attn_last_mean", "tok_cos", "hid_cos", "tok_cos_ln")) { A <- rd_mat(sopt$src_scgpt, case, rr); stopifnot(!is.null(A))
      o <- if (rr == "attn_last_mean") "T" else "sym"; score_all(orient(A, o), gS, gts, mk(cm$key_S, gS, "scgpt_ft_A", rr, o)) } }
  say(sprintf("%s: done, %.1f min", case, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(X, gts); invisible(gc(verbose = FALSE))
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(sopt$out, "metrics_ft.csv"), row.names = FALSE)
say(sprintf("wrote %d rows to %s", nrow(res), file.path(sopt$out, "metrics_ft.csv")))
writeLines(capture.output(sessionInfo()), file.path(sopt$out, "sessionInfo.txt"))
say("done")
