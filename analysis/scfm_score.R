#!/usr/bin/env Rscript
# Preregistered scoring of the learned-attention comparator on the 14 BEELINE cases (docs/preregistration_scfm.md,
# change log entries E1-E13). Stage 1 of 2: every metric, one long table. Stage 2 (analysis/scfm_score_stats.R) turns the table
# into the preregistered interval, decision, secondaries and descriptive tables.
#
# The scoring definitions are EVALUATED from the repository, not copied: analysis/beeline.R (its setup, helpers and
# truths(), everything before its main loop), and score_truth / case_X / edges_to of analysis/beeline_baselines.R (tie-aware
# early precision). Directed TF-node AUPRC ratio and EPR; self-pairs 0; non-finite values 0.
#
# Inputs (no truth file is read before the stage that scores):
#   ~/scfm_models/step4_scoring/<source>/<case>.<readout>.f32   written by analysis/scfm_score_export.py after verifying the
#                                                                SHA-256 manifest of the step-3 freeze (E13)
#   data/beeline/BEELINE-data/.../ExpressionData.csv            (case_X: zero-range genes dropped, as beeline_baselines.R)
#   results/beeline_baselines/<case>/{genie3_scores.rds,grnboost2_edges.tsv,pidc_edges.tsv}   descriptive comparators
#   truths: truths(ds) of analysis/beeline.R (real), or --truth synthetic --synth_dir DIR (smoke test only)
#
# Orientation (primary): score(TF -> g) = A[g, TF], i.e. S = t(A). Other orientations: A and pmax(A, t(A)).
# Gene sets (E12): each contrast is scored on its own G-intersection; the comparators |r|, CLR, CLR-NS are recomputed on it
# from the BEELINE expression (CLR-NS also on the full gene set, reported next to it). Gene-set keys:
#   full = all case genes; A = scPRINT-covered genes; B = arm-B genes (key A when identical); S = scGPT genes; AS = A and S.
#
# Usage (repo root):
#   Rscript analysis/scfm_score.R --out results/scfm_score [--truth real|synthetic --synth_dir DIR]
#       [--cases hESC_500,...] [--threads 8] [--exports ~/scfm_models/step4_scoring] [--att 1]

args0 <- commandArgs(trailingOnly = TRUE)
sopt <- list(out = "results/scfm_score", truth = "real", synth_dir = NA_character_, cases = "all", threads = 8L,
             exports = "~/scfm_models/step4_scoring", att = 1L, seed = 20261021L, resume = 0L)
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(sopt))
  sopt[[k]] <- if (k %in% c("threads", "att", "seed", "resume")) as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
sopt$exports <- path.expand(sopt$exports)

# ---- beeline.R setup and helpers (everything before its main loop), as analysis/beeline_baselines.R does ----------------
src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]
args <- c("--out", sopt$out, "--threads", as.character(sopt$threads))
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  if (identical(deparse(e)[1], "args <- commandArgs(trailingOnly = TRUE)")) next
  eval(e, globalenv())
}
# score_truth (tie-aware EPR), case_X and edges_to from beeline_baselines.R, evaluated from the file
bb <- parse("analysis/beeline_baselines.R")
for (e in bb) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] %in% c("score_truth", "edges_to", "case_X")) eval(e, globalenv())
stopifnot(exists("score_truth"), exists("edges_to"), exists("case_X"))     # case_X: zero-range genes dropped (--context none)
if (sopt$truth == "synthetic") {
  stopifnot(!is.na(sopt$synth_dir))
  truths <- function(ds) { d <- file.path(sopt$synth_dir, ds)
    nm <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", if (ds == "mESC") "LOFGOF", "STRING_regulatory", "STRING_curated")
    stats::setNames(file.path(d, paste0(nm, ".csv")), nm) }
}
base_dir <- "results/beeline_baselines"
say(sprintf("scfm_score: truth=%s out=%s threads=%d att=%d", sopt$truth, sopt$out, sopt$threads, sopt$att))

CASES <- c(t(outer(c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L"), c(500, 1000), paste, sep = "_")))
if (sopt$cases != "all") CASES <- strsplit(sopt$cases, ",")[[1]]

# ---- flat matrix reader -------------------------------------------------------------------------------------------------
rd_mat <- function(src, case, readout) {
  f <- file.path(sopt$exports, src, sprintf("%s.%s.f32", case, readout)); if (!file.exists(f)) return(NULL)
  g <- readLines(file.path(sopt$exports, src, sprintf("%s.genes.txt", case))); G <- length(g)
  m <- matrix(readBin(f, "numeric", n = G * G, size = 4L), G, G); dimnames(m) <- list(g, g)
  stopifnot(all(is.finite(m))); m
}
rd_genes <- function(src, case) { f <- file.path(sopt$exports, src, sprintf("%s.genes.txt", case)); if (file.exists(f)) readLines(f) else NULL }
orient <- function(A, o) switch(o, T = t(A), A = A, max = pmax(A, t(A)), sym = A)
clean <- function(s) { s[!is.finite(s)] <- 0; diag(s) <- 0; s }

# ---- scoring a matrix against all truths -------------------------------------------------------------------------------
rows <- list()
score_all <- function(S, genes, gts, meta) {
  S <- clean(S)
  for (tn in names(gts)) {
    r <- score_truth(S, genes, NULL, gts[[tn]]); if (is.null(r)) next
    rows[[length(rows) + 1L]] <<- data.frame(meta, truth = tn, t(r), row.names = NULL)
  }
}

t_start <- Sys.time(); meta_rows <- list()
partial_dir <- file.path(sopt$out, "partial"); dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
for (ci in seq_along(CASES)) {
  case <- CASES[ci]; ds <- sub("_[0-9]+$", "", case); N <- as.integer(sub("^.*_", "", case))
  pf <- file.path(partial_dir, paste0(case, ".rds"))
  if (sopt$resume == 1L && file.exists(pf)) { rows[[length(rows) + 1L]] <- readRDS(pf); say(case, ": resumed from ", pf); next }
  rows_before <- length(rows); t0 <- Sys.time()
  X <- case_X(ds, N); gfull <- rownames(X)
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  # sources of this case and their gene vectors
  gA <- rd_genes("scprint_medium_A", case); stopifnot(!is.null(gA), all(gA %in% gfull))
  gB <- rd_genes("scprint_medium_B", case); gS <- rd_genes("scgpt_A", case)
  if (!is.null(gB)) stopifnot(all(gB %in% gfull)); if (!is.null(gS)) stopifnot(all(gS %in% gfull))
  gsets <- list(full = gfull, A = gA)
  key_of <- function(g) { for (k in names(gsets)) if (identical(gsets[[k]], g)) return(k); NA_character_ }
  keyB <- NA; keyS <- NA; keyAS <- NA
  if (!is.null(gB)) { keyB <- key_of(gB); if (is.na(keyB)) { gsets$B <- gB; keyB <- "B" } }
  if (!is.null(gS)) { keyS <- key_of(gS); if (is.na(keyS)) { gsets$S <- gS; keyS <- "S" }
    gAS <- gA[gA %in% gS]; keyAS <- key_of(gAS); if (is.na(keyAS)) { gsets$AS <- gAS; keyAS <- "AS" } }
  say(sprintf("%s: %d cells, gene sets %s", case, ncol(X),
              paste(sprintf("%s=%d", names(gsets), lengths(gsets)), collapse = " ")))
  meta_rows[[length(meta_rows) + 1L]] <- data.frame(dataset = ds, N = N, n_cells = ncol(X), n_full = length(gfull), n_A = length(gA),
    n_B = if (is.null(gB)) NA_integer_ else length(gB), n_S = if (is.null(gS)) NA_integer_ else length(gS),
    n_AS = if (is.null(gS)) NA_integer_ else length(gA[gA %in% gS]), key_A = "A", key_B = keyB, key_S = keyS, key_AS = keyAS,
    n_tfs_full = sum(toupper(gfull) %in% select_genes(ds, N)$tfs), n_tfs_A = sum(toupper(gA) %in% select_genes(ds, N)$tfs))
  base_meta <- function(key, source, readout, orient) data.frame(dataset = ds, N = N, gset = key, G = length(gsets[[key]]),
                                                                  source = source, readout = readout, orient = orient)
  # descriptive baselines fitted on the full gene set (cached outputs of results/beeline_baselines), restricted to a gene set
  bd <- file.path(base_dir, case); full_base <- list()
  if (file.exists(file.path(bd, "genie3_scores.rds"))) { W <- readRDS(file.path(bd, "genie3_scores.rds")); full_base$genie3 <- W[gfull, gfull] }
  if (file.exists(file.path(bd, "grnboost2_edges.tsv"))) { e <- utils::read.delim(file.path(bd, "grnboost2_edges.tsv"), stringsAsFactors = FALSE)
    full_base$grnboost2 <- edges_to(e$TF, e$target, e$importance, gfull, FALSE) }
  if (file.exists(file.path(bd, "pidc_edges.tsv"))) { e <- utils::read.delim(file.path(bd, "pidc_edges.tsv"), header = FALSE, stringsAsFactors = FALSE)
    full_base$pidc <- edges_to(e[[1]], e[[2]], e[[3]], gfull, TRUE) }
  # comparators on every gene set: |r|, CLR, CLR-NS recomputed on it (the CLR background is the gene set itself)
  for (key in names(gsets)) {
    g <- gsets[[key]]; Xg <- X[g, , drop = FALSE]
    S <- abs(stats::cor(t(Xg))); score_all(S, g, gts, base_meta(key, "cmp", "abs_pearson", "sym"))
    f <- clr_fit(Xg); S <- f$clr_scores; f$release(FALSE); score_all(S, g, gts, base_meta(key, "cmp", "clr", "sym"))
    Cns <- ns_clr(Xg); score_all(Cns, g, gts, base_meta(key, "cmp", "clr_ns", "sym"))
    if (sopt$att == 1L && key == "A") {          # CLR-NS soft10 attention at its held-out depth (descriptive, preregistered item 5)
      set.seed(sopt$seed + ci)
      t_ns <- heldout_depth(Xg, function(Y) clr:::.softmax_rows(ns_clr(Y), keff = 10, cap = 50L), paste(case, "ns_soft10 on G-A"),
                            rownames(Xg), target_fn = ns_clr)
      Sa <- att_mass(clr:::.softmax_rows(Cns, keff = 10, cap = 50L), t_ns, cols = seq_len(nrow(Xg)))[[1]]
      m <- base_meta(key, "cmp", "clr_ns_att", "sym"); m$depth <- t_ns
      score_all(Sa, g, gts, m); rm(Sa)
    }
    for (b in names(full_base)) score_all(full_base[[b]][g, g], g, gts, base_meta(key, "cmp", b, if (b %in% c("pidc")) "sym" else "dir"))
    rm(Cns, S); invisible(gc(verbose = FALSE))
  }
  rm(full_base)
  # scPRINT sources
  add <- function(srcname, readout, key, orients, file_src = srcname, file_readout = readout, transform = identity) {
    A <- rd_mat(file_src, case, file_readout); if (is.null(A)) return(invisible())
    A <- transform(A)
    for (o in orients) score_all(orient(A, o), gsets[[key]], gts, base_meta(key, srcname, readout, o))
  }
  O3 <- c("T", "A", "max")
  for (rr in c("R0", "R0p", "R2")) {
    add("scprint_medium_A", rr, "A", O3)
    if (!is.na(keyB)) add("scprint_medium_B", rr, keyB, O3)
  }
  add("r1_A", "R1", "A", O3)
  if (!is.na(keyB)) add("r1_B", "R1", keyB, O3)
  add("scprint_large_A", "R0", "A", "T"); add("scprint_small_A", "R0", "A", "T")
  for (sd in 1001:1003) add(paste0("ctrl_perm_s", sd), "R0", "A", "T")
  for (sd in 1:3) { add(paste0("ctrl_rand_esm_s", sd), "R0", "A", "T"); add(paste0("ctrl_rand_emb_s", sd), "R0", "A", "T") }
  add("esm2_cos", "cos", "A", "sym")
  # CLR applied to A (hybrid): CLR calibration of the symmetrised attention (A + t(A)) / 2, arm A, R0
  A <- rd_mat("scprint_medium_A", case, "R0"); Sh <- (A + t(A)) / 2; diag(Sh) <- 0
  score_all(clr_calibrate(Sh, method = "normal", combine = "stouffer"), gsets$A, gts, base_meta("A", "scprint_medium_A", "R0_clr_hybrid", "sym"))
  rm(A, Sh)
  # scGPT (human cases): attention (S = t(A)), token-embedding cosine, hidden-state cosine; exploratory: LayerNorm'd token cosine
  if (!is.na(keyS)) {
    add("scgpt_A", "attn_last_mean", keyS, "T")
    add("scgpt_A", "tok_cos", keyS, "sym"); add("scgpt_A", "hid_cos", keyS, "sym"); add("scgpt_A", "tok_cos_ln", keyS, "sym")
    # scPRINT R0 and ESM2 cosine on the three-way intersection (scPRINT vs scGPT contrast, E12)
    if (keyAS != "A") { A <- rd_mat("scprint_medium_A", case, "R0"); ii <- match(gsets[[keyAS]], gA)
      score_all(t(A[ii, ii]), gsets[[keyAS]], gts, base_meta(keyAS, "scprint_medium_A", "R0", "T")); rm(A) }
    if (keyAS != keyS) for (rr in c("attn_last_mean", "tok_cos", "hid_cos", "tok_cos_ln")) {
      A <- rd_mat("scgpt_A", case, rr); jj <- match(gsets[[keyAS]], gS)
      score_all(orient(A[jj, jj], if (rr == "attn_last_mean") "T" else "sym"), gsets[[keyAS]], gts,
                base_meta(keyAS, "scgpt_A", rr, if (rr == "attn_last_mean") "T" else "sym"))
    }
    # comparators on the three-way set were computed with the gene-set loop above (key AS)
  }
  newrows <- rows[(rows_before + 1L):length(rows)]
  saveRDS(newrows, pf)
  say(sprintf("%s: done, %d score rows, %.1f min", case, sum(vapply(newrows, nrow, 1L)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(X, gts); invisible(gc(verbose = FALSE))
}
res <- do.call(rbind, lapply(rows, function(r) { if (!"depth" %in% names(r)) r$depth <- NA_integer_; r }))
utils::write.csv(res, file.path(sopt$out, "metrics_long.csv"), row.names = FALSE)
if (length(meta_rows)) utils::write.csv(do.call(rbind, meta_rows), file.path(sopt$out, "cases_meta.csv"), row.names = FALSE)
say(sprintf("wrote %d rows to %s in %.1f min", nrow(res), file.path(sopt$out, "metrics_long.csv"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
writeLines(c(capture.output(sessionInfo()), "", paste("seed", sopt$seed)), file.path(sopt$out, "sessionInfo.txt"))
say("done")
