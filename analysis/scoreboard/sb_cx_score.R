#!/usr/bin/env Rscript
# Competitor cells on the complex-truth rows of the raw scoreboard (docs/plans/scoreboard_competitors_plan.md).
# Gates first (a: saved competitor matrices reproduce existing values; b: CLR on each new row reproduces the existing CLR cell;
# c: |r| through the same code path reproduces existing |r| cells, and the yeast |r| recomputed by d5_m_abs_pearson equals its cache).
# Only if EVERY gate of the parts requested passes are the new cells written. Scores are reused (saved matrices), nothing is fitted here.
#   nice -n 19 Rscript analysis/scoreboard/sb_cx_score.R [--parts m3d,ecomac,k562,rpe1,yeast] [--out results/scoreboard/complex_competitors]
args <- commandArgs(trailingOnly = TRUE)
OPT <- list(parts = "m3d,ecomac,k562,rpe1,yeast", out = "results/scoreboard/complex_competitors")
for (i in seq(1, length(args), by = 2)) { k <- sub("^--", "", args[i]); stopifnot(k %in% names(OPT)); OPT[[k]] <- args[i + 1L] }
parts <- strsplit(OPT$parts, ",")[[1]]
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
source("analysis/scoreboard/sb_cx_lib.R"); sb_import_truth_defs()
dir.create(OPT$out, recursive = TRUE, showWarnings = FALSE)
T0 <- Sys.time()
hwm <- function() { l <- readLines("/proc/self/status"); as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", grep("^VmHWM:", l, value = TRUE))) / 1048576 }
G_ <- list(); P_ <- list()
gate <- function(part, gname, item, value, expected, tol = 1e-9) {
  d <- abs(value - expected); ok <- is.finite(d) && d <= tol
  G_[[length(G_) + 1L]] <<- data.frame(part = part, gate = gname, item = item, value = value, expected = expected, abs_diff = d, tol = tol, pass = ok, stringsAsFactors = FALSE)
  sb_say(sprintf("  gate %s | %s | %s: %.12g vs %.12g (diff %.2e) %s", gname, part, item, value, expected, d, if (ok) "PASS" else "FAIL")); invisible(ok) }
cell <- function(ds, truth, method, aupr, p, source, role = "new")
  P_[[length(P_) + 1L]] <<- data.frame(dataset = ds, truth = truth, method = method, aupr = aupr, base_rate = mean(p$label), enrichment = aupr / mean(p$label),
                                       n_pairs = length(p$label), n_pos = sum(p$label), role = role, source = source, stringsAsFactors = FALSE)
score <- function(S, p) { s <- S[cbind(p$i, p$j)]; s[!is.finite(s)] <- NA; ap_fast(s, p$label) }   # the pair-score extraction of clr_logpt_test.R::consume
absr <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; r }
mx <- function(f, tr, m) { x <- utils::read.csv(f, stringsAsFactors = FALSE); x <- x[x$truth == tr & x$method == m, ]; stopifnot(nrow(x) == 1L); x }
sd_cell <- function(ds, tn, m) { x <- utils::read.csv("docs/whitepapers/figures/scoreboard_data.csv", stringsAsFactors = FALSE)
  x <- x[x$dataset == ds & x$truth == tn & x$metric == "pooled_aupr_over_base" & x$scope == "raw" & x$method == m, ]; stopifnot(nrow(x) == 1L); x$enrichment_over_random }

## ---- E. coli: M3D and EcoMAC (saved matrices of analysis/baselines.R) ---------------------------------------------------
for (ds in intersect(parts, c("m3d", "ecomac"))) {
  sb_say("== ", ds); ld <- load_expr(ds); X <- ld$X; genes <- rownames(X); G <- length(genes)
  fp <- readRDS(file.path("results/clr_logpt_test", ds, "splits", "fingerprint.rds"))
  gate(ds, "universe", "gene-set md5 identical to the committed selections' fingerprint (1 = yes)", as.numeric(identical(fp$genes_md5, digest::digest(genes, algo = "md5"))), 1, 0)
  TR <- build_truths_ecoli(ds, genes, NULL, only = c("complexes", "regulons_SC", "regulons_all")); B <- TR$bench
  bdir <- file.path("results", paste0("baselines_", ds)); bm <- file.path(bdir, "metrics.csv")
  rdsm <- list(GENIE3 = "genie3", PIDC = "pidc", GRNBoost2 = "grnboost2")
  S <- lapply(rdsm, function(m) { s <- readRDS(file.path(bdir, paste0(m, "_scores.rds")))
    gate(ds, "universe", paste0(m, " saved matrix: dimnames identical to the gene universe (1 = yes)"), as.numeric(identical(rownames(s), genes) && identical(colnames(s), genes)), 1, 0); s })
  # (a) saved matrices on the existing regulon truths
  for (tn in c("SC", "all")) { b <- B[[c(SC = "regulons_SC", all = "regulons_all")[[tn]]]]
    gate(ds, "a", sprintf("regulons %s base rate vs baselines metrics.csv", tn), mean(b$pairs$label), mx(bm, tn, "genie3")$base_rate)
    for (nm in names(rdsm)) gate(ds, "a", sprintf("%s regulons %s AUPR vs baselines metrics.csv", nm, tn), score(S[[nm]], b$pairs), mx(bm, tn, rdsm[[nm]])$aupr) }
  # (b) CLR on the complexes row, (c) |r|
  ft <- readRDS(file.path("results/clr_logpt_test/truth", ds, "cache", "fit.rds")); C <- ft$C; rm(ft)
  lp <- file.path("results/clr_logpt_test/truth", ds, "metrics.csv"); cx <- B$complexes
  a_clr <- score(C, cx$pairs); rm(C); gate(ds, "b", "CLR complexes AUPR vs clr_logpt_test metrics.csv", a_clr, mx(lp, "complexes", "clr")$aupr)
  gate(ds, "b", "CLR complexes enrichment vs the scoreboard CLR cell", a_clr / mean(cx$pairs$label), sd_cell(ds, "complexes", "CLR"), 1e-9)
  r <- absr(X); a_r <- score(r, cx$pairs); rm(r); gate(ds, "c", "|r| complexes AUPR vs clr_logpt_test metrics.csv", a_r, mx(lp, "complexes", "abs_r")$aupr)
  gate(ds, "c", "|r| complexes enrichment vs the scoreboard |r| cell", a_r / mean(cx$pairs$label), sd_cell(ds, "complexes", "abs_r"), 1e-9)
  cell(ds, "complexes", "CLR", a_clr, cx$pairs, "clr_logpt_test fit.rds (gate b)", "gate"); cell(ds, "complexes", "abs_r", a_r, cx$pairs, "recomputed (gate c)", "gate")
  for (nm in names(rdsm)) cell(ds, "complexes", nm, score(S[[nm]], cx$pairs), cx$pairs, file.path(bdir, paste0(rdsm[[nm]], "_scores.rds")))
  rm(S, X, ld); invisible(gc()); sb_say("peak RSS so far ", round(hwm(), 2), " GB")
}

## ---- K562 CORUM gate (a) only: saved K562 matrices re-scored -----------------------------------------------------------
if ("k562" %in% parts) {
  sb_say("== k562 (gate a only)"); bdir <- "results/baselines_k562ess_truthgenes_genie3"
  genes <- system2("cut", c("-f1", file.path(bdir, "expression.tsv")), stdout = TRUE)[-1]
  mem <- utils::read.csv(file.path(bdir, "regulon_members.csv"), stringsAsFactors = FALSE); mem <- mem[mem$truth == "CORUM", ]
  m <- lapply(strsplit(mem$members, ";", fixed = TRUE), function(g) sort(match(g, genes))); names(m) <- mem$regulon
  pr <- comembership_pairs(list(members = m, class = stats::setNames(rep("complex", length(m)), names(m))), none_op, genes)
  bm <- file.path(bdir, "metrics.csv")
  gate("k562", "a", "CORUM base rate vs baselines metrics.csv", mean(pr$label), mx(bm, "CORUM", "genie3")$base_rate)
  for (nm in c("genie3", "pidc", "grnboost2")) { s <- readRDS(file.path(bdir, paste0(nm, "_scores.rds")))
    gate("k562", "a", sprintf("%s K562 CORUM AUPR vs baselines metrics.csv", nm), score(s, pr), mx(bm, "CORUM", nm)$aupr); rm(s) }
}

## ---- RPE1: gates b and c only (competitors are not run: 2-hour rule, plan section 5) ----------------------------------------
if ("rpe1" %in% parts) {
  sb_say("== rpe1 (gates b, c)"); ld <- load_expr("rpe1"); X <- ld$X; genes <- rownames(X)
  fp <- readRDS("results/clr_logpt_test/rpe1/splits/fingerprint.rds")
  gate("rpe1", "universe", "gene-set md5 identical to the committed selections' fingerprint (1 = yes)", as.numeric(identical(fp$genes_md5, digest::digest(genes, algo = "md5"))), 1, 0)
  TR <- build_truths_human("rpe1", genes, NULL, NULL, NULL); cx <- TR$bench$CORUM_merged; lp <- "results/clr_logpt_test/truth/rpe1/metrics.csv"
  ft <- readRDS("results/clr_logpt_test/truth/rpe1/cache/fit.rds"); C <- ft$C; rm(ft); invisible(gc())
  a_clr <- score(C, cx$pairs); rm(C); gate("rpe1", "b", "CLR CORUM_merged AUPR vs clr_logpt_test metrics.csv", a_clr, mx(lp, "CORUM_merged", "clr")$aupr)
  gate("rpe1", "b", "CLR CORUM_merged enrichment vs the scoreboard CLR cell", a_clr / mean(cx$pairs$label), sd_cell("rpe1", "CORUM_merged", "CLR"), 1e-9)
  r <- absr(X); a_r <- score(r, cx$pairs); rm(r); gate("rpe1", "c", "|r| CORUM_merged AUPR vs clr_logpt_test metrics.csv", a_r, mx(lp, "CORUM_merged", "abs_r")$aupr)
  gate("rpe1", "c", "|r| CORUM_merged enrichment vs the scoreboard |r| cell", a_r / mean(cx$pairs$label), sd_cell("rpe1", "CORUM_merged", "abs_r"), 1e-9)
  cell("rpe1", "CORUM_merged", "CLR", a_clr, cx$pairs, "clr_logpt_test fit.rds (gate b)", "gate"); cell("rpe1", "CORUM_merged", "abs_r", a_r, cx$pairs, "recomputed (gate c)", "gate")
  # competitors (plan addendum 2, section 8): runs on the 8,745 genes left after dropping the duplicated symbols; every score is embedded in the
  # 8,747-gene universe of the CLR cell, pairs of a dropped row get no score (NA, ranked last by ap_fast)
  ww <- "results/scoreboard/complex_competitors/rpe1_work"
  if (file.exists(file.path(ww, "expression.tsv"))) {
    gu <- system2("cut", c("-f1", file.path(ww, "expression.tsv")), stdout = TRUE)[-1]; Gu <- length(gu); pos <- match(gu, genes); G <- length(genes)
    gate("rpe1", "universe", "competitor genes are unique, all in the CLR universe, and 2 rows (the duplicated symbol) are dropped", as.numeric(!anyDuplicated(gu) && !anyNA(pos) && G - Gu == 2L), 1, 0)
    dropped <- setdiff(seq_len(G), pos); npd <- sum(cx$pairs$label == 1 & (cx$pairs$i %in% dropped | cx$pairs$j %in% dropped))
    sb_say("positive CORUM pairs involving a dropped row: ", npd, " (of ", sum(cx$pairs$label), ")")
    gate("rpe1", "universe", "positive truth pairs involving the dropped duplicated-symbol rows (0 = identical to scoring the 8,745-gene truth)", npd, 0, 0)
    embed <- function(M) { E <- matrix(NA_real_, G, G); E[pos, pos] <- M; E }
    rc <- function(nm, M, src) { a <- score(embed(M), cx$pairs); cell("rpe1", "CORUM_merged", nm, a, cx$pairs, src); sb_say("RPE1 ", nm, ": enrichment ", round(a / mean(cx$pairs$label), 4)) }
    f <- file.path(ww, "genie3_W.rds")
    if (file.exists(f)) { W <- readRDS(f); gate("rpe1", "universe", "GENIE3 matrix dimnames equal the competitor gene list (1 = yes)", as.numeric(identical(rownames(W), gu) && identical(colnames(W), gu)), 1, 0)
      rc("GENIE3", sb_sym_max(W), "RPE1 GENIE3 (baselines.R settings), symmetrised by max"); rm(W); invisible(gc()) }
    f <- file.path(ww, "grnboost2_edges.tsv")
    if (file.exists(f)) { e <- utils::read.delim(f, stringsAsFactors = FALSE); M <- matrix(0, Gu, Gu); ia <- match(e$TF, gu); ib <- match(e$target, gu); ok <- !is.na(ia) & !is.na(ib)
      gate("rpe1", "universe", "every GRNBoost2 edge endpoint is in the gene list (1 = yes)", as.numeric(all(ok)), 1, 0); M[cbind(ia[ok], ib[ok])] <- e$importance[ok]; rm(e)
      rc("GRNBoost2", sb_sym_max(M), "RPE1 GRNBoost2 (tools/grnboost2.py, legacy early stopping), symmetrised by max"); rm(M); invisible(gc()) }
    f <- file.path(ww, "pidc.f64")
    if (file.exists(f)) { v <- readBin(f, "double", Gu * (Gu - 1) / 2, size = 8); gate("rpe1", "universe", "PIDC packed file has G(G-1)/2 finite weights (1 = yes)", as.numeric(length(v) == Gu * (Gu - 1) / 2 && all(is.finite(v))), 1, 0)
      M <- matrix(0, Gu, Gu); M[lower.tri(M)] <- v; rm(v); M <- M + t(M)      # row-major upper triangle = column-major lower triangle
      rc("PIDC", M, "RPE1 PIDC (tools/pidc_chunked.jl, bit-identical to the package, gates in pidc_gate/)"); rm(M); invisible(gc()) }
  }
  rm(X, ld); invisible(gc()); sb_say("peak RSS so far ", round(hwm(), 2), " GB")
}

## ---- yeast: DREAM5 network 4 expression, Complex Portal main truth (frozen) --------------------------------------------------
if ("yeast" %in% parts) {
  sb_say("== yeast (DREAM5 net4)"); source("analysis/complex_confirm/cc_lib.R")
  nd <- cc_load_net(); G <- nrow(nd$X); tf <- cc_read_truth_files(); T <- cc_build_truths(nd, tf, ym = FALSE); main <- T$main
  tc <- utils::read.csv(file.path(CC$out, "net4/truth/truth_counts.csv"), stringsAsFactors = FALSE); tc <- tc[tc$truth == "main", ]
  gate("yeast", "b", "truth structure md5 (i, j, label) equals truth_counts.csv (1 = yes)", as.numeric(identical(tc$structure_md5, cc_truth_struct_hash(main))), 1, 0)
  gate("yeast", "b", "truth positive pairs vs truth_counts.csv", sum(main$label), tc$n_pos, 0); gate("yeast", "b", "truth pairs vs truth_counts.csv", length(main$label), tc$n_pairs, 0)
  nclr <- function(S) cc_score_truth(list(x = S), main, NULL, G, 1L)$aupr[[1]]          # the study's own scoring function
  dd <- file.path(D5$out, "net4"); rds <- function(m) d5_read_scores(d5_cache_path(4L, m))
  sha_ok <- function(m) { f <- d5_cache_path(4L, m); sc <- utils::read.delim(paste0(f, ".sha256"), stringsAsFactors = FALSE, colClasses = c(sha256 = "character")); as.numeric(identical(d5_sha256(f), sc$sha256)) }
  for (m in c("genie3_all", "grnboost2_all", "pidc", "abs_pearson", "clr")) gate("yeast", "universe", paste0(m, " cache: SHA-256 equals its sidecar (1 = yes)"), sha_ok(m), 1, 0)
  mm <- utils::read.csv(file.path(CC$out, "net4/truth/metrics.csv"), stringsAsFactors = FALSE); mm <- mm[mm$truth == "main" & mm$scope == "raw" & mm$readout == "clr", ]
  C <- d5_unpack(d5_read_scores(cc_cache_files(CC$out)[1]), G); a_clr <- nclr(C)
  gate("yeast", "b", "CLR (cc_clr.scores) AUPR through cc_score_truth vs study metrics.csv", a_clr, mm$aupr)
  gate("yeast", "b", "CLR enrichment vs the scoreboard CLR cell", a_clr / mean(main$label), sd_cell("dream5_yeast", "complexes_main", "CLR"), 1e-9)
  gate("yeast", "b", "cc_clr.scores equals results/dream5/net4/clr__replicate.scores (max abs diff)", max(abs(d5_read_scores(cc_cache_files(CC$out)[1]) - rds("clr"))), 0, 1e-12)
  gate("yeast", "b", "pair-extraction path (score()) equals cc_score_truth on CLR", score(C, main), a_clr, 1e-15); rm(C)
  # (c) |r|: recomputed by the pipeline's own function vs its cache, then scored
  rr <- d5_m_abs_pearson(nd$X); gate("yeast", "c", "|r| recomputed by d5_m_abs_pearson equals results/dream5/net4/abs_pearson__replicate.scores (max abs diff)",
                                     max(abs(d5_pack(rr) - rds("abs_pearson"))), 0, 1e-12); a_r <- nclr(rr); rm(rr)
  cell("yeast", "complexes_main", "CLR", a_clr, main, "cc_clr.scores (gate b)", "gate")
  cell("yeast", "complexes_main", "abs_r", a_r, main, "results/dream5/net4/abs_pearson__replicate.scores (recomputed, gate c)")
  cell("yeast", "complexes_main", "PIDC", nclr(d5_unpack(rds("pidc"), G)), main, "results/dream5/net4/pidc__replicate.scores")
  dirm <- c(GENIE3 = "genie3_all", GRNBoost2 = "grnboost2_all")
  for (nm in names(dirm)) { m <- dirm[[nm]]; W <- d5_read_scores(d5_cache_path(4L, m)); stopifnot(identical(dim(W), c(G, G)))
    cell("yeast", "complexes_main", nm, nclr(sb_sym_max(W)), main, paste0("results/dream5/net4/", m, "__replicate.scores, symmetrised by max")); rm(W); invisible(gc()) }
  sb_say("peak RSS so far ", round(hwm(), 2), " GB")
}

## ---- write -------------------------------------------------------------------------------------------------------------
gt <- do.call(rbind, G_); utils::write.csv(gt, file.path(OPT$out, "gates.csv"), row.names = FALSE)
nfail <- sum(!gt$pass); sb_say(sprintf("%d gates, %d failed", nrow(gt), nfail))
if (nfail > 0L) { print(gt[!gt$pass, ], row.names = FALSE); stop("GATE FAILURE: no cell written (gates.csv holds the report)") }
pt <- do.call(rbind, P_); utils::write.csv(pt, file.path(OPT$out, "pooled.csv"), row.names = FALSE)
print(pt[pt$role == "new", c("dataset", "truth", "method", "enrichment", "base_rate", "n_pos")], row.names = FALSE, digits = 6)
sb_say(sprintf("done in %.1f min; peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "mins")), hwm()))
