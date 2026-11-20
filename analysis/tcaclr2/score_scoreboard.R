#!/usr/bin/env Rscript
# tcaCLR2 on the scoreboard rows, step 2: every scoreboard metric for tcaCLR2 exactly as for tcaCLR (descriptive point estimates; no interval, no test; exploratory, outside the
# preregistration). One family per call:  m3d | ecomac | k562 | rpe1 | beeline.  Rule of tcaCLR2 (docs/plans/tcaclr2_definition.md) unchanged: depth t2 from the edge TSV of
# edge_depth.R (m3d, k562: results/tcaclr2/) or edge_scoreboard.R (ecomac, rpe1, beeline: results/tcaclr2_scoreboard/); no truth is read to get t2.
# Smoke-test switches (never used for the real run): --t2 N overrides the depth, --quick N (k562) / --genes N (rpe1) take a gene slice, --gate 0 does not stop on the gate, --tag _x renames the outputs.
#   Rscript analysis/tcaclr2/score_scoreboard.R --family m3d [--edge results/tcaclr2/edge_m3d_tf_regulons.tsv] [--out results/tcaclr2_scoreboard] [--threads 8]
# Machinery reused, not re-implemented: operator rz_operator / rz_lazy_P (robust_z_lib.R, checked equal to the package operator as sb_prep does), walk read-out sk_a3r (stack_lib.R; K562: the
# identical a3r of ideas47_lib.R), pooled AUPR with raw and hub_att scopes sb_pooled, per-TF AUROC sb_per_tf_median_auroc, regulon_ap_table, the truth builders of sb_load_ecoli,
# build_truths_ecoli / build_truths_human (clr_logpt_test.R definitions as imported by sb_cx_lib.R), ideas47_lib.R (K562), and for BEELINE beeline.R's set-up, ns_clr and
# beeline_baselines.R::score_truth. Walk: the registered sequential chain (one sparse product per step), read at t = 9 and t = t2 from the same chain.
# Gate (before anything is written): tcaCLR at t = 9 recomputed through this code reproduces every stored scoreboard value of its cells (key a3r of docs/whitepapers/figures/scoreboard_data.csv) to 1e-9.
# Output (--out): cells_<family>.csv (scoreboard_data schema; methods tcaCLR_t9 = the recomputed reference, tcaCLR2), gate_<family>.csv, score_<family>.log. Writes nothing else; refuses frozen paths.
suppressPackageStartupMessages({ source("analysis/tcaclr2/tsb_lib.R") })
topt <- tsb_args(list(family = "", edge = "", edge_dir = TSB$out, out = TSB$out, threads = 8L, t2 = 0L, hub = 1L, tag = "", quick = 0L, genes = 0L, gate = 1L))
FAM <- topt$family; if (!FAM %in% c("m3d", "ecomac", "k562", "rpe1", "beeline")) stop("--family must be m3d, ecomac, k562, rpe1 or beeline")
tsb_refuse_frozen(topt$out); dir.create(topt$out, recursive = TRUE, showWarnings = FALSE)
say <- tsb_say(file.path(topt$out, sprintf("score_%s%s.log", FAM, topt$tag))); say("options: ", paste(names(topt), unlist(topt), sep = "=", collapse = " "))
SBD <- utils::read.csv(TSB$sb_csv, stringsAsFactors = FALSE); PM <- "pooled_aupr_over_base"; T0 <- Sys.time()
SCOPES <- if (topt$hub > 0L) c("raw", "hub_att") else "raw"
edge_file <- function(default) if (nzchar(topt$edge)) topt$edge else default
get_t2 <- function(f) { if (topt$t2 > 0L) { say("t2 overridden on the command line: ", topt$t2, " (smoke test only)"); return(list(t2 = topt$t2, edge = NA_real_)) }
  e <- tsb_edge_read(f); say(sprintf("t2 = %d from %s (edge %.5f)", e$t2, f, e$edge)); e }
rows <- list(); add <- function(r) rows[[length(rows) + 1L]] <<- r
note_for <- function(m, t2, edge) if (m == "tcaCLR2") sprintf("tcaCLR2: depth t2 = %d (edge %.5f; rule of docs/plans/tcaclr2_definition.md); exploratory, outside the preregistration", t2, edge) else "tcaCLR at t = 9 recomputed by this code (gate against the stored scoreboard value)"

# ---- one walk, two read-outs: P -> {tcaCLR_t9, tcaCLR2} --------------------------------------------------------------------------------
readouts <- function(P, t2, a3r_fn) { o <- tsb_walk(P, c(9L, t2), function(t, Pt) a3r_fn(Pt)); list(tcaCLR_t9 = o[["9"]], tcaCLR2 = o[[as.character(t2)]]) }
clean <- function(S) { S[!is.finite(S)] <- 0; diag(S) <- 0; S }

# ---- the E. coli families and the complex-truth families share: pooled raw/hub_att on a bench list -----------------------------------------
pooled_rows <- function(ds, S, bench, truth_name, H_att, m, t2, edge, src) {
  p <- sb_pooled(m, S, bench, H_att, NULL, scopes = SCOPES)
  for (k in seq_len(nrow(p))) add(tsb_row(ds, truth_name[[p$truth[k]]], PM, p$scope[k], m, p$enrichment[k], p$base_rate[k], src, note_for(m, t2, edge)))
}

if (FAM %in% c("m3d", "ecomac")) {
  source("analysis/scoreboard/sb_lib.R"); sb_source(); source("analysis/scoreboard/sb_cx_lib.R"); sb_import_truth_defs()
  ds <- FAM; et <- get_t2(edge_file(if (ds == "m3d") "results/tcaclr2/edge_m3d_tf_regulons.tsv" else file.path(TSB$out, "edge_ecomac.tsv"))); t2 <- et$t2
  D <- sb_load_ecoli(ds); X <- D$X; genes <- D$genes; G <- length(genes); bench <- D$bench; U <- D$U
  TR <- build_truths_ecoli(ds, genes, NULL, only = c("complexes", "regulons_SC", "regulons_all")); bench$complexes <- TR$bench$complexes
  say(sprintf("%s: %d genes x %d samples; truths %s", ds, G, ncol(X), paste(sprintf("%s (%d pairs, %d pos)", names(bench), vapply(bench, function(b) length(b$pairs$label), 0L), vapply(bench, function(b) sum(b$pairs$label), 0)), collapse = "; ")))
  fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = topt$threads); fit$calibrate(method = "normal", combine = "stouffer")
  C <- fit$clr_scores; fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator; fit$release(); rm(fit, X, D, TR)
  A <- rz_operator(C); eq_op <- max(abs(A - A_pkg)); rm(A_pkg, C); say(sprintf("max|A(rz_operator) - package operator| = %.3g", eq_op)); stopifnot(eq_op < 1e-12)
  P <- rz_lazy_P(A); rm(A); q <- rz_stationary(P); H_att <- outer(q, q, "+"); say(sprintf("operator and stationary distribution: %.0f s", as.numeric(difftime(Sys.time(), T0, units = "secs"))))
  SS <- readouts(P, t2, sk_a3r); rm(P); invisible(gc(verbose = FALSE)); say(sprintf("walk to t = 9 and t = %d done; peak RSS %.1f GB", t2, tsb_hwm()))
  src <- sprintf("analysis/tcaclr2/score_scoreboard.R --family %s", ds); tn <- c(SC = "SC", all = "all", complexes = "complexes")
  for (m in names(SS)) { S <- clean(SS[[m]]); SS[[m]] <- NULL
    pooled_rows(ds, S, bench, tn, H_att, m, t2, et$edge, src)
    for (nm in names(U)) { u <- U[[nm]]
      add(tsb_row(ds, nm, "tf_aupr_ratio", "raw", m, ap_fast(S[cbind(u$i, u$j)], u$label) / mean(u$label), NA, src, note_for(m, t2, et$edge)))
      add(tsb_row(ds, nm, "tf_auroc_median", "raw", m, sb_per_tf_median_auroc(S, u, G), NA, src, note_for(m, t2, et$edge))) }
    for (nm in c("SC", "all")) { b <- bench[[nm]]; x <- regulon_ap_table(S, b$reg, b$co, b$operon, genes); y <- x[!is.na(x$ap), ]
      add(tsb_row(ds, nm, "per_regulon_mean_ap_over_base", "raw", m, mean(y$ap) / mean(y$base_rate), mean(y$base_rate), src, note_for(m, t2, et$edge)))
      add(tsb_row(ds, nm, "per_regulon_mean_ap", "raw", m, mean(y$ap), mean(y$base_rate), src, note_for(m, t2, et$edge)))
      add(tsb_row(ds, nm, "per_regulon_mean_ap_norm", "raw", m, mean((y$ap - y$base_rate) / (1 - y$base_rate)), mean(y$base_rate), src, note_for(m, t2, et$edge))) }
    say("scored ", m); rm(S); invisible(gc(verbose = FALSE)) }
}

if (FAM == "k562") {
  et <- get_t2(edge_file("results/tcaclr2/edge_k562_corum.tsv")); t2 <- et$t2
  opt <- list(out = file.path(topt$out, "work_k562"), threads = topt$threads, workers = 4L, jack = 1L, quick = topt$quick)      # as sb_k562.R (no jackknife is used: jack = 1)
  source("analysis/explore_ideas47/ideas47_lib.R"); source("analysis/explore_stack/stack_lib.R"); source("analysis/scoreboard/sb_lib.R"); unlink(file.path(opt$out, "cache"), recursive = TRUE)
  if (topt$quick == 0L) stopifnot(G == 5671L)
  abc <- dorothea_sets(genes, c("A", "B", "C")); bench$DoRothEA_ABC <- list(reg = abc, pairs = comembership_pairs(abc, none, genes))
  rm(halves, M, S9, A); invisible(gc(verbose = FALSE))
  S9a <- a3r(Pt9); Pt <- Pt9; rm(Pt9); invisible(gc(verbose = FALSE))                                                    # t = 9 chain of ideas47_lib.R, continued sequentially to t2
  if (t2 >= 9L) { if (t2 > 9L) for (s in 10:t2) { Pt <- as.matrix(Pt %*% P); if (s %% 8L == 0L) invisible(gc(verbose = FALSE)) }; S2 <- a3r(Pt) } else stop("t2 < 9 not supported for K562 (chain starts at t = 9)")
  rm(Pt, P); invisible(gc(verbose = FALSE)); SS <- list(tcaCLR_t9 = S9a, tcaCLR2 = S2); rm(S9a, S2)
  say(sprintf("walk to t = 9 and t = %d done; peak RSS %.1f GB", t2, tsb_hwm()))
  tn <- c(DoRothEA_AB = "DoRothEA_AB", DoRothEA_ABC = "DoRothEA_ABC", CollecTRI = "CollecTRI", CORUM_all = "CORUM", KD_response = "KD_response", CORUM_noRP = "CORUM_noRP", CORUM_RP = "CORUM_RP",
          CORUM_small = "CORUM_small", CORUM_large = "CORUM_large", CORUM_merged = "CORUM_merged", CORUM_merged_noRP = "CORUM_merged_noRP")
  stopifnot(setequal(names(bench), names(tn)))
  for (m in names(SS)) { S <- clean(SS[[m]]); SS[[m]] <- NULL; pooled_rows("k562", S, bench, tn, H_att, m, t2, et$edge, "analysis/tcaclr2/score_scoreboard.R --family k562"); say("scored ", m); rm(S); invisible(gc(verbose = FALSE)) }
}

if (FAM == "rpe1") {
  source("analysis/scoreboard/sb_lib.R"); sb_source(); source("analysis/scoreboard/sb_cx_lib.R"); sb_import_truth_defs()
  et <- get_t2(edge_file(file.path(TSB$out, "edge_rpe1.tsv"))); t2 <- et$t2
  X <- load_expr("rpe1")$X
  if (topt$genes > 0L && topt$genes < nrow(X)) { set.seed(20261121L); X <- X[sort(sample(nrow(X), topt$genes)), , drop = FALSE] }   # smoke only
  genes <- rownames(X); G <- length(genes)
  TR <- build_truths_human("rpe1", genes, NULL, NULL, NULL); bench <- list(CORUM_merged = TR$bench$CORUM_merged); rm(TR)
  say(sprintf("rpe1: %d genes x %d profiles; truth CORUM_merged %d pairs, %d positive", G, ncol(X), length(bench[[1]]$pairs$label), sum(bench[[1]]$pairs$label)))
  fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = topt$threads); fit$calibrate(method = "normal", combine = "stouffer")
  C <- fit$clr_scores; fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator; fit$release(); rm(fit, X)
  A <- rz_operator(C); eq_op <- max(abs(A - A_pkg)); rm(A_pkg, C); say(sprintf("max|A(rz_operator) - package operator| = %.3g", eq_op)); stopifnot(eq_op < 1e-12)
  P <- rz_lazy_P(A); rm(A); q <- rz_stationary(P); H_att <- outer(q, q, "+")
  SS <- readouts(P, t2, sk_a3r); rm(P); invisible(gc(verbose = FALSE)); say(sprintf("walk to t = 9 and t = %d done; peak RSS %.1f GB", t2, tsb_hwm()))
  for (m in names(SS)) { S <- clean(SS[[m]]); SS[[m]] <- NULL; pooled_rows("rpe1", S, bench, c(CORUM_merged = "CORUM_merged"), H_att, m, t2, et$edge, "analysis/tcaclr2/score_scoreboard.R --family rpe1"); say("scored ", m); rm(S); invisible(gc(verbose = FALSE)) }
}

if (FAM == "beeline") {
  ob <- topt$out; args <- c("--out", ob, "--threads", as.character(topt$threads))                                          # beeline.R's own parser reads `args` (as sb_beeline.R)
  src0 <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src0)[1]
  for (e in parse(text = src0[seq_len(cut_at - 1L)])) { d1 <- deparse(e)[1]
    if (identical(d1, "args <- commandArgs(trailingOnly = TRUE)")) next
    if (startsWith(d1, "script_dir <-")) { script_dir <- "analysis"; next }
    eval(e, globalenv()) }
  bsrc <- readLines("analysis/beeline_baselines.R"); a <- grep("^score_truth <- function", bsrc)[1]; z <- grep("^say\\(\"beeline_baselines: stage", bsrc)[1] - 1L
  eval(parse(text = bsrc[a:z]), globalenv())                                                                              # beeline_baselines.R::score_truth, evaluated from its source
  source("analysis/scoreboard/sb_lib.R"); source("analysis/explore_stack/stack_lib.R"); source("analysis/robust_z/robust_z_lib.R")
  say <- tsb_say(file.path(topt$out, sprintf("score_%s%s.log", FAM, topt$tag)))                                           # beeline.R redefines say
  Ns <- as.integer(strsplit(opt$N, ",")[[1]]); dss <- strsplit(opt$datasets, ",")[[1]]; cases <- expand.grid(N = Ns, ds = dss, stringsAsFactors = FALSE)
  stored <- utils::read.csv("results/scoreboard/beeline/metrics.csv", stringsAsFactors = FALSE); stored <- stored[stored$method == "a3r", ]
  case_X <- function(ds, N) { ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1, check.names = FALSE))
    g <- intersect(select_genes(ds, N)$genes, rownames(ex)); X <- ex[g, , drop = FALSE]; X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE] }
  res <- list(); tcase <- list()
  for (k in seq_len(nrow(cases))) {
    ds <- cases$ds[k]; N <- cases$N[k]; id <- sprintf("beeline_%s_%d", ds, N); t1 <- Sys.time()
    et <- get_t2(edge_file_beeline <- file.path(topt$edge_dir, sprintf("edge_%s.tsv", id))); t2 <- et$t2; tcase[[id]] <- data.frame(case = id, t2 = t2, edge = et$edge)
    X <- case_X(ds, N); genes <- rownames(X); gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE)); Cn <- ns_clr(X); rm(X)
    L <- sk_layer(Cn, ts_out = sort(unique(c(9L, t2)))); SS <- list(tcaCLR_t9 = L$outs[["9"]], tcaCLR2 = L$outs[[as.character(t2)]]); rm(L, Cn)
    for (m in names(SS)) { s <- SS[[m]]; s[!is.finite(s)] <- 0; diag(s) <- 0; SS[[m]] <- s }
    for (tn in names(gts)) for (m in names(SS)) { r <- score_truth(SS[[m]], genes, NULL, gts[[tn]]); if (is.null(r)) next
      res[[length(res) + 1L]] <- data.frame(dataset = ds, N = N, truth = tn, method = m, t(r)) }
    say(sprintf("%s: %d genes, t2 = %d (%.1f s)", id, length(genes), t2, as.numeric(difftime(Sys.time(), t1, units = "secs")))); rm(SS); invisible(gc(verbose = FALSE))
  }
  res <- do.call(rbind, res); utils::write.csv(res, file.path(topt$out, sprintf("beeline_metrics%s.csv", topt$tag)), row.names = FALSE)
  key <- function(x) paste(x$dataset, x$N, x$truth)
  primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF"); okc <- unique(key(stored))                          # the complete case x truth cells of the stored scoreboard run
  for (ts in c("primary", "celltype_ChIP")) { mt <- if (ts == "primary") "beeline_median_auprc_ratio_primary" else "beeline_median_auprc_ratio_celltypeChIP"
    for (m in c("tcaCLR_t9", "tcaCLR2")) { x <- res[res$method == m & key(res) %in% okc & (if (ts == "primary") res$truth %in% primary else res$truth == ts), ]
      add(tsb_row("beeline", "all_14_cases", mt, "raw", m, stats::median(x$auprc_ratio), NA, "analysis/tcaclr2/score_scoreboard.R --family beeline",
        paste0("median over the complete case x truth cells (n = ", nrow(x), "); CLR-NS operator, ", if (m == "tcaCLR2") "t2 per case from the permuted-data edge of the CLR-NS operator (edge_scoreboard.R)" else "t = 9 in every case"))) } }
  utils::write.csv(do.call(rbind, tcase), file.path(topt$out, sprintf("beeline_t2%s.csv", topt$tag)), row.names = FALSE)
}

cells <- do.call(rbind, rows); g <- tsb_gate(cells, SBD, 1e-9, function(s) say(sub("\n$", "", s)))
utils::write.csv(g, file.path(topt$out, sprintf("gate_%s%s.csv", FAM, topt$tag)), row.names = FALSE); print(g, row.names = FALSE, digits = 7)
if (!all(g$pass) && topt$gate > 0L) stop("GATE FAILED: tcaCLR (t = 9) recomputed by this code does not reproduce the stored scoreboard value; no tcaCLR2 cell is written")
utils::write.csv(cells, file.path(topt$out, sprintf("cells_%s%s.csv", FAM, topt$tag)), row.names = FALSE)
say("tcaCLR2 cells:"); print(cells[cells$method == "tcaCLR2", c("dataset", "truth", "metric", "scope", "enrichment_over_random")], row.names = FALSE, digits = 5)
say(sprintf("done: %.1f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), tsb_hwm()))
