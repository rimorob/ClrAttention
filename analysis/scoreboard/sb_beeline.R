#!/usr/bin/env Rscript
# Scoreboard (docs/plans/scoreboard_plan.md), BEELINE job. DESCRIPTIVE: point estimates, no intervals, no tests.
# The 14 cases (7 datasets x TFs+500 / TFs+1000), BEELINE's protocol and truths, scoring by beeline_baselines.R::score_truth (evaluated from that
# file, not copied). Operator = CLR-NS (beeline.R::ns_clr). Attention (the BEELINE method) = soft10 attention on CLR-NS at the held-out depth t*
# of results/beeline_v2/heldout_depths.csv (the existing protocol; attention at t = 9 is also carried). a3r family: t = 9 in every layer, same
# softmax / lazy walk / readout as analysis/explore_stack/stack_lib.R, stacked on CLR-NS; not tuned (design choice recorded in the plan).
# Gate first: |r|, CLR and CLR-NS reproduce results/beeline_baselines/metrics.csv (every truth, AUPRC ratio and EPR, 1e-9); attention on CLR-NS at t*
# reproduces att_clr_absr_blomavg of results/beeline_blom_test/metrics.csv (ap_fast scoring, 1e-9); generic layer-1 vs package operator 1e-12.
# Usage (repo root): Rscript analysis/scoreboard/sb_beeline.R [--out results/scoreboard/beeline] [--threads 4]
args0 <- commandArgs(TRUE); bopt <- list(out = "results/scoreboard/beeline", threads = 4L)
for (k in seq(1, length(args0), by = 2)) { key <- sub("^--", "", args0[k]); bopt[[key]] <- if (key == "threads") as.integer(args0[k + 1]) else args0[k + 1] }
sb_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
# beeline.R's setup and helpers (everything before its main loop), as beeline_baselines.R does; script_dir is pinned to analysis/
src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]
args <- c("--out", bopt$out, "--threads", as.character(bopt$threads))
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  d1 <- deparse(e)[1]
  if (identical(d1, "args <- commandArgs(trailingOnly = TRUE)")) next
  if (startsWith(d1, "script_dir <-")) { script_dir <- "analysis"; next }
  eval(e, globalenv())
}
# beeline_baselines.R::score_truth (tie-aware EPR), evaluated from its source
bsrc <- readLines("analysis/beeline_baselines.R")
a <- grep("^score_truth <- function", bsrc)[1]; z <- grep("^say\\(\"beeline_baselines: stage", bsrc)[1] - 1L
eval(parse(text = bsrc[a:z]), globalenv())
source(file.path(sb_dir, "sb_lib.R")); source("analysis/explore_stack/stack_lib.R"); source("analysis/robust_z/robust_z_lib.R")
sb_say2 <- function(...) sb_say(file.path(bopt$out, "run.log"), ...)
sb_say2("options: ", paste(names(bopt), unlist(bopt), sep = "=", collapse = " "))

Ns_all <- as.integer(strsplit(opt$N, ",")[[1]]); dss <- strsplit(opt$datasets, ",")[[1]]
cases <- expand.grid(N = Ns_all, ds = dss, stringsAsFactors = FALSE)
case_X <- function(ds, N) { ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1, check.names = FALSE))
  g <- intersect(select_genes(ds, N)$genes, rownames(ex)); X <- ex[g, , drop = FALSE]; X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE] }
dep <- utils::read.csv("results/beeline_v2/heldout_depths.csv", stringsAsFactors = FALSE)
# existing attention path of beeline_blom_test.R (soft10 operator, lazy walk, symmetrised mass)
att_blom <- function(S, t) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L)
  P <- Matrix::Diagonal(nrow(A), 1 - 0.5) + 0.5 * Matrix::Matrix(A, sparse = TRUE); Pt <- as.matrix(P)
  if (t > 1) for (s in 2:t) Pt <- as.matrix(Pt %*% P)
  M <- (Pt + t(Pt)) / 2; diag(M) <- 0; M }
edge_sets <- function(genes, gt) { gi <- match(toupper(gt$Gene1), toupper(genes)); gj <- match(toupper(gt$Gene2), toupper(genes))
  ok <- !is.na(gi) & !is.na(gj) & gi != gj; regs <- sort(unique(gi[ok])); if (length(regs) < 1 || sum(ok) < 10) return(NULL)
  G <- length(genes); cand <- cbind(rep(regs, each = G), rep(seq_len(G), times = length(regs))); cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
  list(cand = cand, lab = as.integer((cand[, 1] * (G + 1) + cand[, 2]) %in% unique(gi[ok] * (G + 1) + gj[ok]))) }

rows <- list(); gate_apfast <- list(); op_diff <- 0
for (k in seq_len(nrow(cases))) {
  ds <- cases$ds[k]; N <- cases$N[k]; X <- case_X(ds, N); genes <- rownames(X); t0 <- Sys.time()
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  tstar <- dep$t_soft10[dep$dataset == ds & dep$N == N]; stopifnot(length(tstar) == 1L)
  S <- list(abs_pearson = abs(stats::cor(t(X))))
  f <- clr_fit(X); S$clr <- f$clr_scores; f$release(FALSE)
  S$clr_ns <- ns_clr(X)
  Cn <- S$clr_ns
  S$attention_tstar <- att_blom(Cn, tstar); S$attention_t9 <- att_blom(Cn, 9L)
  L1 <- sk_layer(Cn, ts_out = 9L); S$a3r <- L1$outs[["9"]]; rm(L1)
  L2 <- sk_layer(S$a3r, ts_out = 9L); S$a3r2 <- L2$outs[["9"]]; rm(L2)
  L3 <- sk_layer(S$a3r2, ts_out = 9L); S$a3r3 <- L3$outs[["9"]]; rm(L3)
  Apkg <- { ff <- clr_fit(X); ff$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); o <- ff$operator; ff$release(FALSE); o }   # package operator on the CLR scores (gate (e) form)
  op_diff <- max(op_diff, max(abs(rz_operator(S$clr) - Apkg))); rm(Apkg)
  for (m in names(S)) { s <- S[[m]]; s[!is.finite(s)] <- 0; diag(s) <- 0; S[[m]] <- s }
  for (tn in names(gts)) for (m in names(S)) { r <- score_truth(S[[m]], genes, NULL, gts[[tn]]); if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(dataset = ds, N = N, truth = tn, method = m, t(r)) }
  # ap_fast scoring of attention on CLR-NS at t*, for the gate against beeline_blom_test (that file scores with ap_fast)
  for (tn in names(gts)) { e <- edge_sets(genes, gts[[tn]]); if (is.null(e)) next
    s <- S$attention_tstar[e$cand]; gate_apfast[[length(gate_apfast) + 1L]] <- data.frame(case = sprintf("%s_%d", ds, N), truth = tn, auprc_ratio = ap_fast(s, e$lab) / mean(e$lab)) }
  sb_say2(sprintf("%s TFs+%d: %d genes x %d cells, t* = %d (%.1f min)", ds, N, length(genes), ncol(X), tstar, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(S, X, Cn); invisible(gc(verbose = FALSE))
}
res <- do.call(rbind, rows); ga <- do.call(rbind, gate_apfast)

## ---- equivalence gate -----------------------------------------------------------------------------------------------------
g <- list()
bb <- utils::read.csv("results/beeline_baselines/metrics.csv", stringsAsFactors = FALSE)
for (m in c("abs_pearson", "clr", "clr_ns")) { a <- res[res$method == m, ]; b <- bb[bb$method == m, ]; key <- function(x) paste(x$dataset, x$N, x$truth)
  stopifnot(setequal(key(a), key(b)))
  j <- match(key(a), key(b))
  g <- sb_gate_add(g, sprintf("%s: AUPRC ratio and EPR, %d case x truth cells, vs results/beeline_baselines/metrics.csv", m, nrow(a)), max(abs(a$auprc_ratio - b$auprc_ratio[j]), abs(a$epr - b$epr[j]))) }
bl <- utils::read.csv("results/beeline_blom_test/metrics.csv", stringsAsFactors = FALSE); bl <- bl[bl$method == "att_clr_absr_blomavg" & bl$seed == 1, ]
j <- match(paste(ga$case, ga$truth), paste(bl$case, bl$truth)); stopifnot(!anyNA(j))
g <- sb_gate_add(g, sprintf("attention on CLR-NS at held-out t*: AUPRC ratio (ap_fast), %d case x truth cells, vs beeline_blom_test att_clr_absr_blomavg", nrow(ga)), max(abs(ga$auprc_ratio - bl$auprc_ratio[j])))
g <- sb_gate_add(g, "rz_operator vs package operator on the CLR scores, 14 cases", op_diff, 1e-12)
sb_gate_finish(g, file.path(bopt$out, "equivalence_gate.csv"), sb_say2)

## ---- summaries (as beeline_baselines.R): median over complete primary case x truth cells; median over the 14 cell-type ChIP cases ----------
utils::write.csv(res, file.path(bopt$out, "metrics.csv"), row.names = FALSE)
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
cells <- unique(res[, c("dataset", "N", "truth")]); complete <- merge(cells, stats::aggregate(method ~ dataset + N + truth, res, length))
complete <- complete[complete$method == length(unique(res$method)), c("dataset", "N", "truth")]; rc <- merge(res, complete)
sm <- do.call(rbind, lapply(c("primary", "celltype_ChIP"), function(ts) { x <- rc[if (ts == "primary") rc$truth %in% primary else rc$truth == ts, ]
  do.call(rbind, lapply(split(x, x$method), function(y) data.frame(truths = ts, method = y$method[1], n_cells = nrow(y), median_auprc_ratio = stats::median(y$auprc_ratio), median_epr = stats::median(y$epr), stringsAsFactors = FALSE))) }))
utils::write.csv(sm, file.path(bopt$out, "summary.csv"), row.names = FALSE)
sb_say2("median AUPRC ratio:"); print(sm, digits = 4, row.names = FALSE)
sb_say2("done")
