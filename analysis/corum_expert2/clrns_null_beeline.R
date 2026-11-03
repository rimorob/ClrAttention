#!/usr/bin/env Rscript
# docs/plans/corum_expert2.md, Task B. EXPLORATORY, NOT PREREGISTERED. Nothing from this
# script goes into paper/ or any paper table.
#
# Stage 1 (label-free): for each of the 14 BEELINE cases (7 datasets x TFs+500/1000, the
#   gene sets and expression of analysis/beeline_baselines.R), row diagnostics of |r| of
#   Blom normal scores (ties averaged: CLR-NS's input) and of B-spline MI against a
#   within-gene permutation null (cells permuted independently per gene, one draw; seed
#   fixed), and how much each label-free alternative row calibration changes the ranking
#   relative to CLR-NS. Writes diagnostics.csv and rank_change.csv before any truth is read.
# Stage 2 (truth-based, EXPLORATORY: CLR-NS was itself chosen on these truths, so any
#   comparison here is in-sample): AUPRC ratio of every alternative on BEELINE's primary
#   truths (tie-aware score_truth of beeline_baselines.R), reproduction check of CLR and
#   CLR-NS against results/beeline_baselines/metrics.csv, and the blocked-by-dataset
#   summary (mean log2 ratio to CLR-NS per dataset, exact Wilcoxon, sign count).
#
# Usage (repo root): Rscript analysis/corum_expert2/clrns_null_beeline.R [--out results/corum_expert2/beeline] [--threads 4] [--stage 1|2]
args0 <- commandArgs(trailingOnly = TRUE)
bopt <- list(out = "results/corum_expert2/beeline", threads = 4L, stage = "2", seed = 20261103L)
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i]); bopt[[k]] <- if (k %in% c("threads", "seed")) as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
dir.create(bopt$out, recursive = TRUE, showWarnings = FALSE)

# beeline.R's setup and helpers (everything before its main loop), as beeline_baselines.R does
src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]
args <- c("--out", bopt$out, "--threads", as.character(bopt$threads))
if (!is.null(bopt$datasets)) args <- c(args, "--datasets", bopt$datasets)
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  if (identical(deparse(e)[1], "args <- commandArgs(trailingOnly = TRUE)")) next
  if (grepl("^script_dir <-", deparse(e)[1])) { script_dir <- "analysis"; next }   # helpers live in analysis/
  eval(e, globalenv())
}
source("analysis/corum_expert2/null_diag.R")
# Tie-aware early precision and AUPRC ratio: verbatim from analysis/beeline_baselines.R
score_truth <- function(S, genes, tfs_in, gt) {
  gi <- match(toupper(gt$Gene1), toupper(genes)); gj <- match(toupper(gt$Gene2), toupper(genes))
  ok <- !is.na(gi) & !is.na(gj) & gi != gj
  regs <- sort(unique(gi[ok]))
  if (length(regs) < 1 || sum(ok) < 10) return(NULL)
  G <- length(genes)
  cand <- cbind(rep(regs, each = G), rep(seq_len(G), times = length(regs)))
  cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
  key <- cand[, 1] * (G + 1) + cand[, 2]; pos <- unique(gi[ok] * (G + 1) + gj[ok])
  lab <- as.integer(key %in% pos); s <- S[cand]
  ps <- pr_summary(s, lab); dens <- mean(lab); k <- sum(lab)
  v <- sort(s, decreasing = TRUE)[k]; above <- s > v; tie <- s == v
  prec_k <- (sum(lab[above]) + (k - sum(above)) * sum(lab[tie]) / sum(tie)) / k
  c(n_regulators = length(regs), n_edges = k, density = dens,
    auprc_ratio = ps[["aupr"]] / dens, epr = prec_k / dens, auroc = ps[["auroc"]])
}
case_X <- function(ds, N) {     # as beeline_baselines.R
  ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1, check.names = FALSE))
  g <- intersect(select_genes(ds, N)$genes, rownames(ex)); X <- ex[g, , drop = FALSE]
  X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]
}
Ns_all <- c(500L, 1000L); dss <- strsplit(opt$datasets, ",")[[1]]
cases <- expand.grid(N = Ns_all, ds = dss, stringsAsFactors = FALSE)
mi_of <- function(X) { f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = bopt$threads)
  m <- f$mi; f$release(FALSE); diag(m) <- 0; m }

dg <- list(); rc <- list(); rows <- list(); cache <- list()
for (k in seq_len(nrow(cases))) {
  ds <- cases$ds[k]; N <- cases$N[k]; lab <- sprintf("%s_%d", ds, N); t0 <- Sys.time()
  X <- case_X(ds, N); genes <- rownames(X); G <- length(genes); n <- ncol(X)
  sp <- rowMeans(X == 0)                       # BEELINE inputs are log-normalised: 0 = undetected
  B <- blom_avg(X); R <- stats::cor(t(B)); R[is.na(R)] <- 0; diag(R) <- 0
  set.seed(bopt$seed + k)
  Bp <- t(apply(B, 1, sample)); Rp <- stats::cor(t(Bp)); Rp[is.na(Rp)] <- 0; diag(Rp) <- 0
  Xp <- t(apply(X, 1, sample)); rownames(Xp) <- genes
  M <- mi_of(X); Mp <- mi_of(Xp)
  d1 <- row_null_diag(abs(R), abs(Rp), lab, sp); d1$measure <- "abs_r_blom"
  d2 <- row_null_diag(M, Mp, lab, sp); d2$measure <- "mi_bspline"
  ut <- upper.tri(Rp)
  d1$null_r2_times_nm1 <- mean(Rp[ut]^2) * (n - 1); d2$null_r2_times_nm1 <- NA
  d1$n_cells <- n; d2$n_cells <- n
  # background shift (added after the mDC smoke test showed 99% positive r): global cell-level factor
  det <- log(colSums(X > 0)); Ro <- R; diag(Ro) <- NA
  d1$row_median_signed_r_med <- stats::median(apply(Ro, 1, stats::median, na.rm = TRUE))
  d1$frac_signed_r_pos <- mean(Ro > 0, na.rm = TRUE)
  d1$med_abs_cor_gene_vs_logdetected <- stats::median(abs(stats::cor(t(B), det)))
  d2$row_median_signed_r_med <- NA; d2$frac_signed_r_pos <- NA; d2$med_abs_cor_gene_vs_logdetected <- NA
  dg[[length(dg) + 1L]] <- d1; dg[[length(dg) + 1L]] <- d2
  S <- alt_scores(R, mi = M)
  # label-free input changes that remove the global factor before CLR-NS (not calibration changes)
  Br <- t(stats::resid(stats::lm(t(B) ~ det))); Rr <- stats::cor(t(Br)); Rr[is.na(Rr)] <- 0; diag(Rr) <- 0
  S$ns_detresid <- clr::clr_calibrate(abs(Rr), method = "normal", combine = "stouffer")
  Bs <- (B - rowMeans(B)) / pmax(apply(B, 1, stats::sd), 1e-12); sv <- svd(Bs, nu = 1, nv = 1)
  B1 <- Bs - sv$d[1] * sv$u %*% t(sv$v); R1 <- stats::cor(t(B1)); R1[is.na(R1)] <- 0; diag(R1) <- 0
  S$ns_pc1 <- clr::clr_calibrate(abs(R1), method = "normal", combine = "stouffer")
  d1$frac_signed_r_pos_after_detresid <- mean(Rr[upper.tri(Rr)] > 0); d1$frac_signed_r_pos_after_pc1 <- mean(R1[upper.tri(R1)] > 0)
  d2$frac_signed_r_pos_after_detresid <- NA; d2$frac_signed_r_pos_after_pc1 <- NA
  dg[[length(dg) - 1L]] <- d1; dg[[length(dg)]] <- d2
  rm(Br, Rr, Bs, B1, R1, sv)
  ref <- clr::clr_calibrate(abs(R), method = "normal", combine = "stouffer")
  stopifnot(max(abs(ref - S$clr_ns)) < 1e-8)   # our clr_rowz + stouffer equals the package CLR-NS
  for (m in setdiff(names(S), "clr_ns")) { r <- rank_change(S[[m]], S$clr_ns, lab); r$method <- m; rc[[length(rc) + 1L]] <- r }
  cache[[lab]] <- list(S = S, genes = genes)
  say(sprintf("%s: %d genes x %d cells, diagnostics done (%.1f s)", lab, G, n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
dgd <- do.call(rbind, lapply(dg, function(d) d[, union(c("label", "measure"), names(d))]))
utils::write.csv(dgd, file.path(bopt$out, "diagnostics.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, rc), file.path(bopt$out, "rank_change.csv"), row.names = FALSE)
say("stage 1 (label-free) written")
if (bopt$stage == "1") quit(save = "no")

## ---- stage 2: truths (EXPLORATORY, in-sample) ------------------------------------------
for (k in seq_len(nrow(cases))) {
  ds <- cases$ds[k]; N <- cases$N[k]; lab <- sprintf("%s_%d", ds, N)
  S <- cache[[lab]]$S; genes <- cache[[lab]]$genes
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  for (tn in names(gts)) for (m in names(S)) {
    r <- score_truth(S[[m]], genes, NULL, gts[[tn]]); if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(dataset = ds, N = N, truth = tn, method = m, t(r))
  }
}
res <- do.call(rbind, rows); utils::write.csv(res, file.path(bopt$out, "metrics_EXPLORATORY.csv"), row.names = FALSE)
old <- utils::read.csv("results/beeline_baselines/metrics.csv", stringsAsFactors = FALSE)
chk <- merge(res[res$method %in% c("clr_ns"), ], old[old$method == "clr_ns", ], by = c("dataset", "N", "truth"))
chk2 <- merge(res[res$method == "clr_mi", ], old[old$method == "clr", ], by = c("dataset", "N", "truth"))
chk3 <- merge(res[res$method == "abs_r", ], old[old$method == "abs_pearson", ], by = c("dataset", "N", "truth"))
repro <- data.frame(check = c("CLR-NS AUPRC ratio max |diff| vs results/beeline_baselines", "CLR (MI) AUPRC ratio max |diff|", "|r| Blom vs |r| raw (expected to differ)"),
  n = c(nrow(chk), nrow(chk2), nrow(chk3)),
  max_abs_diff = c(max(abs(chk$auprc_ratio.x - chk$auprc_ratio.y)), max(abs(chk2$auprc_ratio.x - chk2$auprc_ratio.y)), max(abs(chk3$auprc_ratio.x - chk3$auprc_ratio.y))))
utils::write.csv(repro, file.path(bopt$out, "repro_check.csv"), row.names = FALSE); print(repro)
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
summ <- list()
for (tset in c("primary", "celltype_ChIP")) {
  x <- res[if (tset == "primary") res$truth %in% primary else res$truth == "celltype_ChIP", ]
  r <- x[x$method == "clr_ns", c("dataset", "N", "truth", "auprc_ratio")]; names(r)[4] <- "ap_ref"; x <- merge(x, r)
  for (m in setdiff(unique(x$method), "clr_ns")) { y <- x[x$method == m, ]
    lr <- tapply(log2(y$auprc_ratio / y$ap_ref), y$dataset, mean)
    summ[[length(summ) + 1L]] <- data.frame(truths = tset, method = m, n_cells = nrow(y), median_auprc_ratio = stats::median(y$auprc_ratio),
      mean_log2_vs_clrns = mean(lr), min_dataset = min(lr), max_dataset = max(lr), datasets_better = sum(lr > 0), n_datasets = length(lr),
      p_wilcoxon_exact = if (any(lr != 0)) stats::wilcox.test(lr, exact = TRUE)$p.value else NA) } }
summ <- do.call(rbind, summ); summ$median_clrns <- NA
for (tset in c("primary", "celltype_ChIP")) summ$median_clrns[summ$truths == tset] <- stats::median(res$auprc_ratio[res$method == "clr_ns" & (if (tset == "primary") res$truth %in% primary else res$truth == "celltype_ChIP")])
utils::write.csv(summ, file.path(bopt$out, "summary_EXPLORATORY.csv"), row.names = FALSE)
print(summ, digits = 3)
say("done")
