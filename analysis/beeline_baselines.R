#!/usr/bin/env Rscript
# Competitor baselines on the BEELINE single-cell benchmark (paper plan theme 5,
# item 4; approved 2026-10-21): GENIE3, GRNBoost2 and PIDC against CLR, CLR-NS
# (CLR on |r| of Blom normal scores, our best single-cell variant) and |r|, on
# the same 14 cases as analysis/beeline.R: 7 datasets (human hESC, hHep; mouse
# mDC, mESC, mHSC-E/GM/L) x BEELINE gene sets TFs+500 and TFs+1000, BEELINE's
# own protocol (methods fitted on the selected genes only). Same truths and
# scoring as beeline.R (its definitions are evaluated here, not copied):
# directed TF-node AUPRC ratio and EPR per truth (cell-type ChIP-seq,
# non-specific ChIP-seq, STRING, mESC LOF/GOF), plus regulon co-membership.
# GENIE3 and GRNBoost2 are scored on their directed regulator -> target
# weights (every selected gene a candidate regulator, as in analysis/baselines.R);
# PIDC and the correlation/CLR scores are symmetric.
#
# Stages (PIDC's NetworkInference.jl on Julia 1.6 crashes when launched from R,
# so it runs in between, from the job script):
#   --stage export   write <out>/<case>/expression.tsv (genes x cells) per case
#   (job script)     julia +1.6 tools/pidc.jl <case>/expression.tsv <case>/pidc_edges.tsv
#   --stage score    GENIE3, GRNBoost2, read PIDC, score everything
# Usage (repo root):
#   Rscript analysis/beeline_baselines.R --stage export|score [--out results/beeline_baselines]
#       [--threads 16] [--ntrees 1000]

args0 <- commandArgs(trailingOnly = TRUE)
bopt <- list(stage = "score", out = "results/beeline_baselines", threads = 16L, ntrees = 1000L, seed = 20261021L,
             python = "~/venvs/clr-baselines/bin/python")
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i])
  bopt[[k]] <- if (k %in% c("threads", "ntrees", "seed")) as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }

# beeline.R's setup and helpers (everything before its main loop), with our out/threads
src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]
args <- c("--out", bopt$out, "--threads", as.character(bopt$threads))
for (k in intersect(c("datasets", "N"), names(bopt))) args <- c(args, paste0("--", k), bopt[[k]])   # passed through
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  if (identical(deparse(e)[1], "args <- commandArgs(trailingOnly = TRUE)")) next
  eval(e, globalenv())
}
logf <- file.path(bopt$out, "beeline_baselines.log")
# Tie-aware early precision (review 2026-10-21): beeline.R's score_truth breaks ties in the
# top-k by candidate index, which matters for sparse scores (GRNBoost2's zeros) and for
# symmetric scores (both directions of a TF-TF pair tie). Here precision at k is the
# expectation over random tie-breaking: (tp_above + (k - n_above) * tp_tie / n_tie) / k.
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
say("beeline_baselines: stage ", bopt$stage, ", out ", bopt$out, ", threads ", bopt$threads, ", ntrees ", bopt$ntrees)

Ns_all <- as.integer(strsplit(opt$N, ",")[[1]]); dss <- strsplit(opt$datasets, ",")[[1]]
cases <- expand.grid(N = Ns_all, ds = dss, stringsAsFactors = FALSE)
case_dir <- function(ds, N) file.path(bopt$out, paste0(ds, "_", N))
case_X <- function(ds, N) {     # as beeline.R with --context none
  ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1, check.names = FALSE))
  g <- intersect(select_genes(ds, N)$genes, rownames(ex)); X <- ex[g, , drop = FALSE]
  X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]
}

if (bopt$stage == "export") {
  for (k in seq_len(nrow(cases))) {
    ds <- cases$ds[k]; N <- cases$N[k]; d <- case_dir(ds, N); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    X <- case_X(ds, N)
    utils::write.table(data.frame(gene = rownames(X), X, check.names = FALSE), file.path(d, "expression.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    say(sprintf("exported %s TFs+%d: %d genes x %d cells", ds, N, nrow(X), ncol(X)))
  }
  quit(save = "no")
}

## ---- score -----------------------------------------------------------------------
edges_to <- function(a, b, w, genes, sym) { M <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  ia <- match(a, genes); ib <- match(b, genes); ok <- !is.na(ia) & !is.na(ib) & ia != ib
  M[cbind(ia[ok], ib[ok])] <- w[ok]; if (sym) pmax(M, t(M)) else M }
rows <- list(); crow <- list()
for (k in seq_len(nrow(cases))) {
  ds <- cases$ds[k]; N <- cases$N[k]; d <- case_dir(ds, N); X <- case_X(ds, N); genes <- rownames(X)
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  S <- list()
  S$abs_pearson <- abs(stats::cor(t(X)))
  f <- clr_fit(X); S$clr <- f$clr_scores; f$release(FALSE)
  S$clr_ns <- ns_clr(X)
  t0 <- Sys.time(); set.seed(bopt$seed)
  W <- tryCatch(GENIE3::GENIE3(X, nTrees = bopt$ntrees, K = "sqrt", nCores = bopt$threads, verbose = FALSE),
                error = function(e) { say("GENIE3 failed on ", ds, " ", N, ": ", conditionMessage(e)); NULL })
  if (!is.null(W)) S$genie3 <- W[genes, genes]           # directed: rows regulators, columns targets
  say(sprintf("%s TFs+%d: GENIE3 %.1f min", ds, N, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  gf <- file.path(d, "grnboost2_edges.tsv")
  if (!file.exists(gf)) {
    st <- system2(path.expand(bopt$python), c("tools/grnboost2.py", file.path(d, "expression.tsv"), paste0(gf, ".part"),
                  "--workers", bopt$threads, "--seed", bopt$seed), stdout = file.path(d, "grnboost2.log"), stderr = file.path(d, "grnboost2.log"))
    if (st == 0 && file.exists(paste0(gf, ".part"))) file.rename(paste0(gf, ".part"), gf) else say("GRNBoost2 failed on ", ds, " ", N, " (see ", d, "/grnboost2.log)")
  }
  if (file.exists(gf)) { e <- utils::read.delim(gf, stringsAsFactors = FALSE); S$grnboost2 <- edges_to(e$TF, e$target, e$importance, genes, FALSE)
  } else say("GRNBoost2 output missing for ", ds, " ", N)
  pf <- file.path(d, "pidc_edges.tsv")
  if (file.exists(pf)) { e <- utils::read.delim(pf, header = FALSE, stringsAsFactors = FALSE); S$pidc <- edges_to(e[[1]], e[[2]], e[[3]], genes, TRUE)
  } else say("PIDC edges missing for ", ds, " ", N)
  for (m in names(S)) { s <- S[[m]]; s[!is.finite(s)] <- 0; diag(s) <- 0; S[[m]] <- s }
  for (tn in names(gts)) for (m in names(S)) {
    r <- score_truth(S[[m]], genes, NULL, gts[[tn]]); if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(dataset = ds, N = N, truth = tn, method = m, t(r))
  }
  for (m in names(S)) { r <- comembership_truth(S[[m]], genes, gts$celltype_ChIP); if (is.null(r)) next
    crow[[length(crow) + 1L]] <- data.frame(dataset = ds, N = N, method = m, t(r)) }
  cur <- do.call(rbind, rows); cur <- cur[cur$dataset == ds & cur$N == N & cur$truth == "celltype_ChIP", ]
  for (j in seq_len(nrow(cur))) say(sprintf("    %s TFs+%d %-12s AUPRC ratio %.3f  EPR %.3f (cell-type ChIP)", ds, N,
                                            cur$method[j], cur$auprc_ratio[j], cur$epr[j]))
  rm(S, W); invisible(gc())
}
res <- do.call(rbind, rows); utils::write.csv(res, file.path(bopt$out, "metrics.csv"), row.names = FALSE)
if (length(crow)) utils::write.csv(do.call(rbind, crow), file.path(bopt$out, "comembership.csv"), row.names = FALSE)

# Summaries (review 2026-10-21). Primary truths are BEELINE's own (cell-type ChIP,
# non-specific ChIP, STRING, mESC LOF/GOF); the derived STRING variants are secondary.
# Only cells where every method has a result. References: CLR (fixed in advance) and
# CLR-NS (chosen on these truths earlier, beeline_blom_test: disclose). Inference is
# blocked by dataset: per dataset, the mean of log2(method / reference) over gene sets
# and primary truths, then an exact Wilcoxon signed-rank test and a sign test on the 7.
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
cells <- unique(res[, c("dataset", "N", "truth")])
complete <- merge(cells, stats::aggregate(method ~ dataset + N + truth, res, length))
complete <- complete[complete$method == length(unique(res$method)), c("dataset", "N", "truth")]
rc <- merge(res, complete)
summ <- list()
for (refm in c("clr", "clr_ns")) for (tset in c("primary", "celltype_ChIP")) {
  x <- rc[if (tset == "primary") rc$truth %in% primary else rc$truth == "celltype_ChIP", ]
  r <- x[x$method == refm, c("dataset", "N", "truth", "auprc_ratio", "epr")]; names(r)[4:5] <- c("ap_ref", "epr_ref")
  x <- merge(x, r)
  for (m in setdiff(unique(x$method), refm)) { y <- x[x$method == m, ]
    lr <- tapply(log2(y$auprc_ratio / y$ap_ref), y$dataset, mean); le <- tapply(log2(pmax(y$epr, 1e-9) / pmax(y$epr_ref, 1e-9)), y$dataset, mean)
    summ[[length(summ) + 1L]] <- data.frame(reference = refm, truths = tset, method = m, n_cells = nrow(y), n_datasets = length(lr),
      median_auprc_ratio = stats::median(y$auprc_ratio), median_epr = stats::median(y$epr),
      mean_log2_auprc_vs_ref = mean(lr), datasets_better_auprc = sum(lr > 0),
      p_wilcoxon_datasets = if (length(lr) >= 3 && any(lr != 0)) stats::wilcox.test(lr, exact = TRUE)$p.value else NA,
      p_sign_datasets = if (length(lr) >= 1) stats::binom.test(sum(lr > 0), sum(lr != 0))$p.value else NA,
      mean_log2_epr_vs_ref = mean(le), datasets_better_epr = sum(le > 0)) } }
summ <- do.call(rbind, summ); utils::write.csv(summ, file.path(bopt$out, "summary.csv"), row.names = FALSE)
say("blocked by dataset (7 datasets; per dataset the mean log2 ratio over gene sets and truths):")
for (j in seq_len(nrow(summ))) with(summ[j, ], say(sprintf(
  "  vs %-6s %-13s %-12s median AUPRC ratio %.3f  log2 vs ref %+.3f  better in %d/%d datasets  (exact Wilcoxon p = %.3g, sign p = %.3g)  EPR log2 %+.3f",
  reference, truths, method, median_auprc_ratio, mean_log2_auprc_vs_ref, datasets_better_auprc, n_datasets,
  p_wilcoxon_datasets, p_sign_datasets, mean_log2_epr_vs_ref)))
say("done")
