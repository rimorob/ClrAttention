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
  if (!file.exists(gf)) system2(path.expand(bopt$python), c("tools/grnboost2.py", file.path(d, "expression.tsv"), gf,
                                 "--workers", bopt$threads, "--seed", bopt$seed), stdout = FALSE, stderr = FALSE)
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
ref <- res[res$method == "clr_ns", c("dataset", "N", "truth", "auprc_ratio", "epr")]; names(ref)[4:5] <- c("ap_ns", "epr_ns")
w <- merge(res, ref)
say("summary over datasets x gene sets x truths (paired with CLR-NS):")
for (m in unique(w$method)) { x <- w[w$method == m, ]
  p <- if (m != "clr_ns" && sum(x$auprc_ratio != x$ap_ns) >= 5) stats::wilcox.test(x$auprc_ratio, x$ap_ns, paired = TRUE, exact = FALSE)$p.value else NA
  say(sprintf("  %-12s median AUPRC ratio %.3f  median EPR %.3f  beats CLR-NS in %3.0f%%  (Wilcoxon vs CLR-NS p = %.3g, n = %d)",
              m, stats::median(x$auprc_ratio), stats::median(x$epr), 100 * mean(x$auprc_ratio > x$ap_ns), p, nrow(x))) }
say("done")
