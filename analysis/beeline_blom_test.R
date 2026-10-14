#!/usr/bin/env Rscript
# Small test: does a Blom transform with randomly broken ties rescue CLR on
# zero-inflated single-cell data, and does CLR then need MI at all?
#
# scRNA-seq counts are dominated by zeros. On raw values the B-spline MI
# estimator spends its bins on that point mass. Blom scores per gene,
#     b = qnorm((rank - 3/8) / (n + 1/4)), ties broken at random,
# give every gene an exactly Gaussian marginal; the zeros become a randomly
# ordered lower tail, so zero-zero cell pairs carry no ordering signal.
# On Blom data, correlation is the natural dependence measure, and its
# Gaussian MI, -1/2 log(1 - r^2), is a monotone function of |r| -- so the
# assertion to test is that CLR on correlation of Blom data matches or beats
# CLR on B-spline MI.
#
# Methods per case (BEELINE gene sets, no wider context):
#   abs_pearson_raw     |r| on raw values
#   clr_mi_raw          CLR on B-spline MI of raw values (the current method)
#   abs_pearson_blom    |r| on Blom values
#   mi_blom             B-spline MI on Blom values
#   clr_mi_blom         CLR on B-spline MI of Blom values
#   clr_gmi_blom        CLR on Gaussian MI -1/2 log(1 - r^2) of Blom values
#   clr_absr_blom       CLR on |r| of Blom values
#   att_<clr variant>   soft10 attention mass on each CLR variant, at the
#                       held-out t* of results/beeline_v2 for that case
# Blom-based methods are run with `seeds` tie-breaking seeds; the summary
# reports the mean over seeds and the seed-to-seed range. Two tie-handling
# alternatives (added after the first run showed large seed-to-seed swings
# on truths with few edges):
#   *_blomavg   Blom scores of average ranks (ties kept: every zero gets the
#               same score), deterministic
#   *_blomens   CLR averaged over `nens` random tie-breaks (the attention
#               operator is built on the averaged CLR)
#
# Usage (repo root):
#   Rscript analysis/beeline_blom_test.R [--seeds 3] [--N 500,1000] [--workers 6]
#       [--out results/beeline_blom_test]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(beeline = "data/beeline", N = "500,1000",
            datasets = "hESC,hHep,mDC,mESC,mHSC-E,mHSC-GM,mHSC-L",
            out = "results/beeline_blom_test", seeds = 3L, workers = 6L,
            threads = NULL, alpha = 0.5, mem_gb = 1.5,
            depths = "results/beeline_v2/heldout_depths.csv", nens = 10L)
ints <- c("seeds", "workers", "threads", "nens")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key %in% ints) as.integer(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "blom_test.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))       # ap_fast

## ---- data (as beeline_boot.R) ---------------------------------------------------
root <- opt$beeline
expr_dir <- file.path(root, "BEELINE-data", "inputs", "scRNA-Seq")
net_dir <- file.path(root, "Networks")
species <- function(ds) if (ds %in% c("hESC", "hHep")) "human" else "mouse"
tf_file <- function(sp) {
  cand <- c(file.path(root, paste0(sp, "-tfs.csv")), file.path(net_dir, paste0(sp, "-tfs.csv")))
  cand[file.exists(cand)][1]
}
truth_files <- function(ds) {
  sp <- species(ds); d <- file.path(net_dir, sp)
  ct <- switch(ds, hESC = "hESC-ChIP-seq-network.csv", hHep = "HepG2-ChIP-seq-network.csv",
               mDC = "mDC-ChIP-seq-network.csv", mESC = "mESC-ChIP-seq-network.csv",
               "mHSC-ChIP-seq-network.csv")
  ns <- if (sp == "human") "Non-specific-ChIP-seq-network.csv" else "Non-Specific-ChIP-seq-network.csv"
  out <- list(celltype_ChIP = file.path(d, ct), nonspecific_ChIP = file.path(d, ns),
              STRING = file.path(d, "STRING-network.csv"))
  if (ds == "mESC") out$LOFGOF <- file.path(d, "mESC-lofgof-network.csv")
  dd <- file.path(root, "Networks_derived", sp)
  der <- c(STRING_regulatory = "STRING-regulatory.csv", STRING_curated = "STRING-curated.csv")
  for (k in names(der)) if (file.exists(file.path(dd, der[[k]]))) out[[k]] <- file.path(dd, der[[k]])
  out
}
select_genes <- function(ds, N) {
  go <- utils::read.csv(file.path(expr_dir, ds, "GeneOrdering.csv"), stringsAsFactors = FALSE,
                        check.names = FALSE)
  names(go)[1] <- "gene"
  sig <- go[go$VGAMpValue * nrow(go) < 0.01, ]
  sig <- sig[order(-sig$Variance), ]
  tfs <- toupper(utils::read.csv(tf_file(species(ds)), stringsAsFactors = FALSE)$TF)
  union(sig$gene[toupper(sig$gene) %in% tfs], utils::head(sig$gene, N))
}
edge_sets <- function(genes, gt) {
  gi <- match(toupper(gt$Gene1), toupper(genes)); gj <- match(toupper(gt$Gene2), toupper(genes))
  ok <- !is.na(gi) & !is.na(gj) & gi != gj
  regs <- sort(unique(gi[ok]))
  if (length(regs) < 1 || sum(ok) < 10) return(NULL)
  G <- length(genes)
  cand <- cbind(rep(regs, each = G), rep(seq_len(G), times = length(regs)))
  cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
  lab <- as.integer((cand[, 1] * (G + 1) + cand[, 2]) %in% unique(gi[ok] * (G + 1) + gj[ok]))
  list(cand = cand, lab = lab)
}
dep <- utils::read.csv(opt$depths, stringsAsFactors = FALSE)
cases <- list()
for (ds in strsplit(opt$datasets, ",")[[1]]) {
  ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1,
                                  check.names = FALSE))
  gts <- lapply(truth_files(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  for (N in as.integer(strsplit(opt$N, ",")[[1]])) {
    X <- ex[intersect(select_genes(ds, N), rownames(ex)), , drop = FALSE]
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]
    es <- Filter(Negate(is.null), lapply(gts, function(gt) edge_sets(rownames(X), gt)))
    t_soft <- dep$t_soft10[dep$dataset == ds & dep$N == N]
    cases[[sprintf("%s_%d", ds, N)]] <- list(ds = ds, N = N, X = X, es = es,
                                             t = if (length(t_soft)) t_soft else 8L,
                                             zero_frac = mean(X == 0))
  }
  rm(ex, gts)
}
say(sprintf("%d cases; zero fraction %s", length(cases),
            paste(sprintf("%s %.0f%%", names(cases), 100 * vapply(cases, `[[`, 0, "zero_frac")),
                  collapse = ", ")))

## ---- constructions -----------------------------------------------------------------
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
blom <- function(X, ties = "random") {
  n <- ncol(X)
  t(apply(X, 1, function(x) stats::qnorm((rank(x, ties.method = ties) - 3 / 8) / (n + 1 / 4))))
}
abs_cor <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; diag(r) <- 0; r }
bspline <- function(X) { m <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none",
                                                              threads = nthr())$mi
                         diag(m) <- 0; m }
calib <- function(M) { diag(M) <- 0; clr_calibrate(M, method = "normal", combine = "stouffer") }
attention <- function(S, t) {
  A <- clr:::.softmax_rows(S, keff = 10, cap = 50L)
  P <- Matrix::Diagonal(nrow(A), 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- as.matrix(P)
  if (t > 1) for (s in 2:t) Pt <- as.matrix(Pt %*% P)
  M <- (Pt + t(Pt)) / 2; diag(M) <- 0; M
}
score_all <- function(S, es) do.call(rbind, lapply(names(es), function(tn) {
  e <- es[[tn]]; s <- S[e$cand]; s[is.na(s)] <- 0
  dens <- mean(e$lab); k <- sum(e$lab); o <- order(-s)
  data.frame(truth = tn, auprc_ratio = ap_fast(s, e$lab) / dens,
             epr = mean(e$lab[o[seq_len(k)]]) / dens)
}))

run <- function(cn, seed) {
  cs <- cases[[cn]]; out <- list()
  add <- function(S, m) out[[length(out) + 1L]] <<- cbind(case = cn, seed = seed, method = m,
                                                          score_all(S, cs$es))
  if (seed == 1L) {                     # tie-breaking does not touch the raw methods
    add(abs_cor(cs$X), "abs_pearson_raw")
    C <- calib(bspline(cs$X)); add(C, "clr_mi_raw"); add(attention(C, cs$t), "att_clr_mi_raw")
    Ba <- blom(cs$X, "average")
    C <- calib(abs_cor(Ba)); add(C, "clr_absr_blomavg"); add(attention(C, cs$t), "att_clr_absr_blomavg")
    C <- calib(bspline(Ba)); add(C, "clr_mi_blomavg"); add(attention(C, cs$t), "att_clr_mi_blomavg")
    set.seed(777L + match(cn, names(cases)))
    Ce <- Reduce(`+`, lapply(seq_len(opt$nens), function(e) calib(abs_cor(blom(cs$X))))) / opt$nens
    add(Ce, "clr_absr_blomens"); add(attention(Ce, cs$t), "att_clr_absr_blomens")
  }
  set.seed(1000L * seed + match(cn, names(cases)))
  B <- blom(cs$X)
  R <- abs_cor(B); add(R, "abs_pearson_blom")
  M <- bspline(B); add(M, "mi_blom")
  C <- calib(M); add(C, "clr_mi_blom"); add(attention(C, cs$t), "att_clr_mi_blom")
  G <- -0.5 * log(pmax(1 - R^2, 1e-12)); diag(G) <- 0
  C <- calib(G); add(C, "clr_gmi_blom"); add(attention(C, cs$t), "att_clr_gmi_blom")
  C <- calib(R); add(C, "clr_absr_blom"); add(attention(C, cs$t), "att_clr_absr_blom")
  do.call(rbind, out)
}
grid <- expand.grid(seed = seq_len(opt$seeds), case = names(cases), stringsAsFactors = FALSE)
source(file.path(script_dir, "parallel.R"))
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
res <- foreach(k = seq_len(nrow(grid)), .errorhandling = "stop") %dopar%
  run(grid$case[k], grid$seed[k])
par_stop(pc)
res <- do.call(rbind, res)
utils::write.csv(res, file.path(opt$out, "metrics.csv"), row.names = FALSE)

## ---- summary -------------------------------------------------------------------------
# mean over seeds per case x truth x method; seed range of the case-level values
cm <- stats::aggregate(cbind(auprc_ratio, epr) ~ case + truth + method, data = res, FUN = mean)
rg <- stats::aggregate(auprc_ratio ~ case + truth + method, data = res,
                       FUN = function(v) diff(range(v)))
names(rg)[4] <- "seed_range"
cm <- merge(cm, rg)
ref <- cm[cm$method == "clr_mi_raw", c("case", "truth", "auprc_ratio", "epr")]
names(ref)[3:4] <- c("ap_ref", "epr_ref")
m <- merge(cm, ref)
summ <- do.call(rbind, lapply(split(m, m$method), function(x) data.frame(
  method = x$method[1], n = nrow(x),
  median_auprc_ratio = stats::median(x$auprc_ratio), median_epr = stats::median(x$epr),
  median_d_ap_vs_clr_mi_raw = stats::median(x$auprc_ratio - x$ap_ref),
  wins_vs_clr_mi_raw = mean(x$auprc_ratio > x$ap_ref),
  wilcoxon_p = if (all(x$auprc_ratio == x$ap_ref)) NA_real_ else
    suppressWarnings(stats::wilcox.test(x$auprc_ratio, x$ap_ref, paired = TRUE)$p.value),
  max_seed_range = max(x$seed_range))))
summ <- summ[order(-summ$median_auprc_ratio), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
say("all case x truth combinations (mean over tie-breaking seeds), vs CLR on raw MI:")
for (k in seq_len(nrow(summ)))
  say(sprintf("  %-20s median AUPRC ratio %.3f  EPR %.3f | vs clr_mi_raw %+.3f, wins %.0f%%, p = %.2g | max seed range %.3f",
              summ$method[k], summ$median_auprc_ratio[k], summ$median_epr[k],
              summ$median_d_ap_vs_clr_mi_raw[k], 100 * summ$wins_vs_clr_mi_raw[k],
              summ$wilcoxon_p[k], summ$max_seed_range[k]))
bt <- stats::aggregate(auprc_ratio ~ truth + method, data = cm, FUN = stats::median)
bt <- stats::xtabs(auprc_ratio ~ truth + method, data = bt)
utils::write.csv(as.data.frame.matrix(bt), file.path(opt$out, "by_truth.csv"))
say("median AUPRC ratio by truth:")
for (l in utils::capture.output(print(round(bt, 3)))) say("  ", l)
say("done")
