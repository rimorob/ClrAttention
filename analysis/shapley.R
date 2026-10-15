#!/usr/bin/env Rscript
# Which components of completed CLR attention carry the gain? Exact Shapley
# values over six on/off components, on the 466-experiment M3D compendium
# against RegulonDB (co-membership AUPR and median regulon coherence).
#
# Components (on = the most performant configuration; off = plain alternative):
#   Dependence: mutual information   B-spline MI           | |Pearson r|
#   Binning: common HG bins          one joint-histogram count (HG rule) | per-gene Freedman-Diaconis
#                                    (only matters when MI is on)
#   Context normalization (CLR)      row z-scores, Stouffer | raw dependence
#   Softmax with temperature         softmax, keff 10, cap 50 | hard top-10 neighbours, weights
#                                    proportional to scores (only matters with depth)
#   Depth (iterated attention)       attention mass at t = 9 | scores ranked directly
#   Removal of 2 global PCs          top 2 principal components removed | none
#                                    (chosen after inspecting Supplementary Note 7; flagged as such)
# v(S) is the metric with the components in S on. phi_i = sum over S not
# containing i of |S|!(n-|S|-1)!/n! [v(S + i) - v(S)]; the phi sum to
# v(all) - v(none). Drop-one losses v(all) - v(all - i) are reported too.
# Intervals: 100 paired delete-half gene jackknife replicates (AUPR) and 1,000
# paired regulon bootstrap replicates (coherence), shared by all 64 subsets,
# so every Shapley value is recomputed on each replicate.
#
# Usage (repo root):  Rscript analysis/shapley.R [--evidence SC] [--workers 6] [--out results/shapley]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(evidence = "SC", workers = 6L, mem_gb = 6, out = "results/shapley", t = 9L,
            jack = 100L, boot = 1000L, seed = 20261015L, alpha = 0.5, threads = NULL)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% c("workers", "t", "jack", "boot", "seed", "threads")) as.integer(v) else
    if (key %in% c("mem_gb", "alpha")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "shapley.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R")
set.seed(opt$seed)

m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                         check.names = FALSE, stringsAsFactors = FALSE)
X0 <- as.matrix(raw[, -1]); rownames(X0) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X0 <- X0[apply(X0, 1, function(x) diff(range(x))) > 0, ]
genes <- rownames(X0); G <- length(genes)
remove_pcs <- function(Y, k) { Yc <- Y - rowMeans(Y); V <- svd(Yc, nu = 0, nv = k)$v; Y - Yc %*% V %*% t(V) }
Xs <- list(k0 = X0, k2 = remove_pcs(X0, 2L))
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
ck <- if (opt$evidence == "SC") c("C", "S") else c("C", "S", "W", "?")
reg <- make_regulons(rtab, genes, conf_keep = ck)
pairs <- comembership_pairs(reg, operon_of, genes); co <- comember_matrix(reg, G)
halves <- lapply(seq_len(opt$jack), function(r) {          # shared gene half-samples
  keep <- sample(pairs$genes_in, floor(length(pairs$genes_in) / 2))
  inh <- logical(G); inh[keep] <- TRUE; which(inh[pairs$i] & inh[pairs$j]) })

## ---- the 36 distinct score matrices, grouped by base dependence matrix ----------
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
calib <- function(M) { diag(M) <- 0; clr_calibrate(M, method = "normal", combine = "stouffer") }
walk_t <- function(A) {
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- as.matrix(P); for (s in seq_len(opt$t - 1L)) Pt <- as.matrix(Pt %*% P)
  S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S
}
hard_top10 <- function(S) { A <- clr:::.topk_rows(S, 10L); rs <- rowSums(A)
  A[rs > 0, ] <- A[rs > 0, ] / rs[rs > 0]; e <- which(rs == 0); A[cbind(e, e)] <- 1; A }
soft <- function(S) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L); rs <- rowSums(A)
  e <- which(rs == 0); A[cbind(e, e)] <- 1; A }
evaluate <- function(S) {
  diag(S) <- 0; sc <- S[cbind(pairs$i, pairs$j)]
  list(aupr = ap_fast(sc, pairs$label),
       jack = vapply(halves, function(ix) ap_fast(sc[ix], pairs$label[ix]), 0),
       coh = regulon_coherence(S, reg, co, operon_of, genes))
}
groups <- expand.grid(dep = c("mi_hg", "mi_fd", "r"), pc = c("k0", "k2"), stringsAsFactors = FALSE)
source("analysis/parallel.R")
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
ev <- foreach(g = seq_len(nrow(groups)), .errorhandling = "stop") %dopar% {
  X <- Xs[[groups$pc[g]]]; dep <- groups$dep[g]
  D <- if (dep == "r") { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; r } else
    ClrAttention$new(X)$estimate_mi(bins = if (dep == "mi_hg") "hg" else "fd", transform = "none",
                                    threads = nthr())$mi
  diag(D) <- 0
  out <- list()
  for (clr_on in c(FALSE, TRUE)) {
    S <- if (clr_on) calib(D) else D
    key <- function(att) sprintf("%s|%s|clr%d|%s", dep, groups$pc[g], clr_on, att)
    out[[key("none")]] <- evaluate(S)
    out[[key("soft")]] <- evaluate(walk_t(soft(S)))
    out[[key("hard")]] <- evaluate(walk_t(hard_top10(S)))
    par_release()
  }
  out
}
par_stop(pc)
ev <- do.call(c, ev)
say(sprintf("%d distinct configurations evaluated", length(ev)))

## ---- Shapley values -----------------------------------------------------------------
comps <- c("Dependence: mutual information", "Binning: common HG bins", "Context normalization (CLR)",
           "Softmax with temperature", "Depth (iterated attention)", "Removal of 2 global PCs")
n <- length(comps)
subsets <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n))); colnames(subsets) <- comps
key_of <- function(on) {
  dep <- if (!on[1]) "r" else if (on[2]) "mi_hg" else "mi_fd"
  att <- if (!on[5]) "none" else if (on[4]) "soft" else "hard"
  sprintf("%s|%s|clr%d|%s", dep, if (on[6]) "k2" else "k0", on[3], att)
}
keys <- apply(subsets, 1, key_of)
ok_reg <- Reduce(`&`, lapply(ev, function(e) !is.na(e$coh)))
bidx <- replicate(opt$boot, sample(sum(ok_reg), replace = TRUE))
val <- function(metric) {        # rows: subsets; cols: point estimate then replicates
  t(vapply(keys, function(k) { e <- ev[[k]]
    if (metric == "aupr") c(e$aupr, e$jack) else { cc <- e$coh[ok_reg]
      c(stats::median(cc), apply(bidx, 2, function(b) stats::median(cc[b]))) } },
    numeric(1 + if (metric == "aupr") opt$jack else opt$boot)))
}
shap <- function(V) {           # V: 2^n x R values -> n x R Shapley values
  sz <- rowSums(subsets)
  t(vapply(seq_len(n), function(j) {
    phi <- 0
    for (s in which(!subsets[, j])) {
      on <- subsets[s, ]; on[j] <- TRUE
      s2 <- which(apply(subsets, 1, function(r) all(r == on)))
      w <- factorial(sz[s]) * factorial(n - sz[s] - 1) / factorial(n)
      phi <- phi + w * (V[s2, ] - V[s, ])
    }
    phi }, numeric(ncol(V))))
}
full <- which(rowSums(subsets) == n); none <- which(rowSums(subsets) == 0)
res <- list(); cfg <- list()
for (metric in c("aupr", "coherence")) {
  V <- val(if (metric == "aupr") "aupr" else "coh")
  Phi <- shap(V)
  rep_sd <- if (metric == "aupr") function(x) stats::sd(x) else NULL
  ci <- function(x) {          # jackknife: point +/- 1.96 sd (delete-half); bootstrap: percentiles
    if (metric == "aupr") c(x[1] - 1.96 * stats::sd(x[-1]), x[1] + 1.96 * stats::sd(x[-1])) else
      stats::quantile(x[-1], c(0.025, 0.975), names = FALSE) }
  drop <- t(vapply(seq_len(n), function(j) { on <- rep(TRUE, n); on[j] <- FALSE
    s <- which(apply(subsets, 1, function(r) all(r == on))); V[full, ] - V[s, ] }, numeric(ncol(V))))
  for (j in seq_len(n)) {
    res[[length(res) + 1L]] <- data.frame(metric = metric, component = comps[j],
      shapley = Phi[j, 1], shapley_lo = ci(Phi[j, ])[1], shapley_hi = ci(Phi[j, ])[2],
      drop_one = drop[j, 1], drop_lo = ci(drop[j, ])[1], drop_hi = ci(drop[j, ])[2],
      v_full = V[full, 1], v_none = V[none, 1])
  }
  cfg[[metric]] <- data.frame(subsets, metric = metric, value = V[, 1], check.names = FALSE)
}
res <- do.call(rbind, res)
utils::write.csv(res, file.path(opt$out, sprintf("shapley_%s.csv", opt$evidence)), row.names = FALSE)
utils::write.csv(do.call(rbind, cfg), file.path(opt$out, sprintf("configs_%s.csv", opt$evidence)), row.names = FALSE)
for (metric in c("aupr", "coherence")) {
  x <- res[res$metric == metric, ]
  say(sprintf("[%s %s] full %.4f, none %.4f (difference %.4f = sum of Shapley %.4f)", opt$evidence, metric,
              x$v_full[1], x$v_none[1], x$v_full[1] - x$v_none[1], sum(x$shapley)))
  for (k in order(-x$shapley))
    say(sprintf("   %-32s Shapley %+.4f [%+.4f, %+.4f]   drop-one loss %+.4f [%+.4f, %+.4f]", x$component[k],
                x$shapley[k], x$shapley_lo[k], x$shapley_hi[k], x$drop_one[k], x$drop_lo[k], x$drop_hi[k]))
}
say("done")
