#!/usr/bin/env Rscript
# Human perturbation compendium: Perturb-seq pseudobulk profiles as
# "experiments" (Replogle et al. 2022). Co-expression network inference
# assumes the samples explore many network states; unperturbed cells of one
# condition do not (analysis/perturbseq.R is that contrast). Here each sample
# is one knockdown's mean expression change against control cells (z-scored
# per gene; the authors' normalized pseudobulk), the human analogue of the
# M3D experiment compendium.
#
# Design:
#   - knockdowns with a significant transcriptional phenotype (energy test
#     p < --p_max) are split at random into halves A and B;
#   - networks are built on A only;
#   - truths: DoRothEA A-B and A-C, CollecTRI, CORUM (literature), and
#     knockdown-response modules built from B only (no shared samples);
#   - methods: |r|, MI, CLR on each, CLR on |r| of normal scores (CLR-NS), and
#     soft10 attention on every CLR variant at a depth chosen on random halves
#     of A (criterion: top 1% of the same CLR matrix on the other half);
#   - intervals: paired delete-half gene jackknife (AUPR) for attention minus
#     its CLR and for CLR minus |r|.
#
# Usage (repo root):
#   Rscript analysis/perturbseq_bulk.R --h5 data/perturbseq/K562_essential_normalized_bulk_01.h5ad
#       [--out results/perturbseq_bulk_k562ess] [--p_max 1e-3] [--workers 6] [--tmax 32]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(h5 = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad", out = "results/perturbseq_bulk",
            p_max = 1e-3, workers = 6L, mem_gb = 12, tmax = 32L, seed = 20261015L, alpha = 0.5,
            threads = NULL, jack = 100L)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% c("workers", "tmax", "seed", "threads", "jack")) as.integer(v) else
    if (key %in% c("p_max", "mem_gb", "alpha")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "perturbseq_bulk.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
set.seed(opt$seed)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

## ---- compendium -----------------------------------------------------------------------
rd <- function(p) rhdf5::h5read(opt$h5, p)
cat_of <- function(grp, field) {
  x <- rd(paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L]
}
Z <- rd("X"); genes <- cat_of("/var", "gene_name"); rownames(Z) <- genes
pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
pval <- as.numeric(rd("/obs/energy_test_p_value"))
ok <- pval < opt$p_max & !grepl("non-targeting", pert)
K <- pert[ok]; K <- K[sample(length(K))]
A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
X <- Z[, A]; rm(Z); invisible(gc())
# Missing values (present in the genome-wide K562 and RPE1 files): drop genes
# missing in more than 1% of the profiles, set the rest to 0 (no change from
# control cells).
nf <- !is.finite(X); bad_gene <- rowMeans(nf) > 0.01
X <- X[!bad_gene, , drop = FALSE]; n_imp <- sum(!is.finite(X)); X[!is.finite(X)] <- 0
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; genes <- rownames(X); G <- length(genes)
say(sprintf("non-finite values: %d genes dropped (> 1%% missing), %d values set to 0", sum(bad_gene), n_imp))
say(sprintf("%d profiles; %d with a significant phenotype (p < %g): %d for networks (A), %d for response modules (B); %d genes",
            length(pert), length(K), opt$p_max, length(A), length(B), G))

## ---- truths -----------------------------------------------------------------------------
truths <- list(DoRothEA_AB = dorothea_sets(genes), DoRothEA_ABC = dorothea_sets(genes, c("A", "B", "C")),
               CollecTRI = collectri_sets(genes), CORUM = corum_sets(genes),
               KD_response = response_sets(genes, path = opt$h5, use = B))
no_operon <- stats::setNames(character(0), character(0))
bench <- lapply(truths, function(r) list(reg = r, pairs = comembership_pairs(r, no_operon, genes),
                                         co = comember_matrix(r, G)))
for (nm in names(bench)) say(sprintf("truth %-12s %4d sets; %d scored pairs, base rate %.4f", nm,
  length(bench[[nm]]$reg$members), length(bench[[nm]]$pairs$label), mean(bench[[nm]]$pairs$label)))
halves <- lapply(bench, function(b) lapply(seq_len(opt$jack), function(r) {
  keep <- sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2))
  inh <- logical(G); inh[keep] <- TRUE; which(inh[b$pairs$i] & inh[b$pairs$j]) }))

## ---- methods ------------------------------------------------------------------------------
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
blom <- function(M) { n <- ncol(M)
  t(apply(M, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4)))) }
abs_r <- function(M) { r <- abs(stats::cor(t(M))); r[is.na(r)] <- 0; diag(r) <- 0; r }
calib <- function(M) { diag(M) <- 0; clr_calibrate(M, method = "normal", combine = "stouffer") }
mi_of <- function(M) { m <- ClrAttention$new(M)$estimate_mi(bins = "hg", transform = "none", threads = nthr())$mi
                       diag(m) <- 0; m }
base_of <- function(M, kind) switch(kind, absr = abs_r(M), clr_absr = calib(abs_r(M)), mi = mi_of(M),
                                    clr_mi = calib(mi_of(M)), ns_absr = abs_r(blom(M)),
                                    clr_ns = calib(abs_r(blom(M))))
soft <- function(S) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L); e <- which(rowSums(A) == 0)
                      A[cbind(e, e)] <- 1; A }
walk <- function(Am, ts, fn) {
  P <- Matrix::Diagonal(nrow(Am), 1 - opt$alpha) + opt$alpha * Matrix::Matrix(Am, sparse = TRUE); Pt <- NULL
  for (t in seq_len(max(ts))) { Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; fn(S, t) } }
}
heldout_t <- function(M, kind) {
  h <- sample(rep(1:2, length.out = ncol(M))); ut <- which(upper.tri(matrix(0, G, G))); best <- integer(0)
  for (d in 1:2) {
    sv <- base_of(M[, h != d], kind)[ut]; target <- as.integer(sv >= stats::quantile(sv, 0.99)); rm(sv)
    crit <- numeric(opt$tmax)
    walk(soft(base_of(M[, h == d], kind)), seq_len(opt$tmax), function(S, t) crit[t] <<- ap_fast(S[ut], target))
    best <- c(best, which.max(crit))
  }
  as.integer(round(exp(mean(log(best)))))
}
evalS <- function(S) lapply(stats::setNames(names(bench), names(bench)), function(nm) {
  b <- bench[[nm]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]
  list(aupr = ap_fast(sc, b$pairs$label),
       jack = vapply(halves[[nm]], function(ix) ap_fast(sc[ix], b$pairs$label[ix]), 0),
       coherence = stats::median(regulon_coherence(S, b$reg, b$co, no_operon, genes), na.rm = TRUE)) })
kinds <- c("absr", "mi", "ns_absr", "clr_absr", "clr_mi", "clr_ns")
source("analysis/parallel.R")
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
ev <- foreach(kind = kinds, j = seq_along(kinds), .errorhandling = "pass") %dopar% {
  set.seed(opt$seed + j)
  S <- base_of(X, kind); out <- list(); out[[kind]] <- evalS(S); tstar <- NA_integer_
  if (startsWith(kind, "clr")) {
    tstar <- heldout_t(X, kind)
    walk(soft(S), tstar, function(St, t) out[[paste0("att_", kind)]] <<- evalS(St))
  }
  par_release(); list(out = out, tstar = tstar, kind = kind)
}
par_stop(pc)
bad <- vapply(ev, function(x) inherits(x, "error"), NA)
if (any(bad)) say("failed: ", paste(unique(vapply(ev[bad], conditionMessage, "")), collapse = " | "))
ev <- ev[!bad]
tst <- stats::setNames(vapply(ev, `[[`, 0L, "tstar"), vapply(ev, `[[`, "", "kind"))
E <- do.call(c, lapply(ev, `[[`, "out"))
rows <- list()
for (m in names(E)) for (nm in names(bench)) {
  ref <- if (startsWith(m, "att_")) sub("^att_", "", m) else if (startsWith(m, "clr_")) {
    c(clr_absr = "absr", clr_mi = "mi", clr_ns = "ns_absr")[[m]] } else NA
  d <- if (!is.na(ref) && !is.null(E[[ref]])) E[[m]][[nm]]$jack - E[[ref]][[nm]]$jack else NA
  d0 <- if (!is.na(ref) && !is.null(E[[ref]])) E[[m]][[nm]]$aupr - E[[ref]][[nm]]$aupr else NA
  rows[[length(rows) + 1L]] <- data.frame(method = m, truth = nm, aupr = E[[m]][[nm]]$aupr,
    base_rate = mean(bench[[nm]]$pairs$label), coherence = E[[m]][[nm]]$coherence,
    vs = ref, d_aupr = d0, d_lo = d0 - 1.96 * stats::sd(d), d_hi = d0 + 1.96 * stats::sd(d),
    t = if (startsWith(m, "att_")) tst[[sub("^att_", "", m)]] else NA)
}
res <- do.call(rbind, rows); res$aupr_ratio <- res$aupr / res$base_rate
utils::write.csv(res, file.path(opt$out, "metrics.csv"), row.names = FALSE)
for (nm in names(bench)) {
  x <- res[res$truth == nm, ]; x <- x[order(-x$aupr), ]
  say(sprintf("== %s (base rate %.4f)", nm, x$base_rate[1]))
  for (k in seq_len(nrow(x))) say(sprintf("   %-14s %s AUPR %.4f (x%.2f) coherence %.3f%s", x$method[k],
    if (is.na(x$t[k])) "    " else sprintf("t=%2d", x$t[k]), x$aupr[k], x$aupr_ratio[k], x$coherence[k],
    if (is.na(x$vs[k])) "" else sprintf("  | vs %-8s %+.4f [%+.4f, %+.4f]", x$vs[k], x$d_aupr[k], x$d_lo[k], x$d_hi[k])))
}
say("done")
