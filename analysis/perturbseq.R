#!/usr/bin/env Rscript
# Human test of completed CLR attention: K562 essential-gene Perturb-seq
# (Replogle et al. 2022). Networks are built de novo from the NON-TARGETING
# CONTROL cells only and scored against human truths (analysis/human_truths.R):
# DoRothEA A-B and A-C regulons, CollecTRI regulons, CORUM complexes, and
# knockdown-response modules built only from perturbed cells.
#
# Zero inflation and batch (10x lanes = gem groups):
#   pr      analytic Pearson residuals of a negative binomial with the cell's
#           UMI total as offset (Lause et al. 2021, theta = 100), computed WITHIN
#           each lane (depth and lane effects removed together, zeros in deep
#           cells become informative), clipped at +/- sqrt(n_cells)
#   ns      normal scores (Blom; average ranks for ties, so all zeros share one
#           score) computed WITHIN each lane
#   lognorm log1p(counts / total * median total), centred within lane (MI baseline)
# Methods (per representation): |r|, CLR on |r|; on lognorm also B-spline MI and
# CLR on MI; soft10 attention on every CLR variant at its own held-out depth.
# Held-out depth: lanes split at random into two halves (lanes play the role of
# laboratories), the target is the top 1% of the same CLR matrix on the other
# half; truths are never used.
#
# Usage (repo root):  Rscript analysis/perturbseq.R [--min_detect 0.01] [--out results/perturbseq]
#                       [--workers 6] [--tmax 32]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(h5 = "data/perturbseq/K562_essential_raw_singlecell_01.h5ad", out = "results/perturbseq",
            min_detect = 0.01, workers = 6L, mem_gb = 10, tmax = 32L, theta = 100, seed = 20261015L,
            alpha = 0.5, threads = NULL, jack = 100L, boot = 1000L, max_cells = 0L)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% c("workers", "tmax", "seed", "threads", "jack", "boot", "max_cells")) as.integer(v) else
    if (key %in% c("min_detect", "mem_gb", "theta", "alpha")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "perturbseq.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
set.seed(opt$seed)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

## ---- control cells ----------------------------------------------------------------
rd <- function(p) rhdf5::h5read(opt$h5, p)
cat_of <- function(grp, field) {
  x <- rd(paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L]
}
pert <- cat_of("/obs", "gene")
lane <- cat_of("/obs", "gem_group")
genes <- cat_of("/var", "gene_name")
ctrl <- which(pert == "non-targeting")
if (opt$max_cells > 0L && length(ctrl) > opt$max_cells) ctrl <- sort(sample(ctrl, opt$max_cells))   # quick tests only
say(sprintf("%d cells, %d genes; %d non-targeting control cells in %d lanes", length(pert), length(genes),
            length(ctrl), length(unique(lane[ctrl]))))
C <- rhdf5::h5read(opt$h5, "X", index = list(NULL, ctrl))     # genes x control cells (raw UMI)
storage.mode(C) <- "double"; rownames(C) <- genes
lane <- lane[ctrl]
det <- rowMeans(C > 0)
C <- C[det >= opt$min_detect, ]; genes <- rownames(C); G <- length(genes)
tot <- colSums(C)
say(sprintf("%d genes detected in >= %.0f%% of control cells; zero fraction %.1f%%; median UMI/cell %d",
            G, 100 * opt$min_detect, 100 * mean(C == 0), as.integer(stats::median(tot))))

## ---- representations -------------------------------------------------------------------
by_lane <- function(f) { R <- matrix(0, G, ncol(C), dimnames = dimnames(C))
  for (l in unique(lane)) { k <- lane == l; R[, k] <- f(C[, k, drop = FALSE]) }; R }
pearson_res <- function(Cl) {            # analytic NB Pearson residuals within a lane
  n <- colSums(Cl); mu <- outer(rowSums(Cl), n) / sum(Cl)
  z <- (Cl - mu) / sqrt(mu + mu^2 / opt$theta); z[!is.finite(z)] <- 0
  lim <- sqrt(ncol(C)); pmin(pmax(z, -lim), lim)
}
blom_lane <- function(Cl) { n <- ncol(Cl)
  t(apply(Cl, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4)))) }
lognorm_lane <- function(Cl) { L <- log1p(t(t(Cl) / colSums(Cl)) * stats::median(tot)); L - rowMeans(L) }
reps <- list(pr = by_lane(pearson_res), ns = by_lane(blom_lane), lognorm = by_lane(lognorm_lane))
rm(C); invisible(gc())

## ---- truths --------------------------------------------------------------------------
truths <- list(DoRothEA_AB = dorothea_sets(genes), DoRothEA_ABC = dorothea_sets(genes, c("A", "B", "C")),
               CollecTRI = collectri_sets(genes), CORUM = corum_sets(genes), KD_response = response_sets(genes))
no_operon <- stats::setNames(character(0), character(0))
bench <- lapply(truths, function(r) list(reg = r, pairs = comembership_pairs(r, no_operon, genes),
                                         co = comember_matrix(r, G)))
for (nm in names(bench)) say(sprintf("truth %-12s %4d sets; %d scored pairs, base rate %.4f", nm,
  length(bench[[nm]]$reg$members), length(bench[[nm]]$pairs$label), mean(bench[[nm]]$pairs$label)))

## ---- methods ----------------------------------------------------------------------------
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
abs_r <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; diag(r) <- 0; r }
calib <- function(M) { diag(M) <- 0; clr_calibrate(M, method = "normal", combine = "stouffer") }
mi_of <- function(X) { m <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = nthr())$mi
                       diag(m) <- 0; m }
base_of <- function(X, kind) switch(kind, absr = abs_r(X), clr_absr = calib(abs_r(X)), mi = mi_of(X),
                                    clr_mi = calib(mi_of(X)))
soft <- function(S) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L); e <- which(rowSums(A) == 0)
                      A[cbind(e, e)] <- 1; A }
walk <- function(A, ts, fn) {
  P <- Matrix::Diagonal(nrow(A), 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE); Pt <- NULL
  for (t in seq_len(max(ts))) { Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; fn(S, t) } }
}
# held-out depth on random lane halves, criterion = top 1% of the same CLR matrix on the other half
heldout_t <- function(X, kind) {
  ul <- sample(unique(lane)); h <- ifelse(lane %in% ul[seq_len(floor(length(ul) / 2))], 1L, 2L)
  ut <- which(upper.tri(matrix(0, G, G))); best <- integer(0)
  for (d in 1:2) {
    sv <- base_of(X[, h != d], kind)[ut]; target <- as.integer(sv >= stats::quantile(sv, 0.99)); rm(sv)
    crit <- numeric(opt$tmax)
    walk(soft(base_of(X[, h == d], kind)), seq_len(opt$tmax), function(S, t) crit[t] <<- ap_fast(S[ut], target))
    best <- c(best, which.max(crit))
  }
  as.integer(round(exp(mean(log(best)))))
}
tasks <- list(list(rep = "pr", kind = "absr"), list(rep = "pr", kind = "clr_absr"),
              list(rep = "ns", kind = "absr"), list(rep = "ns", kind = "clr_absr"),
              list(rep = "lognorm", kind = "absr"), list(rep = "lognorm", kind = "clr_absr"),
              list(rep = "lognorm", kind = "mi"), list(rep = "lognorm", kind = "clr_mi"))
score <- function(S, label) lapply(names(bench), function(nm) { b <- bench[[nm]]; p <- b$pairs
  sc <- S[cbind(p$i, p$j)]
  data.frame(method = label, truth = nm, aupr = ap_fast(sc, p$label), base_rate = mean(p$label),
             coherence = stats::median(regulon_coherence(S, b$reg, b$co, no_operon, genes), na.rm = TRUE)) })
source("analysis/parallel.R")
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
res <- foreach(tk = tasks, j = seq_along(tasks), .errorhandling = "pass") %dopar% {
  set.seed(opt$seed + j)
  X <- reps[[tk$rep]]; lab <- paste(tk$rep, tk$kind, sep = ":")
  S <- base_of(X, tk$kind); out <- score(S, lab)
  if (startsWith(tk$kind, "clr")) {
    tstar <- heldout_t(X, tk$kind)
    walk(soft(S), tstar, function(St, t) out <<- c(out, score(St, sprintf("%s:attention_t%d", lab, t))))
  }
  par_release(); do.call(rbind, out)
}
par_stop(pc)
bad <- vapply(res, function(x) inherits(x, "error"), NA)
if (any(bad)) say("failed: ", paste(unique(vapply(res[bad], conditionMessage, "")), collapse = " | "))
res <- do.call(rbind, res[!bad])
res$aupr_ratio <- res$aupr / res$base_rate
utils::write.csv(res, file.path(opt$out, "metrics.csv"), row.names = FALSE)
for (nm in names(bench)) {
  x <- res[res$truth == nm, ]; x <- x[order(-x$aupr), ]
  say(sprintf("== %s (base rate %.4f)", nm, x$base_rate[1]))
  for (k in seq_len(nrow(x))) say(sprintf("   %-34s AUPR %.4f (x%.2f)  coherence %.3f", x$method[k], x$aupr[k],
                                          x$aupr_ratio[k], x$coherence[k]))
}
say("done")
