#!/usr/bin/env Rscript
# Is attention depth a resolution parameter? For truths of increasing
# coarseness, score CLR (t = 0) and soft10 attention at every depth t = 1..tmax
# and report the depth at which each truth is recovered best.
#
#   E. coli (M3D, 466 experiments; RegulonDB, all evidence):
#     operons (finest; same-operon pairs are the positives here)
#     TF regulons (single-TF; operon pairs excluded, as in the main benchmark)
#     sigmulons (coarsest; operon pairs excluded; without sigma70 and sigmaS)
#   K562 (Perturb-seq essential compendium, network half A as in
#   perturbseq_bulk.R, same seed):
#     CORUM complexes < knockdown-response modules (half B) < DoRothEA A-B
#     < DoRothEA A-C < CollecTRI
# Prediction if depth sets resolution: argmax t increases from fine to coarse.
#
# Usage (repo root):  Rscript analysis/granularity.R [--tmax 32] [--out results/granularity]
#                       [--reuse_k562 results/granularity_with_sigS]   (E. coli only; K562 rows from that run)

args <- commandArgs(trailingOnly = TRUE)
opt <- list(tmax = 32L, out = "results/granularity", reuse_k562 = "", alpha = 0.5, seed = 20261015L, p_max = 1e-3,
            h5 = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad")
i <- 1L
while (i <= length(args)) { k <- sub("^--", "", args[i]); opt[[k]] <- if (k == "tmax") as.integer(args[i + 1L]) else args[i + 1L]; i <- i + 2L }
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "granularity.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")

soft <- function(S) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L); e <- which(rowSums(A) == 0)
                      A[cbind(e, e)] <- 1; A }
sweep <- function(C, bench, genes, operon_of, label) {
  G <- nrow(C); rows <- list()
  scoreS <- function(S, t) for (nm in names(bench)) { b <- bench[[nm]]; p <- b$pairs
    rows[[length(rows) + 1L]] <<- data.frame(data = label, truth = nm, t = t,
      aupr = ap_fast(S[cbind(p$i, p$j)], p$label), base_rate = mean(p$label),
      coherence = stats::median(regulon_coherence(S, b$reg, b$co, b$operon, genes), na.rm = TRUE)) }
  scoreS(C, 0L)
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(soft(C), sparse = TRUE); Pt <- NULL
  for (t in seq_len(opt$tmax)) { Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    S <- (Pt + t(Pt)) / 2; diag(S) <- 0; scoreS(S, t); say(sprintf("%s t=%d", label, t)) }
  do.call(rbind, rows)
}

## ---- E. coli ------------------------------------------------------------------------
m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; genes <- rownames(X); G <- length(genes)
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
reg_all <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)
pick <- function(keep) list(members = reg_all$members[keep], class = reg_all$class[keep])
tf_reg <- pick(reg_all$class == "TF" & lengths(reg_all$members) <= 500 & lengths(reg_all$members) >= 5)
# Sigmulons without the two global ones: sigma70 (housekeeping, 1,596 genes, 37%) and
# sigmaS/sigma38 (general stress, 331 genes): sigma32, sigma54, sigma24, sigma28, sigma19.
sig_reg <- pick(reg_all$class == "sigma" & lengths(reg_all$members) >= 5 &
                  !names(reg_all$members) %in% c("sigma:sigma70", "sigma:sigma38"))
op <- operon_of[names(operon_of) %in% genes]
op_sets <- lapply(split(match(names(op), genes), op), unique); op_sets <- op_sets[lengths(op_sets) >= 2]
op_reg <- list(members = op_sets, class = stats::setNames(rep("operon", length(op_sets)), names(op_sets)))
none <- stats::setNames(character(0), character(0))
mk <- function(r, o) list(reg = r, pairs = comembership_pairs(r, o, genes), co = comember_matrix(r, G), operon = o)
bench_e <- list(`1 operons` = mk(op_reg, none), `2 TF regulons` = mk(tf_reg, operon_of),
                `3 sigmulons` = mk(sig_reg, operon_of))
for (nm in names(bench_e)) say(sprintf("E. coli %s: %d sets, %d pairs, base rate %.4f", nm,
  length(bench_e[[nm]]$reg$members), length(bench_e[[nm]]$pairs$label), mean(bench_e[[nm]]$pairs$label)))
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")
res_e <- sweep(f$clr_scores, bench_e, genes, operon_of, "E. coli M3D"); f$release(); rm(f, X); invisible(gc())

## ---- K562 compendium -------------------------------------------------------------------
if (nzchar(opt$reuse_k562)) {
  res_h <- utils::read.csv(file.path(opt$reuse_k562, "granularity_curves.csv"), check.names = FALSE)
  res_h <- res_h[startsWith(res_h$data, "K562"), setdiff(names(res_h), "aupr_ratio")]
  say("K562 curves reused from ", opt$reuse_k562, " (", nrow(res_h), " rows)")
} else {
set.seed(opt$seed)                       # same split as perturbseq_bulk.R
rd <- function(p) rhdf5::h5read(opt$h5, p)
cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
Z <- rd("X"); hg <- cat_of("/var", "gene_name"); rownames(Z) <- hg
pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
pval <- as.numeric(rd("/obs/energy_test_p_value"))
K <- pert[pval < opt$p_max & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
Xh <- Z[, A]; rm(Z); Xh[!is.finite(Xh)] <- 0
Xh <- Xh[apply(Xh, 1, function(x) diff(range(x))) > 0, ]; hgenes <- rownames(Xh); H <- length(hgenes)
mkh <- function(r) list(reg = r, pairs = comembership_pairs(r, none, hgenes), co = comember_matrix(r, H), operon = none)
bench_h <- list(`1 CORUM complexes` = mkh(corum_sets(hgenes)),
                `2 KD-response modules` = mkh(response_sets(hgenes, path = opt$h5, use = B)),
                `3 DoRothEA A-B` = mkh(dorothea_sets(hgenes)),
                `4 DoRothEA A-C` = mkh(dorothea_sets(hgenes, c("A", "B", "C"))),
                `5 CollecTRI` = mkh(collectri_sets(hgenes)))
for (nm in names(bench_h)) say(sprintf("K562 %s: %d sets, base rate %.4f", nm,
  length(bench_h[[nm]]$reg$members), mean(bench_h[[nm]]$pairs$label)))
fh <- ClrAttention$new(Xh)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")
res_h <- sweep(fh$clr_scores, bench_h, hgenes, none, "K562 essential compendium")
}

res <- rbind(res_e, res_h); res$aupr_ratio <- res$aupr / res$base_rate
utils::write.csv(res, file.path(opt$out, "granularity_curves.csv"), row.names = FALSE)
best <- do.call(rbind, lapply(split(res, list(res$data, res$truth), drop = TRUE), function(x) data.frame(
  data = x$data[1], truth = x$truth[1], t_best_aupr = x$t[which.max(x$aupr)],
  gain_aupr = max(x$aupr) - x$aupr[x$t == 0], t_best_coherence = x$t[which.max(x$coherence)],
  gain_coherence = max(x$coherence) - x$coherence[x$t == 0])))
best <- best[order(best$data, best$truth), ]
utils::write.csv(best, file.path(opt$out, "granularity_best_t.csv"), row.names = FALSE)
say("best depth per truth (t = 0 is CLR without attention):")
for (k in seq_len(nrow(best))) say(sprintf("  %-26s %-22s AUPR best t = %2d (gain %+.4f) | coherence best t = %2d (gain %+.3f)",
  best$data[k], best$truth[k], best$t_best_aupr[k], best$gain_aupr[k], best$t_best_coherence[k], best$gain_coherence[k]))
say("done")
