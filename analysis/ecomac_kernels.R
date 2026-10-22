#!/usr/bin/env Rscript
# Dependence measure under batch effects: on EcoMAC (2,198 arrays, two
# Affymetrix platforms plus arrays without a GEO record), does CLR on
# B-spline MI beat CLR on correlation? Same benchmark as the M3D study
# (RegulonDB co-membership AUPR and median regulon coherence, strong/confirmed
# and all evidence), on all arrays and within each platform.
#
# Methods: |Pearson r|, raw MI, CLR on MI (HG bins, Stouffer), CLR on |r|,
# CLR on |r| of normal scores (Blom, average ties), and soft10 attention at
# t = 8 on each CLR variant (the M3D held-out depth; not re-selected here).
#
# Subset "all_centered" is the batch-corrected competitor: all arrays, with
# each gene centred and scaled within each platform group before any
# dependence is computed (the simplest location/scale batch correction).
#
# Usage (repo root):  Rscript analysis/ecomac_kernels.R [--subsets all,GPL199,GPL3154] [--threads N]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(ecomac = "data/ecomac", rdb = "data/RegulonDBExtract", out = "results/ecomac_kernels",
            subsets = "all,GPL199,GPL3154,noGEO", threads = NULL, t = 8L, alpha = 0.5)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key %in% c("threads", "t")) as.integer(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "ecomac_kernels.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R")

E <- as.matrix(utils::read.delim(file.path(opt$ecomac, "EcoMAC.tsv"), header = FALSE))
gn <- utils::read.delim(file.path(opt$ecomac, "gene_names.tsv"), header = FALSE, stringsAsFactors = FALSE)
rownames(E) <- gn[[2]]                                   # b-numbers
plat <- utils::read.csv(file.path(opt$ecomac, "array_platform.csv"), stringsAsFactors = FALSE)$platform
stopifnot(length(plat) == ncol(E))
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)

blom_avg <- function(X) { n <- ncol(X)
  t(apply(X, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4)))) }
abs_r <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; diag(r) <- 0; r }
calib <- function(M) { diag(M) <- 0; clr_calibrate(M, method = "normal", combine = "stouffer") }
att <- function(C) {
  A <- clr:::.softmax_rows(C, keff = 10, cap = 50L); G <- nrow(A)
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- as.matrix(P); for (s in seq_len(opt$t - 1L)) Pt <- as.matrix(Pt %*% P)
  S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S
}
sel <- list(all = rep(TRUE, ncol(E)), GPL199 = plat == "Antisense (GPL199)",
            GPL3154 = plat == "Genome 2.0 (GPL3154)", noGEO = plat == "no GEO record",
            all_centered = rep(TRUE, ncol(E)))
center_by <- function(X, grp) {
  for (g in unique(grp)) { k <- grp == g
    mu <- rowMeans(X[, k]); s <- apply(X[, k], 1, stats::sd); s[s <= 0] <- 1
    X[, k] <- (X[, k] - mu) / s }
  X
}
rows <- list()
for (sb in strsplit(opt$subsets, ",")[[1]]) {
  X <- E[, sel[[sb]]]
  if (sb == "all_centered") X <- center_by(X, plat)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  genes <- rownames(X); G <- length(genes)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck)
    list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes), co = comember_matrix(reg, G))
  })
  say(sprintf("%s: %d genes x %d arrays", sb, G, ncol(X)))
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
  Mi <- f$mi; diag(Mi) <- 0; f$release(); rm(f)
  R <- abs_r(X); Rn <- abs_r(blom_avg(X))
  meth <- list(abs_r = R, mi = Mi, clr_mi = calib(Mi), clr_absr = calib(R), clr_ns = calib(Rn))
  rm(Mi, R, Rn); invisible(gc())
  for (m in c("clr_mi", "clr_absr", "clr_ns")) meth[[paste0("att", opt$t, "_", m)]] <- att(meth[[m]])
  for (ev in names(bench)) {
    b <- bench[[ev]]; p <- b$pairs
    for (m in names(meth)) {
      S <- meth[[m]]
      rows[[length(rows) + 1L]] <- data.frame(subset = sb, n_arrays = ncol(X), evidence = ev, method = m,
        aupr = ap_fast(S[cbind(p$i, p$j)], p$label), base_rate = mean(p$label),
        coherence = stats::median(regulon_coherence(S, b$reg, b$co, operon_of, genes), na.rm = TRUE))
      say(sprintf("  [%s %s] %-14s AUPR %.4f  coherence %.3f", sb, ev, m,
                  rows[[length(rows)]]$aupr, rows[[length(rows)]]$coherence))
    }
  }
  rm(meth); invisible(gc())
  utils::write.csv(do.call(rbind, rows), file.path(opt$out, "kernels.csv"), row.names = FALSE)
}
say("done")
