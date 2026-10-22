#!/usr/bin/env Rscript
# Does attention still add to CLR once global expression programs are removed?
#
# Removing the top k principal components before co-expression network
# inference is a standard, label-free way to remove confounders and global
# programs (Parsana et al. 2019, Genome Biology, choosing k with the
# permutation test of sva::num.sv). On EcoMAC it raised CLR's regulon recovery
# sharply, so the question for the paper is whether the attention gain
# survives on PC-corrected data or was partly doing the same job.
#
# Data sets: M3D 466 replicate-averaged experiments; EcoMAC after ComBat by
# platform. For k in --ks and for k_be = sva::num.sv(method = "be") (the
# label-free choice), score CLR(MI), CLR(|r|) and soft10 attention on each
# (t = --t, the held-out depth) against RegulonDB: co-membership AUPR and
# median coherence, strong/confirmed and all evidence. The attention-minus-CLR
# difference gets a paired delete-half gene jackknife interval.
#
# Usage (repo root):  Rscript analysis/pc_removal.R [--ks 0,1,2,3,5,8,12,20] [--t 9]
#                       [--sets m3d,ecomac] [--workers 4] [--out results/pc_removal]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(ks = "0,1,2,3,5,8,12,20", t = 9L, sets = "m3d,ecomac", workers = 4L, mem_gb = 6,
            out = "results/pc_removal", jack = 100L, seed = 20261014L, alpha = 0.5, threads = NULL)
ints <- c("t", "workers", "jack", "seed", "threads")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else if (key %in% c("mem_gb", "alpha")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(file.path(opt$out, "tasks"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "pc_removal.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R")
set.seed(opt$seed)

## ---- data sets -----------------------------------------------------------------------
data_sets <- list()
sets <- strsplit(opt$sets, ",")[[1]]
if ("m3d" %in% sets) {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                           check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
  data_sets$m3d <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
}
if ("ecomac" %in% sets) {
  E <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
  rownames(E) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
  plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
  E <- E[apply(E, 1, function(x) diff(range(x))) > 0, ]
  data_sets$ecomac <- suppressMessages(sva::ComBat(E, batch = factor(plat)))
}
ks <- as.integer(strsplit(opt$ks, ",")[[1]])
kbe <- list()
for (d in names(data_sets)) {
  X <- data_sets[[d]]
  kbe[[d]] <- sva::num.sv(X, mod = matrix(1, ncol(X), 1), method = "be", B = 20, seed = opt$seed)
  say(sprintf("%s: %d genes x %d samples; num.sv (be) k = %d", d, nrow(X), ncol(X), kbe[[d]]))
}
remove_pcs <- function(Y, k) {
  if (k == 0L) return(Y)
  Yc <- Y - rowMeans(Y); V <- svd(Yc, nu = 0, nv = k)$v
  Y - Yc %*% V %*% t(V)
}
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)

## ---- one (data set, k) ------------------------------------------------------------------
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
calib <- function(M) { diag(M) <- 0; clr_calibrate(M, method = "normal", combine = "stouffer") }
att <- function(C) {
  A <- clr:::.softmax_rows(C, keff = 10, cap = 50L); G <- nrow(A)
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- as.matrix(P); for (s in seq_len(opt$t - 1L)) Pt <- as.matrix(Pt %*% P)
  S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S
}
run_one <- function(d, k) {
  X <- remove_pcs(data_sets[[d]], k); genes <- rownames(X); G <- length(genes)
  mi <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = nthr())$mi
  r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0
  S <- list(clr_mi = calib(mi), clr_r = calib(r)); rm(mi, r)
  S$att_clr_mi <- att(S$clr_mi); S$att_clr_r <- att(S$clr_r)
  out <- list()
  for (ev in c("SC", "all")) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck)
    p <- comembership_pairs(reg, operon_of, genes); co <- comember_matrix(reg, G)
    sc <- lapply(S, function(M) M[cbind(p$i, p$j)])
    J <- t(replicate(opt$jack, { keep <- sample(p$genes_in, floor(length(p$genes_in) / 2))
      inh <- logical(G); inh[keep] <- TRUE; sel <- inh[p$i] & inh[p$j]
      vapply(sc, function(s) ap_fast(s[sel], p$label[sel]), 0) }))
    for (m in names(S)) {
      base <- if (startsWith(m, "att_")) sub("^att_", "", m) else m
      dd <- J[, m] - J[, base]
      out[[length(out) + 1L]] <- data.frame(set = d, k = k, evidence = ev, method = m,
        aupr = ap_fast(sc[[m]], p$label),
        coherence = stats::median(regulon_coherence(S[[m]], reg, co, operon_of, genes), na.rm = TRUE),
        d_vs_clr = ap_fast(sc[[m]], p$label) - ap_fast(sc[[base]], p$label),
        d_se = stats::sd(dd))
    }
  }
  do.call(rbind, out)
}
grid <- do.call(rbind, lapply(names(data_sets), function(d)
  data.frame(d = d, k = sort(unique(c(ks, kbe[[d]]))), stringsAsFactors = FALSE)))
source("analysis/parallel.R")
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
res <- foreach(j = seq_len(nrow(grid)), .errorhandling = "pass") %dopar% {
  f <- file.path(opt$out, "tasks", sprintf("%s_k%02d.csv", grid$d[j], grid$k[j]))
  if (file.exists(f)) return(utils::read.csv(f))
  set.seed(opt$seed + j); r <- run_one(grid$d[j], grid$k[j])
  utils::write.csv(r, f, row.names = FALSE); par_release(); r
}
par_stop(pc)
bad <- vapply(res, function(x) inherits(x, "error"), NA)
if (any(bad)) say("failed: ", paste(unique(vapply(res[bad], conditionMessage, "")), collapse = " | "))
res <- do.call(rbind, res[!bad])
res$kbe <- vapply(res$set, function(d) kbe[[d]], 0)
utils::write.csv(res, file.path(opt$out, "pc_removal.csv"), row.names = FALSE)
for (d in names(data_sets)) for (ev in c("SC", "all")) {
  say(sprintf("== %s [%s]  (num.sv k = %d): AUPR / coherence; attention - CLR with 95%% jackknife", d, ev, kbe[[d]]))
  x <- res[res$set == d & res$evidence == ev, ]
  for (k in sort(unique(x$k))) {
    y <- x[x$k == k, ]; g <- function(m, f) y[[f]][y$method == m]
    say(sprintf("  k=%2d%s CLR(MI) %.4f/%.3f  att %.4f/%.3f  d %+.4f [%+.4f, %+.4f] | CLR(|r|) %.4f/%.3f  att %.4f/%.3f  d %+.4f [%+.4f, %+.4f]",
                k, if (k == kbe[[d]]) "*" else " ",
                g("clr_mi", "aupr"), g("clr_mi", "coherence"), g("att_clr_mi", "aupr"), g("att_clr_mi", "coherence"),
                g("att_clr_mi", "d_vs_clr"), g("att_clr_mi", "d_vs_clr") - 1.96 * g("att_clr_mi", "d_se"),
                g("att_clr_mi", "d_vs_clr") + 1.96 * g("att_clr_mi", "d_se"),
                g("clr_r", "aupr"), g("clr_r", "coherence"), g("att_clr_r", "aupr"), g("att_clr_r", "coherence"),
                g("att_clr_r", "d_vs_clr"), g("att_clr_r", "d_vs_clr") - 1.96 * g("att_clr_r", "d_se"),
                g("att_clr_r", "d_vs_clr") + 1.96 * g("att_clr_r", "d_se")))
  }
}
say("done")
