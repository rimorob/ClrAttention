#!/usr/bin/env Rscript
# Label-free choice of batch correction for EcoMAC (2,198 arrays; Antisense
# GPL199, Genome 2.0 GPL3154 and arrays without a GEO record).
#
# Variants: correction in {none, ComBat by platform, ComBat by platform x GEO
# series (series with one array merged within their platform)} crossed with
# removal of the top k principal components (k in --ks) of the corrected,
# gene-centred matrix (an unsupervised stand-in for surrogate variables:
# sva needs a model of interest, which unsupervised inference does not have).
#
# Selection criterion (RegulonDB never used): authors (labs) are split at
# random into two halves of similar array counts. The TARGET is built from
# the held-out labs WITHIN each platform group separately, on uncorrected
# data -- mean of the per-platform CLR(MI) matrices -- and its top 1% of pairs
# are the positives. It contains no between-platform offsets and is the same
# for every variant. Each variant's pooled CLR network on the training labs
# is scored by AUPR for that target (CLR on MI and CLR on |r|). Shared
# platform artefacts therefore cannot raise the criterion; biology that
# reproduces across labs can.
# For information only, every variant is also scored against RegulonDB on all
# arrays (co-membership AUPR and median coherence, strong/confirmed evidence).
#
# Usage (repo root):
#   Rscript analysis/ecomac_batch.R [--nsplit 6] [--ks 0,1,3,5,10,20]
#       [--out results/ecomac_batch] [--workers 5] [--mem_gb 6]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(nsplit = 6L, ks = "0,1,3,5,10,20", out = "results/ecomac_batch", workers = 5L,
            mem_gb = 6, seed = 20261014L, target_frac = 0.01, min_platform = 30L, threads = NULL)
ints <- c("nsplit", "workers", "seed", "min_platform", "threads")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else if (key %in% c("mem_gb", "target_frac")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(file.path(opt$out, "tasks"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "ecomac_batch.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr) })
source("analysis/regulondb.R"); source("analysis/regulons.R")
set.seed(opt$seed)
ks <- as.integer(strsplit(opt$ks, ",")[[1]])

## ---- data ----------------------------------------------------------------------------
E <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
gn <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)
rownames(E) <- gn[[2]]
smp <- utils::read.delim("data/ecomac/samples.tsv", stringsAsFactors = FALSE, check.names = FALSE)
plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
series <- smp[["GEO accession"]]; author <- smp[["Author"]]
stopifnot(length(plat) == ncol(E), length(author) == ncol(E))
E <- E[apply(E, 1, function(x) diff(range(x))) > 0, ]
genes <- rownames(E); G <- length(genes)

## ---- corrected variants ---------------------------------------------------------------
bs <- paste(plat, series, sep = ":")
small <- bs %in% names(which(table(bs) < 2))
bs[small] <- paste(plat[small], "small", sep = ":")
corr <- list(none = E)
say("ComBat by platform (", length(unique(plat)), " batches)")
corr$combat_platform <- suppressMessages(sva::ComBat(E, batch = factor(plat)))
say("ComBat by platform x series (", length(unique(bs)), " batches)")
corr$combat_platform_series <- suppressMessages(sva::ComBat(E, batch = factor(bs)))
remove_pcs <- function(Y, k) {
  if (k == 0L) return(Y)
  Yc <- Y - rowMeans(Y)
  V <- svd(Yc, nu = 0, nv = k)$v
  Y - Yc %*% V %*% t(V)
}
variants <- list()
for (cn in names(corr)) for (k in ks)
  variants[[sprintf("%s_pc%02d", cn, k)]] <- remove_pcs(corr[[cn]], k)
rm(corr); invisible(gc())
say(sprintf("%d variants: %s", length(variants), paste(names(variants), collapse = ", ")))

## ---- splits by author ----------------------------------------------------------------
author_split <- function() {
  sz <- table(author); sz <- sz[sample(length(sz))]   # some author fields are empty
  half <- integer(length(sz)); tot <- c(0, 0)
  for (j in seq_along(sz)) { h <- which.min(tot); half[j] <- h; tot[h] <- tot[h] + sz[[j]] }
  half[match(author, names(sz))]
}
splits <- lapply(seq_len(opt$nsplit), function(s) author_split())

## ---- constructions -------------------------------------------------------------------
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
jit <- function(X) { rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0)) X[rng <= 0, ] <- X[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9); X }
clr_mi <- function(X) ClrAttention$new(jit(X))$estimate_mi(bins = "hg", transform = "none", threads = nthr())$
  calibrate(method = "normal", combine = "stouffer")$clr_scores
clr_r <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; diag(r) <- 0
  clr_calibrate(r, method = "normal", combine = "stouffer") }
ut <- NULL

## ---- tasks: targets first, then (split, variant) -----------------------------------------
source("analysis/parallel.R")
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
tdir <- file.path(opt$out, "tasks")
targets <- foreach(s = seq_len(opt$nsplit), .errorhandling = "stop") %dopar% {
  f <- file.path(tdir, sprintf("target_%02d.rds", s))
  if (file.exists(f)) return(f)
  te <- splits[[s]] == 2L
  groups <- names(which(table(plat[te]) >= opt$min_platform))
  S <- Reduce(`+`, lapply(groups, function(p) clr_mi(E[, te & plat == p]))) / length(groups)
  u <- which(upper.tri(S)); sv <- S[u]
  saveRDS(list(pos = as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)), groups = groups,
               n = vapply(groups, function(p) sum(te & plat == p), 0L)), f)
  f
}
say("targets: ", paste(vapply(targets, function(f) { x <- readRDS(f)
  paste(sprintf("%s=%d", sub(" .*", "", x$groups), x$n), collapse = "/") }, ""), collapse = "; "))
grid <- expand.grid(s = seq_len(opt$nsplit), v = names(variants), stringsAsFactors = FALSE)
res <- foreach(j = seq_len(nrow(grid)), .errorhandling = "pass") %dopar% {
  s <- grid$s[j]; v <- grid$v[j]
  f <- file.path(tdir, sprintf("crit_%02d_%s.csv", s, v))
  if (file.exists(f)) return(utils::read.csv(f))
  tg <- readRDS(targets[[s]])
  X <- variants[[v]][, splits[[s]] == 1L]
  u <- which(upper.tri(matrix(0, G, G)))
  d <- data.frame(split = s, variant = v,
                  crit_clr_mi = ap_fast(clr_mi(X)[u], tg$pos),
                  crit_clr_r = ap_fast(clr_r(X)[u], tg$pos))
  utils::write.csv(d, f, row.names = FALSE); par_release(); d
}
bad <- vapply(res, function(x) inherits(x, "error"), NA)
if (any(bad)) say("tasks failed: ", paste(unique(vapply(res[bad], conditionMessage, "")), collapse = " | "))
crit <- do.call(rbind, res[!bad])
utils::write.csv(crit, file.path(opt$out, "criterion.csv"), row.names = FALSE)

## ---- information only: RegulonDB on all arrays ------------------------------------------
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
reg <- make_regulons(rtab, genes, conf_keep = c("C", "S"))
bench <- list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes), co = comember_matrix(reg, G))
truth <- foreach(v = names(variants), .errorhandling = "pass") %dopar% {
  f <- file.path(tdir, sprintf("truth_%s.csv", v))
  if (file.exists(f)) return(utils::read.csv(f))
  p <- bench$pairs
  d <- do.call(rbind, lapply(c("clr_mi", "clr_r"), function(m) {
    S <- if (m == "clr_mi") clr_mi(variants[[v]]) else clr_r(variants[[v]])
    data.frame(variant = v, method = m, aupr = ap_fast(S[cbind(p$i, p$j)], p$label),
               coherence = stats::median(regulon_coherence(S, bench$reg, bench$co, operon_of, genes), na.rm = TRUE))
  }))
  utils::write.csv(d, f, row.names = FALSE); par_release(); d
}
par_stop(pc)
truth <- do.call(rbind, truth[!vapply(truth, function(x) inherits(x, "error"), NA)])
utils::write.csv(truth, file.path(opt$out, "regulondb_info.csv"), row.names = FALSE)

## ---- summary ----------------------------------------------------------------------------
cs <- stats::aggregate(cbind(crit_clr_mi, crit_clr_r) ~ variant, data = crit, FUN = mean)
csd <- stats::aggregate(cbind(crit_clr_mi, crit_clr_r) ~ variant, data = crit, FUN = stats::sd)
names(csd)[2:3] <- c("sd_clr_mi", "sd_clr_r")
w <- reshape(truth, idvar = "variant", timevar = "method", direction = "wide")
summ <- merge(merge(cs, csd), w)
summ <- summ[order(-summ$crit_clr_mi), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
say(sprintf("label-free criterion (mean over %d author splits) and, for information, RegulonDB (all arrays, SC):", opt$nsplit))
for (k in seq_len(nrow(summ)))
  say(sprintf("  %-30s crit CLR(MI) %.4f (sd %.4f)  CLR(|r|) %.4f | RegulonDB CLR(MI) %.4f / %.3f  CLR(|r|) %.4f / %.3f",
              summ$variant[k], summ$crit_clr_mi[k], summ$sd_clr_mi[k], summ$crit_clr_r[k],
              summ$aupr.clr_mi[k], summ$coherence.clr_mi[k], summ$aupr.clr_r[k], summ$coherence.clr_r[k]))
say("selected by the label-free criterion (CLR(MI)): ", summ$variant[1])
say("done")
