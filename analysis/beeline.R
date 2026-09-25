#!/usr/bin/env Rscript
# BEELINE case study (Pratapa et al., Nature Methods 2020): CLR-as-attention on
# the seven experimental single-cell RNA-seq datasets, with BEELINE's gene
# selection and ground-truth networks, scored with BEELINE's metrics.
#
# Data: Zenodo record 3701939 (BEELINE-data.zip -> inputs/scRNA-Seq/<dataset>/
# ExpressionData.csv (genes x cells, log-normalized) and GeneOrdering.csv
# (VGAM p-value, variance); BEELINE-Networks.zip -> Networks/{human,mouse}/*.csv
# and {human,mouse}-tfs.csv).
#
# Gene selection (BEELINE "TFs + N"): genes with Bonferroni-corrected VGAM
# p < 0.01; all such genes that are TFs, plus the N most variable such genes.
#
# Ground truths per dataset (as in BEELINE): cell-type-specific ChIP-seq,
# non-specific ChIP-seq, STRING; plus the LOF/GOF network for mESC. Gene names
# are matched case-insensitively (mouse networks are upper case).
#
# Evaluation (TF-node, directed): candidate edges = (TF, gene) for every
# regulator of the ground truth present in the data and every other selected
# gene; label = edge in the ground truth. Symmetric scores are used for both
# directions. Metrics: AUPRC ratio (AUPRC / edge density = improvement over a
# random predictor) and early precision ratio (EPR: precision among the top-k
# edges, k = number of true edges, divided by density) -- BEELINE's two
# headline measures. Secondary: regulon co-membership (targets of a TF in the
# cell-type ChIP network as regulons, size 5..500) AUPR ratio.
#
# Methods (fixed in advance, same as the E. coli study): |Pearson|, |Spearman|,
# regularized partial correlation (linear-model baseline, ridge-shrunk
# precision matrix), raw B-spline MI, CLR (HG bins, Stouffer), CLR-attention
# soft10 and |Pearson|-attention pearson_top10 at depths chosen on held-out
# cells (random halves; criterion = predicting the other half's top-1% CLR
# pairs; RegulonDB-style gold standards are never used for the choice).
#
# Usage (repo root):
#   Rscript analysis/beeline.R [--beeline data/beeline] [--N 500,1000]
#       [--datasets hESC,hHep,mDC,mESC,mHSC-E,mHSC-GM,mHSC-L]
#       [--out results/beeline] [--threads N] [--workers W] [--mem_gb 1.5]
# The dataset x gene-set cases run in parallel (foreach over local cores
# minus 2, capped by RAM; analysis/parallel.R), each under its own seed.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(beeline = "data/beeline", N = "500,1000",
            datasets = "hESC,hHep,mDC,mESC,mHSC-E,mHSC-GM,mHSC-L",
            out = "results/beeline", threads = NULL, seed = 20260925L,
            alpha = 0.5, tmax = 64L, workers = NULL, mem_gb = 1.5)
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key %in% c("threads", "seed", "tmax", "workers")) as.integer(args[i + 1L]) else
    if (key == "mem_gb") as.numeric(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "beeline.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))       # pr_summary
source(file.path(script_dir, "regulons.R"))        # regulon_coherence etc.
set.seed(opt$seed)

root <- opt$beeline
expr_dir <- file.path(root, "BEELINE-data", "inputs", "scRNA-Seq")
net_dir <- file.path(root, "Networks")
tf_file <- function(sp) {
  cand <- c(file.path(root, paste0(sp, "-tfs.csv")),
            file.path(net_dir, paste0(sp, "-tfs.csv")))
  cand[file.exists(cand)][1]
}
species <- function(ds) if (ds %in% c("hESC", "hHep")) "human" else "mouse"
truths <- function(ds) {
  sp <- species(ds); d <- file.path(net_dir, sp)
  ct <- switch(ds, hESC = "hESC-ChIP-seq-network.csv",
               hHep = "HepG2-ChIP-seq-network.csv",
               mDC = "mDC-ChIP-seq-network.csv",
               mESC = "mESC-ChIP-seq-network.csv",
               "mHSC-ChIP-seq-network.csv")
  ns <- if (sp == "human") "Non-specific-ChIP-seq-network.csv" else
    "Non-Specific-ChIP-seq-network.csv"
  out <- list(celltype_ChIP = file.path(d, ct), nonspecific_ChIP = file.path(d, ns),
              STRING = file.path(d, "STRING-network.csv"))
  if (ds == "mESC") out$LOFGOF <- file.path(d, "mESC-lofgof-network.csv")
  # Evidence-filtered STRING (analysis/beeline_truths.R), when built.
  dd <- file.path(root, "Networks_derived", sp)
  der <- c(STRING_regulatory = "STRING-regulatory.csv", STRING_curated = "STRING-curated.csv")
  for (k in names(der))
    if (file.exists(file.path(dd, der[[k]]))) out[[k]] <- file.path(dd, der[[k]])
  out
}

# Graded ChIP-seq (tools/get_chipatlas.R): for each TF with a ChIP-Atlas table,
# the Spearman correlation between the TF's row of a score matrix and its
# binding scores (MACS2, +/- 1 kb of the TSS; 0 where no peak), over all other
# selected genes. Two binding sources: experiments in matching cell types
# ("celltype") and the average over all experiments ("all"). Computed within
# each TF, because TFs differ in antibody quality and peak-score scale.
chip_dir <- file.path("data", "chipatlas")
graded_chip <- function(S, genes, ds) {
  g <- if (species(ds) == "human") "hg38" else "mm10"
  key <- if (grepl("^mHSC", ds)) "mHSC" else ds
  up <- toupper(genes); out <- list()
  for (ti in which(file.exists(file.path(chip_dir, g, paste0(up, ".tsv"))))) {
    tb <- utils::read.delim(file.path(chip_dir, g, paste0(up[ti], ".tsv")), stringsAsFactors = FALSE)
    for (src in c("celltype", "all")) {
      col <- if (src == "all") "all" else key
      if (!col %in% names(tb) || all(is.na(tb[[col]]))) next
      b <- numeric(length(genes)); hit <- match(up, tb$gene)
      b[!is.na(hit)] <- tb[[col]][hit[!is.na(hit)]]; b[is.na(b)] <- 0
      keep <- seq_along(genes) != ti
      if (sum(b[keep] > 0) < 20) next
      for (m in names(S))
        out[[length(out) + 1L]] <- data.frame(tf = up[ti], source = src,
          method = sub("_t[0-9]+$", "", m), n_bound = sum(b[keep] > 0),
          rho = suppressWarnings(stats::cor(S[[m]][ti, keep], b[keep], method = "spearman")))
    }
  }
  do.call(rbind, out)
}

select_genes <- function(ds, N) {
  go <- utils::read.csv(file.path(expr_dir, ds, "GeneOrdering.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  names(go)[1] <- "gene"
  sig <- go[go$VGAMpValue * nrow(go) < 0.01, ]
  sig <- sig[order(-sig$Variance), ]
  tfs <- toupper(utils::read.csv(tf_file(species(ds)), stringsAsFactors = FALSE)$TF)
  sel <- union(sig$gene[toupper(sig$gene) %in% tfs], utils::head(sig$gene, N))
  list(genes = sel, tfs = tfs)
}

## ---- helpers (same constructions as the E. coli study) -------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
clr_fit <- function(X) {
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none",
                                       threads = nthr())
  f$calibrate(method = "normal", combine = "stouffer")
  f
}
att_mass <- function(A, ts) {
  G <- nrow(A)
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  out <- list(); Pt <- NULL
  for (t in seq_len(max(ts))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; out[[as.character(t)]] <- S }
  }
  out
}
soft_op <- function(f) { f$build_operator(alpha = opt$alpha, softmax_keff = 10,
                                          softmax_cap = 50L); f$operator }
pear_op <- function(R) { Ar <- abs(R); diag(Ar) <- 0; clr:::.topk_rows(Ar, 10L) }
partial_cor <- function(Z) {
  C <- tcrossprod(Z) / (ncol(Z) - 1)
  lam <- 0.1 * mean(diag(C))                  # fixed ridge (not tuned on truth)
  Om <- solve(C + lam * diag(nrow(C)))
  d <- sqrt(diag(Om))
  P <- -Om / outer(d, d); diag(P) <- 0
  abs(P)
}

# Held-out depth: split cells into random halves; operator on one half,
# criterion = AUPR for the other half's top-1% CLR pairs.
heldout_depth <- function(X, make_A, label) {
  grid <- 2L^(0:floor(log2(opt$tmax)))
  cells <- sample(ncol(X)); h <- split(cells, rep(1:2, length.out = length(cells)))
  best <- integer(0)
  for (d in 1:2) {
    tr <- X[, h[[d]]]; te <- X[, h[[3 - d]]]
    keep <- apply(tr, 1, function(x) diff(range(x))) > 0 &
      apply(te, 1, function(x) diff(range(x))) > 0
    fte <- clr_fit(te[keep, ]); ut <- upper.tri(fte$clr_scores)
    sv <- fte$clr_scores[ut]; target <- as.integer(sv >= stats::quantile(sv, 0.99))
    A <- make_A(tr[keep, ])
    all_t <- sort(unique(c(grid, 3L, 5L, 6L, 7L, 10L, 12L)))
    Sm <- att_mass(A, all_t)
    crit <- vapply(Sm, function(S) pr_summary(S[ut], target)[["aupr"]], numeric(1))
    best <- c(best, as.integer(names(crit)[which.max(crit)]))
  }
  t <- as.integer(round(exp(mean(log(best)))))
  say(sprintf("    held-out depth %-14s directions %s -> t* = %d", label,
              paste(best, collapse = ","), t))
  t
}

score_truth <- function(S, genes, tfs_in, gt) {
  gi <- match(toupper(gt$Gene1), toupper(genes))
  gj <- match(toupper(gt$Gene2), toupper(genes))
  ok <- !is.na(gi) & !is.na(gj) & gi != gj
  regs <- sort(unique(gi[ok]))
  if (length(regs) < 1 || sum(ok) < 10) return(NULL)
  G <- length(genes)
  cand <- cbind(rep(regs, each = G), rep(seq_len(G), times = length(regs)))
  cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
  key <- cand[, 1] * (G + 1) + cand[, 2]
  pos <- unique(gi[ok] * (G + 1) + gj[ok])
  lab <- as.integer(key %in% pos)
  s <- S[cand]
  ps <- pr_summary(s, lab)
  dens <- mean(lab); k <- sum(lab)
  o <- order(-s)
  epr <- mean(lab[o[seq_len(k)]]) / dens
  c(n_regulators = length(regs), n_edges = k, density = dens,
    auprc_ratio = ps[["aupr"]] / dens, epr = epr, auroc = ps[["auroc"]])
}

comembership_truth <- function(S, genes, gt) {
  gi <- match(toupper(gt$Gene1), toupper(genes))
  gj <- match(toupper(gt$Gene2), toupper(genes))
  ok <- !is.na(gi) & !is.na(gj) & gi != gj
  mem <- lapply(split(gj[ok], gi[ok]), unique)
  mem <- mem[lengths(mem) >= 5 & lengths(mem) <= 500]
  if (length(mem) < 3) return(NULL)
  reg <- list(members = mem, class = stats::setNames(rep("TF", length(mem)), names(mem)))
  G <- length(genes)
  pairs <- comembership_pairs(reg, stats::setNames(character(0), character(0)),
                              stats::setNames(genes, genes))
  s <- S[cbind(pairs$i, pairs$j)]
  ps <- pr_summary(s, pairs$label)
  c(n_regulons = length(mem), comember_aupr_ratio = ps[["aupr"]] / mean(pairs$label))
}

## ---- main loop: one task per dataset x gene set, in parallel ----------------
run_case <- function(ds, N) {
  rows <- list(); crow <- list()
  ex <- utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"),
                        row.names = 1, check.names = FALSE)
  ex <- as.matrix(ex)
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  sel <- select_genes(ds, N)
  X <- ex[intersect(sel$genes, rownames(ex)), , drop = FALSE]; rm(ex)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]
  genes <- rownames(X)
  say(sprintf("%s TFs+%d: %d genes (%d TFs) x %d cells", ds, N, nrow(X),
              sum(toupper(genes) %in% sel$tfs), ncol(X)))
  Z <- Zs(X); R <- tcrossprod(Z) / (ncol(Z) - 1)
  Rs <- stats::cor(t(X), method = "spearman")
  f <- clr_fit(X); M <- f$mi; diag(M) <- 0
  t_soft <- heldout_depth(X, function(Y) { g <- clr_fit(Y); A <- soft_op(g); g$release(FALSE); A },
                          paste(ds, N, "soft10"))
  t_pear <- heldout_depth(X, function(Y) pear_op(stats::cor(t(Y))), paste(ds, N, "pearson_top10"))
  depth <- data.frame(dataset = ds, N = N, t_soft10 = t_soft, t_pearson = t_pear)
  S <- list(abs_pearson = abs(R), abs_spearman = abs(Rs),
            partial_cor = partial_cor(Z), mi = M, clr = f$clr_scores)
  S[[sprintf("clr_attention_soft10_t%d", t_soft)]] <- att_mass(soft_op(f), t_soft)[[1]]
  S[[sprintf("pearson_attention_top10_t%d", t_pear)]] <- att_mass(pear_op(R), t_pear)[[1]]
  f$release(FALSE); rm(f, Z, R, Rs, M)
  S <- lapply(S, function(s) { s[is.na(s)] <- 0; diag(s) <- 0; s })
  for (tn in names(gts)) for (m in names(S)) {
    r <- score_truth(S[[m]], genes, NULL, gts[[tn]])
    if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(dataset = ds, N = N, truth = tn,
                                            method = sub("_t[0-9]+$", "", m),
                                            depth = suppressWarnings(as.integer(sub("^.*_t", "", m))),
                                            t(r))
  }
  for (m in names(S)) {
    r <- comembership_truth(S[[m]], genes, gts$celltype_ChIP)
    if (is.null(r)) next
    crow[[length(crow) + 1L]] <- data.frame(dataset = ds, N = N,
                                            method = sub("_t[0-9]+$", "", m), t(r))
  }
  gr <- graded_chip(S, genes, ds)
  if (!is.null(gr)) gr <- cbind(dataset = ds, N = N, gr)
  # Depth diagnostic (never used to choose t): CLR attention at fixed depths.
  drow <- list()
  Sd <- att_mass(soft_op(clr_fit(X)), c(1, 2, 4, 8, 16))
  for (tt in names(Sd)) {
    Sx <- Sd[[tt]]; Sx[is.na(Sx)] <- 0; diag(Sx) <- 0
    for (tn in names(gts)) {
      r <- score_truth(Sx, genes, NULL, gts[[tn]])
      if (!is.null(r)) drow[[length(drow) + 1L]] <- data.frame(dataset = ds, N = N, truth = tn,
                                                               depth = as.integer(tt), t(r))
    }
  }
  rm(Sd)
  res <- do.call(rbind, rows)
  cur <- res[res$truth == "celltype_ChIP", ]
  for (k in seq_len(nrow(cur)))
    say(sprintf("    %s TFs+%d %-28s AUPRC ratio %.3f  EPR %.3f  (celltype ChIP)",
                ds, N, cur$method[k], cur$auprc_ratio[k], cur$epr[k]))
  rm(S); invisible(gc())
  list(rows = res, crow = do.call(rbind, crow), depth = depth, graded = gr,
       diag = do.call(rbind, drow))
}

cases <- expand.grid(N = as.integer(strsplit(opt$N, ",")[[1]]),
                     ds = strsplit(opt$datasets, ",")[[1]], stringsAsFactors = FALSE)
source(file.path(script_dir, "parallel.R"))
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
out <- foreach(ci = seq_len(nrow(cases)), .inorder = TRUE, .errorhandling = "pass") %dopar% {
  if (par_should_stop()) return(NULL)
  set.seed(opt$seed + ci)
  run_case(cases$ds[ci], cases$N[ci])
}
par_stop(pc)
bad <- vapply(out, function(x) inherits(x, "error"), NA)
if (any(bad)) say("cases failed: ", paste(sprintf("%s/%d: %s", cases$ds[bad], cases$N[bad],
                                                  vapply(out[bad], conditionMessage, "")),
                                          collapse = " | "))
out <- out[!bad & !vapply(out, is.null, NA)]
rows <- lapply(out, `[[`, "rows"); crow <- lapply(out, `[[`, "crow")
depths <- lapply(out, `[[`, "depth")
graded <- do.call(rbind, lapply(out, `[[`, "graded"))
diagd <- do.call(rbind, lapply(out, `[[`, "diag"))
if (!is.null(diagd)) utils::write.csv(diagd, file.path(opt$out, "diag_depths.csv"), row.names = FALSE)
if (!is.null(graded)) {
  utils::write.csv(graded, file.path(opt$out, "graded_chip.csv"), row.names = FALSE)
  ref <- graded[graded$method == "clr", c("dataset", "N", "tf", "source", "rho")]
  names(ref)[5] <- "rho_clr"
  gm <- merge(graded, ref)
  for (src in c("celltype", "all")) {
    x0 <- gm[gm$source == src, ]
    if (!nrow(x0)) next
    say(sprintf("graded ChIP-seq (%s experiments): per-TF Spearman rho, median over %d TF x case",
                src, nrow(x0[x0$method == "clr", ])))
    for (m in unique(x0$method)) {
      x <- x0[x0$method == m & !is.na(x0$rho), ]
      say(sprintf("   %-26s median rho %+.3f  rho > 0 in %.0f%%  beats CLR in %.0f%%  (Wilcoxon vs CLR p = %.2g)",
                  m, stats::median(x$rho), 100 * mean(x$rho > 0), 100 * mean(x$rho > x$rho_clr),
                  if (m == "clr") NA_real_ else suppressWarnings(stats::wilcox.test(x$rho, x$rho_clr, paired = TRUE)$p.value)))
    }
  }
}
res <- do.call(rbind, rows); utils::write.csv(res, file.path(opt$out, "tfnode_metrics.csv"), row.names = FALSE)
cm <- do.call(rbind, crow); utils::write.csv(cm, file.path(opt$out, "comembership_metrics.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, depths), file.path(opt$out, "heldout_depths.csv"), row.names = FALSE)

# Summary: per method, median over dataset x N x truth of the ratio to |Pearson|,
# and win counts vs CLR.
ref <- res[res$method == "abs_pearson", c("dataset", "N", "truth", "auprc_ratio", "epr")]
names(ref)[4:5] <- c("ap_ref", "epr_ref")
m <- merge(res, ref)
summ <- do.call(rbind, lapply(split(m, m$method), function(x) data.frame(
  method = x$method[1], n = nrow(x),
  median_auprc_ratio = stats::median(x$auprc_ratio),
  median_epr = stats::median(x$epr),
  median_auprc_vs_pearson = stats::median(x$auprc_ratio / x$ap_ref),
  median_epr_vs_pearson = stats::median(x$epr / x$epr_ref))))
clr_ref <- res[res$method == "clr", c("dataset", "N", "truth", "auprc_ratio")]
names(clr_ref)[4] <- "ap_clr"
w <- merge(res, clr_ref)
summ$frac_beats_clr_auprc <- vapply(summ$method, function(mm) {
  x <- w[w$method == mm, ]; mean(x$auprc_ratio > x$ap_clr) }, numeric(1))
summ <- summ[order(-summ$median_auprc_ratio), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
say("summary over datasets x gene sets x ground truths:")
for (k in seq_len(nrow(summ)))
  say(sprintf("  %-28s median AUPRC ratio %.3f  EPR %.3f  | vs |Pearson| x%.3f / x%.3f | beats CLR in %.0f%%",
              summ$method[k], summ$median_auprc_ratio[k], summ$median_epr[k],
              summ$median_auprc_vs_pearson[k], summ$median_epr_vs_pearson[k],
              100 * summ$frac_beats_clr_auprc[k]))
say("done")
