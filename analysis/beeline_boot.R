#!/usr/bin/env Rscript
# BEELINE uncertainty (docs/COMPUTE_CONSTRAINTS.md items 3 and 8):
#   A. depth stability: nsplit random halvings of the cells per case
#      (dataset x gene set), each used in both directions, with the held-out
#      criterion of beeline.R; reports the distribution of t*;
#   B. nested cell bootstrap: B draws of cells with replacement per case; in
#      every draw the depth is re-selected on one random halving of the
#      draw's cells, then |Pearson|, raw MI, CLR, CLR attention (soft10),
#      |Pearson| attention (pearson_top10), CLR-NS (CLR on |r| of normal
#      scores, see beeline.R) and CLR-NS attention, each attention at the
#      draw's own t*, are scored against every ground truth (AUPRC ratio, EPR).
# Differences against CLR are paired within draw. Per case: percentile
# intervals. Across cases: in every draw index b, the median over cases of
# the difference (each case resampled independently), giving an interval for
# the typical difference per truth type.
#
# Gene selection, truths, constructions and metrics are those of beeline.R.
# --context as in beeline.R: "none" (BEELINE's selected genes only), K (the K
# most variable expressed genes plus the selected ones) or "all"; CLR and the
# operators are built on the context, scoring and the held-out criterion use
# the selected genes only.
#
# Usage (repo root):
#   Rscript analysis/beeline_boot.R [--B 200] [--nsplit 20] [--context none|K|all]
#       [--out results/beeline_boot] [--workers W] [--mem_gb 1.5]
# Each finished task is written to tasks/<case>_<kind><b>.csv; restarts skip them.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(beeline = "data/beeline", N = "500,1000",
            datasets = "hESC,hHep,mDC,mESC,mHSC-E,mHSC-GM,mHSC-L",
            out = "results/beeline_boot", B = 200L, nsplit = 20L, threads = NULL,
            seed = 20261015L, alpha = 0.5, tmax = 64L, workers = NULL, mem_gb = 1.5,
            context = "none")
ints <- c("B", "nsplit", "threads", "seed", "tmax", "workers")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else if (key %in% c("mem_gb", "alpha")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(file.path(opt$out, "tasks"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "beeline_boot.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))       # pr_summary, ap_fast
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

## ---- data, truths and gene selection (as beeline.R) ---------------------------
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
# Candidate (TF, gene) edges and labels per truth, precomputed once per case.
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
context_genes <- function(ex, ds, score_genes) {      # as beeline.R
  expressed <- rownames(ex)[rowMeans(ex > 0) >= 0.01]
  if (opt$context == "all") return(union(score_genes, expressed))
  go <- utils::read.csv(file.path(expr_dir, ds, "GeneOrdering.csv"), stringsAsFactors = FALSE,
                        check.names = FALSE)
  names(go)[1] <- "gene"
  go <- go[go$gene %in% expressed, ]
  union(score_genes, utils::head(go$gene[order(-go$Variance)], as.integer(opt$context)))
}
Ns <- as.integer(strsplit(opt$N, ",")[[1]])
cases <- list()
for (ds in strsplit(opt$datasets, ",")[[1]]) {
  ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1,
                                  check.names = FALSE))
  gts <- lapply(truth_files(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  for (N in Ns) {
    sel <- select_genes(ds, N)
    fg <- if (opt$context == "none") sel else context_genes(ex, ds, sel)
    X <- ex[intersect(fg, rownames(ex)), , drop = FALSE]
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE]
    si <- which(rownames(X) %in% sel)
    es <- Filter(Negate(is.null), lapply(gts, function(gt) edge_sets(rownames(X)[si], gt)))
    cases[[sprintf("%s_%d", ds, N)]] <- list(ds = ds, N = N, X = X, si = si, es = es)
  }
  rm(ex, gts)
}
say(sprintf("%d cases (context %s): %s", length(cases), opt$context,
  paste(sprintf("%s (%d scored / %d fitted genes x %d cells)", names(cases),
  vapply(cases, function(k) length(k$si), 0L), vapply(cases, function(k) nrow(k$X), 0L),
  vapply(cases, function(k) ncol(k$X), 0L)), collapse = ", ")))

## ---- constructions (as beeline.R) -----------------------------------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
jitter_const <- function(X) {
  rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0))
    X[rng <= 0, ] <- X[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9)
  X
}
clr_fit <- function(X) {
  f <- ClrAttention$new(jitter_const(X))$estimate_mi(bins = "hg", transform = "none", threads = nthr())
  f$calibrate(method = "normal", combine = "stouffer")
  f
}
att_walk <- function(A, ts, fn) {
  P <- Matrix::Diagonal(nrow(A), 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- NULL
  for (t in seq_len(max(ts))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; fn(S, t) }
  }
}
soft_op <- function(f) { f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
                         f$operator }
# CLR-NS: CLR on |r| of normal scores (Blom, average ranks for ties); see beeline.R
blom_avg <- function(X) {
  n <- ncol(X)
  t(apply(X, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4))))
}
ns_clr <- function(X) {
  r <- abs(stats::cor(t(blom_avg(X)))); r[is.na(r)] <- 0; diag(r) <- 0
  clr_calibrate(r, method = "normal", combine = "stouffer")
}
pear_op <- function(Y) { Ar <- abs(stats::cor(t(Y))); Ar[is.na(Ar)] <- 0; diag(Ar) <- 0
                         clr:::.topk_rows(Ar, 10L) }
grid_t <- sort(unique(c(2L^(0:floor(log2(opt$tmax))), 3L, 5L, 6L, 7L, 10L, 12L)))
# One random halving of the cells of X, both directions: argmax depth per
# direction for soft10 and pearson_top10 (criterion as in beeline.R).
halving_depths <- function(X, si) {
  cells <- sample(ncol(X)); h <- split(cells, rep(1:2, length.out = length(cells)))
  out <- list()
  for (d in 1:2) {
    tr <- X[, h[[d]]]; te <- X[, h[[3 - d]]]
    keep <- apply(tr, 1, function(x) diff(range(x))) > 0 & apply(te, 1, function(x) diff(range(x))) > 0
    fte <- clr_fit(te[keep, ])
    sk <- which(which(keep) %in% si)                # scored genes among kept
    ut <- upper.tri(matrix(0, length(sk), length(sk)))
    sv <- fte$clr_scores[sk, sk][ut]; target_mi <- as.integer(sv >= stats::quantile(sv, 0.99))
    fte$release(FALSE)
    sv <- ns_clr(te[keep, ])[sk, sk][ut]; target_ns <- as.integer(sv >= stats::quantile(sv, 0.99))
    for (op in c("soft10", "pearson_top10", "ns_soft10")) {
      # each CLR operator is judged against its own CLR in the other half
      target <- if (op == "ns_soft10") target_ns else target_mi
      A <- if (op == "soft10") { g <- clr_fit(tr[keep, ]); a <- soft_op(g); g$release(FALSE); a } else
        if (op == "ns_soft10") clr:::.softmax_rows(ns_clr(tr[keep, ]), keff = 10, cap = 50L) else
        pear_op(tr[keep, ])
      crit <- stats::setNames(numeric(length(grid_t)), grid_t)
      att_walk(A, grid_t, function(S, t) crit[[as.character(t)]] <<- ap_fast(S[sk, sk][ut], target))
      out[[length(out) + 1L]] <- data.frame(direction = d, operator = op,
                                            t_best = as.integer(names(crit)[which.max(crit)]))
    }
  }
  r <- do.call(rbind, out)
  r$t_sel <- ave(r$t_best, r$operator, FUN = function(v) as.integer(round(exp(mean(log(v))))))
  r
}
score_all <- function(S, es) do.call(rbind, lapply(names(es), function(tn) {
  e <- es[[tn]]; s <- S[e$cand]; s[is.na(s)] <- 0
  dens <- mean(e$lab); k <- sum(e$lab)
  o <- order(-s)
  data.frame(truth = tn, auprc_ratio = ap_fast(s, e$lab) / dens,
             epr = mean(e$lab[o[seq_len(k)]]) / dens)
}))

## ---- tasks -------------------------------------------------------------------------
run_task <- function(tk) {
  cs <- cases[[tk$case]]
  if (tk$kind == "split") {
    r <- halving_depths(cs$X, cs$si)
    return(cbind(case = tk$case, kind = "split", b = tk$b, r))
  }
  X <- cs$X[, sample(ncol(cs$X), replace = TRUE)]
  dep <- halving_depths(X, cs$si)
  tsel <- stats::setNames(dep$t_sel[!duplicated(dep$operator)], dep$operator[!duplicated(dep$operator)])
  rows <- list()
  add <- function(S, m, t) { S <- S[cs$si, cs$si]; S[is.na(S)] <- 0; diag(S) <- 0
    rows[[length(rows) + 1L]] <<- cbind(method = m, t = t, score_all(S, cs$es)) }
  Z <- Zs(X)
  add(abs(tcrossprod(Z) / (ncol(Z) - 1)), "abs_pearson", NA); rm(Z)
  f <- clr_fit(X)
  add(f$mi, "mi", NA); add(f$clr_scores, "clr", NA)
  att_walk(soft_op(f), tsel[["soft10"]], function(S, t) add(S, "clr_attention_soft10", t))
  f$release(FALSE)
  att_walk(pear_op(X), tsel[["pearson_top10"]], function(S, t) add(S, "pearson_attention_top10", t))
  Cns <- ns_clr(X); add(Cns, "clr_ns", NA)
  att_walk(clr:::.softmax_rows(Cns, keff = 10, cap = 50L), tsel[["ns_soft10"]],
           function(S, t) add(S, "clr_ns_attention_soft10", t))
  cbind(case = tk$case, kind = "boot", b = tk$b, do.call(rbind, rows))
}
tasks <- list()
for (cn in names(cases)) {
  for (b in seq_len(opt$nsplit)) tasks[[length(tasks) + 1L]] <- list(case = cn, kind = "split", b = b)
  for (b in seq_len(opt$B)) tasks[[length(tasks) + 1L]] <- list(case = cn, kind = "boot", b = b)
}
tname <- function(tk) sprintf("%s_%s%03d", tk$case, tk$kind, tk$b)
done <- sub("\\.csv$", "", list.files(file.path(opt$out, "tasks"), "\\.csv$"))
todo <- Filter(function(tk) !tname(tk) %in% done, tasks)
say(sprintf("%d tasks (%d cases x (%d splits + %d draws)); %d done, %d to run", length(tasks),
            length(cases), opt$nsplit, opt$B, length(done), length(todo)))
source(file.path(script_dir, "parallel.R"))
if (length(todo)) {
  pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
  st <- foreach(tk = todo, .inorder = FALSE, .errorhandling = "pass") %dopar% {
    if (par_should_stop()) return("stopped")
    set.seed(opt$seed + match(tk$case, names(cases)) * 100000L + (tk$kind == "boot") * 10000L + tk$b)
    d <- run_task(tk)
    nm <- tname(tk)
    tmp <- file.path(opt$out, "tasks", paste0(nm, ".csv.part"))
    utils::write.csv(d, tmp, row.names = FALSE)
    file.rename(tmp, file.path(opt$out, "tasks", paste0(nm, ".csv")))
    par_release(); "ok"
  }
  par_stop(pc)
  bad <- vapply(st, function(x) inherits(x, "error"), NA)
  if (any(bad)) say("tasks failed (", sum(bad), "): ",
                    paste(unique(vapply(st[bad], conditionMessage, "")), collapse = " | "))
}

## ---- summaries ---------------------------------------------------------------------
fl <- list.files(file.path(opt$out, "tasks"), "\\.csv$", full.names = TRUE)
sp <- do.call(rbind, lapply(fl[grepl("_split[0-9]+\\.csv$", fl)], utils::read.csv))
bt <- do.call(rbind, lapply(fl[grepl("_boot[0-9]+\\.csv$", fl)], utils::read.csv))
utils::write.csv(sp, file.path(opt$out, "cell_splits.csv"), row.names = FALSE)
utils::write.csv(bt, file.path(opt$out, "boot_draws.csv"), row.names = FALSE)
q <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)
s1 <- unique(sp[, c("case", "b", "operator", "t_sel")])
ts <- do.call(rbind, lapply(split(s1, list(s1$case, s1$operator), drop = TRUE), function(x)
  data.frame(case = x$case[1], operator = x$operator[1], n = nrow(x), t_median = stats::median(x$t_sel),
             t_q25 = q(x$t_sel, .25), t_q75 = q(x$t_sel, .75), t_min = min(x$t_sel), t_max = max(x$t_sel))))
utils::write.csv(ts, file.path(opt$out, "cell_split_tstar.csv"), row.names = FALSE)
say("held-out t* over cell halvings (per case):")
for (k in seq_len(nrow(ts)))
  say(sprintf("  %-14s %-14s median %g  IQR [%g, %g]  range [%d, %d]", ts$case[k], ts$operator[k],
              ts$t_median[k], ts$t_q25[k], ts$t_q75[k], ts$t_min[k], ts$t_max[k]))
ref <- bt[bt$method == "clr", c("case", "b", "truth", "auprc_ratio", "epr")]
names(ref)[4:5] <- c("ap_clr", "epr_clr")
m <- merge(bt, ref)
m$d_ap <- m$auprc_ratio - m$ap_clr; m$d_epr <- m$epr - m$epr_clr
pc <- do.call(rbind, lapply(split(m, list(m$case, m$truth, m$method), drop = TRUE), function(x)
  data.frame(case = x$case[1], truth = x$truth[1], method = x$method[1], n = nrow(x),
             auprc_ratio = mean(x$auprc_ratio), d_ap = mean(x$d_ap), d_ap_lo = q(x$d_ap, .025),
             d_ap_hi = q(x$d_ap, .975), epr = mean(x$epr), d_epr = mean(x$d_epr),
             d_epr_lo = q(x$d_epr, .025), d_epr_hi = q(x$d_epr, .975))))
utils::write.csv(pc, file.path(opt$out, "per_case_intervals.csv"), row.names = FALSE)
m$truth_type <- sub("_(regulatory|curated)$", "_filtered", m$truth)
agg <- stats::aggregate(cbind(d_ap, d_epr) ~ b + truth_type + method, data = m, FUN = stats::median)
ov <- do.call(rbind, lapply(split(agg, list(agg$truth_type, agg$method), drop = TRUE), function(x)
  data.frame(truth_type = x$truth_type[1], method = x$method[1], n_draws = nrow(x),
             median_d_ap = stats::median(x$d_ap), d_ap_lo = q(x$d_ap, .025), d_ap_hi = q(x$d_ap, .975),
             median_d_epr = stats::median(x$d_epr), d_epr_lo = q(x$d_epr, .025), d_epr_hi = q(x$d_epr, .975))))
utils::write.csv(ov, file.path(opt$out, "overall_intervals.csv"), row.names = FALSE)
say("across cases: median over cases of (method - CLR), 2.5-97.5% over draws:")
for (k in seq_len(nrow(ov))) if (ov$method[k] != "clr")
  say(sprintf("  %-18s %-24s dAUPRC ratio %+.3f [%+.3f, %+.3f]  dEPR %+.3f [%+.3f, %+.3f]",
              ov$truth_type[k], ov$method[k], ov$median_d_ap[k], ov$d_ap_lo[k], ov$d_ap_hi[k],
              ov$median_d_epr[k], ov$d_epr_lo[k], ov$d_epr_hi[k]))
pcs <- pc[pc$method == "clr_attention_soft10", ]
say(sprintf("CLR attention vs CLR, per case x truth: AUPRC interval above 0 in %d, below 0 in %d, spans 0 in %d of %d",
            sum(pcs$d_ap_lo > 0), sum(pcs$d_ap_hi < 0), sum(pcs$d_ap_lo <= 0 & pcs$d_ap_hi >= 0), nrow(pcs)))
say("done")
