#!/usr/bin/env Rscript
# Depth selection over many lab splits, the full depth curve, and the softmax
# candidate-cap sweep, on the 466 replicate-averaged M3D experiments
# (docs/COMPUTE_CONSTRAINTS.md items 3, 4 and 10).
#
# followup.R chose t* on ONE balanced split of the 39 labs, by a geometric scan
# plus golden-section search. Here:
#   - nsplit random balanced lab splits (labs shuffled, then each assigned to
#     the half with fewer experiments), each used in both directions;
#   - the held-out criterion (AUPR of attention mass on the train-half
#     operator for the test half's top-1% CLR pairs; RegulonDB never used) is
#     evaluated at EVERY integer depth t = 1..tmax, so the whole curve is kept;
#   - per split, t* = geometric mean of the two directions' argmax; per
#     operator, the reported depth is the median t* over splits;
#   - operators: soft10 at softmax caps 50 (the default), 100, 200 and no cap
#     (every positive CLR score), and the pearson_top10 control;
#   - kernel ablation (added 2026-10-14; select with --ops): the same softmax
#     (keff 10, cap 50) on other score matrices --
#       soft10_clr_absr    CLR (Stouffer) on |Pearson r|
#       soft10_clr_ns      CLR on |r| of normal scores (Blom, average ties)
#       soft10_absr_noclr  raw |r|, no CLR calibration
#       soft10_mi_noclr    raw B-spline MI, no CLR calibration
#     For these, the held-out target is the top 1% of the SAME score matrix
#     computed on the test half (for soft10 and its caps that is CLR on MI,
#     as before; pearson_top10 keeps the CLR-on-MI target of followup.R).
# Then each operator is scored against RegulonDB on the full data at its
# median t* and at the earlier t = 8, with the followup.R intervals: paired
# delete-half jackknife over genes (pooled co-membership AUPR) and paired
# bootstrap over regulons (median coherence), both relative to CLR. A depth
# curve of pooled co-membership AUPR on the full data (diagnostic only) is
# saved for every operator.
#
# Outputs (--out): split_assignments.csv, heldout_curves.csv (split, direction,
# operator, t, criterion), heldout_tstar.csv (per split and direction),
# tstar_summary.csv, ci_comembership.csv, ci_coherence.csv,
# fulldata_depth_curve.csv, lab_splits.log.
#
# Usage (repo root):
#   Rscript analysis/lab_splits.R [--nsplit 50] [--tmax 64]
#       [--ops soft10,soft10_cap100,soft10_cap200,soft10_capall,pearson_top10]
#       [--out results/lab_splits] [--jack 100] [--boot 2000] [--quick 0]
#       [--workers W] [--mem_gb 4]
# Splits run in parallel (analysis/parallel.R); each finished split is written
# to splits/<split>.csv at once, and a restart skips finished splits.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            nsplit = 50L, tmax = 64L,
            ops = "soft10,soft10_cap100,soft10_cap200,soft10_capall,pearson_top10",
            out = "results/lab_splits", jack = 100L, boot = 2000L, quick = 0L,
            seed = 20261014L, alpha = 0.5, target_frac = 0.01, threads = NULL,
            workers = NULL, mem_gb = 4)
ints <- c("nsplit", "tmax", "jack", "boot", "quick", "seed", "threads", "workers")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else
    if (key %in% c("alpha", "target_frac", "mem_gb")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(file.path(opt$out, "splits"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "lab_splits.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
wcsv <- function(x, n) utils::write.csv(x, file.path(opt$out, n), row.names = FALSE)
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))
source(file.path(script_dir, "regulons.R"))
set.seed(opt$seed)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
ops <- strsplit(opt$ops, ",")[[1]]

## ---- data ------------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
raw <- utils::read.delim(file.path(opt$m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                         check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (opt$quick > 0L) {
  v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]
}
genes <- rownames(X); G <- length(genes)
desc <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_descriptions")),
                          stringsAsFactors = FALSE)
exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description)
exper[!grepl("experimenter:", desc$description)] <- NA
lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])
lab <- lab[colnames(X)]
lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab)))
rm(desc, exper)
say(sprintf("compendium: %d genes x %d experiments from %d labs", G, ncol(X),
            length(unique(lab))))

# Random balanced lab splits: shuffle the labs, then give each to the half
# with fewer experiments so far.
lab_split <- function() {
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab])
}
splits <- lapply(seq_len(opt$nsplit), function(s) lab_split())
wcsv(data.frame(experiment = colnames(X), lab = unname(lab),
                do.call(cbind, stats::setNames(splits, sprintf("split_%02d", seq_along(splits))))),
     "split_assignments.csv")

## ---- constructions (as in followup.R / chips_bootstrap.R) --------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
build_fit <- function(Xm) {
  rng <- apply(Xm, 1, function(x) diff(range(x)))
  if (any(rng <= 0))
    Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9)
  f <- ClrAttention$new(Xm)$estimate_mi(bins = "hg", transform = "none", threads = nthr())
  f$calibrate(method = "normal", combine = "stouffer")
  f
}
cap_of <- function(op) switch(op, soft10 = 50L, soft10_cap100 = 100L,
                              soft10_cap200 = 200L, soft10_capall = .Machine$integer.max,
                              NA_integer_)
blom_avg <- function(Xm) {
  n <- ncol(Xm)
  t(apply(Xm, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4))))
}
abs_r <- function(Xm) { r <- abs(stats::cor(t(Xm))); r[is.na(r)] <- 0; diag(r) <- 0; r }
calib <- function(M) clr_calibrate(M, method = "normal", combine = "stouffer")
ablation_ops <- c("soft10_clr_absr", "soft10_clr_ns", "soft10_absr_noclr", "soft10_mi_noclr")
# Score matrix the softmax is applied to (and the held-out target is taken from).
base_scores <- function(f, Xm, op) switch(op,
  soft10_clr_absr = calib(abs_r(Xm)), soft10_clr_ns = calib(abs_r(blom_avg(Xm))),
  soft10_absr_noclr = abs_r(Xm), soft10_mi_noclr = { m <- f$mi; diag(m) <- 0; m },
  f$clr_scores)
# Attention matrix A (rows sum to 1) for an operator on a fit and its |Pearson|.
make_A <- function(f, Xm, op) {
  if (op %in% ablation_ops) {
    A <- clr:::.softmax_rows(base_scores(f, Xm, op), keff = 10, cap = 50L)
    empty <- rowSums(A) == 0; A[cbind(which(empty), which(empty))] <- 1   # self-loop
    return(A)
  }
  if (startsWith(op, "soft10")) {
    f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = cap_of(op))
    return(f$operator)
  }
  if (op == "pearson_top10") {
    Ar <- abs(tcrossprod(Zs(Xm)) / (ncol(Xm) - 1)); diag(Ar) <- 0
    f$set_operator(clr:::.topk_rows(Ar, 10L), alpha = opt$alpha, label = op)
    return(f$operator)
  }
  stop("unknown operator: ", op)
}
# Lazy operator P = (1-a)I + aA: sparse when A is sparse, dense (BLAS) when
# more than 10% of A is non-zero (the uncapped softmax).
lazy_P <- function(A) {
  if (mean(A != 0) > 0.10) { P <- opt$alpha * A; diag(P) <- diag(P) + 1 - opt$alpha; return(P) }
  Matrix::Diagonal(nrow(A), 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
}
# Call fn(S, t) with attention mass S = (P^t + P^t')/2 for t = 1..tmax.
depth_walk <- function(P, tmax, fn) {
  Pt <- NULL
  for (t in seq_len(tmax)) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    S <- Pt + t(Pt); S <- S / 2; diag(S) <- 0
    fn(S, t); rm(S)
    if (t %% 8L == 0L) invisible(gc(verbose = FALSE))
  }
  invisible(NULL)
}

## ---- held-out criterion curves, one task per split ---------------------------
done <- sub("\\.csv$", "", list.files(file.path(opt$out, "splits"), "\\.csv$"))
todo <- setdiff(sprintf("split_%02d", seq_along(splits)), done)
say(sprintf("%d splits already done; %d to run; operators: %s; t = 1..%d",
            length(done), length(todo), paste(ops, collapse = ", "), opt$tmax))
run_split <- function(si) {
  h <- splits[[si]]
  fits <- lapply(1:2, function(d) build_fit(X[, h == d]))
  ut <- which(upper.tri(matrix(0, G, G)))
  rows <- list()
  for (dir in 1:2) {
    tr <- dir; te <- 3L - dir
    sv <- fits[[te]]$clr_scores[ut]
    target_clr <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)); rm(sv)
    for (op in ops) {
      target <- target_clr
      if (op %in% ablation_ops) {
        sv <- base_scores(fits[[te]], X[, h == te], op)[ut]
        target <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)); rm(sv)
      }
      A <- make_A(fits[[tr]], X[, h == tr], op)
      P <- lazy_P(A); rm(A)
      crit <- numeric(opt$tmax)
      depth_walk(P, opt$tmax, function(S, t) crit[t] <<- ap_fast(S[ut], target))
      rm(P)
      rows[[length(rows) + 1L]] <- data.frame(split = si, direction = dir, operator = op,
                                              t = seq_len(opt$tmax), criterion = crit)
      invisible(gc(verbose = FALSE))
    }
  }
  for (f in fits) f$release()
  do.call(rbind, rows)
}
source(file.path(script_dir, "parallel.R"))
if (length(todo)) {
  pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
  st <- foreach(nm = todo, .inorder = FALSE, .errorhandling = "pass") %dopar% {
    if (par_should_stop()) return("stopped")
    si <- as.integer(sub("split_", "", nm))
    set.seed(opt$seed + si)
    t0 <- Sys.time()
    d <- run_split(si)
    tmp <- file.path(opt$out, "splits", paste0(nm, ".csv.part"))
    utils::write.csv(d, tmp, row.names = FALSE)
    file.rename(tmp, file.path(opt$out, "splits", paste0(nm, ".csv")))
    par_release()
    sprintf("%s %.0f s", nm, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  par_stop(pc)
  bad <- vapply(st, function(x) inherits(x, "error"), NA)
  if (any(bad)) say("splits failed: ", paste(vapply(st[bad], conditionMessage, ""), collapse = " | "))
  if (any(unlist(st[!bad]) == "stopped")) say("stop file seen: remaining splits skipped")
  say("split timings: ", paste(unlist(st[!bad]), collapse = "; "))
}
cur <- do.call(rbind, lapply(list.files(file.path(opt$out, "splits"), "\\.csv$", full.names = TRUE),
                             utils::read.csv))
wcsv(cur, "heldout_curves.csv")
ts <- do.call(rbind, lapply(split(cur, list(cur$split, cur$direction, cur$operator), drop = TRUE),
  function(x) { k <- which.max(x$criterion)
    # flatness: depths whose criterion is within 1% of the maximum
    near <- x$t[x$criterion >= 0.99 * x$criterion[k]]
    data.frame(split = x$split[1], direction = x$direction[1], operator = x$operator[1],
               t_best = x$t[k], crit_best = x$criterion[k],
               t_within1pct_lo = min(near), t_within1pct_hi = max(near)) }))
wcsv(ts, "heldout_tstar.csv")
tsplit <- do.call(rbind, lapply(split(ts, list(ts$split, ts$operator), drop = TRUE), function(x)
  data.frame(split = x$split[1], operator = x$operator[1],
             t_star = as.integer(round(exp(mean(log(x$t_best))))))))
tsum <- do.call(rbind, lapply(split(tsplit, tsplit$operator), function(x) {
  y <- ts[ts$operator == x$operator[1], ]
  data.frame(operator = x$operator[1], n_splits = nrow(x),
             t_median = as.integer(stats::median(x$t_star)),
             t_q25 = stats::quantile(x$t_star, 0.25, names = FALSE),
             t_q75 = stats::quantile(x$t_star, 0.75, names = FALSE),
             t_min = min(x$t_star), t_max = max(x$t_star),
             within1pct_lo_median = stats::median(y$t_within1pct_lo),
             within1pct_hi_median = stats::median(y$t_within1pct_hi))
}))
wcsv(tsum, "tstar_summary.csv")
say("held-out t* over splits (t* per split = geometric mean of the two directions):")
for (k in seq_len(nrow(tsum)))
  say(sprintf("  %-15s median t* %d  IQR [%g, %g]  range [%d, %d]  | criterion within 1%% of max: t in [%g, %g] (median over split-directions)",
              tsum$operator[k], tsum$t_median[k], tsum$t_q25[k], tsum$t_q75[k],
              tsum$t_min[k], tsum$t_max[k], tsum$within1pct_lo_median[k],
              tsum$within1pct_hi_median[k]))

## ---- full-data scoring against RegulonDB --------------------------------------
say("full-data scoring at median t* and t = 8")
to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)
bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
  ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
  reg <- make_regulons(rtab, genes, conf_keep = ck)
  list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes),
       co = comember_matrix(reg, G))
})
fit <- build_fit(X)
R0 <- abs(tcrossprod(Zs(X)) / (ncol(X) - 1))
methods <- list(abs_pearson = R0, clr_hg_stouffer = fit$clr_scores)
curve_rows <- list()
for (op in ops) {
  tsel <- tsum$t_median[tsum$operator == op]
  keep_t <- sort(unique(c(8L, tsel)))
  P <- lazy_P(make_A(fit, X, op))
  depth_walk(P, opt$tmax, function(S, t) {
    for (ev in names(bench)) {
      p <- bench[[ev]]$pairs
      curve_rows[[length(curve_rows) + 1L]] <<- data.frame(operator = op, evidence = ev, t = t,
        comember_aupr = ap_fast(S[cbind(p$i, p$j)], p$label))
    }
    if (t %in% keep_t)
      methods[[sprintf("%s:t%03d%s", op, t, if (t == tsel) "*" else "")]] <<- S
  })
  rm(P); invisible(gc(verbose = FALSE))
}
wcsv(do.call(rbind, curve_rows), "fulldata_depth_curve.csv")
fit$release(); rm(fit, R0)
methods <- lapply(methods, function(S) { diag(S) <- 0; S })

# Intervals exactly as in followup.R: delete-half jackknife over genes for
# pooled AUPR, paired bootstrap over regulons for coherence, both vs CLR.
pooled <- list(); coh <- list()
for (ev in names(bench)) {
  p <- bench[[ev]]$pairs
  sc <- lapply(methods, function(S) S[cbind(p$i, p$j)])
  full <- vapply(sc, function(s) pr_summary(s, p$label)[["aupr"]], numeric(1))
  gin <- p$genes_in
  J <- matrix(NA_real_, opt$jack, length(methods), dimnames = list(NULL, names(methods)))
  for (r in seq_len(opt$jack)) {
    keep <- sample(gin, floor(length(gin) / 2))
    inh <- logical(G); inh[keep] <- TRUE
    sel <- inh[p$i] & inh[p$j]
    J[r, ] <- vapply(sc, function(s) ap_fast(s[sel], p$label[sel]), numeric(1))
  }
  se <- apply(J, 2, stats::sd)
  se_d <- apply(J - J[, "clr_hg_stouffer"], 2, stats::sd)
  d0 <- full - full[["clr_hg_stouffer"]]
  pooled[[ev]] <- data.frame(evidence = ev, method = names(methods), aupr = full,
                             se = se, diff_vs_clr = d0, diff_se = se_d,
                             diff_lo = d0 - 1.96 * se_d, diff_hi = d0 + 1.96 * se_d,
                             base_rate = mean(p$label))
  b <- bench[[ev]]
  Cm <- vapply(methods, function(S) regulon_coherence(S, b$reg, b$co, operon_of, genes),
               numeric(length(b$reg$members)))
  ok <- stats::complete.cases(Cm); Cm <- Cm[ok, , drop = FALSE]
  dC <- Cm - Cm[, "clr_hg_stouffer"]
  bt <- replicate(opt$boot, { ix <- sample(nrow(dC), replace = TRUE)
                              apply(dC[ix, , drop = FALSE], 2, stats::median) })
  coh[[ev]] <- data.frame(evidence = ev, method = names(methods),
                          median_auroc = apply(Cm, 2, stats::median),
                          median_diff_vs_clr = apply(dC, 2, stats::median),
                          diff_lo = apply(bt, 1, stats::quantile, 0.025),
                          diff_hi = apply(bt, 1, stats::quantile, 0.975),
                          frac_regulons_improved = colMeans(dC > 0),
                          n_regulons = nrow(dC))
}
pooled <- do.call(rbind, pooled); coh <- do.call(rbind, coh)
wcsv(pooled, "ci_comembership.csv"); wcsv(coh, "ci_coherence.csv")
for (ev in names(bench)) {
  say(sprintf("[%s] co-membership AUPR (vs CLR, 95%% jackknife) | median coherence (vs CLR, 95%% regulon bootstrap):", ev))
  x <- pooled[pooled$evidence == ev, ]; y <- coh[coh$evidence == ev, ]
  for (k in seq_len(nrow(x)))
    say(sprintf("  %-24s AUPR %.4f  %+.4f [%+.4f, %+.4f] | coherence %.3f  %+.3f [%+.3f, %+.3f]",
                x$method[k], x$aupr[k], x$diff_vs_clr[k], x$diff_lo[k], x$diff_hi[k],
                y$median_auroc[k], y$median_diff_vs_clr[k], y$diff_lo[k], y$diff_hi[k]))
}
say("done; outputs in ", normalizePath(opt$out))
