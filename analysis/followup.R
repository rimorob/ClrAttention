#!/usr/bin/env Rscript
# Follow-up analysis on the 466-experiment M3D compendium x RegulonDB:
#   A. Confidence intervals: paired delete-half jackknife over annotated genes
#      for pooled co-membership AUPR (and differences vs single-step CLR), and a
#      paired bootstrap over regulons for per-regulon coherence AUROC.
#   B. Held-out depth selection: experiments are split into two halves by
#      experimenter (lab), so related conditions never straddle the split. The
#      operator is built on the TRAIN half; depth t is chosen to maximize how well
#      attention mass at depth t predicts the TEST half's strongest CLR edges
#      (no RegulonDB involved). The search is a geometric scan (t = 1, 2, 4, ...,
#      tmax) followed by golden-section refinement on integer t inside the
#      bracket around the best scan point. The chosen t* is then scored against
#      RegulonDB on the full data, together with fixed geometric depths.
# Reduced design (user, 2026-09-23): configurations |Pearson|, raw MI,
# parity-2007 CLR, HG-Stouffer CLR; operators soft10, fdr05_top10 and the
# pearson_top10 control; raw values; readout = attention mass.
#
# Usage (repo root):
#   Rscript analysis/followup.R [--rds results/avg/primary_fit.rds]
#       [--m3d data/E_coli_v4_Build_6] [--rdb data/RegulonDBExtract]
#       [--jack 100] [--boot 2000] [--tmax 128] [--Btrain 50] [--threads N]
#       [--out results/avg_followup] [--quick 0]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(rds = "results/avg/primary_fit.rds", m3d = "data/E_coli_v4_Build_6",
            rdb = "data/RegulonDBExtract", jack = 100L, boot = 2000L,
            tmax = 128L, Btrain = 50L, threads = NULL, out = "results/avg_followup",
            quick = 0L, seed = 20260923L, alpha = 0.5, target_frac = 0.01)
ints <- c("jack", "boot", "tmax", "Btrain", "threads", "quick", "seed")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else
    if (key %in% c("alpha", "target_frac")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "followup.log")
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

## ---- data ------------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
raw <- utils::read.delim(file.path(opt$m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                         check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (opt$quick > 0L) {
  v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]
}
genes <- rownames(X); G <- length(genes)
say("compendium: ", G, " genes x ", ncol(X), " experiments")
to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)
bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
  ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
  reg <- make_regulons(rtab, genes, conf_keep = ck)
  list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes),
       co = comember_matrix(reg, G))
})
Zs <- function(M) { Z <- M - rowMeans(M); Z / sqrt(rowSums(Z^2) / (ncol(Z) - 1)) }

## ---- full-data similarity matrices ---------------------------------------
build_fit <- function(Xm, threads = opt$threads) {
  f <- ClrAttention$new(Xm)$estimate_mi(bins = "hg", transform = "none",
                                        threads = threads)
  f$calibrate(method = "normal", combine = "stouffer")
  f
}
fit <- build_fit(X)
Z <- Zs(X); R0 <- tcrossprod(Z) / (ncol(Z) - 1)
if (file.exists(opt$rds) && opt$quick == 0L) {
  rds <- readRDS(opt$rds)
  stopifnot(identical(rds$genes, genes))
  if (!is.null(rds$null)) {
    fit$restore_threshold(rds$null, rds$params$threshold$tau,
                          rds$params$threshold$q %||% 0.05)
    say("restored permutation null from ", opt$rds)
  }
  rm(rds); invisible(gc())
}
if (is.null(tryCatch(fit$null_distribution, error = function(e) NULL))) {
  say("no stored null: running select_threshold(B = ", opt$Btrain, ") on full data")
  fit$select_threshold(B = opt$Btrain, threads = opt$threads)
}
M0 <- fit$mi; diag(M0) <- 0
S_clr <- fit$clr_scores
parity <- ClrAttention$new(X)$estimate_mi(bins = 10, transform = "none",
                                          threads = opt$threads)$
  calibrate(method = "normal", combine = "euclidean")$clr_scores

make_operator <- function(f, R, name) {
  if (name == "soft10") f$build_operator(alpha = opt$alpha, softmax_keff = 10,
                                         softmax_cap = 50L)
  else if (name == "fdr05_top10") { f$reselect(0.05)
    f$build_operator(alpha = opt$alpha, topk_union = 10L) }
  else if (name == "pearson_top10") {
    Ar <- abs(R); diag(Ar) <- 0
    f$set_operator(clr:::.topk_rows(Ar, 10L), alpha = opt$alpha, label = name)
  }
  A <- f$operator
  Matrix::Matrix((1 - opt$alpha) * diag(nrow(A)) + opt$alpha * A, sparse = TRUE)
}
ops <- c("soft10", "fdr05_top10", "pearson_top10")
grid <- 2L^(0:floor(log2(opt$tmax)))

# Attention-mass similarity at a set of depths (sequential sparse steps).
attention_at <- function(P, ts) {
  out <- list(); Pt <- diag(nrow(P))
  for (t in seq_len(max(ts))) {
    Pt <- as.matrix(Pt %*% P)
    if (t %in% ts) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; out[[as.character(t)]] <- S }
  }
  out
}

## ---- B. held-out depth selection -------------------------------------------
say("B. held-out depth selection")
desc <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_descriptions")),
                          stringsAsFactors = FALSE)
exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description)
exper[!grepl("experimenter:", desc$description)] <- NA
lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])
lab <- lab[colnames(X)]
lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab)))
sz <- sort(table(lab), decreasing = TRUE)
half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
split_h <- half[lab]
say(sprintf("split by experimenter: %d labs; half 1 = %d experiments, half 2 = %d",
            length(sz), sum(split_h == 1), sum(split_h == 2)))
wcsv(data.frame(experiment = colnames(X), lab = unname(lab), half = unname(split_h)),
     "heldout_split.csv")

heldout_rows <- list(); tstar <- list()
for (dir in 1:2) {
  tr <- colnames(X)[split_h == dir]; te <- colnames(X)[split_h != dir]
  ftr <- build_fit(X[, tr]); fte <- build_fit(X[, te])
  ut <- upper.tri(fte$clr_scores)
  sv <- fte$clr_scores[ut]
  target <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac))
  say(sprintf("direction %d: train %d / test %d experiments; target = top %.1f%% test-half CLR pairs (%d)",
              dir, length(tr), length(te), 100 * opt$target_frac, sum(target)))
  if ("fdr05_top10" %in% ops) {
    set.seed(opt$seed + dir)
    ftr$select_threshold(B = opt$Btrain, threads = opt$threads)
  }
  Rtr <- tcrossprod(Zs(X[, tr])) / (length(tr) - 1)
  crit <- function(S) pr_summary(S[ut], target)[["aupr"]]
  for (op in ops) {
    P <- make_operator(ftr, Rtr, op)
    memo <- list()
    ck <- list("0" = diag(G))          # checkpoints of P^t (grid depths)
    evalt <- function(ts) {
      for (t in sort(setdiff(ts, as.integer(names(memo))))) {
        t0 <- max(as.integer(names(ck))[as.integer(names(ck)) <= t])
        Pt <- ck[[as.character(t0)]]
        if (t > t0) for (s in seq_len(t - t0)) Pt <- as.matrix(Pt %*% P)
        if (t %in% grid) ck[[as.character(t)]] <<- Pt
        S <- (Pt + t(Pt)) / 2; diag(S) <- 0
        memo[[as.character(t)]] <<- crit(S)
      }
      vapply(as.character(ts), function(k) memo[[k]], numeric(1))
    }
    g <- evalt(grid)
    kbest <- which.max(g)
    lo <- grid[max(1, kbest - 1)]; hi <- grid[min(length(grid), kbest + 1)]
    # golden-section on integers in [lo, hi]
    phi <- (sqrt(5) - 1) / 2
    a <- lo; b <- hi
    while (b - a > 2) {
      c1 <- round(b - phi * (b - a)); c2 <- round(a + phi * (b - a))
      if (c1 == c2) c2 <- c1 + 1L
      f <- evalt(c(c1, c2))
      if (f[1] >= f[2]) b <- c2 else a <- c1
    }
    cand <- a:b; fv <- evalt(cand)
    tb <- cand[which.max(fv)]
    tstar[[op]] <- c(tstar[[op]], tb)
    say(sprintf("  %-14s geometric scan %s -> bracket [%d, %d]; t* = %d (criterion %.4f); %d depths evaluated",
                op, paste(sprintf("%d:%.4f", grid, g), collapse = " "), lo, hi, tb,
                max(fv), length(memo)))
    heldout_rows[[length(heldout_rows) + 1L]] <- data.frame(
      direction = dir, operator = op, t = as.integer(names(memo)),
      criterion = unlist(memo))
    rm(ck); invisible(gc())
  }
  rm(ftr, fte); invisible(gc())
}
ho <- do.call(rbind, heldout_rows); wcsv(ho, "heldout_criterion.csv")
tsel <- vapply(tstar, function(v) as.integer(round(exp(mean(log(v))))), integer(1))
say("held-out t* (geometric mean over the two directions): ",
    paste(names(tsel), tsel, sep = "=", collapse = " "))
wcsv(data.frame(operator = names(tsel), t_dir1 = vapply(tstar, `[`, 0L, 1),
                t_dir2 = vapply(tstar, `[`, 0L, 2), t_star = tsel), "heldout_tstar.csv")

## ---- full-data methods at geometric depths and t* -------------------------
say("scoring methods on the full data")
methods <- list(abs_pearson = abs(R0), mi = M0, clr_parity2007 = parity,
                clr_hg_stouffer = S_clr)
for (op in ops) {
  P <- make_operator(fit, R0, op)
  ts <- sort(unique(c(grid, tsel[[op]])))
  A <- attention_at(P, ts)
  for (k in names(A)) {
    lab_k <- sprintf("%s:t%03d%s", op, as.integer(k),
                     if (as.integer(k) == tsel[[op]]) "*" else "")
    methods[[lab_k]] <- A[[k]]
  }
  rm(A); invisible(gc())
}
# fdr20 with the diffused-profile readout |cor(P^t Z)| (its best readout in
# the main run), at the geometric depths, for the confidence-interval table.
fit$reselect(0.20); fit$build_operator(alpha = opt$alpha)
P20 <- Matrix::Matrix((1 - opt$alpha) * diag(G) + opt$alpha * fit$operator, sparse = TRUE)
Et <- Z
for (t in seq_len(max(grid))) {
  Et <- as.matrix(P20 %*% Et)
  if (t %in% grid) methods[[sprintf("fdr20_abscor:t%03d", t)]] <- abs(stats::cor(t(Et)))
}
diag_zero <- function(S) { diag(S) <- 0; S }
methods <- lapply(methods, diag_zero)
fit$reselect(0.05)
say("methods: ", length(methods))

## ---- A. confidence intervals ----------------------------------------------
say("A. delete-half jackknife over genes (", opt$jack, " half-samples) and regulon bootstrap")
pooled <- list(); jack <- list(); coh <- list()
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
    J[r, ] <- vapply(sc, function(s) pr_summary(s[sel], p$label[sel])[["aupr"]], numeric(1))
  }
  # delete-d jackknife, d = n/2: var(theta_hat) ~ var over half-samples
  se <- apply(J, 2, stats::sd)
  dref <- J - J[, "clr_hg_stouffer"]
  se_d <- apply(dref, 2, stats::sd)
  pooled[[ev]] <- data.frame(evidence = ev, method = names(methods), aupr = full,
                             se = se, lo = full - 1.96 * se, hi = full + 1.96 * se,
                             diff_vs_clr = full - full[["clr_hg_stouffer"]],
                             diff_se = se_d,
                             diff_lo = full - full[["clr_hg_stouffer"]] - 1.96 * se_d,
                             diff_hi = full - full[["clr_hg_stouffer"]] + 1.96 * se_d,
                             base_rate = mean(p$label))
  # per-regulon coherence: paired bootstrap over regulons
  b <- bench[[ev]]
  Cm <- vapply(methods, function(S) regulon_coherence(S, b$reg, b$co, operon_of, genes),
               numeric(length(b$reg$members)))
  ok <- stats::complete.cases(Cm); Cm <- Cm[ok, , drop = FALSE]
  cls <- b$reg$class[ok]
  dC <- Cm - Cm[, "clr_hg_stouffer"]
  bt <- replicate(opt$boot, { ix <- sample(nrow(dC), replace = TRUE)
                              apply(dC[ix, , drop = FALSE], 2, stats::median) })
  coh[[ev]] <- data.frame(evidence = ev, method = names(methods),
                          median_auroc = apply(Cm, 2, stats::median),
                          median_diff_vs_clr = apply(dC, 2, stats::median),
                          diff_lo = apply(bt, 1, stats::quantile, 0.025),
                          diff_hi = apply(bt, 1, stats::quantile, 0.975),
                          frac_regulons_improved = colMeans(dC > 0),
                          wilcoxon_p = apply(dC, 2, function(d)
                            if (all(d == 0)) NA_real_ else
                              suppressWarnings(stats::wilcox.test(d)$p.value)),
                          n_regulons = nrow(dC))
  wcsv(data.frame(regulon = rownames(Cm) %||% names(b$reg$members)[ok],
                  class = cls, size = b$reg$size[ok], Cm, check.names = FALSE),
       sprintf("per_regulon_coherence_%s.csv", ev))
  say(sprintf("  [%s] done", ev))
}
pooled <- do.call(rbind, pooled); coh <- do.call(rbind, coh)
wcsv(pooled, "ci_comembership.csv"); wcsv(coh, "ci_coherence.csv")
fmt <- function(d) sprintf("%-24s AUPR %.4f [%.4f, %.4f]  vs CLR %+.4f [%+.4f, %+.4f]",
                           d$method, d$aupr, d$lo, d$hi, d$diff_vs_clr, d$diff_lo, d$diff_hi)
for (ev in c("SC", "all")) {
  say("pooled co-membership AUPR [", ev, "] (base ", sprintf("%.4f", pooled$base_rate[pooled$evidence == ev][1]), "):")
  for (l in fmt(pooled[pooled$evidence == ev, ])) say("  ", l)
}
say("done; outputs in ", normalizePath(opt$out))
