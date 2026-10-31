#!/usr/bin/env Rscript
# Re-scores one baselines run under three tie conventions for pooled AUPR
# (old = ap_fast, expected under random tie-breaking, optimistic bound, plus the
# negatives-first strict lower bound); see analysis/tie_check/ap_ties.R.
# Data, truths and delete-half gene splits are built by running the SAME code
# as analysis/baselines.R (its lines are parsed and evaluated, not copied), so
# truths and splits are identical to the run being re-scored. GENIE3, PIDC and
# GRNBoost2 scores are read from the run's cached *_scores.rds; |r|, CLR,
# attention and the hub score are not cached by baselines.R, so they are
# recomputed here exactly as baselines.R does (deterministic; the old-convention
# AUPR is checked against the run's metrics.csv).
#
# Usage (repo root):
#   Rscript analysis/tie_check/tie_check.R --run baselines_m3d --dataset m3d [--gene_subset none]
#       [--workers 16] [--out results/tie_check]
args <- commandArgs(trailingOnly = TRUE)
o <- list(run = NULL, dataset = NULL, gene_subset = "none", workers = 16L, out = "results/tie_check",
          t = 9L, alpha = 0.5, seed = 20261016L, jack = 100L, quick = 0L, jack_workers = 4L)
i <- 1L
while (i <= length(args)) { key <- sub("^--", "", args[i]); if (!key %in% names(o)) stop("unknown argument ", args[i])
  v <- args[i + 1L]; o[[key]] <- if (key %in% c("workers", "jack_workers", "t", "seed", "jack", "quick")) as.integer(v) else if (key == "alpha") as.numeric(v) else v; i <- i + 2L }
stopifnot(!is.null(o$run), !is.null(o$dataset))
run_dir <- file.path("results", o$run)
part_dir <- file.path(o$out, "parts", o$run); dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
source("analysis/tie_check/ap_ties.R")

## ---- data, truths, splits: baselines.R lines, evaluated verbatim --------------------------
bl <- readLines("analysis/baselines.R")
i1 <- grep('^say\\(sprintf\\("%s: %d genes x %d samples"', bl)
i2 <- grep("^set.seed\\(opt\\$seed\\)$", bl); i3 <- grep("^# gene-jackknife replicates", bl)
stopifnot(length(i1) == 1, length(i2) == 1, length(i3) == 1, i2 < i3)
bl_args <- c("--dataset", o$dataset, "--gene_subset", o$gene_subset, "--out", file.path(o$out, "scratch", o$run),
             "--seed", o$seed, "--jack", o$jack, "--quick", o$quick, "--workers", o$workers, "--t", o$t, "--alpha", o$alpha)
commandArgs <- function(trailingOnly = FALSE) as.character(bl_args)
eval(parse(text = bl[seq_len(i1)]), envir = globalenv())
eval(parse(text = bl[i2:(i3 - 1L)]), envir = globalenv())
rm(commandArgs)
cat(sprintf("[%s] %s: %d genes, truths: %s\n", format(Sys.time(), "%H:%M:%S"), o$run, G, paste(names(bench), collapse = ",")))

## ---- scoring under the conventions -------------------------------------------------------
convs <- c("old", "expected", "optimistic", "negfirst")
score_all <- function(S, tag) {
  S[!is.finite(S)] <- 0; diag(S) <- 0
  res <- list()
  for (nm in names(bench)) {
    b <- bench[[nm]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]; lab <- b$pairs$label
    pooled <- ap_ties(sc, lab); ts <- tie_summary(sc, lab)
    pp <- ap_prepare(sc, lab); rm(sc); invisible(gc())
    jk <- parallel::mclapply(halves[[nm]], function(ix) { m <- logical(length(pp$o)); m[ix] <- TRUE; ap_ties_subset(pp, m) },
                             mc.cores = o$jack_workers)
    jk <- do.call(rbind, jk); stopifnot(nrow(jk) == o$jack)
    res[[nm]] <- list(pooled = pooled, ties = ts, jack = jk, base_rate = mean(lab)) }
  # matrix-level tie fraction over all off-diagonal gene pairs (upper triangle)
  ut <- S[upper.tri(S)]; os <- sort.int(ut, method = "radix"); last <- c(os[-1L] != os[-length(os)], TRUE)
  m <- diff(c(0L, which(last)))
  res$.matrix <- c(n_pairs = length(ut), n_tied_pairs = sum(m[m > 1]), max_block = max(m),
                   frac_tied = sum(m[m > 1]) / length(ut), frac_largest_block = max(m) / length(ut),
                   largest_block_value = os[which(last)[which.max(m)]])
  res }
hub_stratify <- function(S, H, k = 10L) {      # verbatim from baselines.R
  S[!is.finite(S)] <- 0; ut <- which(upper.tri(S)); h <- H[ut]
  g <- cut(h, unique(stats::quantile(h, seq(0, 1, length.out = k + 1L))), include.lowest = TRUE)
  r <- stats::ave(S[ut], g, FUN = function(x) rank(x) / length(x))
  R <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S)); R[ut] <- r; R + t(R) }
stopifnot(identical(deparse(hub_stratify), deparse(local({ e <- new.env(); eval(parse(text = bl[grep("^hub_stratify <- ", bl):(grep("^hub_stratify <- ", bl) + 5L)]), envir = e); e$hub_stratify }))))

# check against the run's own metrics.csv (old convention must reproduce it)
met <- utils::read.csv(file.path(run_dir, "metrics.csv"), stringsAsFactors = FALSE)
store <- function(S, m, H) {
  t0 <- Sys.time()
  r <- score_all(S, m)
  rh <- if (!is.null(H)) score_all(hub_stratify(S, H), paste0(m, "_hubstrat")) else NULL
  for (v in c(m, if (!is.null(rh)) paste0(m, "_hubstrat"))) {
    rr <- if (v == m) r else rh
    for (nm in names(bench)) { ref <- met$aupr[met$method == v & met$truth == nm]
      if (length(ref) == 1) { d <- abs(rr[[nm]]$pooled[["old"]] - ref)
        cat(sprintf("   check %-24s %-12s old-convention AUPR %.6f vs metrics.csv %.6f (diff %.1e)\n", v, nm, rr[[nm]]$pooled[["old"]], ref, d))
        if (d > 1e-9) warning("MISMATCH vs metrics.csv for ", v, " ", nm) } } }
  saveRDS(list(method = m, plain = r, hubstrat = rh), file.path(part_dir, paste0(m, ".rds")))
  cat(sprintf("[%s] %s scored in %.1f min\n", format(Sys.time(), "%H:%M:%S"), m, as.numeric(difftime(Sys.time(), t0, units = "mins")))) }

## ---- references (as baselines.R) ---------------------------------------------------------
suppressPackageStartupMessages({ library(clr); library(Matrix) })
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = o$workers)$
  calibrate(method = "normal", combine = "stouffer")
clr_s <- f$clr_scores
f$build_operator(alpha = o$alpha, softmax_keff = 10, softmax_cap = 50L)
P <- Matrix::Diagonal(G, 1 - o$alpha) + o$alpha * Matrix::Matrix(f$operator, sparse = TRUE)
Pt <- as.matrix(P); for (s in seq_len(o$t - 1L)) Pt <- as.matrix(Pt %*% P)
St <- (Pt + t(Pt)) / 2; rm(Pt)
pi <- rep(1 / G, G); eps <- 1e-3
for (it in 1:20000) { q <- (1 - eps) * as.numeric(pi %*% P) + eps / G; if (sum(abs(q - pi)) < 1e-13) break; pi <- q }
H <- outer(q, q, "+")
rm(P); f$release(); invisible(gc())
cat(sprintf("[%s] references built\n", format(Sys.time(), "%H:%M:%S")))
todo <- function(m) !file.exists(file.path(part_dir, paste0(m, ".rds")))
if (todo("attention")) store(St, "attention", H)
if (todo("clr")) store(clr_s, "clr", H)
if (todo("hub")) store(H, "hub", NULL)
if (todo("abs_pearson")) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; store(r, "abs_pearson", H); rm(r) }
rm(St, clr_s); invisible(gc())
for (m in c("genie3", "pidc", "grnboost2")) {
  rdsf <- file.path(run_dir, paste0(m, "_scores.rds"))
  if (!file.exists(rdsf) || !todo(m)) next
  S <- readRDS(rdsf); stopifnot(nrow(S) == G); store(S, m, H); rm(S); invisible(gc()) }
cat("done\n")
