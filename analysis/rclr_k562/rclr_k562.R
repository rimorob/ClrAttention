#!/usr/bin/env Rscript
# EXPLORATORY PATHFINDING, NOT PREREGISTERED, NOTHING CONFIRMATORY. Nothing from this
# script goes into paper/ or any paper table.
#
# Question: does rCLR (robust CLR; docs/plans/background_width_exploration.md idea 2 and
# section 4.3; the definition is eb_additive() in analysis/synthetic_width/sim_width.R,
# copied below) change CORUM performance on real K562 essential-gene Perturb-seq the way
# the synthetic work predicted (pairwise score: ribosome up ~4x, small complexes down
# ~0.55x; as attention input: inconsistent)?
#
# Dataset caveat: K562 CORUM was already used in the preregistered keff test
# (analysis/keff_test.R) and in analysis/corum_exploratory/, so this is NOT a held-out
# dataset. Same data (half A, 5,671 genes x 707 knockdown profiles), same CORUM
# (corum_sets), same operator pipeline (cap 50, lazy walk alpha 0.5, t = 9, symmetrised
# S = (P^t + P^t')/2). The cached full-data CLR (results/keff_test/k562ess/cache) is reused
# and checked against a recomputation; MI is recomputed only because the cache stores
# CLR, not MI (rCLR needs MI).
#
# Operators (fixed before running): CLR (reference), rCLR (pairwise), K10 on CLR (the
# reference ratio denominator), K10 on rCLR, K25 on CLR, K25 on rCLR.
# Endpoints (fixed before running; no endpoint was added after seeing results):
#   all CORUM; CORUM excluding ribosome; ribosome alone; small complexes (3-6 genes);
#   large complexes (>= 20 genes). "Ribosome" = the 8 RP-dominated complexes (>= 30% RP
#   genes) of analysis/corum_exploratory/corum_complexes.csv, the same definition as the
#   earlier CORUM_noRP exclusion. Complex size = number of its members present in the data.
#   For each endpoint the pair universe is all gene pairs within the genes of the selected
#   complexes and a pair is positive if it co-occurs in a selected complex
#   (comembership_pairs on the subset, as CORUM_noRP did). Metric: pooled co-membership
#   AUPR and enrichment = AUPR / base rate; reported as ratios to K10 on CLR. A descriptive
#   half-of-genes jackknife (100 reps, seed 20261016) gives a spread for the ratios.
# Label-free side check: fraction of genes with a 'reproducible partner' (shared top-50
#   positive candidates between two halves of the knockdown profiles; the split and rule of
#   analysis/synthetic_width/m3d_label_free.R k562 mode), CLR vs rCLR.
#
# Usage (repo root): Rscript analysis/rclr_k562/rclr_k562.R --out results/rclr_k562 [--quick N] [--threads 8]
args <- commandArgs(TRUE)
opt <- list(out = "results/rclr_k562", quick = 0L, threads = 8L, jack = 100L)
for (k in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[k])]] <- if (is.numeric(opt[[sub("^--", "", args[k])]])) as.integer(args[k + 1]) else args[k + 1]
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
QUICK <- opt$quick > 0L

## ---- data (K562 half A pseudobulk, as keff_test.R / m3d_label_free.R k562 mode) ----------
raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (QUICK) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
genes <- rownames(X); G <- length(genes)
say(sprintf("%d genes x %d profiles (quick = %d)", G, ncol(X), opt$quick))

## ---- operators ------------------------------------------------------------------------
fit_of <- function(Xm) { f <- ClrAttention$new(Xm)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
  f$calibrate(method = "normal", combine = "stouffer"); f }
# rCLR: eb_additive() of analysis/synthetic_width/sim_width.R (U only), G from the matrix.
rclr_of <- function(M) {
  G <- nrow(M); Y <- M^(1 / 3); diag(Y) <- NA
  mu <- stats::median(Y, na.rm = TRUE); a <- rep(0, G)
  for (it in 1:8) { R <- Y - mu - outer(a, a, "+"); da <- apply(R, 1, stats::median, na.rm = TRUE)
    a <- a + da; mu <- mu + stats::median(a); a <- a - stats::median(a) }
  R <- Y - mu - outer(a, a, "+"); s <- apply(R, 1, stats::mad, na.rm = TRUE)
  U <- (sweep(R, 1, s, "/") + sweep(R, 2, s, "/")) / sqrt(2); diag(U) <- NA; U[is.na(U)] <- 0; U
}
selfloop <- function(A) { e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1; A }
lazy_P <- function(A, alpha = 0.5) Matrix::Diagonal(nrow(A), 1 - alpha) + alpha * Matrix::Matrix(A, sparse = TRUE)
walk9 <- function(A) { P <- lazy_P(A); Pt <- as.matrix(P)
  for (t in 2:9) Pt <- as.matrix(Pt %*% P)
  S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }
attn <- function(Cm, k) walk9(selfloop(clr:::.softmax_rows(Cm, keff = k, cap = 50L)))

t0 <- Sys.time()
f <- fit_of(X); M <- f$mi; C <- f$clr_scores; f$release(); rm(f)
say(sprintf("MI + CLR: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cache_f <- "results/keff_test/k562ess/cache/fulldata.rds"
chk <- data.frame(check = character(), value = character())
if (!QUICK) {
  full <- readRDS(cache_f)
  chk <- rbind(chk, data.frame(check = c("genes md5 equals keff_test fingerprint", "max |CLR recomputed - cached CLR|"),
    value = c(as.character(identical(digest::digest(genes, algo = "md5"), full$fingerprint$genes_md5)),
              format(max(abs(C - full$C)), digits = 3))))
  C <- full$C; rm(full); invisible(gc(verbose = FALSE))   # use the cached CLR itself as the reference
}
t0 <- Sys.time(); U <- rclr_of(M); rm(M); invisible(gc(verbose = FALSE))
say(sprintf("rCLR: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
chk <- rbind(chk, data.frame(check = "cor(CLR, rCLR) over upper-triangle pairs", value = format(stats::cor(C[upper.tri(C)], U[upper.tri(U)]), digits = 4)))

## ---- truth: CORUM and the five fixed endpoints ----------------------------------------------
cr <- corum_sets(genes)
cx <- utils::read.csv("analysis/corum_exploratory/corum_complexes.csv", stringsAsFactors = FALSE)
sz <- lengths(cr$members); nm <- names(cr$members)
rp <- nm %in% cx$complex[cx$group == "RP-dominated"]
if (!QUICK) stopifnot(sum(rp) == 8L, all(sz == cx$size[match(nm, cx$complex)]))
sub_of <- function(keep) list(members = cr$members[keep], class = cr$class[keep])
EP <- list(`all CORUM` = rep(TRUE, length(nm)), `CORUM excl. ribosome` = !rp, `ribosome only` = rp,
           `small (3-6 genes)` = sz >= 3 & sz <= 6, `large (>=20 genes)` = sz >= 20)
none <- stats::setNames(character(0), character(0))
EP <- EP[vapply(EP, sum, 0L) > 0L]
bench <- lapply(EP, function(k) { r <- sub_of(k); list(n_cx = sum(k), pairs = comembership_pairs(r, none, genes)) })
for (e in names(bench)) say(sprintf("endpoint %-22s %4d complexes, %8d pairs, base rate %.4f", e, bench[[e]]$n_cx,
  length(bench[[e]]$pairs$label), mean(bench[[e]]$pairs$label)))
big_rp <- sum(rp & sz >= 20)

## ---- score every operator (one dense matrix at a time) --------------------------------------
set.seed(20261016L)                                       # baselines.R / keff_test.R jackknife seed
halves <- lapply(bench, function(b) lapply(seq_len(opt$jack), function(r) {
  keep <- sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2)); inh <- logical(G); inh[keep] <- TRUE
  which(inh[b$pairs$i] & inh[b$pairs$j]) }))
SC <- list()     # operator -> endpoint -> pair scores
grab <- function(op, S) { SC[[op]] <<- lapply(bench, function(b) S[cbind(b$pairs$i, b$pairs$j)]); say("scored ", op) }
grab("CLR", C); grab("rCLR", U)
for (k in c(10, 25)) { grab(sprintf("K%d on CLR", k), attn(C, k)); grab(sprintf("K%d on rCLR", k), attn(U, k)); invisible(gc(verbose = FALSE)) }
ops <- c("CLR", "rCLR", "K10 on CLR", "K10 on rCLR", "K25 on CLR", "K25 on rCLR")
stopifnot(identical(ops, names(SC)[match(ops, names(SC))]))

met <- do.call(rbind, lapply(names(bench), function(e) { b <- bench[[e]]; br <- mean(b$pairs$label)
  do.call(rbind, lapply(ops, function(o) { ap <- ap_fast(SC[[o]][[e]], b$pairs$label)
    data.frame(endpoint = e, n_complexes = b$n_cx, n_genes = length(b$pairs$genes_in), n_pairs = length(b$pairs$label),
               base_rate = br, operator = o, aupr = ap, enrichment = ap / br) })) }))
ref <- met[met$operator == "K10 on CLR", ]
met$ratio_to_K10_CLR <- met$aupr / ref$aupr[match(met$endpoint, ref$endpoint)]
ap_pc <- function(sc, lab) { sc[!is.finite(sc)] <- 0; ap_fast(sc, lab) }
nc <- max(1L, min(opt$threads, parallel::detectCores() - 1L))
jk <- do.call(rbind, lapply(names(bench), function(e) { lab <- bench[[e]]$pairs$label
  A <- do.call(cbind, lapply(ops, function(o) unlist(parallel::mclapply(halves[[e]], function(ix) ap_pc(SC[[o]][[e]][ix], lab[ix]), mc.cores = nc))))
  do.call(rbind, lapply(ops[-3], function(o) { q <- A[, match(o, ops)] / A[, 3]; q <- q[is.finite(q)]
    data.frame(endpoint = e, operator = o, jack_median = stats::median(q), jack_lo = unname(stats::quantile(q, 0.025)),
               jack_hi = unname(stats::quantile(q, 0.975)), n_rep = length(q)) })) }))
say("jackknife done")
met <- merge(met, jk, by = c("endpoint", "operator"), all.x = TRUE, sort = FALSE)
met <- met[order(match(met$endpoint, names(bench)), match(met$operator, ops)), ]
utils::write.csv(met, file.path(opt$out, "endpoint_metrics.csv"), row.names = FALSE)

## ---- label-free side check: reproducible partners, CLR vs rCLR --------------------------------
top50 <- function(Cm) lapply(seq_len(nrow(Cm)), function(i) { r <- Cm[i, ]; r[i] <- -Inf
  o <- order(r, decreasing = TRUE)[1:min(50L, nrow(Cm) - 1L)]; o[r[o] > 0] })
set.seed(20261030L); lab <- colnames(X)                   # m3d_label_free.R k562 mode: each profile its own 'lab'
sz0 <- table(lab); sz0 <- sz0[sample(length(sz0))]
half <- stats::setNames(integer(length(sz0)), names(sz0)); tot <- c(0, 0)
for (l in names(sz0)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz0[[l]] }
h <- unname(half[lab])
fh <- lapply(1:2, function(d) { ff <- fit_of(X[, h == d]); Mh <- ff$mi; Ch <- ff$clr_scores; ff$release(); Uh <- rclr_of(Mh); rm(Mh)
  list(clr = top50(Ch), rclr = top50(Uh)) })
rn <- function(a, b) vapply(seq_len(G), function(i) length(intersect(a[[i]], b[[i]])), 0L)
rep_clr <- rn(fh[[1]]$clr, fh[[2]]$clr); rep_rclr <- rn(fh[[1]]$rclr, fh[[2]]$rclr)
lam <- 50^2 / (G - 1)   # chance overlap; k0 = smallest k with P(Pois(lam) >= k) < 0.05
k0 <- (function() { for (k in 1:100) if (1 - stats::ppois(k - 1, lam) < 0.05) return(k) })()
lf <- do.call(rbind, lapply(list(CLR = rep_clr, rCLR = rep_rclr), function(v)
  data.frame(frac_zero_shared = mean(v == 0), frac_ge_k0 = mean(v >= k0), k0 = k0, median_shared = stats::median(v),
             q10 = unname(stats::quantile(v, .1)), q90 = unname(stats::quantile(v, .9)), mean_shared = mean(v))))
lf <- cbind(operator = rownames(lf), lf, n_genes = G, n_profiles_half = paste(sum(h == 1), sum(h == 2), sep = "/"))
lf$spearman_CLR_vs_rCLR_repro_n <- stats::cor(rep_clr, rep_rclr, method = "spearman")
utils::write.csv(lf, file.path(opt$out, "label_free_repro.csv"), row.names = FALSE)
utils::write.csv(chk, file.path(opt$out, "sanity_checks.csv"), row.names = FALSE)
print(lf, row.names = FALSE)

## ---- summary.md ------------------------------------------------------------------------------------
mt <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
  apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
rt <- do.call(rbind, lapply(names(bench), function(e) { d <- met[met$endpoint == e, ]
  data.frame(endpoint = sprintf("%s (base %.4f)", e, d$base_rate[1]), t(stats::setNames(
    sprintf("%.2f", d$ratio_to_K10_CLR), d$operator))[, ops, drop = FALSE], check.names = FALSE) }))
en <- do.call(rbind, lapply(names(bench), function(e) { d <- met[met$endpoint == e, ]
  data.frame(endpoint = e, t(stats::setNames(sprintf("%.2f", d$enrichment), d$operator))[, ops, drop = FALSE], check.names = FALSE) }))
jr <- met[met$operator != "K10 on CLR", ]; jr$jackknife_ratio <- sprintf("%.2f [%.2f, %.2f]", jr$jack_median, jr$jack_lo, jr$jack_hi)
jt <- do.call(rbind, lapply(names(bench), function(e) { d <- jr[jr$endpoint == e, ]
  data.frame(endpoint = e, t(stats::setNames(d$jackknife_ratio, d$operator))[, ops[-3], drop = FALSE], check.names = FALSE) }))
md <- c("# rCLR on K562 essential Perturb-seq, scored on CORUM (EXPLORATORY)", "",
  "Exploratory pathfinding, not preregistered, not held-out (K562 CORUM was used in the keff test), nothing confirmatory; not for paper/ or paper tables. Script: `analysis/rclr_k562/rclr_k562.R`. Endpoint list and operator list fixed before running.", "",
  sprintf("Data: %d genes x %d knockdown profiles (half A). Ribosome = the 8 RP-dominated CORUM complexes (%d of them have >= 20 genes). Metric: pooled co-membership AUPR; enrichment = AUPR / base rate.", G, ncol(X), big_rp), "",
  "## Ratio of AUPR (= enrichment ratio) to K10 on CLR", "", mt(rt), "",
  "## Enrichment over base rate", "", mt(en), "",
  "## Half-of-genes jackknife of the ratio to K10 on CLR (median [2.5%, 97.5%], 100 reps; descriptive)", "", mt(jt), "",
  "## Label-free: reproducible partners (top-50 positive candidates shared by two halves of the profiles)", "", mt(lf[, c("operator", "frac_zero_shared", "frac_ge_k0", "k0", "median_shared", "q10", "q90")]), "",
  "## Sanity checks", "", mt(chk), "")
writeLines(md, file.path(opt$out, "summary.md"))
say("done")
