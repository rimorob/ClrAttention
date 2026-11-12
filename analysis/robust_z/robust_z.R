#!/usr/bin/env Rscript
# EXPLORATORY, NOT PREREGISTERED, NOTHING CONFIRMATORY. Nothing from this script goes into
# paper/ or any paper table. K562 CORUM is not held out (it was used in the keff test and in
# analysis/corum_exploratory/). No p-values are computed; intervals and point estimates only.
#
# Question (the PI's idea): CLR z-scores each gene's MI row by mean and SD. A gene whose own
# module has many strong partners inflates its row SD, which lowers the z-scores of its true
# partners (may hurt large complexes such as the ribosome). Does a robust row standardisation
# (median and IQR/1.349 instead of mean and SD) change pairwise or attention performance?
#
# Variants (analysis/robust_z/robust_z_lib.R): V0 mean/SD (default CLR; asserted equal to the
# package's CLR, operator and walk mass to 1e-8 before anything else), V0x mean/SD with the
# diagonal excluded from the row (control), V1 median/(IQR/1.349) of raw MI, V2 same on MI^(1/3),
# V3 the existing rCLR (reference only), V4 median/(1.4826 MAD) of raw MI (diagnostic).
# Everything after the row z (truncation at 0, symmetrisation, softmax keff 10 cap 50, lazy walk
# alpha 0.5, S_t = (P^t + P^t')/2) is the default code path.
#
# Datasets and truths as the project's scripts define them:
#   m3d    M3D 466 experiments, RegulonDB strong evidence (SC) and all evidence (analysis/baselines.R)
#   k562   K562 essential Perturb-seq pseudobulk half A, 5,671 truth genes
#          (results/baselines_k562ess_truthgenes_genie3/expression.tsv): CORUM all / excluding
#          ribosome / ribosome only / small (3-6 genes) / large (>= 20 genes) as
#          analysis/rclr_k562/rclr_k562.R, DoRothEA A-B, KD-response modules from half B.
# Read-outs per variant: pairwise score (the symmetrised z) and attention (t = 9, plus the
# label-free held-out depth t* of analysis/lab_splits.R, criterion computed on each variant's own
# score matrix, median over splits). Pooled AUPR / base rate with paired delete-half gene
# jackknife (100 reps, seed 20261016, as baselines.R); per-regulon (macro) AP with delete-20% gene
# jackknife (as baselines.R / macro_aupr.R, first jack_reg reps of the same stream); hub control:
# within-stratum percentile in 10 strata of H_i + H_j with H = pi (stationary distribution of the
# V0 attention walk; the stratifier of the sota table, attention-derived, the same for all
# variants) and, as a method-neutral second stratifier, H = row mean of raw MI; hub-drift
# diagnostic; row-scale diagnostics.
#
# Usage (repo root): Rscript analysis/robust_z/robust_z.R --dataset m3d|k562 --out results/robust_z/<dataset>
#   [--quick N] [--threads 8] [--workers 8] [--nsplit 8] [--tmax 20] [--jack 100] [--jack_reg 100]
args <- commandArgs(TRUE)
opt <- list(dataset = "m3d", out = NULL, quick = 0L, threads = 8L, workers = 8L, nsplit = 8L, tmax = 20L,
            jack = 100L, jack_reg = 100L, t = 9L, alpha = 0.5, seed = 20261016L, depth_seed = 20261014L,
            variants = "V0,V0x,V1,V2,V3,V4", top_pairs = 100000L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); if (!key %in% names(opt)) stop("unknown argument: ", args[k])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[k + 1]) else args[k + 1] }
for (k in c("quick", "threads", "workers", "nsplit", "tmax", "jack", "jack_reg", "t", "seed", "depth_seed", "top_pairs")) opt[[k]] <- as.integer(opt[[k]])
if (is.null(opt$out)) opt$out <- file.path("results/robust_z", opt$dataset)
dir.create(file.path(opt$out, "cache"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "robust_z.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
script_dir <- (function() { f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis/robust_z" })()
source(file.path(script_dir, "robust_z_lib.R"))
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
stopifnot(opt$dataset %in% c("m3d", "k562"))
VARS <- strsplit(opt$variants, ",")[[1]]
none <- stats::setNames(character(0), character(0))
rds_ok <- function(f) file.exists(f)

## ---- data ----------------------------------------------------------------------------------
if (opt$dataset == "m3d") {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
  exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description); exper[!grepl("experimenter:", desc$description)] <- NA
  lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])[colnames(X)]
  lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab))); rm(desc, exper)
} else {
  raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  lab <- colnames(X)                         # each knockdown profile is its own 'lab' (m3d_label_free.R k562 mode)
}
storage.mode(X) <- "double"
if (opt$quick > 0L) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
genes <- rownames(X); G <- length(genes)
say(sprintf("%s: %d genes x %d samples (%d labs)", opt$dataset, G, ncol(X), length(unique(lab))))
UT <- which(upper.tri(matrix(FALSE, G, G)))      # linear indices of the upper triangle

## ---- truths --------------------------------------------------------------------------------
macro_src <- list()      # macro sources: name -> list(reg, co, operon)
if (opt$dataset == "m3d") {
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn); operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck)
    list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes), co = comember_matrix(reg, G), operon = operon_of, macro = ev, sel = NULL) })
  for (nm in names(bench)) macro_src[[nm]] <- bench[[nm]][c("reg", "co", "operon")]
  base_truths <- names(bench); extra_truths <- character(0)
  modsize <- vapply(seq_len(G), function(g) { s <- lengths(bench$SC$reg$members)[vapply(bench$SC$reg$members, function(m) g %in% m, NA)]; if (length(s)) max(s) else NA_real_ }, 0)
  modsize_label <- "largest RegulonDB strong-evidence regulon containing the gene"
} else {
  h5 <- "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"     # half-B knockdown set, as baselines.R
  set.seed(20261015L)
  rd <- function(p) rhdf5::h5read(h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
    if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
    if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  pert <- cat_of("/obs", "gene_transcript"); pval <- as.numeric(rd("/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A_half <- sort(K[seq_len(floor(length(K) / 2))]); B_half <- setdiff(K, A_half)
  rm(rd, cat_of, pert, pval, K)
  tr <- list(DoRothEA_AB = dorothea_sets(genes), CORUM = corum_sets(genes), KD_response = response_sets(genes, path = h5, use = B_half))
  cr <- tr$CORUM; cx <- utils::read.csv("analysis/corum_exploratory/corum_complexes.csv", stringsAsFactors = FALSE)
  sz <- lengths(cr$members); nm <- names(cr$members); rp <- nm %in% cx$complex[cx$group == "RP-dominated"]
  if (opt$quick == 0L) stopifnot(sum(rp) == 8L, all(sz == cx$size[match(nm, cx$complex)]))
  sub_of <- function(keep) list(members = cr$members[keep], class = cr$class[keep])
  mk <- function(reg, macro, sel = NULL) list(reg = reg, pairs = comembership_pairs(reg, none, genes), co = NULL, operon = none, macro = macro, sel = sel)
  bench <- list(DoRothEA_AB = mk(tr$DoRothEA_AB, "DoRothEA_AB"), CORUM_all = mk(cr, "CORUM_all"),
                KD_response = mk(tr$KD_response, "KD_response"),
                CORUM_noRP = mk(sub_of(!rp), "CORUM_all", nm[!rp]), CORUM_RP = mk(sub_of(rp), "CORUM_all", nm[rp]),
                CORUM_small = mk(sub_of(sz >= 3 & sz <= 6), "CORUM_all", nm[sz >= 3 & sz <= 6]),
                CORUM_large = mk(sub_of(sz >= 20), "CORUM_all", nm[sz >= 20]))
  for (src in c("DoRothEA_AB", "CORUM_all", "KD_response")) { r <- bench[[src]]$reg
    macro_src[[src]] <- list(reg = r, co = comember_matrix(r, G), operon = none) }
  base_truths <- c("DoRothEA_AB", "CORUM_all", "KD_response"); extra_truths <- c("CORUM_noRP", "CORUM_RP", "CORUM_small", "CORUM_large")
  mem <- cr$members; modsize <- vapply(seq_len(G), function(g) { s <- sz[vapply(mem, function(m) g %in% m, NA)]; if (length(s)) max(s) else NA_real_ }, 0)
  modsize_label <- "largest CORUM complex containing the gene"
}
for (e in names(bench)) say(sprintf("truth %-12s %5d regulons, %8d pairs, base rate %.4f", e, length(bench[[e]]$reg$members),
  length(bench[[e]]$pairs$label), mean(bench[[e]]$pairs$label)))
base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)

## ---- jackknife streams (seed 20261016, order of analysis/baselines.R) -------------------------
set.seed(opt$seed)
gin_of <- function(reg) sort(unique(unlist(reg$members)))
draw_half <- function(gin) lapply(seq_len(opt$jack), function(r) sort(sample(gin, floor(length(gin) / 2))))
halves <- list()
if (opt$dataset == "m3d") {
  for (e in base_truths) halves[[e]] <- draw_half(bench[[e]]$pairs$genes_in)
} else {          # baselines.R k562ess order: DoRothEA A-B, A-C, CollecTRI, CORUM, KD response (A-C and CollecTRI only advance the stream)
  halves$DoRothEA_AB <- draw_half(bench$DoRothEA_AB$pairs$genes_in)
  invisible(draw_half(gin_of(dorothea_sets(genes, c("A", "B", "C"))))); invisible(draw_half(gin_of(collectri_sets(genes))))
  halves$CORUM_all <- draw_half(bench$CORUM_all$pairs$genes_in); halves$KD_response <- draw_half(bench$KD_response$pairs$genes_in)
}
gene_keep <- list()          # delete-20% gene subsets, first jack_reg of baselines.R's d20 stream
for (r in seq_len(opt$jack_reg)) { k <- logical(G); k[sample(G, G - round(0.2 * G))] <- TRUE; gene_keep[[paste0("d20_", r)]] <- k }
for (e in extra_truths) halves[[e]] <- draw_half(bench[[e]]$pairs$genes_in)
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))

## ---- label-free held-out depth (analysis/lab_splits.R criterion), BEFORE any OpenMP use in the parent ----
lab_split <- function() {
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab]) }
depth_f <- file.path(opt$out, "heldout_curves.csv")
DEPTH_VARS <- setdiff(VARS, "V0x")
if (!file.exists(depth_f) && opt$nsplit > 0L) {
  set.seed(opt$depth_seed); splits <- lapply(seq_len(opt$nsplit), function(s) lab_split())
  depth_split <- function(si) {
    h <- splits[[si]]; set.seed(opt$depth_seed + si); t0 <- Sys.time()
    Ms <- lapply(1:2, function(d) { Xm <- X[, h == d]; rng <- apply(Xm, 1, function(x) diff(range(x)))
      if (any(rng <= 0)) Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9)
      bspline_mi(Xm, bins = "hg", spline_order = 3, threads = 1, transform = "none") })
    rows <- list()
    for (v in DEPTH_VARS) {
      Cs <- lapply(1:2, function(d) rz_scores(Ms[[d]], v))
      for (dir in 1:2) { te <- 3L - dir
        sv <- Cs[[te]][UT]; target <- as.integer(sv >= stats::quantile(sv, 0.99)); rm(sv)
        P <- rz_lazy_P(rz_operator(Cs[[dir]])); crit <- numeric(opt$tmax); Pt <- NULL
        for (t in seq_len(opt$tmax)) { Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
          S <- (Pt + t(Pt)) / 2; crit[t] <- ap_fast(S[UT], target); rm(S) }
        rows[[length(rows) + 1L]] <- data.frame(split = si, direction = dir, variant = v, t = seq_len(opt$tmax), criterion = crit)
        rm(P, Pt, target); invisible(gc(verbose = FALSE)) }
      rm(Cs) }
    attr_t <- as.numeric(difftime(Sys.time(), t0, units = "secs")); out <- do.call(rbind, rows); out$secs <- attr_t; out }
  say(sprintf("depth stage: %d splits, t = 1..%d, variants %s, %d workers", opt$nsplit, opt$tmax, paste(DEPTH_VARS, collapse = ","), min(nc, opt$nsplit)))
  res <- parallel::mclapply(seq_len(opt$nsplit), depth_split, mc.cores = min(nc, opt$nsplit), mc.preschedule = FALSE)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), NA)
  if (any(bad)) stop("depth stage failed for splits: ", paste(which(bad), collapse = ","), " ", paste(unlist(res[bad]), collapse = " | "))
  cur <- do.call(rbind, res); utils::write.csv(cur[, setdiff(names(cur), "secs")], depth_f, row.names = FALSE)
  say(sprintf("depth stage done; per-split seconds: %s", paste(sprintf("%.0f", vapply(res, function(x) x$secs[1], 0)), collapse = " ")))
  rm(res, cur); invisible(gc(verbose = FALSE))
}
tstar <- stats::setNames(rep(opt$t, length(VARS)), VARS)
if (file.exists(depth_f)) {
  cur <- utils::read.csv(depth_f, stringsAsFactors = FALSE)
  ts <- do.call(rbind, lapply(split(cur, list(cur$split, cur$direction, cur$variant), drop = TRUE), function(x) {
    k <- which.max(x$criterion); near <- x$t[x$criterion >= 0.99 * x$criterion[k]]
    data.frame(split = x$split[1], direction = x$direction[1], variant = x$variant[1], t_best = x$t[k], crit_best = x$criterion[k],
               t_within1pct_lo = min(near), t_within1pct_hi = max(near)) }))
  utils::write.csv(ts, file.path(opt$out, "heldout_tstar.csv"), row.names = FALSE)
  tsplit <- do.call(rbind, lapply(split(ts, list(ts$split, ts$variant), drop = TRUE), function(x)
    data.frame(split = x$split[1], variant = x$variant[1], t_star = as.integer(round(exp(mean(log(x$t_best))))))))
  tsum <- do.call(rbind, lapply(split(tsplit, tsplit$variant), function(x) { y <- ts[ts$variant == x$variant[1], ]
    data.frame(variant = x$variant[1], n_splits = nrow(x), t_median = as.integer(floor(stats::median(x$t_star) + 0.5)),
               t_q25 = stats::quantile(x$t_star, 0.25, names = FALSE), t_q75 = stats::quantile(x$t_star, 0.75, names = FALSE),
               t_min = min(x$t_star), t_max = max(x$t_star), tmax_scanned = opt$tmax,
               within1pct_lo_median = stats::median(y$t_within1pct_lo), within1pct_hi_median = stats::median(y$t_within1pct_hi)) }))
  tsum <- tsum[match(DEPTH_VARS, tsum$variant), ]
  utils::write.csv(tsum, file.path(opt$out, "tstar_summary.csv"), row.names = FALSE)
  tstar[tsum$variant] <- tsum$t_median
  say("held-out t* (median over splits): ", paste(sprintf("%s=%d", tsum$variant, tsum$t_median), collapse = " "))
}

## ---- MI and the equivalence check (V0 path = package) -------------------------------------
t0 <- Sys.time()
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
f$calibrate(method = "normal", combine = "stouffer")
M <- f$mi; say(sprintf("MI + CLR: %.1f s; min MI %.3g", as.numeric(difftime(Sys.time(), t0, units = "secs")), min(M)))
C_pkg <- f$clr_scores
C0 <- rz_scores(M, "V0")
f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
A_pkg <- f$operator
P_pkg <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A_pkg, sparse = TRUE)
Pt <- as.matrix(P_pkg); for (s in seq_len(opt$t - 1L)) Pt <- as.matrix(Pt %*% P_pkg)
S_pkg <- (Pt + t(Pt)) / 2; rm(Pt)                 # analysis/baselines.R attention
A0 <- rz_operator(C0); P0 <- rz_lazy_P(A0); S0 <- rz_walk(P0, opt$t)[[as.character(opt$t)]]
eq <- data.frame(check = c("max |C(V0) - package clr_scores|", "max |A(V0) - package operator|", sprintf("max |S_t%d(V0) - package walk mass|", opt$t),
                           "max |C(V0) - C(V0x)| (diagonal convention)", "cor(C(V0), C(V0x)) over upper triangle"),
  value = c(max(abs(C0 - C_pkg)), max(abs(A0 - A_pkg)), max(abs(S0 - S_pkg)), NA, NA), stringsAsFactors = FALSE)
if ("V0x" %in% VARS) { Cx <- rz_scores(M, "V0x"); eq$value[4] <- max(abs(C0 - Cx)); eq$value[5] <- stats::cor(C0[UT], Cx[UT]); rm(Cx) }
utils::write.csv(eq, file.path(opt$out, "equivalence.csv"), row.names = FALSE); print(eq, row.names = FALSE)
stopifnot(all(eq$value[1:3] < 1e-8))
say("EQUIVALENCE PASSED (V0 path reproduces the package CLR, operator and walk mass to < 1e-8)")
rm(C_pkg, A_pkg, P_pkg, S_pkg, S0, A0); f$release(); rm(f); invisible(gc(verbose = FALSE))

## ---- common: stationary distribution of the V0 attention walk, neutral hubness, row-scale ------------
common_f <- file.path(opt$out, "cache", "common.rds")
q_att <- rz_stationary(P0); H_att <- outer(q_att, q_att, "+")
mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+")
rm(C0, P0)
set.seed(7L); samp <- M[sample.int(G * G, min(5e6, G * G))]; thr_glob <- unname(stats::quantile(samp, 0.99)); rm(samp)
cs <- list(V0 = rz_center_scale(M, "V0"), V1 = rz_center_scale(M, "V1"), V2 = rz_center_scale(M, "V2"), V4 = rz_center_scale(M, "V4"))
Moff <- rz_offdiag(M)
rowscale <- data.frame(gene = genes, mean_diag0 = cs$V0$center, sd_diag0 = cs$V0$scale,
  median = cs$V1$center, iqr_scale = cs$V1$scale, mad_scale = cs$V4$scale / 1, median_cube = cs$V2$center, iqr_scale_cube = cs$V2$scale,
  n_strong_global = colSums(Moff >= thr_glob), n_z4_gt3 = colSums((Moff - rep(cs$V4$center, each = G - 1L)) / rep(cs$V4$scale, each = G - 1L) > 3),
  module_size = modsize, stringsAsFactors = FALSE)
rowscale$sd_over_iqr <- rowscale$sd_diag0 / rowscale$iqr_scale; rowscale$sd_over_mad <- rowscale$sd_diag0 / rowscale$mad_scale
rowscale$mean_minus_median_in_iqr <- (rowscale$mean_diag0 - rowscale$median) / rowscale$iqr_scale
rm(Moff, cs)
utils::write.csv(rowscale, file.path(opt$out, "row_scale_genes.csv"), row.names = FALSE)
sp <- function(a, b) { ok <- is.finite(a) & is.finite(b); c(rho = stats::cor(a[ok], b[ok], method = "spearman"), n = sum(ok)) }
rs <- do.call(rbind, lapply(c("mean_diag0", "sd_diag0", "median", "iqr_scale", "mad_scale", "sd_over_iqr", "sd_over_mad", "mean_minus_median_in_iqr"), function(q)
  do.call(rbind, lapply(c("n_strong_global", "n_z4_gt3", "module_size"), function(n) { r <- sp(rowscale[[q]], rowscale[[n]])
    data.frame(row_quantity = q, versus = n, spearman = r[["rho"]], n_genes = r[["n"]]) }))))
rs$module_label <- ifelse(rs$versus == "module_size", modsize_label, NA)
utils::write.csv(rs, file.path(opt$out, "row_scale_spearman.csv"), row.names = FALSE)
dec <- cut(rank(rowscale$n_strong_global, ties.method = "first"), 10, labels = FALSE)
rd10 <- do.call(rbind, lapply(1:10, function(d) { s <- rowscale[dec == d, ]
  data.frame(decile = d, n_genes = nrow(s), n_strong_global_median = stats::median(s$n_strong_global), sd_median = stats::median(s$sd_diag0),
             iqr_scale_median = stats::median(s$iqr_scale), mad_scale_median = stats::median(s$mad_scale),
             sd_over_iqr_median = stats::median(s$sd_over_iqr), sd_over_mad_median = stats::median(s$sd_over_mad)) }))
utils::write.csv(rd10, file.path(opt$out, "row_scale_deciles.csv"), row.names = FALSE)
say(sprintf("row-scale diagnostics written (global strong-partner threshold = MI %.4g, top 1%% of cells)", thr_glob))

## ---- scoring -------------------------------------------------------------------------------
macro_names <- names(macro_src)
jack_one <- function(sc, b, keeps) {
  i <- b$pairs$i; j <- b$pairs$j; lab_ <- b$pairs$label
  unlist(parallel::mclapply(keeps, function(keep) { inh <- logical(G); inh[keep] <- TRUE; ix <- which(inh[i] & inh[j]); ap_fast(sc[ix], lab_[ix]) }, mc.cores = nc)) }
pooled <- function(S) lapply(stats::setNames(names(bench), names(bench)), function(e) { b <- bench[[e]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]
  list(aupr = ap_fast(sc, b$pairs$label), jack = jack_one(sc, b, halves[[e]])) })
evaluate <- function(S, macro = TRUE, hub = TRUE) {
  S[!is.finite(S)] <- 0; diag(S) <- 0
  tt <- Sys.time(); lap <- function(w) { say(sprintf("    evaluate: %s %.0f s", w, as.numeric(difftime(Sys.time(), tt, units = "secs")))); tt <<- Sys.time() }
  out <- list(pooled = pooled(S)); lap("pooled + jackknife")
  if (hub) { out$hub_att <- pooled(rz_hub_stratify(S, H_att)); invisible(gc(verbose = FALSE))
             out$hub_mi <- pooled(rz_hub_stratify(S, H_mi)); invisible(gc(verbose = FALSE)); lap("two hub-stratified scorings") }
  if (macro) out$macro <- lapply(macro_src, function(s) list(
    table = regulon_ap_table(S, s$reg, s$co, s$operon, genes),
    jack = do.call(rbind, parallel::mclapply(names(gene_keep), function(r) { x <- regulon_ap_table(S, s$reg, s$co, s$operon, genes, min_pos = 3L, keep = gene_keep[[r]])
      data.frame(rep = r, regulon = x$regulon, n_pos = x$n_pos, ap = x$ap) }, mc.cores = nc))))
  if (macro) lap("macro tables + gene jackknife")
  # hub drift: the top `top_pairs` pairs of S against two gene-connectedness rankings
  o <- order(S[UT], decreasing = TRUE, method = "radix")[seq_len(opt$top_pairs)]; k <- UT[o] - 1L; gi <- k %% G + 1L; gj <- k %/% G + 1L
  sUT <- S[UT]; sc_top <- sUT[o]
  out$drift <- do.call(rbind, lapply(list(att_pi = q_att, mi_rowmean = mi_rm), function(hh) { r <- rank(hh) / G; top5 <- r > 0.95
    data.frame(share_touching_top5pct = mean(top5[gi] | top5[gj]), mean_endpoint_rank = mean((r[gi] + r[gj]) / 2),
               chance_share = 1 - mean(!top5)^2, min_score_in_top = min(sc_top), frac_positive_score = mean(sUT > 0)) }))
  out$drift$stratifier <- c("att_pi", "mi_rowmean")
  out
}
guard <- list()
say(sprintf("variants: %s; t* used for the depth read-out: %s", paste(VARS, collapse = ","), paste(sprintf("%s=%d", names(tstar), tstar), collapse = " ")))
for (v in VARS) {
  light <- v == "V0x"
  t0 <- Sys.time(); C <- rz_scores(M, v); guard[[v]] <- attr(C, "n_zero_scale"); attr(C, "n_zero_scale") <- NULL
  cf <- file.path(opt$out, "cache", sprintf("%s__pair.rds", v))
  if (!rds_ok(cf)) { saveRDS(evaluate(C, macro = !light, hub = !light), cf); say(sprintf("%s pairwise scored (%.0f s)", v, as.numeric(difftime(Sys.time(), t0, units = "secs")))) }
  ts_need <- unique(c(opt$t, if (!light) tstar[[v]]))
  f9 <- file.path(opt$out, "cache", sprintf("%s__att%d.rds", v, opt$t)); fT <- file.path(opt$out, "cache", sprintf("%s__attT.rds", v))
  if (!rds_ok(f9) || (!light && !rds_ok(fT))) {
    A <- rz_operator(C); P <- rz_lazy_P(A); rm(A); W <- rz_walk(P, ts_need); rm(P)
    if (!rds_ok(f9)) { saveRDS(evaluate(W[[as.character(opt$t)]], macro = !light, hub = !light), f9); say(sprintf("%s attention t=%d scored", v, opt$t)) }
    if (!light && !rds_ok(fT)) { r <- if (tstar[[v]] == opt$t) readRDS(f9) else evaluate(W[[as.character(tstar[[v]])]], macro = FALSE, hub = TRUE)
      r$t <- tstar[[v]]; saveRDS(r, fT); say(sprintf("%s attention t*=%d scored", v, tstar[[v]])) }
    rm(W) }
  rm(C); invisible(gc(verbose = FALSE))
  say(sprintf("%s done (%.0f s); rows hit by the zero-scale guard: %d", v, as.numeric(difftime(Sys.time(), t0, units = "secs")), guard[[v]]))
}
meta <- list(dataset = opt$dataset, G = G, N = ncol(X), base_rate = base_rate, n_pairs = vapply(bench, function(b) length(b$pairs$label), 0L),
             n_regulons = vapply(bench, function(b) length(b$reg$members), 0L), guard = unlist(guard), tstar = tstar, thr_glob = thr_glob, opt = opt,
             macro_sel = lapply(bench, function(b) b$sel), macro_of = vapply(bench, function(b) b$macro %||% NA_character_, ""),
             regulon_size = lapply(macro_src, function(s) stats::setNames(lengths(s$reg$members), names(s$reg$members))))
saveRDS(meta, file.path(opt$out, "cache", "meta.rds"))
say("done")
