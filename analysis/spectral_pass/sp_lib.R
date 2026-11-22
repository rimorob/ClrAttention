# Spectral-embedding exploratory pass (docs/plans/spectral_embedding_pass_plan.md): shared pieces. Source from the repo root.
# Reuses, and does not redefine: the beeline.R set-up (ns_clr, heldout_depth, ...), score_truth (directed TF-node scorer), ap_fast, the complex truth / scorer / planted complexes of
# analysis/scgpt_ext/ext_lib.R, the operator gl_make_P, the symmetriser ad_Psym, the eigen-solver ad_eigs, the permuted-data edge rule ext_edge (ad_permute_genes, ed_t2), and the
# second-CLR read-out sk_c2z of analysis/explore_stack/stack_lib.R. New code here: the unit loader (a copy of the expression construction of score_ext.R, gated against its stored edge
# files), the spectral embedding, the read-outs (cosine, mclust, HDDC, pgmm, ring detector) and the variant registry. No truth is read in this file.
if (!exists("SP_OUT")) SP_OUT <- "results/spectral_pass"
if (!exists("SP_THREADS")) SP_THREADS <- 4L
suppressPackageStartupMessages({ source("analysis/tcaclr2/edge_depth.R"); source("analysis/tcaclr2/tsb_lib.R") })
args <- c("--out", SP_OUT, "--threads", as.character(SP_THREADS))
src0 <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src0)[1]
for (e in parse(text = src0[seq_len(cut_at - 1L)])) { d1 <- deparse(e)[1]
  if (identical(d1, "args <- commandArgs(trailingOnly = TRUE)")) next
  if (startsWith(d1, "script_dir <-")) { script_dir <- "analysis"; next }
  eval(e, globalenv()) }
bsrc <- readLines("analysis/beeline_baselines.R"); a <- grep("^score_truth <- function", bsrc)[1]; z <- grep("^say\\(\"beeline_baselines: stage", bsrc)[1] - 1L
eval(parse(text = bsrc[a:z]), globalenv())
for (e in parse("analysis/beeline_baselines.R")) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] == "case_X") eval(e, globalenv())
source("analysis/scoreboard/sb_lib.R"); source("analysis/explore_stack/stack_lib.R"); source("analysis/robust_z/robust_z_lib.R"); source("analysis/human_truths.R")
source("analysis/scgpt_ext/ext_lib.R")
suppressPackageStartupMessages({ library(mclust) })

## ---- fixed settings (docs/plans/spectral_embedding_pass_plan.md section 3) --------------------------------------------------------------
SP <- list(kcap = 50L, kgrid = c(5L, 10L, 20L, 50L), neig = 60L,        # edge-rule count is taken among the first 60 non-trivial eigenvalues, capped at 50
           mclust_G = 2:30, hddc_K = 2:20, pgmm_G = 2:4, pgmm_q = 1:3, pgmm_kmax = 30L,
           ring_Hmax = 3L, ring_Gmax = 6L, rmix_M = 1:2, rmix_ransac = 6000L, tda = TRUE, ring_cv = 1 / 3, ring_min_n = 8, ring_dbic = 10, ring_restarts = 3L,
           seed = 20261121L)
SP_DEV <- c("hESC", "hHep", "K562", "HepG2", "RPE1")                    # development set; Jurkat is sealed and never listed here

## ---- units: the expression matrix of score_ext.R for one dataset ----------------------------------------------------------------------
# Copy of the construction in score_ext.R (sc: log1p(count / total UMI * median control UMI) on the scGPT gene set, subsample rows; beeline: case_X on the scGPT genes). The scGPT gene
# vector is read from the export (gene list only). Checked by sp_gate_edge (the permuted-data edge must equal results/scgpt_ext/<dataset>/tcaclr2_edge_t2.csv to 1e-9).
cat_of <- function(h5, grp, field) { codes <- rhdf5::h5read(h5, paste0(grp, "/", field)); cats <- rhdf5::h5read(h5, paste0(grp, "/__categories/", field)); as.character(cats)[codes + 1L] }
sp_dataset <- function(dn, expr = "auto") {
  stopifnot(dn %in% SP_DEV)                                                   # Jurkat cannot be opened through this function
  R <- EXT_REG[[dn]]; exports <- path.expand(R$exports); L <- list(dn = dn, R = R)
  if (R$type == "sc") {
    prep <- R$dir
    gs <- utils::read.delim(file.path(prep, "geneset.tsv"), stringsAsFactors = FALSE); sub <- utils::read.delim(file.path(prep, "subsamples.tsv"), stringsAsFactors = FALSE)
    metatxt <- readLines(file.path(prep, "geneset_meta.json"), warn = FALSE)
    MED <- as.numeric(sub("[^0-9.]*([0-9.]+).*", "\\1", grep("median_control_umi", metatxt, value = TRUE)))
    NGS <- as.integer(sub("[^0-9]*([0-9]+).*", "\\1", grep("gene_set_size", metatxt, value = TRUE))); stopifnot(is.finite(MED), nrow(gs) == NGS)
    gs$is_tf_part <- gs$is_tf_part %in% c("True", "TRUE", TRUE)
    EX <- if (expr == "auto") (if (file.exists(R$h5)) "h5" else "files") else expr; stopifnot(EX %in% c("h5", "files"))
    if (EX == "h5") { vg <- cat_of(R$h5, "var", "gene_name"); gi <- match(gs$symbol, vg); stopifnot(!anyNA(gi), !anyDuplicated(gi)); bar <- as.character(rhdf5::h5read(R$h5, "obs/cell_barcode")) }
    else tots <- utils::read.delim(file.path(prep, "cell_totals.tsv"), stringsAsFactors = FALSE)
    L$units <- sort(unique(sub$subsample)); L$unit_case <- function(k) paste0(R$prefix, k); L$expr <- EX
    L$get_X <- function(k, genes) { rows <- sub$file_row[sub$subsample == k]
      if (EX == "h5") {
        C <- rhdf5::h5read(R$h5, "X", index = list(NULL, rows + 1L)); storage.mode(C) <- "double"
        stopifnot(identical(bar[rows + 1L], sub$cell_barcode[sub$subsample == k]))
        Lg <- log1p(t(t(C) / colSums(C)) * MED)[gi, , drop = FALSE]; rownames(Lg) <- gs$symbol
      } else {
        tk <- tots[tots$subsample == k, ]; stopifnot(identical(tk$cell_barcode, sub$cell_barcode[sub$subsample == k]), identical(tk$file_row, rows))
        C <- matrix(readBin(file.path(prep, "expr", sprintf("%s%d.counts.f32", R$prefix, k)), "numeric", n = nrow(gs) * length(rows), size = 4L), nrow(gs), length(rows))
        Lg <- log1p(t(t(C) / tk$total_umi) * MED); rownames(Lg) <- gs$symbol }
      Lg[genes, , drop = FALSE] }
    L$tf_symbols <- gs$symbol[gs$is_tf_part]
  } else {
    L$units <- R$cases; L$unit_case <- function(k) k; L$expr <- "beeline"
    L$get_X <- function(k, genes) { X <- case_X(R$ds, as.integer(sub("^.*_", "", k))); stopifnot(all(genes %in% rownames(X))); X[genes, , drop = FALSE] }
    L$tf_symbols <- NULL
  }
  L$genes_of <- function(k) ext_rd_genes(exports, "scgpt_A", L$unit_case(k))      # gene list only; no score file is opened
  L$unit_X <- function(k) { g <- L$genes_of(k); list(X = L$get_X(k, g), genes = g) }
  L
}
sp_gate_edge <- function(dn, k, e, tol = 1e-9) {      # edge and t2 of the unit vs the stored tcaCLR2 edge file (non-truth output of score_ext.R)
  f <- file.path("results/scgpt_ext", dn, "tcaclr2_edge_t2.csv")
  if (!file.exists(f)) return(list(diff = NA_real_, t2_equal = NA, pass = NA, stored_edge = NA_real_, stored_t2 = NA_integer_))
  s <- utils::read.csv(f, stringsAsFactors = FALSE); s <- s[as.character(s$unit) == as.character(k), ]
  stopifnot(nrow(s) == 1L); d <- abs(s$edge - e$edge)
  list(diff = d, t2_equal = s$t2 == e$t2, pass = d <= tol && s$t2 == e$t2, stored_edge = s$edge, stored_t2 = s$t2)
}

## ---- spectral embedding -----------------------------------------------------------------------------------------------------------------
# Operator: P = gl_make_P(ns_clr(X), keff 10, cap 50, lazy 0.5) (the aCLR operator); symmetrised P_sym = ad_Psym(P) (symmetric, spectrum in [0, 1]); eigenpairs by ad_eigs. The trivial
# eigenpair (value 1) is dropped. Coordinates: Y[, j] = lambda_j^(t2 / 2) * v_j (non-trivial j = 1, 2, ...), so Y Y^T is the rank-k part of P_sym^t2.
# k by the edge rule = number of non-trivial eigenvalues above the permuted-data edge (the same edge as tcaCLR2, ext_edge), among the first SP$neig, capped at SP$kcap.
sp_embed <- function(Cn, edge, t2, neig = SP$neig) {
  P <- gl_make_P(Cn, ED$keff, ED$cap); Ps <- ad_Psym(P); rm(P)
  ev <- ad_eigs(Ps, neig + 1L, vectors = TRUE); stopifnot(abs(ev$values[1L] - 1) < 1e-6)
  lam <- ev$values[-1L]; V <- ev$vectors[, -1L, drop = FALSE]
  n_above <- sum(lam > edge); kE <- min(SP$kcap, n_above)
  list(lam = lam, V = V, edge = edge, t2 = t2, n_above = n_above, k_edge = kE,
       Y = function(k) { k <- min(k, length(lam)); sweep(V[, seq_len(k), drop = FALSE], 2, pmax(lam[seq_len(k)], 0)^(t2 / 2), `*`) })
}
sp_cos <- function(Y) { n <- sqrt(rowSums(Y^2)); n[n == 0] <- 1; Yn <- Y / n; S <- tcrossprod(Yn); diag(S) <- 0; S }

## ---- read-out: mixtures ---------------------------------------------------------------------------------------------------------------
sp_post_score <- function(Z) { Z[!is.finite(Z)] <- 0; S <- tcrossprod(Z); diag(S) <- 0; S }       # pair score = sum_k P(k|i) P(k|j)
# mclust: BIC over G = 1..SP$mclust_G max and all 14 covariance models; the SCORED model is the best BIC with G >= 2 (as specified); `g1_wins` records whether G = 1 had the better BIC.
sp_mclust <- function(Y, G = SP$mclust_G, seed = SP$seed) {
  set.seed(seed); t0 <- Sys.time()
  b <- mclustBIC(Y, G = c(1L, G), verbose = FALSE)                                     # modelNames default = all 14 for d > 1
  ok <- b[, , drop = FALSE]; ok1 <- ok[rownames(ok) == "1", , drop = FALSE]; ok2 <- ok[rownames(ok) != "1", , drop = FALSE]
  best2 <- suppressWarnings(max(ok2, na.rm = TRUE)); best1 <- suppressWarnings(max(ok1, na.rm = TRUE))
  if (!is.finite(best2)) return(list(Z = NULL, G = NA_integer_, model = NA_character_, bic = NA_real_, g1_wins = NA, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), n_feasible = 0L))
  w <- which(ok2 == best2, arr.ind = TRUE)[1, ]; Gb <- as.integer(rownames(ok2)[w[1]]); mn <- colnames(ok2)[w[2]]
  m <- mclust::Mclust(Y, G = Gb, modelNames = mn, verbose = FALSE)
  list(Z = m$z, G = Gb, model = mn, bic = best2, g1_wins = is.finite(best1) && best1 > best2, bic_g1 = best1, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), n_feasible = sum(is.finite(ok)))
}
# HDDC (HDclassif): all 14 models, K = 2..SP$hddc_K, BIC for both the number of components and the intrinsic dimension of each component (d_select = "BIC"), k-means start.
sp_hddc <- function(Y, K = SP$hddc_K, seed = SP$seed) {
  set.seed(seed); t0 <- Sys.time()
  m <- tryCatch(suppressWarnings(HDclassif::hddc(Y, K = K, model = "ALL", criterion = "bic", d_select = "BIC", init = "kmeans", itermax = 60, show = FALSE, nb.rep = 1)), error = function(e) e)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(m, "error") || is.null(m$posterior)) return(list(Z = NULL, K = NA_integer_, model = NA_character_, d = NA_character_, bic = NA_real_, seconds = sec))
  list(Z = m$posterior, K = m$K, model = m$model, d = paste(m$d, collapse = ";"), bic = max(m$BIC, na.rm = TRUE), seconds = sec)
}
# pgmm (mixtures of parsimonious Gaussian mixture / factor analysers): all 8 models, G = 2..SP$pgmm_G, q = 1..SP$pgmm_q (q common to the components), BIC.
sp_pgmm <- function(Y, G = SP$pgmm_G, q = SP$pgmm_q, seed = SP$seed) {
  t0 <- Sys.time(); qmax <- floor((2 * ncol(Y) + 1 - sqrt(8 * ncol(Y) + 1)) / 2 + 1e-9); q <- q[q <= qmax]; if (!length(q)) q <- 1L   # Ledermann bound on q for d columns
  m <- NULL; utils::capture.output(m <- tryCatch(suppressWarnings(suppressMessages(pgmm::pgmmEM(Y, rG = G, rq = q, zstart = 2, loop = 3, seed = seed, relax = FALSE))), error = function(e) e))
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(m, "error") || is.null(m$zhat)) return(list(Z = NULL, G = NA_integer_, q = NA_integer_, model = NA_character_, bic = NA_real_, seconds = sec))
  list(Z = m$zhat, G = m$g, q = m$q, model = m$model, bic = suppressWarnings(max(unlist(m$bic), na.rm = TRUE)), seconds = sec)
}

## ---- read-out: ring detector (polar transform) ------------------------------------------------------------------------------------------
# For one plane (x, y) of the weighted embedding, centre = origin (the coordinates are centred eigenvectors). Polar transform r = sqrt(x^2 + y^2), theta = atan2(y, x).
# RING MODEL (H = 0..SP$ring_Hmax): a mixture of one isotropic Gaussian "bulk" at the origin (density exp(-r^2 / 2 s0^2) / (2 pi s0^2)) and H rings; ring h has a Gaussian radial profile
# N(r; mu_h, s_h) and a uniform angle, i.e. the planar density N(r; mu_h, s_h) / (2 pi r), with the constraint s_h <= mu_h * SP$ring_cv (narrow radial spread) and an effective size of
# at least SP$ring_min_n genes. Fitted by EM (restarts: radius quantile start and SP$ring_restarts seeded random starts). Parameters: 1 (s0) + 3 per ring (mu, s, weight).
# CLUSTER MODEL: mclust BIC (14 covariance models, G = 1..SP$ring_Gmax) on (x, y). Both BICs are 2 logL - npar log n on the same planar likelihood.
# The plane is ACCEPTED as ring-bearing if the best ring model with H >= 1 beats the best cluster model by more than SP$ring_dbic BIC units.
sp_ring_em <- function(r, H, init_mu, min_n = SP$ring_min_n, cv = SP$ring_cv, iters = 200L) {
  n <- length(r); w <- rep(1 / (H + 1), H + 1); s0 <- sqrt(mean(r^2) / 2); mu <- init_mu; s <- pmax(0.1 * mu, 1e-6); ll_old <- -Inf
  for (it in seq_len(iters)) {
    d0 <- w[1] * exp(-r^2 / (2 * s0^2)) / (2 * pi * s0^2)
    dh <- vapply(seq_len(H), function(h) w[h + 1] * stats::dnorm(r, mu[h], s[h]) / (2 * pi * pmax(r, 1e-12)), numeric(n))
    dm <- cbind(d0, matrix(dh, n, H)); tot <- rowSums(dm); tot[tot <= 0] <- .Machine$double.xmin; ll <- sum(log(tot)); g <- dm / tot
    nk <- colSums(g); if (any(nk[-1L] < min_n) && H > 0) return(NULL)
    w <- nk / n; s0 <- sqrt(max(sum(g[, 1] * r^2) / (2 * nk[1]), 1e-12))
    for (h in seq_len(H)) { mu[h] <- sum(g[, h + 1] * r) / nk[h + 1]; s[h] <- sqrt(max(sum(g[, h + 1] * (r - mu[h])^2) / nk[h + 1], 1e-12)); s[h] <- min(max(s[h], 0.01 * mu[h]), cv * mu[h]) }
    if (abs(ll - ll_old) < 1e-8 * abs(ll)) break; ll_old <- ll }
  list(ll = ll, w = w, s0 = s0, mu = mu, s = s, g = g, npar = 1L + 3L * H, H = H)
}
sp_ring_fit <- function(r, H, seed) {
  if (H == 0L) { s0 <- sqrt(mean(r^2) / 2); ll <- sum(log(exp(-r^2 / (2 * s0^2)) / (2 * pi * s0^2))); return(list(ll = ll, w = 1, s0 = s0, mu = numeric(0), s = numeric(0), g = matrix(1, length(r), 1), npar = 1L, H = 0L)) }
  set.seed(seed); best <- NULL; rq <- stats::quantile(r, seq_len(H) / (H + 1)); starts <- list(unname(rq))
  for (b in seq_len(SP$ring_restarts)) starts[[length(starts) + 1L]] <- sort(sample(r[r > stats::median(r)], H))
  for (st in starts) { f <- tryCatch(sp_ring_em(r, H, st), error = function(e) NULL); if (!is.null(f) && is.finite(f$ll) && (is.null(best) || f$ll > best$ll)) best <- f }
  best
}
sp_plane_cluster_bic <- function(x, y) suppressWarnings(max(mclustBIC(cbind(x, y), G = seq_len(SP$ring_Gmax), verbose = FALSE), na.rm = TRUE))
sp_ring_plane <- function(x, y, seed = SP$seed, bc = sp_plane_cluster_bic(x, y)) {
  r <- sqrt(x^2 + y^2); th <- atan2(y, x); n <- length(r)
  fits <- lapply(0:SP$ring_Hmax, function(H) sp_ring_fit(r, H, seed + H)); bics <- vapply(fits, function(f) if (is.null(f)) NA_real_ else 2 * f$ll - f$npar * log(n), 0)
  bR1 <- bics[-1L]; if (all(!is.finite(bR1))) { bR <- -Inf; Hb <- 0L } else { Hb <- which.max(replace(bR1, !is.finite(bR1), -Inf)); bR <- bR1[Hb] }    # best ring model with H >= 1
  acc <- is.finite(bR) && (bR - bc > SP$ring_dbic); f <- fits[[Hb + 1L]]
  list(accepted = acc, H = if (acc) Hb else 0L, bic_ring = bR, bic_cluster = bc, bic_ring_H0 = bics[1L], mu = if (acc) f$mu else numeric(0), s = if (acc) f$s else numeric(0),
       w = if (acc) f$w[-1L] else numeric(0), gamma = if (acc) f$g[, -1L, drop = FALSE] else NULL, theta = th, r = r)
}
# ---- ring family B: mixture of radial Gaussians with a fitted centre ("rmix") ---------------------------------------------------------------
# Planar density of component m with centre c, r = |x - c|:  f_m(x) = exp(-(r - mu_m)^2 / (2 s_m^2)) / Z(mu_m, s_m),  Z = 2 pi [ s^2 exp(-mu^2 / 2 s^2) + mu s sqrt(2 pi) Phi(mu / s) ],
# mu >= 0, s > 0. This is a PROPER density on the plane for every (mu, s) (mu = 0 is an isotropic Gaussian blob; mu > 0 with small s a thin ring; mu > 0 with large s a fat doughnut). It equals the
# PI's N(r; mu, s) / (2 pi r) to first order for thin rings; the literal form is not integrable at r = 0 unless the radial Gaussian vanishes there, so the normalised form above is used.
# Model: sum_m w_m f_m, M = 1 or 2 components, common centre. Parameters counted: 2 (centre) + 2 M (mu, s) + (M - 1) (weights) = 3 M + 1. Fitted by L-BFGS-B on the exact log-likelihood from
# two starts (the centroid; the best circle of SP$rmix_ransac random three-point circles, scored by the number of points within 0.1 RMS radius of the circle). BIC = 2 logL - npar log n, directly comparable with the mclust BIC on the same points.
# A component is a RING component if mu >= s and its effective size is at least SP$ring_min_n; the plane is ACCEPTED if the best model with at least one ring component beats the plane's mclust BIC
# (G = 1..SP$ring_Gmax, 14 models) by more than SP$ring_dbic.
sp_rmix_ll <- function(par, x, y, M, parts = FALSE) {
  cx <- par[1]; cy <- par[2]; r <- sqrt((x - cx)^2 + (y - cy)^2); mu <- par[3 + 2 * (seq_len(M) - 1)]; s <- par[4 + 2 * (seq_len(M) - 1)]
  w <- if (M == 1L) 1 else { a <- c(0, par[3 + 2 * M]); exp(a - max(a)) / sum(exp(a - max(a))) }
  lz <- log(2 * pi) + log(s^2 * exp(-mu^2 / (2 * s^2)) + mu * s * sqrt(2 * pi) * stats::pnorm(mu / s))
  L <- vapply(seq_len(M), function(m) log(w[m]) - (r - mu[m])^2 / (2 * s[m]^2) - lz[m], numeric(length(r)))
  L <- matrix(L, length(r), M); mx <- apply(L, 1, max); ll <- mx + log(rowSums(exp(L - mx)))
  if (parts) return(list(ll = sum(ll), g = exp(L - ll), w = w, mu = mu, s = s, centre = c(cx, cy), r = r))
  -sum(ll)
}
sp_rmix_fit <- function(x, y, M, seed) {
  set.seed(seed); n <- length(x); sc <- sqrt(mean((x - mean(x))^2 + (y - mean(y))^2)); out <- NULL
  ransac <- function() { best <- c(NA, NA, NA); bn <- -1; tol <- 0.1 * sc
    for (b in seq_len(SP$rmix_ransac)) { j <- sample.int(n, 3L); A <- cbind(2 * (x[j[2:3]] - x[j[1]]), 2 * (y[j[2:3]] - y[j[1]])); rhs <- x[j[2:3]]^2 - x[j[1]]^2 + y[j[2:3]]^2 - y[j[1]]^2
      if (abs(det(A)) < 1e-12) next; c0 <- solve(A, rhs); R <- sqrt((x[j[1]] - c0[1])^2 + (y[j[1]] - c0[2])^2); if (!is.finite(R) || R > 4 * sc || R < 0.2 * sc) next
      m <- sum(abs(sqrt((x - c0[1])^2 + (y - c0[2])^2) - R) < tol); if (m > bn) { bn <- m; best <- c(c0, R) } }
    best }
  cs <- list(c(mean(x), mean(y), NA)); rc <- ransac(); if (!anyNA(rc)) cs[[2L]] <- rc
  lo <- c(-3 * sc, -3 * sc, rep(c(0, 0.005 * sc), M), if (M > 1L) -6); hi <- c(3 * sc, 3 * sc, rep(c(4 * sc, 3 * sc), M), if (M > 1L) 6)
  for (c0 in cs) { r0 <- sqrt((x - c0[1])^2 + (y - c0[2])^2); sd0 <- max(stats::mad(r0), 0.02 * sc)
    mu0 <- if (!is.na(c0[3])) { if (M == 1L) c0[3] else c(c0[3], stats::median(r0)) } else { if (M == 1L) stats::median(r0) else stats::quantile(r0, c(0.25, 0.75)) }
    s0 <- if (!is.na(c0[3])) { if (M == 1L) 0.05 * sc else c(0.05 * sc, sd0) } else rep(sd0, M)
    p0 <- c(c0[1:2], as.vector(rbind(pmin(mu0, 4 * sc), s0)), if (M > 1L) 0)
    f <- tryCatch(stats::optim(pmin(pmax(p0, lo), hi), sp_rmix_ll, x = x, y = y, M = M, method = "L-BFGS-B", lower = lo, upper = hi, control = list(maxit = 150, factr = 1e9)), error = function(e) NULL)
    if (!is.null(f) && is.finite(f$value) && (is.null(out) || f$value < out$value)) out <- f }
  if (is.null(out)) return(NULL); d <- sp_rmix_ll(out$par, x, y, M, parts = TRUE); d$npar <- 3L * M + 1L; d$M <- M; d
}
sp_rmix_plane <- function(x, y, seed = SP$seed, bc = sp_plane_cluster_bic(x, y)) {
  n <- length(x); fits <- lapply(SP$rmix_M, function(M) sp_rmix_fit(x, y, M, seed + 10L * M))
  ok <- lapply(fits, function(f) { if (is.null(f)) return(NULL); ring <- f$mu >= f$s & colSums(f$g) >= SP$ring_min_n; if (!any(ring)) return(NULL); f$ring <- ring; f$bic <- 2 * f$ll - f$npar * log(n); f })
  ok <- Filter(Negate(is.null), ok)
  if (!length(ok)) return(list(accepted = FALSE, H = 0L, bic_ring = -Inf, bic_cluster = bc, mu = numeric(0), s = numeric(0), w = numeric(0), gamma = NULL, theta = rep(0, n), centre = c(NA, NA), M = NA_integer_))
  f <- ok[[which.max(vapply(ok, function(z) z$bic, 0))]]; acc <- f$bic - bc > SP$ring_dbic
  list(accepted = acc, H = if (acc) sum(f$ring) else 0L, bic_ring = f$bic, bic_cluster = bc, mu = f$mu[f$ring], s = f$s[f$ring], w = f$w[f$ring], gamma = f$g[, f$ring, drop = FALSE],
       theta = atan2(y - f$centre[2], x - f$centre[1]), centre = f$centre, M = f$M)
}
# persistent-homology VALIDATOR (not a read-out; never enters a pair score): largest H1 bar (death - birth) of the Vietoris-Rips filtration of the plane's points, in units of the plane's RMS radius
# about its mean (TDAstats::calculate_homology, dim 1). Compared with the same statistic on the planes of the per-gene permuted data of the same unit (sp_tda_null).
sp_tda_h1 <- function(x, y) { xy <- cbind(x - mean(x), y - mean(y)); sc <- sqrt(mean(rowSums(xy^2))); h <- TDAstats::calculate_homology(xy / sc, dim = 1); h <- h[h[, "dimension"] == 1, , drop = FALSE]
  if (!nrow(h)) 0 else max(h[, "death"] - h[, "birth"]) }
# all adjacent planes (j, j + 1), j = 1 .. jmax of the weighted embedding Y (computed once: the coordinates do not depend on k). Each plane: list(ring = family A, rmix = family B, h1 = TDA validator).
sp_ring_planes <- function(Y, jmax = ncol(Y) - 1L, seed = SP$seed, tda = SP$tda) lapply(seq_len(min(jmax, ncol(Y) - 1L)), function(j) { x <- Y[, j]; y <- Y[, j + 1L]; bc <- sp_plane_cluster_bic(x, y)
  list(ring = sp_ring_plane(x, y, seed + 100L * j, bc), rmix = sp_rmix_plane(x, y, seed + 100L * j, bc), h1 = if (tda) sp_tda_h1(x, y) else NA_real_) })
# per-gene permuted data of the same unit: the TDA statistic on its planes 1..jmax (the null of the validator). Uses the unit's own t2 and the permuted matrix's own eigenvectors.
sp_tda_null <- function(X, t2, edge, jmax, seed = 90001L) { Xp <- ad_permute_genes(X, seed); emb <- sp_embed(ns_clr(Xp), edge, t2); Y <- emb$Y(jmax + 1L)
  vapply(seq_len(jmax), function(j) sp_tda_h1(Y[, j], Y[, j + 1L]), 0) }
# pair score from the planes j <= k - 1 that were accepted in family `fam` ("ring" or "rmix"): sum over planes and ring components of gamma_ih gamma_jh (1 + cos(theta_i - theta_j)) / 2 (angle about the
# plane's centre: the origin for family A, the fitted centre for family B); zero when no plane is accepted
sp_ring_score <- function(planes, k, G, fam = "ring") {
  S <- matrix(0, G, G); nacc <- 0L; nr <- 0L
  for (j in seq_len(min(length(planes), k - 1L))) { p <- planes[[j]][[fam]]; if (!p$accepted) next
    nacc <- nacc + 1L; nr <- nr + p$H; C <- cos(p$theta); Sn <- sin(p$theta); A <- (1 + outer(C, C) + outer(Sn, Sn)) / 2
    for (h in seq_len(p$H)) S <- S + outer(p$gamma[, h], p$gamma[, h]) * A }
  diag(S) <- 0; list(S = S, planes_accepted = nacc, rings = nr)
}

## ---- variant registry -----------------------------------------------------------------------------------------------------------------
# 5 embeddings (k rule: edge, 5, 10, 20, 50) x 7 read-outs (cos, cos_clr2, mclust, hddc, pgmm, ring = family A, rmix = family B) = 35 variants, M = 35. name = <readout>_k<E|5|10|20|50>
SP_READOUTS <- c("cos", "cos_clr2", "mclust", "hddc", "pgmm", "ring", "rmix")
SP_KRULES <- c("E", "5", "10", "20", "50")
sp_variants <- function() { v <- expand.grid(readout = SP_READOUTS, krule = SP_KRULES, stringsAsFactors = FALSE); v$variant <- paste0(v$readout, "_k", v$krule); v }
sp_kval <- function(emb, krule) if (krule == "E") emb$k_edge else as.integer(krule)

# All 35 pair-score matrices of ONE unit (lazy: returns a list name -> list(S, info)). `only` restricts the read-outs (smoke tests); `say` logs.
sp_unit_scores <- function(emb, say = function(...) invisible(NULL), only = SP_READOUTS, krules = SP_KRULES, ring_planes = NULL, variants = NULL) {   # `variants`: optional character vector of variant names (the frozen list); only those are computed
  G <- nrow(emb$V); out <- list(); info <- list()
  on_of <- function(kr) if (is.null(variants)) only else only[paste0(only, "_k", kr) %in% variants]
  krules <- krules[vapply(krules, function(kr) length(on_of(kr)) > 0L, NA)]
  rk <- krules[vapply(krules, function(kr) any(c("ring", "rmix") %in% on_of(kr)), NA)]; need_ring <- length(rk) > 0L
  kmaxr <- if (need_ring) max(vapply(rk, function(r) sp_kval(emb, r), 0L)) else 0L
  if (need_ring && is.null(ring_planes) && kmaxr >= 2L) { t0 <- Sys.time(); ring_planes <- sp_ring_planes(emb$Y(kmaxr)); say(sprintf("  ring planes 1..%d: %.1f s, accepted ring %d, rmix %d", length(ring_planes), as.numeric(difftime(Sys.time(), t0, units = "secs")), sum(vapply(ring_planes, function(p) p$ring$accepted, NA)), sum(vapply(ring_planes, function(p) p$rmix$accepted, NA)))) }
  for (kr in krules) { k <- sp_kval(emb, kr); on <- on_of(kr)
    if (k < 2L) { for (ro in on) { nm <- paste0(ro, "_k", kr); out[[nm]] <- NULL; info[[nm]] <- data.frame(variant = nm, readout = ro, krule = kr, k = k, status = "not_scored: fewer than 2 eigenvalues above the edge") }; next }
    Y <- emb$Y(k)
    if ("cos" %in% on || "cos_clr2" %in% on) { S <- sp_cos(Y)
      if ("cos" %in% on) { out[[paste0("cos_k", kr)]] <- S; info[[paste0("cos_k", kr)]] <- data.frame(variant = paste0("cos_k", kr), readout = "cos", krule = kr, k = k, status = "ok") }
      if ("cos_clr2" %in% on) { out[[paste0("cos_clr2_k", kr)]] <- sk_c2z(S); info[[paste0("cos_clr2_k", kr)]] <- data.frame(variant = paste0("cos_clr2_k", kr), readout = "cos_clr2", krule = kr, k = k, status = "ok") } }
    if ("mclust" %in% on) { m <- sp_mclust(Y); nm <- paste0("mclust_k", kr)
      if (!is.null(m$Z)) out[[nm]] <- sp_post_score(m$Z)
      info[[nm]] <- data.frame(variant = nm, readout = "mclust", krule = kr, k = k, status = if (is.null(m$Z)) "no feasible model" else "ok", n_comp = m$G, model = m$model, bic = m$bic, g1_wins = m$g1_wins, at_max = isTRUE(m$G == max(SP$mclust_G)), seconds = m$seconds)
      say(sprintf("  mclust k=%d (rule %s): G=%s model=%s G1 better=%s (%.1f s)", k, kr, m$G, m$model, m$g1_wins, m$seconds)) }
    if ("hddc" %in% on) { m <- sp_hddc(Y); nm <- paste0("hddc_k", kr)
      if (!is.null(m$Z)) out[[nm]] <- sp_post_score(m$Z)
      info[[nm]] <- data.frame(variant = nm, readout = "hddc", krule = kr, k = k, status = if (is.null(m$Z)) "no feasible model" else "ok", n_comp = m$K, model = m$model, extra = m$d, bic = m$bic, at_max = isTRUE(m$K == max(SP$hddc_K)), seconds = m$seconds)
      say(sprintf("  hddc k=%d (rule %s): K=%s model=%s d=%s (%.1f s)", k, kr, m$K, m$model, m$d, m$seconds)) }
    if ("pgmm" %in% on && k > SP$pgmm_kmax) { nm <- paste0("pgmm_k", kr); info[[nm]] <- data.frame(variant = nm, readout = "pgmm", krule = kr, k = k, status = sprintf("not_run: pgmm limited to k <= %d (run time; k = 50 ran 43 min on HepG2 and returned no model)", SP$pgmm_kmax)); on <- setdiff(on, "pgmm") }
    if ("pgmm" %in% on) { m <- sp_pgmm(Y); nm <- paste0("pgmm_k", kr)
      if (!is.null(m$Z)) out[[nm]] <- sp_post_score(m$Z)
      info[[nm]] <- data.frame(variant = nm, readout = "pgmm", krule = kr, k = k, status = if (is.null(m$Z)) "no feasible model" else "ok", n_comp = m$G, model = m$model, extra = paste0("q=", m$q), bic = m$bic, at_max = isTRUE(m$G == max(SP$pgmm_G)), seconds = m$seconds)
      say(sprintf("  pgmm k=%d (rule %s): G=%s q=%s model=%s (%.1f s)", k, kr, m$G, m$q, m$model, m$seconds)) }
    for (fam in intersect(c("ring", "rmix"), on)) { rs <- sp_ring_score(ring_planes, k, G, fam); nm <- paste0(fam, "_k", kr); out[[nm]] <- rs$S
      info[[nm]] <- data.frame(variant = nm, readout = fam, krule = kr, k = k, status = "ok", n_comp = rs$rings, extra = paste0("planes_accepted=", rs$planes_accepted)) } }
  list(S = out, info = dplyr_bind(info), ring_planes = ring_planes)
}
dplyr_bind <- function(l) { cols <- unique(unlist(lapply(l, names))); do.call(rbind, lapply(l, function(d) { for (cn in setdiff(cols, names(d))) d[[cn]] <- NA; d[, cols, drop = FALSE] })) }
