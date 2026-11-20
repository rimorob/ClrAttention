#!/usr/bin/env Rscript
# Blocked held-out Gaussian-copula graph-signal likelihood (plan Revision 2, item 3 and 8b; review candidate 1). EXPLORATORY; no truth is read.
#   Rscript analysis/tcaclr_depth/adaptive/graph_signal_lik.R --matrix results/baselines_m3d/expression.tsv --dataset M3D --kind lab --splits 8 --keffs 10 --out results/tcaclr_depth/adaptive_m3d
#
# Model. x_n (one experiment, rank-normal transformed gene by gene, so a Gaussian copula; unit variance per gene) ~ N(0, Sigma),
#   Sigma = U diag(c * lambda_k^t) U' + sigma^2 I,
# U, lambda = the top-r eigenpairs of the symmetrised-by-degree lazy operator P_sym built (grid library operator, keff, cap = 4 keff) from the GRAPH half
# (adaptive_lib.R). The r - G remaining (trailing) modes get the variance sigma^2 + c * lambda_r^t (an upper bound for them). 3 hyperparameters: t (depth-like
# decay, one value per likelihood evaluation), c (amplitude), sigma^2 (noise floor). The floor sigma^2 is also the "mixture floor": a mode with lambda^t -> 0
# (unreachable) keeps the variance sigma^2 and never gives -infinity. The likelihood depends on the data only through the energies e_kn = (u_k' x_n)^2 and
# the tail energy ||x_n||^2 - sum_k e_kn.
# Blocked held-out scheme (no penalty): the graph (and the ranks of the transform) come from half h of the experiments, (c, sigma^2) are fitted by maximum
# likelihood at each t on THAT half's energies, and the log-likelihood is scored on the other half (and the two directions are pooled). Halves are unions of
# resampling units: laboratories (M3D, the stored lab splits) or k-means clusters of the knockdown profiles (K562). The score is per gene and per experiment,
# averaged within a held-out block (laboratory / cluster) first and then over blocks; the reference is the iid model (t = 0, one variance, no graph).
# Selection: the best t over the grid (largest block-mean score) and the ONE-STANDARD-ERROR rule toward the cheaper (smaller) t: the smallest t whose block-mean
# score is within SE(diff) of the best, SE(diff) = sd over blocks of (score(t) - score(best)) / sqrt(#blocks) (paired). A second column,
# t_selected_1se_marginal, uses the SE of the best's own block-mean score (sd over blocks / sqrt(#blocks)), which is much wider; both are reported.
# N_eff = (tr C)^2 / tr(C^2) (experiments) and G_eff (genes) of the matrix and of every held-out half are reported (gl_neff).
# Output (<out>/): gs_surface_<tag>.tsv (split, keff, t [0 = iid reference], ll_per_gene_exp, gain_vs_iid, se_gain, se_diff_best, fitted c / sigma^2 means, flags
#   best / selected), gs_selected_<tag>.tsv (per split and keff), gs_neff_<tag>.tsv, gs_surface_<tag>.png.
.gs_main <- sys.nframe() == 0L
if (!exists("ad_Psym")) source("analysis/tcaclr_depth/adaptive/adaptive_lib.R")

gs_energy <- function(U, Z) { E <- crossprod(U, Z)^2; list(E = E, tail = pmax(colSums(Z^2) - colSums(E), 0)) }

# twice the negative log-likelihood (constants dropped) of experiments with energies E (r x n) and tail energies tl (n), for one (c, s2, t)
gs_nll2 <- function(E, tl, lam, G, c, s2, t) {
  r <- length(lam); v <- s2 + c * lam^t; vt <- s2 + c * lam[r]^t
  sum(log(v)) + (G - r) * log(vt) + colSums(E / v) + tl / vt
}
gs_fit_ct <- function(Ebar, tbar, lam, G, t) {                 # ML (c, s2) at fixed t on the mean energies of the fitting half
  r <- length(lam)
  f <- function(p) { v <- exp(p[2]) + exp(p[1]) * lam^t; vt <- exp(p[2]) + exp(p[1]) * lam[r]^t
    sum(log(v)) + (G - r) * log(vt) + sum(Ebar / v) + tbar / vt }
  s2_0 <- max(tbar / max(G - r, 1), 1e-3); c0 <- max(sum(Ebar) / max(sum(lam^t), 1e-8) / 2, 1e-3)
  o <- stats::optim(c(log(c0), log(s2_0)), f, method = "L-BFGS-B", lower = c(log(1e-6), log(1e-6)), upper = c(log(1e6), log(1e3)))
  c(c = exp(o$par[1]), s2 = exp(o$par[2]))
}

# one direction: graph + fit half Xg, scored half Xh. Returns list(keff -> list(ll = (depths+1) x n_h matrix of per-gene per-experiment log-lik (row 1 = iid), cs = fitted (c, s2) per depth))
gs_direction <- function(Xg, Xh, keffs, depths, r, threads, clipS = TRUE) {
  G <- nrow(Xg); Zg <- ad_rank_normal(Xg); Zh <- ad_rank_normal(Xh)
  f <- ad_fit(Xg, threads); on.exit(f$release())
  tot_g <- colSums(Zg^2); s2_iid <- mean(tot_g) / G; ll_iid <- -0.5 * (G * log(s2_iid) + colSums(Zh^2) / s2_iid) / G
  lapply(stats::setNames(keffs, keffs), function(keff) {
    Ps <- ad_Psym(ad_P_from_fit(f, keff, clipS = clipS)); ev <- ad_eigs(Ps, min(r, G - 5L), vectors = TRUE)
    lam <- pmax(ev$values, 1e-12); eg <- gs_energy(ev$vectors, Zg); eh <- gs_energy(ev$vectors, Zh)
    Ebar <- rowMeans(eg$E); tbar <- mean(eg$tail)
    ll <- matrix(NA_real_, length(depths) + 1L, ncol(Zh)); ll[1L, ] <- ll_iid; cs <- matrix(NA_real_, length(depths), 2L); bnd <- logical(length(depths))
    for (i in seq_along(depths)) { p <- gs_fit_ct(Ebar, tbar, lam, G, depths[i]); cs[i, ] <- p; bnd[i] <- p[["c"]] >= 1e6 * 0.999 || p[["c"]] <= 1e-6 * 1.001 || p[["s2"]] <= 1e-6 * 1.001 || p[["s2"]] >= 1e3 * 0.999   # a fitted parameter sits on its search bound
      ll[i + 1L, ] <- -0.5 * gs_nll2(eh$E, eh$tail, lam, G, p[["c"]], p[["s2"]], depths[i]) / G }
    list(ll = ll, cs = cs, bnd = bnd, lam_r = lam[length(lam)], r = length(lam), tail_share = mean(eh$tail) / mean(colSums(Zh^2)))
  })
}

# one split: both directions pooled; returns the long surface and the selection rows
gs_split <- function(X, h, unit, keffs, depths, r, threads, split = 1L, clipS = TRUE, say = function(...) invisible(NULL)) {
  dirs <- lapply(1:2, function(d) { t0 <- Sys.time(); res <- gs_direction(X[, h == d, drop = FALSE], X[, h != d, drop = FALSE], keffs, depths, r, threads, clipS)
    say(sprintf("split %02d direction %d: %.0f s", split, d, as.numeric(difftime(Sys.time(), t0, units = "secs")))); res })
  uh <- list(unit[h != 1L], unit[h != 2L])                                       # blocks of the held-out half of each direction
  ts <- c(0L, depths); surf <- list(); sel <- list()
  for (kf in as.character(keffs)) {
    B <- do.call(rbind, lapply(1:2, function(d) { ll <- dirs[[d]][[kf]]$ll; b <- tapply(seq_along(uh[[d]]), uh[[d]], function(ix) rowMeans(ll[, ix, drop = FALSE]))
      t(do.call(cbind, b)) }))                                                    # (direction x block) x (t = 0, depths)
    nb <- nrow(B); m <- colMeans(B); i0 <- 1L; cand <- 2:length(ts); ib <- cand[which.max(m[cand])]
    se_d <- vapply(seq_along(ts), function(i) stats::sd(B[, i] - B[, ib]) / sqrt(nb), 0)
    isel <- cand[min(which(m[cand] >= m[ib] - se_d[cand]))]                       # smallest t within one SE of the best (ties at the best itself: se_d = 0)
    se_m <- stats::sd(B[, ib]) / sqrt(nb)                                         # alternative reading: SE of the best's own block-mean score (much wider than the paired SE)
    isel_m <- cand[min(which(m[cand] >= m[ib] - se_m))]
    se_g <- vapply(seq_along(ts), function(i) stats::sd(B[, i] - B[, i0]) / sqrt(nb), 0)
    bf <- c(NA_real_, rowMeans(cbind(dirs[[1]][[kf]]$bnd, dirs[[2]][[kf]]$bnd)))
    cs <- do.call(rbind, lapply(1:2, function(d) dirs[[d]][[kf]]$cs)); cmean <- c(NA_real_, tapply(cs[, 1], rep(seq_along(depths), 2L), mean)); smean <- c(NA_real_, tapply(cs[, 2], rep(seq_along(depths), 2L), mean))
    surf[[kf]] <- data.frame(split = split, keff = as.numeric(kf), t = ts, ll_per_gene_exp = m, gain_vs_iid = m - m[i0], se_gain = se_g, se_diff_best = se_d,
      c_fit = as.numeric(cmean), sigma2_fit = as.numeric(smean), fit_at_bound = bf, best = seq_along(ts) == ib, selected = seq_along(ts) == isel)
    rng <- range(m[cand])
    sel[[kf]] <- data.frame(split = split, keff = as.numeric(kf), blocks = nb, r = dirs[[1]][[kf]]$r, tail_energy_share = mean(c(dirs[[1]][[kf]]$tail_share, dirs[[2]][[kf]]$tail_share)),
      t_best = ts[ib], t_selected_1se = ts[isel], t_selected_1se_marginal = ts[isel_m], se_best_marginal = se_m, gain_best = m[ib] - m[i0], se_gain_best = se_g[ib], ll_range_over_t = rng[2] - rng[1], se_diff_min_nonbest = min(se_d[cand][-(ib - 1L)]))
  }
  list(surface = do.call(rbind, surf), selected = do.call(rbind, sel))
}

# all splits; returns list(surface, selected, neff)
gs_run <- function(X, units, keffs, depths = GL$depths, r = 400L, threads = 8L, workers = 1L, say = function(...) invisible(NULL), clipS = TRUE, neff = TRUE) {
  ns <- length(units$splits)
  out <- ad_lapply(seq_len(ns), function(s) gs_split(X, units$splits[[s]], units$unit, keffs, depths, r, threads, s, clipS, say), workers)
  bad <- vapply(out, function(o) is.character(o), NA); if (any(bad)) stop("split failed: ", out[bad][[1]])
  surf <- do.call(rbind, lapply(out, `[[`, "surface")); sel <- do.call(rbind, lapply(out, `[[`, "selected"))
  ne <- NULL
  if (neff) { a <- gl_neff(X); ne <- data.frame(what = "full matrix", split = 0L, t(a))
    for (s in seq_len(ns)) for (d in 1:2) ne <- rbind(ne, data.frame(what = sprintf("half %d", d), split = s, t(gl_neff(X[, units$splits[[s]] == d, drop = FALSE])))) }
  # pooled over splits: mean surface (descriptive; no SE: the splits share experiments)
  pooled <- stats::aggregate(cbind(ll_per_gene_exp, gain_vs_iid) ~ keff + t, surf, mean)
  list(surface = surf, selected = sel, neff = ne, pooled = pooled)
}

gs_plot <- function(res, file, title = "") {
  s <- res$surface; kf <- sort(unique(s$keff)); cols <- grDevices::hcl.colors(max(length(kf), 2L), "Dark 3")[seq_along(kf)]
  s <- s[s$t > 0, ]; lo <- min(s$gain_vs_iid); hi <- max(s$gain_vs_iid); lo <- min(lo, 0)
  grDevices::png(file, width = 760, height = 520, res = 100); graphics::par(mar = c(4.2, 5, 2.4, 1))
  graphics::plot(range(s$t), c(lo, hi), type = "n", log = "x", xlab = "depth t (log scale)", ylab = "held-out log-likelihood gain over iid model\n(per gene and experiment)", main = title, cex.main = 0.9, xaxt = "n", sub = "thin: splits; thick: mean over splits; dot: 1-SE choice (paired SE)", cex.sub = 0.75)
  graphics::axis(1, at = sort(unique(s$t)), labels = sort(unique(s$t)), cex.axis = 0.7); graphics::abline(h = 0, col = "grey60", lty = 3)
  for (i in seq_along(kf)) {
    ss <- s[s$keff == kf[i], ]; for (sp in unique(ss$split)) { z <- ss[ss$split == sp, ]; graphics::lines(z$t, z$gain_vs_iid, col = grDevices::adjustcolor(cols[i], 0.35), lwd = 1) }
    mm <- stats::aggregate(gain_vs_iid ~ t, ss, mean); graphics::lines(mm$t, mm$gain_vs_iid, col = cols[i], lwd = 2.5)
    se <- res$selected[res$selected$keff == kf[i], ]; graphics::points(se$t_selected_1se, vapply(seq_len(nrow(se)), function(j) ss$gain_vs_iid[ss$split == se$split[j] & ss$t == se$t_selected_1se[j]], 0), pch = 19, col = cols[i], cex = 1.2)
  }
  graphics::legend("bottomright", bty = "n", cex = 0.8, legend = sprintf("keff %g", kf), col = cols, lwd = 2.5)
  grDevices::dev.off(); invisible(file)
}

gs_write <- function(res, out, tag) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE); w <- function(x, nm) utils::write.table(x, file.path(out, sprintf("%s_%s.tsv", nm, tag)), sep = "\t", quote = FALSE, row.names = FALSE)
  w(res$surface, "gs_surface"); w(res$selected, "gs_selected"); if (!is.null(res$neff)) w(res$neff, "gs_neff"); w(res$pooled, "gs_pooled")
  gs_plot(res, file.path(out, sprintf("gs_surface_%s.png", tag)), tag)
}

if (.gs_main) {
  opt <- gl_options(commandArgs(TRUE), list(matrix = "", dataset = "X", kind = "lab", splits = 1L, genes = 0L, keffs = "10", depths = "", r = 400L, clusters = 0L, threads = 8L, workers = 1L,
    seed = AD$seed, out = "", tag = "", lab_splits = "results/lab_splits/split_assignments.csv", clipS = 1L))
  if (!nzchar(opt$matrix) || !nzchar(opt$out)) stop("--matrix and --out are required"); if (!nzchar(opt$tag)) opt$tag <- opt$dataset
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, sprintf("gs_%s.log", opt$tag)))
  say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
  X <- ad_load(opt$matrix, opt$genes, opt$seed); say(sprintf("%s: %d genes x %d experiments", opt$dataset, nrow(X), ncol(X)))
  units <- ad_units(X, opt$kind, opt$splits, opt$seed, opt$lab_splits, opt$clusters); say(sprintf("%d resampling units (%s)", units$k, opt$kind))
  keffs <- as.numeric(strsplit(opt$keffs, ",")[[1]]); depths <- if (nzchar(opt$depths)) as.integer(strsplit(opt$depths, ",")[[1]]) else GL$depths
  t0 <- Sys.time(); res <- gs_run(X, units, keffs, depths, opt$r, opt$threads, opt$workers, say, opt$clipS == 1L); gs_write(res, opt$out, opt$tag)
  print(res$selected, digits = 4, row.names = FALSE); print(res$neff[res$neff$split == 0, ], row.names = FALSE, digits = 4)
  say(sprintf("done: %.1f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), ad_mem()))
}
