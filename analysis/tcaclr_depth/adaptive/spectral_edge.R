#!/usr/bin/env Rscript
# Spectral-edge depth rule (plan Revision 2, item 8a). EXPLORATORY; no truth is read.
#   Rscript analysis/tcaclr_depth/adaptive/spectral_edge.R --matrix results/baselines_m3d/expression.tsv --dataset M3D --keffs 10,28 --perms 5 --out results/tcaclr_depth/adaptive_m3d
#
# For a matrix X (genes x experiments): fit MI + CLR, build the aCLR operator at (keff, cap = gl_cap(keff)) with the grid library, and take its
# symmetrised-by-degree lazy form P_sym (adaptive_lib.R). The NOISE EDGE lambda_edge is the largest eigenvalue other than the trivial 1 of the same
# operator built from PER-GENE PERMUTED data (each gene's values permuted across experiments independently; --perms permutations, mean and range
# reported). Depth rule t* = ln(eps) / ln(lambda_edge) for eps = 0.05 and 0.01 (t* at the mean edge, and at the smallest / largest edge of the range).
# obs_lambda2_perm_p = (1 + #permutations with lambda_2 >= the observed lambda_2) / (P + 1): the permutation p-value of the observed second eigenvalue (needs P >= 19 to reach 0.05).
# Also: the observed top eigenvalues (top 30 stored, --kobs computed), the number of observed non-trivial eigenvalues above the permuted edge
# (mean edge and largest edge), and the second-largest modulus of the NON-symmetric P as a cross-check of the symmetrisation.
# Output (<out>/): spectral_edge_<tag>.tsv (one row per keff), spectral_edge_<tag>_spectrum.tsv (long: keff, rank, observed, perm_mean, perm_min, perm_max,
#   top 30), spectral_edge_<tag>.png (observed vs permuted spectrum, one panel per keff).
.se_main <- sys.nframe() == 0L            # TRUE only under Rscript; the functions below are also sourced by planted_sim.R and negative_control.R
if (!exists("ad_Psym")) source("analysis/tcaclr_depth/adaptive/adaptive_lib.R")

## one matrix -> per-keff top eigenvalues (and the non-symmetric cross-check)
se_spectra <- function(X, keffs, k, threads, clipS = TRUE) {
  f <- ad_fit(X, threads); on.exit(f$release())
  lapply(stats::setNames(keffs, keffs), function(keff) {
    P <- ad_P_from_fit(f, keff, clipS = clipS); Ps <- ad_Psym(P)
    list(values = ad_eigs(Ps, k)$values, nonsym = ad_edge_nonsym(P), nnz = length(P@x), isolated = sum(Matrix::diag(P) > 0.5 + 1e-12 & Matrix::rowSums(P > 0) == 1L))
  })
}

# observed spectrum vs permuted edge. Returns list(summary, spectrum).
se_run <- function(X, keffs, perms = 5L, kobs = 100L, kperm = 30L, threads = 8L, workers = 1L, seed = AD$seed, dataset = "X", clipS = TRUE, say = function(...) invisible(NULL)) {
  t0 <- Sys.time(); G <- nrow(X); N <- ncol(X)
  # task 0 = the observed matrix, tasks 1..perms = per-gene permutations. All clr work happens in the tasks: a process that has run clr (OpenMP) and then forks
  # (mclapply) hangs its children, so the calling process must never fit before it forks (see ad_lapply).
  all <- ad_lapply(0:perms, function(p) se_spectra(if (p == 0L) X else ad_permute_genes(X, seed + 7000L + p), keffs, if (p == 0L) kobs else kperm, threads, clipS), workers)
  bad <- vapply(all, is.character, NA); if (any(bad)) stop("spectrum task failed: ", all[bad][[1]])
  obs <- all[[1L]]; pm <- all[-1L]
  say(sprintf("observed + %d permuted spectra: %.0f s", perms, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rows <- list(); spec <- list()
  for (kf in as.character(keffs)) {
    o <- obs[[kf]]$values; pv <- do.call(cbind, lapply(pm, function(x) x[[kf]]$values[seq_len(kperm)]))     # kperm x perms
    edges <- pv[2L, ]; ns <- vapply(pm, function(x) x[[kf]]$nonsym, 0)                                      # per-permutation largest non-trivial
    em <- mean(edges); elo <- min(edges); ehi <- max(edges); ob <- o[-1L]
    n_unit_o <- sum(o > 1 - AD$unit_tol); n_unit_p <- sum(pv[2L, ] > 1 - AD$unit_tol)
    r <- data.frame(dataset = dataset, G = G, N = N, keff = as.numeric(kf), cap = gl_cap(as.numeric(kf)), perms = perms,
      edge_mean = em, edge_min = elo, edge_max = ehi, edge_nonsym_mean = mean(ns, na.rm = TRUE), obs_lambda2 = o[2L], obs_lambda2_perm_p = (1 + sum(edges >= o[2L])) / (perms + 1), obs_nonsym2 = obs[[kf]]$nonsym,
      kobs = length(o), n_above_mean = sum(ob > em), n_above_max = sum(ob > ehi), n_above_capped = sum(ob > ehi) >= length(ob),
      n_unit_obs = n_unit_o, n_unit_perm = n_unit_p, isolated_obs = obs[[kf]]$isolated, nnz_obs = obs[[kf]]$nnz)
    for (e in AD$eps_edge) { nm <- sprintf("tstar_eps%s", sub("^0\\.", "", format(e)))
      r[[nm]] <- ad_tstar(em, e); r[[paste0(nm, "_lo")]] <- ad_tstar(elo, e); r[[paste0(nm, "_hi")]] <- ad_tstar(ehi, e) }
    top <- o[seq_len(min(30L, length(o)))]; names(top) <- sprintf("lam%02d", seq_along(top)); r <- cbind(r, as.data.frame(t(top)))
    rows[[kf]] <- r
    k <- seq_len(min(30L, length(o)))
    spec[[kf]] <- data.frame(dataset = dataset, keff = as.numeric(kf), rank = k, observed = o[k], perm_mean = rowMeans(pv)[k], perm_min = apply(pv, 1, min)[k], perm_max = apply(pv, 1, max)[k])
  }
  list(summary = do.call(rbind, rows), spectrum = do.call(rbind, spec))
}

se_plot <- function(res, file, title = "") {
  sp <- res$spectrum; kf <- sort(unique(sp$keff)); sm <- res$summary
  nr <- ceiling(length(kf) / min(4L, length(kf))); nc <- min(4L, length(kf))
  grDevices::png(file, width = 340 * nc, height = 300 * nr + 40, res = 100)
  graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2.2, 0.8), oma = c(0, 0, 1.5, 0))
  for (k in kf) {
    s <- sp[sp$keff == k, ]; m <- sm[sm$keff == k, ]
    graphics::plot(s$rank, s$observed, ylim = c(max(0, min(s$observed, s$perm_min) - 0.12), 1.02), type = "n", xlab = "eigenvalue index (1 = trivial)", ylab = "eigenvalue of lazy operator", main = sprintf("keff %g", k), cex.main = 0.95)
    graphics::rect(0.5, m$edge_min, max(s$rank) + 0.5, m$edge_max, col = grDevices::adjustcolor("#d55e00", 0.18), border = NA)
    graphics::abline(h = m$edge_mean, col = "#d55e00", lty = 2)
    graphics::segments(s$rank, s$perm_min, s$rank, s$perm_max, col = "#d55e00")
    graphics::points(s$rank, s$perm_mean, pch = 1, col = "#d55e00", cex = 0.8)
    graphics::points(s$rank, s$observed, pch = 16, col = "#0072b2", cex = 0.8)
    graphics::legend("bottomleft", bty = "n", cex = 0.7, legend = c(sprintf("above edge (max): %d", m$n_above_max), sprintf("edge %.3f (%.3f-%.3f)", m$edge_mean, m$edge_min, m$edge_max), sprintf("t*(0.05) = %.1f", m$tstar_eps05)))
  }
  graphics::mtext(sprintf("%s: observed (blue) vs per-gene permuted (orange; circle = mean, bar = range)", title), outer = TRUE, cex = 0.85)
  grDevices::dev.off(); invisible(file)
}

se_write <- function(res, out, tag) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res$summary, file.path(out, sprintf("spectral_edge_%s.tsv", tag)), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(res$spectrum, file.path(out, sprintf("spectral_edge_%s_spectrum.tsv", tag)), sep = "\t", quote = FALSE, row.names = FALSE)
  se_plot(res, file.path(out, sprintf("spectral_edge_%s.png", tag)), tag)
}

if (.se_main) {
  opt <- gl_options(commandArgs(TRUE), list(matrix = "", dataset = "X", genes = 0L, keffs = "10", perms = 5L, kobs = 100L, threads = 8L, workers = 1L, seed = AD$seed, out = "", tag = "", clipS = 1L))
  if (!nzchar(opt$matrix) || !nzchar(opt$out)) stop("--matrix and --out are required"); if (!nzchar(opt$tag)) opt$tag <- opt$dataset
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, sprintf("spectral_edge_%s.log", opt$tag)))
  say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
  X <- ad_load(opt$matrix, opt$genes, opt$seed); say(sprintf("%s: %d genes x %d experiments", opt$dataset, nrow(X), ncol(X)))
  t0 <- Sys.time(); keffs <- as.numeric(strsplit(opt$keffs, ",")[[1]])
  res <- se_run(X, keffs, opt$perms, opt$kobs, 30L, opt$threads, opt$workers, opt$seed, opt$dataset, opt$clipS == 1L, say)
  se_write(res, opt$out, opt$tag)
  print(res$summary[, c("keff", "edge_mean", "edge_min", "edge_max", "edge_nonsym_mean", "obs_lambda2", "obs_nonsym2", "n_above_mean", "n_above_max", "tstar_eps05", "tstar_eps01")], digits = 4, row.names = FALSE)
  say(sprintf("done: %.1f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), ad_mem()))
}
