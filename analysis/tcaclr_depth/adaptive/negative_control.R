#!/usr/bin/env Rscript
# Negative control for the label-free rules (plan Revision 2, items 9 and 11; review section 4.4a). EXPLORATORY; no truth is read.
#   Rscript analysis/tcaclr_depth/adaptive/negative_control.R --matrix results/baselines_m3d/expression.tsv --dataset M3D --kind lab --perms 4 --keffs 10,28 --out results/tcaclr_depth/adaptive_m3d
#
# --seperms 0 (the job setting) takes the edge from the other controls instead of fresh permutations (see nc_one); --seperms > 0 is the nested version.
# Each control matrix is the real matrix with every gene's values permuted across experiments independently (--perms of them, seeds fixed). Each is run through
# spectral_edge.R (its own --seperms fresh permutations define the edge; observed = the control matrix) and graph_signal_lik.R (same resampling units and splits as
# the real-data run). The rule correctly reports "no structure" for a (control, keff) when
#   edge:  no observed non-trivial eigenvalue above the LARGEST permuted edge (n_above_max == 0; n_above_mean is also reported: with P permutations a fresh null draw
#          exceeds the largest of P with probability 1/(P+1), and the mean edge about half of the time, so these are chance exceedances, not structure), and
#   gs:    the held-out likelihood gain over the iid model at the best t is <= 0 in every split (gain_best <= 0), the surface being flat or falling.
# Alternative (reported alongside, never replacing the above): flag_edge_p = permutation p-value of the control's lambda_2 > 0.05 (needs --seperms >= 19), and
# flag_gs_1se = gain_best <= its standard error (block SE) in every split; flag_both_alt = both.
# Flags are mechanical (the thresholds are the ones written here); nothing is interpreted. Output (<out>/): nc_perm_<tag>.tsv (one row per control and keff),
# nc_summary_<tag>.tsv (share flagged per keff and overall).
.nc_main <- sys.nframe() == 0L
if (!exists("ad_Psym")) source("analysis/tcaclr_depth/adaptive/adaptive_lib.R")
if (!exists("se_run")) source("analysis/tcaclr_depth/adaptive/spectral_edge.R")
if (!exists("gs_run")) source("analysis/tcaclr_depth/adaptive/graph_signal_lik.R")

# One control. seperms > 0: the edge comes from `seperms` fresh permutations of the control (nested, the literal "control through spectral_edge.R").
# seperms = 0: only the control's own spectrum is computed here and the edge is taken, in nc_run, from the OTHER controls (leave-one-out; the controls are exchangeable
# draws from the same permutation null, so this is the same edge at 1/(seperms + 1) of the cost).
nc_one <- function(X, p, units, keffs, depths, r, seperms, kobs, threads, seed, say) {
  t0 <- Sys.time(); Xp <- ad_permute_genes(X, seed + 9000L + p)
  se <- if (seperms > 0L) se_run(Xp, keffs, seperms, kobs, 30L, threads, 1L, seed + 100L * p, paste0("perm", p))$summary else NULL
  obs <- if (seperms > 0L) NULL else se_spectra(Xp, keffs, kobs, threads)
  gs <- gs_run(Xp, units, keffs, depths, r, threads, 1L, neff = FALSE)
  gs$selected$gap <- gs$selected$gain_best - gs$selected$se_gain_best
  gmax <- stats::aggregate(gain_best ~ keff, gs$selected, max); gap <- stats::aggregate(gap ~ keff, gs$selected, max); tsel <- stats::aggregate(t_selected_1se ~ keff, gs$selected, function(x) stats::median(x))
  rng <- stats::aggregate(ll_range_over_t ~ keff, gs$selected, max)
  kf <- gmax$keff; m <- function(v) v[match(kf, gmax$keff)]
  r1 <- data.frame(control = p, keff = kf, gs_gain_best_max_over_splits = gmax$gain_best, gs_gain_minus_se_max = gap$gap[match(kf, gap$keff)], gs_t_selected_median = tsel$t_selected_1se[match(kf, tsel$keff)],
    gs_range_over_t_max = rng$ll_range_over_t[match(kf, rng$keff)])
  if (!is.null(se)) { j <- match(kf, se$keff); r1$edge_mean <- se$edge_mean[j]; r1$edge_max <- se$edge_max[j]; r1$obs_lambda2 <- se$obs_lambda2[j]; r1$edge_perm_p <- se$obs_lambda2_perm_p[j]
    r1$n_above_mean <- se$n_above_mean[j]; r1$n_above_max <- se$n_above_max[j] }
  say(sprintf("control %d: %.0f s", p, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(rows = r1, obs = lapply(obs, function(o) o$values))
}
nc_run <- function(X, units, keffs, depths, perms, r, seperms, kobs, threads, workers, seed, say = function(...) invisible(NULL)) {
  out <- ad_lapply(seq_len(perms), function(p) nc_one(X, p, units, keffs, depths, r, seperms, kobs, threads, seed, say), workers)
  bad <- vapply(out, is.character, NA); if (any(bad)) stop("control failed: ", out[bad][[1]])
  d <- do.call(rbind, lapply(out, `[[`, "rows"))
  if (seperms == 0L) {                                                 # leave-one-out edge from the other controls
    stopifnot(perms >= 3L); add <- list()
    for (kf in as.character(unique(d$keff))) {
      l2 <- vapply(out, function(o) o$obs[[kf]][2L], 0)
      for (p in seq_len(perms)) { oth <- l2[-p]; v <- out[[p]]$obs[[kf]]; ob <- v[-1L]
        add[[length(add) + 1L]] <- data.frame(control = p, keff = as.numeric(kf), edge_mean = mean(oth), edge_max = max(oth), obs_lambda2 = v[2L], edge_perm_p = (1 + sum(oth >= v[2L])) / perms,
          n_above_mean = sum(ob > mean(oth)), n_above_max = sum(ob > max(oth))) }
    }
    d <- merge(d, do.call(rbind, add), by = c("control", "keff"))
  }
  d$flag_edge <- d$n_above_max == 0L; d$flag_gs <- d$gs_gain_best_max_over_splits <= 0; d$flag_both <- d$flag_edge & d$flag_gs
  d$flag_edge_p <- d$edge_perm_p > 0.05; d$flag_gs_1se <- d$gs_gain_minus_se_max <= 0; d$flag_both_alt <- d$flag_edge_p & d$flag_gs_1se
  d <- d[order(d$keff, d$control), ]
  sm <- function(z, keff) data.frame(keff = keff, controls = nrow(z), share_edge_no_structure = mean(z$flag_edge), share_gs_no_structure = mean(z$flag_gs), share_both = mean(z$flag_both),
    share_n_above_mean_zero = mean(z$n_above_mean == 0L), share_edge_p_gt05 = mean(z$flag_edge_p), share_gs_gain_le_se = mean(z$flag_gs_1se), share_both_alt = mean(z$flag_both_alt),
    max_n_above_max = max(z$n_above_max), max_gs_gain = max(z$gs_gain_best_max_over_splits))
  s <- rbind(do.call(rbind, lapply(split(d, d$keff), function(z) sm(z, z$keff[1]))), sm(d, NA_real_))
  list(perm = d, summary = s)
}

if (.nc_main) {
  opt <- gl_options(commandArgs(TRUE), list(matrix = "", dataset = "X", kind = "lab", splits = 1L, genes = 0L, keffs = "10", depths = "", r = 400L, clusters = 0L, perms = 4L, seperms = 3L, kobs = 40L,
    threads = 8L, workers = 1L, seed = AD$seed, out = "", tag = "", lab_splits = "results/lab_splits/split_assignments.csv"))
  if (!nzchar(opt$matrix) || !nzchar(opt$out)) stop("--matrix and --out are required"); if (!nzchar(opt$tag)) opt$tag <- opt$dataset
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, sprintf("nc_%s.log", opt$tag)))
  say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
  X <- ad_load(opt$matrix, opt$genes, opt$seed); say(sprintf("%s: %d genes x %d experiments", opt$dataset, nrow(X), ncol(X)))
  units <- ad_units(X, opt$kind, opt$splits, opt$seed, opt$lab_splits, opt$clusters)
  keffs <- as.numeric(strsplit(opt$keffs, ",")[[1]]); depths <- if (nzchar(opt$depths)) as.integer(strsplit(opt$depths, ",")[[1]]) else GL$depths
  t0 <- Sys.time(); res <- nc_run(X, units, keffs, depths, opt$perms, opt$r, opt$seperms, opt$kobs, opt$threads, opt$workers, opt$seed, say)
  w <- function(x, nm) utils::write.table(x, file.path(opt$out, sprintf("nc_%s_%s.tsv", nm, opt$tag)), sep = "\t", quote = FALSE, row.names = FALSE); w(res$perm, "perm"); w(res$summary, "summary")
  print(res$perm, digits = 4, row.names = FALSE); print(res$summary, digits = 3, row.names = FALSE)
  say(sprintf("done: %.1f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), ad_mem()))
}
