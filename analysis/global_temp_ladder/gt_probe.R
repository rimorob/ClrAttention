#!/usr/bin/env Rscript
# Probe (no truth read): effective-partner distribution of the global-temperature softmax over a wide log2 tau grid, to fix the ladder endpoints in the plan.
#   Rscript analysis/global_temp_ladder/gt_probe.R --ds m3d|ecomac [--threads 8] [--out results/global_temp_ladder/probe]
suppressPackageStartupMessages(source("analysis/global_temp_ladder/gt_lib.R"))
O <- gl_options(commandArgs(TRUE), list(ds = "", threads = 8L, out = file.path(GTL$out, "probe")))
dir.create(O$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(O$out, sprintf("probe_%s.log", O$ds)))
X <- wl_matrix(O$ds); f <- gl_fit(X, O$threads); C <- f$clr_scores; f$release(); rm(f, X)
Tm <- gt_topmat(C); tau2 <- gt_tau_for_eff(Tm); say(sprintf("%s: G = %d, tau(mean eff 10) = %.6f (log2 %.3f); rows with candidates %d", O$ds, nrow(C), tau2, log2(tau2), sum(!is.na(Tm$idx[, 1]))))
grid <- 2^seq(-8, 4, by = 0.5); r <- do.call(rbind, lapply(grid, function(tau) gt_eff_row(Tm, tau, gt_tau_id(tau), O$ds)))
r <- rbind(r, gt_eff_row(Tm, tau2, "V2", O$ds)); utils::write.csv(r, file.path(O$out, sprintf("probe_%s.csv", O$ds)), row.names = FALSE); print(r[, c("op", "tau", "eff_mean", "eff_median", "eff_q10", "eff_q90")], digits = 4, row.names = FALSE)
say(sprintf("peak RSS %.2f GB", gl_vmhwm_gb()))
