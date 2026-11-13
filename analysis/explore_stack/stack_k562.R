#!/usr/bin/env Rscript
# EXPLORATORY, IN-SAMPLE, NOT PREREGISTERED (docs/plans/explore_stack_plan.md). Stacked a3r layers on K562.
# Equivalence gate first (job aborts on failure), then a3r2, a3r3, the t2 grid, controls c1, c2 and layer diagnostics, each scored on the ten
# K562 truths in three scopes (raw, hub pi, hub MI row mean). Set-up, truths, jackknife streams and scoring are those of
# analysis/explore_ideas47/ideas47_lib.R (sourced unchanged).
# Usage (repo root): Rscript analysis/explore_stack/stack_k562.R [--out results/explore_stack/k562] [--threads 8] [--workers 8] [--jack 100] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(out = "results/explore_stack/k562", threads = 8L, workers = 8L, jack = 100L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
stk_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(stk_dir, "../explore_ideas47/ideas47_lib.R"))     # data, truths, jackknife streams, K10 / CLR / a3r, stratifiers, score_readout, e47_gate
source(file.path(stk_dir, "stack_lib.R"))
DS <- "k562"; diag_rows <- list()
say(sprintf("stack k562: G = %d, workers %d, threads %d", G, nc, opt$threads))

## ---- references scored first, then the equivalence gate -----------------------------------------------------
res_clr <- score_readout("CLR", C, macro = FALSE); res_k10 <- score_readout("K10", S9, macro = FALSE)
S1 <- a3r(Pt9); res_a3r <- score_readout("a3r", S1, macro = FALSE)
if (!QUICK) e47_gate(res_clr, res_k10, res_a3r) else say("QUICK run: gate (a)-(c) skipped (gene subset)")
rows <- list()
if (!QUICK) {      # (d) a3r minus K10: point difference and jackknife interval against the full-precision record of corum_analysis.md
  mr <- utils::read.csv("analysis/corum_exploratory/metrics_robust.csv", stringsAsFactors = FALSE); mr <- mr[mr$variant == "a_clr_bilateral_logP_t9_robust", ]
  mp <- c(DoRothEA_AB = "DoRothEA_AB", CollecTRI = "CollecTRI", CORUM = "CORUM_all", KD_response = "KD_response", CORUM_noRP = "CORUM_noRP")
  for (tn in names(mp)) { e <- mp[[tn]]; d <- res_a3r$raw[[e]]$aupr - res_k10$raw[[e]]$aupr; s <- stats::sd(res_a3r$raw[[e]]$jack - res_k10$raw[[e]]$jack); r <- mr[mr$truth == tn, ]
    rows <- sk_gate_add(rows, sprintf("(d) a3r minus K10, %s: d, lo, hi vs metrics_robust.csv (d_vs_att9, att9_lo, att9_hi)", tn),
      max(abs(d - r$d_vs_att9), abs(d - 1.96 * s - r$att9_lo), abs(d + 1.96 * s - r$att9_hi)))
    if (tn == "CORUM") { pr <- sprintf("%+.3f [%+.3f, %+.3f]", d, d - 1.96 * s, d + 1.96 * s)
      rows <- sk_gate_add(rows, sprintf("(d) CORUM a3r minus K10 prints %s, document table +0.059 [+0.033, +0.084]", pr), as.numeric(pr != "+0.059 [+0.033, +0.084]")) } }
}
## (e) generic stacking code with N = 1 equals a3r of the package path
L1 <- sk_layer(C, ts_out = 9L); d_e <- max(abs(L1$outs[["9"]] - S1)); d_A <- max(abs(L1$A - A))
rows <- sk_gate_add(rows, "(e) generic stack N = 1 (operator from CLR scores, walk, a3r) vs a3r of the package path: max |S1' - S1|", d_e, 1e-12)
rows <- sk_gate_add(rows, "(e) generic stack N = 1 operator vs package operator: max |A' - A|", d_A, 1e-12)
sk_gate_finish(rows, file.path(opt$out, "equivalence_gate_stack.csv"), say)
rm(L1); invisible(gc(verbose = FALSE))

## ---- layers ---------------------------------------------------------------------------------------------------
top1_C <- sk_top1(C, UT); top1_1 <- sk_top1(S1, UT)
Lr1 <- sk_layer(A = A, ts_out = 9L)       # layer 1 (reuses the package-path operator): diagnostics
diag_rows[[1]] <- sk_diag(DS, 1L, C, A, Lr1$P, Lr1$reach[["9"]], S1, top1_C, top1_C, mi_rm, q_att, UT); rm(Lr1)
utils::write.csv(do.call(rbind, diag_rows), file.path(opt$out, "diagnostics_layers.csv"), row.names = FALSE)

L2 <- sk_layer(S1, ts_out = c(2L, 4L, 9L)); S2 <- L2$outs[["9"]]
score_readout("a3r2", S2, macro = FALSE); score_readout("a3r2_t2_2", L2$outs[["2"]], macro = FALSE); score_readout("a3r2_t2_4", L2$outs[["4"]], macro = FALSE)
diag_rows[[2]] <- sk_diag(DS, 2L, S1, L2$A, L2$P, L2$reach[["9"]], S2, top1_1, top1_C, mi_rm, q_att, UT)
t2d <- do.call(rbind, lapply(c("2", "4"), function(tt) { x <- L2$outs[[tt]]; data.frame(dataset = DS, layer = 2L, t = as.integer(tt), unreach_share = L2$reach[[tt]][["unreach_share"]],
  rows_lt2_reachable = L2$reach[[tt]][["rows_lt2_reachable"]], S_frac_nonzero = mean(x[UT] > 0), S_rows_all_zero = sum(rowSums(x > 0) == 0L)) }))
utils::write.csv(t2d, file.path(opt$out, "diagnostics_t2_grid.csv"), row.names = FALSE)
rm(L2); utils::write.csv(do.call(rbind, diag_rows), file.path(opt$out, "diagnostics_layers.csv"), row.names = FALSE); invisible(gc(verbose = FALSE))

Sprev <- S2; top1_prev <- sk_top1(S2, UT); rm(S2)
for (k in 3:5) {
  Lk <- sk_layer(Sprev, ts_out = 9L); Sk <- Lk$outs[["9"]]
  if (k == 3L) score_readout("a3r3", Sk, macro = FALSE)       # layers 4 and 5 are label-free diagnostics only: never scored against a truth
  diag_rows[[k]] <- sk_diag(DS, k, Sprev, Lk$A, Lk$P, Lk$reach[["9"]], Sk, top1_prev, top1_C, mi_rm, q_att, UT)
  top1_prev <- sk_top1(Sk, UT); Sprev <- Sk; rm(Lk, Sk); invisible(gc(verbose = FALSE))
  utils::write.csv(do.call(rbind, diag_rows), file.path(opt$out, "diagnostics_layers.csv"), row.names = FALSE)
}
rm(Sprev, top1_prev)

## ---- controls -------------------------------------------------------------------------------------------------
Lc <- sk_layer(A = A, ts_out = c(18L, 27L))
score_readout("c1_t18", Lc$outs[["18"]], macro = FALSE); score_readout("c1_t27", Lc$outs[["27"]], macro = FALSE)
c1d <- do.call(rbind, lapply(c("18", "27"), function(tt) { x <- Lc$outs[[tt]]; data.frame(dataset = DS, t = as.integer(tt), unreach_share = Lc$reach[[tt]][["unreach_share"]],
  rows_lt2_reachable = Lc$reach[[tt]][["rows_lt2_reachable"]], S_frac_nonzero = mean(x[UT] > 0), S_rows_all_zero = sum(rowSums(x > 0) == 0L),
  spearman_Srowmean_mi = stats::cor(rowMeans(x), mi_rm, method = "spearman"), jaccard_top1_a3r = sk_jaccard(sk_top1(x, UT), top1_1)) }))
utils::write.csv(c1d, file.path(opt$out, "diagnostics_c1_depth.csv"), row.names = FALSE); rm(Lc); invisible(gc(verbose = FALSE))
score_readout("c2_z", sk_c2z(S1), macro = FALSE); score_readout("c2_log", sk_c2log(S1), macro = FALSE)
saveRDS(list(base_rate = base_rate, truths = TRUTHS, opt = opt), file.path(opt$out, "cache", "meta.rds"))
say("done")
