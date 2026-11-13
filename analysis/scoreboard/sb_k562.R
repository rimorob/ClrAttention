#!/usr/bin/env Rscript
# Scoreboard (docs/plans/scoreboard_plan.md), K562 job. DESCRIPTIVE: point estimates, no intervals, no tests.
# New cells: DoRothEA A-B-C (raw, hub pi, hub MI) for CLR, attention (K10, t = 9), a3r, a3r2, a3r3. Set-up, truths and stratifiers are those of
# analysis/explore_ideas47/ideas47_lib.R (sourced unchanged; jack = 1 because no jackknife is used here). Gate first: the ten stored truths of
# results/explore_stack/k562/pooled_levels.csv x 3 scopes x 5 readouts are reproduced (1e-9), and DoRothEA ABC attention / CLR reproduce
# results/baselines_k562ess_truthgenes_genie3/metrics.csv (raw) and its hub-stratified rows; the job aborts before writing results on failure.
# Usage (repo root): Rscript analysis/scoreboard/sb_k562.R [--out results/scoreboard/k562] [--threads 8]
args <- commandArgs(TRUE); opt <- list(out = "results/scoreboard/k562", threads = 8L, workers = 4L, jack = 1L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
sb_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(sb_dir, "../explore_ideas47/ideas47_lib.R"))      # data, truths, K10 (S9), CLR (C), a3r, H_att, H_mi, Pt9; defines say, opt$out/cache
source(file.path(sb_dir, "../explore_stack/stack_lib.R")); source(file.path(sb_dir, "sb_lib.R"))
unlink(file.path(opt$out, "cache"), recursive = TRUE)             # not used here
abc <- dorothea_sets(genes, c("A", "B", "C"))
bench$DoRothEA_ABC <- list(reg = abc, pairs = comembership_pairs(abc, none, genes))
say(sprintf("truth DoRothEA_ABC %5d sets, %9d pairs, base rate %.4f", length(abc$members), length(bench$DoRothEA_ABC$pairs$label), mean(bench$DoRothEA_ABC$pairs$label)))
rm(halves); invisible(gc(verbose = FALSE))

S1 <- a3r(Pt9); rm(Pt9); invisible(gc(verbose = FALSE))
L1 <- sk_layer(C, ts_out = 9L); d_e <- max(abs(L1$outs[["9"]] - S1)); rm(L1); say(sprintf("generic stack N = 1 vs a3r: max|diff| = %.3g", d_e)); stopifnot(d_e < 1e-12)
L2 <- sk_layer(S1, ts_out = 9L); S2 <- L2$outs[["9"]]; rm(L2); invisible(gc(verbose = FALSE))
L3 <- sk_layer(S2, ts_out = 9L); S3 <- L3$outs[["9"]]; rm(L3); invisible(gc(verbose = FALSE))
SS <- list(CLR = C, attention = S9, a3r = S1, a3r2 = S2, a3r3 = S3)
pooled <- do.call(rbind, lapply(names(SS), function(m) { r <- sb_pooled(m, SS[[m]], bench, H_att, H_mi); say("scored ", m); r }))

## ---- equivalence gate -------------------------------------------------------------------------------------------------------
rows <- list(); pl <- utils::read.csv("results/explore_stack/k562/pooled_levels.csv", stringsAsFactors = FALSE)
mp2 <- c(attention = "K10", CLR = "CLR", a3r = "a3r", a3r2 = "a3r2", a3r3 = "a3r3"); tr10 <- setdiff(unique(pl$truth), character(0))
for (m in names(mp2)) { d <- 0
  for (sc in c("raw", "hub_att", "hub_mi")) for (tn in tr10) { mine <- pooled$aupr[pooled$readout == m & pooled$truth == tn & pooled$scope == sc]; ref <- pl$aupr[pl$readout == mp2[[m]] & pl$scope == sc & pl$truth == tn]
    stopifnot(length(mine) == 1L, length(ref) == 1L); d <- max(d, abs(mine - ref)) }
  rows <- sb_gate_add(rows, sprintf("%s: AUPR on %d stored truths x 3 scopes vs explore_stack/k562/pooled_levels.csv (%s)", m, length(tr10), mp2[[m]]), d) }
bm <- utils::read.csv("results/baselines_k562ess_truthgenes_genie3/metrics.csv", stringsAsFactors = FALSE)
for (m in c("attention", "CLR")) { key <- if (m == "CLR") "clr" else m
  for (tn in c("DoRothEA_ABC", "DoRothEA_AB")) { bs <- bm[bm$truth == tn, ]; mine <- pooled[pooled$readout == m & pooled$truth == tn, ]
    rows <- sb_gate_add(rows, sprintf("%s %s raw AUPR vs baselines metrics.csv", m, tn), abs(mine$aupr[mine$scope == "raw"] - bs$aupr[bs$method == key]))
    rows <- sb_gate_add(rows, sprintf("%s %s hub_att AUPR vs baselines metrics.csv (%s_hubstrat)", m, tn, key), abs(mine$aupr[mine$scope == "hub_att"] - bs$aupr[bs$method == paste0(key, "_hubstrat")]), 1e-9) } }
sb_gate_finish(rows, file.path(opt$out, "equivalence_gate.csv"), say)
utils::write.csv(pooled, file.path(opt$out, "pooled.csv"), row.names = FALSE)
say("DoRothEA ABC enrichment:"); print(pooled[pooled$truth == "DoRothEA_ABC", ], digits = 4, row.names = FALSE)
say("done")
