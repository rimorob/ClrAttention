#!/usr/bin/env Rscript
# Gate Ga, K562 part (docs/plans/complex_confirm_design.md, Gates, Ga part 2): the code path that the new driver reuses (analysis/explore_stack/
# stack_lib.R::sk_a3r, ideas47_lib.R set-up, truths, jackknife streams and scoring) must reproduce, on K562 and from scratch (no explore_stack cache is
# read), the stored values named in the design. Spent data only; no net4 and no yeast truth is touched.
# Rows: (a) ideas47 gate vs robust_z, baselines, corum_expert2 caches; (d) a3r - K10 point differences and jackknife intervals vs metrics_robust.csv
# (CORUM line prints +0.059 [+0.033, +0.084]; CORUM noRP: d only); (e) generic stack N = 1 equals the package path; and the design's printed values
# (a3r CORUM_raw AUPR 0.105435083 / 7.95x; CLR 0.0847 / 6.39x; attention 0.0465 / 3.51x).
# Usage (repo root): Rscript analysis/complex_confirm/gate_a_k562.R [--out results/complex_confirm/gate_a/k562] [--threads 8] [--workers 8]
args <- commandArgs(TRUE)
opt <- list(out = "results/complex_confirm/gate_a/k562", threads = 8L, workers = 8L, jack = 100L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
source("analysis/explore_ideas47/ideas47_lib.R")   # data, truths, jackknife streams, K10 / CLR / a3r, stratifiers, score_readout, e47_gate
source("analysis/explore_stack/stack_lib.R")
say(sprintf("complex_confirm gate a, K562: G = %d, workers %d, threads %d", G, nc, opt$threads))
res_clr <- score_readout("CLR", C, macro = FALSE); res_k10 <- score_readout("K10", S9, macro = FALSE)
S1 <- a3r(Pt9); res_a3r <- score_readout("a3r", S1, macro = FALSE)
rows <- list()
g1 <- tryCatch({ e47_gate(res_clr, res_k10, res_a3r); NULL }, error = function(e) conditionMessage(e))
if (!is.null(g1)) say("e47_gate reported: ", g1)
gf <- file.path(opt$out, "equivalence_gate.csv")
if (file.exists(gf)) { g <- utils::read.csv(gf, stringsAsFactors = FALSE); for (r in seq_len(nrow(g))) rows[[length(rows) + 1L]] <- data.frame(check = g$check[r], max_abs_diff = g$max_abs_diff[r], pass = g$pass[r]) } else
  rows <- sk_gate_add(rows, "e47_gate did not write its table", Inf)
mr <- utils::read.csv("analysis/corum_exploratory/metrics_robust.csv", stringsAsFactors = FALSE); mr <- mr[mr$variant == "a_clr_bilateral_logP_t9_robust", ]
mp <- c(DoRothEA_AB = "DoRothEA_AB", CollecTRI = "CollecTRI", CORUM = "CORUM_all", KD_response = "KD_response", CORUM_noRP = "CORUM_noRP")
for (tn in names(mp)) { e <- mp[[tn]]; d <- res_a3r$raw[[e]]$aupr - res_k10$raw[[e]]$aupr; s <- stats::sd(res_a3r$raw[[e]]$jack - res_k10$raw[[e]]$jack); r <- mr[mr$truth == tn, ]
  gv <- if (tn == "CORUM_noRP") abs(d - r$d_vs_att9) else max(abs(d - r$d_vs_att9), abs(d - 1.96 * s - r$att9_lo), abs(d + 1.96 * s - r$att9_hi))
  rows <- sk_gate_add(rows, sprintf("(d) a3r minus K10, %s: %s vs metrics_robust.csv", tn, if (tn == "CORUM_noRP") "point difference only" else "d, lo, hi"), gv)
  if (tn == "CORUM") { pr <- sprintf("%+.3f [%+.3f, %+.3f]", d, d - 1.96 * s, d + 1.96 * s)
    rows <- sk_gate_add(rows, sprintf("(d) CORUM a3r minus K10 prints %s, design +0.059 [+0.033, +0.084]", pr), as.numeric(pr != "+0.059 [+0.033, +0.084]")) } }
br <- base_rate[["CORUM_all"]]
a_a3r <- res_a3r$raw$CORUM_all$aupr; a_clr <- res_clr$raw$CORUM_all$aupr; a_k10 <- res_k10$raw$CORUM_all$aupr
rows <- sk_gate_add(rows, sprintf("(design) a3r CORUM_raw AUPR %.9f vs 0.105435083", a_a3r), abs(a_a3r - 0.105435083), 5e-10)
rows <- sk_gate_add(rows, sprintf("(design) a3r CORUM_raw enrichment prints %.2f, design 7.95", a_a3r / br), abs(round(a_a3r / br, 2) - 7.95), 1e-9)
rows <- sk_gate_add(rows, sprintf("(design) CLR CORUM_raw AUPR %.4f vs 0.0847", a_clr), abs(round(a_clr, 4) - 0.0847), 1e-9)
rows <- sk_gate_add(rows, sprintf("(design) CLR CORUM_raw enrichment prints %.2f, design 6.39", a_clr / br), abs(round(a_clr / br, 2) - 6.39), 1e-9)
rows <- sk_gate_add(rows, sprintf("(design) attention CORUM_raw AUPR %.4f vs 0.0465", a_k10), abs(round(a_k10, 4) - 0.0465), 1e-9)
rows <- sk_gate_add(rows, sprintf("(design) attention CORUM_raw enrichment prints %.2f, design 3.51", a_k10 / br), abs(round(a_k10 / br, 2) - 3.51), 1e-9)
L1 <- sk_layer(C, ts_out = 9L); d_e <- max(abs(L1$outs[["9"]] - S1)); d_A <- max(abs(L1$A - A))
rows <- sk_gate_add(rows, "(e) generic stack N = 1 (operator from CLR scores, walk, a3r) vs a3r of the package path: max |S1' - S1|", d_e, 1e-12)
rows <- sk_gate_add(rows, "(e) generic stack N = 1 operator vs package operator: max |A' - A|", d_A, 1e-12)
g <- do.call(rbind, rows); g$part <- "K562"
utils::write.csv(g, file.path(opt$out, "equivalence_gate_k562.csv"), row.names = FALSE); print(g, row.names = FALSE)
if (!all(g$pass)) { say("K562 GATE FAILED"); quit(save = "no", status = 1L) }
say("K562 GATE PASSED")
