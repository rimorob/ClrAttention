#!/usr/bin/env Rscript
# Global-temperature ladder runner (docs/plans/global_temperature_ladder_plan.md). EXPLORATORY, DESCRIPTIVE; spent truths. Run from the repo root.
#   --mode gate  : reference operator (keff 10, cap 50) and V2's own temperature at t = 9 reproduce the stored tables to 1e-9 (real truths; gate only; writes gate.csv)
#   --mode edge  : real-data temperature table (V2's tau, effective partners per tau), permutation edge of the reference operator, V2's tau and every ladder tau; t table
#   --mode score : gate (inline, stops the run on failure), then the ladder (13 taus + V2 reference x depth rules R9 / RE), multi-scale scores MS1-3, comparators
#   --ds m3d|ecomac|k562|yeast|hepg2|kemmeren   [--log2taus -3,0,2] [--perms 10] [--workers 4] [--threads 8] [--synthetic 1] [--tag _smoke] [--out results/global_temp_ladder]
# Smoke switches: --synthetic 1 (random planted truths, no real label scored; the gate is skipped), --log2taus (subset of the ladder), --perms, --tag.
suppressPackageStartupMessages(source("analysis/global_temp_ladder/gt_lib.R"))
O <- gl_options(commandArgs(TRUE), list(mode = "", ds = "", log2taus = "", perms = 10L, workers = 4L, threads = 8L, synthetic = 0L, tag = "", out = GTL$out))
if (!O$mode %in% c("gate", "edge", "score") || !O$ds %in% names(WLD)) stop("--mode gate|edge|score and --ds ", paste(names(WLD), collapse = "|"), " are required")
ds <- O$ds; dir <- file.path(O$out, ds); dir.create(dir, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(dir, sprintf("%s%s.log", O$mode, O$tag)))
say("options: ", paste(names(O), unlist(O), sep = "=", collapse = " ")); T0 <- Sys.time()
if (nzchar(O$log2taus)) { l2 <- as.numeric(strsplit(O$log2taus, ",")[[1]]); stopifnot(all(round(l2, 3) %in% round(log2(GTL$taus), 3))); GTL$taus <- 2^l2 }
f_edge <- file.path(dir, sprintf("edge%s.tsv", O$tag)); f_tau <- file.path(dir, sprintf("tau_table%s.csv", O$tag)); f_gate <- file.path(dir, "gate.csv"); ref_e <- wl_read_t2(ds)
tbl <- function(x, f) { tmp <- paste0(f, ".tmp", Sys.getpid()); utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, f) }
tau_table <- function(Tm, tau2) do.call(rbind, c(list(gt_eff_row(Tm, tau2, "V2", ds)), lapply(GTL$taus, function(tau) gt_eff_row(Tm, tau, gt_tau_id(tau), ds))))

## ---------------------------------------------------------------------------------------------------------------- edge
if (O$mode == "edge") {
  X <- wl_matrix(ds); G <- nrow(X); say(sprintf("%s: %d genes x %d experiments", ds, G, ncol(X)))
  # phase 1 (no truth): real-data fit in a forked child (the parent must not fit before forking), V2's tau and the effective-partner table of every tau
  # (parallel::mcparallel always forks; ad_lapply / mclapply with a single task runs in the parent, after which the later forked workers deadlock)
  ph1 <- parallel::mccollect(parallel::mcparallel(tryCatch({ f <- gl_fit(X, O$threads); C <- f$clr_scores; f$release(); rm(f); Tm <- gt_topmat(C); tau2 <- gt_tau_for_eff(Tm)
    list(tau2 = tau2, tab = tau_table(Tm, tau2), rss = gl_vmhwm_gb()) }, error = function(e) paste("ERROR", conditionMessage(e)))))[[1]]
  if (!is.list(ph1)) stop("real-data phase failed: ", ph1)
  tau2 <- ph1$tau2; wl_save(ph1$tab, f_tau); say(sprintf("V2 temperature of %s: tau = %.8f (log2 %.4f); real-data fit child peak RSS %.1f GB", ds, tau2, log2(tau2), ph1$rss)); print(ph1$tab[, c("op", "tau", "eff_mean", "eff_median", "eff_q10", "eff_q90")], row.names = FALSE, digits = 4)
  ops <- gt_ops(tau2); r <- gt_edges(X, ops, O$perms, O$workers, O$threads, say); tab <- r$tab; tab$dataset <- ds
  wl_save(r$per, file.path(dir, sprintf("edge_perm%s.csv", O$tag))); tbl(tab, f_edge)
  ref <- tab[tab$op == "ref", ]; ok <- abs(ref$edge - ref_e$edge) <= 1e-6 && ref$t == ref_e$t2
  say(sprintf("reference operator: edge %.8f (stored %.8f), t %d (stored t2 %d): %s", ref$edge, ref_e$edge, ref$t, ref_e$t2, if (ok) "PASS" else if (O$perms == ED$perms) "FAIL" else "not comparable (perms != 10)"))
  print(tab[, c("op", "tau", "edge", "tstar", "t", "capped")], row.names = FALSE, digits = 5)
  wl_save(data.frame(dataset = ds, edge = ref$edge, edge_stored = ref_e$edge, t = ref$t, t2_stored = ref_e$t2, perms = O$perms, pass = ok), file.path(dir, sprintf("edge_gate%s.csv", O$tag)))
  say(sprintf("done: %.0f s, peak RSS of the main process %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), gl_vmhwm_gb()))
  if (!ok && O$perms == ED$perms) stop("EDGE GATE FAILED: the reference edge differs from the stored edge file")
  quit(save = "no")
}

## ---------------------------------------------------------------------------------------------------- gate / score set-up
syn <- O$synthetic == 1L; if (O$mode == "gate" && syn) stop("the gate uses real truths")
tt <- if (O$mode == "score") utils::read.delim(f_edge, stringsAsFactors = FALSE) else NULL
D <- wl_load(ds, syn); G <- D$G; U <- D$U; say(sprintf("%s: %d genes x %d experiments; truths %s", ds, G, ncol(D$X),
  paste(sprintf("%s %s (%d pairs, %d pos)", names(D$pooled), if (syn) "[SYNTHETIC]" else "", vapply(D$pooled, function(t) length(t$label), 0L), vapply(D$pooled, function(t) sum(t$label), 0)), collapse = "; ")))
t0 <- Sys.time(); f <- gl_fit(D$X, O$threads); C <- f$clr_scores; f$release(); D$X <- NULL; rm(f); invisible(gc(verbose = FALSE))
say(sprintf("MI + CLR fit: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
Tm <- gt_topmat(C); tau2 <- gt_tau_for_eff(Tm); TT <- tau_table(Tm, tau2); wl_save(TT, if (O$mode == "gate") file.path(dir, sprintf("tau_table_gate%s.csv", O$tag)) else file.path(dir, sprintf("tau_table_score%s.csv", O$tag)))
say(sprintf("V2 temperature: tau = %.8f (log2 %.4f)", tau2, log2(tau2)))
if (O$mode == "score" && !syn) {   # the score job must see the same V2 temperature as the edge job (the permuted operators were built at it)
  te <- utils::read.csv(f_tau, stringsAsFactors = FALSE); if (abs(te$tau[te$op == "V2"] / tau2 - 1) > 1e-9) stop("V2 tau differs from the edge job's: ", te$tau[te$op == "V2"], " vs ", tau2)
  if (!all(abs(te$eff_mean - TT$eff_mean[match(te$op, TT$op)]) < 1e-9)) stop("effective-partner table differs from the edge job's") }
if (D$hub) { Pr <- gl_make_P(C, GTL$ref$keff, GTL$ref$cap); q <- rz_stationary(Pr); D$H <- outer(q, q, "+"); rm(Pr, q) }   # one stratifier per dataset (the registered operator's stationary mass)
rows <- list(); ntime <- list()
addrows <- function(m, rule, score, tau, t) { m$dataset <- ds; m$rule <- rule; m$score <- score; m$tau <- tau; m$t <- t; rows[[length(rows) + 1L]] <<- m }
ev <- function(v) { t1 <- Sys.time(); m <- wl_eval(v, D); ntime[[length(ntime) + 1L]] <<- as.numeric(difftime(Sys.time(), t1, units = "secs")); m }
cells_out <- function() do.call(rbind, rows)

## ------------------------------------------------------------------------------------------------ gate (reference + V2 point)
# gate 1: the registered operator through this code reproduces the stored tcaCLR / tcaCLR2 / CLR cells (the width-ladder gate comparison);
# gate 2: V2's own temperature at t = 9 reproduces results/recipe_ladder/variants_{m3d,ecomac}.csv (V2_global_temp) and the stored tau (variant_info_*.csv). Tolerance 1e-9.
run_gate <- function(res, resV2) {
  SB <- utils::read.csv(TSB$sb_csv, stringsAsFactors = FALSE); R1 <- utils::read.csv(file.path(TSB$out, "cells_reference_tcaclr_all.csv"), stringsAsFactors = FALSE)
  R2 <- utils::read.csv(file.path(TSB$out, "cells_tcaclr2_all.csv"), stringsAsFactors = FALSE); VA <- if (D$hub && ds %in% c("m3d", "ecomac")) utils::read.csv(sprintf("results/direct_edge_diag/variants_%s.csv", ds), stringsAsFactors = FALSE)
  rd <- WLD[[ds]]$ref_ds; g <- list()
  look <- function(src, key, tr, sc, mt) { x <- src$enrichment_over_random[src$dataset == rd & src$method == key & src$truth == tr & src$scope == sc & src$metric == mt]; if (length(x) == 1L) x else NA_real_ }
  for (m in names(res)) { x <- res[[m]]; for (k in seq_len(nrow(x))) { tr <- x$truth[k]; sc <- x$scope[k]; mt <- x$metric[k]
    trn <- if (WLD[[ds]]$kind == "plain") "complexes_main" else tr
    st <- if (mt == "tf_aupr_ratio" && sc == "hub_att") { if (is.null(VA)) NA_real_ else { y <- VA$ratio[VA$dataset == ds & VA$truth == tr & VA$method == m & VA$variant == "d_hubcontrolled"]; if (length(y) == 1L) y else NA_real_ } }
          else look(switch(m, tcaCLR = R1, tcaCLR2 = R2, CLR = SB), switch(m, tcaCLR = "tcaCLR_t9", tcaCLR2 = "tcaCLR2", CLR = "CLR"), trn, sc, mt)
    g[[length(g) + 1L]] <- data.frame(gate = 1L, dataset = ds, score = m, truth = tr, scope = sc, metric = mt, recomputed = x$value[k], stored = st, stringsAsFactors = FALSE) } }
  if (!is.null(VA)) { ge <- utils::read.csv("results/direct_edge_diag/gate_direct_edge.csv", stringsAsFactors = FALSE); ge <- ge[ge$dataset == ds & ge$method %in% names(res), ]
    for (k in seq_len(nrow(ge))) { x <- res[[ge$method[k]]]; v <- x$value[x$truth == ge$truth[k] & x$scope == "raw" & x$metric == ge$metric[k]]; stopifnot(length(v) == 1L)
      g[[length(g) + 1L]] <- data.frame(gate = 1L, dataset = ds, score = ge$method[k], truth = ge$truth[k], scope = "raw", metric = paste0(ge$metric[k], " (gate_direct_edge.csv)"), recomputed = v, stored = ge$stored[k], stringsAsFactors = FALSE) } }
  if (ds %in% c("m3d", "ecomac")) {   # gate 2: V2 (recipe ladder variants file, stage V2_global_temp) and its temperature
    V <- utils::read.csv(sprintf("results/recipe_ladder/variants_%s.csv", ds), stringsAsFactors = FALSE); V <- V[V$stage == "V2_global_temp", ]; vi <- utils::read.csv(sprintf("results/recipe_ladder/variant_info_%s.csv", ds), stringsAsFactors = FALSE)
    for (tr in V$truth) { r <- V[V$truth == tr, ]; get <- function(sc, mt) resV2$value[resV2$truth == tr & resV2$scope == sc & resV2$metric == mt]
      for (z in list(list("raw", "tf_aupr_ratio", r$ratio), list("raw", "tf_auroc_median_tf_rows", r$median_auroc), list("hub_att", "tf_aupr_ratio", r$hubctrl_ratio)))
        g[[length(g) + 1L]] <- data.frame(gate = 2L, dataset = ds, score = "V2", truth = tr, scope = z[[1]], metric = z[[2]], recomputed = get(z[[1]], z[[2]]), stored = z[[3]], stringsAsFactors = FALSE) }
    g[[length(g) + 1L]] <- data.frame(gate = 2L, dataset = ds, score = "V2", truth = "tau", scope = "-", metric = "tau_global (variant_info)", recomputed = tau2, stored = vi$tau_global, stringsAsFactors = FALSE) }
  g <- do.call(rbind, g); g$abs_diff <- abs(g$recomputed - g$stored); g$compared <- is.finite(g$stored); g$pass <- !g$compared | g$abs_diff <= 1e-9
  g$printed_equal <- ifelse(g$compared, signif(g$recomputed, 3) == signif(g$stored, 3), NA)
  say(sprintf("gate 1: %d cells compared, max abs diff %.3g, %d failing", sum(g$compared & g$gate == 1L), max(g$abs_diff[g$compared & g$gate == 1L]), sum(!g$pass & g$gate == 1L)))
  if (any(g$gate == 2L)) { say(sprintf("gate 2 (V2): %d cells compared, max abs diff %.3g, %d failing", sum(g$compared & g$gate == 2L), max(g$abs_diff[g$compared & g$gate == 2L]), sum(!g$pass & g$gate == 2L)))
    print(g[g$gate == 2L, c("truth", "scope", "metric", "recomputed", "stored", "abs_diff")], row.names = FALSE, digits = 8) } else say("gate 2: no stored V2 value for this dataset (recipe ladder ran on M3D and EcoMAC only)")
  print(g[g$compared & g$score == "tcaCLR" & g$gate == 1L, c("score", "truth", "scope", "metric", "recomputed", "stored", "abs_diff")], row.names = FALSE, digits = 7)
  if (all(g$pass) && sum(g$compared & g$gate == 1L) > 0L) { wl_save(g, f_gate); say("GATE PASSED -> ", f_gate) } else { print(g[!g$pass, ], row.names = FALSE); wl_save(g, file.path(dir, "gate_FAILED.csv")); say("GATE FAILED"); quit(status = 1) }
}

## -------------------------------------------------------------------------------------------------------------- score / gate loop
N <- G * (G - 1) / 2; RULES <- c("R9", "RE"); ACC <- lapply(stats::setNames(RULES, RULES), function(r) list(ms1 = NULL, ms2 = NULL, ms3 = NULL, n = 0L, n3 = 0L))
accumulate <- function(rule, tau, v) {   # multi-scale accumulation over the LADDER points only (not the reference V2 point, not the registered reference operator)
  a <- ACC[[rule]]; r <- rank(v, ties.method = "average") / N; z <- (v - mean(v)) / stats::sd(v)
  a$ms1 <- if (is.null(a$ms1)) r else a$ms1 + r; a$ms2 <- if (is.null(a$ms2)) z else a$ms2 + z; a$n <- a$n + 1L
  if (round(log2(tau), 3) %in% round(GTL$ms3_log2, 3)) { a$ms3 <- if (is.null(a$ms3)) r else a$ms3 + r; a$n3 <- a$n3 + 1L }
  ACC[[rule]] <<- a }
ops <- gt_ops(tau2); gate_only <- O$mode == "gate"; resG <- list(); resV2 <- NULL
addrows(ev(C[U$ut]), "both", "CLR", NA, NA); resG$CLR <- rows[[length(rows)]][, c("truth", "scope", "metric", "value")]; say("scored CLR")
for (o in ops) {
  t1 <- Sys.time(); if (gate_only && !o$id %in% c("ref", "V2")) next
  tw <- if (o$id == "ref") ref_e$t2 else if (gate_only) NA_integer_ else tt$t[tt$op == o$id]; stopifnot(length(tw) == 1L)
  P <- gt_P(o, C, Tm); depths <- sort(unique(c(GTL$t_fixed, if (!is.na(tw)) tw)))
  gl_walk(P, depths, function(t, Pt) {
    v <- wl_tca_vec(Pt, U); rl <- c(if (t == GTL$t_fixed) "R9", if (!is.na(tw) && t == tw) "RE"); m <- ev(v)
    if (o$id == "ref") { if (t == GTL$t_fixed) resG$tcaCLR <<- m; if (t == tw) resG$tcaCLR2 <<- m }
    if (o$id == "V2" && t == GTL$t_fixed) { resV2 <<- m
      if (!syn && ds %in% c("m3d", "ecomac")) { source("analysis/direct_edge_diag/lib.R"); S <- wl_unpack(v, D)   # the recipe ladder's per-TF AUROC is tf_rows (TFs with >= 5 targets), not the scoreboard's sb_per_tf_median_auroc (>= 3): gate 2 uses the recipe's own function
        resV2 <<- rbind(m, do.call(rbind, lapply(names(D$edge), function(nm) data.frame(truth = nm, scope = "raw", metric = "tf_auroc_median_tf_rows", value = stats::median(tf_rows(S, D$edge[[nm]], G)$auroc), stringsAsFactors = FALSE)))); rm(S) } }
    for (r in rl) { if (o$id == "ref") addrows(m, r, if (r == "R9") "tcaCLR" else "tcaCLR2", NA, t) else { addrows(m, r, o$id, o$tau, t); if (o$id != "V2") accumulate(r, o$tau, v) } }
    rm(v) })
  rm(P); invisible(gc(verbose = FALSE))
  say(sprintf("operator %s tau %s t_RE %s: %.0f s, peak RSS %.1f GB", o$id, if (is.na(o$tau)) "per-gene" else sprintf("%.5f", o$tau), if (is.na(tw)) "-" else tw, as.numeric(difftime(Sys.time(), t1, units = "secs")), gl_vmhwm_gb()))
  if (o$id == "V2") { if (!syn) run_gate(resG[c("tcaCLR", "tcaCLR2", "CLR")[c("tcaCLR", "tcaCLR2", "CLR") %in% names(resG)]], resV2); if (gate_only) { say(sprintf("done: %.0f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), gl_vmhwm_gb())); quit(save = "no") } }
  wl_save(cells_out(), file.path(dir, sprintf("cells%s%s.csv", if (gate_only) "_gate" else "", O$tag)))
}
rm(C); invisible(gc(verbose = FALSE))
for (r in RULES) { a <- ACC[[r]]
  for (s in c("MS1", "MS2", "MS3")) { vec <- switch(s, MS1 = a$ms1 / a$n, MS2 = a$ms2 / a$n, MS3 = if (a$n3 > 0L) a$ms3 / a$n3 else NULL); if (is.null(vec)) next
    addrows(ev(vec), r, s, NA, NA); say(sprintf("scored %s %s (%d temperatures)", r, s, if (s == "MS3") a$n3 else a$n)); rm(vec) } }
cells <- cells_out(); wl_save(cells, file.path(dir, sprintf("cells%s.csv", O$tag)))
wl_save(data.frame(dataset = ds, G = G, temperatures = length(GTL$taus), synthetic = syn, wall_s = as.numeric(difftime(Sys.time(), T0, units = "secs")), eval_s_mean = mean(unlist(ntime)), n_evals = length(ntime),
  peak_rss_gb = gl_vmhwm_gb(), threads = O$threads), file.path(dir, sprintf("run%s.csv", O$tag)))
say(sprintf("done: %.0f s, %d cells, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), nrow(cells), gl_vmhwm_gb()))
