#!/usr/bin/env Rscript
# Width ladder runner (docs/plans/width_ladder_plan.md). EXPLORATORY, DESCRIPTIVE. Run from the repo root.
#   --mode gate  : reference operator (keff 10, cap 50) only: tcaCLR (t 9), tcaCLR2 (stored t2) and CLR reproduce the registered tables to 1e-9 (real truths; gate only)
#   --mode edge  : permutation edge of the reference operator and of every ladder width (one fit per permutation); t table; reference edge vs the stored edge file
#   --mode score : the ladder (13 widths x depth rules R9 / RE), multi-scale scores MS1-3, comparators; needs the edge table and a passed gate file (not for --synthetic)
#   --ds m3d|ecomac|k562|yeast|hepg2|kemmeren   [--widths 4,16,64] [--perms 10] [--workers 4] [--threads 8] [--synthetic 1] [--tag _smoke] [--out results/width_ladder]
# Smoke switches: --synthetic 1 (random planted truths, no real label scored), --widths, --perms, --tag.
suppressPackageStartupMessages(source("analysis/width_ladder/wl_lib.R"))
O <- gl_options(commandArgs(TRUE), list(mode = "", ds = "", widths = "", perms = 10L, workers = 4L, threads = 8L, synthetic = 0L, tag = "", out = WL$out))
if (!O$mode %in% c("gate", "edge", "score") || !O$ds %in% names(WLD)) stop("--mode gate|edge|score and --ds ", paste(names(WLD), collapse = "|"), " are required")
ds <- O$ds; dir <- file.path(O$out, ds); dir.create(dir, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(dir, sprintf("%s%s.log", O$mode, O$tag)))
say("options: ", paste(names(O), unlist(O), sep = "=", collapse = " ")); T0 <- Sys.time()
widths <- if (nzchar(O$widths)) as.numeric(strsplit(O$widths, ",")[[1]]) else WL$widths; stopifnot(all(widths %in% WL$widths))
f_edge <- file.path(dir, sprintf("edge%s.tsv", O$tag)); f_gate <- file.path(dir, "gate.csv"); ref_e <- wl_read_t2(ds)
tbl <- function(x, f) { tmp <- paste0(f, ".tmp", Sys.getpid()); utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, f) }

## ---------------------------------------------------------------------------------------------------------------- edge
if (O$mode == "edge") {
  X <- wl_matrix(ds); G <- nrow(X); say(sprintf("%s: %d genes x %d experiments", ds, G, ncol(X)))
  ops <- wl_ops(G, widths); r <- wl_edges(X, ops, O$perms, O$workers, O$threads, say); tab <- r$tab; tab$dataset <- ds
  wl_save(r$per, file.path(dir, sprintf("edge_perm%s.csv", O$tag))); tbl(tab, f_edge)
  ref <- tab[tab$op == "ref", ]; ok <- abs(ref$edge - ref_e$edge) <= 1e-6 && ref$t == ref_e$t2
  say(sprintf("reference operator: edge %.8f (stored %.8f), t %d (stored t2 %d): %s", ref$edge, ref_e$edge, ref$t, ref_e$t2, if (ok) "PASS" else if (O$perms == ED$perms) "FAIL" else "not comparable (perms != 10)"))
  print(tab[, c("op", "keff", "cap", "edge", "tstar", "t", "capped")], row.names = FALSE, digits = 5)
  wl_save(data.frame(dataset = ds, edge = ref$edge, edge_stored = ref_e$edge, t = ref$t, t2_stored = ref_e$t2, perms = O$perms, pass = ok), file.path(dir, sprintf("edge_gate%s.csv", O$tag)))
  say(sprintf("done: %.0f s, peak RSS of the main process %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), gl_vmhwm_gb()))
  if (!ok && O$perms == ED$perms) stop("EDGE GATE FAILED: the reference edge differs from the stored edge file")
  quit(save = "no")
}

## ---------------------------------------------------------------------------------------------------- gate / score set-up
syn <- O$synthetic == 1L; if (O$mode == "gate" && syn) stop("the gate uses real truths")
if (O$mode == "score" && !syn && !file.exists(f_gate)) stop("no passed gate file ", f_gate, ": run --mode gate first")
D <- wl_load(ds, syn); G <- D$G; U <- D$U; say(sprintf("%s: %d genes x %d experiments; truths %s", ds, G, ncol(D$X),
  paste(sprintf("%s %s (%d pairs, %d pos)", names(D$pooled), if (syn) "[SYNTHETIC]" else "", vapply(D$pooled, function(t) length(t$label), 0L), vapply(D$pooled, function(t) sum(t$label), 0)), collapse = "; ")))
tt <- if (O$mode == "score") utils::read.delim(f_edge, stringsAsFactors = FALSE) else NULL
t0 <- Sys.time(); f <- gl_fit(D$X, O$threads); C <- f$clr_scores; f$release(); D$X <- NULL; rm(f); invisible(gc(verbose = FALSE))
say(sprintf("MI + CLR fit: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
if (D$hub) { Pr <- gl_make_P(C, WL$ref$keff, WL$ref$cap); q <- rz_stationary(Pr); D$H <- outer(q, q, "+"); rm(Pr, q) }   # one stratifier per dataset (plan, fixed choice 3)
rows <- list(); ntime <- list()
addrows <- function(m, rule, score, keff, cap, t) { m$dataset <- ds; m$rule <- rule; m$score <- score; m$keff <- keff; m$cap <- cap; m$t <- t; rows[[length(rows) + 1L]] <<- m }
ev <- function(v) { t1 <- Sys.time(); m <- wl_eval(v, D); ntime[[length(ntime) + 1L]] <<- as.numeric(difftime(Sys.time(), t1, units = "secs")); m }
cells_out <- function() do.call(rbind, rows)

## --------------------------------------------------------------------------------------------------------------- gate
if (O$mode == "gate") {
  t2 <- ref_e$t2; P <- gl_make_P(C, WL$ref$keff, WL$ref$cap); vv <- list()
  gl_walk(P, sort(unique(c(9L, t2))), function(t, Pt) { v <- wl_tca_vec(Pt, U); if (t == 9L) vv$tcaCLR <<- v; if (t == t2) vv$tcaCLR2 <<- v }); rm(P)
  res <- list(tcaCLR = ev(vv$tcaCLR), tcaCLR2 = ev(vv$tcaCLR2), CLR = ev(C[U$ut]))
  SB <- utils::read.csv(TSB$sb_csv, stringsAsFactors = FALSE); R1 <- utils::read.csv(file.path(TSB$out, "cells_reference_tcaclr_all.csv"), stringsAsFactors = FALSE)
  R2 <- utils::read.csv(file.path(TSB$out, "cells_tcaclr2_all.csv"), stringsAsFactors = FALSE); VA <- if (D$hub && ds %in% c("m3d", "ecomac")) utils::read.csv(sprintf("results/direct_edge_diag/variants_%s.csv", ds), stringsAsFactors = FALSE)
  rd <- WLD[[ds]]$ref_ds; g <- list()
  look <- function(src, key, tr, sc, mt) { x <- src$enrichment_over_random[src$dataset == rd & src$method == key & src$truth == tr & src$scope == sc & src$metric == mt]; if (length(x) == 1L) x else NA_real_ }
  for (m in names(res)) { x <- res[[m]]; for (k in seq_len(nrow(x))) { tr <- x$truth[k]; sc <- x$scope[k]; mt <- x$metric[k]
    trn <- if (WLD[[ds]]$kind == "plain") "complexes_main" else tr
    st <- if (mt == "tf_aupr_ratio" && sc == "hub_att") { if (is.null(VA)) NA_real_ else { y <- VA$ratio[VA$dataset == ds & VA$truth == tr & VA$method == m & VA$variant == "d_hubcontrolled"]; if (length(y) == 1L) y else NA_real_ } }
          else look(switch(m, tcaCLR = R1, tcaCLR2 = R2, CLR = SB), switch(m, tcaCLR = "tcaCLR_t9", tcaCLR2 = "tcaCLR2", CLR = "CLR"), trn, sc, mt)
    g[[length(g) + 1L]] <- data.frame(dataset = ds, score = m, truth = tr, scope = sc, metric = mt, recomputed = x$value[k], stored = st, stringsAsFactors = FALSE) } }
  g <- do.call(rbind, g); g$abs_diff <- abs(g$recomputed - g$stored); g$compared <- is.finite(g$stored); g$pass <- !g$compared | g$abs_diff <= 1e-9
  if (!is.null(VA)) { ge <- utils::read.csv("results/direct_edge_diag/gate_direct_edge.csv", stringsAsFactors = FALSE); ge <- ge[ge$dataset == ds & ge$method %in% names(res), ]   # direct-edge raw numbers of the diagnostic
    for (k in seq_len(nrow(ge))) { x <- res[[ge$method[k]]]; v <- x$value[x$truth == ge$truth[k] & x$scope == "raw" & x$metric == ge$metric[k]]; stopifnot(length(v) == 1L)
      g <- rbind(g, data.frame(dataset = ds, score = ge$method[k], truth = ge$truth[k], scope = "raw", metric = paste0(ge$metric[k], " (gate_direct_edge.csv)"), recomputed = v, stored = ge$stored[k], abs_diff = abs(v - ge$stored[k]), compared = TRUE, pass = abs(v - ge$stored[k]) <= 1e-9)) } }
  say(sprintf("gate: %d cells recomputed, %d compared, max abs diff %.3g, %d failing", nrow(g), sum(g$compared), max(g$abs_diff[g$compared]), sum(!g$pass)))
  print(g[g$compared & g$score == "tcaCLR", c("score", "truth", "scope", "metric", "recomputed", "stored", "abs_diff")], row.names = FALSE, digits = 7)
  if (any(!g$pass)) print(g[!g$pass, ], row.names = FALSE)
  if (all(g$pass) && sum(g$compared) > 0L) { wl_save(g, f_gate); say("GATE PASSED -> ", f_gate) } else { wl_save(g, file.path(dir, "gate_FAILED.csv")); say("GATE FAILED"); quit(status = 1) }
  say(sprintf("done: %.0f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), gl_vmhwm_gb())); quit(save = "no")
}

## -------------------------------------------------------------------------------------------------------------- score
N <- G * (G - 1) / 2; RULES <- c("R9", "RE"); ACC <- lapply(stats::setNames(RULES, RULES), function(r) list(ms1 = NULL, ms2 = NULL, ms3 = NULL, n = 0L, n3 = 0L))
accumulate <- function(rule, keff, v) {
  a <- ACC[[rule]]; r <- rank(v, ties.method = "average") / N; z <- (v - mean(v)) / stats::sd(v)
  a$ms1 <- if (is.null(a$ms1)) r else a$ms1 + r; a$ms2 <- if (is.null(a$ms2)) z else a$ms2 + z; a$n <- a$n + 1L
  if (keff %in% WL$ms3) { a$ms3 <- if (is.null(a$ms3)) r else a$ms3 + r; a$n3 <- a$n3 + 1L }
  ACC[[rule]] <<- a }
addrows(ev(C[U$ut]), "both", "CLR", NA, NA, NA); say("scored CLR")
ops <- wl_ops(G, widths); tw_ref <- ref_e$t2
for (o in ops) {
  t1 <- Sys.time(); tw <- if (o$ref) tw_ref else tt$t[tt$op == o$id]; stopifnot(length(tw) == 1L)
  P <- gl_make_P(C, o$keff, o$cap)
  gl_walk(P, sort(unique(c(WL$t_fixed, tw))), function(t, Pt) {
    v <- wl_tca_vec(Pt, U); rl <- c(if (t == WL$t_fixed) "R9", if (t == tw) "RE"); m <- ev(v)
    for (r in rl) { if (o$ref) addrows(m, r, if (r == "R9") "tcaCLR" else "tcaCLR2", o$keff, o$cap, t) else { addrows(m, r, o$id, o$keff, o$cap, t); accumulate(r, o$keff, v) } }
    rm(v) })
  rm(P); invisible(gc(verbose = FALSE)); wl_save(cells_out(), file.path(dir, sprintf("cells%s.csv", O$tag)))
  say(sprintf("operator %s keff %g cap %d t_RE %d: %.0f s, peak RSS %.1f GB", o$id, o$keff, o$cap, tw, as.numeric(difftime(Sys.time(), t1, units = "secs")), gl_vmhwm_gb()))
}
rm(C); invisible(gc(verbose = FALSE))
for (r in RULES) { a <- ACC[[r]]
  for (s in c("MS1", "MS2", "MS3")) { vec <- switch(s, MS1 = a$ms1 / a$n, MS2 = a$ms2 / a$n, MS3 = if (a$n3 > 0L) a$ms3 / a$n3 else NULL); if (is.null(vec)) next
    addrows(ev(vec), r, s, NA, NA, NA); say(sprintf("scored %s %s (%d widths)", r, s, if (s == "MS3") a$n3 else a$n)); rm(vec) } }
cells <- cells_out(); wl_save(cells, file.path(dir, sprintf("cells%s.csv", O$tag)))
wl_save(data.frame(dataset = ds, G = G, widths = length(widths), synthetic = syn, wall_s = as.numeric(difftime(Sys.time(), T0, units = "secs")), eval_s_mean = mean(unlist(ntime)), n_evals = length(ntime),
  peak_rss_gb = gl_vmhwm_gb(), threads = O$threads), file.path(dir, sprintf("run%s.csv", O$tag)))
say(sprintf("done: %.0f s, %d cells, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), nrow(cells), gl_vmhwm_gb()))
