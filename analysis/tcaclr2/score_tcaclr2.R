#!/usr/bin/env Rscript
# tcaCLR2, step 3: pooled co-membership AUPR of the registered tcaCLR (t = 9, keff 10, cap 50) and tcaCLR2 (depth t2 of edge_depth.R), side by side, on one spent truth, with the
# ratio tcaCLR2 / tcaCLR and the gene delete-half jackknife of the registered scorer for these two cells only. EXPLORATORY; reading rule in docs/plans/tcaclr2_definition.md.
#   Rscript analysis/tcaclr2/score_tcaclr2.R --truth armB_cyc2008 --edge results/tcaclr2/edge_armB_cyc2008.tsv --out results/tcaclr2/score_armB_cyc2008 [--threads 8]
# --truth is a key of GT in analysis/tcaclr_depth/grid_truths.R (k562_corum, net4_cp_yeast, armA_cp_human, armB_cyc2008, m3d_tf_regulons); the matrix and the truths are exactly those of
# stage 2 of the depth grid (stage2_truth_curves.R: gl_truths, gl_ap = ap_fast with the registered tie convention), the scores are the stage-2 tcaCLR read-out (gl_walk + gl_zmat, clipped,
# Stouffer) of P^t, which check_reference.R shows bit-equal (<= 1e-12) to the registered make_scores.R / sk_a3r path. Jackknife: tc_draw_halves(list(main = pp), TC$jack_n, TC$jack_seed) and
# tc_jack_method, tc_interval, exactly as stage2_truth_curves.R does for its two-cell difference (the same halves are used for every scope). Beyond stage 2 this adds the ratio: the
# replicate ratios AP(tcaCLR2) / AP(tcaCLR) go through the same tc_interval, which is an addition of this script (flagged in the README).
# Checks written to checks.csv: for Arm A and Arm B the in-memory tcaCLR (t = 9) against the stored tcaCLR.scores of the confirmation (tolerance 1e-12); with --scores_dir the in-memory
# tcaCLR2 against tcaCLR2.scores written by make_scores_tcaclr2.R (tolerance 1e-12).
# Output (--out): results.csv (one row per scope), checks.csv, score.log. Never writes into results/tcaclr_confirm or docs/tcaclr_confirm_results (refused).
suppressPackageStartupMessages({ source("analysis/tcaclr_depth/grid_lib.R"); source("analysis/tcaclr_depth/grid_truths.R") })   # grid_truths.R sources tc_lib.R
opt <- gl_options(commandArgs(TRUE), list(truth = "", edge = "", out = "", threads = 8L, matrix = "", t2 = 0L, scores_dir = "", genes = 0L))
if (!nzchar(opt$truth) || !nzchar(opt$out)) stop("--truth and --out are required"); if (!nzchar(opt$edge) && opt$t2 < 1L) stop("--edge (or --t2 for tests) is required")
if (grepl("tcaclr_confirm", normalizePath(opt$out, mustWork = FALSE))) stop("refusing to write into the frozen confirmation outputs")
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, "score.log")); say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
t0 <- Sys.time()
if (opt$t2 > 0L) { t2 <- opt$t2; edge <- NA_real_ } else {
  e <- utils::read.delim(opt$edge, stringsAsFactors = FALSE); stopifnot(nrow(e) == 1L, e$t2 >= 1L, e$keff == TC$keff, e$cap == TC$cap, identical(e$dataset, opt$truth)); t2 <- as.integer(e$t2); edge <- e$edge }
tref <- TC$t; say(sprintf("%s: registered tcaCLR t = %d (keff %g, cap %d); tcaCLR2 t2 = %d (edge %.5f)", opt$truth, tref, TC$keff, TC$cap, t2, edge))
X <- gl_truth_matrix(opt$truth, opt$matrix); if (opt$genes > 0L) { set.seed(20261121L); X <- X[sort(sample.int(nrow(X), opt$genes)), , drop = FALSE] }
genes <- rownames(X); G <- length(genes); say(sprintf("%d genes x %d profiles", G, ncol(X)))
truths <- gl_truths(opt$truth, genes)
for (tr in truths) say(sprintf("truth %s: %d pairs, %d positive (base rate %.4f)", tr$scope, tr$n_pairs, tr$n_pos, tr$n_pos / tr$n_pairs))
f <- gl_fit(X, opt$threads); S <- f$clr_scores; f$release(); rm(X, f); invisible(gc(verbose = FALSE)); say(sprintf("MI + CLR fit: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
P <- gl_make_P(S, TC$keff, TC$cap); rm(S)
sc <- list()                                   # sc[[scope]][[as.character(t)]] = tcaCLR scores of the truth pairs (clipped, Stouffer, symmetric read-out)
gl_walk(P, sort(unique(c(tref, t2))), function(t, Pt) {
  Z <- gl_zmat(Pt)
  for (tr in truths) sc[[tr$scope]][[as.character(t)]] <<- (pmax(Z[tr$ix], 0) + pmax(Z[tr$ixT], 0)) / sqrt(2)
  rm(Z); say(sprintf("P^%d read out; peak RSS %.1f GB", t, gl_vmhwm_gb())) })
rm(P); invisible(gc(verbose = FALSE))
chk <- list(); addchk <- function(what, ref, d) { pass <- is.finite(d) && d <= 1e-12; say(sprintf("%s %s: max abs diff %s", if (pass) "PASS" else "FAIL", what, format(d)))
  chk[[length(chk) + 1L]] <<- data.frame(check = what, reference = ref, max_abs_diff = d, pass = pass, stringsAsFactors = FALSE) }
stored_dir <- c(armA_cp_human = "results/tcaclr_confirm/armA/scores", armB_cyc2008 = "results/tcaclr_confirm/armB/scores")[opt$truth]
if (!is.na(stored_dir) && opt$genes == 0L && dir.exists(stored_dir)) {
  st <- tc_scores_open(stored_dir); stopifnot(identical(st$genes, genes)); v <- tc_scores_packed(st, "tcaCLR")
  for (tr in truths) addchk(sprintf("in-memory tcaCLR (t = %d) vs stored tcaCLR.scores, truth pairs", tref), stored_dir, max(abs(sc[[tr$scope]][[as.character(tref)]] - v[tc_pair_index(tr$i, tr$j)]))); rm(v) }
if (nzchar(opt$scores_dir)) {
  st <- tc_scores_open(opt$scores_dir); stopifnot(identical(st$genes, genes)); v <- tc_scores_packed(st, "tcaCLR2")
  for (tr in truths) addchk(sprintf("in-memory tcaCLR2 (t2 = %d) vs tcaCLR2.scores of make_scores_tcaclr2.R, truth pairs", t2), opt$scores_dir, max(abs(sc[[tr$scope]][[as.character(t2)]] - v[tc_pair_index(tr$i, tr$j)]))); rm(v) }
if (length(chk)) { chk <- do.call(rbind, chk); gl_wcsv(chk, file.path(opt$out, "checks.csv")); if (!all(chk$pass)) stop("a check failed") }
res <- list()
for (tr in truths) {
  a <- sc[[tr$scope]][[as.character(t2)]]; b <- sc[[tr$scope]][[as.character(tref)]]
  pp <- list(i = tr$i, j = tr$j, label = tr$label, genes_in = tr$genes_in)
  halves <- tc_draw_halves(list(main = pp), TC$jack_n, TC$jack_seed)$main
  ra <- tc_jack_method(a, pp, halves, G, 1L); rb <- tc_jack_method(b, pp, halves, G, 1L)
  ap2 <- gl_ap(a, tr$label); ap1 <- gl_ap(b, tr$label); d <- ap2 - ap1; ratio <- ap2 / ap1
  iv_d <- tc_interval(d, ra - rb); iv_r <- tc_interval(ratio, ra / rb)
  names(iv_d) <- paste0("diff_", names(iv_d)); names(iv_r) <- paste0("ratio_", names(iv_r))
  res[[tr$scope]] <- cbind(truth = opt$truth, scope = tr$scope, G = G, n_pairs = tr$n_pairs, n_pos = tr$n_pos, edge = edge, t_tcaCLR = tref, t_tcaCLR2 = t2,
    aupr_tcaCLR = ap1, aupr_tcaCLR2 = ap2, ratio = ratio, diff = d, iv_d[, -1L], iv_r[, -1L])
}
res <- do.call(rbind, res); gl_wcsv(res, file.path(opt$out, "results.csv")); print(res[, c("scope", "t_tcaCLR2", "aupr_tcaCLR", "aupr_tcaCLR2", "ratio", "ratio_gene_lo", "ratio_gene_hi", "diff", "diff_gene_lo", "diff_gene_hi")], row.names = FALSE, digits = 4)
say(sprintf("done: %.1f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), gl_vmhwm_gb()))
