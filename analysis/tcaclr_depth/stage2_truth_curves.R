#!/usr/bin/env Rscript
# Stage 2 of docs/plans/tcaclr_depth_grid_plan.md (Revision 2): pooled co-membership AUPR (no jackknife; ap_fast, the registered tie convention) for every
# cell of the grid (aCLR and tcaCLR, both clipS, both clipO) on one spent truth. DESCRIPTIVE; the truth-chosen best cell is an upper bound, never a claim.
#   Rscript analysis/tcaclr_depth/stage2_truth_curves.R --truth armB_cyc2008 --out results/tcaclr_depth/stage2_armB_cyc2008 [--matrix FILE|net4] [--workers 8]
#     truths: k562_corum, net4_cp_yeast, armA_cp_human, armB_cyc2008, m3d_tf_regulons (M3D: scopes SC and all)
# Optional: --lf_cell "keff,cap,clipS,clipO,t" (a tcaCLR cell, e.g. the stage 1 label-free cell) adds the gene delete-half jackknife (500 replicates, seed 20261117,
#   tc_draw_halves / tc_jack_method of the registered scorer) for the difference between that cell and the registered cell, for those two cells only.
# Output: curves.csv (truth, scope, method, keff, cap, clipS, clipO, t, aupr, n_pairs, n_pos), [jackknife.csv], stage2.log.
suppressPackageStartupMessages({ source("analysis/tcaclr_depth/grid_lib.R"); source("analysis/tcaclr_depth/grid_truths.R") })
opt <- gl_options(commandArgs(TRUE), list(truth = "", matrix = "", out = "", workers = 4L, threads = 8L, keffs = "", depths = "", reference = 1L, lf_cell = ""))
if (!nzchar(opt$truth) || !nzchar(opt$out)) stop("--truth and --out are required")
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, "stage2.log")); say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
X <- gl_truth_matrix(opt$truth, opt$matrix); genes <- rownames(X); G <- length(genes); say(sprintf("%s: %d genes x %d profiles", opt$truth, G, ncol(X)))
truths <- gl_truths(opt$truth, genes)
for (tr in truths) say(sprintf("truth %s: %d pairs, %d positive (base rate %.4f)", tr$scope, tr$n_pairs, tr$n_pos, tr$n_pos / tr$n_pairs))
keffs <- if (nzchar(opt$keffs)) as.numeric(strsplit(opt$keffs, ",")[[1]]) else GL$keffs
depths <- if (nzchar(opt$depths)) as.integer(strsplit(opt$depths, ",")[[1]]) else GL$depths
ops <- gl_grid_ops(reference = opt$reference == 1L, keffs = keffs, depths = depths)
lf <- if (nzchar(opt$lf_cell)) { v <- strsplit(opt$lf_cell, ",")[[1]]; list(keff = as.numeric(v[1]), cap = as.integer(v[2]), clipS = as.logical(v[3]), clipO = as.logical(v[4]), t = as.integer(v[5])) } else NULL
cell_id <- function(m, keff, cap, clipS, clipO, t) paste(m, keff, cap, clipS, clipO, t, sep = "|")
keep_ids <- c(cell_id("tcaCLR", GL$ref$keff, GL$ref$cap, TRUE, TRUE, GL$ref$t), if (!is.null(lf)) cell_id("tcaCLR", lf$keff, lf$cap, lf$clipS, lf$clipO, lf$t))
t0 <- Sys.time(); f <- gl_fit(X, opt$threads); say(sprintf("MI + CLR fit: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
Ssc <- list(`TRUE` = f$clr_scores, `FALSE` = gl_unclipped_clr(f$mi)); f$release(); rm(X); invisible(gc(verbose = FALSE))
one_operator <- function(o) {
  t1 <- Sys.time(); rows <- list(); kept <- list()
  P <- gl_make_P(Ssc[[as.character(o$clipS)]], o$keff, o$cap)
  gl_walk(P, o$depths, function(t, Pt) {
    Z <- gl_zmat(Pt)
    for (tr in truths) {
      s <- list(aCLR = (Pt[tr$ix] + Pt[tr$ixT]) / 2); z1 <- Z[tr$ix]; z2 <- Z[tr$ixT]
      s$tca_clip <- (pmax(z1, 0) + pmax(z2, 0)) / sqrt(2); s$tca_unclip <- (z1 + z2) / sqrt(2)
      for (vn in names(s)) {
        rows[[length(rows) + 1L]] <<- data.frame(truth = opt$truth, scope = tr$scope, method = gl_variant_names[[vn]], keff = o$keff, cap = o$cap, clipS = o$clipS,
          clipO = gl_variant_clipO[[vn]], t = t, aupr = gl_ap(s[[vn]], tr$label), n_pairs = tr$n_pairs, n_pos = tr$n_pos, stringsAsFactors = FALSE)
        id <- cell_id(gl_variant_names[[vn]], o$keff, o$cap, o$clipS, gl_variant_clipO[[vn]], t)
        if (id %in% keep_ids) kept[[paste(tr$scope, id)]] <<- s[[vn]] } }
    rm(Z) })
  msg <- sprintf("operator keff %g cap %d clipS %s: %.0f s, peak RSS %.1f GB", o$keff, o$cap, o$clipS, as.numeric(difftime(Sys.time(), t1, units = "secs")), gl_vmhwm_gb())
  list(rows = do.call(rbind, rows), kept = kept, msg = msg)
}
res <- if (opt$workers > 1L) parallel::mclapply(ops, function(o) tryCatch(one_operator(o), error = function(e) list(err = conditionMessage(e))), mc.cores = opt$workers, mc.preschedule = FALSE) else lapply(ops, one_operator)
bad <- vapply(res, function(r) !is.null(r$err) || is.null(r$rows), NA); if (any(bad)) stop("operator failed: ", paste(vapply(res[bad], function(r) as.character(r$err)[1], ""), collapse = " | "))
for (r in res) say(r$msg)
curves <- do.call(rbind, lapply(res, `[[`, "rows")); gl_wcsv(curves, file.path(opt$out, "curves.csv")); say("wrote curves.csv: ", nrow(curves), " cells")
if (!is.null(lf)) {   # gene delete-half jackknife of (label-free cell) - (registered cell), for those two cells only
  kept <- do.call(c, lapply(res, `[[`, "kept")); jr <- list()
  for (tr in truths) { a <- kept[[paste(tr$scope, keep_ids[2])]]; b <- kept[[paste(tr$scope, keep_ids[1])]]; stopifnot(!is.null(a), !is.null(b))
    pp <- list(i = tr$i, j = tr$j, label = tr$label, genes_in = tr$genes_in); halves <- tc_draw_halves(list(main = pp), TC$jack_n, TC$jack_seed)$main
    ra <- tc_jack_method(a, pp, halves, G, 1L); rb <- tc_jack_method(b, pp, halves, G, 1L); d <- gl_ap(a, tr$label) - gl_ap(b, tr$label)
    jr[[tr$scope]] <- cbind(truth = opt$truth, scope = tr$scope, lf_cell = opt$lf_cell, aupr_lf = gl_ap(a, tr$label), aupr_registered = gl_ap(b, tr$label), tc_interval(d, ra - rb)) }
  jr <- do.call(rbind, jr); gl_wcsv(jr, file.path(opt$out, "jackknife.csv")); print(jr, row.names = FALSE)
}
say("done")
