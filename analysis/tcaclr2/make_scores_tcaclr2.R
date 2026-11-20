#!/usr/bin/env Rscript
# tcaCLR2, step 2: the tcaCLR2 score file of one matrix, in exactly the format of analysis/tcaclr_confirm/make_scores.R (genes.txt + packed column-major
# upper-triangle score file, here tcaCLR2.scores). tcaCLR2 = the registered tcaCLR (keff 10, cap 50, lazy 0.5, sk_a3r read-out, clipped, Stouffer, symmetric) at the depth t2 of
# analysis/tcaclr2/edge_depth.R (docs/plans/tcaclr2_definition.md). No truth is read.
#   Rscript analysis/tcaclr2/make_scores_tcaclr2.R --dataset armB_cyc2008 --edge results/tcaclr2/edge_armB_cyc2008.tsv --out results/tcaclr2/scores_armB_cyc2008 [--threads 8]
#   Rscript analysis/tcaclr2/make_scores_tcaclr2.R --expr FILE [--sha FILE] --edge TSV --out DIR
# Self-test (reference check): --t 9 --name tcaCLR_t9 --compare results/tcaclr_confirm/armB/scores   writes tcaCLR_t9.scores to --out and reports the maximum absolute difference to the
#   stored tcaCLR.scores of --compare (the registered Arm B scores; expected 0, tolerance 1e-12) and to the stored genes.txt. With --t the edge file is not read.
# Never writes into results/tcaclr_confirm (refused). The depth chain is the registered one (t - 1 sequential products of the sparse P with the dense power, as make_scores.R).
suppressPackageStartupMessages({ source("analysis/tcaclr_depth/grid_lib.R"); source("analysis/tcaclr_depth/grid_truths.R") })   # grid_truths.R sources tc_lib.R
opt <- tc_options(commandArgs(TRUE), list(dataset = "", expr = "", sha = "", edge = "", t = 0L, name = "tcaCLR2", compare = "", out = "", threads = 8L))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
if (!nzchar(opt$out)) stop("--out is required"); if (grepl("tcaclr_confirm", normalizePath(opt$out, mustWork = FALSE))) stop("refusing to write into the frozen confirmation outputs")
if (!nzchar(opt$dataset) && !nzchar(opt$expr)) stop("--dataset or --expr is required")
t0 <- Sys.time()
if (opt$t > 0L) t2 <- as.integer(opt$t) else {
  if (!nzchar(opt$edge)) stop("--edge (the edge_<dataset>.tsv of edge_depth.R) or --t is required")
  e <- utils::read.delim(opt$edge, stringsAsFactors = FALSE); stopifnot(nrow(e) == 1L, e$t2 >= 1L, is.finite(e$edge), e$keff == TC$keff, e$cap == TC$cap)
  if (nzchar(opt$dataset)) stopifnot(identical(e$dataset, opt$dataset)); t2 <- as.integer(e$t2); say(sprintf("t2 = %d from %s (edge %.5f)", t2, opt$edge, e$edge)) }
if (nzchar(opt$expr)) {
  if (nzchar(opt$sha)) { hs <- tolower(vapply(strsplit(trimws(readLines(opt$sha, warn = FALSE)), "[[:space:]]+"), `[`, "", 1L)); got <- tc_sha256(opt$expr)
    say("expression SHA-256 ", got, if (got %in% hs) " matches the freeze" else " DOES NOT MATCH the freeze"); if (!got %in% hs) stop("the expression file does not match its frozen hash") }
  X <- gl_read_matrix(opt$expr)
} else X <- gl_truth_matrix(opt$dataset)
genes <- rownames(X); G <- length(genes); stopifnot(!anyDuplicated(genes), all(is.finite(X)))
say(sprintf("%s: %d genes x %d profiles; depth %d", if (nzchar(opt$dataset)) opt$dataset else opt$expr, G, ncol(X), t2))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
gf <- file.path(opt$out, "genes.txt"); if (file.exists(gf)) stopifnot(identical(readLines(gf), genes)) else writeLines(genes, gf)
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$calibrate(method = "normal", combine = "stouffer")
f$build_operator(alpha = TC$alpha, softmax_keff = TC$keff, softmax_cap = as.integer(TC$cap))
P <- Matrix::Diagonal(G, 1 - TC$alpha) + TC$alpha * Matrix::Matrix(f$operator, sparse = TRUE); f$release()
Pt <- gl_power_seq(P, t2); rm(P); invisible(gc(verbose = FALSE))
S <- tc_pack(tc_clean_sym(sk_a3r(Pt))); rm(Pt); invisible(gc(verbose = FALSE))
invisible(tc_write_bin(S, file.path(opt$out, paste0(opt$name, ".scores")))); say("wrote ", opt$name, ".scores (", length(S), " pairs)")
if (nzchar(opt$compare)) {
  st <- tc_scores_open(opt$compare); ok <- identical(st$genes, genes); say(sprintf("%s gene order identical to genes.txt of %s", if (ok) "PASS" else "FAIL", opt$compare))
  ref <- tc_scores_packed(st, "tcaCLR"); d <- max(abs(S - ref)); pass <- ok && d <= 1e-12
  say(sprintf("%s: %s (t = %d) vs stored tcaCLR.scores of %s: max abs difference %s over %d pairs (expected 0; tolerance 1e-12)", if (pass) "PASS" else "FAIL", opt$name, t2, opt$compare, format(d), length(S)))
  utils::write.csv(data.frame(comparison = sprintf("%s (t = %d) vs stored tcaCLR", opt$name, t2), reference = opt$compare, n_pairs = length(S), max_abs_diff = d, genes_identical = ok, pass = pass),
    file.path(opt$out, "selftest.csv"), row.names = FALSE)
  if (!pass) { say("SELF-TEST FAILED"); quit(status = 1) } }
say(sprintf("done: %.1f s, peak RSS %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), tc_vmhwm_gb()))
