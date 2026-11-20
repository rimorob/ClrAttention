#!/usr/bin/env Rscript
# tcaCLR2, step 1 (docs/plans/tcaclr2_definition.md): the permuted-data noise edge of the registered aCLR operator and the walk depth t2. EXPLORATORY; no truth is read.
#   Rscript analysis/tcaclr2/edge_depth.R --dataset m3d_tf_regulons [--out results/tcaclr2] [--workers 4] [--threads 2]
#   Rscript analysis/tcaclr2/edge_depth.R --matrix FILE --dataset NAME ...
# --dataset is one of the keys of GT in analysis/tcaclr_depth/grid_truths.R (k562_corum, net4_cp_yeast, armA_cp_human, armB_cyc2008, m3d_tf_regulons); the matrix is the one
# stage 2 of the depth grid scores (net4 from the DREAM5 files). --matrix overrides the file; --genes N takes a fixed random slice (smoke tests only).
# Rule (fixed by the definition note): operator keff 10, cap 50 (NOT 4 x keff), lazy 0.5, clipped CLR score = registered; 10 per-gene permutations (seeds AD$seed + 7000 + p,
# the same permuted matrices as spectral_edge.R); edge = mean over permutations of the largest non-trivial eigenvalue of P_sym (adaptive_lib.R, ad_Psym);
# t2 = max(1, round(ln 0.05 / ln edge)). The operator builder, symmetrisation, eigensolver and permutation are those of analysis/tcaclr_depth/adaptive/ (no second copy).
# Output: <out>/edge_<dataset>.tsv, one row: dataset G N edge edge_min edge_max t2 tstar (unrounded) + keff cap perms eps seed n_constant_genes top_eig_min (sanity: the trivial 1)
#         and <out>/edge_<dataset>.log.
suppressPackageStartupMessages({ source("analysis/tcaclr_depth/adaptive/adaptive_lib.R"); source("analysis/tcaclr_depth/grid_truths.R") })   # grid_truths.R sources tc_lib.R (TC)
ED <- list(keff = TC$keff, cap = TC$cap, alpha = TC$alpha, eps = 0.05, perms = 10L, seed = AD$seed, kvals = 6L)

ed_t2 <- function(edge, eps = ED$eps) { if (!is.finite(edge) || edge >= 1 || edge <= 0) stop("edge ", edge, " outside (0, 1): no depth rule"); max(1L, as.integer(round(log(eps) / log(edge)))) }

# per-permutation edge. All clr work happens inside the tasks (a forked child of a process that has run clr hangs: see ad_lapply).
ed_edge <- function(X, perms = ED$perms, threads = 2L, workers = 1L, seed = ED$seed, say = function(...) invisible(NULL)) {
  t0 <- Sys.time()
  one <- function(p) {
    f <- ad_fit(ad_permute_genes(X, seed + 7000L + p), threads); on.exit(f$release())
    Ps <- ad_Psym(ad_P_from_fit(f, ED$keff, cap = ED$cap, clipS = TRUE)); v <- ad_eigs(Ps, ED$kvals)$values
    c(top = v[1L], edge = v[2L], rss = gl_vmhwm_gb())
  }
  all <- ad_lapply(seq_len(perms), one, workers)
  bad <- vapply(all, function(x) !is.numeric(x), NA); if (any(bad)) stop("permutation task failed: ", all[bad][[1]])
  m <- do.call(rbind, all); say(sprintf("%d permuted spectra: %.0f s (worker peak RSS up to %.1f GB)", perms, as.numeric(difftime(Sys.time(), t0, units = "secs")), max(m[, "rss"])))
  m
}

if (sys.nframe() == 0L) {
  opt <- gl_options(commandArgs(TRUE), list(dataset = "", matrix = "", genes = 0L, perms = ED$perms, workers = 1L, threads = 2L, out = "results/tcaclr2"))
  if (!nzchar(opt$dataset)) stop("--dataset is required")
  if (!nzchar(opt$matrix) && !opt$dataset %in% names(GT)) stop("--dataset must be one of ", paste(names(GT), collapse = ", "), " (or give --matrix)")
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); say <- gl_say(file.path(opt$out, sprintf("edge_%s.log", opt$dataset)))
  say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "), "; rule: keff ", ED$keff, " cap ", ED$cap, " eps ", ED$eps)
  t0 <- Sys.time()
  X <- if (nzchar(opt$matrix)) gl_read_matrix(opt$matrix) else gl_truth_matrix(opt$dataset)
  if (opt$genes > 0L && opt$genes < nrow(X)) { set.seed(ED$seed); X <- X[sort(sample(nrow(X), opt$genes)), , drop = FALSE] }
  G <- nrow(X); N <- ncol(X); ncon <- sum(apply(X, 1, function(x) diff(range(x))) <= 0); say(sprintf("%s: %d genes x %d experiments; %d constant genes", opt$dataset, G, N, ncon))
  m <- ed_edge(X, opt$perms, opt$threads, opt$workers, ED$seed, say)
  edge <- mean(m[, "edge"]); t2 <- ed_t2(edge)
  r <- data.frame(dataset = opt$dataset, G = G, N = N, edge = edge, edge_min = min(m[, "edge"]), edge_max = max(m[, "edge"]), t2 = t2, tstar = log(ED$eps) / log(edge),
    keff = ED$keff, cap = ED$cap, perms = opt$perms, eps = ED$eps, seed = ED$seed, n_constant_genes = ncon, top_eig_min = min(m[, "top"]))
  f <- file.path(opt$out, sprintf("edge_%s.tsv", opt$dataset)); tmp <- paste0(f, ".tmp", Sys.getpid())
  utils::write.table(r, tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE); file.rename(tmp, f)
  say(sprintf("edge %.5f (range %.5f-%.5f); t2 = %d (unrounded %.2f); trivial eigenvalue min %.8f", edge, r$edge_min, r$edge_max, t2, r$tstar, r$top_eig_min))
  say(sprintf("done: %.1f s, peak RSS of the main process %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), gl_vmhwm_gb()))
}
