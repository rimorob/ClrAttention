#!/usr/bin/env Rscript
# tcaCLR2 on the scoreboard rows, step 1: the permuted-data edge and depth t2 (docs/plans/tcaclr2_definition.md, rule unchanged) for the matrices of the scoreboard that
# results/tcaclr2/ does not cover: EcoMAC (ComBat-adjusted, as the scoreboard), RPE1 (half A of the profiles, as clr_logpt_test.R) and the 14 BEELINE cases (CLR-NS operator).
# The edge machinery is edge_depth.R (ed_edge, ed_t2; sourced, not copied) for the MI-CLR matrices. EXPLORATORY; no truth is read.
#   Rscript analysis/tcaclr2/edge_scoreboard.R --id ecomac|rpe1|beeline [--out results/tcaclr2_scoreboard] [--workers 4] [--threads 2] [--genes N] [--cases hESC:500,...]
# --id beeline: for each of the 14 cases (7 datasets x TFs+500 / TFs+1000) the same rule with ONE substitution: the operator is built from CLR-NS (beeline.R::ns_clr, the score the
#   scoreboard's BEELINE cells use for aCLR and tcaCLR) instead of the MI-CLR score, the genes of the case matrix permuted independently (seeds AD$seed + 7000 + p), 10 permutations.
#   The definition note does not say what to do for a CLR-NS operator: this is the literal extension (flagged in the report). gl_make_P is checked equal to rz_operator + rz_lazy_P.
# Output (--out): edge_<id>.tsv (edge_depth.R format; BEELINE: edge_beeline_<dataset>_<N>.tsv per case) and edge_<id>.log.
suppressPackageStartupMessages({ source("analysis/tcaclr2/edge_depth.R"); source("analysis/tcaclr2/tsb_lib.R") })
topt <- tsb_args(list(id = "", out = TSB$out, workers = 4L, threads = 2L, genes = 0L, cases = "", perms = ED$perms))
if (!topt$id %in% c("ecomac", "rpe1", "beeline")) stop("--id must be ecomac, rpe1 or beeline")
tsb_refuse_frozen(topt$out); dir.create(topt$out, recursive = TRUE, showWarnings = FALSE)
say <- tsb_say(file.path(topt$out, sprintf("edge_%s.log", topt$id))); say("options: ", paste(names(topt), unlist(topt), sep = "=", collapse = " "))

write_edge <- function(id, X, m, extra = NULL) {
  edge <- mean(m[, "edge"]); t2 <- ed_t2(edge)
  r <- data.frame(dataset = id, G = nrow(X), N = ncol(X), edge = edge, edge_min = min(m[, "edge"]), edge_max = max(m[, "edge"]), t2 = t2, tstar = log(ED$eps) / log(edge),
    keff = ED$keff, cap = ED$cap, perms = nrow(m), eps = ED$eps, seed = ED$seed, n_constant_genes = sum(apply(X, 1, function(x) diff(range(x))) <= 0), top_eig_min = min(m[, "top"]))
  f <- file.path(topt$out, sprintf("edge_%s.tsv", id)); tmp <- paste0(f, ".tmp", Sys.getpid()); utils::write.table(r, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, f)
  say(sprintf("%s: edge %.5f (range %.5f-%.5f); t2 = %d (unrounded %.2f); trivial eigenvalue min %.8f", id, edge, r$edge_min, r$edge_max, t2, r$tstar, r$top_eig_min)); r
}
t0 <- Sys.time()
if (topt$id %in% c("ecomac", "rpe1")) {
  if (topt$id == "ecomac") { source("analysis/scoreboard/sb_lib.R"); sb_source(); X <- sb_load_ecoli("ecomac")$X }          # the matrix of sb_ecoli.R (ComBat, constant genes dropped)
  else { source("analysis/scoreboard/sb_cx_lib.R"); sb_import_truth_defs(); X <- load_expr("rpe1")$X }                      # the matrix of clr_logpt_test.R / sb_cx_score.R
  if (topt$genes > 0L && topt$genes < nrow(X)) { set.seed(ED$seed); X <- X[sort(sample(nrow(X), topt$genes)), , drop = FALSE] }
  say(sprintf("%s: %d genes x %d experiments; %d constant genes", topt$id, nrow(X), ncol(X), sum(apply(X, 1, function(x) diff(range(x))) <= 0)))
  m <- ed_edge(X, topt$perms, topt$threads, topt$workers, ED$seed, say); write_edge(topt$id, X, m)
} else {
  args <- c("--out", topt$out, "--threads", "1")                                                    # as sb_beeline.R: beeline.R's own option parser reads `args`
  src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]          # beeline.R's set-up and helpers, evaluated as sb_beeline.R does
  for (e in parse(text = src[seq_len(cut_at - 1L)])) { d1 <- deparse(e)[1]
    if (identical(d1, "args <- commandArgs(trailingOnly = TRUE)")) next
    if (startsWith(d1, "script_dir <-")) { script_dir <- "analysis"; next }
    eval(e, globalenv()) }
  source("analysis/robust_z/robust_z_lib.R"); say <- tsb_say(file.path(topt$out, sprintf("edge_%s.log", topt$id)))   # beeline.R redefines say
  Ns <- as.integer(strsplit(opt$N, ",")[[1]]); dss <- strsplit(opt$datasets, ",")[[1]]; cases <- expand.grid(N = Ns, ds = dss, stringsAsFactors = FALSE)
  if (nzchar(topt$cases)) { want <- strsplit(topt$cases, ",")[[1]]; cases <- cases[paste0(cases$ds, ":", cases$N) %in% want, ] }
  case_X <- function(ds, N) { ex <- as.matrix(utils::read.csv(file.path(expr_dir, ds, "ExpressionData.csv"), row.names = 1, check.names = FALSE))
    g <- intersect(select_genes(ds, N)$genes, rownames(ex)); X <- ex[g, , drop = FALSE]; X[apply(X, 1, function(x) diff(range(x))) > 0, , drop = FALSE] }
  for (k in seq_len(nrow(cases))) {
    ds <- cases$ds[k]; N <- cases$N[k]; X <- case_X(ds, N); id <- sprintf("beeline_%s_%d", ds, N); tk <- Sys.time()
    pm <- do.call(rbind, lapply(seq_len(topt$perms), function(p) { Xp <- ad_permute_genes(X, ED$seed + 7000L + p)
      C <- ns_clr(Xp); P <- gl_make_P(C, ED$keff, ED$cap)
      if (p == 1L) { d <- max(abs(P - rz_lazy_P(rz_operator(C)))); say(sprintf("%s: gl_make_P vs rz_lazy_P(rz_operator) max abs diff %.3g", id, d)); stopifnot(d < 1e-12) }
      v <- ad_eigs(ad_Psym(P), ED$kvals)$values; c(top = v[1L], edge = v[2L]) }))
    write_edge(id, X, pm); say(sprintf("%s: %d genes x %d cells, %.1f s", id, nrow(X), ncol(X), as.numeric(difftime(Sys.time(), tk, units = "secs"))))
  }
}
say(sprintf("done: %.1f s, peak RSS of the main process %.2f GB", as.numeric(difftime(Sys.time(), t0, units = "secs")), gl_vmhwm_gb()))
