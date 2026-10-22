#!/usr/bin/env Rscript
# Attention-width (keff) variants, exactly as preregistered in
# docs/preregistration_keff.md (fixed 2026-10-16, commit c0bb150).
#
# Shared: top-50 positive Stouffer CLR candidates per row, lazy walk alpha = 0.5,
# empty rows get a self-loop, S_t = (P^t + P^t')/2 (diagonal zeroed), the same
# seeds and lab splits as lab_splits.R.
#
# Variants (each built on the CLR matrix of the data at hand: a training half
# during selection, the full data for scoring):
#   K3, K5, K10, K25  clr:::.softmax_rows at keff 3, 5, 10 (reference), 25
#   MH     S_t = mean over heads keff 3, 10, 30 of each head's S_t (heads walk
#          separately, in lockstep)
#   O1     one global temperature, w_ij ∝ exp(z_ij / tau) over the cap-50
#          candidates; tau = median over rows of the bisected per-row tau_i of
#          the keff-10 softmax on the FULL-data CLR (one tau per dataset)
#   O1t    (secondary) tau in {0.25, ..., 4} chosen jointly with t (criterion A)
#   O2     per-gene keff_i = clamp(sum over the row's top-50 candidates of
#          (1 - lfdr_ij), 2, 30); lfdr = pi0 f0(z) / f(z), f and f0 kernel
#          densities of z > 0 for the observed CLR and for CLR on permuted data
#          (pooled over pairs); pi0 by Storey (lambda = 0.5) on permutation p;
#          z <= 0 gets lfdr = 1
#   O2e    (secondary) as O2, but lfdr per row by Efron central matching on the
#          unclipped row-standardized MI z
#   O3A    keff in {3,5,7,10,14,20,25,30,40} chosen jointly with t on criterion A
#   O3B    same, criterion B (primary O3)
# Held-out criteria (no truth used):
#   A  AP of S_t (upper triangle) for the other half's top-1% CLR pairs
#      (lab_splits.R criterion)
#   B  mean over genes of the per-row AP of S_t[i, -i] for gene i's top-20
#      held-out CLR partners
# Depth: each variant's own held-out t* (per split-direction argmax over the t
# grid; per split, geometric mean of the two directions; median over splits);
# secondary t = 9. Diagnostic only: full-data truth curves t = 1..32 with the
# truth-chosen best t labelled an upper bound.
# Splits: E. coli M3D 20 random balanced lab splits (lab_splits.R generator
# and seed: the first 20 of its splits), EcoMAC 20 balanced author splits
# (ecomac_batch.R generator), K562 10 random halves of the network-half-A
# perturbations; each split used in both directions.
# Truths (only after <out>/<dataset>/selections.csv is written):
#   E. coli: operons, TF regulons (operon pairs excluded), sigmulons (sigma70,
#   sigmaS excluded) -- definitions of granularity.R (RegulonDB, all evidence).
#   K562: CORUM, knockdown-response modules (half B), DoRothEA A-B, CollecTRI.
# Primary endpoint: pooled co-membership AUPR difference from K10 (each at its
# held-out t*), paired delete-half gene jackknife (100 halves, as
# baselines.R), on M3D TF regulons, EcoMAC TF regulons, K562 CollecTRI, K562
# CORUM; Holm over 8 variants x 4 endpoints = 32 (computed in <out> once the
# per-dataset primary.csv files exist; missing ones count in the 32).
#
# IMPLEMENTATION CHOICES (where the preregistration is silent; most literal reading):
#  1. Truth definitions for E. coli are granularity.R's (all evidence levels,
#     TF regulons 5-500 genes; sigmulons >= 5 genes without sigma70/sigma38);
#     operons are scored with same-operon pairs as the positives.
#  2. O1 primary tau: one number per dataset from the full-data CLR, also used
#     to build O1 on the training halves. tau_i exists only for rows with > 10
#     candidates (other rows are uniform under K10); the median is over those.
#  3. O2: one permutation of the data (every gene permuted independently over
#     samples) per matrix the operator is built on (each training half, and the
#     full data); null = all upper-triangle pairs of that CLR. "Conditioned on
#     z > 0" read literally: lfdr = min(1, pi0 * f0(z | z > 0) / f(z | z > 0)),
#     Gaussian kernel densities (bw.nrd0), evaluated on a 4096-point grid on
#     [0, max z]. Storey pi0 over all upper-triangle pairs, lambda = 0.5,
#     p = (1 + #null >= z) / (1 + #null).
#  4. O2e: per-row Efron central matching = quadratic fit of log KDE (512
#     points, bw.nrd0) over the row's central 50% (quartiles); pi0 capped at 1;
#     if the fitted curvature is not negative, delta0 = median, sigma0 =
#     IQR / 1.349, pi0 = 1. Candidates are the CLR top-50 (shared setting).
#  5. Joint selections (O1t, O3A, O3B): per split-direction joint argmax over
#     (parameter, t); per split, the geometric mean of the two directions for
#     both t (rounded) and the parameter; median over splits; the parameter is
#     snapped to the nearest grid value on the log scale. Ties: first maximum
#     (smaller parameter, smaller t).
#  6. Criterion B: partners = the top 20 POSITIVE held-out CLR scores of row i
#     (fewer if the row has fewer); rows with none are skipped; AP with the
#     repo tie convention (ap_fast).
#  7. K562 "coarse t grid": the same preregistered grid {1,...,32} as E. coli
#     (the only grid specified).
#  8. Row perplexity exp(entropy) is reported on the full-data attention
#     matrix A over rows with candidates (empty rows counted separately); for
#     MH, on the head-averaged A.
#  9. Holm significance is the two-sided z-test of the jackknife difference
#     against 0 (p = 2 Phi(-|d|/se), se = sd of the 100 paired half-sample
#     differences, as baselines.R); "margin >= 0.005" = point estimate beyond
#     +-0.005. A one-sided test against the margin itself is reported as an
#     extra column (p_beyond_margin), not used in the verdict.
# 10. Decision rule: win = Holm p < 0.05 and d >= +0.005; loss = Holm p < 0.05
#     and d <= -0.005; a variant meets the rule with >= 3 wins of 4, no loss,
#     and wins on both M3D and EcoMAC TF regulons. MH rule: win on K562 CORUM
#     and no loss on the regulon endpoints (M3D TF, EcoMAC TF, K562 CollecTRI).
# 11. Coherence: median per-regulon coherence AUROC (regulons.R); vs K10 a
#     paired bootstrap over regulons (B = 1,000) of the median difference.
# 12. Diagnostic (outside the design, from baselines.R): hub-stratified AUPR
#     (10 strata of pi_i + pi_j, pi = stationary distribution of the K10 walk).
# 13. Data loading is baselines.R's (copied: there it is interleaved with truth
#     construction, which must come later here); lab/author labels as
#     lab_splits.R / ecomac_batch.R. Functions lab_split, build_fit, lazy_P
#     (lab_splits.R), author_split (ecomac_batch.R) and hub_stratify
#     (baselines.R) are imported from those scripts by parsing.
# 14. Output files are per dataset in <out>/<dataset>/; the Holm table and
#     verdict in <out>/.
# 15. Truth gate (reviewer, 2026-10-21). Two runs per dataset:
#       selection phase  --stop_after_selection 1: splits, selections.csv, exit
#                        (never reads truths)
#       truth phase      --selections_ref docs/keff_selections/<dataset>_selections.csv:
#                        refused unless <out>/<dataset>/selections.csv existed
#                        when the process started AND is byte-identical to the
#                        committed reference; the selections recomputed from
#                        the saved split curves must also equal it.
#     --quick runs never read truths (they stop after the selections).
#     K562 truths are scored only at each variant's t* and t = 9 (no t = 1-32
#     truth curve); E. coli keeps the diagnostic curve.
#     cache/fulldata.rds and splits/fingerprint.rds carry a fingerprint
#     (dataset, quick, seed, G, md5 of the gene list; for splits also the
#     split assignments and t grid / cap / alpha / topB / target_frac) and
#     are not reused on mismatch.
#
# Usage (repo root):
#   Rscript analysis/keff_test.R --dataset m3d|ecomac|k562ess [--workers 4]
#       [--threads 2] [--quick 0] [--out results/keff_test] [--nsplit 20|10]
#       [--stop_after_selection 0] [--selections_ref <committed csv>] [--hubstrat 1]
#   Rscript analysis/keff_test.R --dataset combine [--out results/keff_test]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(dataset = "m3d", out = "results/keff_test", workers = 4L, threads = 2L, quick = 0L,
            nsplit = NA_integer_, jack = 100L, boot = 1000L, seed = NA_integer_, alpha = 0.5,
            target_frac = 0.01, cap = 50L, topB = 20L, tgrid = "1,2,3,4,5,6,8,10,12,16,20,24,32",
            tmax_curve = 32L, t_secondary = 9L, mem_gb = NA_real_, stop_after_selection = 0L,
            selections_ref = "", hubstrat = 1L, h5 = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad")
ints <- c("workers", "threads", "quick", "nsplit", "jack", "boot", "seed", "cap", "topB", "tmax_curve",
          "t_secondary", "stop_after_selection", "hubstrat")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else if (key %in% c("alpha", "target_frac", "mem_gb")) as.numeric(v) else v
  i <- i + 2L
}
ecoli <- opt$dataset %in% c("m3d", "ecomac")
if (is.na(opt$nsplit)) opt$nsplit <- if (ecoli) 20L else 10L
if (is.na(opt$seed)) opt$seed <- if (ecoli) 20261014L else 20261015L
if (is.na(opt$mem_gb)) opt$mem_gb <- if (ecoli) 3 else 8
tgrid <- as.integer(strsplit(opt$tgrid, ",")[[1]])
kgrid <- c(3, 5, 7, 10, 14, 20, 25, 30, 40)                       # O3 grid
taugrid <- c(0.25, 0.35, 0.5, 0.7, 1, 1.4, 2, 2.8, 4)              # O1-tuned grid
mh_heads <- c(3, 10, 30)
variants <- c("K3", "K5", "K10", "K25", "MH", "O1", "O1t", "O2", "O2e", "O3A", "O3B")
primary_variants <- c("K3", "K5", "K25", "MH", "O1", "O2", "O3A", "O3B")
endpoints <- data.frame(dataset = c("m3d", "ecomac", "k562ess", "k562ess"),
                        truth = c("TF_regulons", "TF_regulons", "CollecTRI", "CORUM"),
                        endpoint = c("M3D TF regulons", "EcoMAC TF regulons", "K562 CollecTRI", "K562 CORUM"),
                        stringsAsFactors = FALSE)
margin <- 0.005

say_to <- function(f) function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = f, append = TRUE) }

## ---- Holm over the 32 primary tests and the decision rule --------------------------------
combine <- function(out, say) {
  fs <- file.path(out, c("m3d", "ecomac", "k562ess"), "primary.csv")
  pr <- do.call(rbind, lapply(fs[file.exists(fs)], utils::read.csv, stringsAsFactors = FALSE))
  grid <- merge(expand.grid(variant = primary_variants, endpoint = endpoints$endpoint, stringsAsFactors = FALSE),
                endpoints, by = "endpoint")
  if (!is.null(pr)) grid <- merge(grid, pr[, c("variant", "endpoint", "t_star", "d", "se", "p")],
                                  by = c("variant", "endpoint"), all.x = TRUE)
  else { grid$t_star <- NA; grid$d <- NA; grid$se <- NA; grid$p <- NA }
  grid$p_holm <- NA_real_
  ok <- !is.na(grid$p)
  grid$p_holm[ok] <- stats::p.adjust(grid$p[ok], "holm", n = 32L)    # missing tests count in the 32
  grid$win <- ok & grid$p_holm < 0.05 & grid$d >= margin
  grid$loss <- ok & grid$p_holm < 0.05 & grid$d <= -margin
  grid <- grid[order(match(grid$variant, primary_variants), match(grid$endpoint, endpoints$endpoint)), ]
  utils::write.csv(grid, file.path(out, "primary_holm.csv"), row.names = FALSE)
  complete <- all(ok)
  v <- do.call(rbind, lapply(primary_variants, function(vn) {
    x <- grid[grid$variant == vn, ]; w <- stats::setNames(x$win, x$endpoint); l <- stats::setNames(x$loss, x$endpoint)
    meets <- sum(w) >= 3 && !any(l) && w[["M3D TF regulons"]] && w[["EcoMAC TF regulons"]]
    mh <- if (vn == "MH") w[["K562 CORUM"]] && !any(l[c("M3D TF regulons", "EcoMAC TF regulons", "K562 CollecTRI")]) else NA
    data.frame(variant = vn, n_tests_present = sum(!is.na(x$p)), wins = sum(w), losses = sum(l),
               wins_both_ecoli = w[["M3D TF regulons"]] && w[["EcoMAC TF regulons"]],
               meets_rule = meets, mh_corum_rule = mh, stringsAsFactors = FALSE) }))
  v$complete <- complete
  if (!is.null(pr) && any(pr$quick != 0)) complete <- FALSE           # --quick smoke runs never give a verdict
  meet <- v$variant[v$meets_rule]
  verdict <- if (!is.null(pr) && any(pr$quick != 0)) "SMOKE TEST (--quick): no verdict" else
    if (!complete) "INCOMPLETE: not all 32 primary tests are present (missing tests counted in Holm n = 32); no verdict" else
    if (length(meet)) paste0("headline changes: ", paste(meet, collapse = ", "), " meet the rule",
                             if (any(meet %in% c("O1", "O2"))) paste0(" (", paste(intersect(meet, c("O1", "O2")), collapse = ", "), " becomes the method)") else "") else
    "headline unchanged (K10): no variant meets the decision rule"
  if (complete && isTRUE(v$mh_corum_rule[v$variant == "MH"])) verdict <- paste0(verdict, "; MH gains on CORUM without losing on regulons: MH becomes the main operator")
  v$verdict <- verdict
  utils::write.csv(v, file.path(out, "verdict.csv"), row.names = FALSE)
  say("primary endpoint, variant - K10 (Holm over 32):")
  for (k in seq_len(nrow(grid))) with(grid[k, ], say(sprintf("  %-4s %-20s d %+.4f se %.4f p %.3g p_holm %.3g%s",
    variant, endpoint, d, se, p, p_holm, if (isTRUE(win)) "  WIN" else if (isTRUE(loss)) "  LOSS" else "")))
  say("verdict: ", verdict)
  invisible(v)
}
if (opt$dataset == "combine") {
  dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
  combine(opt$out, say_to(file.path(opt$out, "keff_test_combine.log"))); quit(save = "no")
}
if (!opt$dataset %in% c("m3d", "ecomac", "k562ess")) stop("unknown dataset: ", opt$dataset)

dout <- file.path(opt$out, opt$dataset)
sel_f <- file.path(dout, "selections.csv")
SEL_EXISTED_AT_START <- file.exists(sel_f)                  # checked before anything is written
truth_run <- opt$quick == 0L && opt$stop_after_selection == 0L
# Truth-phase precondition: the frozen selections existed at start and equal the committed reference.
check_gate <- function() {
  if (!SEL_EXISTED_AT_START) stop("truth phase refused: ", sel_f, " did not exist when this process started ",
    "(run the selection phase with --stop_after_selection 1, commit the file, then rerun with --selections_ref)")
  if (!nzchar(opt$selections_ref) || !file.exists(opt$selections_ref))
    stop("truth phase refused: committed reference --selections_ref '", opt$selections_ref, "' not given or not found")
  rd <- function(f) readBin(f, "raw", file.info(f)$size)
  if (!identical(rd(sel_f), rd(opt$selections_ref)))
    stop("truth phase refused: ", sel_f, " is not byte-identical to ", opt$selections_ref)
  invisible(TRUE)
}
if (truth_run) check_gate()                                  # fail fast, before any compute
dir.create(file.path(dout, "splits"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dout, "cache"), showWarnings = FALSE)
logf <- file.path(dout, "keff_test.log")
say <- say_to(logf)
wcsv <- function(x, n) utils::write.csv(x, file.path(dout, n), row.names = FALSE)
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R")); source(file.path(script_dir, "regulons.R"))
source(file.path(script_dir, "human_truths.R")); source(file.path(script_dir, "parallel.R"))
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

# Import function definitions (only `name <- function(...)` at top level) from
# the analysis scripts, so the shared constructions are not duplicated.
import_defs <- function(file, nms) {
  got <- character(0)
  for (e in parse(file.path(script_dir, file), keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% nms && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      eval(e, globalenv()); got <- c(got, as.character(e[[2]])) }
  if (length(setdiff(nms, got))) stop("not found in ", file, ": ", paste(setdiff(nms, got), collapse = ", "))
}
import_defs("lab_splits.R", c("lab_split", "build_fit", "lazy_P"))   # build_fit uses nthr(); lazy_P uses opt$alpha
import_defs("ecomac_batch.R", "author_split")
import_defs("baselines.R", "hub_stratify")
THREADS <- opt$workers * opt$threads                     # master phases; workers get opt$threads
nthr <- function() THREADS

## ---- data (as baselines.R; no truths here) ----------------------------------------------
if (opt$dataset == "m3d") {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
  desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
  exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description)        # labs as lab_splits.R
  exper[!grepl("experimenter:", desc$description)] <- NA
  lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])
  lab <- lab[colnames(X)]; lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab))); rm(desc, exper)
} else if (opt$dataset == "ecomac") {
  X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
  rownames(X) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
  plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
  author <- utils::read.delim("data/ecomac/samples.tsv", stringsAsFactors = FALSE, check.names = FALSE)[["Author"]]
  stopifnot(length(author) == ncol(X), length(plat) == ncol(X))
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  X <- suppressMessages(sva::ComBat(X, batch = factor(plat)))
} else {
  set.seed(20261015L)                                   # perturbseq_bulk.R split
  rd <- function(p) rhdf5::h5read(opt$h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
    if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
    if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  Z <- rd("X"); rownames(Z) <- cat_of("/var", "gene_name"); pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
  pval <- as.numeric(rd("/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
  X <- Z[, A]; rm(Z); X[!is.finite(X)] <- 0
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  dup <- rownames(X) %in% rownames(X)[duplicated(rownames(X))]
  if (any(dup)) { say(sprintf("dropped %d rows with duplicated gene names", sum(dup))); X <- X[!dup, ] }
  rm(rd, cat_of, pert, pval, K, A, dup)
}
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (opt$quick > 0L) { say(sprintf("--quick %d: top-variance genes of %d", opt$quick, nrow(X)))
  v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
storage.mode(X) <- "double"
genes <- rownames(X); G <- length(genes)
fingerprint <- list(dataset = opt$dataset, quick = opt$quick, seed = opt$seed, G = G,
                    genes_md5 = digest::digest(genes, algo = "md5"))
same_fp <- function(a, b) isTRUE(all.equal(a, b))
say(sprintf("%s: %d genes x %d samples", opt$dataset, G, ncol(X)))

## ---- splits -------------------------------------------------------------------------------
set.seed(opt$seed)
splits <- lapply(seq_len(opt$nsplit), function(s) switch(opt$dataset,
  m3d = lab_split(), ecomac = author_split(), k562ess = sample(rep(1:2, length.out = ncol(X)))))
fp_splits <- c(fingerprint, list(splits_md5 = digest::digest(splits, algo = "md5"), tgrid = tgrid, cap = opt$cap,
                                 alpha = opt$alpha, topB = opt$topB, target_frac = opt$target_frac))
fp_f <- file.path(dout, "splits", "fingerprint.rds")
if (file.exists(fp_f)) { if (!same_fp(readRDS(fp_f), fp_splits)) stop("splits/ fingerprint mismatch (", fp_f, "): not reusing; use a new --out")
} else { if (length(list.files(file.path(dout, "splits"), "\\.csv$"))) stop("splits/ has results without a fingerprint: not reusing")
  saveRDS(fp_splits, fp_f) }
wcsv(data.frame(sample = colnames(X), do.call(cbind, stats::setNames(splits, sprintf("split_%02d", seq_along(splits))))),
     "split_assignments.csv")
say(sprintf("%d %s splits; half sizes (first split): %d / %d", opt$nsplit,
            switch(opt$dataset, m3d = "lab", ecomac = "author", k562ess = "perturbation-half"),
            sum(splits[[1]] == 1L), sum(splits[[1]] == 2L)))

## ---- operators ----------------------------------------------------------------------------
selfloop <- function(A) { e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1; A }
soft_k <- function(C, k) selfloop(clr:::.softmax_rows(C, keff = k, cap = opt$cap))
cand <- function(row, i, n) { row[i] <- -Inf; top <- order(row, decreasing = TRUE)[seq_len(min(opt$cap, n - 1L))]
                              top[row[top] > 0] }
bisect_lt <- function(s, keff) {                       # as clr:::.softmax_rows
  neff <- function(lt) { w <- exp((s - max(s)) / exp(lt)); w <- w / sum(w); exp(-sum(w[w > 0] * log(w[w > 0]))) }
  lo <- log(1e-4); hi <- log(1e4)
  for (it in 1:60) { mid <- (lo + hi) / 2; if (neff(mid) < keff) lo <- mid else hi <- mid }
  (lo + hi) / 2
}
# Softmax over the cap-50 positive candidates with a per-row keff vector (O2,
# O2e; identical to .softmax_rows for a constant vector) or one global tau (O1).
softmax_custom <- function(C, keff = NULL, tau = NULL) {
  n <- nrow(C); A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    top <- cand(C[i, ], i, n); if (!length(top)) next
    s <- C[i, top]
    w <- if (!is.null(tau)) exp((s - max(s)) / tau) else
      if (length(top) <= keff[i]) rep(1, length(top)) else exp((s - max(s)) / exp(bisect_lt(s, keff[i])))
    A[i, top] <- w / sum(w)
  }
  selfloop(A)
}
tau_rows <- function(C, keff = 10) {                    # bisected per-row temperatures (rows with > keff candidates)
  n <- nrow(C); unlist(lapply(seq_len(n), function(i) { top <- cand(C[i, ], i, n)
    if (length(top) <= keff) NULL else exp(bisect_lt(C[i, top], keff)) }))
}
cand_sum <- function(C, one_minus) {                    # per row: sum of (1 - lfdr) over the candidates
  n <- nrow(C); vapply(seq_len(n), function(i) { top <- cand(C[i, ], i, n)
    if (!length(top)) 0 else sum(one_minus(i, top)) }, 0)
}
clamp <- function(x) pmin(pmax(x, 2), 30)
# O2: local FDR from a permutation null (CLR on permuted data), pooled over pairs.
keff_lfdr <- function(C, C0) {
  ut <- upper.tri(C); z <- C[ut]; z0 <- sort(C0[ut]); N0 <- length(z0)
  p <- (1 + N0 - findInterval(z, z0, left.open = TRUE)) / (1 + N0)
  pi0 <- min(1, mean(p > 0.5) / 0.5); rm(p)
  zp <- z[z > 0]; z0p <- z0[z0 > 0]; rm(z, z0)
  hi <- max(zp, z0p)
  d1 <- stats::density(zp, from = 0, to = hi, n = 4096); d0 <- stats::density(z0p, from = 0, to = hi, n = 4096)
  lf <- function(x) { f <- stats::approx(d1$x, d1$y, x, rule = 2)$y; f0 <- stats::approx(d0$x, d0$y, x, rule = 2)$y
    out <- ifelse(f > 0, pmin(1, pi0 * f0 / f), 1); out[x <= 0] <- 1; out }
  keff <- clamp(cand_sum(C, function(i, j) 1 - lf(C[i, j])))
  list(keff = keff, pi0 = pi0)
}
# O2e: per-row Efron central matching on the unclipped row-standardized MI z.
keff_efron <- function(M, C) {
  n <- nrow(M)
  one <- function(i, j) {
    m <- M[i, ]; mu <- mean(m[-i]); sdv <- stats::sd(m[-i]); z <- (m - mu) / sdv
    zz <- z[-i]; d <- stats::density(zz, n = 512); q <- stats::quantile(zz, c(0.25, 0.75), names = FALSE)
    sel <- d$x >= q[1] & d$x <= q[2] & d$y > 0
    b <- if (sum(sel) >= 3) stats::.lm.fit(cbind(1, d$x[sel], d$x[sel]^2), log(d$y[sel]))$coefficients else c(NA, NA, NA)
    if (all(is.finite(b)) && b[3] < 0) { s2 <- -1 / (2 * b[3]); d0 <- b[2] * s2
      pi0 <- min(1, exp(b[1] + d0^2 / (2 * s2) + log(sqrt(2 * pi * s2)))); s0 <- sqrt(s2)
    } else { d0 <- stats::median(zz); s0 <- diff(q) / 1.349; pi0 <- 1 }
    f <- stats::approx(d$x, d$y, z[j], rule = 2)$y
    1 - ifelse(f > 0, pmin(1, pi0 * stats::dnorm(z[j], d0, s0) / f), 1)
  }
  clamp(cand_sum(C, one))
}
perm_rows <- function(Xm) { P <- t(apply(Xm, 1, sample)); dimnames(P) <- dimnames(Xm); P }
# CLR, O2 and O2e keff vectors of one data matrix (a training half or the full data).
prep <- function(Xm, seed) {
  set.seed(seed)
  f <- build_fit(Xm); C <- f$clr_scores; M <- f$mi; f$release(); rm(f)
  ke <- keff_efron(M, C); rm(M); invisible(gc(verbose = FALSE))
  f0 <- build_fit(perm_rows(Xm)); C0 <- f0$clr_scores; f0$release(); rm(f0)
  lf <- keff_lfdr(C, C0); rm(C0); invisible(gc(verbose = FALSE))
  list(C = C, keff_O2 = lf$keff, pi0_O2 = lf$pi0, keff_O2e = ke)
}
# Attention heads of a variant at given parameters; S_t is the mean over heads.
heads_of <- function(v, C, par, keff_O2 = NULL, keff_O2e = NULL) switch(v,
  K3 = list(soft_k(C, 3)), K5 = list(soft_k(C, 5)), K10 = list(soft_k(C, 10)), K25 = list(soft_k(C, 25)),
  MH = lapply(mh_heads, function(k) soft_k(C, k)), O1 = list(softmax_custom(C, tau = par)),
  O1t = list(softmax_custom(C, tau = par)), O2 = list(softmax_custom(C, keff = keff_O2)),
  O2e = list(softmax_custom(C, keff = keff_O2e)), O3A = list(soft_k(C, par)), O3B = list(soft_k(C, par)),
  stop("unknown variant ", v))
# Walk all heads in lockstep; fn(S, t) with S = mean_h (P_h^t + P_h^t')/2 at t in ts.
walk_mass <- function(As, ts, fn) {
  Ps <- lapply(As, function(A) lazy_P(as.matrix(A))); rm(As); H <- length(Ps); Pts <- vector("list", H)
  for (t in seq_len(max(ts))) {
    for (h in seq_len(H)) Pts[[h]] <- if (t == 1L) as.matrix(Ps[[h]]) else as.matrix(Pts[[h]] %*% Ps[[h]])
    if (t %in% ts) { S <- if (H == 1L) Pts[[1]] else Reduce(`+`, Pts) / H
                     S <- (S + t(S)) / 2; diag(S) <- 0; fn(S, t); rm(S) }
    if (t %% 8L == 0L) invisible(gc(verbose = FALSE))
  }
  invisible(NULL)
}
top_partners <- function(C) lapply(seq_len(nrow(C)), function(i) { r <- C[i, ]; r[i] <- -Inf
  o <- order(r, decreasing = TRUE)[seq_len(min(opt$topB, nrow(C) - 1L))]; o[r[o] > 0] })
crit_B <- function(S, tB) {
  n <- nrow(S); lab <- integer(n)
  mean(unlist(lapply(seq_len(n), function(i) { if (!length(tB[[i]])) return(NULL)
    lab[] <- 0L; lab[tB[[i]]] <- 1L; ap_fast(S[i, -i], lab[-i]) })))
}
rowperp <- function(A) { L <- A * log(A); L[A <= 0] <- 0; exp(-rowSums(L)) }

## ---- full-data label-free quantities (cached): CLR, O1 tau, O2 / O2e keff ----------------
cache_f <- file.path(dout, "cache", "fulldata.rds")
if (file.exists(cache_f)) { full <- readRDS(cache_f)
  if (!same_fp(full$fingerprint, fingerprint)) stop("cache fingerprint mismatch (", cache_f, "): not reusing; use a new --out")
  say("full-data CLR and keff vectors reused from ", cache_f)
} else {
  t0 <- Sys.time(); full <- prep(X, opt$seed + 999L)
  full$tau_rows <- tau_rows(full$C, 10); full$tau_O1 <- stats::median(full$tau_rows); full$fingerprint <- fingerprint
  saveRDS(full, cache_f)
  say(sprintf("full-data CLR, permutation null, O2/O2e keff in %.1f min",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
tau_O1 <- full$tau_O1
chk <- max(abs(softmax_custom(full$C, keff = rep(10, G)) - soft_k(full$C, 10)))
if (chk > 1e-12) stop("softmax_custom differs from clr:::.softmax_rows by ", chk)
say(sprintf("O1 width-matched tau = %.4f (median over %d rows of the keff-10 tau_i; IQR %.3f-%.3f)", tau_O1,
            length(full$tau_rows), stats::quantile(full$tau_rows, 0.25), stats::quantile(full$tau_rows, 0.75)))
say(sprintf("O2 full data: pi0 = %.3f, keff median %.1f [%.1f, %.1f]; O2e keff median %.1f [%.1f, %.1f]",
            full$pi0_O2, stats::median(full$keff_O2), min(full$keff_O2), max(full$keff_O2),
            stats::median(full$keff_O2e), min(full$keff_O2e), max(full$keff_O2e)))
rm(full); invisible(gc(verbose = FALSE))

## ---- held-out criterion curves, one task per split (no truths) -----------------------------
run_split <- function(si) {
  h <- splits[[si]]
  hp <- lapply(1:2, function(d) prep(X[, h == d, drop = FALSE], opt$seed + 100L * si + d))
  ut <- which(upper.tri(matrix(0, G, G)))
  rows <- list()
  for (dir in 1:2) {
    tr <- dir; te <- 3L - dir
    sv <- hp[[te]]$C[ut]; tA <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)); rm(sv)
    tB <- top_partners(hp[[te]]$C)
    Cr <- hp[[tr]]$C
    run <- function(variant, param, As, withB = FALSE) {
      As <- lapply(As, function(A) Matrix::Matrix(A, sparse = TRUE))   # drop the dense heads before walking
      cr <- matrix(NA_real_, length(tgrid), 2)
      walk_mass(As, tgrid, function(S, t) { k <- match(t, tgrid)
        cr[k, 1] <<- ap_fast(S[ut], tA); if (withB) cr[k, 2] <<- crit_B(S, tB) })
      rows[[length(rows) + 1L]] <<- data.frame(split = si, direction = dir, variant = variant, param = param,
                                               t = tgrid, crit_A = cr[, 1], crit_B = cr[, 2])
    }
    for (k in kgrid) run("keff", k, list(soft_k(Cr, k)), withB = TRUE)   # K3, K5, K10, K25, O3A, O3B
    run("MH", NA, heads_of("MH", Cr))
    run("O1", tau_O1, heads_of("O1", Cr, tau_O1))
    for (tau in taugrid) run("O1t", tau, heads_of("O1t", Cr, tau))
    run("O2", NA, heads_of("O2", Cr, keff_O2 = hp[[tr]]$keff_O2))
    run("O2e", NA, heads_of("O2e", Cr, keff_O2e = hp[[tr]]$keff_O2e))
  }
  out <- do.call(rbind, rows)
  out$pi0_O2_train <- vapply(out$direction, function(d) hp[[d]]$pi0_O2, 0)
  out$keff_O2_train_median <- vapply(out$direction, function(d) stats::median(hp[[d]]$keff_O2), 0)
  out$keff_O2e_train_median <- vapply(out$direction, function(d) stats::median(hp[[d]]$keff_O2e), 0)
  out
}
done <- sub("\\.csv$", "", list.files(file.path(dout, "splits"), "\\.csv$"))
todo <- setdiff(sprintf("split_%02d", seq_len(opt$nsplit)), done)
say(sprintf("%d splits done, %d to run; t grid %s", length(done), length(todo), opt$tgrid))
reserve_for <- function(W) max(0L, .par_cores() - W * opt$threads)   # workers x threads = the declared cores
if (length(todo)) {
  W <- min(opt$workers, length(todo)); THREADS <- opt$threads
  pc <- par_start(mem_gb = opt$mem_gb, workers = W, reserve = reserve_for(W), say = say,
                  export = setdiff(ls(globalenv()), c("pc")))
  st <- foreach(nm = todo, .inorder = FALSE, .errorhandling = "pass") %dopar% {
    if (par_should_stop()) return("stopped")
    si <- as.integer(sub("split_", "", nm)); t0 <- Sys.time()
    d <- run_split(si)
    tmp <- file.path(dout, "splits", paste0(nm, ".csv.part"))
    utils::write.csv(d, tmp, row.names = FALSE); file.rename(tmp, file.path(dout, "splits", paste0(nm, ".csv")))
    par_release()
    sprintf("%s %.1f min (peak %.0f MB)", nm, as.numeric(difftime(Sys.time(), t0, units = "mins")),
            sum(gc()[, 6]))
  }
  par_stop(pc); THREADS <- opt$workers * opt$threads
  bad <- vapply(st, function(x) inherits(x, "error"), NA)
  if (any(bad)) say("splits failed: ", paste(vapply(st[bad], conditionMessage, ""), collapse = " | "))
  say("split timings: ", paste(unlist(st[!bad]), collapse = "; "))
  if (any(bad) || any(unlist(st[!bad]) == "stopped")) stop("not all splits finished; rerun to resume")
}
cur <- do.call(rbind, lapply(sprintf("%s/splits/split_%02d.csv", dout, seq_len(opt$nsplit)), utils::read.csv))
wcsv(cur, "heldout_curves.csv")

## ---- selections (frozen before any truth is read) --------------------------------------------
snap <- function(x, grid) grid[which.min(abs(log(grid) - log(x)))]
select <- function(x, crit, joint, grid = NULL) {
  per_sd <- do.call(rbind, lapply(split(x, list(x$split, x$direction), drop = TRUE), function(y) {
    y <- y[order(y$param, y$t), ]; k <- which.max(y[[crit]])
    data.frame(split = y$split[1], direction = y$direction[1], t_best = y$t[k], param_best = y$param[k],
               crit_best = y[[crit]][k]) }))
  per_s <- do.call(rbind, lapply(split(per_sd, per_sd$split), function(y)
    data.frame(split = y$split[1], t_split = as.integer(round(exp(mean(log(y$t_best))))),
               param_split = if (joint) exp(mean(log(y$param_best))) else y$param_best[1])))
  list(per_sd = per_sd, per_s = per_s, t_star = as.integer(floor(stats::median(per_s$t_split) + 0.5)),
       param = if (joint) snap(stats::median(per_s$param_split), grid) else per_s$param_split[1])
}
spec <- list(K3 = list("keff", 3, "crit_A", FALSE), K5 = list("keff", 5, "crit_A", FALSE),
             K10 = list("keff", 10, "crit_A", FALSE), K25 = list("keff", 25, "crit_A", FALSE),
             MH = list("MH", NA, "crit_A", FALSE), O1 = list("O1", NA, "crit_A", FALSE),
             O1t = list("O1t", NA, "crit_A", TRUE), O2 = list("O2", NA, "crit_A", FALSE),
             O2e = list("O2e", NA, "crit_A", FALSE), O3A = list("keff", NA, "crit_A", TRUE),
             O3B = list("keff", NA, "crit_B", TRUE))
sels <- list(); persplit <- list()
for (v in variants) {
  s <- spec[[v]]; x <- cur[cur$variant == s[[1]], ]
  if (!is.na(s[[2]])) x <- x[x$param == s[[2]], ]
  r <- select(x, s[[3]], s[[4]], grid = if (v == "O1t") taugrid else kgrid)
  pname <- switch(v, O1 = "tau", O1t = "tau", O3A = "keff", O3B = "keff", K3 = "keff", K5 = "keff", K10 = "keff",
                  K25 = "keff", MH = "keff_heads", O2 = "keff_per_gene", O2e = "keff_per_gene")
  pval <- switch(v, O1 = tau_O1, O1t = r$param, O3A = r$param, O3B = r$param, MH = NA, O2 = NA, O2e = NA, s[[2]])
  sels[[v]] <- data.frame(dataset = opt$dataset, variant = v, criterion = s[[3]], param_name = pname,
    param_value = pval, t_star = r$t_star, t_q25 = stats::quantile(r$per_s$t_split, 0.25, names = FALSE),
    t_q75 = stats::quantile(r$per_s$t_split, 0.75, names = FALSE), t_min = min(r$per_s$t_split),
    t_max = max(r$per_s$t_split), param_split_min = if (s[[4]]) min(r$per_s$param_split) else NA,
    param_split_max = if (s[[4]]) max(r$per_s$param_split) else NA, n_splits = nrow(r$per_s),
    t_secondary = opt$t_secondary, quick = opt$quick, stringsAsFactors = FALSE)
  persplit[[v]] <- cbind(variant = v, r$per_sd)
}
sel <- do.call(rbind, sels)
if (file.exists(sel_f)) {                                       # a frozen (committed) selection exists
  old <- utils::read.csv(sel_f, stringsAsFactors = FALSE)
  same <- nrow(old) == nrow(sel) && all(old$variant == sel$variant) && all(old$t_star == sel$t_star) &&
    isTRUE(all.equal(old$param_value, sel$param_value, tolerance = 1e-8))
  if (!same) stop("recomputed selections differ from the frozen ", sel_f, "; truths not touched")
  say("selections recomputed and identical to the frozen ", sel_f)
} else wcsv(sel, "selections.csv")
wcsv(do.call(rbind, persplit), "selection_per_split_direction.csv")
# O3 width bias: per keff, mean over split-directions of the best criterion over t, and how often it wins
wb <- do.call(rbind, lapply(c("crit_A", "crit_B"), function(cr) {
  x <- cur[cur$variant == "keff", ]
  best <- stats::aggregate(x[[cr]], list(split = x$split, direction = x$direction, keff = x$param), max)
  win <- persplit[[if (cr == "crit_A") "O3A" else "O3B"]]
  data.frame(criterion = cr, keff = kgrid, mean_best_crit = vapply(kgrid, function(k) mean(best$x[best$keff == k]), 0),
             n_split_dirs_selecting = vapply(kgrid, function(k) sum(win$param_best == k), 0L)) }))
wcsv(wb, "o3_width_bias.csv")
say("held-out selections (frozen in ", sel_f, "):")
for (k in seq_len(nrow(sel))) with(sel[k, ], say(sprintf("  %-4s %-6s %-14s %-8s t* %2d (IQR %g-%g, range %d-%d)",
  variant, criterion, param_name, ifelse(is.na(param_value), "", format(signif(param_value, 4))), t_star, t_q25, t_q75, t_min, t_max)))
if (opt$quick > 0L) { say("--quick: stopping after the selections (quick runs never read truths)"); quit(save = "no") }
if (opt$stop_after_selection == 1L) { say("--stop_after_selection 1: selections frozen; stopping before truths"); quit(save = "no") }
check_gate(); TRUTH_GATE <- TRUE                             # re-checked at the gate
say("truth gate passed: ", sel_f, " existed at start and is byte-identical to ", opt$selections_ref)

## ---- full-data operators, row perplexities (label-free) -----------------------------------------
full <- readRDS(cache_f)
selv <- function(v, col) sel[[col]][sel$variant == v]
As <- lapply(stats::setNames(variants, variants), function(v)
  heads_of(v, full$C, selv(v, "param_value"), full$keff_O2, full$keff_O2e))
perp <- do.call(rbind, lapply(variants, function(v) {
  A <- if (length(As[[v]]) == 1L) As[[v]][[1]] else Reduce(`+`, As[[v]]) / length(As[[v]])
  e <- rowSums(A > 0) == 1L & diag(A) == 1; p <- rowperp(A)[!e]; q <- stats::quantile(p, c(0, .05, .25, .5, .75, .95, 1))
  data.frame(variant = v, n_rows = length(p), n_empty_rows = sum(e), mean = mean(p), min = q[1], q05 = q[2], q25 = q[3],
             median = q[4], q75 = q[5], q95 = q[6], max = q[7], frac_below_3 = mean(p < 3), frac_above_30 = mean(p > 30)) }))
wcsv(perp, "row_perplexity.csv")
kq <- function(k, nm) data.frame(variant = nm, t(stats::quantile(k, c(0, .05, .25, .5, .75, .95, 1))), mean = mean(k),
                                 frac_at_2 = mean(k == 2), frac_at_30 = mean(k == 30), check.names = FALSE)
wcsv(rbind(kq(full$keff_O2, "O2"), kq(full$keff_O2e, "O2e")), "o2_keff_distribution.csv")
say("row perplexity (full data): ", paste(sprintf("%s med %.1f <3 %.2f >30 %.2f", perp$variant, perp$median,
                                                  perp$frac_below_3, perp$frac_above_30), collapse = "; "))
As <- lapply(As, function(hl) lapply(hl, function(A) Matrix::Matrix(A, sparse = TRUE)))
# stationary distribution of the K10 walk (teleport 1e-3), for the hub-stratified diagnostic (as baselines.R)
P10 <- lazy_P(as.matrix(As$K10[[1]])); q_pi <- rep(1 / G, G)
for (it in 1:20000) { q <- (1 - 1e-3) * as.numeric(q_pi %*% P10) + 1e-3 / G; if (sum(abs(q - q_pi)) < 1e-13) break; q_pi <- q }
q_pi <- q; rm(P10, full, q); invisible(gc(verbose = FALSE))

## ---- truths: read only now --------------------------------------------------------------------
load_truths <- function() {
  if (!isTRUE(TRUTH_GATE) || !truth_run) stop("truths requested before the selections were frozen")
  check_gate()
  none <- stats::setNames(character(0), character(0))
  if (ecoli) {
    to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
    rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
    operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
    reg_all <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)   # as granularity.R
    pick <- function(keep) list(members = reg_all$members[keep], class = reg_all$class[keep])
    tf_reg <- pick(reg_all$class == "TF" & lengths(reg_all$members) <= 500 & lengths(reg_all$members) >= 5)
    sig_reg <- pick(reg_all$class == "sigma" & lengths(reg_all$members) >= 5 &
                      !names(reg_all$members) %in% c("sigma:sigma70", "sigma:sigma38"))
    op <- operon_of[names(operon_of) %in% genes]
    op_sets <- lapply(split(match(names(op), genes), op), unique); op_sets <- op_sets[lengths(op_sets) >= 2]
    op_reg <- list(members = op_sets, class = stats::setNames(rep("operon", length(op_sets)), names(op_sets)))
    mk <- function(r, o) list(reg = r, pairs = comembership_pairs(r, o, genes), co = comember_matrix(r, G), operon = o)
    list(operons = mk(op_reg, none), TF_regulons = mk(tf_reg, operon_of), sigmulons = mk(sig_reg, operon_of))
  } else {
    mk <- function(r) list(reg = r, pairs = comembership_pairs(r, none, genes), co = comember_matrix(r, G), operon = none)
    list(CORUM = mk(corum_sets(genes)), KD_response = mk(response_sets(genes, path = opt$h5, use = B)),
         DoRothEA_AB = mk(dorothea_sets(genes)), CollecTRI = mk(collectri_sets(genes)))
  }
}
bench <- load_truths()
bench <- bench[vapply(bench, function(b) sum(b$pairs$label) > 0, NA)]
for (nm in names(bench)) say(sprintf("truth %-12s %4d sets, %d pairs, base rate %.4f", nm,
  length(bench[[nm]]$reg$members), length(bench[[nm]]$pairs$label), mean(bench[[nm]]$pairs$label)))
set.seed(20261016L)                                               # baselines.R jackknife seed
halves <- lapply(bench, function(b) lapply(seq_len(opt$jack), function(r) {
  keep <- sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2))
  inh <- logical(G); inh[keep] <- TRUE; which(inh[b$pairs$i] & inh[b$pairs$j]) }))

## ---- full-data scoring, one task per variant ----------------------------------------------------
score_S <- function(S) {
  hs <- if (opt$hubstrat == 1L) hub_stratify(S, outer(q_pi, q_pi, "+")) else NULL
  lapply(stats::setNames(names(bench), names(bench)), function(nm) {
    b <- bench[[nm]]; p <- b$pairs; sc <- S[cbind(p$i, p$j)]
    list(aupr = ap_fast(sc, p$label),
         jack = vapply(halves[[nm]], function(ix) ap_fast(sc[ix], p$label[ix]), 0),
         coh = regulon_coherence(S, b$reg, b$co, b$operon, genes),
         aupr_hubstrat = if (is.null(hs)) NA_real_ else ap_fast(hs[cbind(p$i, p$j)], p$label)) })
}
run_variant <- function(v) {
  ts <- c(heldout = selv(v, "t_star"), t9 = opt$t_secondary)
  curve <- list(); sc <- list()
  steps <- if (ecoli) seq_len(max(opt$tmax_curve, ts)) else sort(unique(ts))   # K562: t* and t = 9 only
  walk_mass(As[[v]], steps, function(S, t) {
    if (ecoli) for (nm in names(bench)) { p <- bench[[nm]]$pairs
      curve[[length(curve) + 1L]] <<- data.frame(variant = v, truth = nm, t = t, aupr = ap_fast(S[cbind(p$i, p$j)], p$label)) }
    for (d in names(ts)[ts == t]) sc[[d]] <<- score_S(S)
  })
  list(variant = v, curve = do.call(rbind, curve), sc = sc)
}
W <- min(opt$workers, length(variants)); THREADS <- opt$threads
pc <- par_start(mem_gb = opt$mem_gb, workers = W, reserve = reserve_for(W), say = say,
                export = setdiff(ls(globalenv()), c("pc", "cur")))
ev <- foreach(v = variants, .errorhandling = "pass") %dopar% { r <- run_variant(v); par_release(); r }
par_stop(pc)
bad <- vapply(ev, function(x) inherits(x, "error"), NA)
if (any(bad)) stop("variants failed: ", paste(vapply(ev[bad], conditionMessage, ""), collapse = " | "))
names(ev) <- vapply(ev, `[[`, "", "variant")

## ---- outputs ---------------------------------------------------------------------------------------
if (ecoli) {
curve <- do.call(rbind, lapply(ev, `[[`, "curve"))
curve$t_heldout <- curve$t == vapply(curve$variant, function(v) selv(v, "t_star"), 0L)
wcsv(curve, "fulldata_depth_curve_DIAGNOSTIC.csv")
ub <- do.call(rbind, lapply(split(curve, list(curve$variant, curve$truth), drop = TRUE), function(x)
  data.frame(variant = x$variant[1], truth = x$truth[1], t_best_truth_chosen_UPPER_BOUND = x$t[which.max(x$aupr)],
             aupr_upper_bound = max(x$aupr), aupr_at_heldout_t = x$aupr[x$t_heldout][1])))
wcsv(ub, "fulldata_truth_chosen_t_UPPER_BOUND.csv")
}
set.seed(opt$seed + 7L)
met <- list(); jr <- list()
for (dp in c("heldout", "t9")) for (nm in names(bench)) {
  ref <- ev$K10$sc[[dp]][[nm]]
  for (v in variants) {
    e <- ev[[v]]$sc[[dp]][[nm]]; dj <- e$jack - ref$jack; d0 <- e$aupr - ref$aupr; se <- stats::sd(dj)
    ok <- !is.na(e$coh) & !is.na(ref$coh); dC <- e$coh[ok] - ref$coh[ok]
    bt <- if (v != "K10" && length(dC) >= 2) replicate(opt$boot, stats::median(dC[sample(length(dC), replace = TRUE)])) else NA
    met[[length(met) + 1L]] <- data.frame(dataset = opt$dataset, variant = v, depth = dp,
      t = if (dp == "t9") opt$t_secondary else selv(v, "t_star"), truth = nm, aupr = e$aupr,
      base_rate = mean(bench[[nm]]$pairs$label), coherence = stats::median(e$coh, na.rm = TRUE),
      d_vs_K10 = d0, se = se, lo = d0 - 1.96 * se, hi = d0 + 1.96 * se,
      p = if (v == "K10" || !(se > 0)) NA_real_ else 2 * stats::pnorm(-abs(d0) / se),
      p_beyond_margin = if (v == "K10" || !(se > 0)) NA_real_ else stats::pnorm(-(abs(d0) - margin) / se),
      coh_median_d_vs_K10 = if (length(dC)) stats::median(dC) else NA_real_,
      coh_d_lo = if (all(is.na(bt))) NA_real_ else stats::quantile(bt, 0.025, names = FALSE),
      coh_d_hi = if (all(is.na(bt))) NA_real_ else stats::quantile(bt, 0.975, names = FALSE),
      n_regulons_coh = length(dC), aupr_hubstrat_DIAGNOSTIC = e$aupr_hubstrat, stringsAsFactors = FALSE)
    jr[[length(jr) + 1L]] <- data.frame(variant = v, depth = dp, truth = nm, rep = seq_along(e$jack), aupr = e$jack)
  }
}
met <- do.call(rbind, met)
met <- merge(met, perp[, c("variant", "frac_below_3", "frac_above_30", "median")], by = "variant", sort = FALSE)
names(met)[names(met) == "median"] <- "row_perplexity_median"
met <- met[order(met$depth, met$truth, match(met$variant, variants)), ]
wcsv(met, "metrics.csv"); wcsv(do.call(rbind, jr), "jack_replicates.csv")
ep <- endpoints[endpoints$dataset == opt$dataset, ]
pri <- met[met$depth == "heldout" & met$truth %in% ep$truth & met$variant %in% primary_variants, ]
pri <- data.frame(dataset = opt$dataset, variant = pri$variant, endpoint = ep$endpoint[match(pri$truth, ep$truth)],
                  t_star = pri$t, d = pri$d_vs_K10, se = pri$se, p = pri$p, quick = opt$quick)
wcsv(pri, "primary.csv")
for (dp in c("heldout", "t9")) for (nm in names(bench)) {
  x <- met[met$depth == dp & met$truth == nm, ]
  say(sprintf("== %s %s at %s (base rate %.4f)", opt$dataset, nm, if (dp == "t9") "t = 9" else "held-out t*", x$base_rate[1]))
  for (k in seq_len(nrow(x))) say(sprintf("   %-4s t %2d AUPR %.4f  vs K10 %+.4f [%+.4f, %+.4f] p %.3g | coherence %.3f (%+.3f [%+.3f, %+.3f])",
    x$variant[k], x$t[k], x$aupr[k], x$d_vs_K10[k], x$lo[k], x$hi[k], x$p[k], x$coherence[k],
    x$coh_median_d_vs_K10[k], x$coh_d_lo[k], x$coh_d_hi[k]))
}
combine(opt$out, say)
say("done; outputs in ", normalizePath(dout))
