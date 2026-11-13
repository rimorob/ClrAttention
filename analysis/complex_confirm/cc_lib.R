# Library for the exploratory estimation of docs/plans/complex_confirm_design.md (rev 2, FIXED 2026-11-13).
# a3r (CLR on log P^t) minus attention (and CLR) on yeast protein complexes in DREAM5 network 4.
# Nothing here is tuned; every constant below is fixed by the design. Sourced by analysis/complex_confirm/cc_run.R
# and by the two gate-a scripts. Source from the repo root.

suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/dream5_lib.R")            # d5_* (net reader, replicate groups, series, refit halves, packing, hashes); regulondb.R (ap_fast); regulons.R
source("analysis/dream5_score_lib.R")      # d5_edges_from_ids, d5_module_setup (Y-M regulons of the DREAM5 gold standard; structure only)
source("analysis/robust_z/robust_z_lib.R") # rz_operator, rz_stationary, rz_hub_stratify
source("analysis/explore_stack/stack_lib.R") # sk_a3r, sk_reach, sk_diag (unchanged since d8fdb58)

## ---- constants fixed by the design ------------------------------------------------------------------------------
CC <- list(
  out = "results/complex_confirm", sel = "docs/complex_confirm_selections", design = "docs/plans/complex_confirm_design.md",
  jack_n = 500L, jack_seed = 20261114L, cjack_n = 500L, cjack_seed = 20261115L, t = 9L,
  exp_n = 100L, group_seed_base = 20263000L, series_seed_base = 20264000L,
  order = c("main", "main_xsub_excl", "main_ribo_joined", "noribo", "nomerge", "size2", "main_paralog_incl",
            "small_3to6", "mid_7to19", "large_ge20", "dream5_ym"),
  ref_zero = 0, ref_min = 0.005, ref_ecoli = 0.027,
  gold = "data/dream5/truth_extracted/Evaluation scripts/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network4.tsv",
  libs = c("analysis/explore_stack/stack_lib.R", "analysis/explore_ideas47/ideas47_lib.R", "analysis/robust_z/robust_z_lib.R"),
  libs_commit = "d8fdb58",
  code = c("analysis/complex_confirm/cc_lib.R", "analysis/complex_confirm/cc_run.R", "analysis/complex_confirm/gate_a_k562.R",
           "analysis/complex_confirm/gate_a_m3d.R", "analysis/explore_stack/stack_lib.R", "analysis/explore_ideas47/ideas47_lib.R",
           "analysis/robust_z/robust_z_lib.R", "analysis/dream5_lib.R", "analysis/dream5_score_lib.R", "analysis/regulondb.R", "analysis/regulons.R"))
stopifnot(D5$refit_seed_base == CC$group_seed_base, D5$series_seed_base == CC$series_seed_base, D5$t == CC$t)
CC$main_names <- CC$order

none_map <- stats::setNames(character(0), character(0))

cc_say <- function(log, ...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)); cat(m, "\n")
  if (!is.null(log)) cat(m, "\n", file = log, append = TRUE)
}
cc_wcsv <- function(x, path) { tmp <- paste0(path, ".tmp", Sys.getpid()); utils::write.csv(x, tmp, row.names = FALSE); file.rename(tmp, path); invisible(path) }
cc_wtsv <- function(x, path) { tmp <- paste0(path, ".tmp", Sys.getpid()); utils::write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, path); invisible(path) }
cc_vmhwm_gb <- function() { l <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  v <- sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", grep("^VmHWM:", l, value = TRUE)); if (length(v)) as.numeric(v) / 1048576 else NA_real_ }
cc_vmrss_gb <- function() { l <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  v <- sub("^VmRSS:[[:space:]]*([0-9]+).*", "\\1", grep("^VmRSS:", l, value = TRUE)); if (length(v)) as.numeric(v) / 1048576 else NA_real_ }

## ---- options ---------------------------------------------------------------------------------------------------
cc_options <- function(args, defaults) {
  opt <- defaults; i <- 1L
  while (i <= length(args)) {
    k <- sub("^--", "", args[i]); if (!k %in% names(opt)) stop("unknown option --", k)
    v <- args[i + 1L]; opt[[k]] <- if (is.character(defaults[[k]])) v else if (is.integer(defaults[[k]])) as.integer(v) else as.numeric(v); i <- i + 2L
  }
  opt
}
# the constants the production run must use (the truth/expjack/summarize guards refuse anything else; synthetic smoke tests are exempt)
cc_check_fixed <- function(opt) {
  if (opt$synthetic == 1L) return(invisible(TRUE))
  bad <- c(opt$quick != 0L, opt$jack != CC$jack_n, opt$t != CC$t, opt$seed != CC$jack_seed, opt$cjack != CC$cjack_n, opt$cseed != CC$cjack_seed, opt$exp_n != CC$exp_n)
  if (any(bad)) stop("production modes refuse non-preregistered settings (--quick, --jack, --t, --seed, --cjack, --cseed, --exp_n, --synthetic)")
  invisible(TRUE)
}

## ---- git guards ------------------------------------------------------------------------------------------------
cc_git_rc <- function(...) { x <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE)); s <- attr(x, "status"); list(out = as.character(x), rc = if (is.null(s)) 0L else s) }
cc_clean_committed_pushed <- function(path, what) {
  if (!file.exists(path)) stop("guard: ", what, " not found: ", path)
  if (cc_git_rc("ls-files", "--error-unmatch", shQuote(path))$rc != 0L) stop("guard: ", what, " is not tracked by git (commit it): ", path)
  st <- cc_git_rc("status", "--porcelain", "--", shQuote(path))$out
  if (length(st) && any(nzchar(st))) stop("guard: ", what, " has uncommitted changes: ", path)
  up <- cc_git_rc("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
  if (up$rc != 0L) stop("guard: no upstream branch, cannot confirm that ", what, " is pushed")
  unp <- cc_git_rc("log", "--format=%h", "@{u}..HEAD", "--", shQuote(path))$out
  if (length(unp) && any(nzchar(unp))) stop("guard: ", what, " has commits that are not pushed (", paste(unp, collapse = ","), ")")
  invisible(TRUE)
}
cc_design_guard <- function() {
  cc_clean_committed_pushed(CC$design, "design file")
  if (!any(grepl("^Status: FIXED", readLines(CC$design, n = 10)))) stop("guard: the design file does not carry a 'Status: FIXED' line")
  invisible(TRUE)
}
cc_code_guard <- function() { for (f in CC$code) cc_clean_committed_pushed(f, "code file"); invisible(TRUE) }
cc_truth_hash_guard <- function() {
  h <- readLines(file.path(CC$sel, "yeast_complexes.sha256")); want <- strsplit(h, "[[:space:]]+")
  for (w in want) { f <- file.path(CC$sel, w[2]); got <- digest::digest(file = f, algo = "sha256")
    if (!identical(tolower(w[1]), got)) stop("guard: truth file SHA-256 mismatch: ", w[2], " ", got, " != ", w[1]) }
  for (f in c("yeast_complexes.sha256", "yeast_complexes_members.tsv", "yeast_complexes_excluded_pairs.tsv", "yeast_complexes_counts.csv", "yeast_complexes_source_log.csv"))
    cc_clean_committed_pushed(file.path(CC$sel, f), "frozen truth file")
  invisible(vapply(want, function(w) w[1], ""))
}
cc_cache_files <- function(root) file.path(root, "net4/labelfree/cache", c("cc_clr.scores", "cc_attention.scores", "cc_a3r.scores", "cc_mi_rowmean.scores"))
# the 6010 caches: manifest committed (tracked, clean, pushed), every file present with exactly the manifest SHA-256
cc_cache_guard <- function(root = CC$out) {
  mf <- file.path(root, "net4/labelfree/score_manifest.tsv"); cc_clean_committed_pushed(mf, "6010 score manifest")
  m <- utils::read.delim(mf, stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
  for (f in cc_cache_files(root)) { r <- m[m$file == f, ]; if (nrow(r) != 1L) stop("guard: the 6010 manifest has no row for ", f)
    if (!identical(d5_sha256(f), r$sha256)) stop("guard: cache changed since the manifest: ", f) }
  invisible(m)
}
cc_lib_identity <- function() {
  rows <- lapply(CC$libs, function(f) { cur <- digest::digest(file = f, algo = "sha256")
    # byte-identical test through git itself (no trailing-newline ambiguity): the blob hash of the working file equals the blob in the commit
    wb <- system2("git", c("hash-object", shQuote(f)), stdout = TRUE); rb <- system2("git", c("rev-parse", shQuote(paste0(CC$libs_commit, ":", f))), stdout = TRUE)
    data.frame(file = f, sha256 = cur, blob_working = wb, blob_commit = rb, identical = identical(wb, rb), stringsAsFactors = FALSE) })
  do.call(rbind, rows)
}

## ---- net4 and the three readouts ------------------------------------------------------------------------------------
cc_load_net <- function(quick = 0L, input = "replicate") {
  nd <- d5_read_net(4L, input = input)
  stopifnot(nrow(nd$X) == 5950L || quick > 0L, !anyDuplicated(nd$gene_names), length(nd$gene_names) == nrow(nd$X))
  if (quick > 0L) { set.seed(77L); idx <- sort(sample(nrow(nd$X), quick)); nd$X <- nd$X[idx, , drop = FALSE]; nd$genes <- nd$genes[idx]
    nd$gene_names <- nd$gene_names[idx]; nd$tf_idx <- integer(0); nd$tf_ids <- character(0); set.seed(NULL) }
  nd
}
# CLR, attention and a3r from ONE fit (exactly dream5_lib.R::d5_m_clr_attention for CLR and attention: same operator, same walk, same
# symmetrisation; a3r = sk_a3r on the very P^9 that attention is read from). base = TRUE keeps the full-data objects.
cc_fit_three <- function(X, workers = 8L, base = FALSE) {
  G <- nrow(X); alpha <- D5$alpha; tt <- D5$t
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = workers)$calibrate(method = "normal", combine = "stouffer")
  out <- list(clr = d5_clean_sym(f$clr_scores))
  if (base) { M <- f$mi; out$mi_rm <- (rowSums(M) - diag(M)) / (G - 1); rm(M) }
  f$build_operator(alpha = alpha, softmax_keff = D5$keff, softmax_cap = as.integer(D5$cap))
  if (base) out$A <- f$operator
  P <- Matrix::Diagonal(G, 1 - alpha) + alpha * Matrix::Matrix(f$operator, sparse = TRUE)
  Pt <- as.matrix(P); for (s in seq_len(tt - 1L)) Pt <- as.matrix(Pt %*% P)
  f$release()
  out$attention <- d5_clean_sym((Pt + t(Pt)) / 2)
  out$a3r <- d5_clean_sym(sk_a3r(Pt))
  if (base) { out$Pt <- Pt; out$P <- P }
  rm(Pt, P); invisible(gc(verbose = FALSE)); out
}
cc_pack_hash <- function(S) d5_hash_scores(d5_pack(S))

## ---- truth (frozen member list; structure from the design) ----------------------------------------------------------------
cc_read_truth_files <- function() {
  list(mem = utils::read.delim(file.path(CC$sel, "yeast_complexes_members.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = ""),
       excl = utils::read.delim(file.path(CC$sel, "yeast_complexes_excluded_pairs.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = ""),
       counts = utils::read.csv(file.path(CC$sel, "yeast_complexes_counts.csv"), stringsAsFactors = FALSE),
       src = utils::read.csv(file.path(CC$sel, "yeast_complexes_source_log.csv"), stringsAsFactors = FALSE))
}
# synthetic truth for smoke tests: random partitions of the (non-decoy) genes with the real complex-size distribution of each variant
cc_synth_files <- function(nd, tf, seed = 5L) {
  pool <- nd$gene_names[!grepl("^decoy", nd$gene_names)]; mem <- list()
  set.seed(seed)
  for (v in unique(tf$mem$variant)) { x <- tf$mem[tf$mem$variant == v, ]; sz <- as.integer(table(x$complex))
    sz <- pmax(2L, pmin(sz, 40L)); sz <- sz[cumsum(sz) <= length(pool) * 0.8]; if (!length(sz)) sz <- 3L
    ord <- sample(pool); cs <- cumsum(sz); st <- c(1L, head(cs, -1L) + 1L)
    mem[[v]] <- do.call(rbind, lapply(seq_along(sz), function(k) data.frame(variant = v, dataset = "dream5_net4", complex = paste0("SYN-", k), orf = ord[st[k]:cs[k]], stringsAsFactors = FALSE))) }
  mem <- do.call(rbind, mem)
  # keep the four ribosome-subunit style variants consistent enough for the builder: exclusions are a few random pairs
  pr <- function(n, scope) { a <- sample(pool, n); b <- sample(pool, n); data.frame(scope = scope, orf_a = pmin(a, b), orf_b = pmax(a, b), reason = "synthetic", stringsAsFactors = FALSE) }
  excl <- rbind(pr(10, "all"), pr(20, "main_xsub_excl"))
  set.seed(NULL); list(mem = mem, excl = excl)
}
cc_excl_keys <- function(excl, scopes, nd) {
  G <- nrow(nd$X); e <- excl[excl$scope %in% scopes, ]
  a <- match(e$orf_a, nd$gene_names); b <- match(e$orf_b, nd$gene_names); ok <- !is.na(a) & !is.na(b)
  unique(pmin(a, b)[ok] * (G + 1) + pmax(a, b)[ok])
}
cc_key <- function(i, j, G) as.numeric(i) * (G + 1) + as.numeric(j)
# complexes (named integer-index vectors) of one variant
cc_sets <- function(mem, variant, nd) {
  x <- mem[mem$variant == variant, ]; if (!nrow(x)) stop("variant not in the member list: ", variant)
  rows <- match(x$orf, nd$gene_names); if (anyNA(rows)) stop("member ORFs outside the net4 universe in variant ", variant)
  lapply(split(rows, x$complex), function(r) sort(unique(r)))
}
# pair universe of a set of complexes: all pairs of genes in at least one complex, minus excluded pairs; label = share a complex
cc_pairs <- function(sets, nd, excl_k = NULL) {
  G <- nrow(nd$X)
  reg <- list(members = sets, class = stats::setNames(rep("complex", length(sets)), names(sets)))
  p <- comembership_pairs(reg, none_map, nd$gene_names)
  k <- cc_key(p$i, p$j, G); rm_ <- if (length(excl_k)) k %in% excl_k else rep(FALSE, length(k))
  list(i = p$i[!rm_], j = p$j[!rm_], label = p$label[!rm_], genes_in = p$genes_in, n_removed = sum(rm_))
}
# union-find id of paralog sets (for the per-complex AP table, which takes an 'operon' map); names = ORFs
cc_paralog_map <- function(excl, nd) {
  e <- excl[excl$scope == "all", ]; orf <- unique(c(e$orf_a, e$orf_b)); id <- stats::setNames(seq_along(orf), orf)
  repeat { ch <- FALSE; for (r in seq_len(nrow(e))) { a <- id[[e$orf_a[r]]]; b <- id[[e$orf_b[r]]]; if (a != b) { id[id == max(a, b)] <- min(a, b); ch <- TRUE } }; if (!ch) break }
  stats::setNames(paste0("para", id), names(id))
}
cc_stratum <- function(main_sets, main_pairs, lo, hi, nd, excl_k) {
  sets_in <- main_sets[lengths(main_sets) >= lo & lengths(main_sets) <= hi]; G <- nrow(nd$X)
  ps <- cc_pairs(sets_in, nd, excl_k); pos_k <- cc_key(ps$i, ps$j, G)[ps$label == 1L]
  lab <- as.integer(cc_key(main_pairs$i, main_pairs$j, G) %in% pos_k)
  keep <- !(main_pairs$label == 1L & lab == 0L)      # co-members only of complexes outside the stratum leave the universe
  list(i = main_pairs$i[keep], j = main_pairs$j[keep], label = lab[keep], genes_in = main_pairs$genes_in, n_removed = 0L, sets = sets_in)
}
cc_ym_pairs <- function(nd, tf = NULL, synthetic = FALSE, seed = 9L) {
  if (synthetic) { set.seed(seed); pool <- which(!grepl("^decoy", nd$gene_names)); sets <- lapply(1:15, function(k) sort(sample(pool, sample(5:40, 1L)))); names(sets) <- paste0("TF", 1:15); set.seed(NULL)
    p <- cc_pairs(sets, nd, NULL); p$sets <- sets; return(p) }
  tab <- utils::read.delim(CC$gold, header = FALSE, stringsAsFactors = FALSE)
  edges <- d5_edges_from_ids(tab, nd); b <- d5_module_setup(edges, nd, 5L, 500L)
  list(i = b$pairs$i, j = b$pairs$j, label = b$pairs$label, genes_in = b$pairs$genes_in, n_removed = 0L, sets = b$reg$members)
}
# all eleven truths of the design's table order (CC$order)
cc_build_truths <- function(nd, tf, synthetic = FALSE, ym = TRUE) {
  mem <- tf$mem; excl <- tf$excl; G <- nrow(nd$X)
  k_all <- cc_excl_keys(excl, "all", nd); k_x <- cc_excl_keys(excl, c("all", "main_xsub_excl"), nd)
  S <- function(v) cc_sets(mem, v, nd)
  T <- list()
  mk <- function(sets, k, role, note = "") { p <- cc_pairs(sets, nd, k); p$sets <- sets; p$role <- role; p$note <- note; p }
  T$main <- mk(S("main"), k_all, "primary")
  T$main_xsub_excl <- mk(S("main"), k_x, "sensitivity", "cross-subunit pairs of both ribosomes removed from the universe")
  T$main_ribo_joined <- mk(S("main_ribo_joined"), k_all, "sensitivity", "40S+60S and 37S+54S joined")
  T$noribo <- mk(S("noribo"), k_all, "sensitivity", "pre-declared two-complex cytoplasmic ribosome list removed")
  T$nomerge <- mk(S("nomerge"), k_all, "sensitivity")
  T$size2 <- mk(S("size2"), k_all, "sensitivity", "minimum size 2")
  T$main_paralog_incl <- mk(S("main"), NULL, "sensitivity", "paralog pairs kept")
  ms <- T$main$sets
  for (st in list(c("small_3to6", 3, 6), c("mid_7to19", 7, 19), c("large_ge20", 20, 100000))) {
    if (!any(lengths(ms) >= as.numeric(st[2]) & lengths(ms) <= as.numeric(st[3]))) { if (!synthetic) stop("stratum ", st[1], " is empty"); next }
    p <- cc_stratum(ms, T$main, as.numeric(st[2]), as.numeric(st[3]), nd, k_all); p$role <- "stratum"; p$note <- "primary universe; positives restricted to the stratum's complexes"; T[[st[1]]] <- p }
  if (ym) { y <- cc_ym_pairs(nd, synthetic = synthetic); y$role <- "descriptive"; y$note <- "touched truth (DREAM5 Y-M regulons, TFs with 5-500 targets); a3r never scored on it; no-harm description"; T$dream5_ym <- y }
  T[intersect(CC$order, names(T))]
}
cc_truth_struct_hash <- function(p) digest::digest(list(p$i, p$j, p$label), algo = "md5")

## ---- scoring ---------------------------------------------------------------------------------------------------------
# paired delete-half gene jackknife of the 26 Oct test (clr_logpt_test.R::jack_dataset): halves are drawn for all truths in table
# order right after set.seed; a half is floor(n/2) of the truth's genes_in; pairs with both genes kept are scored; same halves for all methods
cc_draw_halves <- function(truths, n, seed) {
  set.seed(seed); h <- lapply(truths, function(b) lapply(seq_len(n), function(r) sample(b$genes_in, floor(length(b$genes_in) / 2)))); set.seed(NULL); h
}
# scores: named list of symmetric matrices -> list(aupr = named vector, jack = matrix n x methods)
cc_score_truth <- function(Slist, p, halves, G, nc) {
  sc <- vapply(Slist, function(S) S[cbind(p$i, p$j)], numeric(length(p$i)))
  if (is.null(dim(sc))) sc <- matrix(sc, ncol = length(Slist), dimnames = list(NULL, names(Slist)))
  aupr <- vapply(seq_len(ncol(sc)), function(m) ap_fast(sc[, m], p$label), 0); names(aupr) <- colnames(sc)
  jack <- NULL
  if (!is.null(halves)) {
    one <- function(keep) { inh <- logical(G); inh[keep] <- TRUE; ix <- which(inh[p$i] & inh[p$j]); lab <- p$label[ix]
      vapply(seq_len(ncol(sc)), function(m) ap_fast(sc[ix, m], lab), 0) }
    jl <- if (nc > 1L) parallel::mclapply(halves, one, mc.cores = nc) else lapply(halves, one)
    if (any(vapply(jl, function(x) inherits(x, "try-error") || length(x) != ncol(sc), NA))) stop("a jackknife worker failed")
    jack <- do.call(rbind, jl); colnames(jack) <- colnames(sc)
  }
  list(aupr = aupr, jack = jack)
}
cc_skew <- function(x) { m <- mean(x); mean((x - m)^3) / mean((x - m)^2)^1.5 }
CONTRASTS <- list(P = c("a3r", "attention"), C1 = c("a3r", "clr"), C2 = c("attention", "clr"))
cc_contrast_names <- c(P = "a3r - attention", C1 = "a3r - CLR", C2 = "attention - CLR")
# interval summary of one paired contrast from point estimates and replicate vectors
cc_interval <- function(d, reps) {
  se <- stats::sd(reps); q <- stats::quantile(reps - mean(reps), c(0.025, 0.975), names = FALSE, type = 7)
  data.frame(d = d, se_gene = se, gene_lo = d - 1.96 * se, gene_hi = d + 1.96 * se, pct_lo = d + q[1], pct_hi = d + q[2], skew = cc_skew(reps), n_rep = length(reps))
}

## ---- fixed sentences (design "Estimands and reporting rule"; no p-value anywhere) ----------------------------------------------------
cc_sentences <- function(L, U) {
  c(S_zero = if (L > 0) "The total interval excludes no difference: a3r exceeded attention on yeast complexes (preregistered exploratory estimation; not a confirmatory test)."
             else if (U < 0) "The total interval lies below no difference: a3r was lower than attention on yeast complexes."
             else "The total interval includes no difference.",
    S_min = if (L > CC$ref_min) "The whole interval lies above the smallest effect of interest (+0.005)."
            else if (U < CC$ref_min) "Effects at or above the smallest effect of interest (+0.005) are excluded."
            else "The interval includes the smallest effect of interest (+0.005).",
    S_ecoli = if (U < CC$ref_ecoli) "Effects as large as the E. coli estimate (+0.027) are excluded."
              else if (L > CC$ref_ecoli) "The whole interval lies above the E. coli estimate (+0.027)."
              else "The interval includes the E. coli estimate (+0.027).",
    S_fixed = "Exploratory estimation, not a confirmatory test, a replication or a validation. a3r is not adopted; this study cannot change that by itself (PI decisions 12 and 13 Nov).")
}
CC_NOTE_CLR <- "a3r − CLR has a smaller SE and may exclude 0 where a3r − attention does not; this is never a basis for 'a3r beats CLR' (E. coli S3)."
