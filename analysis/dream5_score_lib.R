# DREAM5 scorers, jackknives, verdict and guards (docs/plans/dream5_design.md, FIXED 2026-11-04).
# Everything here is defined by the design except the TRUTH-STAGE items (cap ties, TF-TF pairs, pair universe), which are
# explicit switches that have NO default (d5_edge_switches()); the lead resolves them at the truth stage by a dated
# change-log entry. Functions that touch a truth refuse to run through d5_truth_guard() unless the score-file checksums
# are committed and unchanged. Tests use decoy truths only (d5_decoy_truth(): random edges, never real truth).
#
# Source from the repo root:  source("analysis/dream5_score_lib.R")   (loads analysis/dream5_lib.R)

if (!exists("D5")) source("analysis/dream5_lib.R")

## ---- truth objects ---------------------------------------------------------------------------------
# edges: data.frame(tf = row index into the TF list, gene = gene index). truth matrix: logical |TF| x G.
d5_truth_matrix <- function(edges, nTF, G) {
  T <- matrix(FALSE, nTF, G); T[cbind(edges$tf, edges$gene)] <- TRUE; T
}
# gene-ID edge table (two columns: regulator ID, target ID) -> index edges; unknown IDs are dropped and counted
d5_edges_from_ids <- function(tab, nd) {
  tf <- match(tab[[1]], nd$tf_ids); ge <- match(tab[[2]], nd$genes)
  ok <- !is.na(tf) & !is.na(ge)
  e <- unique(data.frame(tf = tf[ok], gene = ge[ok])); attr(e, "n_dropped") <- sum(!ok); e
}
# Decoy truth for smoke tests: a made-up network with TF out-degrees from a heavy-tailed distribution and targets drawn
# with unequal in-degree propensities; nothing is derived from any real truth.
d5_decoy_truth <- function(nTF, G, tf_idx, seed, median_out = 8, sd_log = 1.1, max_out = NULL) {
  .with_seed(seed, {
    max_out <- if (is.null(max_out)) max(2L, as.integer(G / 4)) else max_out
    out <- pmin(pmax(round(exp(stats::rnorm(nTF, log(median_out), sd_log))), 1L), max_out)
    prop <- exp(stats::rnorm(G, 0, 1))
    e <- do.call(rbind, lapply(seq_len(nTF), function(a) {
      cand <- setdiff(seq_len(G), tf_idx[a])
      data.frame(tf = a, gene = cand[sample.int(length(cand), min(out[a], length(cand)), prob = prop[cand])])
    }))
    e
  })
}

## ---- module scorer (Y-M and analogues) ------------------------------------------------------------
# Regulons = target sets of each TF with min_size..max_size targets in the expression gene set (make_regulons limits;
# targets only; no operon exclusion). Pairs: all pairs among genes in at least one regulon; positive if they share one.
d5_module_setup <- function(edges, nd, min_size = 5L, max_size = 500L) {
  genes <- nd$genes
  tab <- data.frame(reg_class = "TF", regulator = nd$tf_ids[edges$tf], bnumber = genes[edges$gene], confidence = "C",
                    stringsAsFactors = FALSE)
  reg <- make_regulons(tab, genes, conf_keep = "C", min_size = min_size, max_size = max_size)
  none <- stats::setNames(character(0), character(0))
  pairs <- comembership_pairs(reg, none, genes)
  list(reg = reg, pairs = pairs, n_regulons = length(reg$members), n_pos = sum(pairs$label), n_pairs = length(pairs$label),
       min_size = min_size, max_size = max_size)
}
# Tie-aware (expected over random order within tied blocks) average precision, exact (McSherry and Najork 2008):
# for a block of m items holding p positives, after k0 items and tp0 positives, the item at block position i is positive
# with probability p/m and then has E[precision] = (tp0 + 1 + (i-1)(p-1)/(m-1)) / (k0 + i).
ap_expected <- function(score, label) {
  score[is.na(score)] <- -Inf
  o <- order(score, decreasing = TRUE, method = "radix"); s <- score[o]; l <- label[o]; N <- length(s)
  newb <- c(TRUE, s[-1L] != s[-N]); bid <- cumsum(newb); bstart <- which(newb); m <- diff(c(bstart, N + 1L))
  cl <- cumsum(l); cl0 <- c(0, cl); tp0 <- cl0[bstart]; p <- cl0[bstart + m] - tp0; k0 <- bstart - 1L
  i <- seq_len(N) - k0[bid]; mm <- m[bid]; pp <- p[bid]
  num <- tp0[bid] + 1 + ifelse(mm > 1L, (i - 1) * (pp - 1) / pmax(mm - 1, 1), 0)
  sum((pp / mm) * num / (k0[bid] + i)) / sum(l)
}
# pair scores from a packed upper-triangle vector (column-major, i < j)
d5_pair_scores <- function(v, pairs) v[d5_pair_index(pairs$i, pairs$j)]
# primary Y-M statistic (ties counted at the end of their block, ap_fast); sensitivity: expected-tie AP
d5_module_aupr <- function(v, bench, expected_ties = FALSE, ix = NULL) {
  sc <- d5_pair_scores(v, bench$pairs); lab <- bench$pairs$label
  if (!is.null(ix)) { sc <- sc[ix]; lab <- lab[ix] }
  if (expected_ties) ap_expected(sc, lab) else ap_fast(sc, lab)
}
# Evaluability gate for Y-M (structure only, no score): at least 20 regulons with >= 5 targets and at least 1,000
# positive co-membership pairs in the expression gene set. Regulon count follows the design's wording (>= 5 targets, no
# upper limit); positive pairs are those of the Y-M regulon set (5 to 500 targets). Both counts are reported.
d5_ym_gate <- function(edges, nd) {
  b5 <- d5_module_setup(edges, nd, 5L, .Machine$integer.max); b <- d5_module_setup(edges, nd, 5L, 500L)
  data.frame(n_regulons_ge5 = b5$n_regulons, n_regulons_5_500 = b$n_regulons, n_pos_pairs = b$n_pos,
             n_pos_pairs_ge5_noupper = b5$n_pos, evaluable = b5$n_regulons >= 20L && b$n_pos >= 1000L)
}
# structure-only counts (needs truth; no score)
d5_structure_counts <- function(edges, nd) {
  per_tf <- table(factor(edges$tf, levels = seq_along(nd$tf_ids)))
  g <- d5_ym_gate(edges, nd); b10 <- d5_module_setup(edges, nd, 10L, 500L)
  cbind(data.frame(n_edges = nrow(edges), n_tfs_with_edges = sum(per_tf > 0), median_edges_per_tf = stats::median(per_tf[per_tf > 0]),
                   max_edges_per_tf = max(per_tf), n_tf_tf_edges = sum(edges$gene %in% nd$tf_idx)),
        g, data.frame(n_regulons_ge10 = b10$n_regulons, n_pos_pairs_ge10 = b10$n_pos))
}

## ---- edge scorer (Y-E, E-arm, in silico) -----------------------------------------------------------
# TRUTH-STAGE switches. None has a default; d5_edge_metrics() stops if one that matters is unset. Documented options:
#   ties_at_cap   how the tied block that straddles rank 100,000 is cut:
#                   keep_tied   keep every entry tied with the cap-th score (list longer than the cap)
#                   drop_tied   drop the whole straddling block (list shorter than the cap)
#                   list_order  keep exactly `cap` entries, ties broken by candidate order (column-major cell index)
#   ties_in_list  ties inside the list: end_of_block (positives ranked last in their block, as ap_fast) | expected
#   tftf          TF-TF pairs: keep_both (as submitted) | collapse_max (one orientation per pair, the higher score, ties to the
#                 lower TF-list row; a pair is positive if either orientation is) | drop (TF-TF pairs not scored)
#   universe      all (every candidate pair) | gs (restrict to TFs and genes given in `universe = list(tf =, gene =)`)
#   cap_order     (only if universe == gs) cap_then_universe (cap on all candidates, then restrict) | universe_then_cap
#   jackknife_cap (only with kept subsets) reapply_cap (restrict to kept pairs, then cap) | scaled_cap (same, cap x kept
#                 fraction) | restrict_after_cap (global cap, then restrict to kept pairs)
# Fixed by the paper (not switches): pairs not in the submitted list are randomly ordered at the end (expected AP exactly;
# AUROC with one tied block of mid-ranks).
d5_edge_switches <- function(ties_at_cap = NA, ties_in_list = NA, tftf = NA, universe = NA, cap_order = NA, jackknife_cap = NA)
  list(ties_at_cap = ties_at_cap, ties_in_list = ties_in_list, tftf = tftf, universe = universe, cap_order = cap_order,
       jackknife_cap = jackknife_cap)
.d5_sw_levels <- list(ties_at_cap = c("keep_tied", "drop_tied", "list_order"), ties_in_list = c("end_of_block", "expected"),
                      tftf = c("keep_both", "collapse_max", "drop"), universe = c("all", "gs"),
                      cap_order = c("cap_then_universe", "universe_then_cap"),
                      jackknife_cap = c("reapply_cap", "scaled_cap", "restrict_after_cap"))
.d5_sw_check <- function(sw, need) {
  for (k in need) if (is.null(sw[[k]]) || length(sw[[k]]) != 1L || is.na(sw[[k]]) || !sw[[k]] %in% .d5_sw_levels[[k]])
    stop("TRUTH-STAGE switch '", k, "' is not set to one of: ", paste(.d5_sw_levels[[k]], collapse = ", "),
         " (resolved by the lead at the truth stage; see d5_edge_switches)")
}

# Average precision and AUROC of a ranked list followed by a randomly ordered tail.
#   s, y: score and 0/1 label of the listed pairs; n_pos, N: positives and pairs in the whole universe
d5_list_tail_metrics <- function(s, y, n_pos, N, ties_in_list) {
  K <- length(s); if (n_pos == 0) return(c(aupr = NA_real_, auroc = NA_real_))
  o <- order(s, decreasing = TRUE, method = "radix"); s <- s[o]; y <- y[o]
  tp_L <- sum(y); n0 <- N - n_pos
  # ---- average precision
  if (K > 0L) {
    newb <- c(TRUE, s[-1L] != s[-K]); bid <- cumsum(newb); bstart <- which(newb); m <- diff(c(bstart, K + 1L))
    cl0 <- c(0, cumsum(y)); tp0 <- cl0[bstart]; p <- cl0[bstart + m] - tp0; k0 <- bstart - 1L
    if (ties_in_list == "end_of_block") {
      part <- sum(p * (tp0 + p) / (k0 + m))               # every positive of a block gets precision at the block end
    } else {
      i <- seq_len(K) - k0[bid]; mm <- m[bid]; pp <- p[bid]
      part <- sum((pp / mm) * (tp0[bid] + 1 + ifelse(mm > 1L, (i - 1) * (pp - 1) / pmax(mm - 1, 1), 0)) / (k0[bid] + i))
    }
  } else { part <- 0; bstart <- integer(0); m <- integer(0) }
  Nt <- N - K; n_rest <- n_pos - tp_L
  tail_part <- if (Nt > 0 && n_rest > 0) { i <- seq_len(Nt)
    sum((n_rest / Nt) * (tp_L + 1 + (i - 1) * (n_rest - 1) / max(Nt - 1, 1)) / (K + i)) } else 0
  aupr <- (part + tail_part) / n_pos
  # ---- AUROC (ascending ranks; the tail is one tied block of the lowest ranks)
  rank_sum <- n_rest * (Nt + 1) / 2
  if (K > 0L) {
    bend <- bstart + m - 1L                                  # positions in descending order: block [bstart, bend]
    rank_sum <- rank_sum + if (ties_in_list == "end_of_block")
      sum(p * (N - bend + 1) + p * (p - 1) / 2)             # positives take the lowest ranks of their block
    else sum(p * (N - (bstart + bend) / 2 + 1))              # mid-ranks
  }
  auroc <- (rank_sum - n_pos * (n_pos + 1) / 2) / (as.numeric(n_pos) * n0)
  c(aupr = aupr, auroc = auroc)
}

# M: |TF| x G score matrix (NA on self pairs); truth: logical |TF| x G; tf_idx: gene index of each TF row.
# keep_tf / keep_gene: kept TF rows and gene columns (truth jackknife); universe: list(tf =, gene =) for universe == "gs".
# Returns aupr, auroc, n_pos, N (pairs in the scored universe), K (list length), cap_used.
d5_edge_metrics <- function(M, truth, tf_idx, sw, cap = 100000L, keep_tf = NULL, keep_gene = NULL, universe = NULL) {
  nTF <- nrow(M); G <- ncol(M)
  .d5_sw_check(sw, c("ties_at_cap", "ties_in_list", "tftf", "universe"))
  if (sw$universe == "gs") { .d5_sw_check(sw, "cap_order"); stopifnot(is.list(universe), !is.null(universe$tf), !is.null(universe$gene)) }
  jk <- !is.null(keep_tf) || !is.null(keep_gene)
  if (jk) .d5_sw_check(sw, "jackknife_cap")
  valid <- !is.na(M)
  Y <- truth
  # TF-TF pairs
  if (sw$tftf == "drop") valid[, tf_idx] <- FALSE
  if (sw$tftf == "collapse_max") {
    A <- M[, tf_idx, drop = FALSE]; v <- valid[, tf_idx, drop = FALSE]; yy <- Y[, tf_idx, drop = FALSE]
    ut <- which(upper.tri(A), arr.ind = TRUE)                 # (a, b) with a < b: orientations a->b and b->a
    ab <- cbind(ut[, 1], ut[, 2]); ba <- cbind(ut[, 2], ut[, 1])
    both <- v[ab] & v[ba]
    drop_ab <- both & (A[ab] < A[ba]); drop_ba <- both & !drop_ab        # ties: keep a->b (lower TF-list row first)
    either <- yy[ab] | yy[ba]
    v[ab[drop_ab, , drop = FALSE]] <- FALSE; v[ba[drop_ba, , drop = FALSE]] <- FALSE
    keepo <- both & drop_ba; yy[ab[keepo, , drop = FALSE]] <- either[keepo]
    keepo <- both & drop_ab; yy[ba[keepo, , drop = FALSE]] <- either[keepo]
    valid[, tf_idx] <- v; Y[, tf_idx] <- yy
  }
  masks <- function(rows, cols) { m <- matrix(FALSE, nTF, G); m[rows, cols] <- TRUE; m }
  pre <- valid; post <- matrix(TRUE, nTF, G); cap_eff <- cap
  if (sw$universe == "gs") {
    mu <- masks(universe$tf, universe$gene)
    if (sw$cap_order == "cap_then_universe") post <- post & mu else pre <- pre & mu
  }
  if (jk) {
    mj <- masks(if (is.null(keep_tf)) seq_len(nTF) else keep_tf, if (is.null(keep_gene)) seq_len(G) else keep_gene)
    if (sw$jackknife_cap == "restrict_after_cap") post <- post & mj else {
      pre <- pre & mj
      if (sw$jackknife_cap == "scaled_cap") cap_eff <- max(1L, as.integer(round(cap * sum(mj) / (nTF * G))))
    }
  }
  # the whole scored universe: every restriction applies
  univ <- valid & (if (sw$universe == "gs") masks(universe$tf, universe$gene) else TRUE) & (if (jk) mj else TRUE)
  cells <- which(pre); sc <- M[cells]
  o <- order(sc, decreasing = TRUE, method = "radix"); cells <- cells[o]; sc <- sc[o]
  if (length(cells) > cap_eff) {
    th <- sc[cap_eff]
    nkeep <- switch(sw$ties_at_cap, keep_tied = sum(sc >= th), drop_tied = sum(sc > th), list_order = cap_eff)
    cells <- cells[seq_len(nkeep)]; sc <- sc[seq_len(nkeep)]
  }
  inpost <- post[cells]; cells <- cells[inpost]; sc <- sc[inpost]
  N <- sum(univ); n_pos <- sum(Y[univ])
  r <- d5_list_tail_metrics(sc, as.integer(Y[cells]), n_pos, N, sw$ties_in_list)
  c(r, n_pos = n_pos, N = N, K = length(cells), cap_used = cap_eff)
}

## ---- score retrieval ------------------------------------------------------------------------------
# TF-row matrix from a packed symmetric vector
d5_edge_matrix_sym <- function(v, nd) {
  nTF <- length(nd$tf_idx); G <- length(nd$genes)
  idx <- d5_pair_index(rep(nd$tf_idx, G), rep(seq_len(G), each = nTF))
  M <- matrix(v[idx], nTF, G); M[d5_selfmask(nd)] <- NA; M
}
# score of a method as the object the scorers need: list(edge = |TF| x G matrix | NULL, pair = packed vector | NULL)
d5_method_scores <- function(net, method, nd, input = nd$input, out = D5$out) {
  x <- d5_read_scores(d5_cache_path(net, method, input, out)); kind <- D5$methods[[method]]; G <- length(nd$genes)
  switch(kind,
    sym = list(edge = d5_edge_matrix_sym(x, nd), pair = x),
    tfrow = list(edge = x, pair = NULL),
    dir = { W <- x; E <- W[nd$tf_idx, , drop = FALSE]; W[is.na(W)] <- 0
            list(edge = E, pair = d5_pack(pmax(W, t(W)))) })      # all-genes tree runs: max symmetrisation (analysis/baselines.R)
}

## ---- truth jackknife and experiment jackknife ---------------------------------------------------------------
# SE_J halves, shared by every method: a random half of the TFs and a random half of the genes, per replicate
d5_jack_halves_edge <- function(seed, nTF, G, nrep = 200L)
  .with_seed(seed, lapply(seq_len(nrep), function(r) list(tf = sort(sample.int(nTF, floor(nTF / 2))), gene = sort(sample.int(G, floor(G / 2))))))
d5_jack_halves_module <- function(seed, G, nrep = 200L)
  .with_seed(seed, lapply(seq_len(nrep), function(r) sort(sample.int(G, floor(G / 2)))))

# Y-E truth jackknife: AUPR of every method on every replicate (rows = replicate, columns = method)
d5_jack_edge <- function(Ms, truth, tf_idx, sw, seed, nrep = 200L, cap = 100000L, universe = NULL) {
  h <- d5_jack_halves_edge(seed, nrow(truth), ncol(truth), nrep)
  out <- vapply(seq_len(nrep), function(r) vapply(Ms, function(M) unname(d5_edge_metrics(M, truth, tf_idx, sw, cap, h[[r]]$tf, h[[r]]$gene, universe)["aupr"]), 0), numeric(length(Ms)))
  t(matrix(out, length(Ms), nrep, dimnames = list(names(Ms), NULL)))
}
# Y-M truth jackknife (the kept genes' pairs; primary AP convention)
d5_jack_module <- function(vs, bench, G, seed, nrep = 200L) {
  h <- d5_jack_halves_module(seed, G, nrep); i <- bench$pairs$i; j <- bench$pairs$j
  scs <- lapply(vs, function(v) d5_pair_scores(v, bench$pairs)); lab <- bench$pairs$label
  out <- vapply(seq_len(nrep), function(r) { inh <- logical(G); inh[h[[r]]] <- TRUE; ix <- which(inh[i] & inh[j])
    vapply(scs, function(sc) ap_fast(sc[ix], lab[ix]), 0) }, numeric(length(vs)))
  t(matrix(out, length(vs), nrep, dimnames = list(names(vs), NULL)))
}

# Effect: d = log2(AUPR_a / AUPR_b); paired differences over replicates -> SD
d5_effect <- function(a, b) log2(a / b)
d5_se_paired <- function(rep_a, rep_b) stats::sd(d5_effect(rep_a, rep_b))
d5_combine_se <- function(se_j, se_b) sqrt(se_j^2 + se_b^2)

## ---- verdict (decision table of the design) --------------------------------------------------------------
D5$delta <- log2(1.10); D5$alpha_side <- 0.025
d5_verdict <- function(d, se, delta = D5$delta, alpha = D5$alpha_side) {
  stopifnot(length(d) == 1L, length(se) == 1L)
  if (!is.finite(d) || !is.finite(se) || se < 0) return(data.frame(d = d, SE = se, p_sup = NA_real_, p_harm = NA_real_, p_noninf_lo = NA_real_, p_noninf_hi = NA_real_, verdict = "Inconclusive"))
  z <- if (se == 0) sign(d) * Inf else d / se
  p_sup <- 1 - stats::pnorm(z); p_harm <- stats::pnorm(z)
  p_lo <- if (se == 0) as.numeric(!(d > -delta)) else 1 - stats::pnorm((d + delta) / se)   # test d > -delta
  p_hi <- if (se == 0) as.numeric(!(d < delta)) else 1 - stats::pnorm((delta - d) / se)    # test d < +delta
  v <- if (p_sup < alpha && d >= delta) "Superior"
       else if (p_harm < alpha && d <= -delta) "Inferior"
       else if (p_lo < alpha && p_hi < alpha) "Equivalent"
       else "Inconclusive"
  data.frame(d = d, SE = se, p_sup = p_sup, p_harm = p_harm, p_noninf_lo = p_lo, p_noninf_hi = p_hi, verdict = v)
}

## ---- guard -----------------------------------------------------------------------------------------------
d5_required <- function(nets = D5$nets, methods = D5$primary, input = "replicate")
  expand.grid(net = nets, method = methods, input = input, stringsAsFactors = FALSE)
.git <- function(repo, ...) suppressWarnings(system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = TRUE))
.git_rc <- function(repo, ...) { x <- .git(repo, ...); s <- attr(x, "status"); list(out = as.character(x), rc = if (is.null(s)) 0L else s) }
.git_clean_committed_pushed <- function(repo, path, what) {
  rel <- path
  if (!file.exists(file.path(repo, rel))) stop("DREAM5 truth guard: ", what, " not found: ", rel)
  if (.git_rc(repo, "ls-files", "--error-unmatch", shQuote(rel))$rc != 0L) stop("DREAM5 truth guard: ", what, " is not tracked by git (commit it): ", rel)
  st <- .git_rc(repo, "status", "--porcelain", "--", shQuote(rel))$out
  if (length(st) && any(nzchar(st))) stop("DREAM5 truth guard: ", what, " has uncommitted changes: ", rel)
  up <- .git_rc(repo, "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
  if (up$rc != 0L) stop("DREAM5 truth guard: no upstream branch, cannot confirm that ", what, " is pushed")
  unp <- .git_rc(repo, "log", "--format=%h", "@{u}..HEAD", "--", shQuote(rel))$out
  if (length(unp) && any(nzchar(unp))) stop("DREAM5 truth guard: ", what, " has commits that are not pushed (", paste(unp, collapse = ","), ")")
  invisible(TRUE)
}
# Refuses (stop) unless: the consolidated score manifest (and the refit manifest, if given) is tracked, unchanged vs HEAD and
# pushed; it lists every required score file; and every listed file exists with exactly the manifest's SHA-256.
# `root` is the directory the manifest's file column is relative to (the repo root in production).
d5_truth_guard <- function(manifest = file.path(D5$out, "SCORE_MANIFEST.tsv"), required = d5_required(), repo = ".",
                           refit_manifest = NULL, refit_required = NULL, root = repo, verbose = TRUE) {
  .git_clean_committed_pushed(repo, manifest, "score manifest")
  m <- utils::read.delim(file.path(repo, manifest), stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
  key <- paste(m$net, m$method, m$input); need <- paste(required$net, required$method, required$input)
  miss <- setdiff(need, key)
  if (length(miss)) stop("DREAM5 truth guard: the manifest lacks required score files: ", paste(miss, collapse = "; "))
  chk <- function(df, what) for (k in seq_len(nrow(df))) {
    p <- file.path(root, df$file[k])
    if (!file.exists(p)) stop("DREAM5 truth guard: ", what, " file missing: ", df$file[k])
    if (!identical(d5_sha256(p), df$sha256[k])) stop("DREAM5 truth guard: ", what, " file changed since the manifest: ", df$file[k])
  }
  chk(m[key %in% need, ], "score")
  if (!is.null(refit_manifest)) {
    .git_clean_committed_pushed(repo, refit_manifest, "refit manifest")
    r <- utils::read.delim(file.path(repo, refit_manifest), stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
    if (!is.null(refit_required)) {
      miss <- setdiff(paste(refit_required$net, refit_required$rep, refit_required$method), paste(r$net, r$rep, r$method))
      if (length(miss)) stop("DREAM5 truth guard: the refit manifest lacks ", length(miss), " replicate/method hashes, e.g. ", miss[1])
    }
  }
  if (verbose) message("DREAM5 truth guard: ", sum(key %in% need), " score files verified against the committed manifest")
  invisible(m)
}
# The design must be FIXED, tracked, unchanged and pushed, and the data log committed and pushed (as the truth stage of
# analysis/dream5_download.py requires).
d5_design_guard <- function(design = "docs/plans/dream5_design.md", data_log = "docs/dream5_selections/download_log.tsv", repo = ".") {
  .git_clean_committed_pushed(repo, design, "design file"); .git_clean_committed_pushed(repo, data_log, "data log")
  if (!any(grepl("^Status: FIXED", readLines(file.path(repo, design), n = 10)))) stop("DREAM5 truth guard: the design file does not carry a 'Status: FIXED' line")
  invisible(TRUE)
}
# verify a recomputed refit replicate against the committed hash (used by the truth stage before scoring it)
d5_refit_verify <- function(rep_scores, b, net, refit_manifest) {
  r <- utils::read.delim(refit_manifest, stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
  for (m in names(rep_scores)) {
    h <- r$sha256[r$net == net & r$rep == b & r$method == m]
    if (length(h) != 1L) stop("refit manifest has no hash for net", net, " rep ", b, " ", m)
    if (!identical(d5_hash_scores(rep_scores[[m]]), h)) stop("refit replicate ", b, " (", m, ") differs from the committed hash; not scoring")
  }
  invisible(TRUE)
}

## ---- the truth-phase pipeline (never run before the truth stage; tested on decoys) ---------------------------------
D5$seeds <- list(`4` = c(YE = 20261203L, YM = 20261204L), `3` = c(YE = 20261213L, YM = 20261214L), `1` = c(YE = 20261223L, YM = 20261224L))
# Point estimates, SE_J (200 replicates, design seeds), SE_B (net4: refits recomputed, hash-verified), effects, verdicts.
# The caller (analysis/dream5_truth_run.R) runs d5_truth_guard() first; this function does not read any file except the
# score caches, the refit manifest and what it is handed. `refit_scorer(b)` returns the refit upper-triangle vectors.
d5_truth_pipeline <- function(nd, edges, sw, methods = D5$primary, out = D5$out, nrep = 200L, refit_reps = integer(0),
                              refit_scorer = NULL, refit_verify = NULL, universe = NULL, cap = 100000L,
                              contrast = c("attention", "clr"), seeds = D5$seeds[[as.character(nd$net)]], delta = D5$delta, verbose = TRUE) {
  net <- nd$net; G <- length(nd$genes); nTF <- length(nd$tf_idx); tf_idx <- nd$tf_idx
  say <- function(...) if (verbose) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  truth <- d5_truth_matrix(edges, nTF, G)
  bench <- d5_module_setup(edges, nd, 5L, 500L); bench10 <- d5_module_setup(edges, nd, 10L, 500L)
  gate <- d5_ym_gate(edges, nd)
  sc <- stats::setNames(lapply(methods, function(m) d5_method_scores(net, m, nd, out = out)), methods)
  Ms <- lapply(sc, `[[`, "edge"); vs <- Filter(Negate(is.null), lapply(sc, `[[`, "pair"))
  # ---- point estimates
  pe <- do.call(rbind, lapply(methods, function(m) {
    ye <- d5_edge_metrics(Ms[[m]], truth, tf_idx, sw, cap, universe = universe)
    data.frame(net = net, method = m, endpoint = "YE", aupr = ye[["aupr"]], auroc = ye[["auroc"]], n_pos = ye[["n_pos"]], N = ye[["N"]], aupr_sens = NA_real_, aupr_min10 = NA_real_)
  }))
  if (nrow(gate) && gate$evaluable) pe <- rbind(pe, do.call(rbind, lapply(names(vs), function(m)
    data.frame(net = net, method = m, endpoint = "YM", aupr = d5_module_aupr(vs[[m]], bench), auroc = NA_real_, n_pos = bench$n_pos,
               N = bench$n_pairs, aupr_sens = d5_module_aupr(vs[[m]], bench, expected_ties = TRUE),
               aupr_min10 = if (bench10$n_pos > 0) d5_module_aupr(vs[[m]], bench10) else NA_real_))))
  # ---- SE_J
  say("truth jackknife, ", nrep, " replicates")
  jE <- d5_jack_edge(Ms, truth, tf_idx, sw, seeds[["YE"]], nrep, cap, universe)
  jM <- if (gate$evaluable) d5_jack_module(vs, bench, G, seeds[["YM"]], nrep) else NULL
  # ---- SE_B (refits recomputed and verified against the committed hashes; the truth and gene universe are fixed)
  bE <- bM <- NULL
  if (length(refit_reps)) {
    say("experiment jackknife, ", length(refit_reps), " replicates")
    rr <- lapply(refit_reps, function(b) {
      v <- refit_scorer(b); if (!is.null(refit_verify)) refit_verify(v, b)
      e <- vapply(names(v), function(m) unname(d5_edge_metrics(d5_edge_matrix_sym(v[[m]], nd), truth, tf_idx, sw, cap, universe = universe)["aupr"]), 0)
      mm <- if (gate$evaluable) vapply(names(v), function(m) d5_module_aupr(v[[m]], bench), 0) else NULL
      list(e = e, m = mm) })
    bE <- do.call(rbind, lapply(rr, `[[`, "e")); bM <- if (gate$evaluable) do.call(rbind, lapply(rr, `[[`, "m")) else NULL
  }
  # ---- effects and verdicts
  a <- contrast[1]; b <- contrast[2]; rows <- list()
  for (ep in c("YE", "YM")) {
    if (ep == "YM" && !gate$evaluable) { rows[[ep]] <- data.frame(net = net, endpoint = ep, verdict = "not evaluable"); next }
    J <- if (ep == "YE") jE else jM; Bm <- if (ep == "YE") bE else bM; p <- pe[pe$endpoint == ep, ]
    A <- p$aupr[p$method == a]; B <- p$aupr[p$method == b]; d <- d5_effect(A, B)
    se_j <- d5_se_paired(J[, a], J[, b]); se_b <- if (!is.null(Bm)) d5_se_paired(Bm[, a], Bm[, b]) else NA_real_
    se <- if (is.na(se_b)) NA_real_ else d5_combine_se(se_j, se_b)
    v <- if (net == 4L && !is.na(se)) d5_verdict(d, se, delta) else NULL
    rows[[ep]] <- data.frame(net = net, endpoint = ep, readout = c(YE = "official directed edge AUPR", YM = "regulon co-membership AUPR")[[ep]],
                             contrast = paste(a, "vs", b), aupr_a = A, aupr_b = B, ratio = A / B, difference = A - B, d = d, delta = delta,
                             SE_J = se_j, SE_B = se_b, SE = se, ci_lo = d - 1.96 * se, ci_hi = d + 1.96 * se,
                             ci_J_lo = d - 1.96 * se_j, ci_J_hi = d + 1.96 * se_j,
                             p_sup = if (is.null(v)) NA else v$p_sup, p_harm = if (is.null(v)) NA else v$p_harm,
                             p_noninf_lo = if (is.null(v)) NA else v$p_noninf_lo, p_noninf_hi = if (is.null(v)) NA else v$p_noninf_hi,
                             verdict = if (is.null(v)) "descriptive (no verdict: needs net4 and SE_B)" else v$verdict, stringsAsFactors = FALSE)
  }
  list(point = pe, gate = gate, jack_edge = jE, jack_module = jM, refit_edge = bE, refit_module = bM, rows = rows)
}

# write verdict.csv (the two net4 rows) and the supporting tables
d5_write_results <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  vv <- do.call(rbind, lapply(res$rows, function(r) { miss <- setdiff(names(res$rows$YE), names(r)); for (m in miss) r[[m]] <- NA; r[, names(res$rows$YE)] }))
  utils::write.csv(vv, file.path(dir, "verdict.csv"), row.names = FALSE)
  utils::write.csv(res$point, file.path(dir, "point_estimates.csv"), row.names = FALSE)
  utils::write.csv(res$gate, file.path(dir, "ym_gate.csv"), row.names = FALSE)
  utils::write.csv(cbind(rep = seq_len(nrow(res$jack_edge)), as.data.frame(res$jack_edge)), file.path(dir, "jack_edge.csv"), row.names = FALSE)
  if (!is.null(res$jack_module)) utils::write.csv(cbind(rep = seq_len(nrow(res$jack_module)), as.data.frame(res$jack_module)), file.path(dir, "jack_module.csv"), row.names = FALSE)
  if (!is.null(res$refit_edge)) utils::write.csv(cbind(rep = seq_len(nrow(res$refit_edge)), as.data.frame(res$refit_edge)), file.path(dir, "refit_edge.csv"), row.names = FALSE)
  if (!is.null(res$refit_module)) utils::write.csv(cbind(rep = seq_len(nrow(res$refit_module)), as.data.frame(res$refit_module)), file.path(dir, "refit_module.csv"), row.names = FALSE)
  invisible(vv)
}
