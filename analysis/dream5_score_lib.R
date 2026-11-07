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
# gene-ID edge table -> index edges. The table has exactly two columns (regulator ID, target ID: every row is a positive
# edge) or three, the third a 0/1 label (gold-standard layout): then only the label-1 rows are kept and the label-0 rows are
# dropped and counted; any other third-column value is an error, so a 0-labelled row is never read as a positive. Unknown IDs
# are dropped and counted (n_dropped); self-edges (TF -> itself) are dropped and counted (n_self); repeats are collapsed
# (n_dup). Counts are attributes n_dropped, n_label0, n_self, n_dup of the result.
d5_edges_from_ids <- function(tab, nd) {
  nc <- ncol(tab)
  if (!nc %in% c(2L, 3L)) stop("truth edge table must have exactly two columns (regulator, target) or three (regulator, target, 0/1 label); it has ", nc)
  n_label0 <- 0L
  if (nc == 3L) {
    lab <- suppressWarnings(as.numeric(as.character(tab[[3]])))
    if (anyNA(lab) || !all(lab %in% c(0, 1))) stop("the third column of the truth edge table must hold 0/1 labels only")
    n_label0 <- sum(lab == 0); tab <- tab[lab == 1, 1:2, drop = FALSE]
  }
  tf <- match(tab[[1]], nd$tf_ids); ge <- match(tab[[2]], nd$genes)
  ok <- !is.na(tf) & !is.na(ge)
  selfe <- ok & !is.na(ge) & ge == nd$tf_idx[ifelse(is.na(tf), 1L, tf)]
  keep <- ok & !selfe
  e <- data.frame(tf = tf[keep], gene = ge[keep]); n_dup <- sum(duplicated(e)); e <- unique(e)
  attr(e, "n_dropped") <- sum(!ok); attr(e, "n_label0") <- n_label0; attr(e, "n_self") <- sum(selfe); attr(e, "n_dup") <- n_dup; e
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
d5_ap_expected <- function(score, label) {
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
  if (expected_ties) d5_ap_expected(sc, lab) else ap_fast(sc, lab)
}
# Evaluability gate for Y-M (structure only, no score): at least 20 regulons AND at least 1,000 positive co-membership
# pairs, BOTH counted on the regulon set that the readout scores (5 to 500 targets in the expression gene set; design
# change log 2026-11-05, item 3). The unrestricted count (>= 5 targets, no upper limit) is reported as an extra field,
# n_regulons_ge5 (and n_pos_pairs_ge5_noupper); it does not decide.
d5_ym_gate <- function(edges, nd) {
  b5 <- d5_module_setup(edges, nd, 5L, .Machine$integer.max); b <- d5_module_setup(edges, nd, 5L, 500L)
  data.frame(n_regulons_5_500 = b$n_regulons, n_pos_pairs = b$n_pos,
             n_regulons_ge5 = b5$n_regulons, n_pos_pairs_ge5_noupper = b5$n_pos,
             evaluable = b$n_regulons >= 20L && b$n_pos >= 1000L)
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
#                   drop_tied   drop the whole straddling block (list shorter than the cap); a block straddles only if the entry after
#                               the cap-th has the same score, so a unique cap-th entry is kept (fixed 2026-11-06: it used to be dropped)
#                   list_order  keep exactly `cap` entries, ties broken by candidate order (column-major cell index)
#   ties_in_list  ties inside the list: end_of_block (ap_fast: every positive of a tied block is credited the precision at the END
#                 of the block, (TP before + p)/(rank of the block end); this is higher than ranking the block's positives last,
#                 which would give them precisions (TP before + j)/(rank of block end - p + j), j = 1..p) |
#                 expected (exact mean over random orders inside each tied block)
#   tftf          TF-TF pairs: keep_both (as submitted) | collapse_max (one orientation per pair, the higher score, ties to the
#                 lower TF-list row; a pair is positive if either orientation is) | drop (TF-TF pairs not scored)
#   universe      all (every candidate pair) | gs (restrict to TFs and genes given in `universe = list(tf =, gene =)`)
#   cap_order     (only if universe == gs) cap_then_universe (cap on all candidates, then restrict) | universe_then_cap
#   jackknife_cap (only with kept subsets) reapply_cap (restrict to kept pairs, then cap) | scaled_cap (same, cap x kept
#                 fraction) | restrict_after_cap (global cap, then restrict to kept pairs)
# Fixed by the paper (not switches): pairs not in the submitted list are randomly ordered at the end (expected AP exactly;
# AUROC with one tied block of mid-ranks).
# Y-E decisions use the OFFICIAL trapezoid metrics (d5_official_expected, ties_in_list = expected); the step-wise AP above is a
# labelled sensitivity (D5$ye_readouts). The Y-E primary definition is D5$ye_primary (design change log 2026-11-06).
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

## ---- the official Y-E metrics (DREAM5 Challenge 4 scoring scripts; design change log 2026-11-06, decision 6) --------------
# The official scorer ranks the (capped, universe-restricted) submitted list in file order, appends every unlisted universe
# pair with the expected positive fraction p = (P - TPL)/(T - L) (one fractional step per position), and computes
#   AUROC = trapz(FPR, TPR)                    (DREAM5_Challenge4_Evaluation.m:76)
#   AUPR  = trapz(REC, PREC) / (1 - 1/P)       (:104)
# over the T points (k = 1..T; the point (0, 0) is NOT added, so the first point only starts the first trapezoid).
# These are NOT step-wise average precision (d5_list_tail_metrics, kept as a labelled sensitivity, not for Y-E decisions).
d5_trapz <- function(x, y) { n <- length(x); sum((x[-1L] - x[-n]) * (y[-1L] + y[-n]) / 2) }
# Literal transcription of the official script for a list that is already in rank order (disc: 0/1 discovery vector of the
# listed pairs, P positives and T pairs in the universe). Reference implementation; also used by the file-level scorer.
d5_official_literal <- function(disc, P, T) {
  L <- length(disc); p <- if (L < T) (P - sum(disc)) / (T - L) else 0
  pos <- c(disc, rep(p, T - L)); neg <- c(1 - disc, rep(1 - p, T - L))
  TPk <- cumsum(pos); FPk <- cumsum(neg); K <- seq_len(T)
  TPR <- TPk / P; FPR <- FPk / (T - P); PREC <- TPk / K
  c(aupr = d5_trapz(TPR, PREC) / (1 - 1 / P), auroc = d5_trapz(FPR, TPR))
}
# Order-free expectation of the official statistics when scores are tied (ties_in_list = "expected": the tied entries of a
# block are ranked in uniformly random order; the official scorer's file order is arbitrary for ties).
#   s, y: score and 0/1 label of the listed pairs; P positives and T pairs in the universe.  With no ties this equals
#   d5_official_literal() on the sorted list to rounding error.
# Derivation. Write pos_k, neg_k = 1 - pos_k, TP_k = sum_{j<=k} pos_j. Then, since FPR_{k+1} - FPR_k = neg_{k+1}/N and
# REC_{k+1} - REC_k = pos_{k+1}/P,
#   AUROC = 1/(N P) * sum_{k=2..T} neg_k (TP_{k-1} + pos_k/2)
#   AUPR  = 1/(P (1 - 1/P)) * sum_{k=2..T} pos_k (TP_{k-1}/(k-1) + TP_k/k) / 2,      TP_k = TP_{k-1} + pos_k.
# Both are, term by term, a 0/1 label at position k times an affine function of the labels before k, so the expectation over
# a uniformly random order inside a tied block (m items, p positives, TP0 positives and k0 items before the block, position
# i = 1..m, k = k0 + i) needs only P(y_i = 1) = p/m and E[sum_{j<i} y_j | y_i = 1] = (i-1)(p-1)/(m-1); the order of the other
# blocks does not matter because a block's total positives p, and hence TP0 of every later block, are order-free.
#   AUROC term: neg_k pos_k = 0 for 0/1 labels, so the term is neg_k TP_{k-1}; E[(1-y_i)(TP0 + sum_{j<i} y_j)]
#       = TP0 (1 - p/m) + (i-1) p (m-p) / (m (m-1)); summed over i = 1..m this is (m - p)(TP0 + p/2): the block is a straight
#       segment from (FPR, TPR) at its start to its end (the trapezoid over mid-ranks), independent of k.
#   AUPR term:  E[y_i (TP_{k-1}/(k-1) + (TP_{k-1} + 1)/k)] / 2 = (p/m) (a_i/(k-1) + (a_i + 1)/k) / 2 with
#       a_i = TP0 + (i-1)(p-1)/(m-1) (a_i = TP0 when m = 1); the position k = 1 contributes nothing (no trapezoid ends there).
# The unlisted tail (T - K positions, each pos = p_t = (P - TPL)/(T - K), TP_k = TPL + (k-K) p_t) is deterministic and is
# summed exactly (k >= 2 only: if K = 0 the term k = 1 is dropped, as in the official script).
# Validated by exact enumeration of every within-block order and by Monte Carlo (tests/dream5/test-official.R).
d5_official_expected <- function(s, y, P, T) {
  K <- length(s); NA2 <- c(aupr = NA_real_, auroc = NA_real_)
  if (P < 2 || T - P < 1) return(NA2)
  if (K > T) stop("list longer than the universe")
  o <- order(s, decreasing = TRUE, method = "radix"); s <- s[o]; y <- y[o]
  TPL <- sum(y); Nt <- T - K; pt <- if (Nt > 0) (P - TPL) / Nt else 0
  a_roc <- 0; a_pr <- 0
  if (K > 0L) {
    newb <- c(TRUE, s[-1L] != s[-K]); bid <- cumsum(newb); bstart <- which(newb); m <- diff(c(bstart, K + 1L))
    cl0 <- c(0, cumsum(y)); tp0 <- cl0[bstart]; p <- cl0[bstart + m] - tp0; k0 <- bstart - 1L
    a_roc <- sum((m - p) * (tp0 + p / 2))
    i <- seq_len(K) - k0[bid]; mm <- m[bid]; pp <- p[bid]; k <- seq_len(K)
    a <- tp0[bid] + ifelse(mm > 1L, (i - 1) * (pp - 1) / pmax(mm - 1, 1), 0)
    term <- (pp / mm) * (a / pmax(k - 1, 1) + (a + 1) / k) / 2; term[k == 1L] <- 0
    a_pr <- sum(term)
  }
  if (Nt > 0) {
    a_roc <- a_roc + (1 - pt) * Nt * (TPL + pt * Nt / 2) - (if (K == 0L) (1 - pt) * pt / 2 else 0)
    k <- seq.int(K + 1L, T); tpm <- TPL + (k - K - 1) * pt; tpk <- tpm + pt     # TP_{k-1}, TP_k
    keep <- k >= 2L
    a_pr <- a_pr + sum((pt * (tpm / pmax(k - 1, 1) + tpk / k) / 2)[keep])
  }
  c(aupr = a_pr / P / (1 - 1 / P), auroc = a_roc / ((T - P) * P))
}
# The official p-value: P(X >= x) = sum((X >= x) * Y * dx), dx = X(2) - X(1) (probability.m), from the precomputed null densities
# INPUT/probability_densities/Network{n}_{AUPR,AUROC}.mat, converted once to TSV (analysis/dream5_official_tools.py convert-densities).
# p = 0 exactly above the support top; it is returned as 0 (never floored), and -log10(0) is Inf.
D5$densities_dir <- "data/dream5/truth_extracted_tools/densities"
.d5_dens_cache <- new.env()
d5_read_density <- function(net, metric = c("AUPR", "AUROC"), dir = D5$densities_dir) {
  metric <- match.arg(metric); f <- file.path(dir, sprintf("Network%d_%s.tsv", as.integer(net), metric))
  if (!file.exists(f)) stop("density table not found: ", f, " (convert the .mat files once with analysis/dream5_official_tools.py convert-densities)")
  if (is.null(.d5_dens_cache[[f]])) { d <- utils::read.delim(f, colClasses = "numeric"); .d5_dens_cache[[f]] <- list(X = d$X, Y = d$Y) }
  .d5_dens_cache[[f]]
}
d5_density_p <- function(dens, x) { dx <- dens$X[2] - dens$X[1]; vapply(x, function(v) if (is.na(v)) NA_real_ else sum((dens$X >= v) * dens$Y * dx), 0) }
d5_official_p <- function(net, aupr, auroc, dir = D5$densities_dir)
  c(P_AUPR = d5_density_p(d5_read_density(net, "AUPR", dir), aupr), P_AUROC = d5_density_p(d5_read_density(net, "AUROC", dir), auroc))
# go_all.m:66-68: aupr_score = -mean(log10(P_AUPR)), auroc_score = -mean(log10(P_AUROC)) over the networks, overall = their mean.
# An exact p of 0 gives Inf (not capped).
d5_official_overall <- function(p_aupr, p_auroc) {
  a <- -mean(log10(p_aupr)); b <- -mean(log10(p_auroc)); c(overall = mean(c(a, b)), aupr_score = a, auroc_score = b)
}
# Overall score per method from the per-network official_scores.csv files written by the truth driver (one run per network;
# nets 1, 3, 4): go_all.m's overall = mean(-mean log10 P_AUPR, -mean log10 P_AUROC). Descriptive; Inf if any p is exactly 0.
d5_official_overall_from_dirs <- function(dirs) {
  tb <- do.call(rbind, lapply(dirs, function(d) utils::read.csv(file.path(d, "official_scores.csv"))))
  if (!setequal(tb$net, D5$nets)) stop("official_scores.csv files must cover networks ", paste(D5$nets, collapse = ", "), "; got ", paste(sort(unique(tb$net)), collapse = ", "))
  do.call(rbind, lapply(split(tb, tb$method), function(x) { o <- d5_official_overall(x$P_AUPR, x$P_AUROC); data.frame(method = x$method[1], t(o), row.names = NULL) }))
}
# Gold-derived pair universe (edgelist_2_D5C4_gs.m:3-4): regulators = distinct TFs of the gold edges, targets = every gene in any
# gold edge, either column (so regulators are targets too). `edges` is the index edge table of d5_edges_from_ids(); returns the
# TF rows and gene columns of the net, in TF-list / gene-list order, and the IDs.
d5_gold_universe <- function(edges, nd) {
  tf <- sort(unique(edges$tf)); gene <- sort(unique(c(nd$tf_idx[edges$tf], edges$gene)))
  list(tf = tf, gene = gene, tf_ids = nd$tf_ids[tf], gene_ids = nd$genes[gene])
}
# write the universe as the two ID files the driver reads (one ID per line); returns the paths
d5_write_universe <- function(u, prefix) {
  f <- paste0(prefix, c("_tf.txt", "_gene.txt")); writeLines(u$tf_ids, f[1]); writeLines(u$gene_ids, f[2]); f
}
# check that universe files list exactly the derived universe (sets of IDs)
d5_check_universe_files <- function(u, tf_file, gene_file) {
  a <- readLines(tf_file); b <- readLines(gene_file)
  if (!setequal(a, u$tf_ids) || anyDuplicated(a)) stop("universe TF file ", tf_file, " does not match the TFs of the gold edges")
  if (!setequal(b, u$gene_ids) || anyDuplicated(b)) stop("universe gene file ", gene_file, " does not match the genes of the gold edges")
  TRUE
}
# Y-E readout: the official metrics (primary, decisions) or the step-wise average precision (sensitivity only)
D5$ye_readouts <- c(official = "official Y-E: trapezoid AUPR/(1-1/P) and AUROC with expected-fraction tail (design 2026-11-06, decision 6)",
                    stepwise = "SENSITIVITY ONLY (not for Y-E decisions): step-wise average precision")
.d5_pick <- function(r, metric = c("stepwise", "official"), what = c("aupr", "auroc")) {
  metric <- match.arg(metric); what <- match.arg(what)
  v <- unname(if (metric == "official") r[[paste0(what, "_official")]] else r[[what]])
  if (metric == "official" && is.na(v) && !is.na(r[["n_pos"]]) && r[["n_pos"]] >= 2) stop("official readout needs ties_in_list = 'expected'")
  v
}

# M: |TF| x G score matrix (NA on self pairs); truth: logical |TF| x G; tf_idx: gene index of each TF row.
# keep_tf / keep_gene: kept TF rows and gene columns (truth jackknife); universe: list(tf =, gene =) for universe == "gs".
# Returns aupr, auroc (step-wise AP; sensitivity), n_pos, N (pairs in the scored universe), K (list length), cap_used,
# aupr_official, auroc_official (the official trapezoid metrics; only with ties_in_list == "expected", else NA) and
# n_tied_at_cap (size of the tied block that straddles the cap boundary: the entries sharing the score of the cap-th entry when
# the next entry has the same score; 0 if the cap cuts between two different scores or nothing is cut).
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
  n_tied <- 0L
  if (length(cells) > cap_eff) {
    th <- sc[cap_eff]; straddle <- sc[cap_eff + 1L] == th     # a block straddles the cap only if the next entry ties with the cap-th
    if (straddle) n_tied <- sum(sc == th)
    nkeep <- switch(sw$ties_at_cap, keep_tied = sum(sc >= th), drop_tied = if (straddle) sum(sc > th) else cap_eff, list_order = cap_eff)
    cells <- cells[seq_len(nkeep)]; sc <- sc[seq_len(nkeep)]
  }
  inpost <- post[cells]; cells <- cells[inpost]; sc <- sc[inpost]
  N <- sum(univ); n_pos <- sum(Y[univ]); yl <- as.integer(Y[cells])
  r <- d5_list_tail_metrics(sc, yl, n_pos, N, sw$ties_in_list)
  ro <- if (sw$ties_in_list == "expected") d5_official_expected(sc, yl, n_pos, N) else c(aupr = NA_real_, auroc = NA_real_)
  c(r, n_pos = n_pos, N = N, K = length(cells), cap_used = cap_eff, aupr_official = ro[["aupr"]], auroc_official = ro[["auroc"]], n_tied_at_cap = n_tied)
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
# metric: "official" (Y-E decisions; needs ties_in_list = expected) or "stepwise" (sensitivity)
d5_jack_edge <- function(Ms, truth, tf_idx, sw, seed, nrep = 200L, cap = 100000L, universe = NULL, metric = "stepwise") {
  h <- d5_jack_halves_edge(seed, nrow(truth), ncol(truth), nrep)
  out <- vapply(seq_len(nrep), function(r) vapply(Ms, function(M) .d5_pick(d5_edge_metrics(M, truth, tf_idx, sw, cap, h[[r]]$tf, h[[r]]$gene, universe), metric), 0), numeric(length(Ms)))
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
# The scoring code is part of what is guarded: the score library, the truth driver, the label-free library and the two
# project files the scorers call (ap_fast, make_regulons / comembership_pairs) must be tracked, unchanged and pushed too.
D5$code_files <- c("analysis/dream5_lib.R", "analysis/dream5_score_lib.R", "analysis/dream5_truth_run.R", "analysis/regulondb.R", "analysis/regulons.R")
d5_truth_guard <- function(manifest = file.path(D5$out, "SCORE_MANIFEST.tsv"), required = d5_required(), repo = ".",
                           refit_manifest = NULL, refit_required = NULL, series_manifest = NULL, series_required = NULL, root = repo, verbose = TRUE,
                           code_files = D5$code_files) {
  .git_clean_committed_pushed(repo, manifest, "score manifest")
  for (f in code_files) .git_clean_committed_pushed(repo, f, "scoring code file")
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
  # the supplementary series-level refit manifest (net4): same rule (tracked, unchanged, pushed, complete)
  if (!is.null(series_manifest)) {
    .git_clean_committed_pushed(repo, series_manifest, "series-level refit manifest")
    r <- utils::read.delim(file.path(repo, series_manifest), stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
    if (!is.null(series_required)) {
      miss <- setdiff(paste(series_required$net, series_required$rep, series_required$method), paste(r$net, r$rep, r$method))
      if (length(miss)) stop("DREAM5 truth guard: the series-level refit manifest lacks ", length(miss), " replicate/method hashes, e.g. ", miss[1])
    }
  }
  if (verbose) message("DREAM5 truth guard: ", sum(key %in% need), " score files verified against the committed manifest")
  invisible(m)
}
# Arguments of d5_truth_guard() for the net4 refit manifests. With series_n == 0 the supplementary series-level manifest is neither
# required nor read (the series part is supplementary and cannot block the primary verdict).
d5_truth_guard_args <- function(refit_n, series_n, root_out = D5$out) {
  req <- function(n) expand.grid(net = 4L, rep = seq_len(n), method = D5$refit_methods, stringsAsFactors = FALSE)
  a <- list(refit_manifest = file.path(d5_refit_dir(4L, root_out), "REFIT_MANIFEST.tsv"), refit_required = req(refit_n))
  if (series_n > 0L) { a$series_manifest <- file.path(d5_refit_dir(4L, root_out, unit = "series"), "REFIT_SERIES_MANIFEST.tsv"); a$series_required <- req(series_n) }
  a
}
# The design must be FIXED, tracked, unchanged and pushed, and the data log committed and pushed (as the truth stage of
# analysis/dream5_download.py requires).
d5_design_guard <- function(design = "docs/plans/dream5_design.md", data_log = "docs/dream5_selections/download_log.tsv", repo = ".") {
  .git_clean_committed_pushed(repo, design, "design file"); .git_clean_committed_pushed(repo, data_log, "data log")
  if (!any(grepl("^Status: FIXED", readLines(file.path(repo, design), n = 10)))) stop("DREAM5 truth guard: the design file does not carry a 'Status: FIXED' line")
  invisible(TRUE)
}
# verify a recomputed refit replicate against the committed hash (used by the truth stage before scoring it); `refit_manifest` is
# the group-level REFIT_MANIFEST.tsv or the series-level REFIT_SERIES_MANIFEST.tsv (same layout)
d5_refit_verify <- function(rep_scores, b, net, refit_manifest) {
  r <- utils::read.delim(refit_manifest, stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
  for (m in names(rep_scores)) {
    h <- r$sha256[r$net == net & r$rep == b & r$method == m]
    if (length(h) != 1L) stop("refit manifest has no hash for net", net, " rep ", b, " ", m)
    if (!identical(d5_hash_scores(rep_scores[[m]]), h)) stop("refit replicate ", b, " (", m, ") differs from the committed hash; not scoring")
  }
  invisible(TRUE)
}

## ---- truth driver: command line -> validated options and switches (before any guard, any net or truth file is touched) --------
# Every switch and numeric option is parsed and validated here; a switch or file that the chosen configuration would ignore is an
# error (universe == all with --universe_tf/--universe_gene/--cap_order; --workers for a net other than 4). Nothing is read
# beyond the command line, except that the universe files must exist.
d5_truth_options <- function(args) {
  opt <- list(net = "", truth_edges = "", out_dir = "", ties_at_cap = NA_character_, ties_in_list = NA_character_, tftf = NA_character_,
              universe = NA_character_, universe_tf = "", universe_gene = "", cap_order = NA_character_, jackknife_cap = NA_character_,
              structure_only = 0L, workers = 16L, nrep = 200L, cap = 100000L, refit_n = 100L, series_n = 100L,
              run_kind = "sensitivity", readout = NA_character_)
  given <- character(0); i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i]); if (!grepl("^--", args[i]) || !key %in% names(opt)) stop("unknown argument: ", args[i])
    if (i + 1L > length(args) || grepl("^--", args[i + 1L])) stop("argument ", args[i], " has no value")
    v <- args[i + 1L]
    if (is.integer(opt[[key]])) { x <- suppressWarnings(as.numeric(v)); if (is.na(x) || x != round(x)) stop("argument --", key, " must be an integer, got '", v, "'"); v <- as.integer(x) }
    opt[[key]] <- v; given <- c(given, key); i <- i + 2L
  }
  if (anyDuplicated(given)) stop("argument given twice: --", given[duplicated(given)][1])
  net <- suppressWarnings(as.integer(opt$net))
  if (is.na(net) || !net %in% D5$nets) stop("--net must be one of ", paste(D5$nets, collapse = ", "))
  if (!nzchar(opt$truth_edges) || !nzchar(opt$out_dir)) stop("--truth_edges and --out_dir are required")
  # the switches: all of them are validated first, whatever the mode
  sw <- d5_edge_switches(ties_at_cap = opt$ties_at_cap, ties_in_list = opt$ties_in_list, tftf = opt$tftf, universe = opt$universe,
                         cap_order = opt$cap_order, jackknife_cap = opt$jackknife_cap)
  .d5_sw_check(sw, c("ties_at_cap", "ties_in_list", "tftf", "universe", "jackknife_cap"))   # the truth jackknife always runs
  # run kind: the Y-E PRIMARY run is the design's fixed definition (change log 2026-11-06); a sensitivity run may use any switch
  # values but is labelled so, and its readout is chosen explicitly (default: step-wise AP, never used for Y-E decisions)
  if (!opt$run_kind %in% c("primary", "sensitivity")) stop("--run_kind must be primary or sensitivity")
  primary <- opt$run_kind == "primary"
  if (primary) {
    bad <- names(D5$ye_primary)[!vapply(names(D5$ye_primary), function(k) identical(unname(sw[[k]]), D5$ye_primary[[k]]), NA)]
    if (length(bad)) stop("--run_kind primary (the Y-E primary run) allows only the decisions of the design change log 2026-11-06 (the official scorer's definition): ",
                          paste(sprintf("--%s %s", names(D5$ye_primary), unlist(D5$ye_primary)), collapse = " "), "; got a different value for: ",
                          paste(sprintf("%s = %s", bad, vapply(bad, function(k) as.character(sw[[k]]), "")), collapse = ", "),
                          ". Use --run_kind sensitivity for any other value (such a run is labelled and is not a Y-E decision).")
    if (!is.na(opt$readout) && opt$readout != "official") stop("--run_kind primary uses the official readout only (design 2026-11-06, decision 6); --readout ", opt$readout, " is a sensitivity")
    opt$readout <- "official"
  } else {
    if (is.na(opt$readout)) opt$readout <- "stepwise"
    if (!opt$readout %in% names(D5$ye_readouts)) stop("--readout must be one of ", paste(names(D5$ye_readouts), collapse = ", "))
    if (opt$readout == "official" && sw$ties_in_list != "expected") stop("--readout official needs --ties_in_list expected")
  }
  has_uf <- nzchar(opt$universe_tf) || nzchar(opt$universe_gene)
  if (sw$universe == "gs") {
    .d5_sw_check(sw, "cap_order")
    if (primary) {   # the universe is derived from the gold edges (d5_gold_universe); files are optional and then verified against it
      if (has_uf && !(nzchar(opt$universe_tf) && nzchar(opt$universe_gene))) stop("give both --universe_tf and --universe_gene (they are checked against the gold-derived universe) or neither")
    } else if (!nzchar(opt$universe_tf) || !nzchar(opt$universe_gene)) stop("--universe gs needs both --universe_tf and --universe_gene")
    for (f in c(opt$universe_tf, opt$universe_gene)) if (nzchar(f) && !file.exists(f)) stop("universe file not found: ", f)
  } else {
    if (has_uf) stop("--universe_tf/--universe_gene were given but --universe is 'all': they would be ignored")
    if (!is.na(sw$cap_order)) stop("--cap_order was given but --universe is 'all': it would be ignored")
  }
  if (!opt$structure_only %in% 0:1) stop("--structure_only must be 0 or 1")
  if (net != 4L && "workers" %in% given) stop("--workers was given but net ", net, " has no refits: it would be ignored")
  if (opt$cap < 1L) stop("--cap must be positive")
  if (opt$nrep < 2L) stop("--nrep must be at least 2 (an SD needs two replicates)")
  if (opt$workers < 1L) stop("--workers must be positive")
  if (opt$refit_n < 0L || opt$refit_n > D5$n_refit) stop("--refit_n must be in 0..", D5$n_refit)
  if (opt$series_n < 0L || opt$series_n > D5$n_series_refit) stop("--series_n must be in 0..", D5$n_series_refit)
  if (net == 4L && opt$structure_only == 0L && opt$refit_n < 2L) stop("net 4 needs --refit_n of at least 2 (SE_B)")
  opt$net <- net
  list(opt = opt, sw = sw)
}
# universe files (one gene ID per line) -> index vectors; an ID that is not in the network is an error (never silently NA)
d5_universe_from_files <- function(opt, nd) {
  if (!nzchar(opt$universe_tf)) return(NULL)
  tf <- readLines(opt$universe_tf); ge <- readLines(opt$universe_gene)
  if (any(!tf %in% nd$tf_ids)) stop("universe TF file has ", sum(!tf %in% nd$tf_ids), " IDs that are not in the TF list, e.g. ", tf[!tf %in% nd$tf_ids][1])
  if (any(!ge %in% nd$genes)) stop("universe gene file has ", sum(!ge %in% nd$genes), " IDs that are not in the gene list, e.g. ", ge[!ge %in% nd$genes][1])
  list(tf = match(tf, nd$tf_ids), gene = match(ge, nd$genes))
}

## ---- the truth-phase pipeline (never run before the truth stage; tested on decoys) ---------------------------------
# The design's Y-E primary definition (change log 2026-11-06, decisions 1-5); the readout is the official one (decision 6)
D5$ye_primary <- list(universe = "gs", cap_order = "cap_then_universe", tftf = "keep_both", ties_at_cap = "drop_tied",
                      ties_in_list = "expected", jackknife_cap = "restrict_after_cap")
D5$seeds <- list(`4` = c(YE = 20261203L, YM = 20261204L), `3` = c(YE = 20261213L, YM = 20261214L), `1` = c(YE = 20261223L, YM = 20261224L))
# Point estimates, SE_J (200 replicates, design seeds), SE_B (net4: refits recomputed, hash-verified), effects, verdicts.
# The caller (analysis/dream5_truth_run.R) runs d5_truth_guard() first; this function does not read any file except the
# score caches, the refit manifest and what it is handed. `refit_scorer(b)` returns the refit upper-triangle vectors.
d5_truth_pipeline <- function(nd, edges, sw, methods = D5$primary, out = D5$out, nrep = 200L, refit_reps = integer(0),
                              refit_scorer = NULL, refit_verify = NULL, series_reps = integer(0), series_scorer = NULL, series_verify = NULL,
                              universe = NULL, cap = 100000L, metric = "stepwise", densities = NULL,
                              contrast = c("attention", "clr"), seeds = D5$seeds[[as.character(nd$net)]], delta = D5$delta, verbose = TRUE) {
  if (!metric %in% names(D5$ye_readouts)) stop("metric must be one of ", paste(names(D5$ye_readouts), collapse = ", "))
  if (metric == "official" && !identical(sw$ties_in_list, "expected")) stop("the official Y-E readout needs ties_in_list = 'expected'")
  net <- nd$net; G <- length(nd$genes); nTF <- length(nd$tf_idx); tf_idx <- nd$tf_idx
  say <- function(...) if (verbose) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  truth <- d5_truth_matrix(edges, nTF, G)
  bench <- d5_module_setup(edges, nd, 5L, 500L); bench10 <- d5_module_setup(edges, nd, 10L, 500L)
  gate <- d5_ym_gate(edges, nd)
  sc <- stats::setNames(lapply(methods, function(m) d5_method_scores(net, m, nd, out = out)), methods)
  Ms <- lapply(sc, `[[`, "edge"); vs <- Filter(Negate(is.null), lapply(sc, `[[`, "pair"))
  # ---- point estimates
  yes <- stats::setNames(lapply(methods, function(m) d5_edge_metrics(Ms[[m]], truth, tf_idx, sw, cap, universe = universe)), methods)
  if (metric == "official" && !is.null(universe)) {   # the official T = regulators x targets minus self pairs; a method with NA cells would shrink the universe
    Tg <- length(universe$tf) * length(universe$gene) - sum(nd$tf_idx[universe$tf] %in% universe$gene)
    bad <- methods[vapply(yes, function(r) r[["N"]] != Tg, NA)]
    if (length(bad)) stop("official readout: the scored universe differs from the gold universe size T = ", Tg, " for: ", paste(bad, collapse = ", "), " (non-self cells with NA scores?)")
  }
  pe <- do.call(rbind, lapply(methods, function(m) { ye <- yes[[m]]
    data.frame(net = net, method = m, endpoint = "YE", aupr = .d5_pick(ye, metric, "aupr"), auroc = .d5_pick(ye, metric, "auroc"), n_pos = ye[["n_pos"]], N = ye[["N"]],
               aupr_sens = NA_real_, aupr_min10 = NA_real_, readout = metric, aupr_stepwise = ye[["aupr"]], K_listed = ye[["K"]], n_tied_at_cap = ye[["n_tied_at_cap"]])
  }))
  if (nrow(gate) && gate$evaluable) pe <- rbind(pe, do.call(rbind, lapply(names(vs), function(m)
    data.frame(net = net, method = m, endpoint = "YM", aupr = d5_module_aupr(vs[[m]], bench), auroc = NA_real_, n_pos = bench$n_pos,
               N = bench$n_pairs, aupr_sens = d5_module_aupr(vs[[m]], bench, expected_ties = TRUE),
               aupr_min10 = if (bench10$n_pos > 0) d5_module_aupr(vs[[m]], bench10) else NA_real_,
               readout = "ap_fast (Y-M)", aupr_stepwise = NA_real_, K_listed = NA_real_, n_tied_at_cap = NA_real_))))
  # cap-boundary check (reported, not hidden): entries tied at the cap boundary, per method (expected 0)
  cap_boundary <- data.frame(net = net, method = methods, cap = cap, n_candidates = vapply(Ms, function(M) sum(!is.na(M)), 0),
                             n_tied_at_cap = vapply(yes, function(r) r[["n_tied_at_cap"]], 0), ties_at_cap = sw$ties_at_cap, row.names = NULL)
  if (any(cap_boundary$n_tied_at_cap > 0)) say("CAP BOUNDARY: ties at the cap in ", paste(cap_boundary$method[cap_boundary$n_tied_at_cap > 0], collapse = ", "), " (ties_at_cap = ", sw$ties_at_cap, ")")
  # official p-values and overall-score ingredients (descriptive; lookup in the precomputed null densities; p = 0 stays 0)
  official <- if (metric == "official" && !is.null(densities)) {
    ye <- pe[pe$endpoint == "YE", ]
    pp <- t(vapply(seq_len(nrow(ye)), function(i) d5_official_p(net, ye$aupr[i], ye$auroc[i], densities), c(P_AUPR = 0, P_AUROC = 0)))
    cbind(ye[, c("net", "method", "aupr", "auroc")], pp, neglog10_P_AUPR = -log10(pp[, "P_AUPR"]), neglog10_P_AUROC = -log10(pp[, "P_AUROC"]), row.names = NULL)
  } else NULL
  # ---- SE_J
  say("truth jackknife, ", nrep, " replicates")
  jE <- d5_jack_edge(Ms, truth, tf_idx, sw, seeds[["YE"]], nrep, cap, universe, metric)
  jM <- if (gate$evaluable) d5_jack_module(vs, bench, G, seeds[["YM"]], nrep) else NULL
  # ---- SE_B (refits recomputed and verified against the committed hashes; the truth and gene universe are fixed)
  bE <- bM <- NULL
  if (length(refit_reps)) {
    say("experiment jackknife, ", length(refit_reps), " replicates")
    rr <- lapply(refit_reps, function(b) {
      v <- refit_scorer(b); if (!is.null(refit_verify)) refit_verify(v, b)
      e <- vapply(names(v), function(m) .d5_pick(d5_edge_metrics(d5_edge_matrix_sym(v[[m]], nd), truth, tf_idx, sw, cap, universe = universe), metric), 0)
      mm <- if (gate$evaluable) vapply(names(v), function(m) d5_module_aupr(v[[m]], bench), 0) else NULL
      list(e = e, m = mm) })
    bE <- do.call(rbind, lapply(rr, `[[`, "e")); bM <- if (gate$evaluable) do.call(rbind, lapply(rr, `[[`, "m")) else NULL
  }
  # ---- supplementary series-level jackknife (net4, descriptive; never enters the verdict). The truth driver does NOT use
  # this in-pipeline path: it runs d5_run_series_supplementary() after verdict.csv is written.
  sE <- sM <- NULL
  if (length(series_reps)) {
    sup <- d5_series_supplement(nd, edges, sw, series_reps, series_scorer, series_verify, universe, cap, contrast = contrast, gate = gate, say = say, metric = metric)
    sE <- sup$series_edge; sM <- sup$series_module
  }
  # ---- effects and verdicts
  a <- contrast[1]; b <- contrast[2]; rows <- list()
  for (ep in c("YE", "YM")) {
    if (ep == "YM" && !gate$evaluable) { rows[[ep]] <- data.frame(net = net, endpoint = ep, verdict = "not evaluable"); next }
    J <- if (ep == "YE") jE else jM; Bm <- if (ep == "YE") bE else bM; p <- pe[pe$endpoint == ep, ]
    A <- p$aupr[p$method == a]; B <- p$aupr[p$method == b]; d <- d5_effect(A, B)
    se_j <- d5_se_paired(J[, a], J[, b]); se_b <- if (!is.null(Bm)) d5_se_paired(Bm[, a], Bm[, b]) else NA_real_
    se <- if (is.na(se_b)) NA_real_ else d5_combine_se(se_j, se_b)
    v <- if (net == 4L) d5_verdict(d, se, delta) else NULL
    # net4 with a non-finite d or SE: none of the design table's conditions holds, so the row reads Inconclusive and says why
    note <- if (net == 4L && !(is.finite(d) && is.finite(se)))
      sprintf("Inconclusive because %s not finite (d=%s, SE_J=%s, SE_B=%s): none of the decision table's conditions holds",
              if (!is.finite(d)) "d is" else "SE is", format(d), format(se_j), format(se_b)) else NA_character_
    rows[[ep]] <- data.frame(net = net, endpoint = ep, readout = if (ep == "YE") D5$ye_readouts[[metric]] else "regulon co-membership AUPR (ap_fast)",
                             contrast = paste(a, "vs", b), aupr_a = A, aupr_b = B, ratio = A / B, difference = A - B, d = d, delta = delta,
                             SE_J = se_j, SE_B = se_b, SE = se, ci_lo = d - 1.96 * se, ci_hi = d + 1.96 * se,
                             ci_J_lo = d - 1.96 * se_j, ci_J_hi = d + 1.96 * se_j,
                             p_sup = if (is.null(v)) NA else v$p_sup, p_harm = if (is.null(v)) NA else v$p_harm,
                             p_noninf_lo = if (is.null(v)) NA else v$p_noninf_lo, p_noninf_hi = if (is.null(v)) NA else v$p_noninf_hi,
                             verdict = if (is.null(v)) "descriptive (no verdict: needs net4 and SE_B)" else as.character(v$verdict), verdict_note = note, stringsAsFactors = FALSE)
  }
  list(point = pe, gate = gate, metric = metric, cap_boundary = cap_boundary, official = official, jack_edge = jE, jack_module = jM, refit_edge = bE, refit_module = bM, rows = rows,
       series_edge = sE, series_module = sM, series_line = d5_series_line(net, a, b, sE, sM, bE, bM, gate))
}

# Supplementary series-level jackknife (net4, descriptive): recomputes each series-level refit, verifies it against the
# committed series hashes, scores it. Returns list(series_edge, series_module, series_line); the line carries SE_B of the
# group-level refits (refit_edge / refit_module) for reference only.
d5_series_supplement <- function(nd, edges, sw, series_reps, series_scorer, series_verify = NULL, universe = NULL, cap = 100000L,
                                 contrast = c("attention", "clr"), gate = d5_ym_gate(edges, nd), refit_edge = NULL, refit_module = NULL,
                                 say = function(...) NULL, metric = "stepwise") {
  truth <- d5_truth_matrix(edges, length(nd$tf_idx), length(nd$genes)); bench <- d5_module_setup(edges, nd, 5L, 500L)
  say("series-level supplementary jackknife, ", length(series_reps), " replicates")
  sr <- lapply(series_reps, function(b) {
    v <- series_scorer(b); if (!is.null(series_verify)) series_verify(v, b)
    e <- vapply(names(v), function(m) .d5_pick(d5_edge_metrics(d5_edge_matrix_sym(v[[m]], nd), truth, nd$tf_idx, sw, cap, universe = universe), metric), 0)
    mm <- if (gate$evaluable) vapply(names(v), function(m) d5_module_aupr(v[[m]], bench), 0) else NULL
    list(e = e, m = mm) })
  sE <- do.call(rbind, lapply(sr, `[[`, "e")); sM <- if (gate$evaluable) do.call(rbind, lapply(sr, `[[`, "m")) else NULL
  list(series_edge = sE, series_module = sM,
       series_line = d5_series_line(nd$net, contrast[1], contrast[2], sE, sM, refit_edge, refit_module, gate))
}
# write the supplementary tables
d5_write_series <- function(sup, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(sup$series_edge)) utils::write.csv(cbind(rep = seq_len(nrow(sup$series_edge)), as.data.frame(sup$series_edge)), file.path(dir, "series_edge.csv"), row.names = FALSE)
  if (!is.null(sup$series_module)) utils::write.csv(cbind(rep = seq_len(nrow(sup$series_module)), as.data.frame(sup$series_module)), file.path(dir, "series_module.csv"), row.names = FALSE)
  if (!is.null(sup$series_line)) utils::write.csv(sup$series_line, file.path(dir, "supplementary_series_jackknife.csv"), row.names = FALSE)
  invisible(dir)
}
# Run the supplementary part AFTER the primary outputs (verdict.csv etc.) are on disk: any failure (a refit hash that differs, a
# missing sidecar, a scoring error) is caught, written to SUPPLEMENTARY_SERIES_FAILED.txt and reported as a labelled
# supplementary failure; it cannot change or block the primary verdict. Returns list(ok, sup | error).
d5_run_series_supplementary <- function(res, dir, ...) {
  failf <- file.path(dir, "SUPPLEMENTARY_SERIES_FAILED.txt"); if (file.exists(failf)) file.remove(failf)
  r <- tryCatch(d5_series_supplement(..., refit_edge = res$refit_edge, refit_module = res$refit_module, gate = res$gate), error = function(e) e)
  if (inherits(r, "error")) {
    msg <- paste0("SUPPLEMENTARY series-level jackknife FAILED (the primary verdict in verdict.csv is unaffected): ", conditionMessage(r))
    writeLines(msg, failf); message(msg); return(invisible(list(ok = FALSE, error = conditionMessage(r))))
  }
  d5_write_series(r, dir); invisible(list(ok = TRUE, sup = r))
}

# ONE supplementary line per endpoint: SD of the paired differences d_r = log2(AUPR_a / AUPR_b) over the series-level replicates.
# Labelled series-level, supplementary and descriptive; it is not used in verdict.csv (SE_B of the group-level refits is shown
# next to it for reference only).
d5_series_line <- function(net, a, b, sE, sM, bE = NULL, bM = NULL, gate = NULL) {
  if (is.null(sE)) return(NULL)
  one <- function(ep, S, B) data.frame(net = net, endpoint = ep, contrast = paste(a, "vs", b), n_rep = nrow(S),
    SD_series = d5_se_paired(S[, a], S[, b]), SE_B_group = if (is.null(B)) NA_real_ else d5_se_paired(B[, a], B[, b]),
    label = "series-level supplementary jackknife (59 #Experiment series, delete-half; descriptive; not used in verdict.csv)", stringsAsFactors = FALSE)
  rbind(one("YE", sE, bE), if (!is.null(sM)) one("YM", sM, bM))
}

# write verdict.csv (the two net4 rows) and the supporting tables
d5_write_results <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  vv <- do.call(rbind, lapply(res$rows, function(r) { miss <- setdiff(names(res$rows$YE), names(r)); for (m in miss) r[[m]] <- NA; r[, names(res$rows$YE)] }))
  utils::write.csv(vv, file.path(dir, "verdict.csv"), row.names = FALSE)
  utils::write.csv(res$point, file.path(dir, "point_estimates.csv"), row.names = FALSE)
  if (!is.null(res$cap_boundary)) utils::write.csv(res$cap_boundary, file.path(dir, "cap_boundary.csv"), row.names = FALSE)
  if (!is.null(res$official)) utils::write.csv(res$official, file.path(dir, "official_scores.csv"), row.names = FALSE)
  utils::write.csv(res$gate, file.path(dir, "ym_gate.csv"), row.names = FALSE)
  utils::write.csv(cbind(rep = seq_len(nrow(res$jack_edge)), as.data.frame(res$jack_edge)), file.path(dir, "jack_edge.csv"), row.names = FALSE)
  if (!is.null(res$jack_module)) utils::write.csv(cbind(rep = seq_len(nrow(res$jack_module)), as.data.frame(res$jack_module)), file.path(dir, "jack_module.csv"), row.names = FALSE)
  if (!is.null(res$refit_edge)) utils::write.csv(cbind(rep = seq_len(nrow(res$refit_edge)), as.data.frame(res$refit_edge)), file.path(dir, "refit_edge.csv"), row.names = FALSE)
  if (!is.null(res$refit_module)) utils::write.csv(cbind(rep = seq_len(nrow(res$refit_module)), as.data.frame(res$refit_module)), file.path(dir, "refit_module.csv"), row.names = FALSE)
  if (!is.null(res$series_line)) d5_write_series(res, dir)
  invisible(vv)
}
