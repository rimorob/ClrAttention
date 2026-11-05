setwd(Sys.getenv("D5_ROOT", unset = "."))

## helpers: exact references by enumeration (tiny cases)
perm_all <- function(v) { if (length(v) <= 1) return(list(v)); out <- list(); for (i in seq_along(v)) for (p in perm_all(v[-i])) out[[length(out) + 1]] <- c(v[i], p); out }
ap_of_order <- function(lab) { tp <- cumsum(lab); sum(lab * tp / seq_along(lab)) / sum(lab) }
auc_of_order <- function(lab) { n1 <- sum(lab); n0 <- length(lab) - n1; r <- rank(-seq_along(lab)); (sum(r[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

test_that("d5_ap_expected equals ap_fast without ties and the exact permutation mean with ties", {
  set.seed(1); s <- rnorm(30); y <- rbinom(30, 1, .3); expect_equal(d5_ap_expected(s, y), ap_fast(s, y))
  s <- c(5, 4, 4, 4, 3, 3, 2, 1); y <- c(0, 1, 0, 1, 1, 0, 0, 1)
  nm <- split(seq_along(s), s); blocks <- nm[order(-as.numeric(names(nm)))]
  orders <- list(integer(0)); for (b in blocks) orders <- unlist(lapply(orders, function(a) lapply(perm_all(b), function(p) c(a, p))), recursive = FALSE)
  ex <- mean(vapply(orders, function(o) ap_of_order(y[o]), 0))
  expect_equal(d5_ap_expected(s, y), ex)
  expect_lte(ap_fast(s, y), d5_ap_expected(s, y) + 1e-12)                      # block-end credit lies below the expected-tie AP here
})

test_that("list + random tail: AP and AUROC equal the exact enumeration over tail orders", {
  # universe of 9: list of 4 (tie-free), tail of 5 unlisted pairs with positives n_rest
  s <- c(9, 8, 7, 6); y <- c(1, 0, 1, 0); n_pos <- 4; N <- 9; n_rest <- n_pos - sum(y)     # 2 positives in the tail
  ex_ap <- mean(vapply(perm_all(c(1, 1, 0, 0, 0)), function(t) ap_of_order(c(y, t)), 0))
  ex_auc <- mean(vapply(perm_all(c(1, 1, 0, 0, 0)), function(t) auc_of_order(c(y, t)), 0))
  for (tm in c("end_of_block", "expected")) {
    r <- d5_list_tail_metrics(s, y, n_pos, N, tm)
    expect_equal(unname(r["aupr"]), ex_ap); expect_equal(unname(r["auroc"]), ex_auc)
  }
  # no tail (list = universe), ties inside the list in both conventions
  s <- c(3, 2, 2, 2, 1); y <- c(1, 0, 1, 1, 0)
  ex_e <- mean(vapply(perm_all(2:4), function(p) ap_of_order(y[c(1, p, 5)]), 0))
  expect_equal(unname(d5_list_tail_metrics(s, y, 3, 5, "expected")["aupr"]), ex_e)
  expect_equal(unname(d5_list_tail_metrics(s, y, 3, 5, "end_of_block")["aupr"]), ap_fast(s, y))
  ex_auc_e <- mean(vapply(perm_all(2:4), function(p) auc_of_order(y[c(1, p, 5)]), 0))
  expect_equal(unname(d5_list_tail_metrics(s, y, 3, 5, "expected")["auroc"]), ex_auc_e)
  expect_equal(unname(d5_list_tail_metrics(s, y, 3, 5, "end_of_block")["auroc"]), auc_of_order(c(1, 0, 1, 1, 0)))
  # zero positives
  expect_true(is.na(unname(d5_list_tail_metrics(c(2, 1), c(0, 0), 0, 5, "expected")["aupr"])))
})

test_that("edge metrics: switches have no defaults; perfect and random scores behave", {
  set.seed(2); nTF <- 6; G <- 30; tf_idx <- 1:6
  tr <- d5_decoy_truth(nTF, G, tf_idx, 11, median_out = 4); T <- d5_truth_matrix(tr, nTF, G)
  M <- matrix(runif(nTF * G), nTF, G); M[cbind(1:nTF, tf_idx)] <- NA
  expect_error(d5_edge_metrics(M, T, tf_idx, d5_edge_switches()), "TRUTH-STAGE switch")
  expect_error(d5_edge_metrics(M, T, tf_idx, d5_edge_switches(ties_at_cap = "bogus", ties_in_list = "expected", tftf = "keep_both", universe = "all")), "ties_at_cap")
  expect_error(d5_edge_metrics(M, T, tf_idx, TESTSW(universe = "gs"), universe = list(tf = 1:3, gene = 1:10)), "cap_order")
  expect_error(d5_edge_metrics(M, T, tf_idx, TESTSW(), keep_tf = 1:3), "jackknife_cap")
  # perfect scores, cap above the number of positives
  Tn <- T; Tn[cbind(1:nTF, tf_idx)] <- FALSE
  P <- ifelse(Tn, 1, 0) + matrix(runif(nTF * G, 0, .01), nTF, G); P[cbind(1:nTF, tf_idx)] <- NA
  r <- d5_edge_metrics(P, T, tf_idx, TESTSW(), cap = sum(Tn)); expect_equal(unname(r["aupr"]), 1); expect_equal(unname(r["auroc"]), 1)
  # capped perfect list: K = cap positives at the top, the rest randomly ordered
  r2 <- d5_edge_metrics(P, T, tf_idx, TESTSW(), cap = 10)
  expect_equal(unname(r2["K"]), 10); expect_lt(unname(r2["aupr"]), 1); expect_gt(unname(r2["aupr"]), 10 / sum(Tn) * 0.9)
  n_rest <- sum(Tn) - 10; Nt <- sum(!is.na(M)) - 10
  expect_equal(unname(r2["aupr"]), (10 + sum((n_rest / Nt) * (10 + 1 + (seq_len(Nt) - 1) * (n_rest - 1) / (Nt - 1)) / (10 + seq_len(Nt)))) / sum(Tn))
  # random scores: AUPR near the base rate (average over draws), AUROC near 0.5
  rr <- replicate(40, { Mr <- matrix(runif(nTF * G), nTF, G); Mr[cbind(1:nTF, tf_idx)] <- NA; d5_edge_metrics(Mr, T, tf_idx, TESTSW(), cap = 1e6)[c("aupr", "auroc")] })
  expect_lt(abs(mean(rr[2, ]) - 0.5), 0.03); expect_lt(abs(mean(rr[1, ]) - sum(Tn) / (nTF * G - nTF)), 0.03)
  # no cap and tie-free: equal to ap_fast on the candidate list
  ok <- !is.na(M); expect_equal(unname(d5_edge_metrics(M, Tn, tf_idx, TESTSW(), cap = 1e6)["aupr"]), ap_fast(M[ok], as.integer(Tn[ok])))
})

test_that("each truth-stage switch does what its documentation says", {
  nTF <- 4; G <- 10; tf_idx <- 1:4
  T <- matrix(FALSE, nTF, G); T[cbind(c(1, 2, 3, 4, 2), c(6, 7, 8, 2, 1))] <- TRUE        # edges 4->2 and 2->1 are TF-TF
  M <- matrix(0, nTF, G); M[cbind(1:4, tf_idx)] <- NA
  vals <- c(10, 9, 8, 8, 8, 8, 7, 6, 5, 4, 3, 2, 1); M[!is.na(M)][1:13] <- NA                # placeholder (overwritten below)
  M <- matrix(NA_real_, nTF, G); cells <- which(row(M) != matrix(tf_idx, nTF, G, byrow = FALSE)[, 1][row(M)] | col(M) != tf_idx[row(M)])
  M[!(col(M) == tf_idx[row(M)])] <- rev(seq_len(sum(!(col(M) == tf_idx[row(M)])))) / 10    # distinct scores
  # ties straddling the cap: set cells 20..24 (in descending order) equal
  o <- order(M, decreasing = TRUE, na.last = NA); M[o[5:9]] <- M[o[5]]
  cap <- 7
  r_keep <- d5_edge_metrics(M, T, tf_idx, TESTSW(ties_at_cap = "keep_tied"), cap); r_drop <- d5_edge_metrics(M, T, tf_idx, TESTSW(ties_at_cap = "drop_tied"), cap)
  r_ord <- d5_edge_metrics(M, T, tf_idx, TESTSW(ties_at_cap = "list_order"), cap)
  expect_equal(unname(c(r_keep["K"], r_drop["K"], r_ord["K"])), c(9, 4, 7))
  # TF-TF handling
  Nall <- sum(!is.na(M))
  expect_equal(unname(d5_edge_metrics(M, T, tf_idx, TESTSW(tftf = "keep_both"), 1e6)["N"]), Nall)
  expect_equal(unname(d5_edge_metrics(M, T, tf_idx, TESTSW(tftf = "drop"), 1e6)["N"]), Nall - nTF * (nTF - 1))
  rc <- d5_edge_metrics(M, T, tf_idx, TESTSW(tftf = "collapse_max"), 1e6)
  expect_equal(unname(rc["N"]), Nall - nTF * (nTF - 1) / 2)
  expect_equal(unname(d5_edge_metrics(M, T, tf_idx, TESTSW(tftf = "keep_both"), 1e6)["n_pos"]), 5)
  expect_equal(unname(d5_edge_metrics(M, T, tf_idx, TESTSW(tftf = "drop"), 1e6)["n_pos"]), 3)
  expect_equal(unname(rc["n_pos"]), 5)         # 4->2 and 2->1 are different TF pairs; each is positive in its kept orientation
  T2 <- T; T2[2, 4] <- TRUE                    # both orientations of the pair (2,4) are edges: one positive after collapsing
  expect_equal(unname(d5_edge_metrics(M, T2, tf_idx, TESTSW(tftf = "collapse_max"), 1e6)["n_pos"]), 5)
  expect_equal(unname(d5_edge_metrics(M, T2, tf_idx, TESTSW(tftf = "keep_both"), 1e6)["n_pos"]), 6)
  # collapse keeps the higher-scoring orientation: make 3->1 score high, 1->3 low; the kept cell is 3->1
  M3 <- M; M3[3, 1] <- 100; M3[1, 3] <- -1
  T3 <- matrix(FALSE, nTF, G); T3[3, 1] <- TRUE
  expect_equal(unname(d5_edge_metrics(M3, T3, tf_idx, TESTSW(tftf = "collapse_max"), 1)["aupr"]), 1)    # cap 1: the top cell is the positive one
  # pair universe and cap order
  u <- list(tf = 1:2, gene = 5:10)
  a <- d5_edge_metrics(M, T, tf_idx, TESTSW(universe = "gs", cap_order = "cap_then_universe"), cap = 20, universe = u)
  b <- d5_edge_metrics(M, T, tf_idx, TESTSW(universe = "gs", cap_order = "universe_then_cap"), cap = 20, universe = u)
  nu <- sum(!is.na(M[1:2, 5:10])); expect_equal(unname(c(a["N"], b["N"])), c(nu, nu))
  expect_equal(unname(b["K"]), min(20, nu)); expect_lt(unname(a["K"]), unname(b["K"]) + 1)
  expect_equal(unname(a["n_pos"]), sum(T[1:2, 5:10]))
  # truth jackknife restriction
  kt <- c(1, 3); kg <- c(2, 6, 7, 8)
  nk <- sum(!is.na(M[kt, kg]))
  ra <- d5_edge_metrics(M, T, tf_idx, TESTSW(jackknife_cap = "reapply_cap"), 12, kt, kg)
  rs <- d5_edge_metrics(M, T, tf_idx, TESTSW(jackknife_cap = "scaled_cap"), 12, kt, kg)
  rg <- d5_edge_metrics(M, T, tf_idx, TESTSW(jackknife_cap = "restrict_after_cap"), 12, kt, kg)
  expect_equal(unname(c(ra["N"], rs["N"], rg["N"])), c(nk, nk, nk))
  expect_equal(unname(ra["cap_used"]), 12); expect_equal(unname(rs["cap_used"]), max(1, round(12 * 8 / 40)))
  expect_lte(unname(rg["K"]), 12); expect_gte(unname(ra["K"]), unname(rg["K"]))
})

test_that("module scorer: regulons, pairs and evaluability gate on decoy truths", {
  nd <- d5_synth_net(120, 10, 30, 5); G <- 120
  e <- d5_decoy_truth(10, G, nd$tf_idx, 3, median_out = 9)
  expect_false(any(e$gene == nd$tf_idx[e$tf]))                                          # no self edges
  expect_identical(e, d5_decoy_truth(10, G, nd$tf_idx, 3, median_out = 9)); expect_false(identical(e, d5_decoy_truth(10, G, nd$tf_idx, 4, median_out = 9)))
  b <- d5_module_setup(e, nd)
  sz <- table(e$tf); keep <- names(sz)[sz >= 5 & sz <= 500]
  expect_equal(b$n_regulons, length(keep))
  # brute force positive pairs
  genes_in <- sort(unique(e$gene[e$tf %in% as.integer(keep)])); pos <- 0
  mem <- lapply(as.integer(keep), function(a) e$gene[e$tf == a]); pr <- combn(genes_in, 2)
  pos <- sum(apply(pr, 2, function(p) any(vapply(mem, function(m) all(p %in% m), NA)))); expect_equal(b$n_pos, pos)
  expect_equal(b$n_pairs, choose(length(genes_in), 2))
  # a perfect scorer reaches AUPR 1 under both tie conventions; scores via the packed upper triangle
  S <- matrix(0, G, G); S[cbind(b$pairs$i, b$pairs$j)] <- b$pairs$label + runif(length(b$pairs$i), 0, 1e-3); S <- pmax(S, t(S))
  v <- d5_pack(S); expect_equal(d5_module_aupr(v, b), 1); expect_equal(d5_module_aupr(v, b, expected_ties = TRUE), 1)
  # limits: >= 10 targets gives a subset
  b10 <- d5_module_setup(e, nd, 10L, 500L); expect_lte(b10$n_regulons, b$n_regulons)
  # gate
  g <- d5_ym_gate(e, nd); expect_equal(g$n_regulons_5_500, length(keep)); expect_equal(g$n_regulons_ge5, length(keep))   # no regulon exceeds 500 here
  expect_equal(g$evaluable, length(keep) >= 20 && g$n_pos_pairs >= 1000)
  e2 <- d5_decoy_truth(10, G, nd$tf_idx, 3, median_out = 2, sd_log = 0.1); expect_false(d5_ym_gate(e2, nd)$evaluable)
  # a case that is evaluable
  nd2 <- d5_synth_net(300, 40, 20, 6); e3 <- d5_decoy_truth(40, 300, nd2$tf_idx, 9, median_out = 12, sd_log = 0.4)
  g3 <- d5_ym_gate(e3, nd2); expect_gte(g3$n_regulons_ge5, 20); expect_true(g3$evaluable)
  sc <- d5_structure_counts(e3, nd2); expect_equal(sc$n_edges, nrow(e3))
})

test_that("Y-M gate: both counts on the scored 5-500 regulon set; the unrestricted count is an extra, non-deciding field (change log 2026-11-05)", {
  nd <- d5_synth_net(1300, 60, 20, 8); G <- 1300; set.seed(3)
  targets <- function(a, n) sample(setdiff(seq_len(G), nd$tf_idx[a]), n)
  mk <- function(sizes) do.call(rbind, lapply(seq_along(sizes), function(a) data.frame(tf = a, gene = targets(a, sizes[a]))))
  # (a) the unrestricted count passes (25 regulons with >= 5 targets) but only 15 regulons are in 5-500 -> NOT evaluable
  e <- mk(c(rep(30, 15), rep(600, 10))); g <- d5_ym_gate(e, nd)
  expect_equal(g$n_regulons_ge5, 25L); expect_equal(g$n_regulons_5_500, 15L)
  expect_gte(g$n_pos_pairs, 1000); expect_false(g$evaluable)                      # the old (unrestricted) rule would have passed
  expect_true(g$n_regulons_ge5 >= 20 && g$n_pos_pairs >= 1000)
  # (b) 20 regulons in 5-500 (so the regulon count passes) but < 1,000 positive pairs there, while the two huge regulons make
  # the unrestricted pair count large: NOT evaluable; the unrestricted pair count is reported but does not decide
  e <- mk(c(rep(5, 20), rep(600, 2))); g <- d5_ym_gate(e, nd)
  expect_equal(g$n_regulons_5_500, 20L); expect_lt(g$n_pos_pairs, 1000); expect_gt(g$n_pos_pairs_ge5_noupper, 1000); expect_false(g$evaluable)
  # (c) both counts pass on the 5-500 set -> evaluable, whatever sits above 500
  e <- mk(c(rep(30, 20), 600)); g <- d5_ym_gate(e, nd)
  expect_equal(g$n_regulons_5_500, 20L); expect_gte(g$n_pos_pairs, 1000); expect_true(g$evaluable); expect_equal(g$n_regulons_ge5, 21L)
  # (d) exactly at the boundary: 19 regulons of 5-500 plus one of 600 is not evaluable (unrestricted count 20); 20 are
  e <- mk(c(rep(30, 19), 600)); g <- d5_ym_gate(e, nd); expect_equal(g$n_regulons_ge5, 20L); expect_false(g$evaluable)
  # the gate columns are all reported and the structure counts carry them
  expect_true(all(c("n_regulons_5_500", "n_pos_pairs", "n_regulons_ge5", "n_pos_pairs_ge5_noupper", "evaluable") %in% names(g)))
  expect_true(all(c("n_regulons_ge5", "evaluable") %in% names(d5_structure_counts(e, nd))))
})

test_that("truth jackknife: design seeds, halves, shared across methods, reproducible, caller RNG untouched", {
  h <- d5_jack_halves_edge(20261203L, 333, 5950, 200L)
  expect_length(h, 200); expect_length(h[[1]]$tf, 166); expect_length(h[[1]]$gene, 2975)
  expect_identical(h, d5_jack_halves_edge(20261203L, 333, 5950, 200L)); expect_false(identical(h[[1]], d5_jack_halves_edge(20261204L, 333, 5950, 200L)[[1]]))
  set.seed(20261203L); t1 <- sort(sample.int(333, 166)); g1 <- sort(sample.int(5950, 2975)); expect_identical(h[[1]]$tf, t1); expect_identical(h[[1]]$gene, g1)   # draw order: TFs then genes
  set.seed(9); a <- runif(1); set.seed(9); invisible(d5_jack_halves_module(20261204L, 5950, 200L)); expect_equal(runif(1), a)
  expect_identical(d5_jack_halves_module(20261204L, 100L, 3L), d5_jack_halves_module(20261204L, 100L, 3L))
  expect_identical(D5$seeds[["4"]], c(YE = 20261203L, YM = 20261204L)); expect_identical(D5$seeds[["3"]], c(YE = 20261213L, YM = 20261214L))
  expect_identical(D5$seeds[["1"]], c(YE = 20261223L, YM = 20261224L))
  # replicates reproducible; same halves for every method
  nTF <- 6; G <- 40; tf_idx <- 1:6; T <- d5_truth_matrix(d5_decoy_truth(nTF, G, tf_idx, 2, 5), nTF, G)
  mk <- function(s) { set.seed(s); M <- matrix(runif(nTF * G), nTF, G); M[cbind(1:nTF, tf_idx)] <- NA; M }
  Ms <- list(a = mk(1), b = mk(2), a2 = mk(1)); sw <- TESTSW(jackknife_cap = "reapply_cap")
  j1 <- d5_jack_edge(Ms, T, tf_idx, sw, 77L, 8L, cap = 30); j2 <- d5_jack_edge(Ms, T, tf_idx, sw, 77L, 8L, cap = 30)
  expect_identical(j1, j2); expect_equal(dim(j1), c(8, 3)); expect_identical(j1[, "a"], j1[, "a2"]); expect_false(identical(j1[, "a"], j1[, "b"]))
  expect_false(identical(j1, d5_jack_edge(Ms, T, tf_idx, sw, 78L, 8L, cap = 30)))
  # module jackknife: same halves, restricted to kept genes' pairs
  nd <- d5_synth_net(80, 8, 25, 7); e <- d5_decoy_truth(8, 80, nd$tf_idx, 5, 8); b <- d5_module_setup(e, nd)
  vs <- list(x = d5_pack(abs(cor(t(nd$X)))), y = d5_pack(abs(cor(t(nd$X)))))
  jm <- d5_jack_module(vs, b, 80, 20261204L, 6L); expect_identical(jm[, "x"], jm[, "y"]); expect_identical(jm, d5_jack_module(vs, b, 80, 20261204L, 6L))
  h2 <- d5_jack_halves_module(20261204L, 80, 6L); inh <- logical(80); inh[h2[[1]]] <- TRUE; ix <- which(inh[b$pairs$i] & inh[b$pairs$j])
  expect_equal(unname(jm[1, "x"]), ap_fast(d5_pair_scores(vs$x, b$pairs)[ix], b$pairs$label[ix]))
  expect_equal(d5_se_paired(c(1, 2, 3), c(1, 1, 1)), sd(log2(c(1, 2, 3)))); expect_equal(d5_combine_se(3, 4), 5)
})
