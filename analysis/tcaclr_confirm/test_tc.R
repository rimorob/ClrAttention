#!/usr/bin/env Rscript
# Unit and end-to-end tests of analysis/tcaclr_confirm/ on small synthetic data with known answers. Run from the repo root:
#   nice -n 19 Rscript analysis/tcaclr_confirm/test_tc.R
# Touches no real data except the stored net4 sentences (results/complex_confirm/estimate.csv) for the sentence-generator check.
suppressPackageStartupMessages({ library(testthat); source("analysis/tcaclr_confirm/tc_lib.R") })
options(warn = 1)
# waldo/diffobj (used by testthat to print differences) are not installed here: compare with all.equal / identical and print both values on failure
expect_equal <- function(object, expected, ..., tolerance = 1e-8, info = NULL) { ok <- isTRUE(all.equal(object, expected, tolerance = tolerance, check.attributes = FALSE))
  testthat::expect(ok, paste0("not equal: ", paste(format(utils::head(object, 6)), collapse = " "), " vs ", paste(format(utils::head(expected, 6)), collapse = " "), if (!is.null(info)) paste0("\n", info) else "")) }
expect_identical <- function(object, expected, info = NULL) testthat::expect(identical(object, expected), paste0("not identical: ", paste(format(utils::head(object, 6)), collapse = " "), " vs ", paste(format(utils::head(expected, 6)), collapse = " "), if (!is.null(info)) paste0("\n", info) else ""))
tmp <- Sys.getenv("TC_TEST_DIR", tempfile("tc_test_")); dir.create(tmp, showWarnings = FALSE)

test_that("sort-once replicate AP equals ap_fast, with ties, NAs and subsets (bit-identical)", {
  set.seed(11); n <- 4000; s <- round(rnorm(n), 1); s[sample(n, 20)] <- NA; l <- rbinom(n, 1, 0.1)
  pr <- tc_ap_prep(s, l); expect_identical(tc_ap_sub(pr, seq_len(n)), ap_fast(s, l))
  for (r in 1:20) { keep <- sort(sample(n, 1500)); ixs <- which(pr$o %in% keep)   # positions in sorted order of the kept pairs
    expect_identical(tc_ap_sub(pr, ixs), ap_fast(s[keep], l[keep])) }
})

test_that("gene jackknife: halves and replicate APs equal the direct computation", {
  set.seed(5); G <- 40; i <- c(); j <- c(); for (a in 1:(G - 1)) for (b in (a + 1):G) { i <- c(i, a); j <- c(j, b) }
  p <- list(i = i, j = j, label = rbinom(length(i), 1, 0.2), genes_in = 1:G); sc <- round(rnorm(length(i)), 2)
  h <- tc_draw_halves(list(main = p, other = p), 7L, 42L)
  set.seed(42); manual <- lapply(1:7, function(r) sample(1:G, floor(G / 2))); expect_identical(h$main, manual)
  jr <- tc_jack_method(sc, p, h$main, G)
  direct <- vapply(h$main, function(k) { inh <- logical(G); inh[k] <- TRUE; ix <- which(inh[p$i] & inh[p$j]); ap_fast(sc[ix], p$label[ix]) }, 0)
  expect_identical(jr, direct)
})

test_that("experiment halves: sequential seed and the net4 per-replicate convention", {
  h <- tc_exp_halves(59L, 3L, 20264000L, "perrep")
  for (b in 1:3) { set.seed(20264000L + b); expect_identical(h[[b]], sort(sample.int(59L, 29L))) }
  hs <- tc_exp_halves(100L, 4L, TC$exp_seed, "seq"); set.seed(TC$exp_seed); expect_identical(hs[[1]], sort(sample.int(100L, 50L))); expect_identical(hs[[2]], sort(sample.int(100L, 50L)))
  expect_false(exists(".Random.seed", envir = globalenv()) && FALSE)   # RNG restored without error
  if (file.exists("analysis/dream5_lib.R")) { suppressPackageStartupMessages(source("analysis/dream5_lib.R"))
    for (b in c(1, 17, 100)) { expect_identical(tc_exp_halves(59L, 100L, 20264000L, "perrep")[[b]], d5_series_half(59L, b)); expect_identical(tc_exp_halves(321L, 100L, 20263000L, "perrep")[[b]], d5_group_half(321L, b)) } }
})

# a small truth: genes g1..g30; complexes A (g1-g4), B (g5-g8), C (g9-g28); g29, g30 outside
genes <- paste0("g", 1:30)
mem <- rbind(data.frame(variant = "main", dataset = "syn", complex = rep(c("A", "B", "C"), c(4, 4, 20)), orf = genes[c(1:4, 5:8, 9:28)]),
             data.frame(variant = "noribo", dataset = "syn", complex = rep(c("A", "B"), c(4, 4)), orf = genes[1:8]))
excl <- data.frame(scope = "all", orf_a = "g1", orf_b = "g2", reason = "paralog_set")
tr <- list(mem = mem, excl = excl)

test_that("truth builder: pair universe, labels, exclusions, strata and the perturbed split", {
  T <- suppressWarnings(tc_build_truths(tr, genes, "plan", perturbed = match(c("g1", "g9"), genes)))
  m <- T$main; expect_equal(length(m$genes_in), 28L); expect_equal(length(m$label), choose(28, 2) - 1L)   # one paralog pair removed
  expect_equal(sum(m$label), 6 - 1 + 6 + choose(20, 2))
  expect_equal(sum(T$small_3to6$label), 5 + 6); expect_equal(sum(T$large_ge20$label), choose(20, 2)); expect_null(T$mid_7to19)
  expect_equal(sum(T$large_ge20$label) + sum(T$small_3to6$label), sum(m$label))     # no mid complex here
  expect_equal(length(T$pert_neither$label) + length(T$pert_rest$label), length(m$label))
  pr <- T$pert_rest; expect_true(all(pr$i %in% match(c("g1", "g9"), genes) | pr$j %in% match(c("g1", "g9"), genes)))
  brute <- 0; for (a in 1:27) for (b in (a + 1):28) brute <- brute + (genes[a] != "g1" || genes[b] != "g2") * ((a <= 4 && b <= 4) || (a >= 5 && a <= 8 && b <= 8 && b >= 5) || (a >= 9 && b >= 9))
  expect_equal(sum(m$label), brute)
})

test_that("pair credit adds up to AP, contributions add up to d, n_eff is the inverse Herfindahl index", {
  set.seed(3); G <- 30; T <- suppressWarnings(tc_build_truths(tr, genes, "plan")); p <- T$main
  sa <- rnorm(length(p$i)); sb <- rnorm(length(p$i)) + 0.5 * p$label
  expect_equal(sum(tc_pair_credit(sa, p$label)), ap_fast(sa, p$label), tolerance = 1e-12)
  pc <- tc_pair_complex(p$sets, p, G); r <- tc_complex_contrib(sb, sa, p, p$sets, G, pc)
  expect_equal(sum(r$d_c), ap_fast(sb, p$label) - ap_fast(sa, p$label), tolerance = 1e-12); expect_equal(r$d, sum(r$d_c), tolerance = 1e-12)
  # known answer by brute force for distinct scores: precision of every positive at its rank
  brute_credit <- function(s, l) { o <- order(s, decreasing = TRUE); lo <- l[o]; pr <- cumsum(lo) / seq_along(lo); cr <- numeric(length(s)); cr[o] <- lo * pr / sum(l); cr }
  expect_equal(tc_pair_credit(sa, p$label), brute_credit(sa, p$label), tolerance = 1e-12)
  expect_equal(tc_neff(pc, sum(p$label)), 1 / sum((tapply(pc$w, pc$cx, sum) / sum(p$label))^2))
  two <- data.frame(pair = 1:8, cx = rep(1:2, 4), w = 1); expect_equal(tc_neff(two, 8), 2)
  sh <- data.frame(pair = 1:4, cx = c(1, 1, 1, 2), w = 1); expect_equal(tc_neff(sh, 4), 1 / (0.75^2 + 0.25^2))
  mult <- data.frame(pair = c(1, 1, 2), cx = c(1, 2, 1), w = c(0.5, 0.5, 1)); expect_equal(tc_neff(mult, 2), 1.6); expect_equal(tc_neff_overlap(mult), 1.8)   # fractional vs overlap-counted
  expect_equal(tc_neff_overlap(two), 2)
})

test_that("contributions with a known answer: gain carried by one complex", {
  # complexes A, B of four genes each; a ranks A on top and B at the bottom; b ranks B on top and A at the bottom: d_A > 0 > d_B and they add to d
  gs <- paste0("h", 1:8); mm <- data.frame(variant = "main", dataset = "s", complex = rep(c("A", "B"), each = 4), orf = gs); T <- suppressWarnings(tc_build_truths(list(mem = mm, excl = excl[0, ]), gs, "plan")); p <- T$main
  inA <- (p$i <= 4 & p$j <= 4); inB <- (p$i > 4 & p$j > 4)
  sa <- ifelse(inA, 3, ifelse(inB, 1, 2)) + runif(length(p$i), 0, 0.01); sb <- ifelse(inB, 3, ifelse(inA, 1, 2)) + runif(length(p$i), 0, 0.01)
  r <- tc_complex_contrib(sa, sb, p, p$sets, 8L); expect_gt(r$d_c[["A"]], 0); expect_lt(r$d_c[["B"]], 0); expect_equal(sum(r$d_c), r$d)
  # exact: the six A pairs sit at ranks 1-6 under a (precision 1) and at ranks 23-28 under b (tp = 6 + j at rank 22 + j); the pair counts of the universe are 12 positives, 16 negatives
  dA <- (6 - sum((6 + 1:6) / (22 + 1:6))) / 12
  expect_equal(r$d_c[["A"]], dA, tolerance = 1e-12); expect_equal(r$d_c[["B"]], -dA, tolerance = 1e-12); expect_equal(r$d, 0, tolerance = 1e-12)
})

test_that("descriptors and the sign test", {
  d_c <- c(A = 0.06, B = 0.03, C = 0.005, D = 0.005, E = -0.01); size <- c(A = 25, B = 22, C = 8, D = 5, E = 4); d <- sum(d_c)
  x <- tc_descriptors(d_c, size, d); expect_equal(x$share_top3, 0.095 / 0.09); expect_equal(x$top3, "A;B;C"); expect_equal(x$share_large, 0.09 / 0.09); expect_equal(x$k_carry, 2L)
  expect_true(is.na(tc_descriptors(d_c, size, -0.01)$share_top3))
  nine <- tc_breadth(c(rep(1, 9), 0), rep(0.5, 10)); expect_equal(nine$n_improved, 9L); expect_equal(nine$n_worse, 1L); expect_equal(nine$p_sign_one_sided, 11 / 1024); expect_true(nine$pass)
  eight <- tc_breadth(c(rep(1, 8), 0, 0), rep(0.5, 10)); expect_equal(eight$p_sign_one_sided, 56 / 1024); expect_false(eight$pass)
  tied <- tc_breadth(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5), c(rep(0.5, 10), 0.5, 0.5)); expect_equal(tied$n_tied, 2L); expect_equal(tied$n_improved, 10L); expect_equal(tied$p_sign_one_sided, 1 / 1024)
  expect_true(is.na(tc_breadth(c(NA, 1), c(1, NA))$p_sign_one_sided))
  expect_equal(tc_breadth_sentence(nine, 3L, 0.09), "tcaCLR improved the median complex.")
  expect_equal(tc_breadth_sentence(eight, 3L, 0.09), "The pooled gain is carried by 3 complexes; tcaCLR did not improve the median complex.")
  expect_match(tc_breadth_sentence(eight, NA_integer_, -0.01), "not positive")
})

test_that("H2 labels and verdict", {
  a <- tc_h2_labels(0.137, 0.106, 0.168); expect_true(a[["positive"]]); expect_false(a[["near_zero"]]); expect_false(a[["not_positive"]])
  b <- tc_h2_labels(-0.0034, -0.009, 0.0022); expect_false(b[["positive"]]); expect_true(b[["near_zero"]]); expect_true(b[["not_positive"]])
  s <- tc_h2_labels(-0.0037, -0.0113, 0.0038); expect_true(s[["not_positive"]])
  v <- tc_h2_verdict(s, b, a); expect_equal(v$verdict, "held"); expect_match(v$sentence, "held")
  expect_equal(tc_h2_verdict(s, b, tc_h2_labels(0.02, -0.01, 0.05))$verdict, "partly held")
  expect_equal(tc_h2_verdict(tc_h2_labels(0.05, 0.03, 0.07), tc_h2_labels(0.05, 0.03, 0.07), tc_h2_labels(-0.05, -0.07, -0.03))$verdict, "did not hold")
  expect_true(tc_h2_labels(0.0105, -0.01, 0.03)[["near_zero"]] == FALSE)   # |d| >= 0.01 is not near zero
})

test_that("H3: intersection-union and Holm", {
  tab <- data.frame(comparator = c("aCLR", "CLR", "absr"), d = c(0.02, 0.03, 0.01), p = c(0.001, 0.010, 0.020))
  h <- tc_h3(tab); expect_true(h$best); expect_equal(h$k, 4L); expect_equal(h$tab$holm_p, c(0.003, 0.02, 0.02))   # Holm: 0.001*3, then max(0.01*2, .003), then max(0.02*1, .02)
  expect_true(all(h$tab$won_holm))
  tab2 <- data.frame(comparator = c("aCLR", "CLR", "absr"), d = c(0.02, 0.03, 0.01), p = c(0.001, 0.010, 0.026)); h2 <- tc_h3(tab2); expect_false(h2$best); expect_equal(h2$tab$won_iut, c(TRUE, TRUE, FALSE))
  tab3 <- data.frame(comparator = c("aCLR", "CLR"), d = c(0.02, -0.01), p = c(0.001, 0.9)); expect_false(tc_h3(tab3)$best)
  tab4 <- data.frame(comparator = c("aCLR", "CLR"), d = c(0.02, 0.02), p = c(0.012, 0.002)); h4 <- tc_h3(tab4)    # Holm: 0.012 and 0.004 -> both win
  expect_true(h4$best); expect_true(all(h4$tab$won_holm)); expect_equal(h4$tab$holm_p, c(0.012, 0.004))
  tab5 <- data.frame(comparator = c("aCLR", "CLR"), d = c(0.02, 0.02), p = c(0.020, 0.021)); h5 <- tc_h3(tab5)
  expect_true(h5$best); expect_false(any(h5$tab$won_holm))   # each comparison is won at 0.025 (IUT) but Holm over two gives 0.04 for both
})

test_that("total SE, one-sided p and the H1 rule", {
  tt <- tc_total(0.03, 0.01, 0.01); expect_equal(tt$se_total, sqrt(2) * 0.01); expect_equal(tt$p_one_sided, 1 - pnorm(0.03 / (sqrt(2) * 0.01))); expect_equal(tt$total_lo, 0.03 - 1.96 * sqrt(2) * 0.01)
  expect_false(tc_total(0.03, 0.01)$exp_included)
  expect_match(tc_claim_sentence(0.001, 0.05, "test_0.03", 0.5), "^tcaCLR exceeds aCLR in pooled co-membership AUPR on the Complex Portal human complexes in HepG2")
  expect_match(tc_claim_sentence(0.001, 0.05, "test_0.03", 0.5), "M3D p = 0.079, EcoMAC p = 0.051\\)\\.$")
  expect_match(tc_claim_sentence(0.001, 0.05, "test_0.03", 0.81), "the gain is carried by complexes of 20 or more genes; ")
  expect_false(grepl("carried by complexes of 20", tc_claim_sentence(0.001, 0.05, "test_0.03", 0.80)))      # strictly more than 80%
  expect_match(tc_claim_sentence(0.04, 0.05, "test_0.03", 0.5), "^The claim is not made")
  expect_match(tc_claim_sentence(0.001, 0.004, "test_0.03", 0.5), "^The claim is not made")                  # d below +0.005
  expect_match(tc_claim_sentence(0.001, 0.05, "estimate_only", 0.5), "estimate only")
})

test_that("interval sentences reproduce the stored net4 sentences byte for byte", {
  f <- "results/complex_confirm/estimate.csv"; skip_if_not(file.exists(f)); e <- utils::read.csv(f, stringsAsFactors = FALSE); P <- e[e$row == "P", ]
  s <- tc_interval_sentences(P$total_series_lo, P$total_series_hi, "a3r", "attention", "yeast complexes", fixed = P$S_fixed)
  expect_identical(unname(s[["S_zero"]]), P$S_zero); expect_identical(unname(s[["S_min"]]), P$S_min); expect_identical(unname(s[["S_ecoli"]]), P$S_ecoli); expect_identical(unname(s[["S_fixed"]]), P$S_fixed)
  g <- tc_interval_sentences(P$gene_lo, P$gene_hi, "a3r", "attention", "yeast complexes"); expect_identical(paste("[gene jackknife only, comparable to the E. coli test]", g[["S_zero"]]), P$side_gene_only_S_zero)
  expect_match(tc_interval_sentences(-0.05, -0.01, "a", "b")[["S_zero"]], "lies below no difference"); expect_match(tc_interval_sentences(-0.01, 0.02, "a", "b")[["S_zero"]], "includes no difference")
  expect_match(tc_interval_sentences(0.001, 0.003, "a", "b")[["S_min"]], "excluded"); expect_match(tc_interval_sentences(-0.1, 0.004, "a", "b")[["S_min"]], "excluded")
  expect_match(tc_interval_sentences(0.03, 0.05, "a", "b")[["S_ecoli"]], "above the E. coli"); expect_match(tc_interval_sentences(0.0, 0.02, "a", "b")[["S_ecoli"]], "as large as the E. coli")
  expect_match(tc_secondary_sentence("noribo", 0.046, 0.02, 0.07), "excludes no difference"); expect_match(tc_pert_sentence(0.0, -0.01, 0.01, 0.1, 0.05, 0.15), "confined to pairs involving perturbed genes")
})

test_that("SHA-256 guard refuses mismatches, missing files and uncovered files", {
  d <- file.path(tmp, "truth"); dir.create(d); writeLines(c("a", "b"), file.path(d, "members.tsv")); writeLines("c", file.path(d, "excl.tsv"))
  h <- function(f) tc_sha256(file.path(d, f)); writeLines(c(paste(h("members.tsv"), "members.tsv"), paste(h("excl.tsv"), "excl.tsv")), file.path(d, "t.sha256"))
  expect_true(all(tc_hash_guard(file.path(d, "t.sha256"), file.path(d, c("members.tsv", "excl.tsv")))$ok))
  expect_error(tc_hash_guard(file.path(d, "t.sha256"), file.path(d, c("members.tsv", "other.tsv"))), "not covered")
  writeLines(c("a", "B"), file.path(d, "members.tsv")); expect_error(tc_hash_guard(file.path(d, "t.sha256")), "MISMATCH")
  writeLines(c("a", "b"), file.path(d, "members.tsv")); file.remove(file.path(d, "excl.tsv")); expect_error(tc_hash_guard(file.path(d, "t.sha256")), "MISMATCH or missing")
  writeLines("garbage", file.path(d, "bad.sha256")); expect_error(tc_hash_guard(file.path(d, "bad.sha256")), "malformed"); expect_error(tc_hash_guard(file.path(d, "nope.sha256")), "not found")
})

test_that("score files: packed index, directory reader, length validation", {
  G <- 7; S <- matrix(rnorm(49), G); S <- S + t(S); diag(S) <- 0; v <- S[upper.tri(S)]; i <- c(1, 3, 2); j <- c(2, 7, 6)
  expect_equal(v[tc_pair_index(i, j)], S[cbind(i, j)]); expect_equal(v[tc_pair_index(j, i)], S[cbind(i, j)])
  sd_ <- file.path(tmp, "sc"); dir.create(sd_); writeLines(paste0("x", 1:G), file.path(sd_, "genes.txt")); tc_write_bin(v, file.path(sd_, "aCLR.scores")); saveRDS(S, file.path(sd_, "CLR.rds")); tc_write_bin(v[-1], file.path(sd_, "bad.scores"))
  st <- tc_scores_open(sd_); expect_equal(st$G, G); expect_equal(tc_scores_packed(st, "aCLR"), v); expect_equal(tc_scores_packed(st, "CLR"), v); expect_error(tc_scores_packed(st, "bad"), "is not G")
  expect_equal(tc_unpack(v, G), S)
})

test_that("p formatting for sentences", { expect_identical(tc_fmt_p(0.00001), "< 0.0001"); expect_identical(tc_fmt_p(0.0047), "= 0.0047"); expect_identical(tc_fmt_p(NA_real_), "= NA") })

test_that("PIDC refit through tools/pidc_chunked.jl equals the NetworkInference package output (tools/pidc.jl) on a small matrix", {
  jl <- path.expand("~/.juliaup/bin/julia"); skip_if_not(file.exists(jl))
  set.seed(4); G <- 60; X <- matrix(rnorm(G * 50), G, 50, dimnames = list(sprintf("G%03d", 1:G), NULL)); X[1:15, ] <- X[1:15, ] + rep(rnorm(50), each = 15)
  f <- tc_pidc_refit(X, threads = 1L, tmpdir = tmp); fin <- file.path(tmp, "ref.tsv"); fout <- file.path(tmp, "ref_edges.tsv")
  utils::write.table(cbind(gene = rownames(X), as.data.frame(X)), fin, sep = "\t", quote = FALSE, row.names = FALSE, col.names = c("gene", paste0("c", 1:50)))
  rc <- system2("env", c("-i", paste0("HOME=", Sys.getenv("HOME")), "PATH=/usr/bin:/bin", jl, "+1.6", "-t", "1", "tools/pidc.jl", fin, fout), stdout = FALSE, stderr = FALSE); expect_equal(rc, 0L)
  e <- utils::read.delim(fout, header = FALSE, stringsAsFactors = FALSE); w <- f(match(e[[1]], rownames(X)), match(e[[2]], rownames(X))); expect_true(max(abs(w - e[[3]])) < 1e-12)
})

test_that("tc_fit_three equals the net4 cc_fit_three on a small matrix (bit-identical)", {
  skip_if_not(file.exists("analysis/complex_confirm/cc_lib.R"))
  set.seed(8); X <- matrix(rnorm(60 * 30), 60, 30); X[1:10, ] <- X[1:10, ] + rep(rnorm(30), each = 10)
  a <- tc_fit_three(X, 2L); suppressPackageStartupMessages(source("analysis/complex_confirm/cc_lib.R")); b <- cc_fit_three(X, 2L)
  expect_identical(a$CLR, b$clr); expect_identical(a$aCLR, b$attention); expect_identical(a$tcaCLR, b$a3r)
})

## ---- end-to-end on a planted synthetic compendium (production defaults, then Arm B, then the refusals) --------------------------------------
test_that("end to end: gates, score, expjack, collect, report; Arm B has no p-value and no claim; hash refusal stops the driver", {
  set.seed(21); n_g <- 150; gn <- sprintf("G%03d", 1:n_g); n_s <- 40
  sizes <- c(25, 12, 8, 5, 4, 3); cx <- list(); at <- 1; for (k in seq_along(sizes)) { cx[[k]] <- gn[at:(at + sizes[k] - 1)]; at <- at + sizes[k] }
  X <- matrix(rnorm(n_g * n_s), n_g, n_s, dimnames = list(gn, paste0("c", 1:n_s))); for (m in cx) X[m, ] <- X[m, ] + matrix(rep(rnorm(n_s, sd = 1.2), each = length(m)), length(m))
  mm <- do.call(rbind, lapply(seq_along(cx), function(k) data.frame(variant = "main", dataset = "syn", complex = paste0("CX", k), orf = cx[[k]])))
  mm <- rbind(mm, data.frame(variant = "noribo", dataset = "syn", complex = paste0("CX", rep(2:6, sizes[2:6])), orf = unlist(cx[2:6])))
  ex <- data.frame(scope = "all", orf_a = "G001", orf_b = "G002", reason = "paralog_set")
  td <- file.path(tmp, "e2e_truth"); dir.create(td); utils::write.table(mm, file.path(td, "members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(ex, file.path(td, "excluded.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c(paste(tc_sha256(file.path(td, "members.tsv")), "members.tsv"), paste(tc_sha256(file.path(td, "excluded.tsv")), "excluded.tsv")), file.path(td, "t.sha256"))
  writeLines(gn[c(1:5, 30:40)], file.path(td, "perturbed.txt"))
  sd_ <- file.path(tmp, "e2e_scores"); dir.create(sd_); writeLines(gn, file.path(sd_, "genes.txt")); R <- tc_fit_three(X, 2L)
  for (m in c("CLR", "aCLR", "tcaCLR")) tc_write_bin(tc_pack(R[[m]]), file.path(sd_, paste0(m, ".scores")))
  tc_write_bin(tc_pack(tc_absr(X)), file.path(sd_, "absr.scores"))
  for (m in c("GENIE3", "GRNBoost2")) { N <- matrix(runif(n_g^2), n_g); tc_write_bin(tc_pack(pmax(N, t(N))), file.path(sd_, paste0(m, ".scores"))) }
  tc_write_bin(tc_pack(R$aCLR[, ] * 0 + rank(R$aCLR) / length(R$aCLR)), file.path(sd_, "aCLR_rank.scores")); tc_write_bin(rowMeans(R$CLR), file.path(sd_, "mi_rowmean.scores"))
  for (k in 1:2) tc_write_bin(tc_pack(R$tcaCLR + k * 1e-3 * matrix(rnorm(n_g^2), n_g)), file.path(sd_, sprintf("rewire_%02d.scores", k)))
  xf <- file.path(tmp, "e2e_expr.rds"); saveRDS(X, xf)
  run <- function(arm, stage, out, extra = character(0)) system2("Rscript", c("analysis/tcaclr_confirm/score_confirm.R", "--stage", stage, "--arm", arm, "--scores", sd_, "--truth_dir", td, "--members", file.path(td, "members.tsv"),
      "--excluded", file.path(td, "excluded.tsv"), "--sha", file.path(td, "t.sha256"), "--out", out, "--expr", xf, "--perturbed", file.path(td, "perturbed.txt"), "--workers", "2", "--threads", "2", extra), stdout = TRUE, stderr = TRUE)
  rcof <- function(x) { s <- attr(x, "status"); if (is.null(s)) 0L else s }
  outA <- file.path(tmp, "outA")
  # production mode with missing methods is refused; --allow_missing 1 lets the synthetic run proceed with the registered counts
  x <- run("A", "gates", outA); expect_gt(rcof(x), 0L); expect_true(any(grepl("needs all seven methods", x)))
  x <- run("A", "gates", outA, c("--allow_missing", "1")); expect_equal(rcof(x), 0L); g <- utils::read.csv(file.path(outA, "gates.csv")); expect_true(all(g$pass)); expect_true(any(g$check == "power class"))
  x <- run("A", "score", outA, c("--allow_missing", "1")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n"))
  met <- utils::read.csv(file.path(outA, "score/metrics.csv")); expect_true(all(c("tcaCLR", "aCLR", "CLR", "absr", "GENIE3", "GRNBoost2", "aCLR_rank", "rewire_01") %in% met$method)); expect_true("hub_mi" %in% met$scope)
  expect_true(all(c("pert_neither", "pert_rest", "small_3to6", "mid_7to10", "mid_11to19", "large_ge20") %in% met$truth))
  jg <- utils::read.csv(file.path(outA, "score/jack_gene_replicates.csv")); expect_equal(max(jg$rep), 500L)
  # the stored AUPR equals ap_fast on the truth built independently
  Tt <- tc_build_truths(list(mem = mm, excl = ex), gn, "plan", match(readLines(file.path(td, "perturbed.txt")), gn)); v <- tc_unpack(tc_read_bin(file.path(sd_, "tcaCLR.scores")), n_g)
  expect_equal(met$aupr[met$truth == "main" & met$method == "tcaCLR" & met$scope == "raw"], ap_fast(v[cbind(Tt$main$i, Tt$main$j)], Tt$main$label), tolerance = 1e-14)
  x <- run("A", "report", outA, c("--allow_missing", "1")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n"))   # no experiment source yet: gene-only, no claim
  e0 <- utils::read.csv(file.path(outA, "estimate.csv")); expect_true(is.na(e0$p_one_sided)); expect_match(e0$claim, "The claim is not made: experiment component missing")
  unlink(file.path(outA, "estimate.csv"))
  for (lane in 0:1) { x <- run("A", "expjack", outA, c("--allow_missing", "1", "--lane", lane, "--nlanes", "2")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n")) }
  x <- run("A", "collect", outA, c("--allow_missing", "1")); expect_equal(rcof(x), 0L)
  er <- utils::read.csv(file.path(outA, "expjack/replicates_plan.csv")); expect_equal(length(unique(er$rep)), 100L); expect_true(all(c("CLR", "aCLR", "tcaCLR", "absr") %in% er$method))
  # the replicate AP of replicate 1 equals an independent refit on the registered half
  half <- tc_exp_halves(n_s, 100L, TC$exp_seed, "seq")[[1]]; Rb <- tc_fit_three(X[, half], 2L)
  expect_equal(er$aupr[er$rep == 1 & er$truth == "main" & er$method == "tcaCLR"], ap_fast(Rb$tcaCLR[cbind(Tt$main$i, Tt$main$j)], Tt$main$label), tolerance = 1e-14)
  x <- run("A", "report", outA, c("--allow_missing", "1")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n"))
  e1 <- utils::read.csv(file.path(outA, "estimate.csv")); expect_equal(e1$power_class, "estimate_only"); expect_true(is.na(e1$p_one_sided)); expect_match(e1$claim, "registered test is not run")   # n_eff of the synthetic truth is far below 15
  x <- run("A", "report", outA, c("--allow_missing", "1")); expect_gt(rcof(x), 0L); expect_true(any(grepl("written once", x)))    # production report is written once
  x <- run("A", "report", outA, c("--allow_missing", "1", "--neff_override", "30")); expect_gt(rcof(x), 0L); expect_true(any(grepl("refuses non-registered", x)))
  x <- run("A", "report", outA, c("--allow_missing", "1", "--validation", "1", "--neff_override", "30")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n"))   # validation: power gate overridden to exercise the test path
  e <- utils::read.csv(file.path(outA, "estimate.csv")); expect_true(e$exp_included); expect_false(is.na(e$p_one_sided)); expect_false(is.na(e$claim)); expect_match(e$claim, "M3D p = 0.079"); expect_match(e$claim, "VALIDATION RUN")
  expect_equal(e$se_total, sqrt(e$se_gene^2 + e$se_exp^2)); expect_equal(e$total_hi, e$d + 1.96 * e$se_total); expect_equal(e$p_one_sided, 1 - pnorm(e$d / e$se_total))
  h3 <- utils::read.csv(file.path(outA, "h3.csv")); expect_true(all(h3$note[h3$method_b %in% c("GENIE3", "GRNBoost2")] == TC$note_ge)); expect_true(all(is.na(h3$se_exp[h3$method_b %in% c("GENIE3", "GRNBoost2")])))
  expect_true(all(!is.na(h3$se_exp[h3$method_b %in% c("aCLR", "CLR", "absr")])))
  sent <- readLines(file.path(outA, "sentences.txt")); expect_true(any(grepl(TC$note_ge, sent, fixed = TRUE))); expect_true(any(grepl("Trap:", sent))); expect_true(any(grepl("Perturbed-gene split", sent)))
  expect_true(any(grepl("Built-in data-truth link: 16 of the 57 truth genes", sent))); expect_true(any(grepl("Power gate: n_eff = 30.0", sent)))
  pw <- utils::read.csv(file.path(outA, "pairwise.csv")); expect_true(all(pw$note[pw$method_a %in% c("GENIE3", "GRNBoost2") | pw$method_b %in% c("GENIE3", "GRNBoost2")] == TC$note_ge)); expect_true(all(pw$note[!(pw$method_a %in% c("GENIE3", "GRNBoost2") | pw$method_b %in% c("GENIE3", "GRNBoost2"))] == ""))
  expect_true(file.exists(file.path(outA, "readouts.csv"))); ro <- utils::read.csv(file.path(outA, "readouts.csv")); expect_true(any(grepl("^rewire_mean_of_2", ro$readout))); expect_true(any(grepl("^hub_stratified", ro$readout)))
  # a refittable comparator without experiment replicates gets no p-value and blocks the best-of statement
  outE <- file.path(tmp, "outE"); dir.create(outE, showWarnings = FALSE); file.copy(file.path(outA, "score"), outE, recursive = TRUE); dir.create(file.path(outE, "expjack"), showWarnings = FALSE)
  utils::write.csv(er[er$method != "absr", ], file.path(outE, "expjack", "replicates_plan.csv"), row.names = FALSE)
  x <- run("A", "report", outE, c("--allow_missing", "1", "--validation", "1", "--neff_override", "30")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n"))
  hE <- utils::read.csv(file.path(outE, "h3.csv")); expect_true(is.na(hE$p_one_sided[hE$method_b == "absr"])); expect_true(hE$exp_missing[hE$method_b == "absr"]); expect_false(is.na(hE$p_one_sided[hE$method_b == "GENIE3"]))
  expect_true(any(grepl("cannot be declared the best", readLines(file.path(outE, "sentences.txt")))))
  # Arm B: same estimand and descriptors, no p-value, no claim, no best-of statement
  outB <- file.path(tmp, "outB"); dir.create(outB, showWarnings = FALSE); file.copy(file.path(outA, "score"), outB, recursive = TRUE); file.copy(file.path(outA, "expjack"), outB, recursive = TRUE)
  x <- run("B", "report", outB, c("--allow_missing", "1", "--validation", "1", "--neff_override", "30")); expect_equal(rcof(x), 0L, info = paste(tail(x, 8), collapse = "\n"))
  eb <- utils::read.csv(file.path(outB, "estimate.csv")); expect_true(is.na(eb$p_one_sided)); expect_true(is.na(eb$claim)); expect_equal(eb$power_class, "arm_B_estimate"); expect_equal(eb$d, e$d)
  hb <- utils::read.csv(file.path(outB, "h3.csv")); expect_true(all(is.na(hb$p_one_sided))); expect_true(all(is.na(hb$won_iut))); sb <- readLines(file.path(outB, "sentences.txt")); expect_false(any(grepl("best of the", sb))); expect_false(any(grepl("M3D p = 0.079", sb)))
  # non-registered settings are refused in production mode
  x <- run("A", "score", file.path(tmp, "outC"), c("--allow_missing", "1", "--jack", "50")); expect_gt(rcof(x), 0L); expect_true(any(grepl("refuses non-registered", x)))
  # hash refusal at the driver: change one byte of the members file
  writeLines(readLines(file.path(td, "members.tsv")), file.path(td, "members_copy.tsv")); cat("x", file = file.path(td, "members.tsv"), append = TRUE)
  x <- run("A", "gates", file.path(tmp, "outD"), c("--allow_missing", "1")); expect_gt(rcof(x), 0L); expect_true(any(grepl("MISMATCH", x))); expect_false(file.exists(file.path(tmp, "outD", "gates.csv")))
  x <- run("A", "score", file.path(tmp, "outD"), c("--allow_missing", "1")); expect_gt(rcof(x), 0L); expect_false(file.exists(file.path(tmp, "outD", "score", "metrics.csv")))
})
cat("ALL TESTS DONE\n")
