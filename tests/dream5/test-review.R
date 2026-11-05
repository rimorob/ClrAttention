setwd(Sys.getenv("D5_ROOT", unset = "."))
# Tests for the step-6 review fixes (S1 truth edge tables, S2 option validation, N7 ap_expected naming). Decoy and synthetic data only.

test_that("S1: d5_edges_from_ids needs two columns, or three with a 0/1 label; label-0 rows are never positives; self-edges are dropped", {
  nd <- d5_synth_net(40, 6, 12, 3, net = 4L)               # TF IDs G1..G6, genes G1..G40
  two <- data.frame(a = c("G1", "G1", "G2", "G2", "G9", "G3", "G3"), b = c("G10", "G10", "G11", "G2", "G12", "G13", "NOPE"), stringsAsFactors = FALSE)
  e <- d5_edges_from_ids(two, nd)
  # G1->G10 (twice: one kept), G2->G11, G2->G2 self (dropped), G9->G12 (G9 is not a TF: unknown regulator), G3->G13, G3->NOPE (unknown)
  expect_equal(nrow(e), 3); expect_equal(attr(e, "n_self"), 1); expect_equal(attr(e, "n_dup"), 1); expect_equal(attr(e, "n_dropped"), 2); expect_equal(attr(e, "n_label0"), 0)
  expect_equal(sort(paste(nd$tf_ids[e$tf], nd$genes[e$gene])), sort(c("G1 G10", "G2 G11", "G3 G13")))
  # three columns: label 1 kept, label 0 dropped and counted, never positive
  three <- data.frame(a = c("G1", "G2", "G3", "G4"), b = c("G10", "G11", "G12", "G13"), l = c(1, 0, 1, 0))
  e3 <- d5_edges_from_ids(three, nd); expect_equal(nrow(e3), 2); expect_equal(attr(e3, "n_label0"), 2)
  expect_equal(sort(paste(nd$tf_ids[e3$tf], nd$genes[e3$gene])), c("G1 G10", "G3 G12"))
  e3s <- d5_edges_from_ids(transform(three, l = as.character(l)), nd); expect_identical(e3s, e3)       # "0"/"1" as text is read too
  # a label-0-only table yields no positives at all
  e0 <- d5_edges_from_ids(transform(three, l = 0), nd); expect_equal(nrow(e0), 0); expect_equal(attr(e0, "n_label0"), 4)
  # anything else is an error
  expect_error(d5_edges_from_ids(three[, 1, drop = FALSE], nd), "exactly two columns")
  expect_error(d5_edges_from_ids(cbind(three, x = 1), nd), "exactly two columns")
  expect_error(d5_edges_from_ids(transform(three, l = c(1, 2, 0, 1)), nd), "0/1 labels only")
  expect_error(d5_edges_from_ids(transform(three, l = c("1", "x", "0", "1")), nd), "0/1 labels only")
  expect_error(d5_edges_from_ids(transform(three, l = c(1, NA, 0, 1)), nd), "0/1 labels only")
})

test_that("S2: d5_truth_options validates every switch and flags what would be ignored", {
  base <- c("--net", "4", "--truth_edges", "t.tsv", "--out_dir", "o")
  swa <- c("--ties_at_cap", "keep_tied", "--ties_in_list", "end_of_block", "--tftf", "keep_both", "--universe", "all", "--jackknife_cap", "reapply_cap")
  o <- d5_truth_options(c(base, swa)); expect_equal(o$opt$net, 4L); expect_equal(o$sw$jackknife_cap, "reapply_cap"); expect_equal(o$opt$series_n, 100L)
  for (k in c("--ties_at_cap", "--ties_in_list", "--tftf", "--universe", "--jackknife_cap")) {          # each switch is required
    i <- match(k, swa); expect_error(d5_truth_options(c(base, swa[-c(i, i + 1)])), "TRUTH-STAGE switch", info = k)
    bad <- swa; bad[i + 1] <- "nonsense"; expect_error(d5_truth_options(c(base, bad)), "TRUTH-STAGE switch", info = k)
  }
  expect_error(d5_truth_options(c(base, swa, "--cap_order", "cap_then_universe")), "--cap_order was given but --universe is 'all'")
  expect_error(d5_truth_options(c(base, swa, "--universe_tf", "a.txt")), "would be ignored")
  expect_error(d5_truth_options(c(base, swa, "--universe_gene", "a.txt")), "would be ignored")
  gs <- replace(swa, match("all", swa), "gs")
  expect_error(d5_truth_options(c(base, gs)), "TRUTH-STAGE switch 'cap_order'")
  expect_error(d5_truth_options(c(base, gs, "--cap_order", "universe_then_cap")), "needs both --universe_tf and --universe_gene")
  expect_error(d5_truth_options(c(base, gs, "--cap_order", "universe_then_cap", "--universe_tf", "nofile1", "--universe_gene", "nofile2")), "universe file not found")
  f1 <- tempfile(); f2 <- tempfile(); writeLines("G1", f1); writeLines("G2", f2)
  og <- d5_truth_options(c(base, gs, "--cap_order", "universe_then_cap", "--universe_tf", f1, "--universe_gene", f2)); expect_equal(og$sw$cap_order, "universe_then_cap")
  nd <- d5_synth_net(40, 6, 12, 3, net = 4L); u <- d5_universe_from_files(og$opt, nd); expect_equal(u$tf, 1); expect_equal(u$gene, 2)
  writeLines("ZZZ", f2); expect_error(d5_universe_from_files(og$opt, nd), "not in the gene list")
  expect_null(d5_universe_from_files(o$opt, nd))
  # other options, and jackknife_cap is checked up front (not only when the jackknife starts)
  expect_error(d5_truth_options(c(base, swa, "--nrep", "1")), "at least 2"); expect_error(d5_truth_options(c(base, swa, "--series_n", "101")), "series_n")
  expect_error(d5_truth_options(c(base, swa, "--series_n", "x")), "must be an integer"); expect_error(d5_truth_options(c(base, swa, "--workers")), "no value")
  expect_error(d5_truth_options(c(base, swa, "--bogus", "1")), "unknown argument"); expect_error(d5_truth_options(c(base, swa, "--cap", "0")), "positive")
  expect_error(d5_truth_options(c(base, swa, "--structure_only", "2")), "0 or 1"); expect_error(d5_truth_options(c(base, swa, "--nrep", "5", "--nrep", "6")), "twice")
  expect_equal(d5_truth_options(c(base, swa, "--series_n", "0"))$opt$series_n, 0L)
  expect_error(d5_truth_options(c("--net", "3", "--truth_edges", "t", "--out_dir", "o", swa, "--workers", "8")), "would be ignored")
  expect_error(d5_truth_options(c("--net", "2", "--truth_edges", "t", "--out_dir", "o", swa)), "--net must be")
  expect_equal(d5_truth_options(c("--net", "3", "--truth_edges", "t", "--out_dir", "o", swa))$opt$net, 3L)
})

test_that("N7: the dream5 expected-tie AP is d5_ap_expected; ap_expected stays the regulons.R definition", {
  expect_false(any(grepl("^ap_expected <-", readLines("analysis/dream5_score_lib.R"))))
  e <- new.env(); sys.source("analysis/regulons.R", envir = e)
  expect_true(identical(deparse(ap_expected), deparse(e$ap_expected)))                  # not shadowed after both libraries are sourced
  expect_true(exists("d5_ap_expected"))
})
