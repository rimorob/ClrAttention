# Official Y-E metrics (trapezoid AUPR/(1-1/P), AUROC, expected-fraction tail), p-value lookup, gold-derived universe, run kinds, scorer gate.
# Only synthetic data and the organisers' shipped SAMPLE files (their demo predictions 'myteam' and OUTPUT/myteam.txt) are used; no score of ours.
setwd(Sys.getenv("D5_ROOT", unset = "."))   # repo root: testthat runs each file in tests/dream5
source("analysis/dream5_official.R")
SCR <- "data/dream5/truth_extracted/Evaluation scripts"
have_demo <- file.exists(file.path(SCR, "OUTPUT", "myteam.txt")) && file.exists(file.path(D5$densities_dir, "Network4_AUROC.tsv"))
have_py <- nzchar(Sys.which("python3")) && system2("python3", c("-I", "-c", shQuote("import scipy, numpy")), stdout = FALSE, stderr = FALSE) == 0L

perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(r) c(v[i], r))))
# every within-block order of a listed vector (blocks = equal scores, s sorted descending)
all_orders <- function(s) {
  bl <- split(seq_along(s), factor(s, levels = unique(s))); combos <- list(integer(0))
  for (b in bl) { pp <- perms(b); combos <- do.call(c, lapply(combos, function(c0) lapply(pp, function(p) c(c0, p)))) }
  combos
}

test_that("official trapezoid: literal transcription equals the MATLAB loop forms; hand-computed case", {
  # hand case: P = 2, T = 4, list = (1, 0); tail has 2 slots with p = (2 - 1)/(4 - 2) = 0.5
  r <- d5_official_literal(c(1, 0), 2, 4)
  TPk <- c(1, 1, 1.5, 2); FPk <- c(0, 1, 1.5, 2); TPR <- TPk / 2; FPR <- FPk / 2; PREC <- TPk / 1:4
  auroc <- 0; for (i in 1:3) auroc <- auroc + (FPR[i + 1] - FPR[i]) * (TPR[i + 1] + TPR[i]) / 2        # the commented loop of the official script
  aupr <- 0;  for (i in 1:3) aupr <- aupr + (TPR[i + 1] - TPR[i]) * (PREC[i + 1] + PREC[i]) / 2
  expect_equal(unname(r["auroc"]), auroc, tolerance = 1e-14); expect_equal(unname(r["aupr"]), aupr / (1 - 1 / 2), tolerance = 1e-14)
  # random larger lists against the loop form
  set.seed(5); disc <- rbinom(60, 1, .3); P <- sum(disc) + 7; T <- 200
  r <- d5_official_literal(disc, P, T); p <- 7 / 140; pos <- c(disc, rep(p, 140)); neg <- 1 - pos
  TPk <- cumsum(pos); FPk <- cumsum(neg); TPR <- TPk / P; FPR <- FPk / (T - P); PREC <- TPk / 1:T
  lo <- 0; for (i in 1:(T - 1)) lo <- lo + (FPR[i + 1] - FPR[i]) * (TPR[i + 1] + TPR[i]) / 2
  lp <- 0; for (i in 1:(T - 1)) lp <- lp + (TPR[i + 1] - TPR[i]) * (PREC[i + 1] + PREC[i]) / 2
  expect_equal(unname(r["auroc"]), lo, tolerance = 1e-12); expect_equal(unname(r["aupr"]), lp / (1 - 1 / P), tolerance = 1e-12)
})

test_that("expected official statistics: no ties equal the literal list; exact enumeration of all within-block orders", {
  set.seed(11)
  for (rep in 1:25) {                                           # no ties: strictly decreasing scores
    K <- sample(0:150, 1); T <- K + sample(5:400, 1); y <- rbinom(K, 1, .2); P <- sum(y) + sample(2:30, 1); P <- min(P, T - 1)
    if (P < sum(y)) next
    expect_equal(d5_official_expected(rev(seq_len(K)), y, P, T), d5_official_literal(y, P, T), tolerance = 1e-9, info = rep)
  }
  maxerr <- 0; ncase <- 0
  for (rep in 1:80) {                                           # ties: exhaustive enumeration (K <= 8, so <= 8! orders)
    K <- sample(0:8, 1); T <- K + sample(0:6, 1); if (T < 4) next
    P <- sample(2:(T - 1), 1); s <- sort(sample(1:3, K, replace = TRUE), decreasing = TRUE); y <- rbinom(K, 1, .5)
    if (sum(y) > P) next
    ex <- rowMeans(sapply(all_orders(s), function(idx) d5_official_literal(y[idx], P, T)))
    maxerr <- max(maxerr, abs(ex - d5_official_expected(s, y, P, T))); ncase <- ncase + 1
  }
  expect_gt(ncase, 40); expect_lt(maxerr, 1e-12)
  # a case with one big tie covering the whole list and one with the tail empty
  s <- rep(1, 7); y <- c(1, 0, 1, 1, 0, 0, 1)
  ex <- rowMeans(sapply(all_orders(s), function(idx) d5_official_literal(y[idx], 4, 7))); expect_equal(d5_official_expected(s, y, 4, 7), ex, tolerance = 1e-12)
  expect_true(all(is.na(d5_official_expected(c(2, 1), c(1, 0), 1, 10))))     # P < 2: AUPR normalisation undefined
})

test_that("expected official statistics: Monte Carlo on a mid-size case with large tied blocks and a tail", {
  set.seed(2026)
  K <- 600; s <- sort(sample(seq_len(25), K, replace = TRUE), decreasing = TRUE); y <- rbinom(K, 1, ifelse(s > 18, .25, .05)); P <- sum(y) + 40; T <- 4000
  blocks <- split(seq_len(K), factor(s, levels = unique(s)))
  draw <- function() { idx <- unlist(lapply(blocks, function(b) b[sample.int(length(b))]), use.names = FALSE); d5_official_literal(y[idx], P, T) }
  sims <- t(replicate(2500, draw())); ex <- d5_official_expected(s, y, P, T)
  for (m in c("aupr", "auroc")) {
    se <- stats::sd(sims[, m]) / sqrt(nrow(sims)); expect_lt(abs(mean(sims[, m]) - ex[[m]]), 4 * se, label = m)
  }
  expect_gt(max(lengths(blocks)), 25)                                         # genuinely tied blocks
})

test_that("official AUROC with expected ties equals the mid-rank AUROC of d5_list_tail_metrics; official AUPR differs from step-wise AP", {
  set.seed(8); s <- sort(sample(1:6, 40, replace = TRUE), decreasing = TRUE); y <- rbinom(40, 1, .3); P <- sum(y) + 5; T <- 300
  o <- d5_official_expected(s, y, P, T); st <- d5_list_tail_metrics(s, y, P, T, "expected")
  expect_equal(unname(o["auroc"]), unname(st["auroc"]), tolerance = 1e-12)
  expect_gt(abs(unname(o["aupr"]) - unname(st["aupr"])), 1e-6)
})

test_that("d5_edge_metrics: official fields, cap-boundary tie count, drop_tied fix, readout selection", {
  set.seed(31); nTF <- 12; G <- 80; tf_idx <- 1:12
  tr <- d5_decoy_truth(nTF, G, tf_idx, 7, median_out = 6); T <- d5_truth_matrix(tr, nTF, G)
  M <- matrix(rnorm(nTF * G), nTF, G); M[cbind(1:nTF, tf_idx)] <- NA
  u <- list(tf = 1:7, gene = 10:60); sw <- do.call(d5_edge_switches, D5$ye_primary); cap <- 150
  r <- d5_edge_metrics(M, T, tf_idx, sw, cap, universe = u)
  # independent construction: cap on all candidates, then restrict to the universe, official literal on the kept list
  o <- order(M, decreasing = TRUE, na.last = NA)[seq_len(cap)]; inu <- (row(M)[o] %in% u$tf) & (col(M)[o] %in% u$gene)
  uni <- matrix(FALSE, nTF, G); uni[u$tf, u$gene] <- TRUE; uni[cbind(1:nTF, tf_idx)] <- FALSE
  lit <- d5_official_literal(as.numeric(T[o][inu]), sum(T[uni]), sum(uni))
  expect_equal(r[["aupr_official"]], lit[["aupr"]], tolerance = 1e-12); expect_equal(r[["auroc_official"]], lit[["auroc"]], tolerance = 1e-12)
  expect_equal(r[["auroc_official"]], r[["auroc"]], tolerance = 1e-12)         # AUROC: official trapezoid = mid-rank AUROC
  expect_equal(r[["n_tied_at_cap"]], 0); expect_equal(r[["N"]], sum(uni)); expect_equal(r[["K"]], sum(inu))
  # end_of_block has no official counterpart
  e <- d5_edge_metrics(M, T, tf_idx, d5_edge_switches("drop_tied", "end_of_block", "keep_both", "all"), cap)
  expect_true(is.na(e[["aupr_official"]])); expect_error(.d5_pick(e, "official"), "ties_in_list")
  expect_equal(.d5_pick(r, "official"), r[["aupr_official"]]); expect_equal(.d5_pick(r, "stepwise"), r[["aupr"]])
  # cap boundary: a block of 4 tied scores straddling the cap is counted; a unique cap-th score is not a tie
  M2 <- M; ord <- order(M2, decreasing = TRUE, na.last = NA); M2[ord[148:151]] <- M2[ord[148]]
  swa <- d5_edge_switches("drop_tied", "expected", "keep_both", "all")
  d <- d5_edge_metrics(M2, T, tf_idx, swa, 150); k <- d5_edge_metrics(M2, T, tf_idx, d5_edge_switches("keep_tied", "expected", "keep_both", "all"), 150)
  expect_equal(d[["n_tied_at_cap"]], 4); expect_equal(c(d[["K"]], k[["K"]]), c(147, 151))
  expect_equal(d5_edge_metrics(M, T, tf_idx, swa, 150)[["K"]], 150)             # drop_tied does not drop a lone cap-th entry (no straddling block)
  M4 <- M; M4[ord[149]] <- M4[ord[150]]                                         # a block that ENDS exactly at the cap does not straddle it
  d4 <- d5_edge_metrics(M4, T, tf_idx, swa, 150); expect_equal(c(d4[["n_tied_at_cap"]], d4[["K"]]), c(0, 150))
  # nothing cut: no boundary
  expect_equal(d5_edge_metrics(M, T, tf_idx, swa, 1e6)[["n_tied_at_cap"]], 0)
})

test_that("official p-value lookup, overall score and the shipped OUTPUT/myteam.txt (printed precision), against the Python port", {
  skip_if_not(have_demo)
  files <- function(n) c(gold = sprintf("%s/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network%d.tsv", SCR, n),
                         pred = sprintf("%s/INPUT/predictions/myteam/DREAM5_NetworkInference_myteam_Network%d.txt", SCR, n))
  r <- d5_off_score_team(files, D5$densities_dir); tb <- r$table
  ref <- readLines(file.path(SCR, "OUTPUT", "myteam.txt")); num <- function(i) as.numeric(strsplit(trimws(ref[i]), "\t")[[1]])
  expect_equal(round(r$overall[["overall"]], 3), num(2)[1]); expect_equal(round(r$overall[["aupr_score"]], 3), num(2)[2]); expect_equal(round(r$overall[["auroc_score"]], 3), num(2)[3])
  expect_equal(as.numeric(sprintf("%.3e", tb$P_AUPR)), num(5)); expect_equal(as.numeric(sprintf("%.3e", tb$P_AUROC)), num(8))
  expect_equal(round(tb$AUPR, 3), num(11)); expect_equal(round(tb$AUROC, 3), num(14))
  # the values quoted in the truth-stage log (6 digits)
  expect_equal(round(tb$AUPR, 6), c(0.027128, 0.027424, 0.019338)); expect_equal(round(tb$AUROC, 6), c(0.681817, 0.599660, 0.515037))
  expect_equal(tb$P, c(4012, 2066, 3940)); expect_equal(tb$T, c(278392, 152280, 227202)); expect_equal(tb$L, c(86225, 12032, 12305))
  # p-value lookup: below the support the p is the total mass (<1), above the support top it is exactly 0
  dn <- d5_read_density(4, "AUPR"); expect_lt(abs(d5_density_p(dn, -1) - sum(dn$Y) * (dn$X[2] - dn$X[1])), 1e-12)
  expect_equal(d5_density_p(dn, 0.1218), 0); expect_gt(d5_density_p(dn, 0.12), 0)   # the density underflows to 0 just below its support top (official behaviour)
  expect_equal(d5_density_p(d5_read_density(4, "AUROC"), 0.6175), 0)
  expect_true(all(diff(d5_density_p(dn, seq(0, 0.2, by = 0.01))) <= 0))
  # overall score: an exact p of 0 gives Inf, never capped
  o <- d5_official_overall(c(0.5, 0, 0.1), c(0.2, 0.2, 0.2)); expect_true(is.infinite(o[["overall"]]) && is.infinite(o[["aupr_score"]]) && is.finite(o[["auroc_score"]]))
  expect_equal(unname(d5_official_overall(c(.1, .01, .001), c(.1, .1, .1))), c(mean(c(2, 1)), 2, 1))
  skip_if_not(have_py)
  pf <- tempfile(fileext = ".tsv")
  system2("python3", c("-I", "analysis/dream5_official_tools.py", "port-run", shQuote("data/dream5/truth_extracted_tools/d5_official_port.py"),
          shQuote(file.path(SCR, "INPUT/gold_standard_edges_only")), shQuote(file.path(SCR, "INPUT/probability_densities")), shQuote(file.path(SCR, "INPUT/predictions")), pf), stdout = FALSE)
  pp <- utils::read.delim(pf)
  for (m in c("AUPR", "AUROC", "P_AUPR", "P_AUROC")) expect_lt(max(abs(pp[[m]] - tb[[m]])), 1e-9, label = m)
  # the matrix pipeline code on the demo lists (strictly decreasing scores): equal to the port to 1e-9
  for (n in c(3L, 4L)) {
    nd <- d5_read_net(n); gold <- utils::read.delim(files(n)[["gold"]], header = FALSE, stringsAsFactors = FALSE)
    mp <- d5_off_matrix_path(d5_off_read_edges(files(n)[["pred"]]), gold, nd); i <- match(n, tb$net)
    expect_lt(abs(mp[["AUPR"]] - pp$AUPR[i]), 1e-9); expect_lt(abs(mp[["AUROC"]] - pp$AUROC[i]), 1e-9); expect_equal(mp[["N"]], pp$T[i]); expect_equal(mp[["n_pos"]], pp$P[i])
  }
})

test_that("gold-derived universe: function output equals the extraction's universe files; sizes match the official T", {
  skip_if_not(file.exists(file.path(SCR, "INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network4.tsv")) && file.exists("data/dream5/truth_extracted_universe/net4_universe_tf.txt"))
  for (n in c(1L, 3L, 4L)) {
    nd <- d5_read_net(n); gold <- utils::read.delim(sprintf("%s/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network%d.tsv", SCR, n), header = FALSE, stringsAsFactors = FALSE)
    edges <- d5_edges_from_ids(gold, nd); u <- d5_gold_universe(edges, nd)
    ft <- sprintf("data/dream5/truth_extracted_universe/net%d_universe_tf.txt", n); fg <- sprintf("data/dream5/truth_extracted_universe/net%d_universe_gene.txt", n)
    expect_true(d5_check_universe_files(u, ft, fg)); expect_setequal(u$tf_ids, readLines(ft)); expect_setequal(u$gene_ids, readLines(fg))
    expect_equal(nrow(edges), c(`1` = 4012, `3` = 2066, `4` = 3940)[[as.character(n)]]); expect_equal(attr(edges, "n_dropped"), 0)
    expect_equal(c(length(u$tf), length(u$gene)), list(`1` = c(178, 1565), `3` = c(141, 1081), `4` = c(114, 1994))[[as.character(n)]])
    nself <- length(intersect(nd$tf_idx[u$tf], u$gene)); T <- length(u$tf) * length(u$gene) - nself
    expect_equal(T, c(`1` = 278392, `3` = 152280, `4` = 227202)[[as.character(n)]])
    # the universe the metrics see: sum of the masked valid cells
    expect_equal(sum(outer(seq_along(nd$tf_idx), seq_along(nd$genes), function(i, j) i %in% u$tf & j %in% u$gene & nd$tf_idx[i] != j)), T)
  }
  tmp <- tempfile(); fs <- d5_write_universe(u, tmp); expect_true(d5_check_universe_files(u, fs[1], fs[2]))
  writeLines(c(readLines(fs[1]), "G1"), fs[1]); expect_error(d5_check_universe_files(u, fs[1], fs[2]), "does not match|G1")
  writeLines(readLines(fs[2])[-1], fs[2]); writeLines(u$tf_ids, fs[1]); expect_error(d5_check_universe_files(u, fs[1], fs[2]), "does not match the genes")
})

test_that("driver options: the Y-E primary run allows only the design's decisions; sensitivity runs are labelled; official readout needs expected ties", {
  base <- c("--net", "4", "--truth_edges", "t.tsv", "--out_dir", "o")
  prim <- c("--run_kind", "primary", "--universe", "gs", "--cap_order", "cap_then_universe", "--tftf", "keep_both", "--ties_at_cap", "drop_tied",
            "--ties_in_list", "expected", "--jackknife_cap", "restrict_after_cap")
  o <- d5_truth_options(c(base, prim)); expect_equal(o$opt$run_kind, "primary"); expect_equal(o$opt$readout, "official")
  for (k in names(D5$ye_primary)) expect_identical(o$sw[[k]], D5$ye_primary[[k]])
  alt <- list(`--universe` = "all", `--cap_order` = "universe_then_cap", `--tftf` = "collapse_max", `--ties_at_cap` = "keep_tied", `--ties_in_list` = "end_of_block", `--jackknife_cap` = "reapply_cap")
  for (k in names(alt)) {
    a <- prim; a[match(k, a) + 1L] <- alt[[k]]
    if (k == "--universe") a <- a[!a %in% c("--cap_order", "cap_then_universe")]
    expect_error(d5_truth_options(c(base, a)), "allows only the decisions of the design change log 2026-11-06", info = k)
  }
  expect_error(d5_truth_options(c(base, prim, "--readout", "stepwise")), "official readout only")
  expect_equal(d5_truth_options(c(base, prim, "--readout", "official"))$opt$readout, "official")
  expect_error(d5_truth_options(c(base, prim, "--universe_tf", "a.txt")), "both --universe_tf and --universe_gene")
  expect_error(d5_truth_options(c(base, prim, "--universe_tf", "nofile", "--universe_gene", "nofile2")), "universe file not found")
  # the switches are still required in a primary run
  expect_error(d5_truth_options(c(base, prim[-(11:12)])), "TRUTH-STAGE switch")
  swa <- c("--ties_at_cap", "keep_tied", "--ties_in_list", "end_of_block", "--tftf", "keep_both", "--universe", "all", "--jackknife_cap", "reapply_cap")
  expect_error(d5_truth_options(c(base, swa, "--run_kind", "bogus")), "primary or sensitivity")
  # sensitivity (default): any switch values; stepwise readout by default; official only with expected ties
  s <- d5_truth_options(c(base, swa)); expect_equal(c(s$opt$run_kind, s$opt$readout), c("sensitivity", "stepwise"))
  expect_error(d5_truth_options(c(base, swa, "--readout", "official")), "needs --ties_in_list expected")
  swb <- swa; swb[4] <- "expected"; expect_equal(d5_truth_options(c(base, swb, "--readout", "official"))$opt$readout, "official")
  expect_error(d5_truth_options(c(base, swa, "--readout", "nonsense")), "--readout must be one of")
  # the driver script wires the primary pieces
  src <- paste(readLines("analysis/dream5_truth_run.R"), collapse = "\n")
  expect_match(src, "d5_gold_universe\\(edges, nd\\)"); expect_match(src, "metric = metric"); expect_match(src, "CAP-BOUNDARY CHECK"); expect_match(src, "cap_boundary")
})

test_that("pipeline with the official readout on a synthetic net: labelled readout, cap-boundary table, official p-values, stepwise sensitivity column", {
  skip_if_not(file.exists(file.path(D5$densities_dir, "Network3_AUPR.tsv")))
  old <- D5$ntrees; D5$ntrees <- 20L; on.exit(D5$ntrees <- old, add = TRUE)
  nd <- d5_synth_net(300, 40, 24, 6, net = 3L); out <- tempfile(); meths <- c("clr", "attention")
  for (m in meths) suppressMessages(d5_run_method(nd, m, "replicate", 2L, out))
  edges <- d5_decoy_truth(40, 300, nd$tf_idx, 91, median_out = 12, sd_log = 0.4)
  u <- d5_gold_universe(edges, nd); sw <- do.call(d5_edge_switches, D5$ye_primary)
  expect_error(d5_truth_pipeline(nd, edges, TESTSW(), methods = meths, out = out, nrep = 4L, cap = 600L, metric = "official", verbose = FALSE), "ties_in_list")
  r <- d5_truth_pipeline(nd, edges, sw, methods = meths, out = out, nrep = 6L, cap = 600L, universe = u, metric = "official", densities = D5$densities_dir, verbose = FALSE)
  expect_equal(r$metric, "official"); expect_true(all(r$point$readout[r$point$endpoint == "YE"] == "official"))
  expect_true(all(r$point$aupr_stepwise[r$point$endpoint == "YE"] > 0)); expect_false(isTRUE(all.equal(r$point$aupr_stepwise[1:2], r$point$aupr[1:2])))
  expect_equal(r$cap_boundary$n_tied_at_cap, c(0, 0)); expect_equal(r$cap_boundary$method, meths)
  expect_equal(nrow(r$official), 2); expect_true(all(r$official$P_AUPR >= 0 & r$official$P_AUPR <= 1)); expect_equal(r$official$neglog10_P_AUPR, -log10(r$official$P_AUPR))
  expect_match(r$rows$YE$readout, "official Y-E")
  dd <- tempfile(); d5_write_results(r, dd); expect_true(all(file.exists(file.path(dd, c("cap_boundary.csv", "official_scores.csv", "verdict.csv")))))
  rs <- d5_truth_pipeline(nd, edges, sw, methods = meths, out = out, nrep = 6L, cap = 600L, universe = u, metric = "stepwise", verbose = FALSE)
  expect_match(rs$rows$YE$readout, "SENSITIVITY ONLY"); expect_equal(rs$point$aupr[1:2], r$point$aupr_stepwise[1:2]); expect_null(rs$official)
})

PUBTAB <- "analysis/dream5_published_table.py"
# a synthetic export of the organisers' sheet: same layout as r-survey4.xls.csv (title rows, category rows, rows with the within-category number in column 2)
synth_sheet <- function(path, drop = NULL, counts = c(Regression = 8, `Mutual information` = 5, Correlation = 3, `Bayesian networks` = 6, Other = 8, Meta = 5), bad = NULL) {
  ln <- c(strrep(",", 46), ",Method,Score,,,AUPR,,,AUROC,,,\"Net2 (S. aureus, not used for scoring)\",,,,\"Survey\"", ",,Overall,AUPR,AUROC,Net1,Net3,Net4,Net1,Net3,Net4,AUPR,AUROC,,,1: yes")
  k <- 0
  for (cat in names(counts)) { ln <- c(ln, sprintf(",%s,,,,,,,,,,,,,%s", cat, cat))
    for (i in seq_len(counts[[cat]])) { k <- k + 1; v <- c(k + 0.25, k + 0.2, k + 0.3, k / 100, k / 200, k / 300, 0.5 + k / 100, 0.6, 0.7, 0.08, 0.55)
      if (identical(bad, c(cat, i))) v[3] <- NA
      ln <- c(ln, paste0(",", i, ",", paste(ifelse(is.na(v), "", v), collapse = ","), ",,", i, ", ,*")) } }
  ln <- c(ln, ",Community,,,,,,,,,,,,,Community", ",RankSum,52.4271,58.5993,46.2549,0.327,0.09,0.022,0.809,0.65,0.528,0.0971,0.5845", strrep(",", 12))
  if (!is.null(drop)) ln <- ln[-grep(drop, ln)[1]]
  writeLines(ln, path); path
}
test_that("published-table converter: tidy output, team naming, row-count assertions fail loudly", {
  skip_if_not(nzchar(Sys.which("python3")))
  td <- tempfile(); dir.create(td); src <- synth_sheet(file.path(td, "s.csv")); out <- file.path(td, "t.csv")
  st <- system2("python3", c("-I", PUBTAB, "tidy", shQuote(src), shQuote(out)), stdout = TRUE, stderr = TRUE); expect_null(attr(st, "status"))
  tb <- utils::read.csv(out, stringsAsFactors = FALSE)
  expect_equal(names(tb), c("team", "score_overall", "score_aupr", "score_auroc", "aupr1", "aupr3", "aupr4", "auroc1", "auroc3", "auroc4")); expect_equal(nrow(tb), 36L)
  expect_equal(tb$team[c(1, 8, 9, 13, 14, 16, 17, 22, 36)], c("Regression1", "Regression8", "MI1", "MI5", "Correlation1", "Correlation3", "Bayesian1", "Bayesian6", "Community"))
  expect_equal(tb$team[nrow(tb)], "Community"); expect_equal(sum(grepl("^Other", tb$team)), 8L); expect_equal(sum(grepl("^Meta", tb$team)), 5L)
  expect_equal(tb$aupr1[tb$team == "Regression3"], 0.03); expect_equal(tb$score_overall[tb$team == "Community"], 52.4271); expect_equal(tb$auroc4[tb$team == "Community"], 0.528)
  # a missing row, an extra row, a missing number: nonzero status and a clear message; nothing is written
  chk <- function(srcf, pat) { o <- file.path(td, "x.csv"); r <- suppressWarnings(system2("python3", c("-I", PUBTAB, "tidy", shQuote(srcf), shQuote(o)), stdout = TRUE, stderr = TRUE))
    expect_true(!is.null(attr(r, "status")) && attr(r, "status") != 0L); expect_match(paste(r, collapse = " "), pat); expect_false(file.exists(o)) }
  chk(synth_sheet(file.path(td, "d.csv"), drop = "^,Community,,"), "row label|unexpected row counts")
  chk(synth_sheet(file.path(td, "d2.csv"), drop = "^,6,"), "unexpected row counts")
  chk(synth_sheet(file.path(td, "e.csv"), counts = c(Regression = 9, `Mutual information` = 5, Correlation = 3, `Bayesian networks` = 6, Other = 8, Meta = 5)), "unexpected row counts")
  chk(synth_sheet(file.path(td, "m.csv"), bad = c("Other", 2)), "column 4 \\(score_auroc\\) is empty")
})
test_that("published-table converter on the real sheet (when present): 36 named rows, values as printed", {
  f <- "data/dream5/published_extracted/results/r-survey4.xls.csv"; skip_if_not(file.exists(f) && nzchar(Sys.which("python3")))
  out <- tempfile(fileext = ".csv"); st <- system2("python3", c("-I", PUBTAB, "tidy", shQuote(f), shQuote(out)), stdout = TRUE, stderr = TRUE); expect_null(attr(st, "status"))
  tb <- utils::read.csv(out, stringsAsFactors = FALSE); expect_equal(nrow(tb), 36L); expect_equal(anyDuplicated(tb$team), 0L)
  expect_true(all(c("Regression8", "MI5", "Correlation3", "Bayesian6", "Other8", "Meta5", "Community") %in% tb$team))
  expect_equal(abs(tb$score_overall - (tb$score_aupr + tb$score_auroc) / 2) <= 2e-3, rep(TRUE, 36))
})
test_that("flat prediction directory: symlinks for all folders, files untouched, refuses incomplete teams", {
  skip_if_not(nzchar(Sys.which("python3")))
  td <- tempfile(); root <- file.path(td, "Network_predictions")
  mk <- function(sub, t) { d <- file.path(root, sub, t); dir.create(d, recursive = TRUE); for (n in 1:4) writeLines(paste("G1", "G2", n), file.path(d, sprintf("DREAM5_NetworkInference_%s_Network%d.txt", t, n))) }
  mk("Challenge participants", "Regression1"); mk("Challenge participants", "Meta2"); mk("Off-the-shelf methods", "MI1")
  dir.create(file.path(root, "Community integration"), recursive = TRUE); for (n in 1:4) writeLines("x", file.path(root, "Community integration", sprintf("DREAM5_NetworkInference_Community_Network%d.txt", n)))
  before <- tools::md5sum(list.files(root, recursive = TRUE, full.names = TRUE))
  flat <- file.path(td, "flat"); st <- system2("python3", c("-I", PUBTAB, "flat-pred", shQuote(root), shQuote(flat)), stdout = TRUE, stderr = TRUE); expect_null(attr(st, "status"))
  expect_equal(sort(list.files(flat)), c("Community", "MI1", "Meta2", "Regression1"))
  for (t in c("Community", "MI1", "Meta2", "Regression1")) { f <- file.path(flat, t, sprintf("DREAM5_NetworkInference_%s_Network3.txt", t)); expect_true(nzchar(Sys.readlink(f))); expect_true(file.exists(f)) }
  expect_equal(readLines(file.path(flat, "MI1", "DREAM5_NetworkInference_MI1_Network3.txt")), "G1 G2 3")
  expect_equal(tools::md5sum(list.files(root, recursive = TRUE, full.names = TRUE)), before)          # the extracted files are never modified
  st2 <- system2("python3", c("-I", PUBTAB, "flat-pred", shQuote(root), shQuote(flat)), stdout = TRUE, stderr = TRUE); expect_null(attr(st2, "status"))   # idempotent
  unlink(file.path(root, "Challenge participants", "Meta2", "DREAM5_NetworkInference_Meta2_Network4.txt"))
  st3 <- suppressWarnings(system2("python3", c("-I", PUBTAB, "flat-pred", shQuote(root), shQuote(file.path(td, "flat2"))), stdout = TRUE, stderr = TRUE)); expect_true(!is.null(attr(st3, "status"))); expect_match(paste(st3, collapse = " "), "missing prediction file")
})

test_that("scorer gate: refuses without both data inputs; passes on the shipped demo (three renamed copies); fails on a wrong published table; R, port and matrix path agree", {
  skip_if_not(have_demo && have_py)
  source("analysis/dream5_scorer_gate.R", local = TRUE)
  o0 <- d5_gate_args(c("--out_dir", tempfile()))
  r <- d5_gate_run(o0, say = function(...) NULL); expect_equal(r$status, 2L); expect_match(r$msg, "REFUSED.*--pred_dir.*--results_csv")
  tdir <- tempfile(); pred <- file.path(tdir, "pred"); teams <- c("myteam", "teamb", "teamc")
  for (t in teams) { dir.create(file.path(pred, t), recursive = TRUE)
    for (n in c(1, 3, 4)) file.copy(sprintf("%s/INPUT/predictions/myteam/DREAM5_NetworkInference_myteam_Network%d.txt", SCR, n), file.path(pred, t, sprintf("DREAM5_NetworkInference_%s_Network%d.txt", t, n))) }
  ref <- readLines(file.path(SCR, "OUTPUT", "myteam.txt")); num <- function(i) as.numeric(strsplit(trimws(ref[i]), "\t")[[1]])
  # the shipped demo output: line 2 = overall / AUPR score / AUROC score, line 11 = AUPR nets, line 14 = AUROC nets
  mk <- function(f, fudge = 0, score_fudge = 0, names = teams) { tb <- data.frame(team = names, score_overall = num(2)[1] + score_fudge, score_aupr = num(2)[2], score_auroc = num(2)[3], check.names = FALSE)
    for (k in 1:3) { n <- c(1, 3, 4)[k]; tb[[sprintf("aupr%d", n)]] <- num(11)[k] + fudge; tb[[sprintf("auroc%d", n)]] <- num(14)[k] }
    utils::write.csv(tb, f, row.names = FALSE) }
  csv <- file.path(tdir, "results.csv"); mk(csv)
  o <- d5_gate_args(c("--pred_dir", pred, "--results_csv", csv, "--out_dir", file.path(tdir, "out")))
  res <- d5_gate_run(o, say = function(...) NULL)
  expect_equal(res$status, 0L); expect_true(res$pass); expect_equal(sort(res$ok_teams), sort(teams)); expect_lt(res$d_rp, 1e-9); expect_lt(res$d_rm, 1e-9)
  expect_true(all(c("R", "port", "matrix") %in% res$comparison$impl)); expect_true(all(res$comparison$pass)); expect_lte(max(res$comparison$abs_diff), 1e-3)
  expect_setequal(unique(res$comparison$quantity), c("AUPR", "AUROC", "score_overall", "score_aupr", "score_auroc")); expect_equal(nrow(res$comparison), 3 * 3 * (6 + 3))
  expect_true(all(file.exists(file.path(tdir, "out", c("gate_report.md", "gate_comparison.csv"))))); expect_match(readLines(file.path(tdir, "out", "gate_report.md"))[1], "PASS")
  # fewer teams than required, a wrong published AUPR, a wrong published aggregate score, a missing column, a name that matches no directory
  o2 <- d5_gate_args(c("--pred_dir", pred, "--results_csv", csv, "--out_dir", file.path(tdir, "o2"), "--teams", "myteam", "--matrix", "0"))
  expect_equal(d5_gate_run(o2, say = function(...) NULL)$status, 1L)                                     # one team < min_teams 3
  csv2 <- file.path(tdir, "bad.csv"); mk(csv2, fudge = 0.01)
  r3 <- d5_gate_run(d5_gate_args(c("--pred_dir", pred, "--results_csv", csv2, "--out_dir", file.path(tdir, "o3"), "--matrix", "0")), say = function(...) NULL); expect_equal(r3$status, 1L); expect_length(r3$ok_teams, 0)
  csv2b <- file.path(tdir, "bad2.csv"); mk(csv2b, score_fudge = 0.01)
  r3b <- d5_gate_run(d5_gate_args(c("--pred_dir", pred, "--results_csv", csv2b, "--out_dir", file.path(tdir, "o3b"), "--matrix", "0")), say = function(...) NULL); expect_equal(r3b$status, 1L)
  expect_true(all(r3b$comparison$pass[r3b$comparison$quantity != "score_overall"])); expect_false(any(r3b$comparison$pass[r3b$comparison$quantity == "score_overall"]))
  tb <- utils::read.csv(csv, check.names = FALSE); utils::write.csv(tb[, !grepl("auroc3", names(tb))], file.path(tdir, "nocol.csv"), row.names = FALSE)
  expect_error(d5_gate_run(d5_gate_args(c("--pred_dir", pred, "--results_csv", file.path(tdir, "nocol.csv"), "--out_dir", file.path(tdir, "o4"))), say = function(...) NULL), "lacks the columns auroc3")
  csv5 <- file.path(tdir, "other_names.csv"); mk(csv5, names = c("Alpha", "Beta", "Gamma"))
  expect_equal(d5_gate_run(d5_gate_args(c("--pred_dir", pred, "--results_csv", csv5, "--out_dir", file.path(tdir, "o5"))), say = function(...) NULL)$status, 2L)   # matched by name, never by value
  # the command line script refuses (status 2) with a missing input
  st <- suppressWarnings(system2("Rscript", c("analysis/dream5_scorer_gate.R", "--pred_dir", shQuote(file.path(tdir, "nonexistent")), "--results_csv", shQuote(csv)), stdout = FALSE, stderr = FALSE))
  expect_equal(st, 2L)
})

test_that("demo evidence: the aggregate scores of go_all.m are mean(-log10 p) per metric and overall = their mean (shipped OUTPUT/myteam.txt)", {
  skip_if_not(have_demo)
  ref <- readLines(file.path(SCR, "OUTPUT", "myteam.txt")); num <- function(i) as.numeric(strsplit(trimws(ref[i]), "\t")[[1]])
  a <- d5_gate_aggregates(num(5), num(8))                       # the printed p-values (4 significant digits)
  expect_lte(abs(a[["score_overall"]] - num(2)[1]), 1e-3); expect_lte(abs(a[["score_aupr"]] - num(2)[2]), 1e-3); expect_lte(abs(a[["score_auroc"]] - num(2)[3]), 1e-3)
  expect_equal(a[["score_overall"]], (a[["score_aupr"]] + a[["score_auroc"]]) / 2)
})

test_that("gate helpers: tolerance rule, Inf handling (no floor or cap), name matching, aggregates", {
  expect_equal(d5_gate_close(c(0.1, 0.1, 0.1, Inf, Inf, NA, 5), c(0.1005, 0.102, NA, Inf, 50, 1, 5.0011), 1e-3), c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
  a <- d5_gate_aggregates(c(0.1, 0.01, 0), c(0.5, 0.5, 0.5)); expect_true(is.infinite(a[["score_overall"]]) && is.infinite(a[["score_aupr"]]) && is.finite(a[["score_auroc"]]))
  b <- d5_gate_aggregates(c(1e-300, 1e-300, 1e-300), c(0.1, 0.1, 0.1)); expect_equal(b[["score_aupr"]], 300)          # tiny but nonzero p is not floored either
  pub <- data.frame(team = c("Regression1", "MI2"), score_overall = c(5, 1), score_aupr = c(6, 1), score_auroc = c(4, 1), aupr1 = .2, aupr3 = .1, aupr4 = .05, auroc1 = .7, auroc3 = .6, auroc4 = .55)
  d <- data.frame(team = "regression1", net = c(1, 3, 4), AUPR = c(.2, .1, .05) + 1e-4, AUROC = c(.7, .6, .55), P_AUPR = 10^-6, P_AUROC = 10^-4)
  cm <- d5_gate_compare(list(R = d), pub, 1e-3); expect_equal(nrow(cm), 3 * 2 + 3); expect_true(all(cm$team == "regression1"))              # matched by name, case-insensitively
  expect_equal(cm$ours[cm$quantity == "score_aupr"], 6); expect_equal(cm$ours[cm$quantity == "score_auroc"], 4); expect_true(all(cm$pass))
  d2 <- d[1:2, ]; expect_false(any(is.na(d5_gate_compare(list(matrix = d2), pub)$net) ))                                                         # a team lacking a network gets no aggregate rows
  expect_equal(nrow(d5_gate_compare(list(R = transform(d, team = "Unknown")), pub)), 0L)                                                        # unknown names are never compared
  tb <- data.frame(team = c("a", "b"), score_overall = 1, score_aupr = 1, score_auroc = 1, aupr1 = 1, aupr3 = 1, aupr4 = 1, auroc1 = 1, auroc3 = 1, auroc4 = "x")
  f <- tempfile(fileext = ".csv"); utils::write.csv(tb, f, row.names = FALSE); expect_error(d5_gate_read_results(f), "non-numeric or empty values in column auroc4")
  tb$auroc4 <- 1; tb$team <- "a"; utils::write.csv(tb, f, row.names = FALSE); expect_error(d5_gate_read_results(f), "duplicated team names")
})

test_that("overall official score combines the three per-network tables; Inf is not capped", {
  mk <- function(net, pa, po) { d <- tempfile(); dir.create(d); write.csv(data.frame(net = net, method = c("a", "b"), P_AUPR = pa, P_AUROC = po), file.path(d, "official_scores.csv"), row.names = FALSE); d }
  ds <- c(mk(1, c(.1, 0), c(.1, .01)), mk(3, c(.01, .5), c(.1, .1)), mk(4, c(.001, .5), c(.1, .1)))
  o <- d5_official_overall_from_dirs(ds)
  expect_equal(o$method, c("a", "b")); expect_equal(unname(unlist(o[1, c("overall", "aupr_score", "auroc_score")])), unname(d5_official_overall(c(.1, .01, .001), c(.1, .1, .1))))
  expect_true(is.infinite(o$overall[2]) && is.infinite(o$aupr_score[2]) && is.finite(o$auroc_score[2]))
  expect_error(d5_official_overall_from_dirs(ds[1:2]), "must cover networks")
})
