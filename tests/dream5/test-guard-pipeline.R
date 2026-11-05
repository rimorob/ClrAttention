setwd(Sys.getenv("D5_ROOT", unset = "."))

## ---- the truth guard (temporary git repos) -----------------------------------------------------------
git <- function(repo, ...) { x <- suppressWarnings(system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = TRUE)); invisible(x) }
make_repo <- function() {
  bare <- tempfile("bare"); system2("git", c("init", "--bare", "-q", "-b", "main", shQuote(bare)))
  repo <- tempfile("repo"); dir.create(repo); git(repo, "init", "-q", "-b", "main"); git(repo, "config", "user.email", "t@t"); git(repo, "config", "user.name", "t")
  git(repo, "remote", "add", "origin", shQuote(bare))
  dir.create(file.path(repo, "results/dream5/net1"), recursive = TRUE); dir.create(file.path(repo, "results/dream5/net4/refit"), recursive = TRUE)
  repo
}
guard_fixture <- function() {
  repo <- make_repo(); req <- d5_required(c(1L, 4L), c("clr", "attention"))
  files <- character(0)
  for (k in seq_len(nrow(req))) { f <- sprintf("results/dream5/net%d/%s__replicate.scores", req$net[k], req$method[k]); writeLines(paste("x", k), file.path(repo, f)); files <- c(files, f) }
  m <- data.frame(net = req$net, method = req$method, input = req$input, file = files, bytes = 1,
                  sha256 = vapply(files, function(f) d5_sha256(file.path(repo, f)), ""), stringsAsFactors = FALSE)
  write.table(m, file.path(repo, "results/dream5/SCORE_MANIFEST.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rf <- expand.grid(net = 4L, rep = 1:2, method = c("clr", "attention"), stringsAsFactors = FALSE); rf$sha256 <- paste0(rep("a", 64), collapse = "")
  write.table(rf, file.path(repo, "results/dream5/net4/refit/REFIT_MANIFEST.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(repo = repo, req = req, rf = rf)
}
commit_push <- function(repo, files) { git(repo, "add", "-f", files); git(repo, "commit", "-q", "-m", "x"); git(repo, "push", "-q", "-u", "origin", "main") }
M <- "results/dream5/SCORE_MANIFEST.tsv"; RM <- "results/dream5/net4/refit/REFIT_MANIFEST.tsv"

test_that("guard passes only for a committed, pushed, unchanged manifest whose files still match", {
  fx <- guard_fixture(); repo <- fx$repo
  expect_error(d5_truth_guard(M, fx$req, repo, verbose = FALSE), "not tracked")                  # untracked
  commit_push(repo, c(M, RM, list.files(file.path(repo, "results"), recursive = TRUE, full.names = FALSE) |> (\(x) file.path("results", x))()))
  expect_silent(d5_truth_guard(M, fx$req, repo, verbose = FALSE))
  rr <- data.frame(net = 4L, rep = 1:2, method = c("clr", "attention"))
  expect_silent(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, refit_required = expand.grid(net = 4L, rep = 1:2, method = c("clr", "attention"), stringsAsFactors = FALSE), verbose = FALSE))
  expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, refit_required = expand.grid(net = 4L, rep = 1:3, method = c("clr", "attention"), stringsAsFactors = FALSE), verbose = FALSE), "lacks")
  # the supplementary series-level manifest is part of the required set: absent -> refused; incomplete -> refused; then committed -> passes
  SM <- "results/dream5/net4/refit_series/REFIT_SERIES_MANIFEST.tsv"; sreq <- expand.grid(net = 4L, rep = 1:2, method = c("clr", "attention"), stringsAsFactors = FALSE)
  expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, series_manifest = SM, series_required = sreq, verbose = FALSE), "series-level refit manifest not found")
  dir.create(file.path(repo, dirname(SM)), recursive = TRUE); sf <- sreq; sf$sha256 <- strrep("b", 64)
  write.table(sf, file.path(repo, SM), sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, series_manifest = SM, series_required = sreq, verbose = FALSE), "not tracked")
  git(repo, "add", "-f", SM); git(repo, "commit", "-q", "-m", "series"); expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, series_manifest = SM, series_required = sreq, verbose = FALSE), "not pushed")
  git(repo, "push", "-q")
  expect_silent(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, series_manifest = SM, series_required = sreq, verbose = FALSE))
  expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, series_manifest = SM, series_required = expand.grid(net = 4L, rep = 1:3, method = c("clr", "attention"), stringsAsFactors = FALSE), verbose = FALSE), "series-level refit manifest lacks")
  writeLines(c(readLines(file.path(repo, SM)), ""), file.path(repo, SM))
  expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, series_manifest = SM, series_required = sreq, verbose = FALSE), "uncommitted"); git(repo, "checkout", "--", SM)
  # a required file missing from the manifest
  expect_error(d5_truth_guard(M, d5_required(c(1L, 4L), c("clr", "attention", "pidc")), repo, verbose = FALSE), "lacks required")
  # a score file changed after the manifest was committed
  f <- file.path(repo, "results/dream5/net1/clr__replicate.scores"); old <- readLines(f); writeLines("tampered", f)
  expect_error(d5_truth_guard(M, fx$req, repo, verbose = FALSE), "changed since the manifest"); writeLines(old, f)
  expect_silent(d5_truth_guard(M, fx$req, repo, verbose = FALSE))
  # a missing score file
  file.rename(f, paste0(f, ".bak")); expect_error(d5_truth_guard(M, fx$req, repo, verbose = FALSE), "missing"); file.rename(paste0(f, ".bak"), f)
  # manifest modified in the working tree
  writeLines(c(readLines(file.path(repo, M)), ""), file.path(repo, M)); expect_error(d5_truth_guard(M, fx$req, repo, verbose = FALSE), "uncommitted")
  git(repo, "checkout", "--", M); expect_silent(d5_truth_guard(M, fx$req, repo, verbose = FALSE))
  # committed but not pushed
  txt <- readLines(file.path(repo, M)); writeLines(txt, file.path(repo, M)); git(repo, "commit", "-q", "--allow-empty", "-m", "noop")
  m2 <- read.delim(file.path(repo, M), colClasses = c(sha256 = "character")); m2$bytes <- 2
  write.table(m2, file.path(repo, M), sep = "\t", quote = FALSE, row.names = FALSE); git(repo, "add", M); git(repo, "commit", "-q", "-m", "m2")
  expect_error(d5_truth_guard(M, fx$req, repo, verbose = FALSE), "not pushed")
  git(repo, "push", "-q"); expect_silent(d5_truth_guard(M, fx$req, repo, verbose = FALSE))
  # a tampered refit manifest in the working tree is refused too
  writeLines(c(readLines(file.path(repo, RM)), ""), file.path(repo, RM))
  expect_error(d5_truth_guard(M, fx$req, repo, refit_manifest = RM, verbose = FALSE), "uncommitted")
})

test_that("guard refuses without an upstream branch, and the design guard needs FIXED + committed", {
  fx <- guard_fixture(); repo <- fx$repo
  git(repo, "add", "-f", "results"); git(repo, "commit", "-q", "-m", "x")          # committed, never pushed, no upstream
  expect_error(d5_truth_guard(M, fx$req, repo, verbose = FALSE), "no upstream")
  dir.create(file.path(repo, "docs/plans"), recursive = TRUE); dir.create(file.path(repo, "docs/dream5_selections"), recursive = TRUE)
  writeLines(c("# t", "", "Status: DRAFT"), file.path(repo, "docs/plans/dream5_design.md")); writeLines("x", file.path(repo, "docs/dream5_selections/download_log.tsv"))
  git(repo, "add", "docs"); git(repo, "commit", "-q", "-m", "d"); git(repo, "push", "-q", "-u", "origin", "main")
  expect_error(d5_design_guard(repo = repo), "FIXED")
  writeLines(c("# t", "", "Status: FIXED 2026-11-04"), file.path(repo, "docs/plans/dream5_design.md"))
  expect_error(d5_design_guard(repo = repo), "uncommitted")
  git(repo, "commit", "-q", "-am", "fix"); expect_error(d5_design_guard(repo = repo), "not pushed"); git(repo, "push", "-q")
  expect_true(d5_design_guard(repo = repo))
})

test_that("the truth driver calls both guards before it opens the truth file; no code reads the sealed archive", {
  src <- paste(readLines("analysis/dream5_truth_run.R"), collapse = "\n")
  gpos <- regexpr("d5_truth_guard\\(", src); dpos <- regexpr("d5_design_guard\\(\\)", src); rpos <- regexpr("read\\.delim\\(opt\\$truth_edges", src)
  expect_true(gpos > 0 && dpos > 0 && rpos > 0 && gpos < rpos && dpos < rpos)
  expect_match(src, "series_manifest = series_manifest"); expect_match(src, "REFIT_SERIES_MANIFEST.tsv"); expect_match(src, "d5_refit_verify\\(v, b, 4L, series_manifest\\)")
  # label-free files never mention the scoring archive or gold-standard file names
  fs <- c("analysis/dream5_lib.R", "analysis/dream5_score_lib.R", "analysis/dream5_run.R", "analysis/dream5_truth_run.R", "tools/grnboost2.py")
  for (f in fs) expect_false(any(grepl("Evaluation scripts|gold_standard|synapse/scoring|\\.zip", readLines(f))), info = f)
})

## ---- the truth pipeline on a synthetic network with decoy truths ----------------------------------------------
test_that("pipeline end to end on decoys: guards-free core, jackknives, refit hash check, verdict.csv", {
  old <- D5$ntrees; D5$ntrees <- 20L; on.exit(D5$ntrees <- old, add = TRUE)
  nd <- d5_synth_net(300, 40, 24, 6, net = 4L); out <- tempfile()
  meths <- c("abs_pearson", "clr", "clrns", "attention", "genie3_tf", "grnboost2_tf")
  for (m in meths) suppressMessages(d5_run_method(nd, m, "replicate", 2L, out))
  for (b in 1:4) d5_refit_run(nd, b, 2L, out)
  for (b in 1:3) d5_refit_run(nd, b, 2L, out, unit = "series")
  sman <- d5_refit_manifest(4, out, "series"); sfile <- d5_refit_manifest_file(4, out, "series")
  rman <- d5_refit_manifest(4, out); rfile <- file.path(d5_refit_dir(4, out), "REFIT_MANIFEST.tsv")
  edges <- d5_decoy_truth(40, 300, nd$tf_idx, 91, median_out = 12, sd_log = 0.4)
  sw <- TESTSW(jackknife_cap = "reapply_cap")
  run <- function(...) d5_truth_pipeline(nd, edges, sw, methods = meths, out = out, nrep = 10L, refit_reps = 1:4, cap = 600L,
         refit_scorer = function(b) d5_refit_scores(nd, b, 2L), refit_verify = function(v, b) d5_refit_verify(v, b, 4L, rfile), verbose = FALSE, ...)
  r1 <- run(); r2 <- run()
  # supplementary series-level jackknife: scored, verified against its own committed hashes, never in the verdict
  runs <- function(...) run(series_reps = 1:3, series_scorer = function(b) d5_refit_scores(nd, b, 2L, unit = "series"),
                            series_verify = function(v, b) d5_refit_verify(v, b, 4L, sfile), ...)
  rs <- runs(); expect_identical(rs$rows, r1$rows); expect_identical(rs$refit_edge, r1$refit_edge); expect_identical(rs$jack_edge, r1$jack_edge)  # verdict inputs unchanged
  expect_null(r1$series_edge); expect_null(r1$series_line)
  expect_equal(dim(rs$series_edge), c(3, 4)); expect_equal(dim(rs$series_module), c(3, 4))
  sl <- rs$series_line; expect_equal(sl$endpoint, c("YE", "YM")); expect_equal(sl$n_rep, c(3, 3))
  expect_equal(sl$SD_series[1], sd(log2(rs$series_edge[, "attention"] / rs$series_edge[, "clr"])))
  expect_equal(sl$SD_series[2], sd(log2(rs$series_module[, "attention"] / rs$series_module[, "clr"])))
  expect_equal(sl$SE_B_group, c(rs$rows$YE$SE_B, rs$rows$YM$SE_B)); expect_true(all(grepl("series-level supplementary", sl$label)) && all(grepl("not used in verdict", sl$label)))
  dd <- tempfile(); d5_write_results(rs, dd); vcs <- read.csv(file.path(dd, "verdict.csv")); expect_false(any(grepl("series|SD_series", names(vcs))))
  expect_true(all(file.exists(file.path(dd, c("series_edge.csv", "series_module.csv", "supplementary_series_jackknife.csv")))))
  expect_equal(read.csv(file.path(dd, "supplementary_series_jackknife.csv"))$SD_series, sl$SD_series)
  d0 <- tempfile(); d5_write_results(r1, d0); expect_false(file.exists(file.path(d0, "supplementary_series_jackknife.csv")))
  expect_identical(read.csv(file.path(d0, "verdict.csv")), vcs)                                   # verdict.csv byte-for-byte the same content with or without the series replicates
  badS <- read.delim(sfile, colClasses = c(sha256 = "character")); badS$sha256[1] <- strrep("0", 64); sbf <- tempfile(); write.table(badS, sbf, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(run(series_reps = 1L, series_scorer = function(b) d5_refit_scores(nd, b, 2L, unit = "series"), series_verify = function(v, b) d5_refit_verify(v, b, 4L, sbf)), "differs from the committed hash")
  expect_error(run(series_reps = 4L, series_scorer = function(b) d5_refit_scores(nd, b, 2L, unit = "series"), series_verify = function(v, b) d5_refit_verify(v, b, 4L, sfile)), "no hash")
  expect_identical(r1$rows, r2$rows); expect_identical(r1$jack_edge, r2$jack_edge); expect_identical(r1$refit_edge, r2$refit_edge)       # fixed seeds
  expect_true(r1$gate$evaluable)
  expect_equal(dim(r1$jack_edge), c(10, 6)); expect_equal(dim(r1$refit_edge), c(4, 4)); expect_equal(dim(r1$jack_module), c(10, 4)); expect_equal(dim(r1$refit_module), c(4, 4))
  v <- r1$rows$YE
  expect_equal(v$d, log2(v$aupr_a / v$aupr_b)); expect_equal(v$SE_J, sd(log2(r1$jack_edge[, "attention"] / r1$jack_edge[, "clr"])))
  expect_equal(v$SE_B, sd(log2(r1$refit_edge[, "attention"] / r1$refit_edge[, "clr"]))); expect_equal(v$SE, sqrt(v$SE_J^2 + v$SE_B^2))
  expect_equal(v$ci_lo, v$d - 1.96 * v$SE); expect_equal(v$ci_hi, v$d + 1.96 * v$SE)
  vm <- r1$rows$YM; expect_equal(vm$SE, sqrt(vm$SE_J^2 + vm$SE_B^2))
  expect_true(all(c(v$verdict, vm$verdict) %in% c("Superior", "Inferior", "Equivalent", "Inconclusive")))
  expect_equal(v$verdict, d5_verdict(v$d, v$SE)$verdict)
  # written verdict.csv
  d <- tempfile(); d5_write_results(r1, d); vc <- read.csv(file.path(d, "verdict.csv"))
  expect_equal(vc$endpoint, c("YE", "YM")); expect_true(all(c("d", "SE_J", "SE_B", "SE", "ci_lo", "ci_hi", "verdict", "p_sup", "p_harm") %in% names(vc)))
  expect_true(all(file.exists(file.path(d, c("point_estimates.csv", "jack_edge.csv", "jack_module.csv", "refit_edge.csv", "refit_module.csv", "ym_gate.csv")))))
  # point estimates: |r| and the others finite; module AUPR above the base rate for CLR on planted structure is not asserted (decoy truth)
  expect_true(all(is.finite(r1$point$aupr)))
  ym <- r1$point[r1$point$endpoint == "YM", ]; expect_true(all(is.finite(ym$aupr_min10))); expect_true(all(is.finite(ym$aupr_sens))); expect_true(all(ym$aupr_sens >= ym$aupr - 1e-12))
  # tampered hash in the refit manifest: the truth stage refuses to score
  bad <- read.delim(rfile, colClasses = c(sha256 = "character")); bad$sha256[1] <- strrep("0", 64); bf <- tempfile(); write.table(bad, bf, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(d5_truth_pipeline(nd, edges, sw, methods = meths, out = out, nrep = 3L, refit_reps = 1L, cap = 600L, refit_scorer = function(b) d5_refit_scores(nd, b, 2L),
               refit_verify = function(v, b) d5_refit_verify(v, b, 4L, bf), verbose = FALSE), "differs from the committed hash")
  # net other than 4: effects reported with SE_J only, no verdict
  nd3 <- nd; nd3$net <- 3L
  for (m in c("clr", "attention")) file.copy(d5_cache_path(4, m, "replicate", out), { dir.create(d5_dir(3, out), showWarnings = FALSE); d5_cache_path(3, m, "replicate", out) })
  r3 <- d5_truth_pipeline(nd3, edges, sw, methods = c("clr", "attention"), out = out, nrep = 5L, cap = 600L, verbose = FALSE)
  expect_match(r3$rows$YE$verdict, "descriptive"); expect_true(is.na(r3$rows$YE$SE_B)); expect_true(is.finite(r3$rows$YE$SE_J))
  # Y-M not evaluable: the row says so and nothing is scored for it
  e_small <- d5_decoy_truth(40, 300, nd$tf_idx, 5, median_out = 2, sd_log = 0.1)
  r4 <- d5_truth_pipeline(nd, e_small, sw, methods = meths, out = out, nrep = 3L, cap = 600L, verbose = FALSE)
  expect_equal(r4$rows$YM$verdict, "not evaluable"); expect_null(r4$jack_module); expect_false("YM" %in% r4$point$endpoint)
})

test_that("the truth driver refuses to run in the current repository state when a guard fails (nonexistent truth path, never read)", {
  out <- suppressWarnings(system2("Rscript", c("analysis/dream5_truth_run.R", "--net", "4", "--truth_edges", "/nonexistent/never_read.tsv", "--out_dir", tempfile()), stdout = TRUE, stderr = TRUE))
  expect_false(is.null(attr(out, "status")))
  expect_true(any(grepl("truth guard|cannot open", out)))
  expect_false(any(grepl("truth: [0-9]+ edges kept", out)))
})
