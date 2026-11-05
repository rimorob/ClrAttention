#!/usr/bin/env Rscript
# Run from the repo root:  nice -n 19 Rscript tests/dream5/run_tests.R
suppressPackageStartupMessages(library(testthat))
Sys.setenv(D5_ROOT = normalizePath("."))
res <- as.data.frame(testthat::test_dir("tests/dream5", reporter = "summary", stop_on_failure = FALSE, env = new.env(parent = globalenv())))
cat(sprintf("\nDREAM5 R tests: %d expectations passed, %d failed, %d errors, %d skipped (%d test blocks)\n",
            sum(res$passed), sum(res$failed), sum(res$error), sum(res$skipped), nrow(res)))
quit(status = as.integer(sum(res$failed) + sum(res$error) > 0))
