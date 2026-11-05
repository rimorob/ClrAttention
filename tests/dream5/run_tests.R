#!/usr/bin/env Rscript
# ONE command runs every DREAM5 test (R testthat files and the GRNBoost2 --regulators Python test).
# Run from the repo root:  nice -n 19 Rscript tests/dream5/run_tests.R
# The Python test uses the clr-baselines virtualenv (override with D5_PYTHON=<python>); --skip_python 1 skips it.
suppressPackageStartupMessages(library(testthat))
Sys.setenv(D5_ROOT = normalizePath("."))
skip_py <- identical(commandArgs(trailingOnly = TRUE), c("--skip_python", "1"))
res <- as.data.frame(testthat::test_dir("tests/dream5", reporter = "summary", stop_on_failure = FALSE, env = new.env(parent = globalenv())))
cat(sprintf("\nDREAM5 R tests: %d expectations passed, %d failed, %d errors, %d skipped (%d test blocks)\n",
            sum(res$passed), sum(res$failed), sum(res$error), sum(res$skipped), nrow(res)))
py_fail <- 0L
if (!skip_py) {
  py <- Sys.getenv("D5_PYTHON", unset = path.expand("~/venvs/clr-baselines/bin/python"))
  out <- suppressWarnings(system2(py, "tests/dream5/test_grnboost2_regulators.py", stdout = TRUE, stderr = TRUE)); st <- attr(out, "status")
  npass <- sum(grepl("^PASS ", out)); nfail <- sum(grepl("^FAIL ", out)) + as.integer(!is.null(st) && st != 0L && !any(grepl("^FAIL ", out)))
  cat(out, sep = "\n"); cat(sprintf("DREAM5 Python tests (GRNBoost2 --regulators): %d checks passed, %d failed\n", npass, nfail)); py_fail <- nfail
}
quit(status = as.integer(sum(res$failed) + sum(res$error) + py_fail > 0))
