#!/usr/bin/env Rscript
# Reproduction gate for the GRNBoost2 rerun (27 Oct 2026): every row that does NOT depend on GRNBoost2 must be
# unchanged between the superseded metrics file and the rerun's.
#   Rscript analysis/rerun_gate.R old.csv new.csv key1,key2,... col1,col2,... [soft_methods=genie3]
# Rows are matched on the keys (which must include "method"); methods starting with "grnboost2" are skipped.
# Methods listed in soft_methods are reported but do not fail the gate (e.g. GENIE3 recomputed with parallel RNG).
# Exit 1 if any other method differs by more than 1e-6 or is missing.
a <- commandArgs(trailingOnly = TRUE)
old <- utils::read.csv(a[1], stringsAsFactors = FALSE); new <- utils::read.csv(a[2], stringsAsFactors = FALSE)
keys <- strsplit(a[3], ",")[[1]]; cols <- strsplit(a[4], ",")[[1]]
soft <- if (length(a) >= 5) strsplit(a[5], ",")[[1]] else character(0)
old <- old[!startsWith(old$method, "grnboost2"), ]; new <- new[!startsWith(new$method, "grnboost2"), ]
m <- merge(old, new, by = keys, suffixes = c(".old", ".new"), all = TRUE)
bad <- FALSE
miss <- m[is.na(m[[paste0(cols[1], ".old")]]) | is.na(m[[paste0(cols[1], ".new")]]), ]
if (nrow(miss)) { cat("rows present in only one file:", nrow(miss), "\n"); print(utils::head(miss[, keys], 10)); bad <- TRUE }
for (meth in unique(m$method)) {
  x <- m[m$method == meth, ]
  d <- vapply(cols, function(cn) max(abs(x[[paste0(cn, ".old")]] - x[[paste0(cn, ".new")]]), na.rm = TRUE), 0)
  ok <- all(d <= 1e-6); flag <- if (ok) "OK" else if (meth %in% soft) "DIFFERS (soft)" else "FAILS"
  cat(sprintf("  %-26s max |old - new| %s  %s\n", meth, paste(sprintf("%s=%.2e", cols, d), collapse = " "), flag))
  if (!ok && !(meth %in% soft)) bad <- TRUE
}
if (bad) { cat("GATE FAILED\n"); quit(status = 1) } else cat("GATE PASSED\n")
