# Shared helpers for the complex-truth competitor cells of the scoreboard (docs/plans/scoreboard_competitors_plan.md).
# Truth construction is NOT re-implemented: the function definitions of analysis/clr_logpt_test.R (truth builders for M3D,
# EcoMAC and RPE1 CORUM, data loaders) are evaluated verbatim from that file (top-level function / constant definitions only;
# none of its modes, gates or outputs run). Source from the repo root.
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
SYN <- FALSE; script_dir <- "analysis"; opt <- list(expr = ""); SELDIR <- "docs/clr_logpt_selections"
sb_import_truth_defs <- function(file = "analysis/clr_logpt_test.R",
    want = c("merge_sets", "h5_cat", "half_split", "load_expr", "mk_truth", "as_reg", "ecoli_regtab", "ecoli_members",
             "complexes_with_pos", "build_truths_ecoli", "build_truths_human", "build_truths", "none_op", "RP_RE")) {
  ex <- parse(file, keep.source = FALSE); got <- character(0)
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% want) {
    eval(e, globalenv()); got <- c(got, as.character(e[[2]])) }
  stopifnot(setequal(got, want)); invisible(got)
}
sb_say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n") }
sb_sym_max <- function(W) { W[is.na(W)] <- 0; W <- pmax(W, t(W)); diag(W) <- 0; W }   # baselines.R sym_max (+ zero diagonal)
