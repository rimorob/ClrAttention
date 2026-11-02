#!/usr/bin/env Rscript
# EXPLORATORY (not preregistered) checks on the scored BEELINE table; reads only results/scfm_score/{metrics_long,primary_cells}.csv.
#  (a) the mHSC-E_500 organism-head flag (QC: 0.67 of its cells classified human): primary contrast with mHSC-E taken from TFs+1000 only
#  (b) leave-one-dataset-out: the same exact interval on the 6 remaining d_i (C by the same rule, n = 6)
#  (c) the per-truth d_i (primary contrast restricted to one truth), n = 7
# Usage: Rscript analysis/scfm_score_exploratory.R --dir results/scfm_score
a <- commandArgs(TRUE); dir <- if (length(a) >= 2) a[2] else "results/scfm_score"
DELTA <- 0.15
src <- readLines("analysis/scfm_score_stats.R")           # exact-interval helpers of the stats script (walsh, c_crit, hl, p_ni)
for (e in parse(text = src)) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] %in% c("walsh", "c_crit", "hl", "p_ni")) eval(e)
DELTA <- 0.15
p <- read.csv(file.path(dir, "primary_cells.csv"))
DS <- c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L")
dv <- function(x) { v <- tapply(x$lr, x$dataset, mean); v[DS[DS %in% names(v)]] }
out <- character()
say <- function(...) { m <- paste0(...); cat(m, "\n"); out <<- c(out, m) }
d <- dv(p); h <- hl(d)
say("EXPLORATORY CHECKS (not preregistered)")
say(sprintf("primary (reference): HL %+.3f, [%+.3f, %+.3f]", h$hl, h$L, h$U))
q <- p[!(p$dataset == "mHSC-E" & p$N == 500), ]; d2 <- dv(q); h2 <- hl(d2)
say(sprintf("(a) mHSC-E from TFs+1000 only (mHSC-E_500 organism-head flag): d_mHSC-E %+.3f (was %+.3f); HL %+.3f, [%+.3f, %+.3f]", d2["mHSC-E"], d["mHSC-E"], h2$hl, h2$L, h2$U))
say("(b) leave-one-dataset-out (n = 6; coverage ", sprintf("%.3f", hl(d[-1])$coverage), "):")
for (k in seq_along(d)) { hh <- hl(d[-k]); say(sprintf("    without %-8s HL %+.3f, [%+.3f, %+.3f], L > -delta: %s", names(d)[k], hh$hl, hh$L, hh$U, hh$L > -DELTA)) }
say("(c) per-truth d_i (n = 7 except LOFGOF, mESC only):")
for (tn in c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")) { v <- dv(p[p$truth == tn, ]); hh <- hl(v)
  say(sprintf("    %-17s n=%d  %s  HL %+.3f%s", tn, length(v), paste(sprintf("%s %+.2f", names(v), v), collapse = ", "), hh$hl,
              if (is.finite(hh$L)) sprintf(" [%+.3f, %+.3f]", hh$L, hh$U) else "")) }
writeLines(out, file.path(dir, "exploratory.txt"))
