#!/usr/bin/env Rscript
# Method comparison ("SOTA") table: pooled co-membership AUPR over base rate for
# each dataset x truth, raw and hub-stratified. A method is bold when it ties
# the best method in its row: the paired difference best - method, over the
# shared delete-half gene jackknife replicates (jack_replicates.csv from
# analysis/baselines.R; SE = sd of the replicate differences, delete-d factor 1),
# has a 95% interval that includes 0.
#
# Usage (repo root):  Rscript analysis/sota_table.R [--out docs/whitepapers/sota_table.md]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(out = "docs/whitepapers/sota_table.md")
if (length(args) >= 2 && args[1] == "--out") opt$out <- args[2]
runs <- c(M3D = "results/baselines_m3d", EcoMAC = "results/baselines_ecomac",
          K562 = "results/baselines_k562ess_truthgenes_genie3")
methods <- c(attention = "Attention", clr = "CLR", genie3 = "GENIE3", pidc = "PIDC",
             grnboost2 = "GRNBoost2", abs_pearson = "|r|", hub = "Hub only")
hub_methods <- c("attention", "clr", "genie3", "pidc")

lines <- c("# Method comparison", "",
  sprintf("*Generated %s by `analysis/sota_table.R`.* Pooled co-membership AUPR divided by the base rate (1.00 = random), CLR attention at t = 9. **Bold**: best in the row or tied with it (95%% interval of the paired jackknife difference includes 0). \"-\": not run.", format(Sys.time(), "%Y-%m-%d %H:%M")), "")

table_for <- function(hubstrat) {
  ms <- if (hubstrat) hub_methods else names(methods)
  out <- c(paste0("| Dataset | Truth | ", paste(methods[ms], collapse = " | "), " |"),
           paste0("|---|---|", paste(rep("---", length(ms)), collapse = "|"), "|"))
  for (d in names(runs)) {
    mf <- file.path(runs[[d]], "metrics.csv"); jf <- file.path(runs[[d]], "jack_replicates.csv")
    if (!file.exists(mf)) next
    x <- utils::read.csv(mf, stringsAsFactors = FALSE)
    j <- if (file.exists(jf)) utils::read.csv(jf, stringsAsFactors = FALSE) else NULL
    for (tr in unique(x$truth)) {
      y <- x[x$truth == tr, ]; key <- if (hubstrat) paste0(ms, "_hubstrat") else ms
      a <- stats::setNames(y$aupr[match(key, y$method)], ms); b <- y$base_rate[1]
      best <- names(which.max(a))
      tie <- vapply(ms, function(m) {
        if (is.na(a[[m]])) return(NA)
        if (m == best) return(TRUE)
        if (is.null(j)) return(NA)
        jb <- j[j$truth == tr & j$method == key[ms == best], ]; jm <- j[j$truth == tr & j$method == key[ms == m], ]
        if (!nrow(jb) || !nrow(jm)) return(NA)
        dd <- jb$aupr[order(jb$rep)] - jm$aupr[order(jm$rep)]
        (a[[best]] - a[[m]]) < 1.96 * stats::sd(dd) }, NA)
      cell <- vapply(ms, function(m) {
        if (is.na(a[[m]])) return("-")
        v <- sprintf("%.2f", a[[m]] / b)
        if (isTRUE(tie[[m]])) paste0("**", v, "**") else if (is.na(tie[[m]])) paste0(v, "?") else v }, "")
      out <- c(out, sprintf("| %s | %s | %s |", d, gsub("_", " ", tr), paste(cell, collapse = " | ")))
    }
  }
  out
}
lines <- c(lines, "## Raw scores", "", table_for(FALSE), "",
           "## Hub-controlled (within-stratum percentile in 10 strata of pi_i + pi_j)", "", table_for(TRUE), "",
           "\"?\": tie test unavailable (no jackknife replicates saved for this run yet).")
writeLines(lines, opt$out)
cat(lines, sep = "\n")
