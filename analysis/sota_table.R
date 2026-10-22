#!/usr/bin/env Rscript
# Method comparison ("SOTA") table: pooled co-membership AUPR over base rate for
# each dataset x truth, raw and hub-stratified. A method is bold when it ties
# the best method in its row: the paired difference best - method, over the
# shared delete-half gene jackknife replicates (jack_replicates.csv from
# analysis/baselines.R, merged on the replicate index; SE = sd of the replicate
# differences, delete-d factor 1), has a Bonferroni-corrected 95% interval
# (over the comparisons in the row) that includes 0. A dagger marks cells tied
# only under the correction (uncorrected interval excludes 0).
#
# Usage (repo root):  Rscript analysis/sota_table.R [--out docs/whitepapers/sota_table.md]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(out = "docs/whitepapers/sota_table.md")
if (length(args) >= 2 && args[1] == "--out") opt$out <- args[2]
runs <- c(M3D = "results/baselines_m3d", EcoMAC = "results/baselines_ecomac",
          K562 = "results/baselines_k562ess_truthgenes_genie3")
methods <- c(attention = "Attention", clr = "CLR", genie3 = "GENIE3", pidc = "PIDC",
             grnboost2 = "GRNBoost2", abs_pearson = "\\|r\\|", hub = "Hub only")
hub_methods <- c("attention", "clr", "genie3", "pidc", "grnboost2", "abs_pearson")

lines <- c("# Method comparison", "",
  sprintf("*Generated %s by `analysis/sota_table.R`.* Pooled co-membership AUPR divided by the base rate (1.00 = random); CLR attention at t = 9. **Bold**: best in the row or tied with it (Bonferroni-corrected 95%% interval of the paired jackknife difference includes 0); **x**† tied only after the correction. \"-\": not run.", format(Sys.time(), "%Y-%m-%d %H:%M")), "")

table_for <- function(hubstrat) {
  ms <- if (hubstrat) hub_methods else names(methods)
  out <- c(paste0("| Dataset | Truth | Base rate | ", paste(methods[ms], collapse = " | "), " |"),
           paste0("|---|---|---|", paste(rep("---", length(ms)), collapse = "|"), "|"))
  for (d in names(runs)) {
    mf <- file.path(runs[[d]], "metrics.csv"); jf <- file.path(runs[[d]], "jack_replicates.csv")
    if (!file.exists(mf)) next
    x <- utils::read.csv(mf, stringsAsFactors = FALSE)
    j <- if (file.exists(jf)) utils::read.csv(jf, stringsAsFactors = FALSE) else NULL
    for (tr in unique(x$truth)) {
      y <- x[x$truth == tr, ]; key <- if (hubstrat) paste0(ms, "_hubstrat") else ms
      a <- stats::setNames(y$aupr[match(key, y$method)], ms); b <- y$base_rate[1]
      best <- names(which.max(a)); ncomp <- sum(!is.na(a)) - 1L
      zc <- stats::qnorm(1 - 0.025 / max(1L, ncomp))           # Bonferroni over the row's comparisons
      z <- vapply(ms, function(m) {                              # (best - m) / SE of the paired difference
        if (is.na(a[[m]]) || m == best || is.null(j)) return(NA_real_)
        jb <- j[j$truth == tr & j$method == key[ms == best], c("rep", "aupr")]
        jm <- j[j$truth == tr & j$method == key[ms == m], c("rep", "aupr")]
        mm <- merge(jb, jm, by = "rep")
        if (nrow(mm) != nrow(jb) || nrow(mm) != nrow(jm)) stop("jackknife replicates do not match for ", d, " ", tr, " ", m)
        (a[[best]] - a[[m]]) / stats::sd(mm$aupr.x - mm$aupr.y) }, 0)
      cell <- vapply(ms, function(m) {
        if (is.na(a[[m]])) return("-")
        v <- sprintf("%.2f", a[[m]] / b)
        if (m == best || (!is.na(z[[m]]) && z[[m]] < 1.96)) paste0("**", v, "**")
        else if (!is.na(z[[m]]) && z[[m]] < zc) paste0("**", v, "**†") else v }, "")
      out <- c(out, sprintf("| %s | %s | %.3f | %s |", d, gsub("_", " ", tr), b, paste(cell, collapse = " | ")))
    }
  }
  out
}
lines <- c(lines, "## Raw scores", "", table_for(FALSE), "",
           "## Hub-controlled (within-stratum percentile in 10 strata of pi_i + pi_j)", "", table_for(TRUE), "",
           "**Notes.**",
           "- Intervals resample genes (delete-half jackknife, 100 replicates) with the fitted networks held fixed; they do not resample experiments or cells, or vary the random seeds of GENIE3/GRNBoost2.",
           "- Depth t = 9 was chosen on M3D with a label-free held-out criterion (lab splits; no truths used), then transferred unchanged to EcoMAC and K562. M3D is therefore not an untouched test set for the depth.",
           "- Hub control: each method's scores are converted to percentiles within 10 strata of pi_i + pi_j, where pi is the stationary distribution of the attention walk (an attention-derived stratifier; a method-neutral one is planned). \"Hub only\" ranks pairs by pi_i + pi_j alone (a degree baseline).",
           "- K562: the Perturb-seq essential compendium, network built from knockdown half A on the 5,671 genes in at least one literature truth; the knockdown-response truth comes from the other half (B). All K562 TF-truth rows are within 10% of random.",
           "- EcoMAC may share arrays with M3D, so it is not yet an independent replication.",
           "- GRNBoost2 returns a sparse network (most pairs tied at zero), which lowers its pooled AUPR.")
writeLines(lines, opt$out)
cat(lines, sep = "\n")
