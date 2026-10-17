#!/usr/bin/env Rscript
# Macro-averaged per-regulon AP: every regulon counts once, next to the pooled
# (micro) co-membership AUPR that large regulons dominate. Reads the
# per_regulon.csv tables written by analysis/baselines.R. Definition fixed
# with the statistics reviewer before scoring (docs/paper_plan.md, theme 5):
#   - candidates as regulon_coherence(); AP with ties broken in expectation;
#     regulons with >= 10 positive (non-same-operon) pairs; the exclusions
#     depend on the truth only, so every method is scored on the same regulons
#   - primary: per-regulon paired difference dAP = AP(method) - AP(attention);
#     mean (regulon bootstrap, 10,000, BCa) and median (Wilcoxon signed rank),
#     plus the share of regulons where attention wins
#   - levels reported as (AP - b) / (1 - b), b = the regulon's base rate
#   - strata: size 5-9 / 10-24 / 25-99 / >= 100; class (TF vs sigma)
#   - macro-to-micro curve: regulons weighted by n_pos^gamma, gamma 0..1
#   - sensitivity: regulons with > 50% same-operon member pairs removed;
#     leave-one-regulon-out jackknife SE of the mean
#   - Holm families declared in advance: E. coli (M3D, EcoMAC x SC, all x
#     comparators x mean, median) and K562 (truths x comparators x mean, median)
#
# Usage (repo root):  Rscript analysis/macro_aupr.R [--dirs results/baselines_m3d,results/baselines_ecomac,...]
#                       [--out results/macro_aupr] [--B 10000]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(dirs = "results/baselines_m3d,results/baselines_ecomac,results/baselines_k562ess_genie3",
            out = "results/macro_aupr", B = 10000L, seed = 20261017L)
i <- 1L
while (i <= length(args)) { k <- sub("^--", "", args[i]); opt[[k]] <- if (k %in% c("B", "seed")) as.integer(args[i + 1L]) else args[i + 1L]; i <- i + 2L }
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "macro_aupr.log")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages(library(boot))
set.seed(opt$seed)

files <- file.path(strsplit(opt$dirs, ",")[[1]], "per_regulon.csv")
files <- files[file.exists(files)]
if (!length(files)) stop("no per_regulon.csv found")
pr <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
pr <- pr[!is.na(pr$ap), ]
pr$ap_norm <- (pr$ap - pr$base_rate) / (1 - pr$base_rate)
pr$size_bin <- cut(pr$size, c(0, 9, 24, 99, Inf), labels = c("5-9", "10-24", "25-99", ">=100"))
pr$family <- ifelse(pr$dataset %in% c("m3d", "ecomac"), "E. coli", "K562")
say("per-regulon rows: ", nrow(pr), " from ", paste(files, collapse = ", "))

boot_mean <- function(d) {
  if (length(d) < 5) return(c(lo = NA, hi = NA, p = NA))
  b <- boot::boot(d, function(x, i) mean(x[i]), R = opt$B)
  ci <- tryCatch(boot::boot.ci(b, type = "bca")$bca[4:5], error = function(e) stats::quantile(b$t, c(0.025, 0.975)))
  c(lo = ci[1], hi = ci[2], p = min(1, max(1 / opt$B, 2 * min(mean(b$t <= 0), mean(b$t >= 0)))))   # floored at 1/B
}
loro_se <- function(d) { n <- length(d); if (n < 3) return(NA); m <- (sum(d) - d) / (n - 1); sqrt((n - 1) / n * sum((m - mean(m))^2)) }

rows <- list(); strata <- list(); gam <- list(); levels <- list()
for (key in unique(paste(pr$dataset, pr$truth, sep = "|"))) {
  x <- pr[paste(pr$dataset, pr$truth, sep = "|") == key, ]
  ds <- x$dataset[1]; tr <- x$truth[1]
  for (m in unique(x$method)) { y <- x[x$method == m, ]
    levels[[length(levels) + 1L]] <- data.frame(dataset = ds, truth = tr, method = m, n_regulons = nrow(y),
      mean_ap_norm = mean(y$ap_norm), median_ap_norm = stats::median(y$ap_norm),
      median_fold = stats::median(y$ap / y$base_rate), median_nonzero = stats::median(y$nonzero)) }
  att <- x[x$method == "attention", ]
  for (m in setdiff(unique(x$method), "attention")) {
    y <- x[x$method == m, ]; j <- match(att$regulon, y$regulon); ok <- !is.na(j)
    a <- att[ok, ]; d <- y$ap[j[ok]] - a$ap                     # comparator - attention
    bm <- boot_mean(d)
    wil <- if (sum(d != 0) >= 5) stats::wilcox.test(d, exact = FALSE)$p.value else NA
    noop <- a$frac_same_operon <= 0.5
    rows[[length(rows) + 1L]] <- data.frame(family = x$family[1], dataset = ds, truth = tr, comparator = m,
      n_regulons = length(d), mean_d = mean(d), mean_lo = bm[["lo"]], mean_hi = bm[["hi"]], p_mean = bm[["p"]],
      loro_se = loro_se(d), median_d = stats::median(d), p_median = wil, attention_wins = mean(d < 0),
      mean_d_no_operon_dominated = mean(d[noop]), median_d_no_operon_dominated = stats::median(d[noop]))
    for (bn in levels(pr$size_bin)) { s <- a$size_bin == bn
      if (any(s)) strata[[length(strata) + 1L]] <- data.frame(dataset = ds, truth = tr, comparator = m, stratum = paste("size", bn),
        n = sum(s), mean_d = mean(d[s]), median_d = stats::median(d[s]), attention_wins = mean(d[s] < 0)) }
    for (cl in unique(a$class)) { s <- a$class %in% cl
      strata[[length(strata) + 1L]] <- data.frame(dataset = ds, truth = tr, comparator = m, stratum = paste("class", cl),
        n = sum(s), mean_d = mean(d[s]), median_d = stats::median(d[s]), attention_wins = mean(d[s] < 0)) }
    for (g in c(0, 0.25, 0.5, 0.75, 1)) { w <- a$n_pos^g
      gam[[length(gam) + 1L]] <- data.frame(dataset = ds, truth = tr, comparator = m, gamma = g, weighted_mean_d = sum(w * d) / sum(w)) }
  }
}
res <- do.call(rbind, rows)
long <- rbind(data.frame(i = seq_len(nrow(res)), stat = "mean", p = res$p_mean, family = res$family),
              data.frame(i = seq_len(nrow(res)), stat = "median", p = res$p_median, family = res$family))
long$p_holm <- NA
for (f in unique(long$family)) { s <- long$family == f; long$p_holm[s] <- stats::p.adjust(long$p[s], "holm") }
res$p_mean_holm <- long$p_holm[long$stat == "mean"]; res$p_median_holm <- long$p_holm[long$stat == "median"]
utils::write.csv(res, file.path(opt$out, "macro_vs_attention.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, levels), file.path(opt$out, "macro_levels.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, strata), file.path(opt$out, "macro_strata.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, gam), file.path(opt$out, "macro_gamma.csv"), row.names = FALSE)

say("comparator - attention, per-regulon AP (negative = attention better); Holm within family")
for (k in seq_len(nrow(res))) with(res[k, ], say(sprintf(
  "%-7s %-14s %-10s n=%3d  mean %+.4f [%+.4f, %+.4f] p_holm %.3g | median %+.4f p_holm %.3g | attention wins %.0f%%",
  dataset, truth, comparator, n_regulons, mean_d, mean_lo, mean_hi, p_mean_holm, median_d, p_median_holm, 100 * attention_wins)))
say("done")
