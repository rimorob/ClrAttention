#!/usr/bin/env Rscript
# Macro-averaged per-regulon AP: every regulon counts once, next to the pooled
# (micro) co-membership AUPR that large regulons dominate. Reads the
# per_regulon.csv and regulon_members.csv tables written by analysis/baselines.R.
# Definition fixed with the statistics reviewer before scoring (docs/paper_plan.md, theme 5):
#   - candidates as regulon_coherence(); AP with ties broken in expectation;
#     regulons with >= 10 positive (non-same-operon) pairs (so CORUM complexes
#     of 3-4 genes drop out); exclusions depend on the truth only
#   - comparisons only within one run (one results directory), so every
#     method is paired with the attention matrix of the same run; duplicated
#     (run, method, truth, regulon) rows stop the script
#   - primary: per-regulon paired difference dAP = AP(comparator) - AP(attention);
#     mean (regulon bootstrap, BCa interval, p by inverting the BCa interval)
#     and median (Wilcoxon signed rank); share of regulons where attention
#     wins, and share tied
#   - overlap-aware sensitivity: bootstrap over connected groups of regulons
#     (linked when member Jaccard > 0.2), percentile interval
#   - levels as (AP - b) / (1 - b), b = the regulon's base rate
#   - strata: size 5-9 / 10-24 / 25-99 / >= 100; class (TF vs sigma)
#   - pair-weighted means of per-regulon AP, weights n_pos^gamma, gamma 0..1
#     (gamma = 1 is NOT the pooled AUPR, which ranks all pairs together)
#   - operon-dominated regulons (> 50% same-operon member pairs) removed, with
#     interval and test
#   - Holm families fixed in advance, with fixed sizes (missing entries count):
#       E. coli: 5 comparators (abs_pearson, clr, genie3, grnboost2, pidc) x
#                runs baselines_m3d, baselines_ecomac x evidence SC, all x
#                mean, median = 40 tests
#       K562 primary: run baselines_k562ess_truthgenes_genie3 x 5 comparators x
#                5 truths x mean, median = 50 tests
#       K562 secondary (all genes; shares most genes with the primary run, not
#                independent): baselines_k562ess_genie3, 40 tests (no PIDC)
#
# Usage (repo root):  Rscript analysis/macro_aupr.R [--dirs dir1,dir2,...] [--out results/macro_aupr] [--B 10000]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(dirs = paste("results/baselines_m3d", "results/baselines_ecomac", "results/baselines_k562ess_genie3",
                         "results/baselines_k562ess_truthgenes_genie3", sep = ","),
            out = "results/macro_aupr", B = 10000L, seed = 20261017L, jaccard = 0.2)
i <- 1L
while (i <= length(args)) { k <- sub("^--", "", args[i])
  opt[[k]] <- if (k %in% c("B", "seed")) as.integer(args[i + 1L]) else if (k == "jaccard") as.numeric(args[i + 1L]) else args[i + 1L]
  i <- i + 2L }
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "macro_aupr.log")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages(library(boot))
set.seed(opt$seed)

families <- list(
  "E. coli"        = list(runs = c("baselines_m3d", "baselines_ecomac"), n = 40L),
  "K562 primary"   = list(runs = "baselines_k562ess_truthgenes_genie3", n = 50L),
  "K562 secondary" = list(runs = "baselines_k562ess_genie3", n = 40L))
family_of <- function(run) { f <- names(families)[vapply(families, function(x) run %in% x$runs, NA)]
  if (length(f)) f else "other (no Holm family)" }

dirs <- strsplit(opt$dirs, ",")[[1]]
f1 <- file.path(dirs, "per_regulon.csv"); f2 <- file.path(dirs, "regulon_members.csv")
miss <- c(f1, f2)[!file.exists(c(f1, f2))]
if (length(miss)) stop("missing input: ", paste(miss, collapse = ", "))
pr <- do.call(rbind, lapply(f1, utils::read.csv, stringsAsFactors = FALSE))
mem <- do.call(rbind, lapply(f2, utils::read.csv, stringsAsFactors = FALSE))
key <- paste(pr$run, pr$method, pr$truth, pr$regulon)
if (anyDuplicated(key)) stop("duplicated (run, method, truth, regulon) rows: ", paste(utils::head(unique(key[duplicated(key)])), collapse = "; "))
pr <- pr[!is.na(pr$ap), ]
pr$ap_norm <- (pr$ap - pr$base_rate) / (1 - pr$base_rate)
pr$size_bin <- cut(pr$size, c(0, 9, 24, 99, Inf), labels = c("5-9", "10-24", "25-99", ">=100"))
say("per-regulon rows: ", nrow(pr), " from ", paste(dirs, collapse = ", "))

boot_mean <- function(d) {
  if (length(d) < 5 || all(d == d[1])) return(c(lo = NA, hi = NA, p = NA))
  b <- boot::boot(d, function(x, i) mean(x[i]), R = opt$B)
  ci <- function(a) boot::boot.ci(b, conf = 1 - a, type = "bca")$bca[4:5]
  lohi <- ci(0.05)
  excl <- function(a) { r <- ci(a); r[1] > 0 || r[2] < 0 }       # BCa interval at level 1 - a excludes 0
  if (!excl(1 - 1e-9)) p <- 1 else if (excl(1 / opt$B)) p <- 1 / opt$B else {
    lo <- 1 / opt$B; hi <- 1 - 1e-9                                 # smallest a at which it excludes 0
    for (it in 1:30) { mid <- sqrt(lo * hi); if (excl(mid)) hi <- mid else lo <- mid }
    p <- hi }
  c(lo = lohi[1], hi = lohi[2], p = p)
}
groups_of <- function(run, truth, regs) {                        # connected components, member Jaccard > threshold
  m <- mem[mem$run == run & mem$truth == truth, ]; m <- m[match(regs, m$regulon), ]
  sets <- strsplit(m$members, ";", fixed = TRUE); n <- length(sets); g <- seq_len(n)
  find <- function(x) { while (g[x] != x) x <- g[x]; x }
  for (a in seq_len(n - 1L)) for (b in (a + 1L):n) {
    inter <- length(intersect(sets[[a]], sets[[b]]))
    if (inter && inter / length(union(sets[[a]], sets[[b]])) > opt$jaccard) { ra <- find(a); rb <- find(b); if (ra != rb) g[rb] <- ra }
  }
  vapply(seq_len(n), find, 0L)
}
cluster_ci <- function(d, grp) {
  if (length(unique(grp)) < 5) return(c(NA, NA))
  sp <- split(d, grp); s <- vapply(sp, sum, 0); n <- lengths(sp)
  bt <- replicate(opt$B, { k <- sample(length(sp), replace = TRUE); sum(s[k]) / sum(n[k]) })
  stats::quantile(bt, c(0.025, 0.975), names = FALSE)
}
wil <- function(d) if (sum(d != 0) >= 5) stats::wilcox.test(d, exact = FALSE)$p.value else NA_real_

rows <- list(); strata <- list(); gam <- list(); levels <- list()
for (k in unique(paste(pr$run, pr$truth, sep = "|"))) {
  x <- pr[paste(pr$run, pr$truth, sep = "|") == k, ]
  rn <- x$run[1]; ds <- x$dataset[1]; tr <- x$truth[1]
  for (m in unique(x$method)) { y <- x[x$method == m, ]
    levels[[length(levels) + 1L]] <- data.frame(run = rn, dataset = ds, truth = tr, method = m, n_regulons = nrow(y),
      mean_ap_norm = mean(y$ap_norm), median_ap_norm = stats::median(y$ap_norm),
      median_fold = stats::median(y$ap / y$base_rate), median_nonzero = stats::median(y$nonzero)) }
  att <- x[x$method == "attention", ]
  if (!nrow(att)) next
  grp_all <- groups_of(rn, tr, att$regulon)
  for (m in setdiff(unique(x$method), "attention")) {
    y <- x[x$method == m, ]; j <- match(att$regulon, y$regulon); ok <- !is.na(j)
    a <- att[ok, ]; d <- y$ap[j[ok]] - a$ap                     # comparator - attention
    bm <- boot_mean(d); cc <- cluster_ci(d, grp_all[ok])
    noop <- a$frac_same_operon <= 0.5; bn <- boot_mean(d[noop])
    rows[[length(rows) + 1L]] <- data.frame(family = family_of(rn), run = rn, dataset = ds, truth = tr, comparator = m,
      n_regulons = length(d), n_overlap_groups = length(unique(grp_all[ok])),
      mean_d = mean(d), mean_lo = bm[["lo"]], mean_hi = bm[["hi"]], p_mean = bm[["p"]],
      mean_lo_overlap = cc[1], mean_hi_overlap = cc[2],
      median_d = stats::median(d), p_median = wil(d), attention_wins = mean(d < 0), tied = mean(d == 0),
      n_no_operon_dominated = sum(noop), mean_d_noop = mean(d[noop]), mean_lo_noop = bn[["lo"]], mean_hi_noop = bn[["hi"]],
      median_d_noop = stats::median(d[noop]), p_median_noop = wil(d[noop]))
    for (b in levels(pr$size_bin)) { s <- a$size_bin == b
      if (any(s)) strata[[length(strata) + 1L]] <- data.frame(run = rn, truth = tr, comparator = m, stratum = paste("size", b),
        n = sum(s), mean_d = mean(d[s]), median_d = stats::median(d[s]), attention_wins = mean(d[s] < 0), tied = mean(d[s] == 0)) }
    for (cl in unique(a$class)) { s <- a$class %in% cl
      strata[[length(strata) + 1L]] <- data.frame(run = rn, truth = tr, comparator = m, stratum = paste("class", cl),
        n = sum(s), mean_d = mean(d[s]), median_d = stats::median(d[s]), attention_wins = mean(d[s] < 0), tied = mean(d[s] == 0)) }
    for (g in c(0, 0.25, 0.5, 0.75, 1)) { w <- a$n_pos^g
      gam[[length(gam) + 1L]] <- data.frame(run = rn, truth = tr, comparator = m, gamma = g,
        pair_weighted_mean_d = sum(w * d) / sum(w)) }
  }
}
res <- do.call(rbind, rows)
res$p_mean_holm <- NA_real_; res$p_median_holm <- NA_real_
for (f in names(families)) { s <- res$family == f
  if (!any(s)) next
  ph <- stats::p.adjust(c(res$p_mean[s], res$p_median[s]), "holm", n = families[[f]]$n)
  res$p_mean_holm[s] <- ph[seq_len(sum(s))]; res$p_median_holm[s] <- ph[sum(s) + seq_len(sum(s))] }
utils::write.csv(res, file.path(opt$out, "macro_vs_attention.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, levels), file.path(opt$out, "macro_levels.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, strata), file.path(opt$out, "macro_strata.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, gam), file.path(opt$out, "macro_pair_weighted.csv"), row.names = FALSE)

say("comparator - attention, per-regulon AP (negative = attention better); Holm within the fixed families")
for (k in seq_len(nrow(res))) with(res[k, ], say(sprintf(
  "%-36s %-12s %-10s n=%3d  mean %+.4f [%+.4f, %+.4f] (overlap groups [%+.4f, %+.4f]) p_holm %.3g | median %+.4f p_holm %.3g | attention wins %.0f%%, tied %.0f%%",
  run, truth, comparator, n_regulons, mean_d, mean_lo, mean_hi, mean_lo_overlap, mean_hi_overlap, p_mean_holm,
  median_d, p_median_holm, 100 * attention_wins, 100 * tied)))
say("done")
