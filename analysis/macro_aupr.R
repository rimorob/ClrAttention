#!/usr/bin/env Rscript
# Macro-averaged per-regulon AP: every regulon counts once, next to the pooled
# (micro) co-membership AUPR that large regulons dominate. Reads per_regulon.csv
# and per_regulon_jack.csv written by analysis/baselines.R. Definition fixed
# with the statistics reviewer before scoring (docs/paper_plan.md, theme 5):
#   - candidates as regulon_coherence(); AP with ties broken in expectation;
#     regulons with >= 10 positive (non-same-operon) pairs in the full data
#     (so CORUM complexes of 3-4 genes drop out); the set depends on the truth
#     only. A gene in several regulons counts in each (no merging).
#   - comparisons only within one run (one results directory); duplicated
#     (run, method, truth, regulon) rows or missing inputs stop the script
#   - per-regulon paired difference dAP = AP(comparator) - AP(attention);
#     mean and median over regulons. Two sources of uncertainty, both required:
#       regulons: BCa bootstrap over regulons for the mean (p by inverting the
#                 interval); Wilcoxon signed rank for the median
#       genes:    delete-d gene jackknife, d = 20% of genes (primary; random
#                 subsets, variance factor (n - d)/d = 4; Shao & Tu 1995 5.2),
#                 200 replicates; the regulon set is fixed per design to the
#                 regulons keeping >= 3 positive pairs in every replicate, and
#                 the gene-interval estimate is computed on that same set
#                 (number excluded reported); d = 50% as a stability check,
#                 flagged when the 95% interval of log(SE ratio) lies wholly
#                 outside +-log(1.25)
#     A difference counts only if both intervals exclude 0; Holm on the larger p.
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
            out = "results/macro_aupr", B = 10000L, seed = 20261017L)
i <- 1L
while (i <= length(args)) { k <- sub("^--", "", args[i])
  opt[[k]] <- if (k %in% c("B", "seed")) as.integer(args[i + 1L]) else args[i + 1L]; i <- i + 2L }
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
f1 <- file.path(dirs, "per_regulon.csv"); f2 <- file.path(dirs, "per_regulon_jack.csv")
miss <- c(f1, f2)[!file.exists(c(f1, f2))]
if (length(miss)) stop("missing input: ", paste(miss, collapse = ", "))
pr <- do.call(rbind, lapply(f1, utils::read.csv, stringsAsFactors = FALSE))
pj <- do.call(rbind, lapply(f2, utils::read.csv, stringsAsFactors = FALSE))
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
wil <- function(d) if (sum(d != 0) >= 5) stats::wilcox.test(d, exact = FALSE)$p.value else NA_real_
# delete-d gene jackknife of a statistic f over the per-regulon differences;
# regulon set fixed to `regs`; returns SE and the share of (regulon, replicate) cells skipped
gene_jack <- function(jc, ja, d_full, regs, design, f) {
  factor <- c(d20 = 0.8 / 0.2, d50 = 1)[[design]]
  jc <- jc[jc$design == design & jc$regulon %in% regs, ]; ja <- ja[ja$design == design & ja$regulon %in% regs, ]
  m <- merge(jc[, c("rep", "regulon", "ap")], ja[, c("rep", "regulon", "ap")], by = c("rep", "regulon"), suffixes = c(".c", ".a"))
  R <- length(unique(jc$rep))
  okc <- tapply(!is.na(m$ap.c) & !is.na(m$ap.a), m$regulon, sum)
  fixed <- names(okc)[okc == R]                                   # scored in every replicate
  m <- m[m$regulon %in% fixed, ]
  th <- vapply(split(m$ap.c - m$ap.a, m$rep), f, 0)
  list(se = sqrt(factor * mean((th - mean(th))^2)), est = f(d_full[regs %in% fixed]),
       n = length(fixed), excluded = length(regs) - length(fixed), R = R)
}
unstable <- function(j50, j20) { if (!(j50$se > 0 && j20$se > 0)) return(NA)
  l <- log(j50$se / j20$se); s <- sqrt(1 / j50$R + 1 / j20$R) / sqrt(2)
  abs(l) - 1.96 * s > log(1.25) }
zp <- function(est, se) if (is.finite(se) && se > 0) 2 * stats::pnorm(-abs(est) / se) else NA_real_

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
  pjr <- pj[pj$run == rn & pj$truth == tr, ]; ja <- pjr[pjr$method == "attention", ]
  for (m in setdiff(unique(x$method), "attention")) {
    y <- x[x$method == m, ]; j <- match(att$regulon, y$regulon); ok <- !is.na(j)
    a <- att[ok, ]; d <- y$ap[j[ok]] - a$ap                     # comparator - attention
    jc <- pjr[pjr$method == m, ]
    bm <- boot_mean(d); pw <- wil(d)
    jm20 <- gene_jack(jc, ja, d, a$regulon, "d20", mean); jm50 <- gene_jack(jc, ja, d, a$regulon, "d50", mean)
    jd20 <- gene_jack(jc, ja, d, a$regulon, "d20", stats::median)
    pjm <- zp(jm20$est, jm20$se); pjd <- zp(jd20$est, jd20$se)
    noop <- a$frac_same_operon <= 0.5; bn <- boot_mean(d[noop])
    rows[[length(rows) + 1L]] <- data.frame(family = family_of(rn), run = rn, dataset = ds, truth = tr, comparator = m,
      n_regulons = length(d),
      mean_d = mean(d), mean_lo_regulons = bm[["lo"]], mean_hi_regulons = bm[["hi"]], p_mean_regulons = bm[["p"]],
      n_regulons_genes = jm20$n, excluded_genes = jm20$excluded, mean_d_genes = jm20$est,
      mean_lo_genes = jm20$est - 1.96 * jm20$se, mean_hi_genes = jm20$est + 1.96 * jm20$se, p_mean_genes = pjm,
      p_mean = max(bm[["p"]], pjm),
      median_d = stats::median(d), p_median_regulons = pw,
      median_d_genes = jd20$est, median_lo_genes = jd20$est - 1.96 * jd20$se, median_hi_genes = jd20$est + 1.96 * jd20$se,
      p_median_genes = pjd, p_median = max(pw, pjd),
      jack_se_ratio_d50_d20 = jm50$se / jm20$se, jack_unstable = unstable(jm50, jm20), excluded_genes_d50 = jm50$excluded,
      attention_wins = mean(d < 0), tied = mean(d == 0),
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

say("comparator - attention, per-regulon AP (negative = attention better); both regulon and gene intervals shown; Holm on the larger p")
for (k in seq_len(nrow(res))) with(res[k, ], say(sprintf(
  "%-36s %-12s %-10s n=%3d  mean %+.4f regulons [%+.4f, %+.4f] genes [%+.4f, %+.4f] p_holm %.3g | median %+.4f p_holm %.3g | attention wins %.0f%%, tied %.0f%%%s",
  run, truth, comparator, n_regulons, mean_d, mean_lo_regulons, mean_hi_regulons, mean_lo_genes, mean_hi_genes, p_mean_holm,
  median_d, p_median_holm, 100 * attention_wins, 100 * tied, if (isTRUE(jack_unstable)) "  [gene jackknife unstable: d50/d20 SE ratio outside +-25%]" else "")))
say("done")
