#!/usr/bin/env Rscript
# Sensitivity analysis (NOT preregistered): collects results/ecomac_nonoverlap/<subset>/{metrics,jack_replicates}.csv
# into one table. Intervals as baselines.R: point +- 1.96 sd over the 100 shared delete-half gene replicates;
# the ratio attention / CLR uses the same replicates (ratio of the replicate AUPRs).
subs <- c("A", "B", "C", "Bg3154", "Bmix", "R1370", "A_tstar", "B_tstar", "Bg3154_tstar", "Bmix_tstar")
od <- "results/ecomac_nonoverlap"
rows <- list()
for (s in subs) {
  f <- file.path(od, s, "metrics.csv"); if (!file.exists(f)) next
  m <- utils::read.csv(f); jr <- utils::read.csv(file.path(od, s, "jack_replicates.csv"))
  for (tr in c("SC", "all")) for (hs in c(FALSE, TRUE)) {
    a <- paste0("attention", if (hs) "_hubstrat" else ""); c0 <- paste0("clr", if (hs) "_hubstrat" else "")
    g <- function(meth, col) m[[col]][m$method == meth & m$truth == tr]
    ja <- jr$aupr[jr$method == a & jr$truth == tr]; jc <- jr$aupr[jr$method == c0 & jr$truth == tr]
    rat <- g(a, "aupr") / g(c0, "aupr"); jrat <- ja / jc
    rows[[length(rows) + 1L]] <- data.frame(subset = s, truth = tr, hub_controlled = hs, base_rate = g(a, "base_rate"),
      aupr_abs_pearson = if (hs) NA else g("abs_pearson", "aupr"), aupr_clr = g(c0, "aupr"), aupr_attention = g(a, "aupr"),
      aupr_genie3 = if (hs || !"genie3" %in% m$method) NA else g("genie3", "aupr"),
      diff_att_minus_clr = g(a, "d_vs_clr"), diff_lo = g(a, "d_clr_lo"), diff_hi = g(a, "d_clr_hi"),
      ratio_att_over_clr = rat, ratio_lo = rat - 1.96 * stats::sd(jrat), ratio_hi = rat + 1.96 * stats::sd(jrat),
      tf_aupr_ratio_clr = if (hs) NA else g("clr", "tf_aupr_ratio"), tf_aupr_ratio_att = if (hs) NA else g("attention", "tf_aupr_ratio"))
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(od, "summary_metrics.csv"), row.names = FALSE)
options(width = 250); print(format(res, digits = 3), row.names = FALSE)

## paired contrasts of the attention - CLR AUPR difference between disjoint subsets (same 100 gene halves,
## identical gene sets): this resamples genes only, not arrays, so it is descriptive.
contrast <- function(s1, s2) do.call(rbind, lapply(c("SC", "all"), function(tr) do.call(rbind, lapply(c(FALSE, TRUE), function(hs) {
  a <- paste0("attention", if (hs) "_hubstrat" else ""); c0 <- paste0("clr", if (hs) "_hubstrat" else "")
  dd <- function(s) { m <- utils::read.csv(file.path(od, s, "metrics.csv")); jr <- utils::read.csv(file.path(od, s, "jack_replicates.csv"))
    list(d = m$aupr[m$method == a & m$truth == tr] - m$aupr[m$method == c0 & m$truth == tr],
         j = jr$aupr[jr$method == a & jr$truth == tr] - jr$aupr[jr$method == c0 & jr$truth == tr]) }
  x <- dd(s1); y <- dd(s2); d <- x$d - y$d; se <- stats::sd(x$j - y$j)
  data.frame(contrast = paste(s1, "-", s2), truth = tr, hub_controlled = hs, d1 = x$d, d2 = y$d, diff = d,
             lo = d - 1.96 * se, hi = d + 1.96 * se) }))))
ct <- rbind(contrast("B", "C"), contrast("Bg3154", "Bmix"), contrast("B", "A"), contrast("B", "R1370"))
utils::write.csv(ct, file.path(od, "summary_contrasts.csv"), row.names = FALSE)
print(format(ct, digits = 3), row.names = FALSE)
