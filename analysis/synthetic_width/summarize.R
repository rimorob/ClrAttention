#!/usr/bin/env Rscript
# Summaries for sim_width.R: mean (SD) over seeds of each metric, and paired
# ratios to the method (K10, attention readout, t = 9) with the number of
# seeds in which the operator beats it. SYNTHETIC ONLY.
OUT <- "analysis/synthetic_width/out"
mets <- c("cpx_x", "cpx_noribo_x", "s5_x", "s20_x", "s80_x", "small_x", "mid_x", "reg_x", "sat_auc")
for (scen in c("messy", "clean")) {
  f <- file.path(OUT, scen, "metrics_by_seed.csv"); if (!file.exists(f)) next
  r <- utils::read.csv(f)
  agg <- do.call(rbind, lapply(split(r, factor(r$op, unique(r$op))), function(d)
    data.frame(op = d$op[1], n_seeds = nrow(d), t(sapply(mets, function(v)
      if (v == "sat_auc") sprintf("%.3f (%.3f)", mean(d[[v]]), stats::sd(d[[v]])) else
        sprintf("%.1f (%.1f)", mean(d[[v]]), stats::sd(d[[v]])))))))
  utils::write.csv(agg, file.path(OUT, sprintf("summary_%s.csv", scen)), row.names = FALSE)
  ref <- r[r$op == "K10|att_t9", ]
  pr <- do.call(rbind, lapply(split(r, factor(r$op, unique(r$op))), function(d) {
    d <- d[match(ref$seed, d$seed), ]
    data.frame(op = d$op[1], t(sapply(setdiff(mets, "sat_auc"), function(v) {
      q <- d[[v]] / ref[[v]]; sprintf("%.2f (%d/%d)", exp(mean(log(q))), sum(q > 1), length(q)) })))
  }))
  utils::write.csv(pr, file.path(OUT, sprintf("paired_vs_K10_%s.csv", scen)), row.names = FALSE)
  cat("\n==== ", scen, ": mean (SD) over seeds, x base rate ====\n", sep = ""); print(agg, row.names = FALSE)
  cat("\n==== ", scen, ": geometric-mean ratio to K10|att_t9 (seeds better) ====\n", sep = ""); print(pr, row.names = FALSE)
}
