#!/usr/bin/env Rscript
# Summaries of sim_asymmetry.R outputs (synthetic only): means over seeds, and
# paired differences against a3r at t = 9 (same simulated data per seed).
# Writes out/summary_<scenario>.csv, out/paired_vs_a3r.csv, out/per_clique_wide.csv.
D <- "analysis/synthetic_asymmetry/out"
ops <- c("absr", "CLR", "att_t9", "a3r_t9", "C1_a3c_t9", "C2_sink_att_t9", "C2_sink_a3r_t9",
         "C3_lcclr", "C3_lcclr_t4", "C4a_resid_att_t9", "C4b_units_first",
         "att_t4", "a3r_t4", "C1_a3c_t4", "C2_sink_a3r_t4")
met <- c("cpx_x", "cpx_noribo_x", "reg_x", "edge_x", "sat_auc", "hub_auc", "sat_in_topK", "dir_target_to_tf")
options(width = 220)
pair <- list(); pcw <- list()
for (sc in c("base", "global", "tf_matched")) {
  m <- utils::read.csv(file.path(D, sc, "metrics_by_seed.csv"))
  mm <- stats::aggregate(m[met], by = list(op = m$op), FUN = mean)
  mm <- mm[match(ops, mm$op), ]
  utils::write.csv(mm, file.path(D, paste0("summary_", sc, ".csv")), row.names = FALSE)
  cat("\n==== scenario", sc, "(mean over", length(unique(m$seed)), "seeds) ====\n")
  print(format(mm, digits = 3), row.names = FALSE)
  ref <- m[m$op == "a3r_t9", ]
  for (o in setdiff(ops, "a3r_t9")) { x <- m[m$op == o, ]; x <- x[match(ref$seed, x$seed), ]
    for (v in c("cpx_x", "cpx_noribo_x", "reg_x", "edge_x", "sat_auc", "hub_auc")) {
      d <- x[[v]] / ref[[v]]
      pair[[length(pair) + 1]] <- data.frame(scenario = sc, op = o, metric = v, ratio_mean = mean(d),
        ratio_min = min(d), ratio_max = max(d), seeds_better = sum(d > 1), seeds = length(d)) } }
  p <- utils::read.csv(file.path(D, sc, "per_clique_ap.csv"))
  p <- p[p$op %in% c("absr", "CLR", "att_t9", "a3r_t9", "C1_a3c_t9", "C2_sink_a3r_t9", "C3_lcclr_t4", "C4b_units_first"), ]
  p$cell <- sprintf("L%.2f_n%s", p$loading, sub(" .*", "", p$size_class))
  w <- stats::reshape(p[, c("op", "cell", "ap_clique")], idvar = "op", timevar = "cell", direction = "wide")
  names(w) <- sub("ap_clique.", "", names(w), fixed = TRUE)
  cells <- unique(p$cell[order(p$loading, as.numeric(sub(" .*", "", p$size_class)))])
  w <- w[, c("op", cells)]; pcw[[sc]] <- cbind(scenario = sc, w)
  if (sc == "base") { cat("\nper-clique AP (base):\n"); print(format(w, digits = 2), row.names = FALSE) }
}
pair <- do.call(rbind, pair)
utils::write.csv(pair, file.path(D, "paired_vs_a3r.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, pcw), file.path(D, "per_clique_wide.csv"), row.names = FALSE)
cat("\npaired ratio to a3r_t9 (mean [min, max] over seeds, seeds better):\n")
pw <- within(pair[pair$op %in% c("CLR", "att_t9", "C1_a3c_t9", "C2_sink_a3r_t9", "C2_sink_att_t9", "C3_lcclr", "C3_lcclr_t4", "C4b_units_first"), ],
             txt <- sprintf("%.3f [%.3f, %.3f] %d/%d", ratio_mean, ratio_min, ratio_max, seeds_better, seeds))
print(stats::reshape(pw[, c("scenario", "op", "metric", "txt")], idvar = c("scenario", "op"), timevar = "metric", direction = "wide"),
      row.names = FALSE)
