# method_global_tau.png: ONE temperature for all genes (exploratory variant) versus the per-gene temperature, on real M3D rows (00_compute_m3d.R).
# Run:  nice -n 19 Rscript analysis/method_figs/14_method_global_tau.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); ri <- ex$row_info; H <- 7.2
bars2 <- function(ax, v, col, w = 0.7) for (q in seq_along(v)) polyline(ax$X(c(q, q)), c(ax$Y(0), ax$Y(v[q])), col, 1.8)
draw <- function() {
  txt("Global temperature: one tau for every gene (exploratory variant)", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  badge_exploratory(MF_W - 0.1, Y(0.16))
  # strip: what is shared, what changes
  steps <- c("CLR scores", "top 50 positive\nper gene", "softmax", "lazy walk P^9", "second CLR")
  xs <- 0.1 + (0:4) * 1.38
  for (k in 1:5) { hl <- k == 3; rect(xs[k], Y(1.05), 1.22, 0.55, fill = if (hl) "#c9c2f4" else "#f1f0ed", col = if (hl) VIT else NA, lwd = 1.2)
    lines_txt(strsplit(steps[k], "\n")[[1]], xs[k] + 0.61, Y(0.78 + (length(strsplit(steps[k], "\n")[[1]]) == 2) * -0.0) + 0.0 - 0.0, 0.15, 7, if (hl) VIT else ink[["secondary"]], "centre", if (hl) "bold" else "plain") }
  for (k in 1:4) harrow(xs[k] + 1.23, Y(0.775), xs[k + 1] - 0.01, Y(0.775), ink[["muted"]], 1)
  lines_txt(c(sprintf("per-gene (aCLR): tau solved for each row so that exp(entropy) = 10   |   global: ONE tau = %.4f, solved so that the MEAN", ex$tau_v2), "effective partners over all genes = 10; everything else is unchanged (same candidates, same walk, same second CLR)"), 0.1, Y(1.3), 0.15, 7, ink[["secondary"]])
  # panel A: weights on four real rows
  ya <- 1.85; txt("Same four genes, same scores, two temperature rules", 0.1, Y(ya), 8, face = "bold", just = "left")
  txt("weight on each of the 50 candidates, by rank (1 to 50)", 0.1, Y(ya + 0.18), 7, ink[["secondary"]], just = "left")
  ord <- names(ri)[order(vapply(ri, `[[`, 0, "tau"))]; mh <- 0.55
  rowspec <- list(list(top = 0.6, key = "w_loc", col = "#8c80e0", lab = "per-gene tau"), list(top = 1.95, key = "w_glob", col = "#3a2aa0", lab = "one global tau"))
  for (rs in rowspec) txt(rs$lab, 0.13, Y(ya + rs$top + mh / 2 + 0.1), 7, if (rs$key == "w_loc") VIT else "#3a2aa0", rot = 90, face = "bold")
  for (k in seq_along(ord)) { r <- ri[[ord[k]]]; xx <- 0.5 + (k - 1) * 0.82
    txt(r$gene, xx + 0.25, Y(ya + 0.42), 7.5, face = "bold")
    for (rs in rowspec) { yb_ <- Y(ya + rs$top + mh)
      a1 <- axes_box(xx, xx + 0.55, yb_, mh, c(0, 51), c(0, 1), c(1, 50), c(0, 0.5, 1), yfmt = function(v) sub("^0", "", sprintf("%g", v)))
      bars2(a1, r[[rs$key]], rs$col)
      tau_ <- if (rs$key == "w_loc") r$tau else ex$tau_v2; e_ <- if (rs$key == "w_loc") r$eff_loc else r$eff_glob
      txt(sprintf("tau %.2f", tau_), xx + 0.27, yb_ - 0.3, 7, if (rs$key == "w_loc") VIT else "#3a2aa0"); txt(sprintf("exp(H) %.1f", e_), xx + 0.27, yb_ - 0.45, 7, if (rs$key == "w_loc") VIT else "#3a2aa0") } }
  # panel B: ECDF of effective partners
  xb <- 3.95; txt("Effective partners per gene, all genes", xb, Y(ya), 8, face = "bold", just = "left")
  txt("exp(entropy) of each gene's attention row", xb, Y(ya + 0.18), 7, ink[["secondary"]], just = "left")
  ax <- axes_box(xb + 0.5, xb + 2.9, Y(ya + 0.5 + 1.7), 1.7, log10(c(1, 50)), c(0, 1), log10(c(1, 2, 5, 10, 20, 50)), c(0, 0.25, 0.5, 0.75, 1), xfmt = function(v) sprintf("%g", round(10^v)), yfmt = function(v) sub("^0", "", sprintf("%g", v)), xlab = "effective partners", ylab = "share of genes at or below")
  eg <- sort(ex$eff_glob[!is.na(ex$eff_glob)]); eg <- pmin(eg, 50); yy <- seq_along(eg) / length(eg)
  polyline(ax$X(log10(pmax(eg, 1))), ax$Y(yy), "#3a2aa0", 1.6)
  el <- ex$eff_loc; polyline(ax$X(log10(c(min(el), median(el), median(el), max(el)))), ax$Y(c(0, 0, 1, 1)), "#8c80e0", 2.2)
  lines_txt(c("per-gene tau:", "all genes at 10"), ax$X(log10(10)) + 0.1, ax$Y(0.3), 0.15, 7, VIT)
  lines_txt(c("one global tau:", "genes range", "widely"), ax$X(log10(2)) - 0.05, ax$Y(0.93), 0.15, 7, "#3a2aa0")
  qg <- quantile(ex$eff_glob, c(.1, .5, .9), na.rm = TRUE)
  lines_txt(c(sprintf("global: 10th / 50th / 90th pct %.1f / %.1f / %.1f;", qg[1], qg[2], qg[3]), sprintf("range %.1f to %.1f; mean %.1f", min(eg), max(eg), mean(ex$eff_glob, na.rm = TRUE)),
              sprintf("per-gene: all %d rows at %.1f", ex$G, mean(el))), xb, Y(ya + 0.5 + 1.7 + 0.55), 0.15, 7, ink[["secondary"]])
  # panel C: distribution of per-gene tau vs global
  yc <- 5.25; txt("One global value versus the per-gene values", 0.1, Y(yc), 8, face = "bold", just = "left")
  lt <- log2(ex$tau_loc); lo <- -2; hi <- 3
  hp <- hist_panel(0.7, 3.4, Y(yc + 0.3 + 1.1), 1.0, lt, c(lo, hi), seq(lo, hi, by = 0.25), fill = "#c9c2f4", ticks = -2:3, tick_fmt = function(v) sprintf("%d", v), xlab = "log2 of the temperature")
  seg(hp$X(log2(ex$tau_v2)), Y(yc + 0.3 + 1.1), hp$X(log2(ex$tau_v2)), Y(yc + 0.3 + 1.1) + 1.1, "#3a2aa0", 1.6)
  txt(sprintf("global tau %.4f", ex$tau_v2), hp$X(log2(ex$tau_v2)) - 0.06, Y(yc + 0.3 + 1.1) + 1.12, 7, "#3a2aa0", just = "right")
  txt(sprintf("%d per-gene taus", sum(!is.na(ex$tau_loc))), hp$X(1.2), Y(yc + 0.3 + 1.1) + 0.75, 7, VIT, just = "left")
  e <- ex$eco_v2
  lines_txt(c("EcoMAC, from the stored table (not redrawn):", sprintf("global tau %.4f; effective partners", e$tau), sprintf("10th / 50th / 90th pct %.1f / %.1f / %.1f,", e$eff_q10, e$eff_median, e$eff_q90), sprintf("range %.1f to %.1f.", e$eff_min, e$eff_max),
              "Descriptive: no comparison of network", "quality is shown here."), 3.95, Y(yc + 0.4), 0.16, 7, ink[["secondary"]])
}
mf_save("global_tau", H, draw)
cat(sprintf("tau_v2 %.4f; global eff quantiles %s\n", ex$tau_v2, paste(round(quantile(ex$eff_glob, c(0, .1, .5, .9, 1), na.rm = TRUE), 2), collapse = " ")))
