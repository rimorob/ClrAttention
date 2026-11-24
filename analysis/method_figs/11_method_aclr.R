# method_aclr.png: aCLR = CLR scores -> top-50 -> softmax with a per-gene temperature (keff 10) -> row-stochastic A -> lazy walk -> symmetric read-out.
# Real M3D numbers (00_compute_m3d.R). Run:  nice -n 19 Rscript analysis/method_figs/11_method_aclr.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); nm <- ex$names; n <- length(nm); ri <- ex$row_info
H <- 6.95; CE <- 0.145
bars <- function(ax, v, col, w = 0.7) for (k in seq_along(v)) polyline(ax$X(c(k, k)), c(ax$Y(0), ax$Y(v[k])), col, 1.7)
draw <- function() {
  txt("aCLR: CLR scores as attention, one temperature per gene", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  badge_example(MF_W - 0.1, Y(0.16), sprintf("Example: E. coli M3D (%d genes x %d profiles); query gene lexA", ex$G, ex$N))
  top <- 0.6; L <- ri[["lexA"]]; yb <- Y(top + 1.55)
  x1 <- 0.1; panel_head(1, "Top 50", x1, Y(top), "CLR score by rank")
  ax <- axes_box(x1 + 0.35, x1 + 1.2, yb, 0.95, c(0, 51), c(0, 16), c(1, 25, 50), c(0, 5, 10, 15), xlab = "rank of partner")
  bars(ax, L$score, "#f6b595")
  lines_txt(c("lexA row: the 50 highest", sprintf("positive of %d scores", ex$G - 1)), x1, Y(top + 2.07), 0.15, 7, ink[["secondary"]])
  x2 <- 1.47; panel_head(2, "Softmax", x2, Y(top), "weight by rank")
  ax2 <- axes_box(x2 + 0.35, x2 + 1.2, yb, 0.95, c(0, 51), c(0, 0.5), c(1, 25, 50), c(0, 0.25, 0.5), yfmt = function(v) sub("^0", "", sprintf("%g", v)), xlab = "rank of partner")
  bars(ax2, L$w_loc, "#8c80e0")
  lines_txt(c(sprintf("tau %.2f is solved so that", L$tau), sprintf("exp(entropy) = %.1f", L$eff_loc)), x2, Y(top + 2.07), 0.15, 7, VIT)
  x3 <- 2.82; panel_head(3, "Attention A", x3, Y(top), "row-stochastic, x 100")
  mat_draw(x3 + 0.33, Y(top + 0.62), CE, ex$A_sub, nm, ramp = RAMP$violet, vmax = 0.45, gamma = 0.8, fmt = function(v) sprintf("%d", round(v * 100)))
  txt(sprintf("%d of %d genes shown", n, ex$G), x3 + 0.33, Y(top + 0.62 + n * CE + 0.14), 7, ink[["secondary"]], just = "left")
  x4 <- 4.9; panel_head(4, "Walk, symmetric", x4, Y(top), "(P^9 + P^9') / 2, x 100")
  mat_draw(x4 + 0.33, Y(top + 0.62), CE, ex$aclr_sub * 100, nm, ramp = RAMP$magenta, vmax = 12, gamma = 0.8, fmt = function(v) sprintf("%d", round(v)))
  txt("P = 0.5 I + 0.5 A, depth t = 9", x4 + 0.33, Y(top + 0.62 + n * CE + 0.14), 7, MGT, just = "left")
  txt("fixed: effective partners keff 10, cap 50 candidates per gene, lazy weight 0.5, walk depth t = 9", 0.1, Y(top + 2.6), 7, ink[["secondary"]], just = "left")
  rb <- 3.4
  txt("Local temperature: rows with different score spreads get different tau, the same number of effective partners", 0.1, Y(rb), 8, face = "bold", just = "left")
  ord <- names(ri)[order(vapply(ri, `[[`, 0, "tau"))]
  yb1 <- Y(rb + 1.3); yb2 <- Y(rb + 2.55); mh <- 0.75
  for (k in seq_along(ord)) { r <- ri[[ord[k]]]; xx <- 0.1 + (k - 1) * 1.72
    txt(sprintf("%s: tau %.2f", r$gene, r$tau), xx + 0.4, Y(rb + 0.35), 7.5, VIT, "bold", just = "left")
    a1 <- axes_box(xx + 0.4, xx + 1.55, yb1, mh, c(0, 51), c(0, 16), c(1, 25, 50), c(0, 5, 10, 15), ylab = if (k == 1) "CLR score" else NULL)
    for (q in seq_along(r$score)) polyline(a1$X(c(q, q)), c(a1$Y(0), a1$Y(r$score[q])), "#f6b595", 1.6)
    a2 <- axes_box(xx + 0.4, xx + 1.55, yb2, mh, c(0, 51), c(0, 0.5), c(1, 25, 50), c(0, 0.25, 0.5), yfmt = function(v) sub("^0", "", sprintf("%g", v)), xlab = "rank", ylab = if (k == 1) "weight" else NULL)
    for (q in seq_along(r$w_loc)) polyline(a2$X(c(q, q)), c(a2$Y(0), a2$Y(r$w_loc[q])), "#8c80e0", 1.6)
    txt(sprintf("weight on rank 1: %.2f", r$w_loc[1]), xx + 1.55, a2$Y(0.46), 7, ink[["secondary"]], just = "right")
    txt(sprintf("exp(entropy) = %.1f", r$eff_loc), xx + 1.55, a2$Y(0.37), 7, VIT, just = "right") }
  lines_txt(c(sprintf("rlmL, flgJ and napF sit at the 10th, 50th and 90th percentile of tau over all %d genes (%.2f / %.2f / %.2f); lexA is shown for continuity.", ex$G, quantile(ex$tau_loc, .1), quantile(ex$tau_loc, .5), quantile(ex$tau_loc, .9)),
              "tau follows the spread of the scores: a flat row gets a small tau, a steep row a large one, so each row attends to about 10 partners."), 0.1, Y(rb + 3.2), 0.16, 7, ink[["secondary"]])
}
mf_save("aclr", H, draw)
cat(sprintf("eff_loc: min %.4f max %.4f; tau range %.3f-%.3f\n", min(ex$eff_loc), max(ex$eff_loc), min(ex$tau_loc, na.rm = TRUE), max(ex$tau_loc, na.rm = TRUE)))
