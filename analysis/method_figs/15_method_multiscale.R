# method_multiscale.png: multi-scale scores MS1/MS2/MS3 over the 13-temperature ladder, on real M3D pairs (00_compute_m3d.R: walk depth 9, second CLR at every temperature).
# Run:  nice -n 19 Rscript analysis/method_figs/15_method_multiscale.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); nm <- ex$names; n <- length(nm); lg <- ex$ladder_log2; H <- 6.85
stopifnot(length(lg) == 13, abs(ex$ms1[1, 11] - mean(ex$pct[1, 11, ])) < 1e-12, abs(ex$ms1_8[1, 11] - mean(ex$pct[1, 11, round(lg, 2) %in% round(ex$ms3_log2, 2)])) < 1e-12)
draw <- function() {
  txt("Multi-scale: average over temperatures, no tuning", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  badge_exploratory(MF_W - 0.1, Y(0.16))
  # ---- ladder strip ----
  panel_head(1, "Thirteen temperatures", 0.1, Y(0.55), "log2 tau from -4 to +2 in half octaves; at each, the same walk (t = 9) and second CLR give one score per pair")
  xs <- function(v) 1.6 + (v + 4) / 6 * 5.2; yl <- Y(1.3)
  seg(xs(-4), yl, xs(2), yl, ink[["secondary"]], 0.9)
  for (k in seq_along(lg)) { dot(xs(lg[k]), yl, "#8c80e0", 6); if (lg[k] == round(lg[k])) txt(sprintf("%+d", as.integer(lg[k])), xs(lg[k]), yl + 0.14, 7, ink[["secondary"]]) else seg(xs(lg[k]), yl, xs(lg[k]), yl + 0.04, ink[["secondary"]], 0.7)
    txt(sprintf("%.1f", ex$eff_tau[k]), xs(lg[k]), yl - 0.17, 7, VIT) }
  txt("log2 tau", 1.3, yl + 0.14, 7, ink[["secondary"]], just = "right"); txt("mean effective", 1.3, yl - 0.1, 7, VIT, just = "right"); txt("partners per gene", 1.3, yl - 0.24, 7, VIT, just = "right")
  m3 <- range(ex$ms3_log2); seg(xs(m3[1]), yl - 0.42, xs(m3[2]), yl - 0.42, MGT, 1.4); seg(xs(m3[1]), yl - 0.42, xs(m3[1]), yl - 0.36, MGT, 1.4); seg(xs(m3[2]), yl - 0.42, xs(m3[2]), yl - 0.36, MGT, 1.4)
  txt("MS3 uses these 8 temperatures", xs(mean(m3)), yl - 0.56, 7, MGT)
  # ---- pairs across the ladder ----
  yt <- 2.35; panel_head(2, "Each pair across the ladder", 0.1, Y(yt), "rank percentile of the pair's score among all 9.2 million pairs")
  a <- axes_box(0.75, 3.5, Y(yt + 0.55 + 1.65), 1.65, c(-4, 2), c(0, 1), -4:2, c(0, 0.25, 0.5, 0.75, 1), xfmt = function(v) sprintf("%+d", v), yfmt = function(v) sub("^0", "", sprintf("%g", v)), xlab = "log2 tau", ylab = "percentile")
  upp <- which(upper.tri(ex$ms1), arr.ind = TRUE)
  for (q in seq_len(nrow(upp))) polyline(a$X(lg), a$Y(ex$pct[upp[q, 1], upp[q, 2], ]), "#dcdbd7", 0.8)
  sel <- list(list(j = 3, col = "black", lw = 1.8, who = "lexA-dinI"), list(j = 10, col = "#3a2aa0", lw = 1.8, who = "lexA-rpoN"), list(j = 11, col = "#8c80e0", lw = 1.8, who = "lexA-yeaR"))
  for (s in sel) { polyline(a$X(lg), a$Y(ex$pct[1, s$j, ]), s$col, s$lw); m <- ex$ms1[1, s$j]; seg(a$X(2) , a$Y(m), a$X(2) + 0.06, a$Y(m), s$col, 1.5); dot(a$X(2) + 0.12, a$Y(m), s$col, 5, 18) }
  txt("lexA-dinI", a$X(-3.9), a$Y(1) - 0.1, 7, "black", just = "left")
  txt(sprintf("lexA-%s", nm[10]), a$X(-3.9), a$Y(ex$pct[1, 10, 1]) - 0.12, 7, "#3a2aa0", just = "left")
  txt(sprintf("lexA-%s", nm[11]), a$X(-3.9), a$Y(ex$pct[1, 11, 1]) - 0.12, 7, "#6250d6", just = "left")
  txt("grey: the other pairs of the 11 genes", a$X(2), a$Y(0.04), 7, ink[["muted"]], just = "right")
  txt("diamond at right: MS1 (mean of the 13)", 0.75, Y(yt + 0.55 + 1.65 + 0.5), 7, ink[["secondary"]], just = "left")
  # ---- MS1 matrix ----
  xm <- 3.95; panel_head(3, "MS1 for the 11 genes", xm, Y(yt), "mean of the 13 percentiles")
  ms <- ex$ms1; diag(ms) <- 0
  mat_draw(xm + 0.36, Y(yt + 0.85), 0.215, ms, nm, ramp = RAMP$violet, vmax = 1, gamma = 1.4, d = 2, thresh = 0.7)
  # ---- definitions ----
  yd <- 6.1; seg(0.1, Y(yd - 0.25), 6.9, Y(yd - 0.25), ink[["grid"]], 0.8)
  defs <- c("MS1: mean over the 13 temperatures of each pair's rank percentile (average-tie rank / number of pairs, within that temperature).",
            "MS2: mean over the 13 temperatures of each pair's score standardised within the temperature (z over all pairs).",
            "MS3: as MS1 over the 8 temperatures log2 tau = -2.5 to +1.  Equal weights; no temperature is selected, tuned or dropped.")
  lines_txt(defs, 0.1, Y(yd), 0.17, 7, ink[["secondary"]])
}
mf_save("multiscale", H, draw)
cat(sprintf("MS1 lexA-%s %.3f lexA-%s %.3f; MS3 %.3f %.3f; MS2 %.2f %.2f\n", nm[10], ex$ms1[1, 10], nm[11], ex$ms1[1, 11], ex$ms1_8[1, 10], ex$ms1_8[1, 11], ex$ms2[1, 10], ex$ms2[1, 11]))
