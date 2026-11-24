# method_tcaclr.png: tcaCLR = aCLR walk P^9, then a second CLR on log P^9 over reachable pairs. Real eleven-gene M3D neighbourhood (00_compute_m3d.R).
# Run:  nice -n 19 Rscript analysis/method_figs/12_method_tcaclr.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); nm <- ex$names; n <- length(nm); J <- 3L; BG <- n
zi <- ex$Z2_sub[1, J]; zj <- ex$Z2_sub[J, 1]; sij <- ex$S2_sub[1, J]
stopifnot(abs((max(zi, 0) + max(zj, 0)) / sqrt(2) - sij) < 1e-9)
# z of the traced entries recomputed from the stored log P^9, mu and sd of the rows
stopifnot(abs((ex$logP_sub[1, J] - ex$mu2[1]) / ex$sd2[1] - zi) < 1e-9)
CELL <- 0.225; H <- 8.3
draw <- function() {
  txt("tcaCLR: a second CLR on the walk matrix", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  badge_example(MF_W - 0.1, Y(0.16), sprintf("Example: %d genes of E. coli M3D (%d genes x %d profiles)", n, ex$G, ex$N))
  xl <- 0.1; xr <- 3.65; r1 <- 0.55; r2 <- r1 + 3.75
  panel_head(1, "Walk matrix P^9 (aCLR operator)", xl, Y(r1), "entry (i, j): chance a 9-step lazy walk from i is at j, x 100")
  g1 <- mat_draw(xl + 0.42, Y(r1 + 0.82), CELL, ex$P9_sub * 100, nm, ramp = RAMP$magenta, vmax = 18, gamma = 0.8, d = 1, diag_grey = FALSE)
  lines_txt(c("rows sum to 1 over all 4297 genes"), xl + 0.42, Y(r1 + 0.82 + n * CELL + 0.14), 0.16, 7, ink[["secondary"]])
  panel_head(2, "Each gene against its own walk background", xr, Y(r1), "row z-score of log P^9 over the pairs the walk reaches")
  v <- ex$logp_row_i[is.finite(ex$logp_row_i)]; lo <- floor(min(v)); hi <- ceiling(max(v))
  xr0 <- xr + 0.55; xr1 <- xr + 3.3; yb <- Y(r1 + 2.35)
  hp <- hist_panel(xr0, xr1, yb, 1.05, v, c(lo, hi), seq(lo, hi, by = 0.5), ticks = seq(ceiling(lo / 4) * 4, hi, by = 4), tick_fmt = function(t) sprintf("%d", t), xlab = "log P^9 from lexA to each reached gene")
  txt("number of genes", xr0 - 0.12, yb + 0.5, 7, ink[["secondary"]], rot = 90)
  mu <- ex$mu2[1]; s <- ex$sd2[1]
  seg(hp$X(mu), yb, hp$X(mu), yb + 1.2, ink[["secondary"]], 1); txt(sprintf("mean %.1f", mu), hp$X(mu) - 0.05, yb + 1.24, 7, ink[["secondary"]], just = "right")
  seg(hp$X(mu), yb + 0.7, hp$X(mu + s), yb + 0.7, ink[["primary"]], 1.4); txt(sprintf("sd %.1f", s), hp$X(mu + s) + 0.05, yb + 0.7, 7, ink[["primary"]], just = "left")
  lv <- ex$logP_sub[1, J]; seg(hp$X(lv), yb, hp$X(lv), yb + 0.3, "black", 1.1); dot(hp$X(lv), yb, "black", 5)
  lines_txt(c(sprintf("%s:", nm[J]), sprintf("log P^9 %.1f, z = %.1f", lv, ex$Z2_sub[1, J])), hp$X(lv) + 0.08, yb + 0.55, 0.15, 7, "black", just = "right")
  lv <- ex$logP_sub[1, BG]; seg(hp$X(lv), yb, hp$X(lv), yb + 1.42, ink[["secondary"]], 1.1, 3); dot(hp$X(lv), yb, ink[["secondary"]], 5)
  txt(sprintf("%s (dotted): log P^9 %.1f, z = %.1f", nm[BG], lv, ex$Z2_sub[1, BG]), xr1, yb + 1.5, 7, ink[["secondary"]], just = "right")
  lines_txt(c(sprintf("lexA reaches %d of %d other genes;", ex$reach_i, ex$G - 1),
              "unreached pairs (P^9 = 0) stay out of the", "background and score 0"), xr + 0.55, Y(r1 + 3.0), 0.16, 7, ink[["secondary"]])
  panel_head(3, "Clip at 0", xl, Y(r2), "as in CLR: z below the row mean is set to 0")
  g3 <- mat_draw(xl + 0.42, Y(r2 + 0.82), CELL, ex$Z2_sub, nm, ramp = RAMP$orange, vmax = 6, gamma = 0.9, d = 1, neg_grey = TRUE)
  cell_box(g3, 1, J, "black", 1.6, 1); cell_box(g3, J, 1, "black", 1.6, 2); cell_box(g3, 1, BG, ink[["secondary"]], 1.6, 3)
  txt(sprintf("grey = z < 0, clipped to 0 (%d of %d off-diagonal cells)", sum(ex$Z2_sub[row(ex$Z2_sub) != col(ex$Z2_sub)] < 0), n * (n - 1)), xl + 0.42, Y(r2 + 0.82 + n * CELL + 0.14), 7, ink[["secondary"]], just = "left")
  panel_head(4, "Symmetric Stouffer score", xr, Y(r2), "S2_ij = (max(z_ij, 0) + max(z_ji, 0)) / sqrt(2)")
  g4 <- mat_draw(xr + 0.42, Y(r2 + 0.82), CELL, ex$S2_sub, nm, ramp = RAMP$orange, vmax = 8, gamma = 0.9, d = 1)
  cell_box(g4, 1, J, "black", 1.6, 1); cell_box(g4, J, 1, "black", 1.6, 2); cell_box(g4, 1, BG, ink[["secondary"]], 1.6, 3)
  lines_txt(c(sprintf("lexA-%s (solid): (%.1f + %.1f) / sqrt(2) = %.1f", nm[J], zi, zj, sij),
              sprintf("lexA-%s (dotted): (%.1f + %.1f) / sqrt(2) = %.1f", nm[BG], max(ex$Z2_sub[1, BG], 0), max(ex$Z2_sub[BG, 1], 0), ex$S2_sub[1, BG]),
              "symmetric: this is the tcaCLR score"), xr + 0.42, Y(r2 + 0.82 + n * CELL + 0.14), 0.16, 7, ink[["secondary"]])
  harrow(xl + 3.45, Y(r1 + 1.6), xr - 0.02, Y(r1 + 1.6), ink[["muted"]], 1.2)
  harrow(xl + 3.45, Y(r2 + 1.6), xr - 0.02, Y(r2 + 1.6), ink[["muted"]], 1.2)
}
mf_save("tcaclr", H, draw)
