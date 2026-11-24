# method_clr.png: CLR in four steps on a real eleven-gene neighbourhood of M3D (numbers from analysis/method_figs/00_compute_m3d.R).
# Run:  nice -n 19 Rscript analysis/method_figs/10_method_clr.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); nm <- ex$names; n <- length(nm)
J <- 3L                                   # traced pair: lexA (row 1) and dinI
BG <- n                                   # a background gene (clipping example)
zi <- ex$Zraw_sub[1, J]; zj <- ex$Zraw_sub[J, 1]; sij <- ex$S_sub[1, J]
stopifnot(abs((max(zi, 0) + max(zj, 0)) / sqrt(2) - sij) < 1e-9)
CELL <- 0.225; H <- 8.3
draw <- function() {
  txt("CLR: from mutual information to a context-normalised score", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  badge_example(MF_W - 0.1, Y(0.16), sprintf("Example: %d genes of E. coli M3D (%d genes x %d profiles)", n, ex$G, ex$N))
  lab <- function(r, cc) NULL
  xl <- 0.1; xr <- 3.65; r1 <- 0.55; r2 <- r1 + 3.75
  # ---- panel 1: MI ----
  panel_head(1, "B-spline mutual information", xl, Y(r1), "B-spline estimator, one value per gene pair")
  g1 <- mat_draw(xl + 0.42, Y(r1 + 0.82), CELL, ex$D_sub, nm, ramp = RAMP$blue, vmax = 0.6, gamma = 1, d = 2)
  txt("grey diagonal: not used", xl + 0.42, Y(r1 + 0.82 + n * CELL + 0.14), 7, ink[["muted"]], just = "left")
  # ---- panel 2: background distribution of one row ----
  panel_head(2, "Each gene against its own background", xr, Y(r1), "row z-score: (MI - row mean) / row sd, over all genes")
  v <- ex$mi_row_i[-ex$i]; xr0 <- xr + 0.55; xr1 <- xr + 3.3; yb <- Y(r1 + 2.35)
  hp <- hist_panel(xr0, xr1, yb, 1.05, v, c(0, 0.6), seq(0, 0.6, by = 0.015), ticks = c(0, 0.2, 0.4, 0.6), tick_fmt = function(t) sub("^0\\.", ".", sprintf("%g", t)),
                   xlab = sprintf("MI of lexA with each of the other %d genes", ex$G - 1))
  txt("number of genes", xr0 - 0.12, yb + 0.5, 7, ink[["secondary"]], rot = 90)
  bw <- 0.015; xs <- seq(0, 0.6, length.out = 200); cnt <- length(v) * bw * dnorm(xs, ex$mu_i, ex$sd_i)
  polyline(hp$X(xs), hp$Y(cnt), BLT, 1.1)
  mu <- ex$mu_i; s <- ex$sd_i; mi_j <- ex$D_sub[1, J]; mi_b <- ex$D_sub[1, BG]
  seg(hp$X(mu), yb, hp$X(mu), yb + 1.2, ink[["secondary"]], 1)
  txt(sprintf("mean %s", fnum(mu, 3)), hp$X(mu) + 0.05, yb + 1.24, 7, ink[["secondary"]], just = "left")
  seg(hp$X(mu), yb + 0.7, hp$X(mu + s), yb + 0.7, ink[["primary"]], 1.4); txt(sprintf("sd %s", fnum(s, 3)), hp$X(mu + s) + 0.05, yb + 0.7, 7, ink[["primary"]], just = "left")
  for (q in list(list(v = mi_j, who = nm[J], col = "black", ty = 1, ytop = 0.5), list(v = mi_b, who = nm[BG], col = ink[["secondary"]], ty = 3, ytop = 1.5))) {
    seg(hp$X(q$v), yb, hp$X(q$v), yb + q$ytop - 0.05, q$col, 1.1, q$ty); dot(hp$X(q$v), yb, q$col, 5)
    lb <- sprintf("%s: MI %s, z = %.1f", q$who, fnum(q$v, 3), (q$v - mu) / s)
    if (q$who == nm[J]) txt(lb, hp$X(q$v) - 0.06, yb + 0.62, 7, q$col, just = "right") else txt(lb, hp$X(q$v) - 0.02, yb + q$ytop + 0.03, 7, q$col, just = "left") }
  lines_txt(c(sprintf("z_ij = (MI_ij - mean_i) / sd_i  (row i against its own row)"),
              sprintf("lexA row: %d values, normal approximation", ex$G - 1)),
            xr + 0.55, Y(r1 + 3.0), 0.16, 7, ink[["secondary"]])
  # ---- panel 3: clip ----
  panel_head(3, "Clip at 0", xl, Y(r2), "negative z (below the row mean) carries no evidence of dependence: set to 0")
  g3 <- mat_draw(xl + 0.42, Y(r2 + 0.82), CELL, ex$Zraw_sub, nm, ramp = RAMP$orange, vmax = 15, gamma = 0.9, d = 1, neg_grey = TRUE, diag_grey = TRUE)
  cell_box(g3, 1, J, "black", 1.6, 1); cell_box(g3, J, 1, "black", 1.6, 2)
  cell_box(g3, 1, BG, ink[["secondary"]], 1.6, 3)
  txt("row i = the gene whose background is used", xl + 0.42, Y(r2 + 0.82 + n * CELL + 0.14), 7, ink[["secondary"]], just = "left")
  txt(sprintf("grey = z < 0, clipped to 0 (%d of %d off-diagonal cells)", sum(ex$Zraw_sub[row(ex$Zraw_sub) != col(ex$Zraw_sub)] < 0), n * (n - 1)), xl + 0.42, Y(r2 + 0.82 + n * CELL + 0.3), 7, ink[["secondary"]], just = "left")
  # ---- panel 4: Stouffer ----
  panel_head(4, "Symmetric Stouffer score", xr, Y(r2), "S_ij = (max(z_ij, 0) + max(z_ji, 0)) / sqrt(2)")
  g4 <- mat_draw(xr + 0.42, Y(r2 + 0.82), CELL, ex$S_sub, nm, ramp = RAMP$orange, vmax = 18, gamma = 0.9, d = 1)
  cell_box(g4, 1, J, "black", 1.6, 1); cell_box(g4, J, 1, "black", 1.6, 2); cell_box(g4, 1, BG, ink[["secondary"]], 1.6, 3)
  zb <- ex$Zraw_sub[1, BG]; zb2 <- ex$Zraw_sub[BG, 1]
  lines_txt(c(sprintf("lexA-%s (solid): (%.1f + %.1f) / sqrt(2) = %.1f", nm[J], zi, zj, sij),
              sprintf("lexA-%s (dotted): (%.1f + %.1f) / sqrt(2) = %.1f", nm[BG], max(zb, 0), max(zb2, 0), ex$S_sub[1, BG]),
              "symmetric: S_ij = S_ji; this is the CLR score"), xr + 0.42, Y(r2 + 0.82 + n * CELL + 0.14), 0.16, 7, ink[["secondary"]])
  # arrows between panels
  harrow(xl + 3.45, Y(r1 + 1.6), xr - 0.02, Y(r1 + 1.6), ink[["muted"]], 1.2)
  harrow(xl + 3.45, Y(r2 + 1.6), xr - 0.02, Y(r2 + 1.6), ink[["muted"]], 1.2)
}
mf_save("clr", H, draw)
cat(sprintf("pair lexA-%s: MI %.4f z_ij %.3f z_ji %.3f S %.3f; background gene %s z_ij %.3f z_ji %.3f S %.3f\n", nm[J], ex$D_sub[1, J], zi, zj, sij, nm[BG], ex$Zraw_sub[1, BG], ex$Zraw_sub[BG, 1], ex$S_sub[1, BG]))
