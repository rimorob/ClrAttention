# method_tcaclr2.png: data-derived walk depth. Observed and permuted spectra of M3D (00_compute_m3d.R: 10 per-gene permutations, seeds as in the stored edge files);
# the five-dataset table is read from the stored edge files. Run:  nice -n 19 Rscript analysis/method_figs/13_method_tcaclr2.R
source("analysis/method_figs/mf_lib.R")
ex <- readRDS("analysis/method_figs/cache/m3d_example.rds"); et <- ex$edge_tab
H <- 5.05
box <- function(num, title, lines, x, y, w = 1.28, h = 0.95, col = ink[["primary"]]) {
  rect(x, y - h, w, h, fill = "#f7f6f3", col = NA); rect(x, y - h, 0.035, h, fill = MGT)
  chip(num, x + 0.2, y - 0.17, 0.095); txt(title, x + 0.34, y - 0.17, 7, face = "bold", just = "left")
  lines_txt(lines, x + 0.1, y - 0.42, 0.135, 7, ink[["secondary"]]) }
draw <- function() {
  txt("tcaCLR2: the walk depth from the data, no truth read", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  badge_example(MF_W - 0.1, Y(0.16), sprintf("Example: E. coli M3D (%d genes x %d profiles)", ex$G, ex$N))
  y0 <- Y(0.4); xs <- 0.1 + (0:4) * 1.38
  box(1, "Permute", c("each gene's profile", "shuffled on its own,", "10 times"), xs[1], y0)
  box(2, "Refit", c("CLR scores and the", "aCLR operator (keff 10,", "cap 50, lazy 0.5)"), xs[2], y0)
  box(3, "Eigenvalues", c("largest non-trivial", "modulus of the", "symmetrised operator"), xs[3], y0)
  box(4, "Edge", c("mean over the 10", "permutations", "(noise level of the walk)"), xs[4], y0)
  box(5, "Depth", c("t2 = max(1, round(", "ln 0.05 / ln edge))", "the 0.05 convention"), xs[5], y0)
  for (k in 1:4) harrow(xs[k] + 1.285, y0 - 0.48, xs[k + 1] - 0.005, y0 - 0.48, ink[["muted"]], 1)
  # ---- spectrum ----
  yt <- 1.65; hh <- 2.0
  txt("Eigenvalue spectrum of the symmetrised operator", 0.1, Y(yt), 8, face = "bold", just = "left")
  K <- 40; yr <- c(floor(min(ex$ev_obs, ex$perm_ev) * 20) / 20, 1)
  a <- axes_box(0.62, 3.1, Y(yt + 0.35 + hh), hh, c(0, K + 1), yr, c(1, 10, 20, 30, 40), seq(yr[1], 1, by = 0.05), yfmt = function(v) sub("^0", "", sprintf("%.2f", v)), xlab = "eigenvalue index (largest first)", ylab = "eigenvalue")
  for (p in 1:ncol(ex$perm_ev)) polyline(a$X(1:K), a$Y(ex$perm_ev[, p]), "#b9b8b3", 0.9)
  polyline(a$X(1:K), a$Y(ex$ev_obs), BLT, 1.6)
  seg(a$X(2), a$Y(ex$edge), a$X(K + 1), a$Y(ex$edge), MGT, 1.2, 2)
  dot(a$X(2), a$Y(ex$edge), MGT, 5)
  txt(sprintf("edge %.3f", ex$edge), a$X(K + 0.5), a$Y(ex$edge) + 0.1, 7, MGT, just = "right")
  txt("observed data", a$X(14), a$Y(ex$ev_obs[14]) + 0.14, 7, BLT, just = "left")
  txt("10 permutations", a$X(14), a$Y(mean(ex$perm_ev[14, ])) - 0.16, 7, ink[["secondary"]], just = "left")
  txt("trivial eigenvalue 1", a$X(2.5), a$Y(1) - 0.03, 7, ink[["secondary"]], just = "left")
  # ---- decay ----
  txt("Noise decay and the depth", 3.45, Y(yt), 8, face = "bold", just = "left")
  b <- axes_box(3.95, 5.55, Y(yt + 0.35 + hh), hh, c(0, 45), c(0, 1), c(0, 10, 20, 30, 40), c(0, 0.05, 0.25, 0.5, 0.75, 1), yfmt = function(v) sub("^0", "", sprintf("%g", v)), xlab = "walk depth t", ylab = "edge^t")
  tt <- seq(0, 45, length.out = 200); polyline(b$X(tt), b$Y(ex$edge^tt), MGT, 1.6)
  seg(b$X(0), b$Y(0.05), b$X(45), b$Y(0.05), ink[["secondary"]], 0.9, 3)
  seg(b$X(ex$t2), b$Y(0.05), b$X(ex$t2), b$Y(0), ink[["primary"]], 1)
  dot(b$X(ex$t2), b$Y(ex$edge^ex$t2), "black", 5)
  lines_txt(c(sprintf("t2 = %d", ex$t2), sprintf("(t* = %.2f)", ex$tstar)), b$X(15), b$Y(0.42), 0.15, 7, ink[["primary"]]); seg(b$X(ex$t2) - 0.02, b$Y(0.05) + 0.06, b$X(20), b$Y(0.33), ink[["muted"]], 0.7)
  txt("0.05", b$X(45) , b$Y(0.05) + 0.1, 7, ink[["secondary"]], just = "right")
  # ---- table ----
  tx <- 5.85; txt("Five datasets", tx, Y(yt), 8, face = "bold", just = "left")
  nmap <- c(m3d = "M3D", ecomac = "EcoMAC", k562 = "K562", hepg2 = "HepG2", kemmeren = "Kemmeren"); rows <- et[et$dataset %in% names(nmap), ]; rows <- rows[match(names(nmap), rows$dataset), ]
  txt("dataset", tx, Y(yt + 0.38), 7, ink[["secondary"]], just = "left"); txt("edge", tx + 0.78, Y(yt + 0.38), 7, ink[["secondary"]], just = "right"); txt("t2", tx + 1.05, Y(yt + 0.38), 7, ink[["secondary"]], just = "right")
  seg(tx, Y(yt + 0.46), tx + 1.1, Y(yt + 0.46), ink[["muted"]], 0.7)
  for (k in seq_len(nrow(rows))) { yy <- Y(yt + 0.62 + (k - 1) * 0.2)
    txt(nmap[[rows$dataset[k]]], tx, yy, 7, just = "left"); txt(sprintf("%.3f", rows$edge[k]), tx + 0.78, yy, 7, just = "right"); txt(sprintf("%d", rows$t2[k]), tx + 1.05, yy, 7, face = "bold", just = "right") }
  lines_txt(c("edge: mean over 10", "permutations; values from", "the stored runs. The M3D", "panels redraw the same", "permutations and agree."), tx, Y(yt + 1.75), 0.14, 7, ink[["secondary"]])
  seg(0.1, Y(4.5), 6.9, Y(4.5), ink[["grid"]], 0.8)
  txt("Not adapted: keff 10, cap 50, lazy weight 0.5. Not read: any regulon, complex or other truth, any held-out score. Applies to tcaCLR (second CLR on log P^t).", 0.1, Y(4.72), 7, ink[["secondary"]], just = "left")
}
mf_save("tcaclr2", H, draw)
cat(sprintf("M3D edge %.6f t2 %d tstar %.3f; table: %s\n", ex$edge, ex$t2, ex$tstar, paste(et$dataset, et$t2, collapse = " ")))
