# method_baselines.png: the four comparison methods, one card each. Neutral wording, no results. Schematics are drawings (no numbers); the |r| card shows the real M3D
# profiles of lexA and recA (analysis/method_figs/01_compute_m3d_variants.R). Run:  nice -n 19 Rscript analysis/method_figs/17_method_baselines.R
source("analysis/method_figs/mf_lib.R")
v <- readRDS("analysis/method_figs/cache/m3d_variants.rds"); H <- 4.4
GREY <- ink[["muted"]]; DK <- ink[["secondary"]]
node <- function(x, y, r = 0.045, fill = GREY) grid.circle(x, y, r, default.units = "native", gp = gpar(fill = fill, col = NA))
tree <- function(x, y, s = 0.11, col = DK) { seg(x, y, x - s, y - s, col, 0.9); seg(x, y, x + s, y - s, col, 0.9); seg(x - s, y - s, x - 1.5 * s, y - 2 * s, col, 0.9); seg(x - s, y - s, x - 0.5 * s, y - 2 * s, col, 0.9)
  node(x, y, 0.03, col); node(x - s, y - s, 0.03, col); node(x + s, y - s, 0.03, col); node(x - 1.5 * s, y - 2 * s, 0.025, GREY); node(x - 0.5 * s, y - 2 * s, 0.025, GREY) }
draw <- function() {
  txt("Comparison methods: what each one scores", 0.1, Y(0.16), 8.5, face = "bold", just = "left")
  pill("schematics are drawings, not data; the |r| panel is real", MF_W - 0.1, Y(0.16), 0.2, 7, "#f1f0ed", ink[["secondary"]], just = "right")
  cw <- 1.55; xs <- 0.12 + (0:3) * (cw + 0.19); top <- 0.45; ch <- 3.6
  cards <- list(
    list(t = "GENIE3", s = "tree ensemble per target", b = "For each target gene, a random-forest-type ensemble predicts its expression from all other genes. A candidate regulator's score for that target is its importance in the ensemble. One model per gene; scores are directed."),
    list(t = "GRNBoost2", s = "boosted trees per target", b = "The same per-target scheme with gradient-boosted trees and early stopping. A candidate regulator's score for that target is its importance in the boosted ensemble. One model per gene; scores are directed."),
    list(t = "PIDC", s = "partial information decomposition", b = "Expression is discretised and the mutual information of each pair is decomposed against every third gene. A pair's score aggregates the share of its information that is unique to the pair. Scores are undirected."),
    list(t = "|r|", s = "absolute Pearson correlation", b = "The absolute Pearson correlation between the two genes' profiles across samples. One value per pair; undirected; no normalisation against a background."))
  for (k in 1:4) { x <- xs[k]; cd <- cards[[k]]
    rect(x, Y(top + ch), cw, ch, fill = "#f7f6f3"); rect(x, Y(top) - 0.035, cw, 0.035, fill = GREY)
    txt(cd$t, x + 0.1, Y(top + 0.2), 9, face = "bold", just = "left"); txt(cd$s, x + 0.1, Y(top + 0.4), 7, ink[["secondary"]], just = "left", face = "italic")
    cy <- Y(top + 1.15)   # schematic centre height
    if (k == 1) { for (q in 1:5) { node(x + 0.2, cy + (q - 3) * 0.17); seg(x + 0.25, cy + (q - 3) * 0.17, x + 0.52, cy + ((q - 3) * 0.17) * 0.4, GREY, 0.5) }
      for (q in -1:1) tree(x + 0.72, cy + q * 0.3 + 0.1); harrow(x + 0.95, cy, x + 1.2, cy, DK, 0.9); node(x + 1.35, cy, 0.075, DK)
      txt("predictors", x + 0.1, cy - 0.58, 7, DK, just = "left"); txt("target", x + 1.35, cy - 0.2, 7, DK) }
    if (k == 2) { for (q in 0:2) { tree(x + 0.3 + q * 0.4, cy + 0.2, 0.09); if (q < 2) txt("+", x + 0.5 + q * 0.4, cy + 0.12, 8, DK) }
      harrow(x + 0.3, cy - 0.28, x + 1.2, cy - 0.28, DK, 0.9); txt("trees added in sequence,", x + cw / 2, cy - 0.5, 7, DK); txt("each fitted to the previous error", x + cw / 2, cy - 0.64, 7, DK) }
    if (k == 3) { px <- c(x + 0.45, x + 1.15, x + 0.8); py <- c(cy - 0.2, cy - 0.2, cy + 0.38)
      seg(px[1], py[1], px[2], py[2], DK, 1.6); seg(px[1], py[1], px[3], py[3], GREY, 0.8); seg(px[2], py[2], px[3], py[3], GREY, 0.8)
      for (q in 1:3) node(px[q], py[q], 0.06, if (q < 3) DK else GREY); txt("X", px[1] - 0.14, py[1], 7, DK); txt("Y", px[2] + 0.14, py[2], 7, DK); txt("Z", px[3], py[3] + 0.15, 7, DK)
      txt("X-Y information", x + cw / 2, cy - 0.45, 7, DK); txt("not carried by Z", x + cw / 2, cy - 0.59, 7, DK) }
    if (k == 4) { bx <- x + 0.3; bw <- 1.1; bh <- 0.95; by <- cy - 0.5
      rect(bx, by, bw, bh, fill = NA, col = ink[["grid"]], lwd = 0.8); xr <- range(v$x); yr <- range(v$y)
      px <- bx + (v$x - xr[1]) / diff(xr) * bw; py <- by + (v$y - yr[1]) / diff(yr) * bh
      grid.points(px, py, default.units = "native", pch = 19, size = unit(1.6, "pt"), gp = gpar(col = "#52514e", alpha = 0.5))
      txt(sprintf("%s", v$pair[1]), bx + bw / 2, by - 0.1, 7, DK, face = "italic"); txt(v$pair[2], bx - 0.07, by + bh / 2, 7, DK, rot = 90, face = "italic")
      txt(sprintf("r = %.2f, %d profiles", v$r, v$n), x + cw / 2, cy - 0.72, 7, DK) }
    lines_txt(strwrap(cd$b, 28), x + 0.1, Y(top + 2.15), 0.15, 7, ink[["secondary"]]) }
  txt("Each method turns the same genes-by-profiles matrix into one score per gene pair (or per directed pair), which is then ranked.", 0.1, Y(4.2), 7, ink[["secondary"]], just = "left")
}
mf_save("baselines", H, draw)
