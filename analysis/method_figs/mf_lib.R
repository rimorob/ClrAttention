# Shared drawing helpers for the eight method figures (analysis/method_figs/). Same canvas conventions and colour key as analysis/fig_method.R and
# analysis/plot_theme.R: blue = dependence (MI), orange = CLR normalisation, violet = softmax, magenta = depth / walk, grey = neutral, aqua = reserved for the
# full method of the paper (not used here). Canvas in inches, y measured UP; Y(d) converts a distance from the top. All text >= 7 pt at final size (7.0 in wide).
suppressPackageStartupMessages({ library(grid); library(Matrix) })
source("analysis/plot_theme.R")
FONT <- "Liberation Sans"
MF_W <- 7.0
ORT <- "#b5461a"; VIT <- "#4a3bb0"; MGT <- "#a8345f"; BLT <- "#184f95"   # darker tints of the key colours, for TEXT and thin lines only
ramp_of <- function(hi, mid) colorRampPalette(c("#fcfcfb", mid, hi))(100)
RAMP <- list(blue = ramp_of("#184f95", "#9ec5f4"), orange = ramp_of("#a8390f", "#f6b595"), violet = ramp_of("#3a2aa0", "#c9c2f4"),
             magenta = ramp_of("#a8345f", "#f5bfd3"), grey = colorRampPalette(c("#fcfcfb", "#b9b8b3", "#3d3c3a"))(100))

mf_state <- new.env(); mf_state$H <- 5
Y <- function(d) mf_state$H - d
txt <- function(label, x, y, size = 7, col = ink[["primary"]], face = "plain", just = "centre", rot = 0)
  grid.text(label, x, y, default.units = "native", just = just, rot = rot, gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1, arrow = NULL)
  grid.segments(x0, y0, x1, y1, default.units = "native", arrow = arrow, gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt", fill = col))
rect <- function(x, y, w, h, fill = NA, col = NA, lwd = 1, alpha = 1, lty = 1, just = c("left", "bottom"))
  grid.rect(x, y, w, h, default.units = "native", just = just, gp = gpar(fill = fill, col = col, lwd = lwd, alpha = alpha, lty = lty, linejoin = "mitre"))
harrow <- function(x0, y0, x1, y1, col = ink[["secondary"]], lwd = 1) seg(x0, y0, x1, y1, col, lwd, arrow = arrow(length = unit(0.075, "inches"), type = "closed"))
tw <- function(label, size = 7, face = "plain") convertWidth(grobWidth(textGrob(label, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face))), "native", TRUE)
pill <- function(label, x, y, h = 0.2, size = 7, fill = "#f1f0ed", col = ink[["primary"]], face = "plain", border = NA, just = "centre", padx = 0.2) {
  w <- tw(label, size, face) + padx
  x0 <- switch(just, centre = x, left = x + w / 2, right = x - w / 2)
  grid.roundrect(x0, y, w, h, default.units = "native", r = unit(h / 2, "inches"), gp = gpar(fill = fill, col = border, lwd = 0.8))
  txt(label, x0, y, size, col, face); invisible(w)
}
chip <- function(num, x, y, r = 0.105) { grid.circle(x, y, r, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA)); txt(num, x, y - 0.002, 8, "white", "bold") }
panel_head <- function(num, title, x, y, sub = NULL, col = NULL) {
  chip(num, x + 0.105, y); txt(title, x + 0.27, y, 8, face = "bold", just = "left")
  if (!is.null(sub)) txt(sub, x + 0.27, y - 0.17, 7, ink[["secondary"]], just = "left")
  if (!is.null(col)) rect(x, y - 0.3, 0.035, 0.0, fill = col)
}
badge_exploratory <- function(x, y) pill("EXPLORATORY: not preregistered", x, y, h = 0.2, size = 7, fill = ink[["primary"]], col = "white", face = "bold", just = "right")
badge_example <- function(x, y, label) pill(label, x, y, h = 0.2, size = 7, fill = "#f1f0ed", col = ink[["secondary"]], just = "right")

# number format: drop the leading zero (.07), keep sign
fnum <- function(v, d = 2) { s <- formatC(v, format = "f", digits = d); s <- sub("^(-?)0\\.", "\\1.", s); s[as.numeric(s) == 0] <- "0"; s }
# matrix with numbers; top-left (x0, ytop); returns geometry. neg: cells with value < 0 are filled grey (clipped / not used). labels on top (rotated) and left.
mat_draw <- function(x0, ytop, cell, M, rn, cn = rn, ramp = RAMP$grey, vmax = NULL, gamma = 0.8, numbers = TRUE, d = 2, fmt = NULL, neg_grey = FALSE,
                     diag_grey = TRUE, top_labels = TRUE, left_labels = TRUE, lab_size = 7, num_size = 7, thresh = 0.62, hatch_neg = FALSE) {
  n <- nrow(M); m <- ncol(M); off <- M; if (diag_grey && n == m) diag(off) <- NA
  if (is.null(vmax)) vmax <- max(off[is.finite(off)], na.rm = TRUE)
  for (r in seq_len(n)) for (cc in seq_len(m)) {
    v <- M[r, cc]; isd <- diag_grey && n == m && r == cc
    v01 <- min(max(v / vmax, 0), 1)^gamma
    fl <- if (isd) "#ececea" else if (neg_grey && v < 0) "#dcdbd7" else ramp[max(1, ceiling(v01 * 100))]
    rect(x0 + (cc - 1) * cell, ytop - r * cell, cell, cell, fill = fl, col = "white", lwd = 0.5)
    if (numbers && !isd) { lab <- if (is.null(fmt)) fnum(v, d) else fmt(v)
      tc <- if (neg_grey && v < 0) ink[["muted"]] else if (v01 > thresh) "white" else ink[["primary"]]
      txt(lab, x0 + (cc - 0.5) * cell, ytop - (r - 0.5) * cell, num_size, tc) }
  }
  if (left_labels) for (r in seq_len(n)) txt(rn[r], x0 - 0.04, ytop - (r - 0.5) * cell, lab_size, ink[["secondary"]], just = "right", face = "italic")
  if (top_labels) for (cc in seq_len(m)) txt(cn[cc], x0 + (cc - 0.5) * cell + 0.0, ytop + 0.04, lab_size, ink[["secondary"]], just = "left", rot = 90, face = "italic")
  invisible(list(x0 = x0, ytop = ytop, cell = cell))
}
cell_box <- function(g, r, cc, col = "black", lwd = 1.6, lty = 1) rect(g$x0 + (cc - 1) * g$cell, g$ytop - r * g$cell, g$cell, g$cell, col = col, lwd = lwd, lty = lty)
row_band <- function(g, r, n, fill = ink[["secondary"]], alpha = 0.18) rect(g$x0 - 0.02, g$ytop - r * g$cell, n * g$cell + 0.04, g$cell, fill = fill, alpha = alpha)

# histogram panel (counts or densities). x range xr, bars from breaks; returns mapping functions
hist_panel <- function(x0, x1, ybase, h, vals, xr, breaks, fill = "#9ec5f4", ticks = NULL, tick_fmt = function(v) sprintf("%g", v), xlab = NULL, ymax = NULL, border = NA) {
  hh <- hist(pmin(pmax(vals, xr[1]), xr[2]), breaks = breaks, plot = FALSE); cnt <- hh$counts; if (is.null(ymax)) ymax <- max(cnt)
  X_ <- function(v) x0 + (v - xr[1]) / (xr[2] - xr[1]) * (x1 - x0); Y_ <- function(c) ybase + c / ymax * h
  for (k in seq_along(cnt)) rect(X_(hh$breaks[k]), ybase, X_(hh$breaks[k + 1]) - X_(hh$breaks[k]), Y_(cnt[k]) - ybase, fill = fill, col = "white", lwd = 0.3)
  seg(x0, ybase, x1, ybase, ink[["secondary"]], 0.8)
  if (!is.null(ticks)) for (tk in ticks) { seg(X_(tk), ybase, X_(tk), ybase - 0.04, ink[["secondary"]], 0.7); txt(tick_fmt(tk), X_(tk), ybase - 0.13, 7, ink[["secondary"]]) }
  if (!is.null(xlab)) txt(xlab, (x0 + x1) / 2, ybase - 0.3, 7, ink[["secondary"]])
  list(X = X_, Y = Y_, counts = cnt, breaks = hh$breaks, ymax = ymax)
}
# generic scatter/line axes: a small plot box with ticks; returns X_, Y_
axes_box <- function(x0, x1, ybase, h, xr, yr, xticks, yticks, xfmt = function(v) sprintf("%g", v), yfmt = function(v) sprintf("%g", v), xlab = NULL, ylab = NULL, grid = TRUE) {
  X_ <- function(v) x0 + (v - xr[1]) / (xr[2] - xr[1]) * (x1 - x0); Y_ <- function(v) ybase + (v - yr[1]) / (yr[2] - yr[1]) * h
  if (grid) for (tk in yticks) seg(x0, Y_(tk), x1, Y_(tk), ink[["grid"]], 0.5)
  seg(x0, ybase, x1, ybase, ink[["secondary"]], 0.8); seg(x0, ybase, x0, ybase + h, ink[["secondary"]], 0.8)
  for (tk in xticks) { seg(X_(tk), ybase, X_(tk), ybase - 0.04, ink[["secondary"]], 0.7); txt(xfmt(tk), X_(tk), ybase - 0.13, 7, ink[["secondary"]]) }
  for (tk in yticks) { seg(x0 - 0.04, Y_(tk), x0, Y_(tk), ink[["secondary"]], 0.7); txt(yfmt(tk), x0 - 0.07, Y_(tk), 7, ink[["secondary"]], just = "right") }
  if (!is.null(xlab)) txt(xlab, (x0 + x1) / 2, ybase - 0.3, 7, ink[["secondary"]])
  if (!is.null(ylab)) txt(ylab, x0 - 0.42, ybase + h / 2, 7, ink[["secondary"]], rot = 90)
  list(X = X_, Y = Y_)
}
polyline <- function(x, y, col, lwd = 1, lty = 1) grid.lines(x, y, default.units = "native", gp = gpar(col = col, lwd = lwd, lty = lty))
dot <- function(x, y, col, size = 5, pch = 19) grid.points(x, y, default.units = "native", pch = pch, size = unit(size, "pt"), gp = gpar(col = col, fill = col))
# wrapped text lines (explicit vector of lines)
lines_txt <- function(v, x, y, dy = 0.15, size = 7, col = ink[["secondary"]], just = "left", face = "plain") for (q in seq_along(v)) txt(v[q], x, y - (q - 1) * dy, size, col, face, just)

mf_save <- function(name, H, draw, out_dir = "docs/checkpoint_collins_src/figures") {
  mf_state$H <- H; f <- file.path(out_dir, paste0("method_", name, ".png"))
  png(f, width = MF_W * 300, height = H * 300, res = 300, type = "cairo", bg = "white")
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA)); pushViewport(viewport(xscale = c(0, MF_W), yscale = c(0, H), clip = "on")); draw(); popViewport()
  invisible(dev.off()); cat("wrote", f, sprintf("(%.1f x %.1f in)\n", MF_W, H))
}
# neighbourhood names with the roles used in all figures
nb_roles <- function(ex) c("query gene", rep("top CLR partner", 5), rep("found by the walk", 3))
