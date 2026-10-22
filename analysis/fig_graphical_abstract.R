# Graphical abstract (16:9), drawn entirely in code with grid (no generative imagery).
# Brief: docs/whitepapers/graphical_abstract_prompt.md, revised after PI review:
# published CLR (Faith et al. 2007) is ONLY the per-gene context normalization;
# the softmax and the multi-step walk are this paper's additions.
# Palette: analysis/plot_theme.R (blue = MI, orange = CLR normalization, violet = softmax,
# magenta = depth/walk, aqua = the full method of this work, grey = competitors,
# red = loss marker only). Network drawings and the heat map are schematic (seeded
# layouts, not data). Bars use real numbers from docs/whitepapers/sota_table.md
# (raw pooled co-membership AUPR / base rate; this work at t = 9).
#
# Output: docs/whitepapers/figures/graphical_abstract.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")

out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"   # Arial-metric sans

# ---- numbers (docs/whitepapers/sota_table.md, raw scores) -------------------
this_lab <- "CLR + softmax + walk (this work)"
strip_a <- data.frame(   # M3D, strong-evidence (SC) RegulonDB regulons
  label = c(this_lab, "GENIE3", "PIDC", "CLR"),
  value = c(1.68, 1.51, 1.45, 1.44),
  col = c(pal[["aqua"]], ink[["muted"]], ink[["muted"]], pal[["orange"]]),
  stringsAsFactors = FALSE)
strip_b <- data.frame(   # K562, CORUM complexes
  label = c(this_lab, "CLR", "GENIE3", "|r| (correlation)"),
  value = c(3.51, 6.39, 6.32, 7.00),
  col = c(pal[["aqua"]], pal[["orange"]], ink[["muted"]], ink[["muted"]]),
  stringsAsFactors = FALSE)

# ---- helpers ----------------------------------------------------------------
txt <- function(label, x, y, size = 13, col = ink[["primary"]], face = "plain",
                just = "centre", rot = 0, lineheight = 1.05)
  grid.text(label, x, y, default.units = "native", just = just, rot = rot,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face,
                      lineheight = lineheight))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, alpha = 1, lty = 1, arrow = NULL)
  grid.segments(x0, y0, x1, y1, default.units = "native", arrow = arrow,
                gp = gpar(col = col, lwd = lwd, alpha = alpha, lty = lty, lineend = "round",
                          fill = col))
node <- function(x, y, r = 0.075, fill = "white", col = ink[["secondary"]], lwd = 1.2, alpha = 1)
  grid.circle(x, y, r, default.units = "native",
              gp = gpar(fill = fill, col = col, lwd = lwd, alpha = alpha))
query <- function(x, y, r = 0.15) {      # query gene i: neutral black
  node(x, y, r, fill = ink[["primary"]], col = ink[["primary"]])
  txt("i", x, y, 12, "white", "bold.italic")
}
ring <- function(n, cx, cy, r, phase = 0)
  data.frame(x = cx + r * cos(phase + 2 * pi * (seq_len(n) - 1) / n),
             y = cy + r * sin(phase + 2 * pi * (seq_len(n) - 1) / n))
trim <- function(x0, y0, x1, y1, r0, r1) {
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2); ux <- (x1 - x0) / d; uy <- (y1 - y0) / d
  list(x0 = x0 + ux * r0, y0 = y0 + uy * r0, x1 = x1 - ux * r1, y1 = y1 - uy * r1)
}
card <- function(x, y, w, h, num, title, size = 12.5) {
  grid.roundrect(x, y, w, h, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.1, "inches"), gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.4))
  grid.circle(x + 0.27, y + h - 0.29, 0.16, default.units = "native",
              gp = gpar(fill = ink[["primary"]], col = NA))
  txt(num, x + 0.27, y + h - 0.295, 11.5, "white", "bold")
  grid.text(title, x + 0.53, y + h - 0.17, default.units = "native", just = c("left", "top"),
            gp = gpar(fontfamily = FONT, fontsize = size, fontface = "bold",
                      col = ink[["primary"]], lineheight = 1.0))
}
chevron <- function(x, y, col = ink[["muted"]])
  grid.polygon(c(x - 0.055, x + 0.075, x - 0.055), c(y + 0.13, y, y - 0.13),
               default.units = "native", gp = gpar(fill = col, col = NA))
pill <- function(label, x, y, w, h = 0.36, fill = "white", col = ink[["primary"]], size = 12,
                 face = "bold", tcol = ink[["primary"]], lwd = 1.3) {
  grid.roundrect(x, y, w, h, default.units = "native", r = unit(h / 2, "inches"),
                 gp = gpar(fill = fill, col = col, lwd = lwd))
  txt(label, x, y, size, tcol, face)
}
band <- function(x0, x1, y, fill, label) {   # labelled band + bracket over a group of panels
  grid.roundrect(x0, y, x1 - x0, 0.46, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = fill, col = NA))
  txt(label, x0 + 0.18, y + 0.23, 14, face = "bold", just = "left")
  by <- y - 0.1
  grid.lines(c(x0 + 0.05, x0 + 0.05, x1 - 0.05, x1 - 0.05), c(by - 0.09, by, by, by - 0.09),
             default.units = "native", gp = gpar(col = ink[["secondary"]], lwd = 1.3))
}

# bar strip with its own axis; returns nothing
bar_strip <- function(d, x0, w, y_top, xmax, ticks, loss = FALSE) {
  lab_x <- x0 + 0.25; bx0 <- x0 + 2.85; bx1 <- x0 + w - 0.6
  sc <- (bx1 - bx0) / xmax
  ry <- y_top - (0:3) * 0.31; bh <- 0.21
  # grid + random line
  for (t in ticks) seg(bx0 + t * sc, ry[1] + 0.17, bx0 + t * sc, ry[4] - 0.17, ink[["grid"]], 0.9)
  for (k in 1:4) {
    txt(d$label[k], lab_x, ry[k], 11, face = if (k == 1) "bold" else "plain", just = "left")
    grid.rect(bx0, ry[k], d$value[k] * sc, bh, default.units = "native",
              just = c("left", "centre"), gp = gpar(fill = d$col[k], col = NA))
    if (loss && k == 1) {
      txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc - 0.05, ry[k], 10.5, "white", "bold",
          just = "right")
    } else {
      txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc + 0.05, ry[k], 10.5,
          face = if (k == 1) "bold" else "plain", just = "left")
    }
  }
  seg(bx0 + sc, ry[1] + 0.26, bx0 + sc, ry[4] - 0.17, ink[["primary"]], 1.1, lty = 2)
  txt("random = 1", bx0 + sc + 0.04, ry[1] + 0.3, 9.5, ink[["secondary"]], just = "left")
  seg(bx0, ry[1] + 0.17, bx0, ry[4] - 0.17, ink[["secondary"]], 1)
  ay <- ry[4] - 0.17
  seg(bx0, ay, bx1, ay, ink[["secondary"]], 1)
  for (t in ticks) { seg(bx0 + t * sc, ay, bx0 + t * sc, ay - 0.05, ink[["secondary"]], 1)
    txt(as.character(t), bx0 + t * sc, ay - 0.16, 9.5, ink[["secondary"]]) }
  txt("precision–recall area as a multiple of random", (bx0 + bx1) / 2, ay - 0.38, 9.5,
      ink[["secondary"]])
  if (loss) {   # red loss marker: from the CLR value back to this work's value
    ly <- ry[1]
    seg(bx0 + d$value[2] * sc, ly, bx0 + d$value[1] * sc + 0.04, ly, pal[["red"]], 2.4,
        arrow = arrow(length = unit(0.08, "inches"), type = "closed"))
    seg(bx0 + d$value[2] * sc, ly - 0.09, bx0 + d$value[2] * sc, ly + 0.09, pal[["red"]], 2.4)
    txt(sprintf("%.1f× lower", d$value[2] / d$value[1]),
        bx0 + d$value[2] * sc + 0.07, ly, 10, face = "bold", just = "left")
  }
}

# ---- drawing ----------------------------------------------------------------
draw_abstract <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))

  # headline: the finding, with its cost
  grid.text(expression(bold("Adding a softmax and a multi-step walk to CLR's normalization recovers large ") *
                         bolditalic("E. coli") * bold(" regulons better")),
            0.3, 8.62, default.units = "native", just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 18, col = ink[["primary"]]))
  txt("(16% over CLR), but loses resolution on tight protein complexes",
      0.3, 8.26, 18, face = "bold", just = "left")

  W <- 2.4; gap <- 0.2; Y <- 1.5; H <- 5.5
  X <- 0.3 + (0:3) * (W + gap)
  xr <- X[4] + W                      # right edge of the method panels (10.5)
  band(X[1], X[2] + W, 7.38, "#ecebe7", "CLR (Faith et al. 2007)")
  band(X[3], xr, 7.38, "#e6e3fa", "This work: steps taken from attention")

  card(X[1], Y, W, H, 1, "Expression\ncompendium")
  card(X[2], Y, W, H, 2, "Mutual information,\nz-scored against\neach gene's own\nbackground")
  card(X[3], Y, W, H, 3, "Softmax with\nper-gene\ntemperature")
  card(X[4], Y, W, H, 4, "Multi-step walk,\nnumber of steps\nchosen without\nlabels")
  for (i in 1:3) chevron(X[i] + W + gap / 2, Y + H / 2)

  # ---------------- panel 1: heat-map motif + organisms ----------------------
  x0 <- X[1]
  set.seed(7)
  ng <- 12; ne <- 14
  mod <- rep(1:3, c(5, 4, 3))
  prof <- matrix(rnorm(3 * ne), 3)
  M <- prof[mod, ] * 0.95 + matrix(rnorm(ng * ne, sd = 0.55), ng)
  M <- (M - min(M)) / diff(range(M))
  ramp <- colorRampPalette(seq_blue[1:7])(101)
  hx0 <- x0 + 0.48; hx1 <- x0 + W - 0.2; hy0 <- 3.75; hy1 <- 5.85
  cw <- (hx1 - hx0) / ne; ch <- (hy1 - hy0) / ng
  for (g in 1:ng) for (e in 1:ne)
    grid.rect(hx0 + (e - 1) * cw, hy1 - g * ch, cw * 0.9, ch * 0.88, default.units = "native",
              just = c("left", "bottom"), gp = gpar(fill = ramp[1 + round(100 * M[g, e])], col = NA))
  txt("genes", x0 + 0.3, (hy0 + hy1) / 2, 11.5, ink[["secondary"]], rot = 90)
  txt("experiments", hx0, hy0 - 0.22, 11.5, ink[["secondary"]], just = "left")
  seg(hx0 + 1.0, hy0 - 0.22, hx1, hy0 - 0.22, ink[["muted"]], 1.2,
      arrow = arrow(length = unit(0.06, "inches"), type = "closed"))
  # bacterium (E. coli) and human cell (K562)
  bx <- x0 + 0.62; by <- 2.75
  grid.curve(bx + 0.3, by, bx + 0.62, by - 0.04, default.units = "native", curvature = -0.6,
             square = FALSE, ncp = 8, gp = gpar(col = ink[["secondary"]], lwd = 1.2))
  grid.roundrect(bx, by, 0.66, 0.3, default.units = "native", r = unit(0.15, "inches"),
                 gp = gpar(fill = seq_blue[3], col = ink[["secondary"]], lwd = 1.3))
  txt("E. coli", bx + 0.05, by - 0.45, 12.5, face = "italic")
  txt("M3D, EcoMAC\nmicroarrays", bx + 0.05, by - 0.88, 10.5, ink[["secondary"]])
  hx <- x0 + 1.78; hy <- 2.75
  grid.circle(hx, hy, 0.25, default.units = "native",
              gp = gpar(fill = seq_blue[2], col = ink[["secondary"]], lwd = 1.3))
  grid.circle(hx + 0.04, hy + 0.01, 0.1, default.units = "native",
              gp = gpar(fill = seq_blue[5], col = ink[["secondary"]], lwd = 1))
  txt("human", hx, hy - 0.45, 12.5)
  txt("K562\nPerturb-seq", hx, hy - 0.88, 10.5, ink[["secondary"]])

  # ---------------- panel 2: CLR = MI z-scored per gene ----------------------
  x0 <- X[2]; cx <- x0 + W / 2; cy <- 4.4
  pill("existing method", x0 + 1.02, 5.88, 1.5, 0.3, fill = "#ecebe7", col = NA, size = 10.5,
       face = "plain", tcol = ink[["secondary"]])
  nb <- ring(10, cx, cy, 0.78, phase = 0.25)
  mi <- c(0.95, 0.25, 0.4, 0.15, 0.8, 0.3, 0.55, 0.12, 0.7, 0.2)
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.15, 0.07)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["blue"]], 0.8 + 6 * mi[k], alpha = 0.35 + 0.6 * mi[k])
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.07)
  query(cx, cy)
  txt("line width = MI", cx, 3.36, 10.5, ink[["secondary"]])
  # this gene's own background (orange = CLR normalization)
  gx0 <- x0 + 0.95; gx1 <- x0 + W - 0.15; gy0 <- 2.05; gh <- 0.72
  zz <- seq(-3.2, 3.2, length.out = 120)
  gxs <- gx0 + (zz + 3.2) / 6.4 * (gx1 - gx0)
  grid.polygon(c(gxs, rev(gxs)), c(gy0 + gh * exp(-zz^2 / 2), rep(gy0, length(zz))),
               default.units = "native", gp = gpar(fill = pal[["orange"]], alpha = 0.22, col = NA))
  grid.lines(gxs, gy0 + gh * exp(-zz^2 / 2), default.units = "native",
             gp = gpar(col = pal[["orange"]], lwd = 2.2))
  seg(gx0, gy0, gx1, gy0, ink[["muted"]], 1)
  zt <- c(-1.6, -1.1, -0.7, -0.4, -0.15, 0.1, 0.35, 0.6, 1.0, 1.4, 2.0, 2.75)
  for (z in zt) { xx <- gx0 + (z + 3.2) / 6.4 * (gx1 - gx0)
    seg(xx, gy0 - 0.1, xx, gy0 + 0.02, pal[["blue"]], ifelse(z > 2.5, 2.8, 1.3)) }
  txt("z", gx0 + (2.75 + 3.2) / 6.4 * (gx1 - gx0), gy0 + 0.27, 12.5, face = "bold.italic")
  txt("gene i's\nown\nbackground", x0 + 0.17, gy0 + 0.35, 10.5, ink[["secondary"]], just = "left")

  # ---------------- panel 3: softmax ------------------------------------------
  x0 <- X[3]; cx <- x0 + W / 2; cy <- 4.35
  nb <- ring(20, cx, cy, 0.86, phase = 0.1)
  top <- c(1, 3, 5, 6, 9, 11, 13, 15, 17, 19)   # ~10 partners
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.15, 0.075)
    if (k %in% top) seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 3,
                        arrow = arrow(length = unit(0.065, "inches"), type = "closed"))
    else seg(s$x0, s$y0, s$x1, s$y1, ink[["muted"]], 0.8, alpha = 0.45)
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.06,
                                    col = if (k %in% top) ink[["secondary"]] else ink[["muted"]],
                                    alpha = if (k %in% top) 1 else 0.7)
  query(cx, cy)
  txt("~10 effective partners\nper gene, set in advance", cx, 2.55, 11, ink[["secondary"]])

  # ---------------- panel 4: multi-step walk ----------------------------------
  x0 <- X[4]; hcx <- x0 + W / 2; hcy <- 4.35
  set.seed(11)
  radii <- c(0.3, 0.6, 0.9)
  for (j in 3:1) grid.circle(hcx, hcy, radii[j], default.units = "native",
                             gp = gpar(fill = pal[["magenta"]], alpha = c(0.32, 0.2, 0.11)[j], col = NA))
  for (j in 1:3) grid.circle(hcx, hcy, radii[j], default.units = "native",
                             gp = gpar(fill = NA, col = pal[["magenta"]], lwd = 1.3, lty = 2))
  ang <- runif(28, 0, 2 * pi); rr <- sqrt(runif(28)) * 0.85 + 0.12
  pts <- data.frame(x = c(hcx, hcx + rr * cos(ang)), y = c(hcy, hcy + rr * sin(ang)))
  D <- as.matrix(dist(pts))
  for (a in 1:28) for (b in (a + 1):29) if (D[a, b] < 0.36)
    seg(pts$x[a], pts$y[a], pts$x[b], pts$y[b], ink[["secondary"]], 0.8, alpha = 0.55)
  for (k in 2:29) node(pts$x[k], pts$y[k], 0.04, col = ink[["secondary"]], lwd = 0.9)
  query(hcx, hcy, 0.12)
  # steps scale (log3): t = 1, 3, 9, 27
  sx0 <- x0 + 0.4; sx1 <- x0 + W - 0.4; sy <- 2.6
  tx <- sx0 + (0:3) / 3 * (sx1 - sx0)
  txt("steps t (log scale)", (sx0 + sx1) / 2, sy + 0.3, 10.5, ink[["secondary"]])
  seg(sx0, sy, sx1, sy, pal[["magenta"]], 4)
  for (k in 1:4) {
    seg(tx[k], sy - 0.09, tx[k], sy + 0.09, ink[["secondary"]], 1.2)
    txt(c("1", "3", "9", "27")[k], tx[k], sy - 0.26, 11, ink[["secondary"]],
        face = if (k == 3) "bold" else "plain")
  }
  grid.circle(tx[3], sy, 0.07, default.units = "native", gp = gpar(fill = ink[["primary"]], col = "white"))
  txt("t = 9: chosen on\nheld-out labs", tx[3], 1.92, 10.5, ink[["secondary"]])

  # ---------------- under the bands ---------------------------------------------
  txt("Published CLR stops here: a score for each gene pair", (X[1] + X[2] + W) / 2, 1.17,
      12, ink[["secondary"]], face = "italic")
  pill("No learned parameters", (X[3] + xr) / 2, 1.17, 2.7, 0.38, size = 12.5)
  # transformer analogy strip spanning both bands
  grid.roundrect(X[1], 0.3, xr - X[1], 0.46, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "#f3f3f1", col = ink[["grid"]]))
  lab <- "Transformer analogy:"
  txt(lab, X[1] + 0.18, 0.53, 11.5, face = "bold", just = "left")
  lw <- convertWidth(grobWidth(textGrob(lab, gp = gpar(fontfamily = FONT, fontsize = 11.5,
                                                         fontface = "bold"))), "native", TRUE)
  txt("MI ≈ query·key score;  CLR z-score ≈ scaled logit;  softmax and walk ≈ attention and stacked layers",
      X[1] + 0.18 + lw + 0.1, 0.53, 11.5, ink[["secondary"]], just = "left")

  # ---------------- panel 5: results -------------------------------------------
  x5 <- xr + 0.25; W5 <- 15.7 - x5
  grid.roundrect(x5, 0.3, W5, 7.54, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.1, "inches"), gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.4))
  grid.circle(x5 + 0.27, 7.55, 0.16, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
  txt("5", x5 + 0.27, 7.545, 11.5, "white", "bold")
  txt("What it recovers", x5 + 0.53, 7.55, 14, face = "bold", just = "left")

  # (a) regulons: icon + title + definition + bars + footnote
  ix <- x5 + W5 - 0.55; iy <- 6.92
  tg <- ring(9, ix, iy, 0.3, phase = 0.2)
  for (k in 1:9) seg(ix, iy, tg$x[k], tg$y[k], pal[["aqua"]], 1.2)
  for (k in 1:9) node(tg$x[k], tg$y[k], 0.04, fill = "#d1efe4", col = pal[["aqua"]], lwd = 1)
  node(ix, iy, 0.09, fill = pal[["aqua"]], col = "#0d6e4c", lwd = 1)
  txt("E. coli regulons (experimentally\nsupported, RegulonDB)", x5 + 0.25, 7.0, 12.5,
      face = "bold", just = "left", lineheight = 1.0)
  txt("regulon = genes controlled by one transcription factor (TF)", x5 + 0.25, 6.5, 10,
      ink[["secondary"]], just = "left")
  bar_strip(strip_a, x5, W5, 6.0, 2, c(0, 0.5, 1, 1.5, 2))
  txt("per regulon: tie with CLR and GENIE3;\nahead in 1,000/1,000 resamplings vs CLR",
      x5 + 0.25, 4.1, 10, ink[["secondary"]], face = "italic", just = "left")

  seg(x5 + 0.25, 3.78, x5 + W5 - 0.25, 3.78, ink[["grid"]], 1.2)

  # (b) complexes: icon + title + loss marker + bars
  cxi <- x5 + W5 - 0.55; cyi <- 3.3
  cl <- ring(5, cxi, cyi, 0.22, phase = pi / 2)
  for (a in 1:4) for (b in (a + 1):5) seg(cl$x[a], cl$y[a], cl$x[b], cl$y[b], ink[["muted"]], 1.3)
  for (k in 1:5) node(cl$x[k], cl$y[k], 0.05, fill = "#e9e8e4", col = ink[["muted"]], lwd = 1)
  txt("Protein complexes (K562 CORUM;\nmostly ribosome): this work loses 1.8×", x5 + 0.25, 3.33,
      12.5, face = "bold", just = "left", lineheight = 1.0)
  txt("complex = proteins that bind as one machine", x5 + 0.25, 2.85, 10, ink[["secondary"]],
      just = "left")
  bar_strip(strip_b, x5, W5, 2.36, 8, c(0, 2, 4, 6, 8), loss = TRUE)

  popViewport()
}

# ---- render -------------------------------------------------------------------
png(file.path(out_dir, "graphical_abstract.png"), width = 2400, height = 1350, res = 150,
    type = "cairo", bg = "white")
draw_abstract(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract.pdf"), width = 16, height = 9, bg = "white")
draw_abstract(); invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract.{png,pdf}"), "\n")
