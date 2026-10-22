# Graphical abstract (16:9), drawn entirely in code with grid (no generative imagery).
# Brief: docs/whitepapers/graphical_abstract_prompt.md. Palette: analysis/plot_theme.R.
# Quantitative elements use real numbers from docs/whitepapers/sota_table.md
# (pooled co-membership AUPR / base rate; CLR attention at t = 9). Network drawings
# and the heat-map motif are schematic (seeded layouts, not data).
#
# Output: docs/whitepapers/figures/graphical_abstract.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")

out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"   # Arial-metric sans; falls back to the device sans if absent

# ---- numbers (docs/whitepapers/sota_table.md, raw scores) -------------------
bars <- data.frame(
  label = c("CLR attention", "GENIE3", "PIDC", "CLR"),
  value = c(1.68, 1.51, 1.45, 1.44),
  col   = c(pal[["aqua"]], ink[["muted"]], ink[["muted"]], pal[["orange"]]),
  stringsAsFactors = FALSE)
corum <- c(attention = 3.51, clr = 6.39)   # K562 CORUM, raw

# ---- helpers ----------------------------------------------------------------
gp_txt <- function(size, col = ink[["primary"]], face = "plain", ...)
  gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face, ...)
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
ring <- function(n, cx, cy, r, phase = 0)
  data.frame(x = cx + r * cos(phase + 2 * pi * (seq_len(n) - 1) / n),
             y = cy + r * sin(phase + 2 * pi * (seq_len(n) - 1) / n))
# shorten a segment so it stops at node boundaries
trim <- function(x0, y0, x1, y1, r0, r1) {
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2); ux <- (x1 - x0) / d; uy <- (y1 - y0) / d
  list(x0 = x0 + ux * r0, y0 = y0 + uy * r0, x1 = x1 - ux * r1, y1 = y1 - uy * r1)
}
card <- function(x, y, w, h, num, title) {
  grid.roundrect(x, y, w, h, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.12, "inches"),
                 gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.4))
  grid.circle(x + 0.32, y + h - 0.34, 0.19, default.units = "native",
              gp = gpar(fill = ink[["primary"]], col = NA))
  txt(num, x + 0.32, y + h - 0.345, 13, "white", "bold")
  txt(title, x + 0.62, y + h - 0.34, 15.5, face = "bold", just = "left", lineheight = 0.95)
}
chevron <- function(x, y)
  grid.polygon(c(x - 0.07, x + 0.09, x - 0.07), c(y + 0.16, y, y - 0.16),
               default.units = "native", gp = gpar(fill = ink[["muted"]], col = NA))
check_mark <- function(x, y, s = 0.16, col = pal[["aqua"]]) {
  grid.circle(x, y, s * 1.35, default.units = "native", gp = gpar(fill = col, col = NA))
  grid.lines(c(x - s * 0.6, x - s * 0.12, x + s * 0.65), c(y + s * 0.02, y - s * 0.5, y + s * 0.55),
             default.units = "native",
             gp = gpar(col = "white", lwd = 3.2, lineend = "round", linejoin = "round"))
}
neutral_mark <- function(x, y, s = 0.16) {
  grid.circle(x, y, s * 1.35, default.units = "native", gp = gpar(fill = ink[["muted"]], col = NA))
  seg(x - s * 0.6, y, x + s * 0.6, y, "white", 3.2)
}

# ---- drawing ----------------------------------------------------------------
draw_abstract <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))

  # title (the message in one line) + subline
  txt("CLR (2007) is context-normalized attention; completing it, with no learning, recovers large regulatory programs",
      0.32, 8.56, 19.5, face = "bold", just = "left")
  txt("The gain requires context normalization; depth trades pair-level resolution for module-level recovery.",
      0.32, 8.14, 14, ink[["secondary"]], just = "left")

  W <- 3.72; H <- 7.45; Y <- 0.3; gap <- 0.21
  X <- 0.3 + (0:3) * (W + gap)
  titles <- c("Expression compendium", "CLR = context-normalized\nattention",
              "Completing the attention", "What it recovers")
  for (i in 1:4) card(X[i], Y, W, H, i, titles[i])
  for (i in 1:3) chevron(X[i] + W + gap / 2, Y + H / 2)

  # ---------------- panel 1: heat-map motif + organisms ----------------------
  x0 <- X[1]
  set.seed(7)
  ng <- 16; ne <- 20
  mod <- rep(1:3, c(6, 5, 5))
  prof <- matrix(rnorm(3 * ne), 3)
  M <- prof[mod, ] * 0.95 + matrix(rnorm(ng * ne, sd = 0.55), ng)
  M <- (M - min(M)) / diff(range(M))
  ramp <- colorRampPalette(seq_blue)(101)
  hx0 <- x0 + 0.62; hx1 <- x0 + W - 0.28; hy0 <- 3.55; hy1 <- 6.95
  cw <- (hx1 - hx0) / ne; ch <- (hy1 - hy0) / ng
  for (g in 1:ng) for (e in 1:ne)
    grid.rect(hx0 + (e - 1) * cw, hy1 - g * ch, cw * 0.92, ch * 0.88, default.units = "native",
              just = c("left", "bottom"),
              gp = gpar(fill = ramp[1 + round(100 * M[g, e])], col = NA))
  txt("genes", x0 + 0.38, (hy0 + hy1) / 2, 13, ink[["secondary"]], rot = 90)
  txt("experiments", (hx0 + hx1) / 2, hy0 - 0.25, 13, ink[["secondary"]])
  seg(hx0 + 1.95, hy0 - 0.25, hx1, hy0 - 0.25, ink[["muted"]], 1.3,
      arrow = arrow(length = unit(0.07, "inches"), type = "closed"))

  # bacterium (E. coli)
  bx <- x0 + 0.98; by <- 2.25
  grid.curve(bx + 0.42, by, bx + 0.95, by - 0.05, default.units = "native", curvature = -0.6,
             square = FALSE, ncp = 8, gp = gpar(col = ink[["secondary"]], lwd = 1.3))
  grid.roundrect(bx, by, 0.95, 0.42, default.units = "native", r = unit(0.21, "inches"),
                 gp = gpar(fill = seq_blue[3], col = ink[["secondary"]], lwd = 1.4))
  txt("E. coli", bx + 0.08, by - 0.62, 14, face = "italic")
  txt("M3D · EcoMAC\nmicroarrays", bx + 0.08, by - 1.15, 12, ink[["secondary"]])
  # human cell (K562)
  hx <- x0 + 2.75; hy <- 2.25
  grid.circle(hx, hy, 0.36, default.units = "native",
              gp = gpar(fill = seq_blue[2], col = ink[["secondary"]], lwd = 1.4))
  grid.circle(hx + 0.05, hy + 0.02, 0.15, default.units = "native",
              gp = gpar(fill = seq_blue[5], col = ink[["secondary"]], lwd = 1.1))
  txt("human", hx, hy - 0.62, 14)
  txt("K562 Perturb-seq\npseudobulk", hx, hy - 1.15, 12, ink[["secondary"]])

  # ---------------- panel 2: CLR as context-normalized attention -------------
  x0 <- X[2]; cx <- x0 + W / 2; cy <- 5.45
  set.seed(3)
  nb <- ring(12, cx, cy, 1.18, phase = 0.2)
  mi <- c(0.95, 0.25, 0.4, 0.15, 0.8, 0.3, 0.2, 0.55, 0.12, 0.35, 0.7, 0.18)
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.2, 0.09)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["blue"]], 0.8 + 7 * mi[k], alpha = 0.35 + 0.6 * mi[k])
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.09)
  node(cx, cy, 0.2, fill = "white", col = pal[["orange"]], lwd = 3.2)
  txt("i", cx, cy, 14, face = "bold.italic")
  txt("edge width = mutual information", cx, 4.02, 12, ink[["secondary"]])

  # background curve: this gene's own MI distribution, standardized
  gx0 <- x0 + 1.3; gx1 <- x0 + W - 0.3; gy0 <- 2.45; gh <- 0.95
  zz <- seq(-3.2, 3.2, length.out = 120)
  gxs <- gx0 + (zz + 3.2) / 6.4 * (gx1 - gx0)
  grid.polygon(c(gxs, rev(gxs)), c(gy0 + gh * exp(-zz^2 / 2), rep(gy0, length(zz))),
               default.units = "native",
               gp = gpar(fill = pal[["orange"]], alpha = 0.22, col = NA))
  grid.lines(gxs, gy0 + gh * exp(-zz^2 / 2), default.units = "native",
             gp = gpar(col = pal[["orange"]], lwd = 2.4))
  seg(gx0, gy0, gx1, gy0, ink[["muted"]], 1)
  zt <- c(-1.6, -1.1, -0.7, -0.4, -0.15, 0.1, 0.35, 0.6, 1.0, 1.4, 2.0, 2.75)
  for (z in zt) seg(gx0 + (z + 3.2) / 6.4 * (gx1 - gx0), gy0 - 0.12,
                    gx0 + (z + 3.2) / 6.4 * (gx1 - gx0), gy0 + 0.02,
                    pal[["blue"]], ifelse(z > 2.5, 3, 1.5))
  zx <- gx0 + (2.75 + 3.2) / 6.4 * (gx1 - gx0)
  txt("z", zx, gy0 + 0.32, 14, face = "bold.italic")
  txt("gene i's own\nbackground", x0 + 0.28, gy0 + 0.42, 12, ink[["secondary"]], just = "left")

  # analogy box
  bxl <- x0 + 0.25; bxw <- W - 0.5
  grid.roundrect(bxl, 1.02, bxw, 0.82, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.06, "inches"), gp = gpar(fill = "#f3f3f1", col = ink[["grid"]]))
  txt("dependence ≈ query · key", cx, 1.62, 13)
  txt("context normalization ≈ scaled logits", cx, 1.24, 13)
  txt("CLR (2007): attention before transformers", cx, 0.62, 13, ink[["secondary"]],
      face = "italic")

  # ---------------- panel 3: softmax + depth ---------------------------------
  x0 <- X[3]; cx <- x0 + 1.25; cy <- 5.55
  nb <- ring(22, cx, cy, 0.95, phase = 0.1)
  top <- c(1, 3, 5, 6, 9, 11, 13, 16, 18, 21)   # ~10 partners
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.17, 0.08)
    if (k %in% top) seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 3.4,
                        arrow = arrow(length = unit(0.07, "inches"), type = "closed"))
    else seg(s$x0, s$y0, s$x1, s$y1, ink[["muted"]], 0.8, alpha = 0.45)
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.065,
                                    col = if (k %in% top) ink[["secondary"]] else ink[["muted"]],
                                    alpha = if (k %in% top) 1 else 0.7)
  node(cx, cy, 0.17, col = pal[["orange"]], lwd = 3)
  txt("i", cx, cy, 13, face = "bold.italic")
  txt("softmax", x0 + 2.42, 6.0, 14, face = "bold", just = "left")
  txt("weight on\n~10 partners", x0 + 2.42, 5.42, 12.5, ink[["secondary"]], just = "left")

  # multi-step walk halos over a small network
  set.seed(11)
  hcx <- x0 + 1.25; hcy <- 2.55
  radii <- c(0.38, 0.72, 1.06)
  for (j in 3:1) grid.circle(hcx, hcy, radii[j], default.units = "native",
                             gp = gpar(fill = pal[["magenta"]], alpha = c(0.32, 0.2, 0.11)[j],
                                       col = NA))
  for (j in 1:3) grid.circle(hcx, hcy, radii[j], default.units = "native",
                             gp = gpar(fill = NA, col = pal[["magenta"]], lwd = 1.4, lty = 2))
  ang <- runif(30, 0, 2 * pi); rr <- sqrt(runif(30)) * 1.0 + 0.12
  pts <- data.frame(x = hcx + rr * cos(ang), y = hcy + rr * sin(ang))
  D <- as.matrix(dist(rbind(c(hcx, hcy), pts)))
  for (a in 1:30) for (b in (a + 1):31) if (b <= 31 && D[a, b] < 0.42)
    seg(c(hcx, pts$x)[a], c(hcy, pts$y)[a], c(hcx, pts$x)[b], c(hcy, pts$y)[b],
        ink[["secondary"]], 0.9, alpha = 0.55)
  for (k in 1:30) node(pts$x[k], pts$y[k], 0.045, col = ink[["secondary"]], lwd = 1)
  node(hcx, hcy, 0.1, col = pal[["orange"]], lwd = 2.6)

  txt("multi-step\nwalk", x0 + 2.42, 3.25, 14, face = "bold", just = "left")
  # depth scale (log3): t = 1, 3, 9, 27
  sx0 <- x0 + 2.5; sx1 <- x0 + W - 0.28; sy <- 2.15
  tx <- sx0 + (0:3) / 3 * (sx1 - sx0)
  seg(sx0, sy, sx1, sy, pal[["magenta"]], 4)
  for (k in 1:4) {
    seg(tx[k], sy - 0.1, tx[k], sy + 0.1, ink[["secondary"]], 1.3)
    txt(c("1", "3", "9", "27")[k], tx[k], sy - 0.3, 12.5, ink[["secondary"]])
  }
  txt("depth t (log)", sx0 - 0.05, sy + 0.32, 12.5, ink[["secondary"]], just = "left")

  # pill: no learned parameters
  grid.roundrect(x0 + W / 2, 0.72, 2.9, 0.46, default.units = "native",
                 r = unit(0.23, "inches"), gp = gpar(fill = "white", col = ink[["primary"]],
                                                     lwd = 1.4))
  txt("no learned parameters", x0 + W / 2, 0.72, 13.5, face = "bold")

  # ---------------- panel 4: recovery + bar strip -----------------------------
  x0 <- X[4]
  # regulon (left)
  rcx <- x0 + 0.98; rcy <- 5.75
  tg <- ring(13, rcx, rcy, 0.68, phase = 0.3)
  tg$x <- tg$x + c(0.08, -0.04, 0.06, -0.08, 0.05, 0.02, -0.06, 0.07, -0.03, 0.04, -0.07, 0.03, 0)
  for (k in 1:13) { s <- trim(rcx, rcy, tg$x[k], tg$y[k], 0.15, 0.06)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["aqua"]], 1.8) }
  for (k in 1:13) node(tg$x[k], tg$y[k], 0.06, fill = "#d1efe4", col = pal[["aqua"]], lwd = 1.3)
  node(rcx, rcy, 0.15, fill = pal[["aqua"]], col = "#0d6e4c", lwd = 1.5)
  txt("TF", rcx, rcy, 11, "white", "bold")
  check_mark(rcx + 0.72, rcy + 0.78)
  txt("large regulatory\nprograms ↑", rcx, 4.62, 13, face = "bold")
  txt("E. coli regulons", rcx, 4.1, 12, ink[["secondary"]], face = "italic")

  # tight complex (right)
  ccx <- x0 + 2.78; ccy <- 5.75
  cl <- ring(5, ccx, ccy, 0.34, phase = pi / 2)
  for (a in 1:4) for (b in (a + 1):5) seg(cl$x[a], cl$y[a], cl$x[b], cl$y[b], ink[["muted"]], 1.8)
  for (k in 1:5) node(cl$x[k], cl$y[k], 0.075, fill = "#e9e8e4", col = ink[["muted"]], lwd = 1.3)
  neutral_mark(ccx + 0.52, ccy + 0.78)
  txt("tight complexes:\nCLR better", ccx, 4.62, 13, face = "bold")
  txt(sprintf("K562 CORUM:\nattention %.2f · CLR %.2f", corum[["attention"]], corum[["clr"]]),
      ccx, 3.98, 11.5, ink[["secondary"]], lineheight = 1.0)

  # bar strip: M3D strong-evidence regulons, AUPR / base rate
  seg(x0 + 0.25, 3.62, x0 + W - 0.25, 3.62, ink[["grid"]], 1.2)
  txt("E. coli (M3D), strong-evidence regulons", x0 + 0.25, 3.38, 12.5, face = "bold", just = "left")
  txt("AUPR / base rate (1 = random); t = 9", x0 + 0.25, 3.08, 12, ink[["secondary"]], just = "left")
  lab_x <- x0 + 0.25; bx0 <- x0 + 1.95; bxs <- (x0 + W - 0.6 - bx0) / 1.8   # scale 0..1.8
  rows_y <- 2.68 - (0:4) * 0.39; bh <- 0.24
  for (k in 1:4) {
    txt(bars$label[k], lab_x, rows_y[k], 11.5, face = if (k == 1) "bold" else "plain",
        just = "left")
    grid.rect(bx0, rows_y[k], bars$value[k] * bxs, bh, default.units = "native",
              just = c("left", "centre"), gp = gpar(fill = bars$col[k], col = NA))
    txt(sprintf("%.2f", bars$value[k]), bx0 + bars$value[k] * bxs + 0.06, rows_y[k], 11.5,
        just = "left", face = if (k == 1) "bold" else "plain")
  }
  # placeholder: learned attention (scFM), pending
  txt("learned attn. (scFM)", lab_x, rows_y[5], 11, ink[["secondary"]], just = "left")
  grid.rect(bx0, rows_y[5], 1.8 * bxs, bh, default.units = "native", just = c("left", "centre"),
            gp = gpar(fill = NA, col = ink[["muted"]], lty = 2, lwd = 1.3))
  txt("pending", bx0 + 0.9 * bxs, rows_y[5], 11, ink[["muted"]], face = "italic")
  # random line at 1
  seg(bx0 + bxs, rows_y[1] + 0.25, bx0 + bxs, rows_y[4] - 0.17, ink[["secondary"]], 1, lty = 3)
  seg(bx0, rows_y[1] + 0.2, bx0, rows_y[5] - 0.2, ink[["muted"]], 1)

  popViewport()
}

# ---- render -------------------------------------------------------------------
png(file.path(out_dir, "graphical_abstract.png"), width = 2400, height = 1350, res = 150,
    type = "cairo", bg = "white")
draw_abstract(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract.pdf"), width = 16, height = 9, bg = "white")
draw_abstract(); invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract.{png,pdf}"), "\n")
