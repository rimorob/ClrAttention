# Graphical abstract (PORTRAIT, W_CAN x H_CAN in, printed 7.0 in wide), drawn entirely in code with grid (no generative imagery).
# Brief: docs/whitepapers/graphical_abstract_prompt.md, revised after PI review:
# published CLR (Faith et al. 2007) is ONLY the per-gene context normalization;
# the softmax and the multi-step walk are this paper's additions. Panel 2 (PI wording, 19 Nov): the 2007 paper (Faith et al.)
# describes only the Euclidean combination sqrt(zi^2 + zj^2); the released 2007 software offered both the Euclidean and the
# Stouffer combination (zi + zj)/sqrt(2); the Stouffer combination is the one used here. No statement about who introduced
# Stouffer or when. (The B-spline MI estimator is the same family in both, Daub et al. 2004.)
# Names (19 Nov): aCLR = the adopted method (earlier text and file names: "attention"), tcaCLR = the log + second CLR read-out
# (earlier: "a3r"). The stacked variants of tcaCLR (two and three layers) are not shown on this figure.
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

out_dir <- Sys.getenv("SB_OUT", "docs/whitepapers/figures"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
FONT <- "Liberation Sans"   # Arial-metric sans

# ---- numbers (docs/whitepapers/sota_table.md, raw scores) -------------------
this_lab <- "aCLR (adopted)"
# The tcaCLR row is a value of results/explore_stack/{m3d,k562}/pooled_levels.csv
# (raw scope; same base rate as the other bars). Not adopted: 26 Oct preregistered test; K562 is in-sample for it.
# ---- tcaCLR2 (exploratory): same metric and base-rate scaling as the other bars (raw pooled AUPR / base rate), from results/tcaclr2_scoreboard/cells_tcaclr2_all.csv ----
T2_NAME <- "tcaCLR2"
t2c <- utils::read.csv(Sys.getenv("SB_T2", "results/tcaclr2_scoreboard/cells_tcaclr2_all.csv"), stringsAsFactors = FALSE)
t2e <- function(ds, truth) { x <- t2c[t2c$dataset == ds & t2c$truth == truth & t2c$metric == "pooled_aupr_over_base" & t2c$scope == "raw" & t2c$method == "tcaCLR2", ]; stopifnot(nrow(x) == 1L); x$enrichment_over_random }
T2_M3D <- t2e("m3d", "SC"); T2_K562 <- t2e("k562", "CORUM")
a3r_lab <- "tcaCLR (lead candidate)"
strip_a <- data.frame(   # M3D, strong-evidence (SC) RegulonDB regulons
  label = c(this_lab, a3r_lab, "tcaCLR2 (exploratory)", "GENIE3", "PIDC", "CLR"),
  value = c(1.68, 1.72, T2_M3D, 1.51, 1.45, 1.44),
  col = c(pal[["aqua"]], "#bfe6d8", "#bfe6d8", ink[["muted"]], ink[["muted"]], pal[["orange"]]),
  stringsAsFactors = FALSE)
strip_b <- data.frame(   # K562, CORUM complexes
  label = c("aCLR", "tcaCLR (picked on K562)", "tcaCLR2 (exploratory)", "CLR", "GENIE3", "|r| (correlation)"),
  value = c(3.51, 7.95, T2_K562, 6.39, 6.32, 7.00),
  col = c(pal[["aqua"]], "#bfe6d8", "#bfe6d8", pal[["orange"]], ink[["muted"]], ink[["muted"]]),
  stringsAsFactors = FALSE)

EXPL <- "#bfe6d8"
strip_a$bold <- c(TRUE, rep(FALSE, 5)); strip_b$bold <- c(TRUE, rep(FALSE, 5))

# ---- PORTRAIT layout: canvas W_CAN x H inches (H from the layout below), printed 7.0 in wide: scale 7.0 / W_CAN = 0.843;
# every font size here is >= 8.5 pt, i.e. >= 7.1 pt printed. Panels 1-3 sit in one row, panel 4 (the walk) is a full-width row,
# the transformer analogy follows, and panel 5 (what it recovers) closes the figure.
W_CAN <- 8.3
FS_S <- 8.5                                                  # smallest text (canvas pt)
txt <- function(label, x, y, size = 13, col = ink[["primary"]], face = "plain",
                just = "centre", rot = 0, lineheight = 1.0)
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
  txt("i", x, y, 11, "white", "bold.italic")
}
ring <- function(n, cx, cy, r, phase = 0)
  data.frame(x = cx + r * cos(phase + 2 * pi * (seq_len(n) - 1) / n),
             y = cy + r * sin(phase + 2 * pi * (seq_len(n) - 1) / n))
trim <- function(x0, y0, x1, y1, r0, r1) {
  d <- sqrt((x1 - x0)^2 + (y1 - y0)^2); ux <- (x1 - x0) / d; uy <- (y1 - y0) / d
  list(x0 = x0 + ux * r0, y0 = y0 + uy * r0, x1 = x1 - ux * r1, y1 = y1 - uy * r1)
}
text_w <- function(s, size, face = "plain") convertWidth(grobWidth(textGrob(s, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face))), "inches", TRUE)
wrap_to <- function(s, size, width, face = "plain") {     # word wrap (never inside a word) to a width in inches
  words <- strsplit(s, " ", fixed = TRUE)[[1]]; lines <- character(); cur <- ""
  for (w in words) { cand <- if (nzchar(cur)) paste(cur, w) else w
    if (!nzchar(cur) || text_w(cand, size, face) <= width) cur <- cand else { lines <- c(lines, cur); cur <- w } }
  c(lines, cur) }
# card: (x, yt) = top-left corner, w, h; number circle and title (top-left)
card <- function(x, yt, w, h, num, title, size = 10.5) {
  grid.roundrect(x, yt - h, w, h, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.2))
  grid.circle(x + 0.21, yt - 0.23, 0.13, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
  txt(num, x + 0.21, yt - 0.235, 9.5, "white", "bold")
  grid.text(title, x + 0.43, yt - 0.1, default.units = "native", just = c("left", "top"),
            gp = gpar(fontfamily = FONT, fontsize = size, fontface = "bold", col = ink[["primary"]], lineheight = 1.0))
}
chevron <- function(x, y, col = ink[["muted"]])
  grid.polygon(c(x - 0.04, x + 0.06, x - 0.04), c(y + 0.11, y, y - 0.11),
               default.units = "native", gp = gpar(fill = col, col = NA))
pill <- function(label, x, y, w, h = 0.3, fill = "white", col = ink[["primary"]], size = 10,
                 face = "bold", tcol = ink[["primary"]], lwd = 1.2) {
  grid.roundrect(x, y, w, h, default.units = "native", r = unit(h / 2, "inches"),
                 gp = gpar(fill = fill, col = col, lwd = lwd))
  txt(label, x, y, size, tcol, face)
}
band <- function(x0, x1, ytop, fill, label, hb = 0.36, size = 11.5) {   # labelled band + bracket over a group of panels; (ytop = top edge)
  grid.roundrect(x0, ytop - hb, x1 - x0, hb, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.07, "inches"), gp = gpar(fill = fill, col = NA))
  txt(label, x0 + 0.14, ytop - hb / 2, size, face = "bold", just = "left", lineheight = 0.95)
  by <- ytop - hb - 0.09
  grid.lines(c(x0 + 0.04, x0 + 0.04, x1 - 0.04, x1 - 0.04), c(by - 0.07, by, by, by - 0.07),
             default.units = "native", gp = gpar(col = ink[["secondary"]], lwd = 1.1))
}
# dotted green outline (new in this work); (x, yt) = top-left corner
newbox <- function(x, yt, w, h, r = 0.07)
  grid.roundrect(x, yt - h, w, h, default.units = "native", just = c("left", "bottom"),
                 r = unit(r, "inches"), gp = DOT)

# bar strip: rows (label only; the lead rows are dashed); yt = first row centre
bar_strip <- function(d, x0, w, y_first, xmax, ticks, loss = FALSE, pitch = 0.19) {
  lab_x <- x0 + 0.1; bx0 <- x0 + 1.75; bx1 <- x0 + w - 0.45
  sc <- (bx1 - bx0) / xmax
  n <- nrow(d); ry <- y_first - (seq_len(n) - 1) * pitch; bh <- 0.14
  for (t in ticks) seg(bx0 + t * sc, ry[1] + 0.15, bx0 + t * sc, ry[n] - 0.15, ink[["grid"]], 0.8)
  for (k in seq_len(n)) {
    txt(d$label[k], lab_x, ry[k], 9.3, face = if (d$bold[k]) "bold" else "plain", just = "left")
    lead_row <- d$col[k] == EXPL
    grid.rect(bx0, ry[k], d$value[k] * sc, bh, default.units = "native", just = c("left", "centre"),
              gp = gpar(fill = d$col[k], col = if (lead_row) pal[["aqua"]] else NA, lwd = 1.0,
                        lty = if (lead_row) 2 else 1))
    if (loss && k == 1) txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc - 0.04, ry[k], 9, "white", "bold", just = "right")
    else txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc + 0.04, ry[k], 9.3,
        face = if (d$bold[k]) "bold" else "plain", just = "left")
  }
  seg(bx0 + sc, ry[1] + 0.2, bx0 + sc, ry[n] - 0.15, ink[["primary"]], 1.0, lty = 2)
  txt("random = 1", bx0 + sc + 0.04, ry[1] + 0.22, FS_S, ink[["secondary"]], just = "left")
  seg(bx0, ry[1] + 0.15, bx0, ry[n] - 0.15, ink[["secondary"]], 0.9)
  ay <- ry[n] - 0.15
  seg(bx0, ay, bx1, ay, ink[["secondary"]], 0.9)
  for (t in ticks) { seg(bx0 + t * sc, ay, bx0 + t * sc, ay - 0.04, ink[["secondary"]], 0.9)
    txt(as.character(t), bx0 + t * sc, ay - 0.13, FS_S, ink[["secondary"]]) }
  txt("precision\u2013recall area as a multiple of random", (bx0 + bx1) / 2 - 0.1, ay - 0.29, FS_S, ink[["secondary"]])
  if (loss) {   # red loss marker: from the CLR value back to this work's value
    ly <- ry[1]; ci <- which(d$label == "CLR")
    seg(bx0 + d$value[ci] * sc, ly, bx0 + d$value[1] * sc + 0.04, ly, pal[["red"]], 2.2,
        arrow = arrow(length = unit(0.07, "inches"), type = "closed"))
    seg(bx0 + d$value[ci] * sc, ly - 0.07, bx0 + d$value[ci] * sc, ly + 0.07, pal[["red"]], 2.2)
    txt(sprintf("%.1f\u00d7 lower", d$value[ci] / d$value[1]), bx0 + d$value[ci] * sc + 0.06, ly, 9, face = "bold", just = "left")
  }
  ay - 0.29      # y of the axis title
}

# ---- vertical layout (distance from the top of the canvas) --------------------------------------------------------------------
D_BAND <- 0.82; HB <- 0.46; D_R1 <- 1.52; H1 <- 4.1       # band top, band height, row 1 top, row 1 height
D_NOTE <- D_R1 + H1 + 0.14
D_R2 <- D_NOTE + 0.17; H2 <- 1.95                 # row 2 (the walk)
D_ANA <- D_R2 + H2 + 0.07; HANA <- 0.3           # transformer analogy
D_R3 <- D_ANA + HANA + 0.09; H3 <- 3.35           # results panel
H_CAN <- D_R3 + H3 + 0.1

# ---- drawing ----------------------------------------------------------------
draw_abstract <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  H <- H_CAN
  pushViewport(viewport(xscale = c(0, W_CAN), yscale = c(0, H)))
  T <- H
  MG <- 0.12; XR <- W_CAN - 0.12; TW <- XR - MG

  # headline: the finding, with its cost
  grid.text(expression(bold("Adding a softmax and a multi-step walk to CLR's normalization recovers large ") *
                         bolditalic("E. coli")),
            MG, T - 0.27, default.units = "native", just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 13.5, col = ink[["primary"]]))
  txt("regulons better (16% over CLR), but loses resolution on tight protein complexes",
      MG, T - 0.52, 13.5, face = "bold", just = "left")

  # ---- row 1: panels 1-3 ----
  gap <- 0.12; Wv <- c(1.8, 3.0, 3.0)
  X <- MG + c(0, cumsum(Wv[1:2] + gap))
  yt <- T - D_R1
  band(X[1], X[2] + Wv[2], T - D_BAND, "#ecebe7", "CLR (Faith et al. 2007)", HB)
  band(X[3], X[3] + Wv[3], T - D_BAND, "#e6e3fa", "aCLR (this work):\nsteps taken from attention", HB, 10.5)
  card(X[1], yt, Wv[1], H1, 1, "Expression\ncompendium")
  card(X[2], yt, Wv[2], H1, 2, "Mutual information,\nz-scored against each\ngene's own background")
  card(X[3], yt, Wv[3], H1, 3, "Softmax with\nper-gene\ntemperature")
  for (i in 1:2) chevron(X[i] + Wv[i] + gap / 2, yt - H1 / 2)

  # ---------------- panel 1: heat-map motif + organisms ----------------------
  x0 <- X[1]
  set.seed(7)
  ng <- 12; ne <- 14
  mod <- rep(1:3, c(5, 4, 3))
  prof <- matrix(rnorm(3 * ne), 3)
  M <- prof[mod, ] * 0.95 + matrix(rnorm(ng * ne, sd = 0.55), ng)
  M <- (M - min(M)) / diff(range(M))
  ramp <- colorRampPalette(seq_blue[1:7])(101)
  hx0 <- x0 + 0.44; hx1 <- x0 + Wv[1] - 0.12; hy1 <- yt - 0.7; hy0 <- yt - 1.95
  cw <- (hx1 - hx0) / ne; ch <- (hy1 - hy0) / ng
  for (g in 1:ng) for (e in 1:ne)
    grid.rect(hx0 + (e - 1) * cw, hy1 - g * ch, cw * 0.9, ch * 0.88, default.units = "native",
              just = c("left", "bottom"), gp = gpar(fill = ramp[1 + round(100 * M[g, e])], col = NA))
  txt("genes", x0 + 0.27, (hy0 + hy1) / 2, 9.5, ink[["secondary"]], rot = 90)
  txt("experiments", hx0, hy0 - 0.17, 9.5, ink[["secondary"]], just = "left")
  seg(hx0 + 0.78, hy0 - 0.17, hx1, hy0 - 0.17, ink[["muted"]], 1.0,
      arrow = arrow(length = unit(0.05, "inches"), type = "closed"))
  by <- yt - 2.58; bx <- x0 + 0.4
  grid.curve(bx + 0.26, by, bx + 0.44, by - 0.04, default.units = "native", curvature = -0.6,
             square = FALSE, ncp = 8, gp = gpar(col = ink[["secondary"]], lwd = 1.1))
  grid.roundrect(bx, by, 0.52, 0.24, default.units = "native", r = unit(0.12, "inches"),
                 gp = gpar(fill = seq_blue[3], col = ink[["secondary"]], lwd = 1.2))
  txt("E. coli", bx + 0.04, by - 0.36, 10.5, face = "italic")
  txt("M3D, EcoMAC\nmicroarrays", bx + 0.04, by - 0.72, FS_S, ink[["secondary"]])
  hx <- x0 + 1.37; hy <- by
  grid.circle(hx, hy, 0.21, default.units = "native", gp = gpar(fill = seq_blue[2], col = ink[["secondary"]], lwd = 1.2))
  grid.circle(hx + 0.04, hy + 0.01, 0.08, default.units = "native", gp = gpar(fill = seq_blue[5], col = ink[["secondary"]], lwd = 1))
  txt("human", hx, hy - 0.36, 10.5)
  txt("K562\nPerturb-seq", hx, hy - 0.72, FS_S, ink[["secondary"]])

  # ---------------- panel 2: CLR criss-cross (after Faith et al. 2007, Fig. 2A) ----
  # MI matrix with gene i's row and gene j's column; each is z-scored against its own
  # background; the pair score combines the two (Stouffer (zi + zj)/sqrt(2); Euclidean sqrt(zi^2 + zj^2) ranks nearly identically; negatives clipped to 0).
  x0 <- X[2]; W2 <- Wv[2]
  pill("existing method", x0 + 1.0, yt - 0.74, 1.5, 0.24, fill = "#ecebe7", col = NA, size = 9.3,
       face = "plain", tcol = ink[["secondary"]])
  set.seed(5)
  n <- 8; cs <- 0.115; mx0 <- x0 + 0.95; my1 <- yt - 1.1; my0 <- my1 - n * cs; mx1 <- mx0 + n * cs
  A <- matrix(runif(n * n)^1.6, n); A <- (A + t(A)) / 2
  ri <- 3; cj <- 7; A[ri, cj] <- A[cj, ri] <- 0.97
  ramp <- colorRampPalette(seq_blue[1:7])(101)
  for (r in 1:n) for (c in 1:n)
    grid.rect(mx0 + (c - 1) * cs, my1 - r * cs, cs * 0.9, cs * 0.9, default.units = "native",
              just = c("left", "bottom"),
              gp = gpar(fill = if (r == c) ink[["grid"]] else ramp[1 + round(100 * A[r, c])], col = NA))
  rowy <- my1 - (ri - 0.5) * cs - 0.005; colx <- mx0 + (cj - 0.5) * cs - 0.005
  grid.rect(mx0 - 0.025, rowy, n * cs + 0.035, cs + 0.04, default.units = "native",
            just = c("left", "centre"), gp = gpar(fill = NA, col = pal[["orange"]], lwd = 2.0))
  grid.rect(colx, my0 - 0.025, cs + 0.04, n * cs + 0.035, default.units = "native",
            just = c("centre", "bottom"), gp = gpar(fill = NA, col = pal[["orange"]], lwd = 2.0))
  grid.rect(colx, rowy, cs * 0.95, cs * 0.95, default.units = "native",
            gp = gpar(fill = NA, col = ink[["primary"]], lwd = 2.2))
  txt("i", mx0 - 0.12, rowy, 11, face = "bold.italic")
  txt("j", colx, my1 + 0.1, 11, face = "bold.italic")
  txt("MI of\nevery\ngene\npair", mx1 + 0.1, (my0 + my1) / 2, 8.8, ink[["secondary"]], just = "left")
  # the two backgrounds: row i (z_i) and column j (z_j)
  cb <- yt - 2.5; chh <- 0.3
  curves <- list(c(x0 + 0.15, x0 + 1.3, 2.3), c(colx - 0.55, colx + 0.55, 1.6))
  zz <- seq(-3.2, 3.2, length.out = 100)
  zt <- c(-1.7, -1.2, -0.8, -0.5, -0.2, 0.05, 0.3, 0.6, 0.95, 1.3)
  for (k in 1:2) {
    a <- curves[[k]][1]; b <- curves[[k]][2]; zm <- curves[[k]][3]
    xs <- a + (zz + 3.2) / 6.4 * (b - a)
    grid.polygon(c(xs, rev(xs)), c(cb + chh * exp(-zz^2 / 2), rep(cb, length(zz))),
                 default.units = "native", gp = gpar(fill = pal[["orange"]], alpha = 0.22, col = NA))
    grid.lines(xs, cb + chh * exp(-zz^2 / 2), default.units = "native",
               gp = gpar(col = pal[["orange"]], lwd = 1.8))
    seg(a, cb, b, cb, ink[["muted"]], 0.9)
    for (z in zt) { xx <- a + (z + 3.2) / 6.4 * (b - a); seg(xx, cb - 0.07, xx, cb + 0.01, pal[["blue"]], 1.0) }
    xm <- a + (zm + 3.2) / 6.4 * (b - a)
    seg(xm, cb - 0.1, xm, cb + 0.12, ink[["primary"]], 2.3)
  }
  lc <- mean(curves[[1]][1:2]); rc <- mean(curves[[2]][1:2])
  grid.text(expression(z[i] * "  (row " * italic(i) * ")"), lc, cb - 0.2, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 9.3, col = ink[["primary"]]))
  grid.text(expression(z[j] * "  (column " * italic(j) * ")"), rc, cb - 0.2, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 9.3, col = ink[["primary"]]))
  # connectors from the highlighted row / column to their backgrounds
  ar <- arrow(length = unit(0.05, "inches"), type = "closed")
  grid.lines(c(mx0 - 0.2, lc, lc), c(rowy, rowy, cb + chh + 0.07), default.units = "native",
             arrow = ar, gp = gpar(col = ink[["muted"]], lwd = 1.1, lty = 3, fill = ink[["muted"]]))
  seg(colx, my0 - 0.04, rc, cb + chh + 0.07, ink[["muted"]], 1.1, lty = 3, arrow = ar)
  # combination into one score
  fx <- x0 + W2 / 2
  seg(lc, cb - 0.3, fx - 0.3, yt - 3.0, ink[["secondary"]], 1.1, arrow = ar)
  seg(rc, cb - 0.3, fx + 0.3, yt - 3.0, ink[["secondary"]], 1.1, arrow = ar)
  grid.text(expression((z[i] + z[j]) / sqrt(2)), fx, yt - 3.15, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 12, col = ink[["primary"]]))
  txt("Stouffer combination (used here); z clipped at 0", fx, yt - 3.4, FS_S, ink[["secondary"]])
  grid.text(expression("2007 paper: Euclidean " * sqrt(z[i]^2 + z[j]^2) * "; 2007 software: both"), fx, yt - 3.6, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = FS_S, col = ink[["secondary"]]))
  grid.roundrect(fx, yt - 3.85, 2.35, 0.26, default.units = "native", r = unit(0.05, "inches"),
                 gp = gpar(fill = "#fbe3d9", col = pal[["orange"]], lwd = 1.3))
  txt("CLR score for the pair i–j", fx, yt - 3.85, 9.8, face = "bold")

  # ---------------- panel 3: softmax ------------------------------------------
  x0 <- X[3]; cx <- x0 + Wv[3] / 2; cy <- yt - 1.95
  nb <- ring(20, cx, cy, 0.85, phase = 0.1)
  top <- c(1, 3, 5, 6, 9, 11, 13, 15, 17, 19)   # ~10 partners
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.13, 0.065)
    if (k %in% top) seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 2.6,
                        arrow = arrow(length = unit(0.06, "inches"), type = "closed"))
    else seg(s$x0, s$y0, s$x1, s$y1, ink[["muted"]], 0.7, alpha = 0.45)
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.055,
                                    col = if (k %in% top) ink[["secondary"]] else ink[["muted"]],
                                    alpha = if (k %in% top) 1 else 0.7)
  query(cx, cy, 0.13)
  txt("~10 effective partners\nper gene, set in advance", cx, yt - 3.35, 9.8, ink[["secondary"]])

  # "Published CLR stops here" under panels 1-2
  txt("Published CLR stops here: a score for each gene pair", (X[1] + X[2] + Wv[2]) / 2, T - D_NOTE, 9.8, ink[["secondary"]], face = "italic")

  # ---------------- row 2, panel 4: multi-step walk (full width) ---------------------------
  yq <- T - D_R2; x0 <- MG; W4 <- TW
  card(x0, yq, W4, H2, 4, "Multi-step walk,\nnumber of steps\nchosen without\nlabels")
  hcx <- x0 + 3.3; hcy <- yq - 1.0
  set.seed(11)
  radii <- c(0.3, 0.6, 0.9)
  for (j in 3:1) grid.circle(hcx, hcy, radii[j], default.units = "native",
                             gp = gpar(fill = pal[["magenta"]], alpha = c(0.32, 0.2, 0.11)[j], col = NA))
  for (j in 1:3) grid.circle(hcx, hcy, radii[j], default.units = "native",
                             gp = gpar(fill = NA, col = pal[["magenta"]], lwd = 1.2, lty = 2))
  ang <- runif(28, 0, 2 * pi); rr <- sqrt(runif(28)) * 0.85 + 0.12
  pts <- data.frame(x = c(hcx, hcx + rr * cos(ang)), y = c(hcy, hcy + rr * sin(ang)))
  D <- as.matrix(dist(pts))
  for (a in 1:28) for (b in (a + 1):29) if (D[a, b] < 0.36)
    seg(pts$x[a], pts$y[a], pts$x[b], pts$y[b], ink[["secondary"]], 0.7, alpha = 0.55)
  for (k in 2:29) node(pts$x[k], pts$y[k], 0.035, col = ink[["secondary"]], lwd = 0.8)
  query(hcx, hcy, 0.11)
  # steps scale (log3): t = 1, 3, 9, 27
  sx0 <- x0 + 4.2; sx1 <- x0 + 5.95; sy <- yq - 1.05
  tx <- sx0 + (0:3) / 3 * (sx1 - sx0)
  txt("steps t (log scale)", (sx0 + sx1) / 2, sy + 0.3, FS_S, ink[["secondary"]])
  seg(sx0, sy, sx1, sy, pal[["magenta"]], 3.5)
  for (k in 1:4) {
    seg(tx[k], sy - 0.07, tx[k], sy + 0.07, ink[["secondary"]], 1.1)
    txt(c("1", "3", "9", "27")[k], tx[k], sy - 0.2, 9, ink[["secondary"]], face = if (k == 3) "bold" else "plain")
  }
  grid.circle(tx[3], sy, 0.055, default.units = "native", gp = gpar(fill = ink[["primary"]], col = "white"))
  txt("t = 9: chosen on\nheld-out labs", tx[3], sy - 0.55, FS_S, ink[["secondary"]])
  pill("No learned parameters", x0 + W4 - 0.95, yq - H2 / 2, 1.75, 0.3, size = 9.8)

  # transformer analogy strip
  ya <- T - D_ANA
  grid.roundrect(MG, ya - HANA, TW, HANA, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.07, "inches"), gp = gpar(fill = "#f3f3f1", col = ink[["grid"]]))
  lab <- "Transformer analogy:"
  txt(lab, MG + 0.14, ya - HANA / 2, 9.5, face = "bold", just = "left")
  lw <- text_w(lab, 9.5, "bold")
  txt("MI ≈ query·key score;  CLR z-score ≈ scaled logit;  softmax ≈ softmax;  walk ≈ repeated attention with residual",
      MG + 0.14 + lw + 0.08, ya - HANA / 2, 9, ink[["secondary"]], just = "left")

  # ---------------- panel 5: results (full width, two blocks side by side) -----------------------------------------
  y5 <- T - D_R3; x5 <- MG; W5 <- TW
  grid.roundrect(x5, y5 - H3, W5, H3, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.2))
  grid.circle(x5 + 0.21, y5 - 0.23, 0.13, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
  txt("5", x5 + 0.21, y5 - 0.235, 9.5, "white", "bold")
  txt("What it recovers", x5 + 0.43, y5 - 0.23, 11.5, face = "bold", just = "left")
  d0 <- y5 - 0.42
  seg(x5 + 4.02, d0 - 0.05, x5 + 4.02, y5 - H3 + 0.12, ink[["grid"]], 1.0)
  # (a) regulons
  bxa <- x5 + 0.1; wa <- 3.85
  ix <- bxa + wa - 0.35; iy <- d0 - 0.3
  tg <- ring(9, ix, iy, 0.2, phase = 0.2)
  for (k in 1:9) seg(ix, iy, tg$x[k], tg$y[k], pal[["aqua"]], 1.0)
  for (k in 1:9) node(tg$x[k], tg$y[k], 0.03, fill = "#d1efe4", col = pal[["aqua"]], lwd = 0.9)
  node(ix, iy, 0.065, fill = pal[["aqua"]], col = "#0d6e4c", lwd = 0.9)
  txt("E. coli regulons (experimentally\nsupported, RegulonDB)", bxa + 0.1, d0 - 0.2, 10.5, face = "bold", just = "left", lineheight = 1.0)
  txt("regulon = genes controlled by one transcription factor (TF)", bxa + 0.1, d0 - 0.6, FS_S, ink[["secondary"]], just = "left")
  ya_ax <- bar_strip(strip_a, bxa, wa, d0 - 1.0, 2, c(0, 0.5, 1, 1.5, 2))
  txt("per regulon: tie with CLR and GENIE3;\nahead in 1,000/1,000 resamplings vs CLR", bxa + 0.1, ya_ax - 0.3, FS_S, ink[["secondary"]], face = "italic", just = "left")
  # (b) complexes
  bxb <- x5 + 4.1; wb <- 3.85
  cxi <- bxb + wb - 0.35; cyi <- d0 - 0.3
  cl <- ring(5, cxi, cyi, 0.17, phase = pi / 2)
  for (a in 1:4) for (b in (a + 1):5) seg(cl$x[a], cl$y[a], cl$x[b], cl$y[b], ink[["muted"]], 1.1)
  for (k in 1:5) node(cl$x[k], cl$y[k], 0.035, fill = "#e9e8e4", col = ink[["muted"]], lwd = 0.9)
  txt("Protein complexes (K562 CORUM;\nmostly ribosome): this work loses 1.8×", bxb + 0.1, d0 - 0.2,
      10.5, face = "bold", just = "left", lineheight = 1.0)
  txt("complex = proteins that bind as one machine", bxb + 0.1, d0 - 0.6, FS_S, ink[["secondary"]], just = "left")
  bar_strip(strip_b, bxb, wb, d0 - 1.0, 8, c(0, 2, 4, 6, 8), loss = TRUE)

  popViewport()
}

# ---- render -------------------------------------------------------------------
png(tempfile(fileext = ".png"), width = 800, height = 800, type = "cairo")   # measuring device is open while text widths are taken
png(file.path(out_dir, "graphical_abstract.png"), width = round(W_CAN * 300), height = round(H_CAN * 300), res = 300,
    type = "cairo", bg = "white")
draw_abstract(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract.pdf"), width = W_CAN, height = H_CAN, bg = "white")
draw_abstract(); invisible(dev.off())
invisible(dev.off())
cat(sprintf("p1 long: canvas %.2f x %.2f in; printed 7.00 x %.2f in (scale %.3f); min font %.1f pt -> %.1f pt printed\n", W_CAN, H_CAN, H_CAN * 7 / W_CAN, 7 / W_CAN, FS_S, FS_S * 7 / W_CAN))
cat("wrote", file.path(out_dir, "graphical_abstract.{png,pdf}"), "\n")
