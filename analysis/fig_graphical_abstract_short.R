# Graphical abstract for the SHORT pre-read, page 1 (16:9), drawn entirely in code with grid.
# Copy of analysis/fig_graphical_abstract.R, revised 15 Nov (PI design):
#  - signage: a DOTTED GREEN outline marks every element that is new in this work, including every
#    element using human data (the original CLR paper evaluated no human data); historical CLR parts are
#    labelled "Faith et al. 2007". Legend entry: "dotted green = new in this work".
#  - panel 2: the 2007 packages we found combine the clipped row z-scores as sqrt(zi^2 + zj^2) (Euclidean);
#    the version used here combines them as (zi + zj)/sqrt(2) (Stouffer). NO statement about where the
#    Stouffer combination came from, and it is not marked as new (project memory clr-2007-history).
#  - panels 3 and 4: how partners are chosen (top 50 positive CLR scores, ranked), the effective number of
#    partners (~10, per-gene temperature) and what one walk step does (rows of A; A is not symmetric).
#  - panel 5: CLR (historical), CLR attention, and the lead method; the stacked variants are dropped.
#  - the lead method's name is NOT decided: it is defined once, in LEAD_NAME below (a rename is a one-line
#    change and a re-run).
# Palette: analysis/plot_theme.R. Network drawings, heat map and the weights in panel 3 are schematic (seeded
# layouts; the weights are computed with the method's own rule from made-up scores), not data. Bars use real
# numbers from docs/whitepapers/sota_table.md and results/explore_stack (raw pooled co-membership AUPR /
# base rate; t = 9).
#
# Output: docs/whitepapers/figures/graphical_abstract_short_p1.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_short.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")

out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"   # Arial-metric sans

# ---- the lead method's name: DEFINED HERE ONLY (rename = change this line and re-run) ----------------------------------
LEAD_NAME <- "tcaCLR"                     # DEFINED HERE ONLY (temporary name; internal name a3r)
ACLR_NAME <- "aCLR"                       # DEFINED HERE ONLY (temporary name)
NAMES_NOTE <- "Method names are temporary."
LEAD_STATUS <- "lead candidate"      # honest status label used in the legend and panel 5
GREEN <- pal[["green"]]                           # dotted green = new in this work
DOT <- gpar(col = GREEN, lty = "11", lwd = 2.6, fill = NA)

# ---- numbers (docs/whitepapers/sota_table.md; results/explore_stack/*/pooled_levels.csv; raw scores) ---------------------
EXPL <- "#bfe6d8"
strip_a <- data.frame(   # M3D, strong-evidence (SC) RegulonDB regulons
  label = c("CLR (historical)", ACLR_NAME, LEAD_NAME, "GENIE3", "PIDC"),
  sub = c("", "adopted", "lead candidate", "", ""),
  value = c(1.44, 1.68, 1.72, 1.51, 1.45),
  col = c(pal[["orange"]], pal[["aqua"]], EXPL, ink[["muted"]], ink[["muted"]]),
  new = c(FALSE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
strip_b <- data.frame(   # K562, CORUM complexes
  label = c("CLR (historical)", ACLR_NAME, LEAD_NAME, "GENIE3", "|r| (correlation)"),
  sub = c("", "adopted", "lead; picked on K562", "", ""),
  value = c(6.39, 3.51, 7.95, 6.32, 7.00),
  col = c(pal[["orange"]], pal[["aqua"]], EXPL, ink[["muted"]], ink[["muted"]]),
  new = c(FALSE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)

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

# dotted green outline (new in this work); (x, y) = bottom-left corner
newbox <- function(x, y, w, h, r = 0.08)
  grid.roundrect(x, y, w, h, default.units = "native", just = c("left", "bottom"),
                 r = unit(r, "inches"), gp = DOT)

# bar strip: rows with an optional second line; dotted green around the rows that are new
bar_strip <- function(d, x0, w, y_top, xmax, ticks, loss = FALSE, pitch = 0.34) {
  lab_x <- x0 + 0.25; bx0 <- x0 + 2.0; bx1 <- x0 + w - 0.45
  sc <- (bx1 - bx0) / xmax
  n <- nrow(d); ry <- y_top - (seq_len(n) - 1) * pitch; bh <- 0.18
  for (t in ticks) seg(bx0 + t * sc, ry[1] + 0.2, bx0 + t * sc, ry[n] - 0.2, ink[["grid"]], 0.9)
  for (k in seq_len(n)) {
    has_sub <- nzchar(d$sub[k]); dy <- if (has_sub) 0.06 else 0
    if (d$new[k]) grid.roundrect(x0 + 0.12, ry[k], w - 0.24, pitch - 0.03, default.units = "native",
                                 just = c("left", "centre"), r = unit(0.05, "inches"), gp = DOT)
    txt(d$label[k], lab_x, ry[k] + dy, 11, face = if (d$new[k]) "bold" else "plain", just = "left")
    if (has_sub) txt(d$sub[k], lab_x, ry[k] - 0.09, 8.5, ink[["secondary"]], "italic", just = "left")
    lead_row <- d$col[k] == EXPL
    grid.rect(bx0, ry[k], d$value[k] * sc, bh, default.units = "native", just = c("left", "centre"),
              gp = gpar(fill = d$col[k], col = if (lead_row) pal[["aqua"]] else NA, lwd = 1.1,
                        lty = if (lead_row) 2 else 1))
    if (loss && k == 2) txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc - 0.05, ry[k], 10, "white", "bold", just = "right")
    else txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc + 0.05, ry[k], 10.5,
        face = if (d$new[k]) "bold" else "plain", just = "left")
  }
  seg(bx0 + sc, ry[1] + 0.3, bx0 + sc, ry[n] - 0.2, ink[["primary"]], 1.1, lty = 2)
  txt("random = 1", bx0 + sc + 0.04, ry[1] + 0.33, 9, ink[["secondary"]], just = "left")
  seg(bx0, ry[1] + 0.2, bx0, ry[n] - 0.2, ink[["secondary"]], 1)
  ay <- ry[n] - 0.2
  seg(bx0, ay, bx1, ay, ink[["secondary"]], 1)
  for (t in ticks) { seg(bx0 + t * sc, ay, bx0 + t * sc, ay - 0.05, ink[["secondary"]], 1)
    txt(as.character(t), bx0 + t * sc, ay - 0.16, 9.5, ink[["secondary"]]) }
  txt("precision-recall area, multiple of random", (bx0 + bx1) / 2 - 0.3, ay - 0.34, 9,
      ink[["secondary"]])
  if (loss) {   # red loss marker: from the CLR value back to CLR attention's value (row 2)
    ly <- ry[2]
    seg(bx0 + d$value[1] * sc, ly, bx0 + d$value[2] * sc + 0.04, ly, pal[["red"]], 2.4,
        arrow = arrow(length = unit(0.08, "inches"), type = "closed"))
    seg(bx0 + d$value[1] * sc, ly - 0.08, bx0 + d$value[1] * sc, ly + 0.08, pal[["red"]], 2.4)
  }
}

# schematic weights for panel 3: made-up CLR scores for one gene's 50 candidates, then the method's own rule
# (softmax(score / T), T by bisection so that exp(entropy) = 10)
sm_weights <- function(s, keff = 10) {
  f <- function(lt) { w <- exp((s - max(s)) / exp(lt)); w <- w / sum(w); exp(-sum(w * log(w))) }
  lo <- -6; hi <- 6
  for (it in 1:80) { mid <- (lo + hi) / 2; if (f(mid) < keff) lo <- mid else hi <- mid }
  w <- exp((s - max(s)) / exp((lo + hi) / 2)); w / sum(w)
}

# ---- drawing ----------------------------------------------------------------
draw_abstract <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))

  grid.text(expression(bold("Adding a softmax and a multi-step walk to CLR recovers large ") *
                         bolditalic("E. coli") * bold(" regulons better")),
            0.3, 8.62, default.units = "native", just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 18, col = ink[["primary"]]))
  txt("(16% over CLR), but does poorly on tight correlation blocks",
      0.3, 8.26, 18, face = "bold", just = "left")

  txt(paste0(LEAD_NAME, " tentatively overcomes this problem. ", NAMES_NOTE), 0.3, 7.97, 12,
      ink[["secondary"]], "bold.italic", just = "left")
  # legend (signage convention)
  lx <- 11.84
  grid.roundrect(lx, 8.2, 0.42, 0.2, default.units = "native", just = c("left", "centre"),
                 r = unit(0.04, "inches"), gp = DOT)
  txt("dotted green = new in this work", lx + 0.52, 8.2, 10.5, face = "bold", just = "left")
  grid.rect(lx, 7.94, 0.42, 0.13, default.units = "native", just = c("left", "centre"),
            gp = gpar(fill = EXPL, col = pal[["aqua"]], lwd = 1.1, lty = 2))
  txt(paste0("dashed bar = ", LEAD_NAME, ", ", LEAD_STATUS), lx + 0.52, 7.94, 9.5, ink[["secondary"]], just = "left")

  gap <- 0.18; Y <- 1.7; H <- 5.3
  Wv <- c(1.8, 3.0, 3.0, 3.0)
  X <- 0.3 + c(0, cumsum(Wv[1:3] + gap))
  xr <- X[4] + Wv[4]
  band(X[1], X[2] + Wv[2], 7.38, "#ecebe7", "CLR: Faith et al. 2007 (historical)")
  band(X[3], xr, 7.38, "#e6e3fa", paste0(ACLR_NAME, ": new in this work"))
  newbox(X[3], 7.38, xr - X[3], 0.46)

  card(X[1], Y, Wv[1], H, 1, "Expression\ncompendium", 12)
  card(X[2], Y, Wv[2], H, 2, "Mutual information,\nz-scored against each\ngene's own background", 12)
  card(X[3], Y, Wv[3], H, 3, "Softmax with a\nper-gene temperature\nover the top 50 partners", 12)
  card(X[4], Y, Wv[4], H, 4, "Multi-step walk;\nsteps chosen\nwithout labels", 12)
  newbox(X[3], Y, Wv[3], H, 0.1); newbox(X[4], Y, Wv[4], H, 0.1)
  for (i in 1:3) chevron(X[i] + Wv[i] + gap / 2, Y + H / 2)

  # ---------------- panel 1: heat-map motif + organisms ----------------------
  x0 <- X[1]
  set.seed(7)
  ng <- 12; ne <- 14
  mod <- rep(1:3, c(5, 4, 3))
  prof <- matrix(rnorm(3 * ne), 3)
  M <- prof[mod, ] * 0.95 + matrix(rnorm(ng * ne, sd = 0.55), ng)
  M <- (M - min(M)) / diff(range(M))
  ramp <- colorRampPalette(seq_blue[1:7])(101)
  hx0 <- x0 + 0.48; hx1 <- x0 + Wv[1] - 0.17; hy0 <- 4.2; hy1 <- 6.1
  cw <- (hx1 - hx0) / ne; ch <- (hy1 - hy0) / ng
  for (g in 1:ng) for (e in 1:ne)
    grid.rect(hx0 + (e - 1) * cw, hy1 - g * ch, cw * 0.9, ch * 0.88, default.units = "native",
              just = c("left", "bottom"), gp = gpar(fill = ramp[1 + round(100 * M[g, e])], col = NA))
  txt("genes", x0 + 0.3, (hy0 + hy1) / 2, 11.5, ink[["secondary"]], rot = 90)
  txt("experiments", hx0, hy0 - 0.22, 11.5, ink[["secondary"]], just = "left")
  seg(hx0 + 0.95, hy0 - 0.22, hx1, hy0 - 0.22, ink[["muted"]], 1.2,
      arrow = arrow(length = unit(0.06, "inches"), type = "closed"))
  # bacterium (E. coli) and human cell (K562; human data is new in this work)
  bx <- x0 + 0.42; by <- 3.25
  grid.curve(bx + 0.26, by, bx + 0.46, by - 0.04, default.units = "native", curvature = -0.6,
             square = FALSE, ncp = 8, gp = gpar(col = ink[["secondary"]], lwd = 1.2))
  grid.roundrect(bx, by, 0.55, 0.26, default.units = "native", r = unit(0.13, "inches"),
                 gp = gpar(fill = seq_blue[3], col = ink[["secondary"]], lwd = 1.3))
  txt("E. coli", bx + 0.04, by - 0.42, 12, face = "italic")
  txt("M3D, EcoMAC\nmicroarrays", bx + 0.04, by - 0.85, 9, ink[["secondary"]])
  hx <- x0 + 1.4; hy <- 3.25
  grid.circle(hx, hy, 0.23, default.units = "native",
              gp = gpar(fill = seq_blue[2], col = ink[["secondary"]], lwd = 1.3))
  grid.circle(hx + 0.04, hy + 0.01, 0.09, default.units = "native",
              gp = gpar(fill = seq_blue[5], col = ink[["secondary"]], lwd = 1))
  txt("human", hx, hy - 0.42, 12)
  txt("K562\nPerturb-seq", hx, hy - 0.85, 9, ink[["secondary"]])
  newbox(x0 + 1.04, 2.1, 0.7, 1.65, 0.06)      # all human data is new in this work

  # ---------------- panel 2: CLR criss-cross (after Faith et al. 2007, Fig. 2A) ----
  x0 <- X[2]; W2 <- Wv[2]
  pill("CLR as originally published", x0 + W2 / 2, 5.98, 2.35, 0.26, fill = "#ecebe7", col = NA, size = 10.5,
       face = "plain", tcol = ink[["secondary"]])
  set.seed(5)
  n <- 8; cs <- 0.135; mx0 <- x0 + 0.95; my1 <- 5.52; my0 <- my1 - n * cs; mx1 <- mx0 + n * cs
  A <- matrix(runif(n * n)^1.6, n); A <- (A + t(A)) / 2
  ri <- 3; cj <- 7; A[ri, cj] <- A[cj, ri] <- 0.97
  ramp <- colorRampPalette(seq_blue[1:7])(101)
  for (r in 1:n) for (c in 1:n)
    grid.rect(mx0 + (c - 1) * cs, my1 - r * cs, cs * 0.9, cs * 0.9, default.units = "native",
              just = c("left", "bottom"),
              gp = gpar(fill = if (r == c) ink[["grid"]] else ramp[1 + round(100 * A[r, c])], col = NA))
  rowy <- my1 - (ri - 0.5) * cs - 0.0068; colx <- mx0 + (cj - 0.5) * cs - 0.0068
  grid.rect(mx0 - 0.03, rowy, n * cs + 0.045, cs + 0.05, default.units = "native",
            just = c("left", "centre"), gp = gpar(fill = NA, col = pal[["orange"]], lwd = 2.2))
  grid.rect(colx, my0 - 0.03, cs + 0.05, n * cs + 0.045, default.units = "native",
            just = c("centre", "bottom"), gp = gpar(fill = NA, col = pal[["orange"]], lwd = 2.2))
  grid.rect(colx, rowy, cs * 0.95, cs * 0.95, default.units = "native",
            gp = gpar(fill = NA, col = ink[["primary"]], lwd = 2.4))
  txt("i", mx0 - 0.13, rowy, 12.5, face = "bold.italic")
  txt("j", colx, my1 + 0.14, 12.5, face = "bold.italic")
  txt("MI of\nevery\ngene\npair", mx1 + 0.12, (my0 + my1) / 2, 10, ink[["secondary"]], just = "left")
  cb <- 3.78; chh <- 0.36
  curves <- list(c(x0 + 0.15, x0 + 1.3, 2.3), c(colx - 0.55, colx + 0.55, 1.6))
  zz <- seq(-3.2, 3.2, length.out = 100)
  zt <- c(-1.7, -1.2, -0.8, -0.5, -0.2, 0.05, 0.3, 0.6, 0.95, 1.3)
  for (k in 1:2) {
    a <- curves[[k]][1]; b <- curves[[k]][2]; zm <- curves[[k]][3]
    xs <- a + (zz + 3.2) / 6.4 * (b - a)
    grid.polygon(c(xs, rev(xs)), c(cb + chh * exp(-zz^2 / 2), rep(cb, length(zz))),
                 default.units = "native", gp = gpar(fill = pal[["orange"]], alpha = 0.22, col = NA))
    grid.lines(xs, cb + chh * exp(-zz^2 / 2), default.units = "native",
               gp = gpar(col = pal[["orange"]], lwd = 2))
    seg(a, cb, b, cb, ink[["muted"]], 1)
    for (z in zt) { xx <- a + (z + 3.2) / 6.4 * (b - a); seg(xx, cb - 0.08, xx, cb + 0.01, pal[["blue"]], 1.1) }
    xm <- a + (zm + 3.2) / 6.4 * (b - a)
    seg(xm, cb - 0.11, xm, cb + 0.14, ink[["primary"]], 2.6)
  }
  lc <- mean(curves[[1]][1:2]); rc <- mean(curves[[2]][1:2])
  grid.text(expression(z[i] * "  (row " * italic(i) * ")"), lc, cb - 0.27, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 10.5, col = ink[["primary"]]))
  grid.text(expression(z[j] * "  (column " * italic(j) * ")"), rc, cb - 0.27, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 10.5, col = ink[["primary"]]))
  ar <- arrow(length = unit(0.055, "inches"), type = "closed")
  grid.lines(c(mx0 - 0.22, lc, lc), c(rowy, rowy, cb + chh + 0.08), default.units = "native",
             arrow = ar, gp = gpar(col = ink[["muted"]], lwd = 1.2, lty = 3, fill = ink[["muted"]]))
  seg(colx, my0 - 0.05, rc, cb + chh + 0.08, ink[["muted"]], 1.2, lty = 3, arrow = ar)
  txt("z-scores clipped at 0, then combined:", x0 + W2 / 2, 3.2, 9, ink[["secondary"]], face = "italic")
  # the two combinations (no statement about where the Stouffer combination came from)
  bxw <- W2 - 0.3; bxx <- x0 + 0.15
  grid.roundrect(bxx, 2.62, bxw, 0.46, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.06, "inches"), gp = gpar(fill = "#f1f0ec", col = ink[["muted"]], lwd = 1.1))
  txt("Euclidean combination", bxx + 0.1, 2.93, 10, face = "bold", just = "left")
  txt("(2007 packages we found)", bxx + 0.1, 2.74, 9, ink[["secondary"]], just = "left")
  grid.text(expression(sqrt(z[i]^2 + z[j]^2)), bxx + bxw - 0.1, 2.85, default.units = "native", just = "right",
            gp = gpar(fontfamily = FONT, fontsize = 12, col = ink[["primary"]]))
  grid.roundrect(bxx, 2.07, bxw, 0.46, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.06, "inches"), gp = gpar(fill = "#fbe3d9", col = pal[["orange"]], lwd = 1.4))
  txt("Stouffer combination", bxx + 0.1, 2.38, 10, face = "bold", just = "left")
  txt("(used here)", bxx + 0.1, 2.19, 9, ink[["secondary"]], just = "left")
  grid.text(expression((z[i] + z[j]) / sqrt(2)), bxx + bxw - 0.1, 2.3, default.units = "native", just = "right",
            gp = gpar(fontfamily = FONT, fontsize = 12, col = ink[["primary"]]))
  txt("= the CLR score for the pair i-j", x0 + W2 / 2, 1.88, 9.5, ink[["secondary"]], face = "italic")

  # ---------------- panel 3: softmax over the top 50 ---------------------------
  x0 <- X[3]; cx <- x0 + Wv[3] / 2; cy <- 5.02
  nb <- ring(20, cx, cy, 0.6, phase = 0.1)
  top <- c(1, 3, 5, 6, 9, 11, 13, 15, 17, 19)
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.13, 0.06)
    if (k %in% top) seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 2.6,
                        arrow = arrow(length = unit(0.055, "inches"), type = "closed"))
    else seg(s$x0, s$y0, s$x1, s$y1, ink[["muted"]], 0.8, alpha = 0.45)
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.05,
                                    col = if (k %in% top) ink[["secondary"]] else ink[["muted"]],
                                    alpha = if (k %in% top) 1 else 0.7)
  query(cx, cy, 0.12)
  txt("gene i attends over its own top 50\npositive CLR scores (20 drawn)", cx, 4.2, 9.5, ink[["secondary"]])
  # ranked scores and weights
  k50 <- 1:50; sc50 <- 2 + 4 * exp(-(k50 - 1) / 11); w50 <- sm_weights(sc50)
  px0 <- x0 + 0.4; px1 <- x0 + Wv[3] - 0.2; pb <- 3.0; ph <- 0.85; bw <- (px1 - px0) / 50
  for (k in k50) grid.rect(px0 + (k - 1) * bw, pb, bw * 0.8, w50[k] / max(w50) * ph, default.units = "native",
                           just = c("left", "bottom"), gp = gpar(fill = pal[["violet"]], col = NA))
  grid.lines(px0 + (k50 - 0.4) * bw, pb + (sc50 - 1.5) / (6 - 1.5) * ph, default.units = "native",
             gp = gpar(col = ink[["muted"]], lwd = 1.6, lty = 1))
  seg(px0 - 0.03, pb, px1, pb, ink[["secondary"]], 1)
  for (k in c(1, 10, 50)) { seg(px0 + (k - 0.6) * bw, pb, px0 + (k - 0.6) * bw, pb - 0.05, ink[["secondary"]], 1)
    txt(as.character(k), px0 + (k - 0.6) * bw, pb - 0.16, 9, ink[["secondary"]]) }
  txt("partner rank by CLR score", (px0 + px1) / 2, pb - 0.33, 9, ink[["secondary"]])
  txt("attention weight", px0 + 1.45, pb + ph * 0.78, 9, pal[["violet"]], "bold", just = "left")
  txt("CLR score", px0 + 1.45, pb + ph * 0.55, 9, ink[["secondary"]], just = "left")
  txt("Weights fall with rank, so partners are\ncontiguous in CLR rank. Up to 50 partners (fewer\nif fewer positive scores; 10 or fewer: equal\nweights). Temperature T_i is set per gene so\nthat ~10 partners count.", cx, 2.2, 8.5, ink[["secondary"]], lineheight = 1.0)

  # ---------------- panel 4: multi-step walk ----------------------------------
  x0 <- X[4]
  ix <- x0 + 0.4; iy <- 4.8
  l1x <- x0 + 1.4; l1y <- c(5.45, 4.8, 4.15)
  l2x <- x0 + 2.45; l2y <- c(5.7, 5.35, 5.0, 4.65, 4.3, 3.95)
  edges2 <- list(c(1, 1), c(1, 2), c(1, 3), c(2, 3), c(2, 4), c(2, 5), c(3, 5), c(3, 6))
  for (e in edges2) { s <- trim(l1x, l1y[e[1]], l2x, l2y[e[2]], 0.07, 0.05)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 1.2, alpha = 0.7,
        arrow = arrow(length = unit(0.045, "inches"), type = "closed")) }
  for (k in 1:3) { s <- trim(ix, iy, l1x, l1y[k], 0.13, 0.07)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], c(3.2, 2.4, 1.6)[k],
        arrow = arrow(length = unit(0.06, "inches"), type = "closed")) }
  for (k in 1:3) node(l1x, l1y[k], 0.07, fill = "#f6c9dc", col = pal[["magenta"]], lwd = 1.3)
  for (k in 1:6) node(l2x, l2y[k], 0.05, fill = "#fbe6ef", col = pal[["magenta"]], lwd = 1)
  query(ix, iy, 0.12)
  txt("step 1", l1x, 5.85, 9.5, ink[["secondary"]], "bold"); txt("step 2", l2x, 6.0, 9.5, ink[["secondary"]], "bold")
  txt("gene i", ix, 5.35, 9.5, ink[["secondary"]], "bold")
  txt("Each step moves a gene's mass to its weighted\npartners, then on to those partners' own partners\n(rows of the attention matrix A). A is not\nsymmetric. Half the mass stays at each step.",
      x0 + Wv[4] / 2, 3.3, 8.5, ink[["secondary"]], lineheight = 1.0)
  sx0 <- x0 + 0.4; sx1 <- x0 + Wv[4] - 0.4; sy <- 2.3
  tx <- sx0 + (0:3) / 3 * (sx1 - sx0)
  txt("steps t (log scale)", (sx0 + sx1) / 2, sy + 0.22, 9.5, ink[["secondary"]])
  seg(sx0, sy, sx1, sy, pal[["magenta"]], 4)
  for (k in 1:4) {
    seg(tx[k], sy - 0.08, tx[k], sy + 0.08, ink[["secondary"]], 1.2)
    txt(c("1", "3", "9", "27")[k], tx[k], sy - 0.2, 10, ink[["secondary"]], face = if (k == 3) "bold" else "plain")
  }
  grid.circle(tx[3], sy, 0.065, default.units = "native", gp = gpar(fill = ink[["primary"]], col = "white"))
  txt("t = 9, chosen on held-out labs", tx[3] + 0.1, sy - 0.4, 9, ink[["secondary"]])

  # ---------------- lead-method strip, under panels 3-4 ----------------------------------------
  txt("Published CLR stops here: a score for each gene pair", (X[1] + X[2] + Wv[2]) / 2, 1.35,
      11.5, ink[["secondary"]], face = "italic")
  sx <- X[3]; sw <- xr - X[3]; sy0 <- 0.9; sh <- 0.62
  grid.roundrect(sx, sy0, sw, sh, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "#f6fcf9", col = NA))
  newbox(sx, sy0, sw, sh)
  txt(paste0(LEAD_NAME, " adds two steps after the walk:"), sx + 0.15, sy0 + sh - 0.17, 11, face = "bold", just = "left")
  pill("log of the walk P^9", sx + 0.95, sy0 + 0.2, 1.6, 0.27, fill = pal[["magenta"]], col = NA, size = 9.5, tcol = "white")
  grid.polygon(c(sx + 1.85, sx + 1.98, sx + 1.85), c(sy0 + 0.1, sy0 + 0.2, sy0 + 0.3), default.units = "native",
               gp = gpar(fill = ink[["muted"]], col = NA))
  pill("second CLR: row z over reachable pairs, clip at 0, Stouffer", sx + 4.1, sy0 + 0.2, 3.95, 0.27, fill = pal[["orange"]], col = NA, size = 9.5, tcol = "white")
  txt("no learned parameters", sx + sw - 0.15, sy0 + sh - 0.17, 10, ink[["secondary"]], "italic", just = "right")
  # transformer analogy strip spanning both bands
  grid.roundrect(X[1], 0.3, xr - X[1], 0.46, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "#f3f3f1", col = ink[["grid"]]))
  lab <- "Transformer analogy:"
  txt(lab, X[1] + 0.18, 0.53, 11.5, face = "bold", just = "left")
  lw <- convertWidth(grobWidth(textGrob(lab, gp = gpar(fontfamily = FONT, fontsize = 11.5,
                                                         fontface = "bold"))), "native", TRUE)
  txt("MI = query-key score;  CLR z-score = scaled logit;  softmax = softmax;  walk = repeated attention with residual",
      X[1] + 0.18 + lw + 0.1, 0.53, 11, ink[["secondary"]], just = "left")

  # ---------------- panel 5: results -------------------------------------------
  x5 <- xr + 0.2; W5 <- 15.75 - x5
  grid.roundrect(x5, 0.3, W5, 7.54, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.1, "inches"), gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.4))
  grid.circle(x5 + 0.27, 7.55, 0.16, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
  txt("5", x5 + 0.27, 7.545, 11.5, "white", "bold")
  txt("What it recovers", x5 + 0.53, 7.55, 14, face = "bold", just = "left")

  ix <- x5 + W5 - 0.45; iy <- 6.9
  tg <- ring(9, ix, iy, 0.26, phase = 0.2)
  for (k in 1:9) seg(ix, iy, tg$x[k], tg$y[k], pal[["aqua"]], 1.2)
  for (k in 1:9) node(tg$x[k], tg$y[k], 0.035, fill = "#d1efe4", col = pal[["aqua"]], lwd = 1)
  node(ix, iy, 0.08, fill = pal[["aqua"]], col = "#0d6e4c", lwd = 1)
  txt("E. coli regulons (experimentally\nsupported, RegulonDB)", x5 + 0.25, 6.95, 12,
      face = "bold", just = "left", lineheight = 1.0)
  txt("regulon = genes controlled by one\ntranscription factor (TF)", x5 + 0.25, 6.45, 9.5,
      ink[["secondary"]], just = "left")
  bar_strip(strip_a, x5, W5, 5.95, 2, c(0, 0.5, 1, 1.5, 2))
  txt(paste0("per regulon: ", ACLR_NAME, " ties CLR and GENIE3"),
      x5 + 0.25, 3.8, 8.5, ink[["secondary"]], face = "italic", just = "left")

  # (b) complexes: human data, new in this work
  cxi <- x5 + W5 - 0.45; cyi <- 3.3
  cl <- ring(5, cxi, cyi, 0.2, phase = pi / 2)
  for (a in 1:4) for (b in (a + 1):5) seg(cl$x[a], cl$y[a], cl$x[b], cl$y[b], ink[["muted"]], 1.3)
  for (k in 1:5) node(cl$x[k], cl$y[k], 0.045, fill = "#e9e8e4", col = ink[["muted"]], lwd = 1)
  txt(paste0("Human K562 protein complexes\n(CORUM; mostly ribosome): ", ACLR_NAME, "\nscores 1.8x lower than CLR"), x5 + 0.25, 3.28,
      11, face = "bold", just = "left", lineheight = 1.0)
  txt("complex = proteins that bind as one machine", x5 + 0.25, 2.93, 9, ink[["secondary"]], just = "left")
  bar_strip(strip_b, x5, W5, 2.38, 9, c(0, 2, 4, 6, 8), loss = TRUE)
  newbox(x5 + 0.1, 0.38, W5 - 0.2, 3.27, 0.08)

  popViewport()
}

# ---- render -------------------------------------------------------------------
png(file.path(out_dir, "graphical_abstract_short_p1.png"), width = 2400, height = 1350, res = 150,
    type = "cairo", bg = "white")
draw_abstract(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_short_p1.pdf"), width = 16, height = 9, bg = "white")
draw_abstract(); invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract_short_p1.{png,pdf}"), "\n")
