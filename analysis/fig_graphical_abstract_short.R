# Graphical abstract for the SHORT pre-read, page 1 (portrait, 8.3 in wide canvas, printed 7.0 in wide), drawn entirely in code with grid.
# Copy of analysis/fig_graphical_abstract.R, revised 15 Nov (PI design):
#  - signage: a DOTTED GREEN outline marks every element that is new in this work, including every
#    element using human data (the original CLR paper evaluated no human data); historical CLR parts are
#    labelled "Faith et al. 2007". Legend entry: "dotted green = new in this work".
#  - panel 2 (PI wording, 19 Nov): the 2007 paper (Faith et al.) describes only the Euclidean combination sqrt(zi^2 + zj^2);
#    the released 2007 software offered both the Euclidean and the Stouffer combination (zi + zj)/sqrt(2); the Stouffer
#    combination is the one used here. NO statement about who introduced Stouffer or when, and it is not marked as new.
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

out_dir <- Sys.getenv("SB_OUT", "docs/whitepapers/figures"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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
# ---- tcaCLR2 (exploratory): same metric and base-rate scaling as the other bars (raw pooled AUPR / base rate), from results/tcaclr2_scoreboard/cells_tcaclr2_all.csv ----
T2_NAME <- "tcaCLR2"
t2c <- utils::read.csv(Sys.getenv("SB_T2", "results/tcaclr2_scoreboard/cells_tcaclr2_all.csv"), stringsAsFactors = FALSE)
t2e <- function(ds, truth) { x <- t2c[t2c$dataset == ds & t2c$truth == truth & t2c$metric == "pooled_aupr_over_base" & t2c$scope == "raw" & t2c$method == "tcaCLR2", ]; stopifnot(nrow(x) == 1L); x$enrichment_over_random }
T2_M3D <- t2e("m3d", "SC"); T2_K562 <- t2e("k562", "CORUM")
strip_a <- data.frame(   # M3D, strong-evidence (SC) RegulonDB regulons
  label = c("CLR (historical)", ACLR_NAME, LEAD_NAME, T2_NAME, "GENIE3", "PIDC"),
  sub = c("", "adopted", "lead candidate", "exploratory", "", ""),
  value = c(1.44, 1.68, 1.72, T2_M3D, 1.51, 1.45),
  col = c(pal[["orange"]], pal[["aqua"]], EXPL, EXPL, ink[["muted"]], ink[["muted"]]),
  new = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
strip_b <- data.frame(   # K562, CORUM complexes
  label = c("CLR (historical)", ACLR_NAME, LEAD_NAME, T2_NAME, "GENIE3", "|r| (correlation)"),
  sub = c("", "adopted", "lead; picked on K562", "exploratory", "", ""),
  value = c(6.39, 3.51, 7.95, T2_K562, 6.32, 7.00),
  col = c(pal[["orange"]], pal[["aqua"]], EXPL, EXPL, ink[["muted"]], ink[["muted"]]),
  new = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)

# ---- PORTRAIT layout: canvas W_CAN x H inches (H from the layout below), printed 7.0 in wide: scale 7.0 / W_CAN = 0.843;
# every font size here is >= 8.4 pt, i.e. >= 7.1 pt printed. Panels 1-3 sit in one row, panel 4 (the walk) is a full-width row,
# the lead-method strip and the transformer analogy follow, and panel 5 (what it recovers) closes the figure.
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
band <- function(x0, x1, ytop, fill, label) {   # labelled band + bracket over a group of panels; (ytop = top edge)
  grid.roundrect(x0, ytop - 0.36, x1 - x0, 0.36, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.07, "inches"), gp = gpar(fill = fill, col = NA))
  txt(label, x0 + 0.14, ytop - 0.18, 11.5, face = "bold", just = "left")
  by <- ytop - 0.45
  grid.lines(c(x0 + 0.04, x0 + 0.04, x1 - 0.04, x1 - 0.04), c(by - 0.07, by, by, by - 0.07),
             default.units = "native", gp = gpar(col = ink[["secondary"]], lwd = 1.1))
}
# dotted green outline (new in this work); (x, yt) = top-left corner
newbox <- function(x, yt, w, h, r = 0.07)
  grid.roundrect(x, yt - h, w, h, default.units = "native", just = c("left", "bottom"),
                 r = unit(r, "inches"), gp = DOT)

# bar strip: rows (label, optional status text on the same line); dotted green around the rows that are new; yt = first row centre
bar_strip <- function(d, x0, w, y_first, xmax, ticks, loss = FALSE, pitch = 0.19) {
  lab_x <- x0 + 0.1; bx0 <- x0 + 2.0; bx1 <- x0 + w - 0.4
  sc <- (bx1 - bx0) / xmax
  n <- nrow(d); ry <- y_first - (seq_len(n) - 1) * pitch; bh <- 0.14
  for (t in ticks) seg(bx0 + t * sc, ry[1] + 0.15, bx0 + t * sc, ry[n] - 0.15, ink[["grid"]], 0.8)
  for (k in seq_len(n)) {
    has_sub <- nzchar(d$sub[k])
    if (d$new[k]) grid.roundrect(x0 + 0.03, ry[k], w - 0.06, pitch - 0.02, default.units = "native",
                                 just = c("left", "centre"), r = unit(0.04, "inches"), gp = gpar(col = GREEN, lty = "11", lwd = 2.0, fill = NA))
    txt(d$label[k], lab_x, ry[k], 9.3, face = if (d$new[k]) "bold" else "plain", just = "left")
    if (has_sub) { lw <- text_w(d$label[k], 9.3, if (d$new[k]) "bold" else "plain")
      txt(d$sub[k], lab_x + lw + 0.06, ry[k], FS_S, ink[["secondary"]], "italic", just = "left") }
    lead_row <- d$col[k] == EXPL
    grid.rect(bx0, ry[k], d$value[k] * sc, bh, default.units = "native", just = c("left", "centre"),
              gp = gpar(fill = d$col[k], col = if (lead_row) pal[["aqua"]] else NA, lwd = 1.0,
                        lty = if (lead_row) 2 else 1))
    if (loss && k == 2) txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc - 0.04, ry[k], 9, "white", "bold", just = "right")
    else txt(sprintf("%.2f", d$value[k]), bx0 + d$value[k] * sc + 0.04, ry[k], 9.3,
        face = if (d$new[k]) "bold" else "plain", just = "left")
  }
  seg(bx0 + sc, ry[1] + 0.2, bx0 + sc, ry[n] - 0.15, ink[["primary"]], 1.0, lty = 2)
  txt("random = 1", bx0 + sc + 0.04, ry[1] + 0.22, FS_S, ink[["secondary"]], just = "left")
  seg(bx0, ry[1] + 0.15, bx0, ry[n] - 0.15, ink[["secondary"]], 0.9)
  ay <- ry[n] - 0.15
  seg(bx0, ay, bx1, ay, ink[["secondary"]], 0.9)
  for (t in ticks) { seg(bx0 + t * sc, ay, bx0 + t * sc, ay - 0.04, ink[["secondary"]], 0.9)
    txt(as.character(t), bx0 + t * sc, ay - 0.13, FS_S, ink[["secondary"]]) }
  txt("precision-recall area, multiple of random", (bx0 + bx1) / 2 - 0.3, ay - 0.29, FS_S, ink[["secondary"]])
  if (loss) {   # red loss marker: from the CLR value back to CLR attention's value (row 2)
    ly <- ry[2]
    seg(bx0 + d$value[1] * sc, ly, bx0 + d$value[2] * sc + 0.04, ly, pal[["red"]], 2.2,
        arrow = arrow(length = unit(0.07, "inches"), type = "closed"))
    seg(bx0 + d$value[1] * sc, ly - 0.07, bx0 + d$value[1] * sc, ly + 0.07, pal[["red"]], 2.2)
  }
  ay - 0.29      # y of the axis title
}

# schematic weights for panel 3: made-up CLR scores for one gene's 50 candidates, then the method's own rule
# (softmax(score / T), T by bisection so that exp(entropy) = 10)
sm_weights <- function(s, keff = 10) {
  f <- function(lt) { w <- exp((s - max(s)) / exp(lt)); w <- w / sum(w); exp(-sum(w * log(w))) }
  lo <- -6; hi <- 6
  for (it in 1:80) { mid <- (lo + hi) / 2; if (f(mid) < keff) lo <- mid else hi <- mid }
  w <- exp((s - max(s)) / exp((lo + hi) / 2)); w / sum(w)
}

# ---- vertical layout (distance from the top of the canvas) --------------------------------------------------------------------
D_BAND <- 1.12; D_R1 <- 1.7; H1 <- 3.94          # band top, row 1 top, row 1 height
D_NOTE <- D_R1 + H1 + 0.12
D_R2 <- D_NOTE + 0.15; H2 <- 1.45                 # row 2 (the walk)
D_STR <- D_R2 + H2 + 0.07; HSTR <- 0.56           # lead-method strip
D_ANA <- D_STR + HSTR + 0.07; HANA <- 0.3        # transformer analogy
D_R3 <- D_ANA + HANA + 0.09; H3 <- 3.2          # results panel
H_CAN <- D_R3 + H3 + 0.1

# ---- drawing ----------------------------------------------------------------
draw_abstract <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  H <- H_CAN
  pushViewport(viewport(xscale = c(0, W_CAN), yscale = c(0, H)))
  T <- H
  MG <- 0.12; XR <- W_CAN - 0.12; TW <- XR - MG

  grid.text(expression(bold("Adding a softmax and a multi-step walk to CLR recovers large ") *
                         bolditalic("E. coli") * bold(" regulons")),
            MG, T - 0.27, default.units = "native", just = "left",
            gp = gpar(fontfamily = FONT, fontsize = 14.5, col = ink[["primary"]]))
  txt("better (16% over CLR), but does poorly on tight correlation blocks",
      MG, T - 0.51, 14.5, face = "bold", just = "left")
  txt(paste0(LEAD_NAME, " recovers the K562 complex loss (a selected result; see text). ", NAMES_NOTE), MG, T - 0.77, 10.5,
      ink[["secondary"]], "bold.italic", just = "left")
  # legend (signage convention)
  ly <- T - 1.0
  grid.roundrect(MG, ly, 0.34, 0.16, default.units = "native", just = c("left", "centre"),
                 r = unit(0.03, "inches"), gp = DOT)
  txt("dotted green = new in this work", MG + 0.43, ly, 9.5, face = "bold", just = "left")
  lx2 <- 3.35
  grid.rect(lx2, ly, 0.34, 0.11, default.units = "native", just = c("left", "centre"),
            gp = gpar(fill = EXPL, col = pal[["aqua"]], lwd = 1.0, lty = 2))
  txt(paste0("dashed bars = ", LEAD_NAME, " (", LEAD_STATUS, "), ", T2_NAME, " (exploratory)"), lx2 + 0.43, ly, FS_S, ink[["secondary"]], just = "left")

  # ---- row 1: panels 1-3 ----
  gap <- 0.12; Wv <- c(1.8, 3.0, 3.0)
  X <- MG + c(0, cumsum(Wv[1:2] + gap))
  yt <- T - D_R1; ybm <- yt - H1
  band(X[1], X[2] + Wv[2], T - D_BAND, "#ecebe7", "CLR: Faith et al. 2007 (historical)")
  band(X[3], X[3] + Wv[3], T - D_BAND, "#e6e3fa", paste0(ACLR_NAME, ": new in this work"))
  newbox(X[3], T - D_BAND, Wv[3], 0.36)
  card(X[1], yt, Wv[1], H1, 1, "Expression\ncompendium")
  card(X[2], yt, Wv[2], H1, 2, "Mutual information,\nz-scored against each\ngene's own background")
  card(X[3], yt, Wv[3], H1, 3, "Softmax with a\nper-gene temperature\nover the top 50 partners")
  newbox(X[3], yt, Wv[3], H1, 0.08)
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
  # bacterium (E. coli) and human cell (K562; human data is new in this work)
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
  newbox(x0 + 1.0, hy + 0.32, 0.74, 1.3, 0.05)      # all human data is new in this work

  # ---------------- panel 2: CLR criss-cross (after Faith et al. 2007, Fig. 2A) ----
  x0 <- X[2]; W2 <- Wv[2]
  pill("CLR as originally published", x0 + W2 / 2, yt - 0.74, 2.3, 0.22, fill = "#ecebe7", col = NA, size = 9.3,
       face = "plain", tcol = ink[["secondary"]])
  set.seed(5)
  n <- 8; cs <- 0.1; mx0 <- x0 + 0.95; my1 <- yt - 1.1; my0 <- my1 - n * cs; mx1 <- mx0 + n * cs
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
  cb <- yt - 2.4; chh <- 0.28
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
  grid.text(expression(z[i] * "  (row " * italic(i) * ")"), lc, cb - 0.22, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 9.3, col = ink[["primary"]]))
  grid.text(expression(z[j] * "  (column " * italic(j) * ")"), rc, cb - 0.22, default.units = "native",
            gp = gpar(fontfamily = FONT, fontsize = 9.3, col = ink[["primary"]]))
  ar <- arrow(length = unit(0.05, "inches"), type = "closed")
  grid.lines(c(mx0 - 0.2, lc, lc), c(rowy, rowy, cb + chh + 0.07), default.units = "native",
             arrow = ar, gp = gpar(col = ink[["muted"]], lwd = 1.1, lty = 3, fill = ink[["muted"]]))
  seg(colx, my0 - 0.04, rc, cb + chh + 0.07, ink[["muted"]], 1.1, lty = 3, arrow = ar)
  txt("z-scores clipped at 0, then combined:", x0 + W2 / 2, yt - 2.88, FS_S, ink[["secondary"]], face = "italic")
  # the two combinations (no statement about who introduced the Stouffer combination or when)
  bxw <- W2 - 0.3; bxx <- x0 + 0.15
  grid.roundrect(bxx, yt - 3.4, bxw, 0.4, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.05, "inches"), gp = gpar(fill = "#f1f0ec", col = ink[["muted"]], lwd = 1.0))
  txt("Euclidean combination", bxx + 0.09, yt - 3.1, 9.3, face = "bold", just = "left")
  txt("(2007 paper and 2007 software)", bxx + 0.09, yt - 3.28, FS_S, ink[["secondary"]], just = "left")
  grid.text(expression(sqrt(z[i]^2 + z[j]^2)), bxx + bxw - 0.09, yt - 3.19, default.units = "native", just = "right",
            gp = gpar(fontfamily = FONT, fontsize = 10.5, col = ink[["primary"]]))
  grid.roundrect(bxx, yt - 3.86, bxw, 0.4, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.05, "inches"), gp = gpar(fill = "#fbe3d9", col = pal[["orange"]], lwd = 1.3))
  txt("Stouffer combination", bxx + 0.09, yt - 3.56, 9.3, face = "bold", just = "left")
  txt("(2007 software; used here)", bxx + 0.09, yt - 3.74, FS_S, ink[["secondary"]], just = "left")
  grid.text(expression((z[i] + z[j]) / sqrt(2)), bxx + bxw - 0.09, yt - 3.65, default.units = "native", just = "right",
            gp = gpar(fontfamily = FONT, fontsize = 10.5, col = ink[["primary"]]))

  # ---------------- panel 3: softmax over the top 50 ---------------------------
  x0 <- X[3]; cx <- x0 + Wv[3] / 2; cy <- yt - 1.3
  nb <- ring(20, cx, cy, 0.5, phase = 0.1)
  top <- c(1, 3, 5, 6, 9, 11, 13, 15, 17, 19)
  for (k in seq_len(nrow(nb))) {
    s <- trim(cx, cy, nb$x[k], nb$y[k], 0.12, 0.055)
    if (k %in% top) seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 2.3,
                        arrow = arrow(length = unit(0.05, "inches"), type = "closed"))
    else seg(s$x0, s$y0, s$x1, s$y1, ink[["muted"]], 0.7, alpha = 0.45)
  }
  for (k in seq_len(nrow(nb))) node(nb$x[k], nb$y[k], 0.045,
                                    col = if (k %in% top) ink[["secondary"]] else ink[["muted"]],
                                    alpha = if (k %in% top) 1 else 0.7)
  query(cx, cy, 0.11)
  txt("gene i attends over its own top 50\npositive CLR scores (20 drawn)", cx, yt - 2.0, 8.8, ink[["secondary"]])
  # ranked scores and weights
  k50 <- 1:50; sc50 <- 2 + 4 * exp(-(k50 - 1) / 11); w50 <- sm_weights(sc50)
  px0 <- x0 + 0.35; px1 <- x0 + Wv[3] - 0.2; pb <- yt - 2.95; ph <- 0.6; bw <- (px1 - px0) / 50
  for (k in k50) grid.rect(px0 + (k - 1) * bw, pb, bw * 0.8, w50[k] / max(w50) * ph, default.units = "native",
                           just = c("left", "bottom"), gp = gpar(fill = pal[["violet"]], col = NA))
  grid.lines(px0 + (k50 - 0.4) * bw, pb + (sc50 - 1.5) / (6 - 1.5) * ph, default.units = "native",
             gp = gpar(col = ink[["muted"]], lwd = 1.4, lty = 1))
  seg(px0 - 0.03, pb, px1, pb, ink[["secondary"]], 0.9)
  for (k in c(1, 10, 50)) { seg(px0 + (k - 0.6) * bw, pb, px0 + (k - 0.6) * bw, pb - 0.04, ink[["secondary"]], 0.9)
    txt(as.character(k), px0 + (k - 0.6) * bw, pb - 0.13, FS_S, ink[["secondary"]]) }
  txt("partner rank by CLR score", (px0 + px1) / 2, pb - 0.28, FS_S, ink[["secondary"]])
  txt("attention weight", px0 + 1.25, pb + ph * 0.8, FS_S, pal[["violet"]], "bold", just = "left")
  txt("CLR score", px0 + 1.25, pb + ph * 0.5, FS_S, ink[["secondary"]], just = "left")
  wl <- wrap_to("Weights fall with rank, so partners are contiguous in CLR rank. Up to 50 partners (fewer if fewer positive scores; 10 or fewer: equal weights). Temperature T_i is set per gene so that ~10 partners count.", FS_S, 2.75)
  for (k in seq_along(wl)) txt(wl[k], cx, yt - 3.4 - (k - 1) * 0.125, FS_S, ink[["secondary"]])

  # "Published CLR stops here" under panels 1-2
  txt("Published CLR stops here: a score for each gene pair", (X[1] + X[2] + Wv[2]) / 2, T - D_NOTE, 9.8, ink[["secondary"]], face = "italic")

  # ---------------- row 2, panel 4: multi-step walk (full width) ---------------------------
  yq <- T - D_R2; x0 <- MG; W4 <- TW
  card(x0, yq, W4, H2, 4, "Multi-step walk;\nsteps chosen\nwithout labels")
  newbox(x0, yq, W4, H2, 0.08)
  ix <- x0 + 2.05; iy <- yq - 0.88
  l1x <- x0 + 2.95; l1y <- iy + c(0.3, 0, -0.3)
  l2x <- x0 + 3.95; l2y <- iy + c(0.45, 0.27, 0.09, -0.09, -0.27, -0.45)
  edges2 <- list(c(1, 1), c(1, 2), c(1, 3), c(2, 3), c(2, 4), c(2, 5), c(3, 5), c(3, 6))
  for (e in edges2) { s <- trim(l1x, l1y[e[1]], l2x, l2y[e[2]], 0.07, 0.05)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], 1.1, alpha = 0.7,
        arrow = arrow(length = unit(0.045, "inches"), type = "closed")) }
  for (k in 1:3) { s <- trim(ix, iy, l1x, l1y[k], 0.12, 0.07)
    seg(s$x0, s$y0, s$x1, s$y1, pal[["violet"]], c(3.0, 2.2, 1.5)[k],
        arrow = arrow(length = unit(0.055, "inches"), type = "closed")) }
  for (k in 1:3) node(l1x, l1y[k], 0.065, fill = "#f6c9dc", col = pal[["magenta"]], lwd = 1.2)
  for (k in 1:6) node(l2x, l2y[k], 0.045, fill = "#fbe6ef", col = pal[["magenta"]], lwd = 0.9)
  query(ix, iy, 0.11)
  txt("step 1", l1x, yq - 0.22, FS_S, ink[["secondary"]], "bold"); txt("step 2", l2x, yq - 0.12, FS_S, ink[["secondary"]], "bold")
  txt("gene i", ix, iy + 0.25, FS_S, ink[["secondary"]], "bold")
  wl <- wrap_to("Each step moves a gene's mass to its weighted partners, then on to those partners' own partners (rows of the attention matrix A). A is not symmetric. Half the mass stays at each step.", FS_S, 2.2)
  for (k in seq_along(wl)) txt(wl[k], x0 + 4.45, iy + ((length(wl) - 1) / 2 - (k - 1)) * 0.125, FS_S, ink[["secondary"]], just = "left")
  sx0 <- x0 + 6.75; sx1 <- x0 + W4 - 0.3; sy <- yq - 0.85
  tx <- sx0 + (0:3) / 3 * (sx1 - sx0)
  txt("steps t (log scale)", (sx0 + sx1) / 2, sy + 0.28, FS_S, ink[["secondary"]])
  seg(sx0, sy, sx1, sy, pal[["magenta"]], 3.5)
  for (k in 1:4) {
    seg(tx[k], sy - 0.07, tx[k], sy + 0.07, ink[["secondary"]], 1.1)
    txt(c("1", "3", "9", "27")[k], tx[k], sy - 0.19, 9, ink[["secondary"]], face = if (k == 3) "bold" else "plain")
  }
  grid.circle(tx[3], sy, 0.055, default.units = "native", gp = gpar(fill = ink[["primary"]], col = "white"))
  txt("t = 9, chosen on\nheld-out labs", tx[3], sy - 0.45, FS_S, ink[["secondary"]])

  # ---------------- lead-method strip -------------------------------------------
  ys <- T - D_STR
  grid.roundrect(MG, ys - HSTR, TW, HSTR, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.07, "inches"), gp = gpar(fill = "#f6fcf9", col = NA))
  newbox(MG, ys, TW, HSTR)
  txt(paste0(LEAD_NAME, " adds two steps after the walk:"), MG + 0.14, ys - 0.15, 10, face = "bold", just = "left")
  txt("no learned parameters", MG + TW - 0.14, ys - 0.15, 9.3, ink[["secondary"]], "italic", just = "right")
  pill("log of the walk P^9", MG + 0.95, ys - 0.43, 1.6, 0.25, fill = pal[["magenta"]], col = NA, size = 9, tcol = "white")
  grid.polygon(c(MG + 1.9, MG + 2.02, MG + 1.9), c(ys - 0.35, ys - 0.43, ys - 0.51), default.units = "native",
               gp = gpar(fill = ink[["muted"]], col = NA))
  pill("second CLR: row z over reachable pairs, clip at 0, Stouffer", MG + 4.2, ys - 0.43, 4.0, 0.25, fill = pal[["orange"]], col = NA, size = 9, tcol = "white")
  # transformer analogy strip
  ya <- T - D_ANA
  grid.roundrect(MG, ya - HANA, TW, HANA, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.07, "inches"), gp = gpar(fill = "#f3f3f1", col = ink[["grid"]]))
  lab <- "Transformer analogy:"
  txt(lab, MG + 0.14, ya - HANA / 2, 9.5, face = "bold", just = "left")
  lw <- text_w(lab, 9.5, "bold")
  txt("MI = query-key score;  CLR z-score = scaled logit;  softmax = softmax;  walk = repeated attention with residual",
      MG + 0.14 + lw + 0.08, ya - HANA / 2, 9, ink[["secondary"]], just = "left")

  # ---------------- panel 5: results (full width, two blocks side by side) -----------------------------------------
  y5 <- T - D_R3; x5 <- MG; W5 <- TW
  grid.roundrect(x5, y5 - H3, W5, H3, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "white", col = ink[["grid"]], lwd = 1.2))
  grid.circle(x5 + 0.21, y5 - 0.23, 0.13, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
  txt("5", x5 + 0.21, y5 - 0.235, 9.5, "white", "bold")
  txt("What it recovers", x5 + 0.43, y5 - 0.23, 11.5, face = "bold", just = "left")
  d0 <- y5 - 0.42
  # (a) regulons
  bxa <- x5 + 0.1; wa <- 3.85
  ix <- bxa + wa - 0.35; iy <- d0 - 0.3
  tg <- ring(9, ix, iy, 0.2, phase = 0.2)
  for (k in 1:9) seg(ix, iy, tg$x[k], tg$y[k], pal[["aqua"]], 1.0)
  for (k in 1:9) node(tg$x[k], tg$y[k], 0.03, fill = "#d1efe4", col = pal[["aqua"]], lwd = 0.9)
  node(ix, iy, 0.065, fill = pal[["aqua"]], col = "#0d6e4c", lwd = 0.9)
  txt("E. coli regulons (experimentally\nsupported, RegulonDB)", bxa + 0.1, d0 - 0.2, 10.5, face = "bold", just = "left", lineheight = 1.0)
  txt("regulon = genes controlled by one\ntranscription factor (TF)", bxa + 0.1, d0 - 0.58, FS_S, ink[["secondary"]], just = "left")
  ya_ax <- bar_strip(strip_a, bxa, wa, d0 - 1.0, 2, c(0, 0.5, 1, 1.5, 2))
  txt(paste0("per regulon: ", ACLR_NAME, " ties CLR and GENIE3"), bxa + 0.1, ya_ax - 0.2, FS_S, ink[["secondary"]], face = "italic", just = "left")
  # (b) complexes: human data, new in this work
  bxb <- x5 + 4.1; wb <- 3.85
  cxi <- bxb + wb - 0.35; cyi <- d0 - 0.3
  cl <- ring(5, cxi, cyi, 0.17, phase = pi / 2)
  for (a in 1:4) for (b in (a + 1):5) seg(cl$x[a], cl$y[a], cl$x[b], cl$y[b], ink[["muted"]], 1.1)
  for (k in 1:5) node(cl$x[k], cl$y[k], 0.035, fill = "#e9e8e4", col = ink[["muted"]], lwd = 0.9)
  txt(paste0("Human K562 protein complexes\n(CORUM; mostly ribosome): ", ACLR_NAME, "\nscores 1.8x lower than CLR"), bxb + 0.1, d0 - 0.26,
      10.5, face = "bold", just = "left", lineheight = 1.0)
  txt("complex = proteins that bind as one machine", bxb + 0.1, d0 - 0.66, FS_S, ink[["secondary"]], just = "left")
  yb_ax <- bar_strip(strip_b, bxb, wb, d0 - 1.0, 9, c(0, 2, 4, 6, 8), loss = TRUE)
  newbox(bxb - 0.05, d0 + 0.04, wb + 0.1, d0 + 0.04 - (yb_ax - 0.16), 0.07)

  popViewport()
}

# ---- render -------------------------------------------------------------------
png(tempfile(fileext = ".png"), width = 800, height = 800, type = "cairo")   # measuring device is open while text widths are taken
png(file.path(out_dir, "graphical_abstract_short_p1.png"), width = round(W_CAN * 300), height = round(H_CAN * 300), res = 300,
    type = "cairo", bg = "white")
draw_abstract(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_short_p1.pdf"), width = W_CAN, height = H_CAN, bg = "white")
draw_abstract(); invisible(dev.off())
invisible(dev.off())
cat(sprintf("p1: canvas %.2f x %.2f in; printed 7.00 x %.2f in (scale %.3f); min font %.1f pt -> %.1f pt printed\n", W_CAN, H_CAN, H_CAN * 7 / W_CAN, 7 / W_CAN, FS_S, FS_S * 7 / W_CAN))
cat("wrote", file.path(out_dir, "graphical_abstract_short_p1.{png,pdf}"), "\n")
