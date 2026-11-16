# Graphical abstract for the SHORT pre-read, page 2 (16:9): the versions compared, one row each, drawn in code with grid.
# Copy of analysis/fig_graphical_abstract_p2.R, revised 15 Nov: only CLR, CLR attention and the lead method, plus two
# of the original ablation rows (attention without CLR normalization). Dropped rows: CLR-NS, |r| (correlation),
# attention on CLR of |r|, attention on CLR-NS, and the stacked variants (a3r2, a3r3).
# Dotted green = new in this work (softmax, walk, log + second CLR; all human-data columns).
# Numbers, all copied from existing result files:
#   E. coli SC AUPR: Table 1 of docs/whitepapers/paper.md (466 experiments, strong-evidence regulons,
#     each operator at its own held-out depth t*; CLR = scores ranked directly)
#   lead method: E. coli SC AUPR 0.1622 and K562 CORUM 7.95 from results/explore_stack/{m3d,k562}/pooled_levels.csv;
#     K562 CORUM for CLR 6.39 and CLR attention 3.51 from docs/whitepapers/sota_table.md
#   BEELINE median AUPRC ratio (44 case-truth cells): CLR 1.48 (plain CLR); CLR attention 1.54 and the lead method 1.73
#     use the CLR-NS operator (lead at t = 9, not tuned), docs/whitepapers/figures/scoreboard_data.csv
# The lead method's name is defined once, in LEAD_NAME (a rename is a one-line change and a re-run).
#
# Output: docs/whitepapers/figures/graphical_abstract_short_p2.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_short_p2.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")

out_dir <- "docs/whitepapers/figures"
LEAD_NAME <- "tcaCLR"                     # DEFINED HERE ONLY (temporary name; internal name a3r)
ACLR_NAME <- "aCLR"                       # DEFINED HERE ONLY (temporary name)
NAMES_NOTE <- "Method names are temporary."
LEAD_STATUS <- "main lead, not yet adopted"
GREEN <- pal[["green"]]
DOT <- gpar(col = GREEN, lty = "11", lwd = 2.4, fill = NA)
FONT <- "Liberation Sans"

txt <- function(label, x, y, size = 12, col = ink[["primary"]], face = "plain",
                just = "centre", lineheight = 1.05)
  grid.text(label, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face,
                      lineheight = lineheight))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1)
  grid.segments(x0, y0, x1, y1, default.units = "native",
                gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round"))
pill <- function(label, x, y, w, h = 0.34, fill = "white", col = ink[["primary"]], size = 11,
                 face = "bold", tcol = ink[["primary"]], lwd = 1.3, lty = 1)
{
  grid.roundrect(x, y, w, h, default.units = "native", r = unit(h / 2, "inches"),
                 gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
  txt(label, x, y, size, tcol, face)
}

# component chips: filled = present in this version, dashed grey dash = absent
dep_cols <- c(MI = pal[["blue"]], r = ink[["muted"]], NS = ink[["secondary"]])
dep_lab  <- c(MI = "mutual information", r = "|r|", NS = "|r| of normal scores")
chip <- function(present, x, y, w, label, col, tcol = "white", new = FALSE) {
  if (present) {
    pill(label, x, y, w, fill = col, col = col, tcol = tcol, size = 10.5)
    if (new) grid.roundrect(x, y, w + 0.12, 0.34 + 0.1, default.units = "native", r = unit(0.22, "inches"), gp = DOT)
  } else seg(x - 0.14, y, x + 0.14, y, ink[["grid"]], 2)
}

# ---- the versions ----------------------------------------------------------------
v <- data.frame(
  name = c("CLR (historical)", "Attention on |r|", "Attention on mutual information",
           ACLR_NAME, LEAD_NAME),
  sub = c("Faith et al. 2007: context normalization, scores ranked directly",
          "softmax + walk, no CLR normalization",
          "softmax + walk, no CLR normalization",
          "MI, CLR normalization, softmax + walk (adopted)",
          paste0(ACLR_NAME, ", then log + second CLR (main lead)")),
  dep = c("MI", "r", "MI", "MI", "MI"),
  cn  = c(TRUE, FALSE, FALSE, TRUE, TRUE),
  sm  = c(FALSE, TRUE, TRUE, TRUE, TRUE),
  lc  = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  walklab = c(NA, "t* = 8", "t* = 6", "t* = 9", "t = 9"),
  ecoli = c(0.136, 0.143, 0.143, 0.158, 0.1622),
  corum = c("6.39", "", "", "3.51", "7.95"),
  beeline = c("1.48", "", "", "1.54", "1.73"),
  stringsAsFactors = FALSE)

draw_p2 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))

  txt(paste0(paste0("The versions compared: CLR, ", ACLR_NAME, " and "), LEAD_NAME),
      0.5, 8.55, 22, face = "bold", just = "left")
  txt("E. coli: strong-evidence regulons, 466 experiments (RegulonDB).  K562: CORUM complexes.  Single cells: BEELINE, median over 44 case-truth cells (14 dataset-and-gene-set cases).",
      0.5, 8.1, 11.5, ink[["secondary"]], just = "left")
  grid.roundrect(0.5, 7.78, 0.42, 0.2, default.units = "native", just = c("left", "centre"), r = unit(0.04, "inches"), gp = DOT)
  txt("dotted green = new in this work (the softmax, the walk, the log + second CLR, and all human data)", 1.02, 7.78, 10.5, ink[["secondary"]], just = "left")

  txt(NAMES_NOTE, 15.7, 7.78, 10.5, ink[["secondary"]], "italic", just = "right")
  xn <- 0.4; xd <- 4.95; xc <- 6.75; xs <- 8.0; xw <- 9.1; xl <- 10.6; xe0 <- 11.3; xe1 <- 12.9; xk <- 13.25; xb <- 14.5
  hy <- 7.2
  for (h in list(list("Version", xn, "left"), list("Dependence", xd, "centre"),
                 list("Context\nnormalization", xc, "centre"), list("Softmax", xs, "centre"),
                 list("Walk", xw, "centre"), list("Log +\nsecond CLR", xl, "centre"),
                 list("E. coli\nAUPR (SC)", (xe0 + xe1) / 2, "centre"),
                 list("K562 CORUM\n(x random)", xk, "left"),
                 list("BEELINE\nmedian ratio", xb, "left")))
    txt(h[[1]], h[[2]], hy, 11, ink[["secondary"]], "bold", just = h[[3]], lineheight = 0.95)
  for (xx in c(xk, xb)) grid.roundrect(xx - 0.1, hy, 1.15, 0.55, default.units = "native", just = c("left", "centre"),
                                       r = unit(0.05, "inches"), gp = DOT)
  seg(0.4, hy - 0.34, 15.7, hy - 0.34, ink[["secondary"]], 1.2)

  ys <- seq(6.35, by = -0.85, length.out = nrow(v))
  ax <- function(a) xe0 + (a - 0.13) / (0.17 - 0.13) * (xe1 - xe0)
  for (a in c(0.13, 0.14, 0.15, 0.16, 0.17)) {
    seg(ax(a), ys[1] + 0.42, ax(a), ys[length(ys)] - 0.42, ink[["grid"]], 0.9)
    txt(sprintf("%.2f", a), ax(a), ys[length(ys)] - 0.5, 9, ink[["secondary"]])
  }
  txt("axis starts at 0.13", (xe0 + xe1) / 2, ys[length(ys)] - 0.7, 9, ink[["secondary"]], "italic")

  for (k in seq_len(nrow(v))) {
    y <- ys[k]; adopted <- k == 4; lead <- k == 5
    if (adopted || lead) grid.roundrect(0.3, y, 15.5, 0.74, default.units = "native", just = c("left", "centre"),
                             r = unit(0.08, "inches"),
                             gp = gpar(fill = if (adopted) "#e3f5ee" else "#f4fbf8", col = pal[["aqua"]],
                                       lwd = if (adopted) 1.4 else 1, lty = if (adopted) 1 else 2))
    txt(v$name[k], xn, y + 0.13, 12.5, face = "bold", just = "left")
    txt(v$sub[k], xn, y - 0.17, 8.5, ink[["secondary"]], just = "left")
    chip(TRUE, xd, y, 2.0, dep_lab[[v$dep[k]]], dep_cols[[v$dep[k]]])
    chip(v$cn[k], xc, y, 1.15, "CLR z", pal[["orange"]])
    chip(v$sm[k], xs, y, 1.0, "softmax", pal[["violet"]], new = TRUE)
    chip(v$sm[k], xw, y, 1.0, if (is.na(v$walklab[k])) "walk" else v$walklab[k], pal[["magenta"]], new = TRUE)
    chip(v$lc[k], xl, y, 1.45, "log + CLR", pal[["aqua"]], new = TRUE)
    grid.circle(ax(v$ecoli[k]), y, 0.08, default.units = "native",
                gp = gpar(fill = if (adopted || lead) pal[["aqua"]] else ink[["secondary"]], col = NA))
    txt(sprintf("%.3f", v$ecoli[k]), ax(v$ecoli[k]) + 0.14, y + 0.2, 9.5, just = "left", face = if (adopted || lead) "bold" else "plain")
    if (nzchar(v$corum[k])) txt(v$corum[k], xk, y, 12, just = "left", face = if (adopted || lead) "bold" else "plain")
    if (nzchar(v$beeline[k])) txt(v$beeline[k], xb, y, 12, just = "left", face = if (adopted || lead) "bold" else "plain")
  }
  fy0 <- ys[length(ys)] - 0.95
  txt("t* = held-out depth, chosen per version; the method comparison itself fixes t = 9. The E. coli column is the component ablation (each operator at its own t*), except for the lead method (t = 9).",
      0.4, fy0, 9, ink[["secondary"]], just = "left")
  txt(paste0(LEAD_NAME, " was picked among about 20 read-outs by its K562 scores (a selection, not a fit: it has no fitted parameters). It is not yet adopted under the preregistered rule."),
      0.4, fy0 - 0.25, 9, ink[["secondary"]], just = "left")
  txt(paste0("BEELINE: CLR is plain CLR; ", ACLR_NAME, " and the lead method use the CLR-NS operator (CLR on |r| of normal scores), the lead at t = 9, not tuned (exploratory)."),
      0.4, fy0 - 0.5, 9, ink[["secondary"]], just = "left")

  fy <- 0.4
  grid.roundrect(0.3, fy, 15.5, 0.5, default.units = "native", just = c("left", "centre"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "#f1f0ec", col = NA))
  txt("Compared with (BEELINE median ratio):", 0.55, fy, 11, face = "bold", just = "left")
  txt("GENIE3 1.69    PIDC 1.56    GRNBoost2 1.47", 5.15, fy, 11, just = "left")
  txt("scPRINT zero-shot attention: CLR-NS ahead in 7 of 7 datasets (+1.03 log2 AUPR)",
      9.3, fy, 11, just = "left")

  popViewport()
}

png(file.path(out_dir, "graphical_abstract_short_p2.png"), width = 2400, height = 1350, res = 150,
    type = "cairo", bg = "white")
draw_p2(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_short_p2.pdf"), width = 16, height = 9, bg = "white")
draw_p2(); invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract_short_p2.{png,pdf}"), "\n")
