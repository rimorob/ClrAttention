# Graphical abstract, page 2 (16:9): the versions of CLR we tested, one row each, drawn in code with grid.
# Companion of analysis/fig_graphical_abstract.R (page 1); same palette and conventions.
# Numbers, all copied from existing result files:
#   E. coli SC AUPR: Table 1 of docs/whitepapers/paper.md (466 experiments, strong-evidence regulons,
#     each operator at its own held-out depth t*; "CLR scores ranked directly" = CLR without softmax)
#   BEELINE median AUPRC ratio over 14 cases: results/beeline_baselines/summary.csv
#     (|r| 1.20, CLR 1.48, CLR-NS 1.60, GENIE3 1.69, PIDC 1.56, GRNBoost2 1.47)
#   scFM: results/scfm_score (BEELINE, preregistered primary: CLR-NS minus scPRINT R0, log2 AUPRC, +1.03)
# Versions without a number in those files show a dash; nothing is interpolated.
#
# Output: docs/whitepapers/figures/graphical_abstract_p2.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_p2.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")

out_dir <- "docs/whitepapers/figures"
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
chip <- function(present, x, y, w, label, col, tcol = "white") {
  if (present) pill(label, x, y, w, fill = col, col = col, tcol = tcol, size = 10.5)
  else seg(x - 0.14, y, x + 0.14, y, ink[["grid"]], 2)
}

# ---- the versions ----------------------------------------------------------------
# dep: dependence measure; cn: context normalization; sm: softmax; wk: walk
v <- data.frame(
  name = c("CLR", "CLR-NS", "|r| (correlation)",
           "Attention on |r|", "Attention on mutual information",
           "Attention on CLR of |r|", "Attention on CLR-NS", "CLR attention (this work)"),
  sub = c("context normalization, scores ranked directly",
          "CLR on |r| of normal scores (Blom)",
          "no normalization, scores ranked directly",
          "softmax + walk, no CLR normalization",
          "softmax + walk, no CLR normalization",
          "CLR normalization, then softmax + walk",
          "CLR-NS scores, then softmax + walk",
          "MI, CLR normalization, softmax + walk"),
  dep = c("MI", "NS", "r", "r", "MI", "r", "NS", "MI"),
  cn  = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  sm  = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  tstar = c(NA, NA, NA, 8, 6, 7, 5, 9),
  ecoli = c(0.136, NA, NA, 0.143, 0.143, 0.153, 0.151, 0.158),
  beeline = c("1.48", "1.60", "1.20", "", "", "", "no gain over CLR-NS", ""),
  stringsAsFactors = FALSE)

draw_p2 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))

  txt("The eight versions of CLR we tested, from raw dependence to the full method",
      0.5, 8.55, 20, face = "bold", just = "left")
  txt("E. coli: strong-evidence regulons, 466 experiments (RegulonDB).  Single cells: BEELINE, median over 14 dataset-and-gene-set cases.",
      0.5, 8.1, 11.5, ink[["secondary"]], just = "left")

  # column x positions
  xn <- 0.5; xd <- 5.4; xc <- 7.45; xs <- 8.95; xw <- 10.2; xe0 <- 11.45; xe1 <- 13.6; xb <- 14.05
  hy <- 7.55
  for (h in list(list("Version", xn, "left"), list("Dependence", xd, "centre"),
                 list("Context\nnormalization", xc, "centre"), list("Softmax", xs, "centre"),
                 list("Walk", xw, "centre"), list("E. coli\nAUPR (SC)", (xe0 + xe1) / 2, "centre"),
                 list("BEELINE\nmedian ratio", xb, "left")))
    txt(h[[1]], h[[2]], hy, 11, ink[["secondary"]], "bold", just = h[[3]], lineheight = 0.95)
  seg(0.5, hy - 0.34, 15.5, hy - 0.34, ink[["secondary"]], 1.2)

  ys <- seq(6.78, by = -0.7, length.out = nrow(v))
  # E. coli dot axis 0.13 .. 0.16 (deliberately not from zero; ticks labelled)
  ax <- function(a) xe0 + (a - 0.13) / (0.16 - 0.13) * (xe1 - xe0)
  for (a in c(0.13, 0.14, 0.15, 0.16)) {
    seg(ax(a), ys[1] + 0.3, ax(a), ys[length(ys)] - 0.3, ink[["grid"]], 0.9)
    txt(sprintf("%.2f", a), ax(a), ys[length(ys)] - 0.52, 9.5, ink[["secondary"]])
  }
  txt("axis starts at 0.13", (xe0 + xe1) / 2, ys[length(ys)] - 0.8, 9.5, ink[["secondary"]], "italic")

  for (k in seq_len(nrow(v))) {
    y <- ys[k]; last <- k == nrow(v)
    if (last) grid.roundrect(0.35, y, 15.3, 0.64, default.units = "native", just = c("left", "centre"),
                             r = unit(0.08, "inches"),
                             gp = gpar(fill = "#e3f5ee", col = pal[["aqua"]], lwd = 1.4))
    txt(v$name[k], xn, y + 0.12, 13, face = "bold", just = "left")
    txt(v$sub[k], xn, y - 0.17, 10, ink[["secondary"]], just = "left")
    chip(TRUE, xd, y, 2.3, dep_lab[[v$dep[k]]], dep_cols[[v$dep[k]]])
    chip(v$cn[k], xc, y, 1.35, "CLR z", pal[["orange"]])
    chip(v$sm[k], xs, y, 1.1, "softmax", pal[["violet"]])
    chip(v$sm[k], xw, y, 1.0, if (is.na(v$tstar[k])) "walk" else sprintf("t* = %d", v$tstar[k]),
         pal[["magenta"]])
    if (!is.na(v$ecoli[k])) {
      grid.circle(ax(v$ecoli[k]), y, 0.085, default.units = "native",
                  gp = gpar(fill = if (last) pal[["aqua"]] else ink[["secondary"]], col = NA))
      txt(sprintf("%.3f", v$ecoli[k]), ax(v$ecoli[k]) + 0.16, y + 0.2, 9.5, just = "left",
          face = if (last) "bold" else "plain")
    } else txt("not in the ablation", (xe0 + xe1) / 2, y, 9.5, ink[["muted"]], "italic")
    if (nzchar(v$beeline[k])) txt(v$beeline[k], xb, y, 11, just = "left",
                                  face = if (k == 2) "bold" else "plain")
  }
  txt("t* = held-out depth, chosen per version; the method comparison itself fixes t = 9.",
      0.5, ys[length(ys)] - 0.52, 9.5, ink[["secondary"]], just = "left")

  # footer strip: what CLR was compared with
  fy <- 0.62
  grid.roundrect(0.35, fy, 15.3, 0.62, default.units = "native", just = c("left", "centre"),
                 r = unit(0.08, "inches"), gp = gpar(fill = "#f1f0ec", col = NA))
  txt("Compared with (BEELINE median ratio):", 0.55, fy, 11, face = "bold", just = "left")
  txt("GENIE3 1.69    PIDC 1.56    GRNBoost2 1.47", 5.15, fy, 11, just = "left")
  txt("scPRINT zero-shot attention: CLR-NS ahead in 7 of 7 datasets (+1.03 log2 AUPR)",
      9.3, fy, 11, just = "left")

  popViewport()
}

png(file.path(out_dir, "graphical_abstract_p2.png"), width = 2400, height = 1350, res = 150,
    type = "cairo", bg = "white")
draw_p2(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_p2.pdf"), width = 16, height = 9, bg = "white")
draw_p2(); invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract_p2.{png,pdf}"), "\n")
