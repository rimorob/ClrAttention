# Graphical abstract, page 2 (16:9): the versions of CLR we tested, one row each, drawn in code with grid.
# Companion of analysis/fig_graphical_abstract.R (page 1); same palette and conventions.
# Numbers, all copied from existing result files:
#   E. coli SC AUPR: Table 1 of docs/whitepapers/paper.md (466 experiments, strong-evidence regulons,
#     each operator at its own held-out depth t*; "CLR scores ranked directly" = CLR without softmax)
#   BEELINE median AUPRC ratio over 14 cases: results/beeline_baselines/summary.csv
#     (|r| 1.20, CLR 1.48, CLR-NS 1.60, GENIE3 1.69, PIDC 1.56, GRNBoost2 1.47)
#   scFM: results/scfm_score (BEELINE, preregistered primary: CLR-NS minus scPRINT R0, log2 AUPRC, +1.03)
#   a3r rows (added 13 Nov; exploratory, not adopted): E. coli SC AUPR at t = 9 in every layer from
#     results/explore_stack/m3d/pooled_levels.csv (a3r 0.1622, a3r2 0.1629, a3r3 0.1566); K562 CORUM enrichment
#     from results/explore_stack/k562/pooled_levels.csv (7.95, 8.19, 8.15); K562 CORUM for CLR 6.39, |r| 7.00,
#     CLR attention 3.51 from docs/whitepapers/sota_table.md. a3r was picked among ~20 readouts on K562 truths.
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
           "Attention on CLR of |r|", "Attention on CLR-NS", "CLR attention (adopted)",
           "a3r (exploratory)", "a3r2: two layers (exploratory)", "a3r3: three layers (exploratory)"),
  sub = c("context normalization, scores ranked directly",
          "CLR on |r| of normal scores (Blom)",
          "no normalization, scores ranked directly",
          "softmax + walk, no CLR normalization",
          "softmax + walk, no CLR normalization",
          "CLR normalization, then softmax + walk",
          "CLR-NS scores, then softmax + walk",
          "MI, CLR normalization, softmax + walk",
          "attention, then log of the walk + a second CLR",
          "a3r fed back as the score of a second layer",
          "a3r fed back twice"),
  dep = c("MI", "NS", "r", "r", "MI", "r", "NS", "MI", "MI", "MI", "MI"),
  cn  = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  sm  = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  lc  = c(rep(FALSE, 8), TRUE, TRUE, TRUE),
  walklab = c(NA, NA, NA, "t* = 8", "t* = 6", "t* = 7", "t* = 5", "t* = 9", "t = 9", "t = 9 ×2", "t = 9 ×3"),
  ecoli = c(0.136, NA, NA, 0.143, 0.143, 0.153, 0.151, 0.158, 0.1622, 0.1629, 0.1566),
  corum = c("6.39", "", "7.00", "", "", "", "", "3.51", "7.95", "8.19", "8.15"),
  beeline = c("1.48", "1.60", "1.20", "", "", "", "no gain over CLR-NS", "", "", "", ""),
  stringsAsFactors = FALSE)

draw_p2 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))

  txt("The versions of CLR we tested, from raw dependence to the full method and its exploratory extensions",
      0.5, 8.55, 20, face = "bold", just = "left")
  txt("E. coli: strong-evidence regulons, 466 experiments (RegulonDB).  K562: CORUM complexes.  Single cells: BEELINE, median over 44 case-truth cells (14 dataset-and-gene-set cases).",
      0.5, 8.1, 11.5, ink[["secondary"]], just = "left")

  # column x positions
  xn <- 0.4; xd <- 4.95; xc <- 6.75; xs <- 8.0; xw <- 9.1; xl <- 10.6; xe0 <- 11.3; xe1 <- 12.9; xk <- 13.25; xb <- 14.5
  hy <- 7.6
  for (h in list(list("Version", xn, "left"), list("Dependence", xd, "centre"),
                 list("Context\nnormalization", xc, "centre"), list("Softmax", xs, "centre"),
                 list("Walk", xw, "centre"), list("Log +\nsecond CLR", xl, "centre"),
                 list("E. coli\nAUPR (SC)", (xe0 + xe1) / 2, "centre"),
                 list("K562 CORUM\n(× random)", xk, "left"),
                 list("BEELINE\nmedian ratio", xb, "left")))
    txt(h[[1]], h[[2]], hy, 11, ink[["secondary"]], "bold", just = h[[3]], lineheight = 0.95)
  seg(0.4, hy - 0.34, 15.7, hy - 0.34, ink[["secondary"]], 1.2)

  ys <- seq(6.85, by = -0.52, length.out = nrow(v))
  # E. coli dot axis 0.13 .. 0.17 (deliberately not from zero; ticks labelled)
  ax <- function(a) xe0 + (a - 0.13) / (0.17 - 0.13) * (xe1 - xe0)
  for (a in c(0.13, 0.14, 0.15, 0.16, 0.17)) {
    seg(ax(a), ys[1] + 0.3, ax(a), ys[length(ys)] - 0.3, ink[["grid"]], 0.9)
    txt(sprintf("%.2f", a), ax(a), ys[length(ys)] - 0.34, 9, ink[["secondary"]])
  }
  txt("axis starts at 0.13", (xe0 + xe1) / 2, ys[length(ys)] - 0.5, 9, ink[["secondary"]], "italic")

  for (k in seq_len(nrow(v))) {
    y <- ys[k]; last <- k == 8; expl <- k > 8
    if (last || expl) grid.roundrect(0.3, y, 15.5, 0.5, default.units = "native", just = c("left", "centre"),
                             r = unit(0.08, "inches"),
                             gp = gpar(fill = if (last) "#e3f5ee" else "#f4fbf8", col = pal[["aqua"]],
                                       lwd = if (last) 1.4 else 1, lty = if (last) 1 else 2))
    txt(v$name[k], xn, y + 0.1, 12, face = "bold", just = "left")
    txt(v$sub[k], xn, y - 0.14, 9, ink[["secondary"]], just = "left")
    chip(TRUE, xd, y, 2.0, dep_lab[[v$dep[k]]], dep_cols[[v$dep[k]]])
    chip(v$cn[k], xc, y, 1.15, "CLR z", pal[["orange"]])
    chip(v$sm[k], xs, y, 1.0, "softmax", pal[["violet"]])
    chip(v$sm[k], xw, y, 1.0, if (is.na(v$walklab[k])) "walk" else v$walklab[k], pal[["magenta"]])
    chip(v$lc[k], xl, y, 1.45, "log + CLR", pal[["aqua"]])
    if (!is.na(v$ecoli[k])) {
      grid.circle(ax(v$ecoli[k]), y, 0.075, default.units = "native",
                  gp = gpar(fill = if (last || expl) pal[["aqua"]] else ink[["secondary"]], col = NA))
      txt(sprintf("%.3f", v$ecoli[k]), ax(v$ecoli[k]) + 0.13, y + 0.17, 9, just = "left",
          face = if (last) "bold" else "plain")
    } else txt("not in the ablation", (xe0 + xe1) / 2, y, 9, ink[["muted"]], "italic")
    if (nzchar(v$corum[k])) txt(v$corum[k], xk, y, 11, just = "left",
                                face = if (last || expl) "bold" else "plain")
    if (nzchar(v$beeline[k])) txt(v$beeline[k], xb, y, if (k == 7) 9 else 11, just = "left",
                                  face = if (k == 2) "bold" else "plain")
  }
  txt("t* = held-out depth, chosen per version; the method comparison itself fixes t = 9. Exploratory rows: t = 9 in every layer, not re-selected; none was adopted (26 Oct preregistered test, 12 Nov stacked-layer run).",
      0.4, ys[length(ys)] - 0.5, 9, ink[["secondary"]], just = "left")
  txt("a3r was picked among about 20 read-outs by its K562 scores (a selection, not a fit: it has no fitted parameters). Gaps between a3r, a3r2 and a3r3 are inside the noise. BEELINE: not run for a3r.",
      0.4, ys[length(ys)] - 0.75, 9, ink[["secondary"]], just = "left")

  # footer strip: what CLR was compared with
  fy <- 0.4
  grid.roundrect(0.3, fy, 15.5, 0.5, default.units = "native", just = c("left", "centre"),
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
