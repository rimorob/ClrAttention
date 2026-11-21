# Graphical abstract for the SHORT pre-read, page 2 (portrait, 8.3 in wide canvas, printed 7.0 in wide): the versions compared, one row each, drawn in code with grid.
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

out_dir <- Sys.getenv("SB_OUT", "docs/whitepapers/figures"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
LEAD_NAME <- "tcaCLR"                     # DEFINED HERE ONLY (temporary name; internal name a3r)
ACLR_NAME <- "aCLR"                       # DEFINED HERE ONLY (temporary name)
NAMES_NOTE <- "Method names are temporary."
LEAD_STATUS <- "lead candidate"
GREEN <- pal[["green"]]
DOT <- gpar(col = GREEN, lty = "11", lwd = 2.4, fill = NA)
FONT <- "Liberation Sans"

# ---- tcaCLR2 (exploratory): values from results/tcaclr2_scoreboard/cells_tcaclr2_all.csv (analysis/tcaclr2/assemble_tcaclr2_scoreboard.R) --------------------
T2_NAME <- "tcaCLR2"
t2c <- utils::read.csv(Sys.getenv("SB_T2", "results/tcaclr2_scoreboard/cells_tcaclr2_all.csv"), stringsAsFactors = FALSE)
t2v <- function(ds, truth, metric, scope = "raw") t2c[t2c$dataset == ds & t2c$truth == truth & t2c$metric == metric & t2c$scope == scope & t2c$method == "tcaCLR2", ]
t2_ecoli <- with(t2v("m3d", "SC", "pooled_aupr_over_base"), enrichment_over_random * base_rate)           # pooled AUPR, strong-evidence regulons (the E. coli column)
t2_corum <- sprintf("%.2f", t2v("k562", "CORUM", "pooled_aupr_over_base")$enrichment_over_random)          # K562 CORUM, x random
t2_beel  <- sprintf("%.2f", t2v("beeline", "all_14_cases", "beeline_median_auprc_ratio_primary")$enrichment_over_random)
T2D <- utils::read.csv(Sys.getenv("SB_T2D", "results/tcaclr2_scoreboard/t2_by_dataset.csv"), stringsAsFactors = FALSE)
T2_NOTE <- paste0(T2_NAME, ": exploratory, outside the preregistration; ", LEAD_NAME, " with the walk depth set by the data (t2 = M3D ", T2D$t2[T2D$label == "M3D"], ", K562 ", T2D$t2[T2D$label == "K562"], "; BEELINE per case ", T2D$t2[grepl("BEELINE", T2D$label)], "), not t = 9. E. coli: the M3D strong-evidence AUPR at t2; BEELINE: CLR-NS operator, t2 per case.")

# ---- PORTRAIT layout: canvas W_CAN x H inches, printed 6.9 in wide (scale 6.9 / W_CAN = 0.831); every font here is >= 8.6 pt -> >= 7.1 pt printed
W_CAN <- 8.3
txt <- function(label, x, y, size = 12, col = ink[["primary"]], face = "plain",
                just = "centre", lineheight = 1.0)
  grid.text(label, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face,
                      lineheight = lineheight))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1)
  grid.segments(x0, y0, x1, y1, default.units = "native",
                gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "round"))
FS <- 8.8; FS_N <- 8.6        # chip text, small text (canvas pt)
pill <- function(label, x, y, w, h = 0.3, fill = "white", col = ink[["primary"]], size = FS,
                 face = "bold", tcol = ink[["primary"]], lwd = 1.1, lty = 1)
{
  grid.roundrect(x, y, w, h, default.units = "native", r = unit(h / 2, "inches"),
                 gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
  txt(label, x, y, size, tcol, face)
}
text_w <- function(s, size, face = "plain") convertWidth(grobWidth(textGrob(s, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face))), "inches", TRUE)
wrap_to <- function(s, size, width, face = "plain") {     # word wrap (never inside a word) to a width in inches
  words <- strsplit(s, " ", fixed = TRUE)[[1]]; lines <- character(); cur <- ""
  for (w in words) { cand <- if (nzchar(cur)) paste(cur, w) else w
    if (!nzchar(cur) || text_w(cand, size, face) <= width) cur <- cand else { lines <- c(lines, cur); cur <- w } }
  c(lines, cur) }

# component chips: filled = present in this version, dashed grey dash = absent
dep_cols <- c(MI = pal[["blue"]], r = ink[["muted"]], NS = ink[["secondary"]])
dep_lab  <- c(MI = "mutual information", r = "|r|", NS = "|r| of normal scores")
chip <- function(present, x, y, w, label, col, tcol = "white", new = FALSE) {
  if (present) {
    pill(label, x, y, w, fill = col, col = col, tcol = tcol)
    if (new) grid.roundrect(x, y, w + 0.1, 0.3 + 0.08, default.units = "native", r = unit(0.19, "inches"), gp = DOT)
  } else seg(x - 0.11, y, x + 0.11, y, ink[["grid"]], 2)
}

# ---- the versions ----------------------------------------------------------------
v <- data.frame(
  name = c("CLR (historical)", "Attention on |r|", "Attention on mutual information",
           ACLR_NAME, LEAD_NAME, T2_NAME),
  sub = c("Faith et al. 2007: context normalization, scores ranked directly",
          "softmax + walk, no CLR normalization",
          "softmax + walk, no CLR normalization",
          "MI, CLR normalization, softmax + walk (adopted)",
          paste0(ACLR_NAME, ", then log + second CLR (lead candidate)"),
          paste0(LEAD_NAME, " with the walk depth set by the data (exploratory)")),
  dep = c("MI", "r", "MI", "MI", "MI", "MI"),
  cn  = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
  sm  = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  lc  = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  walklab = c(NA, "t* = 8", "t* = 6", "t* = 9", "t = 9", "t = t2"),
  ecoli = c(0.136, 0.143, 0.143, 0.158, 0.1622, t2_ecoli),
  corum = c("6.39", "", "", "3.51", "7.95", t2_corum),
  beeline = c("1.48", "", "", "1.54", "1.73", t2_beel),
  stringsAsFactors = FALSE)

# column geometry: each row is two lines, name + description on top, the chips and numbers below
XL <- 0.15; XR <- W_CAN - 0.15
xd <- 0.97; wd <- 1.5; xc <- 2.17; wc <- 0.66; xs <- 3.0; ws <- 0.8; xw <- 3.9; ww <- 0.7; xl <- 4.85; wl <- 1.0
xe0 <- 5.6; xe1 <- 6.65; xk <- 6.8; xb <- 7.64
RH2 <- 0.66

# Two passes (as in the scoreboards): pass 1 on a throw-away device measures the height, pass 2 draws on a canvas of exactly that height.
draw_p2 <- function(H) {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, W_CAN), yscale = c(0, H)))
  T <- H
  ttl <- wrap_to(paste0("The methods compared: CLR, ", ACLR_NAME, ", ", LEAD_NAME, " and ", T2_NAME), 15, XR - XL, "bold")
  for (k in seq_along(ttl)) txt(ttl[k], XL, T - 0.3 - (k - 1) * 0.26, 15, face = "bold", just = "left")
  y <- T - 0.3 - (length(ttl) - 1) * 0.26 - 0.3
  sl <- wrap_to("E. coli: strong-evidence regulons, 466 experiments (RegulonDB).  K562: CORUM complexes.  Single cells: BEELINE, median over 44 case-truth cells (14 dataset-and-gene-set cases).", 9.4, XR - XL)
  for (k in seq_along(sl)) { txt(sl[k], XL, y, 9.4, ink[["secondary"]], just = "left"); y <- y - 0.17 }
  y <- y - 0.1
  grid.roundrect(XL, y, 0.34, 0.16, default.units = "native", just = c("left", "centre"), r = unit(0.03, "inches"), gp = DOT)
  lg <- wrap_to("dotted green = new in this work (the softmax, the walk, the log + second CLR, and all human data)", 9.2, 5.6)
  for (k in seq_along(lg)) txt(lg[k], XL + 0.43, y - (k - 1) * 0.16, 9.2, ink[["secondary"]], just = "left")
  txt(NAMES_NOTE, XR, y, 9.2, ink[["secondary"]], "italic", just = "right")
  y <- y - (length(lg) - 1) * 0.16 - 0.4
  # headers
  hy <- y
  for (h in list(list("Method", XL, "left", 0.17), list("Dependence", xd, "centre", -0.1),
                 list("Context\nnormalization", xc + 0.0, "centre"), list("Softmax", xs, "centre"),
                 list("Walk", xw, "centre"), list("Log +\nsecond CLR", xl, "centre"),
                 list("E. coli\nAUPR (SC)", (xe0 + xe1) / 2, "centre"),
                 list("K562\nCORUM\n(x random)", xk, "left"),
                 list("BEELINE\nmedian\nratio", xb, "left")))
    txt(h[[1]], h[[2]], hy + (if (length(h) > 3) h[[4]] else 0), FS_N, ink[["secondary"]], "bold", just = h[[3]], lineheight = 0.95)
  for (xx in c(xk, xb)) grid.roundrect(xx - 0.06, hy, if (xx == xk) 0.76 else 0.68, 0.56, default.units = "native", just = c("left", "centre"),
                                       r = unit(0.04, "inches"), gp = DOT)
  y_rule <- hy - 0.36
  seg(XL, y_rule, XR, y_rule, ink[["secondary"]], 1.1)
  n <- nrow(v)
  ys <- y_rule - 0.06 - RH2 * (seq_len(n) - 0.5)         # row centres (block of two lines)
  ax <- function(a) xe0 + (a - 0.13) / (0.17 - 0.13) * (xe1 - xe0)
  y_end <- ys[n] - RH2 / 2
  for (a in c(0.13, 0.14, 0.15, 0.16, 0.17)) {
    seg(ax(a), ys[1] + RH2 / 2 - 0.2, ax(a), y_end + 0.03, ink[["grid"]], 0.8)
    txt(sprintf("%.2f", a), ax(a), y_end - 0.12, FS_N, ink[["secondary"]])
  }
  txt("axis starts at 0.13", (xe0 + xe1) / 2, y_end - 0.29, FS_N, ink[["secondary"]], "italic")
  for (k in seq_len(n)) {
    yc <- ys[k]; adopted <- k == 4; lead <- k >= 5
    if (adopted || lead) grid.roundrect(XL - 0.08, yc, XR - XL + 0.16, RH2 - 0.04, default.units = "native", just = c("left", "centre"),
                             r = unit(0.07, "inches"),
                             gp = gpar(fill = if (adopted) "#e3f5ee" else "#f4fbf8", col = pal[["aqua"]],
                                       lwd = if (adopted) 1.3 else 0.9, lty = if (adopted) 1 else 2))
    yt <- yc + 0.17; yb <- yc - 0.12
    txt(v$name[k], XL, yt, 10.3, face = "bold", just = "left")
    nw <- text_w(v$name[k], 10.3, "bold")
    txt(v$sub[k], XL + nw + 0.12, yt - 0.005, FS_N, ink[["secondary"]], just = "left")
    chip(TRUE, xd, yb, wd, dep_lab[[v$dep[k]]], dep_cols[[v$dep[k]]])
    chip(v$cn[k], xc, yb, wc, "CLR z", pal[["orange"]])
    chip(v$sm[k], xs, yb, ws, "softmax", pal[["violet"]], new = TRUE)
    chip(v$sm[k], xw, yb, ww, if (is.na(v$walklab[k])) "walk" else v$walklab[k], pal[["magenta"]], new = TRUE)
    chip(v$lc[k], xl, yb, wl, "log + CLR", pal[["aqua"]], new = TRUE)
    grid.circle(ax(v$ecoli[k]), yb, 0.07, default.units = "native",
                gp = gpar(fill = if (adopted || lead) pal[["aqua"]] else ink[["secondary"]], col = NA))
    right <- v$ecoli[k] > 0.16
    txt(sprintf("%.3f", v$ecoli[k]), ax(v$ecoli[k]) + (if (right) -0.1 else 0.1), yb + 0.14, FS_N, just = if (right) "right" else "left", face = if (adopted || lead) "bold" else "plain")
    if (nzchar(v$corum[k])) txt(v$corum[k], xk, yb, 10, just = "left", face = if (adopted || lead) "bold" else "plain")
    if (nzchar(v$beeline[k])) txt(v$beeline[k], xb, yb, 10, just = "left", face = if (adopted || lead) "bold" else "plain")
  }
  yf <- y_end - 0.5
  notes <- c("t* = held-out depth, chosen per method; the method comparison itself fixes t = 9. The E. coli column is the component ablation (each operator at its own t*), except for the lead method (t = 9).",
    paste0(LEAD_NAME, " was picked among about 20 read-outs by its K562 scores (a selection, not a fit: it has no fitted parameters). Its preregistered E. coli test did not reach the preregistered threshold."),
    paste0("BEELINE: CLR is plain CLR; ", ACLR_NAME, " and the lead method use the CLR-NS operator (CLR on |r| of normal scores), the lead at t = 9, not tuned (exploratory)."),
    T2_NOTE)
  for (nt in notes) { ln <- wrap_to(nt, FS_N, XR - XL)
    for (m in seq_along(ln)) { txt(ln[m], XL, yf, FS_N, ink[["secondary"]], just = "left"); yf <- yf - 0.135 }
    yf <- yf - 0.03 }
  yf <- yf - 0.1
  grid.roundrect(XL - 0.08, yf - 0.38, XR - XL + 0.16, 0.5, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.07, "inches"), gp = gpar(fill = "#f1f0ec", col = NA))
  txt("Compared with (BEELINE median ratio):", XL + 0.05, yf - 0.03, 9.2, face = "bold", just = "left")
  txt("GENIE3 1.69    PIDC 1.56    GRNBoost2 1.47", XL + 3.05, yf - 0.03, 9.2, just = "left")
  txt("scPRINT zero-shot attention: CLR-NS ahead in 7 of 7 datasets (+1.03 log2 AUPR)", XL + 0.05, yf - 0.24, 9.2, just = "left")
  yf <- yf - 0.38
  popViewport()
  H - yf + 0.12
}

png(tempfile(fileext = ".png"), width = 800, height = 800, type = "cairo")
H_CAN <- draw_p2(40); invisible(dev.off())
png(file.path(out_dir, "graphical_abstract_short_p2.png"), width = round(W_CAN * 300), height = round(H_CAN * 300), res = 300,
    type = "cairo", bg = "white")
invisible(draw_p2(H_CAN)); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_short_p2.pdf"), width = W_CAN, height = H_CAN, bg = "white")
invisible(draw_p2(H_CAN)); invisible(dev.off())
cat(sprintf("p2: canvas %.2f x %.2f in; printed %.2f x %.2f in (scale %.3f); min font 8.6 pt -> %.1f pt printed\n", W_CAN, H_CAN, 6.9, H_CAN * 6.9 / W_CAN, 6.9 / W_CAN, 8.6 * 6.9 / W_CAN))
cat("wrote", file.path(out_dir, "graphical_abstract_short_p2.{png,pdf}"), "\n")
