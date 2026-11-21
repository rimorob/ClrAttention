# Graphical abstract, page 2 (portrait, 8.3 in wide canvas, printed 7.0 in wide): the versions of CLR we tested, one row each, drawn in code with grid.
# Companion of analysis/fig_graphical_abstract.R (page 1); same palette and conventions.
# Numbers, all copied from existing result files:
#   E. coli SC AUPR: Table 1 of docs/whitepapers/paper.md (466 experiments, strong-evidence regulons,
#     each operator at its own held-out depth t*; "CLR scores ranked directly" = CLR without softmax)
#   BEELINE median AUPRC ratio over 14 cases: results/beeline_baselines/summary.csv
#     (|r| 1.20, CLR 1.48, CLR-NS 1.60, GENIE3 1.69, PIDC 1.56, GRNBoost2 1.47)
#   scFM: results/scfm_score (BEELINE, preregistered primary: CLR-NS minus scPRINT R0, log2 AUPRC, +1.03)
#   a3r row (added 13 Nov; exploratory, not adopted): E. coli SC AUPR at t = 9 from
#     results/explore_stack/m3d/pooled_levels.csv (a3r 0.1622); K562 CORUM enrichment
#     from results/explore_stack/k562/pooled_levels.csv (7.95); K562 CORUM for CLR 6.39, |r| 7.00,
#     CLR attention 3.51 from docs/whitepapers/sota_table.md. a3r was picked among ~20 readouts on K562 truths.
# Names (19 Nov): aCLR = the adopted CLR attention (earlier text "attention"), tcaCLR = the log + second CLR read-out
# (earlier text "a3r"). The stacked variants of tcaCLR are not shown. Defined below.
# BEELINE for the adopted method and the tcaCLR family: median AUPRC ratio over the 44 case-truth cells of the primary truths,
#   from docs/whitepapers/figures/scoreboard_data.csv (CLR-NS operator; the tcaCLR family at t = 9, not tuned): 1.54, 1.73.
# Versions without a number in those files show a dash; nothing is interpolated.
#
# Output: docs/whitepapers/figures/graphical_abstract_p2.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_p2.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")

out_dir <- Sys.getenv("SB_OUT", "docs/whitepapers/figures"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
FONT <- "Liberation Sans"
ACLR_NAME <- "aCLR"       # DEFINED HERE ONLY (temporary name; earlier text "attention")
LEAD_NAME <- "tcaCLR"     # DEFINED HERE ONLY (temporary name; earlier text "a3r")

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
chip <- function(present, x, y, w, label, col, tcol = "white") {
  if (present) pill(label, x, y, w, fill = col, col = col, tcol = tcol) else seg(x - 0.11, y, x + 0.11, y, ink[["grid"]], 2)
}

# ---- the versions ----------------------------------------------------------------
# dep: dependence measure; cn: context normalization; sm: softmax; wk: walk
v <- data.frame(
  name = c("CLR", "CLR-NS", "|r| (correlation)",
           "Attention on |r|", "Attention on mutual information",
           "Attention on CLR of |r|", "Attention on CLR-NS", paste0(ACLR_NAME, " (adopted)"),
           paste0(LEAD_NAME, " (lead candidate)"), paste0(T2_NAME, " (exploratory)")),
  sub = c("context normalization, scores ranked directly",
          "CLR on |r| of normal scores (Blom)",
          "no normalization, scores ranked directly",
          "softmax + walk, no CLR normalization",
          "softmax + walk, no CLR normalization",
          "CLR normalization, then softmax + walk",
          "CLR-NS scores, then softmax + walk",
          "MI, CLR normalization, softmax + walk",
          paste0(ACLR_NAME, ", then log of the walk + a second CLR"),
          paste0(LEAD_NAME, " with the walk depth set by the data, not t = 9")),
  dep = c("MI", "NS", "r", "r", "MI", "r", "NS", "MI", "MI", "MI"),
  cn  = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  sm  = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  lc  = c(rep(FALSE, 8), TRUE, TRUE),
  walklab = c(NA, NA, NA, "t* = 8", "t* = 6", "t* = 7", "t* = 5", "t* = 9", "t = 9", "t = t2"),
  ecoli = c(0.136, NA, NA, 0.143, 0.143, 0.153, 0.151, 0.158, 0.1622, t2_ecoli),
  corum = c("6.39", "", "7.00", "", "", "", "", "3.51", "7.95", t2_corum),
  beeline = c("1.48", "1.60", "1.20", "", "", "", "no gain over CLR-NS", "1.54", "1.73", t2_beel),
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
  ttl <- wrap_to("The versions of CLR we tested, from raw dependence to the full method and its exploratory extensions", 14.5, XR - XL, "bold")
  for (k in seq_along(ttl)) txt(ttl[k], XL, T - 0.3 - (k - 1) * 0.25, 14.5, face = "bold", just = "left")
  y <- T - 0.3 - (length(ttl) - 1) * 0.25 - 0.3
  sl <- wrap_to("E. coli: strong-evidence regulons, 466 experiments (RegulonDB).  K562: CORUM complexes.  Single cells: BEELINE, median over 44 case-truth cells (14 dataset-and-gene-set cases).", 9.4, XR - XL)
  for (k in seq_along(sl)) { txt(sl[k], XL, y, 9.4, ink[["secondary"]], just = "left"); y <- y - 0.17 }
  y <- y - 0.38
  # headers
  hy <- y
  for (h in list(list("Version", XL, "left", 0.17), list("Dependence", xd, "centre", -0.1),
                 list("Context\nnormalization", xc + 0.0, "centre"), list("Softmax", xs, "centre"),
                 list("Walk", xw, "centre"), list("Log +\nsecond CLR", xl, "centre"),
                 list("E. coli\nAUPR (SC)", (xe0 + xe1) / 2, "centre"),
                 list("K562\nCORUM\n(\u00d7 random)", xk, "left"),
                 list("BEELINE\nmedian\nratio", xb, "left")))
    txt(h[[1]], h[[2]], hy + (if (length(h) > 3) h[[4]] else 0), FS_N, ink[["secondary"]], "bold", just = h[[3]], lineheight = 0.95)
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
    yc <- ys[k]; adopted <- k == 8; lead <- k > 8
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
    chip(v$sm[k], xs, yb, ws, "softmax", pal[["violet"]])
    chip(v$sm[k], xw, yb, ww, if (is.na(v$walklab[k])) "walk" else v$walklab[k], pal[["magenta"]])
    chip(v$lc[k], xl, yb, wl, "log + CLR", pal[["aqua"]])
    if (!is.na(v$ecoli[k])) {
      grid.circle(ax(v$ecoli[k]), yb, 0.07, default.units = "native",
                  gp = gpar(fill = if (adopted || lead) pal[["aqua"]] else ink[["secondary"]], col = NA))
      right <- v$ecoli[k] > 0.16
      txt(sprintf("%.3f", v$ecoli[k]), ax(v$ecoli[k]) + (if (right) -0.1 else 0.1), yb + 0.14, FS_N, just = if (right) "right" else "left", face = if (adopted) "bold" else "plain")
    } else txt("not in the ablation", (xe0 + xe1) / 2, yb, FS_N, ink[["muted"]], "italic")
    if (nzchar(v$corum[k])) txt(v$corum[k], xk, yb, 10, just = "left", face = if (adopted || lead) "bold" else "plain")
    if (nzchar(v$beeline[k])) txt(if (k == 7) "no gain over\nCLR-NS" else v$beeline[k], if (k == 7) xk else xb, yb, if (k == 7) FS_N else 10, just = "left", face = if (k == 2) "bold" else "plain", lineheight = 0.95)
  }
  yf <- y_end - 0.5
  notes <- c("t* = held-out depth, chosen per version; the method comparison itself fixes t = 9.",
    paste0(LEAD_NAME, ": lead candidate; its 26 Oct preregistered E. coli test did not reach the threshold. It was picked among about 20 read-outs by its K562 scores (a selection, not a fit: it has no fitted parameters)."),
    paste0("BEELINE: CLR is plain CLR; the ", ACLR_NAME, " and ", LEAD_NAME, " rows use the CLR-NS operator (", LEAD_NAME, " family at t = 9, not tuned)."),
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
png(file.path(out_dir, "graphical_abstract_p2.png"), width = round(W_CAN * 300), height = round(H_CAN * 300), res = 300,
    type = "cairo", bg = "white")
invisible(draw_p2(H_CAN)); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_p2.pdf"), width = W_CAN, height = H_CAN, bg = "white")
invisible(draw_p2(H_CAN)); invisible(dev.off())
cat(sprintf("p2: canvas %.2f x %.2f in; printed %.2f x %.2f in (scale %.3f); min font 8.6 pt -> %.1f pt printed\n", W_CAN, H_CAN, 6.9, H_CAN * 6.9 / W_CAN, 6.9 / W_CAN, 8.6 * 6.9 / W_CAN))
cat("wrote", file.path(out_dir, "graphical_abstract_p2.{png,pdf}"), "\n")
