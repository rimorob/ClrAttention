# Graphical abstract, page 3 (16:9): "Every metric, wins and losses". A descriptive scoreboard, drawn in code with grid.
# Companion of analysis/fig_graphical_abstract.R (page 1) and _p2.R; same palette and conventions.
# Numbers: docs/whitepapers/figures/scoreboard_data.csv (built by analysis/scoreboard_matrix.R from existing result files and the
# scoreboard jobs of docs/plans/scoreboard_plan.md). Nothing is a test: no intervals, no p-values, no winner declared.
# Two pages are written because one page would force type below 9 pt:
#   graphical_abstract_p3   raw scores
#   graphical_abstract_p3b  hub-controlled versions (within-stratum percentile in 10 strata of pi_i + pi_j) of the rows that have them
# Cells: the number, and a bar from random (tick = baseline, 1.00 except the AUROC row, chance 0.50) with length proportional to
# (value - baseline), scaled within the ROW. Highest value of the row: bold, dot. Ties are not marked here. Missing cell: grey dash ("not run").
#
# Output: docs/whitepapers/figures/graphical_abstract_p3{,b}.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_p3.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")
out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"
D <- utils::read.csv(file.path(out_dir, "scoreboard_data.csv"), stringsAsFactors = FALSE)

txt <- function(label, x, y, size = 10, col = ink[["primary"]], face = "plain", just = "centre", lineheight = 1.0)
  grid.text(label, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face, lineheight = lineheight))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1)
  grid.segments(x0, y0, x1, y1, default.units = "native", gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))

# ---- columns ----------------------------------------------------------------------------------------------------------
meths <- c("Attention", "CLR", "a3r", "a3r2", "a3r3", "GENIE3", "PIDC", "GRNBoost2", "abs_r")
head1 <- c("Attention", "CLR", "a3r", "a3r2", "a3r3", "GENIE3", "PIDC", "GRNBoost2", "|r|")
head2 <- c("adopted", "", "exploratory", "exploratory", "exploratory", "", "", "", "")
EXPL <- "#bfe6d8"
fillc <- c(Attention = pal[["aqua"]], CLR = pal[["orange"]], a3r = EXPL, a3r2 = EXPL, a3r3 = EXPL,
           GENIE3 = ink[["muted"]], PIDC = ink[["muted"]], GRNBoost2 = ink[["muted"]], abs_r = ink[["muted"]])
X_LAB <- 0.3                    # row labels start
X_C0 <- 4.75                    # first method column starts
CW <- (15.75 - X_C0) / 9        # column width
NUM_X <- 0.16                   # number, left edge inside a cell
ORG_X <- 0.62                   # bar origin (the 'random' tick) inside a cell
BAR_W <- CW - ORG_X - 0.06      # longest bar
cx0 <- function(j) X_C0 + (j - 1) * CW

# ---- row specs ---------------------------------------------------------------------------------------------------------
# list(label, dataset, truth, metric, scope, base (baseline of the bar), digits)
R <- function(label, ds, truth, metric, scope = "raw", base = 1, digits = 2) list(label = label, ds = ds, truth = truth, metric = metric, scope = scope, base = base, digits = digits)
PM <- "pooled_aupr_over_base"
raw_rows <- list(
  list(group = "E. coli, M3D compendium (466 experiments, RegulonDB)"),
  R("Regulons, strong evidence", "m3d", "SC", PM),
  R("Regulons, all evidence", "m3d", "all", PM),
  R("Direct TF-to-target edges (AUPR ratio)", "m3d", "SC", "tf_aupr_ratio"),
  R("Direct edges: median per-TF AUROC (chance 0.50)", "m3d", "SC", "tf_auroc_median", base = 0.5, digits = 3),
  R("Per-regulon mean AP, strong evidence", "m3d", "SC", "per_regulon_mean_ap_over_base"),
  R("Per-regulon mean AP, all evidence", "m3d", "all", "per_regulon_mean_ap_over_base"),
  R("Protein complexes (Complex Portal)", "m3d", "complexes", PM),
  list(group = "E. coli, EcoMAC compendium (2,198 arrays)"),
  R("Regulons, strong evidence", "ecomac", "SC", PM),
  R("Regulons, all evidence", "ecomac", "all", PM),
  R("Protein complexes (Complex Portal)", "ecomac", "complexes", PM),
  list(group = "Human K562 Perturb-seq (essential genes)"),
  R("DoRothEA A-B (TF targets)", "k562", "DoRothEA_AB", PM),
  R("DoRothEA A-B-C (TF targets)", "k562", "DoRothEA_ABC", PM),
  R("CollecTRI (TF targets)", "k562", "CollecTRI", PM),
  R("CORUM protein complexes", "k562", "CORUM", PM),
  R("Knockdown-response modules", "k562", "KD_response", PM),
  list(group = "Human RPE1 Perturb-seq"),
  R("CORUM protein complexes (merged)", "rpe1", "CORUM_merged", PM),
  list(group = "Single cells (BEELINE: 14 dataset-and-gene-set cases, 44 case-truth cells)"),
  R("Median AUPRC ratio, primary truths", "beeline", "all_14_cases", "beeline_median_auprc_ratio_primary"),
  R("Median AUPRC ratio, cell-type ChIP truth", "beeline", "all_14_cases", "beeline_median_auprc_ratio_celltypeChIP"))
hub_rows <- list(
  list(group = "E. coli, M3D compendium"),
  R("Regulons, strong evidence", "m3d", "SC", PM, "hub_att"),
  R("Regulons, all evidence", "m3d", "all", PM, "hub_att"),
  R("Protein complexes (Complex Portal)", "m3d", "complexes", PM, "hub_att"),
  list(group = "E. coli, EcoMAC compendium"),
  R("Regulons, strong evidence", "ecomac", "SC", PM, "hub_att"),
  R("Regulons, all evidence", "ecomac", "all", PM, "hub_att"),
  R("Protein complexes (Complex Portal)", "ecomac", "complexes", PM, "hub_att"),
  list(group = "Human K562 Perturb-seq"),
  R("DoRothEA A-B (TF targets)", "k562", "DoRothEA_AB", PM, "hub_att"),
  R("DoRothEA A-B-C (TF targets)", "k562", "DoRothEA_ABC", PM, "hub_att"),
  R("CollecTRI (TF targets)", "k562", "CollecTRI", PM, "hub_att"),
  R("CORUM protein complexes", "k562", "CORUM", PM, "hub_att"),
  R("Knockdown-response modules", "k562", "KD_response", PM, "hub_att"),
  list(group = "Human RPE1 Perturb-seq"),
  R("CORUM protein complexes (merged)", "rpe1", "CORUM_merged", PM, "hub_att"))

cell_vals <- function(r) { x <- D[D$dataset == r$ds & D$truth == r$truth & D$metric == r$metric & D$scope == r$scope, ]
  v <- stats::setNames(rep(NA_real_, length(meths)), meths)
  ok <- x$status != "not_run" & !is.na(x$enrichment_over_random) & x$method %in% meths
  v[x$method[ok]] <- x$enrichment_over_random[ok]; v }
fmt <- function(v, digits) if (digits == 2 && v >= 10) sprintf("%.1f", v) else sprintf(paste0("%.", digits, "f"), v)

# ---- one page ----------------------------------------------------------------------------------------------------------
draw_page <- function(rows, sub, y_top, footer) {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))
  txt("Every metric, wins and losses", 0.3, 8.65, 21, face = "bold", just = "left")
  txt(sub, 0.3, 8.28, 11, ink[["secondary"]], just = "left")
  nr <- sum(vapply(rows, function(r) is.null(r$group), NA)); ng <- length(rows) - nr
  y_bot <- 1.12; GH <- 0.74                                   # group-header height relative to a row
  pitch <- (y_top - y_bot) / (nr + ng * GH)
  # column bands (adopted and exploratory columns) and headers
  hy <- y_top + 0.5
  grid.rect(cx0(1), y_bot - 0.05, CW, y_top - y_bot + 0.78, default.units = "native", just = c("left", "bottom"), gp = gpar(fill = "#eaf7f1", col = NA))
  grid.rect(cx0(3), y_bot - 0.05, 3 * CW, y_top - y_bot + 0.78, default.units = "native", just = c("left", "bottom"), gp = gpar(fill = "#f6fcf9", col = pal[["aqua"]], lty = 2, lwd = 1))
  for (j in seq_along(meths)) {
    x <- cx0(j) + CW / 2
    txt(head1[j], x, hy + 0.12, 11, face = "bold")
    if (nzchar(head2[j])) txt(head2[j], x, hy - 0.12, 9.5, if (j == 1) "#0d6e4c" else ink[["secondary"]], "italic")
  }
  seg(X_LAB, y_top + 0.1, 15.75, y_top + 0.1, ink[["secondary"]], 1)
  y <- y_top
  for (r in rows) {
    if (!is.null(r$group)) {
      y <- y - GH * pitch / 2
      txt(r$group, X_LAB, y, 10.5, ink[["primary"]], "bold", just = "left")
      seg(X_LAB, y - GH * pitch / 2 + 0.005, 15.75, y - GH * pitch / 2 + 0.005, ink[["grid"]], 0.9)
      y <- y - GH * pitch / 2; next }
    y <- y - pitch / 2
    v <- cell_vals(r); best <- if (all(is.na(v))) NA else which.max(v)
    span <- max(v - r$base, na.rm = TRUE); if (!is.finite(span) || span <= 0) span <- 1
    txt(r$label, X_LAB + 0.12, y, 10, just = "left")
    for (j in seq_along(meths)) {
      x0 <- cx0(j)
      seg(x0 + ORG_X, y - pitch * 0.42, x0 + ORG_X, y + pitch * 0.42, ink[["secondary"]], 0.8)          # random
      if (is.na(v[j])) { txt("–", x0 + NUM_X + 0.1, y, 11, ink[["muted"]]); next }
      frac <- (v[j] - r$base) / span; len <- if (frac >= 0) frac * BAR_W else max(frac, -0.08 / BAR_W) * BAR_W
      expl <- meths[j] %in% c("a3r", "a3r2", "a3r3")
      bh <- pitch * 0.5
      grid.rect(x0 + ORG_X, y, abs(len), bh, default.units = "native", just = c(if (len >= 0) "left" else "right", "centre"),
                gp = gpar(fill = fillc[[meths[j]]], col = if (expl) pal[["aqua"]] else NA, lwd = 1.1, lty = if (expl) 2 else 1))
      isb <- !is.na(best) && j == best
      if (isb) grid.circle(x0 + 0.07, y, 0.035, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
      txt(fmt(v[j], r$digits), x0 + NUM_X + 0.02, y, 10, face = if (isb) "bold" else "plain", just = "left")
    }
    y <- y - pitch / 2
  }
  # footer
  for (k in seq_along(footer)) txt(footer[[k]]$text, 0.3, y_bot - 0.22 - (k - 1) * 0.185, footer[[k]]$size, footer[[k]]$col, footer[[k]]$face, just = "left")
  popViewport()
}

foot_common <- list(
  list(text = "Bars start at random (tick; 1.00, or 0.50 for the AUROC row) and are scaled within each row; the value, not the bar, carries the comparison across rows. Bold with a dot: highest value in the row.",
       size = 9, col = ink[["secondary"]], face = "plain"),
  list(text = "Ties are not marked here (see the method table, docs/whitepapers/sota_table.md, for Bonferroni-tied cells). No intervals, no tests. A grey dash is a cell that was not run.",
       size = 9, col = ink[["secondary"]], face = "plain"),
  list(text = "Pooled AUPR / base rate (1.00 = random), raw scores unless marked hub-controlled; BEELINE median over the 44 case-truth cells of the primary truths (14 cases); a3r was picked among about 20 read-outs by its K562 scores (selected, not fitted);",
       size = 9, col = ink[["secondary"]], face = "plain"),
  list(text = "K562 and M3D are tuning sets; attention and a3r settings were developed on M3D; DREAM5 yeast is separate (page 2 / section 8) and was not run for a3r.",
       size = 9, col = ink[["secondary"]], face = "plain"))
foot_raw <- c(foot_common, list(
  list(text = "BEELINE: scores of CLR-NS for the attention and a3r columns (attention at its held-out depth, a3r family at t = 9, not tuned); CLR is plain CLR; median over 14 cases x BEELINE primary truths (44 case-truth cells), second row cell-type ChIP only.",
       size = 9, col = ink[["secondary"]], face = "plain")))
foot_hub <- c(foot_common, list(
  list(text = "Hub-controlled: each method's scores are converted to percentiles within 10 strata of pi_i + pi_j (pi of the attention walk), the stratifier of the method table. Not run for direct edges, per-regulon AP and BEELINE.",
       size = 9, col = ink[["secondary"]], face = "plain")))

png(file.path(out_dir, "graphical_abstract_p3.png"), width = 2400, height = 1350, res = 150, type = "cairo", bg = "white")
draw_page(raw_rows, "Raw scores. Adopted method (attention, t = 9) and CLR, then three exploratory a3r variants (dashed), then competitors.", 7.3, foot_raw); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_p3.pdf"), width = 16, height = 9, bg = "white")
draw_page(raw_rows, "Raw scores. Adopted method (attention, t = 9) and CLR, then three exploratory a3r variants (dashed), then competitors.", 7.3, foot_raw); invisible(dev.off())
png(file.path(out_dir, "graphical_abstract_p3b.png"), width = 2400, height = 1350, res = 150, type = "cairo", bg = "white")
draw_page(hub_rows, "Hub-controlled versions of the same rows (page 3 shows the raw scores): controlling for highly connected (hub) genes.", 7.3, foot_hub); invisible(dev.off())
cairo_pdf(file.path(out_dir, "graphical_abstract_p3b.pdf"), width = 16, height = 9, bg = "white")
draw_page(hub_rows, "Hub-controlled versions of the same rows (page 3 shows the raw scores): controlling for highly connected (hub) genes.", 7.3, foot_hub); invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract_p3{,b}.{png,pdf}"), "\n")
