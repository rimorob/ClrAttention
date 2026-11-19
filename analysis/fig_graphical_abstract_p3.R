# Graphical abstract, page 3 (16:9) of the LONG checkpoint: "Every metric, wins and losses". A descriptive scoreboard, drawn in code with grid.
# Companion of analysis/fig_graphical_abstract.R (page 1) and _p2.R; same palette and conventions. Rebuilt 19 Nov from the
# SHORT pre-read's scoreboard code (analysis/fig_graphical_abstract_short_p3.R), with the long figure's extra columns kept:
# the stacked variants of the lead method (CSV keys a3r2, a3r3: two and three layers; exploratory, not adopted).
# Names: aCLR (CSV key "Attention", earlier text "attention") and tcaCLR (CSV key "a3r", earlier text "a3r"); both defined below.
# Rows: the 19 rows of the earlier long figure plus the two rows of the later preregistered runs (HepG2 and Kemmeren),
# read from docs/tcaclr_confirm_results/{armA,armB}/metrics.csv (seven methods: CLR, aCLR, tcaCLR, GENIE3, PIDC, GRNBoost2, |r|).
# The stacked variants were not run on the new rows: dashes there.
# Numbers: docs/whitepapers/figures/scoreboard_data.csv (built by analysis/scoreboard_matrix.R); nothing is a test:
# no intervals, no p-values, no statistical winner. A bold value with a dot is the highest printed value of its row
# (equal printed values share the mark and the tally). Two pages:
#   graphical_abstract_p3   raw scores, with a tally row
#   graphical_abstract_p3b  hub-controlled versions (within-stratum percentile in 10 strata of pi_i + pi_j) of the rows that have them
#
# Output: docs/whitepapers/figures/graphical_abstract_p3{,b}.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_p3.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")
out_dir <- Sys.getenv("SB_OUT", "docs/whitepapers/figures")           # SB_OUT / SB_CSV / SB_VERBOSE: only for re-running on another copy of the data (checks)
FONT <- "Liberation Sans"
LEAD_NAME <- "tcaCLR"                     # DEFINED HERE ONLY (temporary name; internal name a3r, earlier text "a3r")
ACLR_NAME <- "aCLR"                       # DEFINED HERE ONLY (temporary name; internal name Attention, earlier text "attention")
NAMES_NOTE <- "Method names are temporary."
LEAD_KEY <- "a3r"
LEAD_STATUS <- "lead candidate"
GREEN <- pal[["green"]]
D <- utils::read.csv(Sys.getenv("SB_CSV", file.path(out_dir, "scoreboard_data.csv")), stringsAsFactors = FALSE)
# Two rows added after the tcaCLR confirmatory runs (docs/tcaclr_confirm_results/{armA,armB}/metrics.csv: pooled AUPR and
# base rate of the primary truth, raw scores). Scaled like the other complex rows: pooled AUPR / base rate. Not in scoreboard_data.csv.
CONF_DIR <- Sys.getenv("SB_CONF", "docs/tcaclr_confirm_results")
conf_rows <- function(arm, ds) {
  m <- utils::read.csv(file.path(CONF_DIR, arm, "metrics.csv"), stringsAsFactors = FALSE)
  m <- m[m$truth == "main" & m$scope == "raw", ]
  key <- c(tcaCLR = "a3r", aCLR = "Attention", CLR = "CLR", GENIE3 = "GENIE3", PIDC = "PIDC", GRNBoost2 = "GRNBoost2", absr = "abs_r")
  data.frame(dataset = ds, truth = "complexes_main", metric = "pooled_aupr_over_base", scope = "raw", method = unname(key[m$method]),
             enrichment_over_random = m$aupr / m$base_rate, base_rate = m$base_rate,
             source_file = file.path(CONF_DIR, arm, "metrics.csv"), status = "from_existing", note = "tcaCLR confirmatory run", stringsAsFactors = FALSE)
}
D <- rbind(D, conf_rows("armA", "hepg2"), conf_rows("armB", "kemmeren"))

txt <- function(label, x, y, size = 10, col = ink[["primary"]], face = "plain", just = "centre", lineheight = 1.0)
  grid.text(label, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face, lineheight = lineheight))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1)
  grid.segments(x0, y0, x1, y1, default.units = "native", gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))

# ---- columns ----------------------------------------------------------------------------------------------------------
meths <- c("Attention", "CLR", LEAD_KEY, "a3r2", "a3r3", "GENIE3", "PIDC", "GRNBoost2", "abs_r")
head1 <- c(ACLR_NAME, "CLR", LEAD_NAME, LEAD_NAME, LEAD_NAME, "GENIE3", "PIDC", "GRNBoost2", "|r|")
head2 <- c("adopted", "", "lead candidate", "stacked, 2 layers", "stacked, 3 layers", "", "", "", "")
head3 <- c("", "", "", "exploratory", "exploratory", "", "", "", "")
EXPL <- "#bfe6d8"
fillc <- c(Attention = pal[["aqua"]], CLR = pal[["orange"]], EXPL, a3r2 = EXPL, a3r3 = EXPL,
           GENIE3 = ink[["muted"]], PIDC = ink[["muted"]], GRNBoost2 = ink[["muted"]], abs_r = ink[["muted"]])
names(fillc)[3] <- LEAD_KEY
EXPL_KEYS <- c(LEAD_KEY, "a3r2", "a3r3")
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
  R("Median AUPRC ratio, cell-type ChIP truth", "beeline", "all_14_cases", "beeline_median_auprc_ratio_celltypeChIP"),
  list(group = "Yeast, DREAM5 net4 expression (536 chips)"),
  R("Protein complexes (Complex Portal, exploratory estimation)", "dream5_yeast", "complexes_main", PM),
  list(group = "Human HepG2 Perturb-seq pseudobulk (later preregistered run, Complex Portal human complexes)"),
  R("Complex Portal human complexes", "hepg2", "complexes_main", PM),
  list(group = "Yeast, Kemmeren deletion compendium (1,283 profiles; estimate only)"),
  R("CYC2008 complexes", "kemmeren", "complexes_main", PM))
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
draw_page <- function(rows, sub, y_top, footer, y_bot = 1.43, title = TITLE, tally = FALSE) {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, 16), yscale = c(0, 9)))
  txt(title, 0.3, 8.65, 21, face = "bold", just = "left")
  txt(sub, 0.3, 8.28, 11, ink[["secondary"]], just = "left")
  nr <- sum(vapply(rows, function(r) is.null(r$group), NA)); ng <- length(rows) - nr
  nrows <- nr; if (tally) nr <- nr + 1          # the tally row takes one row pitch
  firsts <- stats::setNames(integer(length(meths)), meths)
  GH <- 0.74                                   # group-header height relative to a row
  pitch <- (y_top - y_bot) / (nr + ng * GH)
  # column bands (adopted column; the lead method and its stacked variants) and headers
  hy <- y_top + 0.5
  grid.rect(cx0(1), y_bot - 0.05, CW, y_top - y_bot + 0.78, default.units = "native", just = c("left", "bottom"), gp = gpar(fill = "#eaf7f1", col = NA))
  grid.rect(cx0(3), y_bot - 0.05, 3 * CW, y_top - y_bot + 0.78, default.units = "native", just = c("left", "bottom"), gp = gpar(fill = "#f6fcf9", col = pal[["aqua"]], lty = 2, lwd = 1))
  for (j in seq_along(meths)) {
    x <- cx0(j) + CW / 2
    txt(head1[j], x, hy + 0.1, 11, face = "bold")
    if (nzchar(head2[j])) txt(head2[j], x, hy - 0.1, 9.5, if (j == 1) "#0d6e4c" else ink[["secondary"]], "italic")
    if (nzchar(head3[j])) txt(head3[j], x, hy - 0.28, 9.5, ink[["secondary"]], "italic")
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
    v <- cell_vals(r)
    # first or shared first: highest value, or tied with it at the printed precision
    best <- if (all(is.na(v))) integer(0) else which(!is.na(v) & vapply(v, function(z) !is.na(z) && fmt(z, r$digits) == fmt(max(v, na.rm = TRUE), r$digits), NA))
    firsts[best] <- firsts[best] + 1L
    if (nzchar(Sys.getenv("SB_VERBOSE")) && tally) cat(sprintf("ROW %-62s | %s\n", paste(r$ds, r$truth, r$metric), paste(names(v)[best], collapse = ",")))
    span <- max(v - r$base, na.rm = TRUE); if (!is.finite(span) || span <= 0) span <- 1
    txt(r$label, X_LAB + 0.12, y, 10, just = "left")
    for (j in seq_along(meths)) {
      x0 <- cx0(j)
      seg(x0 + ORG_X, y - pitch * 0.42, x0 + ORG_X, y + pitch * 0.42, ink[["secondary"]], 0.8)          # random
      if (is.na(v[j])) { txt("–", x0 + NUM_X + 0.1, y, 11, ink[["muted"]]); next }
      frac <- (v[j] - r$base) / span; len <- if (frac >= 0) frac * BAR_W else max(frac, -0.08 / BAR_W) * BAR_W
      expl <- meths[j] %in% EXPL_KEYS
      bh <- pitch * 0.5
      grid.rect(x0 + ORG_X, y, abs(len), bh, default.units = "native", just = c(if (len >= 0) "left" else "right", "centre"),
                gp = gpar(fill = fillc[[meths[j]]], col = if (expl) pal[["aqua"]] else NA, lwd = 1.1, lty = if (expl) 2 else 1))
      isb <- j %in% best
      if (isb) grid.circle(x0 + 0.07, y, 0.035, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
      txt(fmt(v[j], r$digits), x0 + NUM_X + 0.02, y, 10, face = if (isb) "bold" else "plain", just = "left")
    }
    y <- y - pitch / 2
  }
  if (tally) {
    seg(X_LAB, y, 15.75, y, ink[["secondary"]], 1.2)
    y <- y - pitch / 2
    txt(paste0("Rows with the highest value or tied for it (of ", nrows, ")"), X_LAB + 0.12, y, 10, face = "bold", just = "left")
    win <- firsts == max(firsts)
    if (nzchar(Sys.getenv("SB_VERBOSE"))) cat("TALLY", paste(names(firsts), firsts, sep = "=", collapse = " "), "of", nrows, "\n")
    for (j in seq_along(meths)) {
      x0 <- cx0(j)
      if (win[j]) grid.roundrect(x0 + NUM_X - 0.1, y, 1.0, pitch * 0.8, default.units = "native", just = c("left", "centre"),
                                 r = unit(0.05, "inches"), gp = gpar(fill = NA, col = ink[["primary"]], lwd = 1.4))
      txt(sprintf("%d/%d", firsts[[j]], nrows), x0 + NUM_X + 0.02, y, 11, face = if (win[j]) "bold" else "plain", just = "left")
    }
    y <- y - pitch / 2
  }
  # footer
  for (k in seq_along(footer)) txt(footer[[k]]$text, 0.3, y_bot - 0.22 - (k - 1) * (if (length(footer) > 7) 0.178 else 0.185), footer[[k]]$size, footer[[k]]$col, footer[[k]]$face, just = "left")
  popViewport()
}

fl <- function(t) list(text = t, size = 9, col = ink[["secondary"]], face = "plain")
foot_common <- list(
  fl("Bars start at random (tick; 1.00, or 0.50 for the AUROC row) and are scaled within each row; the value, not the bar, carries the comparison across rows."),
  fl("A bold value with a dot is the highest printed value in its row (equal printed values share the mark and the tally), not a statistical winner; no intervals or tests are shown. A grey dash is a cell that was not run, which cannot win."),
  fl(paste0("Pooled AUPR / base rate (1.00 = random), raw scores unless marked hub-controlled; ", LEAD_NAME, " was picked among about 20 read-outs by its K562 scores (selected, not fitted); K562 and M3D are tuning sets;")),
  fl(paste0(ACLR_NAME, " and ", LEAD_NAME, " settings were developed on M3D. ", LEAD_NAME, " is the lead candidate: its preregistered E. coli test did not reach the threshold; its preregistered HepG2 test (a separate result on a new truth) did. The stacked variants (2 and 3 layers) are exploratory and not adopted.")))
foot_raw <- c(foot_common, list(
  fl(paste0("DREAM5 yeast row: one exploratory estimation run (plan fixed in advance), not a confirmatory test; competitors from saved score matrices; ", LEAD_NAME, " was not run on the DREAM5 regulon or edge truths (except a descriptive line); the stacked variants were not run on it.")),
  fl(paste0("HepG2 and Kemmeren rows: later preregistered runs, one pass each, ", ACLR_NAME, " and ", LEAD_NAME, " not tuned on them; Kemmeren is estimate only; stacked variants not run. HepG2 shares its lab and screen family with K562; there ", LEAD_NAME, " is not shown better than CLR, |r|, GENIE3 or PIDC.")),
  fl("RPE1 row: GENIE3, PIDC and GRNBoost2 were run on its 8,745 genes (two duplicated-symbol rows dropped) and scored on the same truth. E. coli complex rows: competitors re-scored from the saved score matrices of the earlier comparison."),
  fl(paste0("BEELINE: median over the 44 case-truth cells of the primary truths (14 dataset-and-gene-set cases); ", ACLR_NAME, " and the ", LEAD_NAME, " columns use the CLR-NS operator (", ACLR_NAME, " at its held-out depth, ", LEAD_NAME, " family at t = 9, not tuned); CLR is plain CLR; second row cell-type ChIP only."))))
foot_hub <- c(foot_common, list(
  fl("The yeast complex row of the raw scoreboard is left out here: the study stored its hub-controlled scores only with the MI row-mean stratifier, not the stratifier used on this page, so they are not comparable."),
  fl(paste0("Hub-controlled: each method's scores are converted to percentiles within 10 strata of pi_i + pi_j (pi = stationary distribution of the ", ACLR_NAME, " walk), the stratifier of the method table. Not run for direct edges, per-regulon AP and BEELINE.")),
  fl("The competitors were scored on the complex rows (raw scores only) and not hub-controlled, so those cells are dashes here. The HepG2 and Kemmeren rows of the raw scoreboard (later runs) have no hub-controlled version and are not shown.")))

SUB_RAW <- paste0("Raw scores. ", ACLR_NAME, " (adopted, t = 9) and CLR, then the lead method ", LEAD_NAME, " and two exploratory stacked variants (dashed), then competitors.")
SUB_HUB <- "Hub-controlled versions of the same rows (page 3 shows the raw scores): each method's scores are compared only among gene pairs of similar hub degree."
TITLE <- "Every metric, wins and losses"
for (nm in c("", "b")) {
  rows <- if (nm == "") raw_rows else hub_rows
  sub <- paste(if (nm == "") SUB_RAW else SUB_HUB, NAMES_NOTE)
  ft <- if (nm == "") foot_raw else foot_hub
  png(file.path(out_dir, paste0("graphical_abstract_p3", nm, ".png")), width = 2400, height = 1350, res = 150, type = "cairo", bg = "white")
  draw_page(rows, sub, 7.3, ft, y_bot = if (nm == "") 1.58 else 1.43, title = if (nm == "") TITLE else "Hub-controlled scoreboard", tally = nm == ""); invisible(dev.off())
  cairo_pdf(file.path(out_dir, paste0("graphical_abstract_p3", nm, ".pdf")), width = 16, height = 9, bg = "white")
  draw_page(rows, sub, 7.3, ft, y_bot = if (nm == "") 1.58 else 1.43, title = if (nm == "") TITLE else "Hub-controlled scoreboard", tally = nm == ""); invisible(dev.off())
}
cat("wrote", file.path(out_dir, "graphical_abstract_p3{,b}.{png,pdf}"), "\n")
