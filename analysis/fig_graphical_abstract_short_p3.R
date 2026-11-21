# Scoreboard for the SHORT pre-read (portrait, 8.3 in wide canvas, printed 7.0 in wide): "Every metric, wins and losses". A descriptive scoreboard, drawn in code with grid.
# Copy of analysis/fig_graphical_abstract_p3.R, revised 15 Nov. Columns kept: CLR (historical), CLR attention, the lead
# method, GENIE3, PIDC, GRNBoost2, |r|. Columns dropped: a3r2 and a3r3 (the stacked variants). All metric rows kept
# (no method or truth is dropped other than those two columns). The best value of a row is recomputed over the kept columns.
# Numbers: docs/whitepapers/figures/scoreboard_data.csv (built by analysis/scoreboard_matrix.R); nothing is a test.
# Two pages: graphical_abstract_short_p3 (raw, main body) and graphical_abstract_short_p3b (hub-controlled, appendix).
# The lead method's name is defined once, in LEAD_NAME (the CSV keeps its internal key "a3r"; a rename is a one-line change here).
# Dotted green = new in this work (human-data groups).
#
# Output: docs/whitepapers/figures/graphical_abstract_short_p3{,b}.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_graphical_abstract_short_p3.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")
out_dir <- Sys.getenv("SB_OUT", "docs/whitepapers/figures"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)           # SB_OUT / SB_CSV / SB_VERBOSE: only for re-running on another copy of the data (checks)
FONT <- "Liberation Sans"
LEAD_NAME <- "tcaCLR"                     # DEFINED HERE ONLY (temporary name; internal name a3r)
ACLR_NAME <- "aCLR"                       # DEFINED HERE ONLY (temporary name)
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
# tcaCLR2 cells (exploratory): same schema as scoreboard_data.csv, method key "tcaCLR2"
T2_KEY <- "tcaCLR2"; T2_NAME <- "tcaCLR2"
t2c <- utils::read.csv(Sys.getenv("SB_T2", "results/tcaclr2_scoreboard/cells_tcaclr2_all.csv"), stringsAsFactors = FALSE)
D <- rbind(D, t2c[t2c$method == T2_KEY, names(D)])
T2_DEPTHS <- utils::read.csv(Sys.getenv("SB_T2D", "results/tcaclr2_scoreboard/t2_by_dataset.csv"), stringsAsFactors = FALSE)   # label, t2 (assemble_tcaclr2_scoreboard.R)
T2_NOTE1 <- paste0(T2_NAME, " is exploratory and outside the preregistration: ", LEAD_NAME, " with the depth set by the data (t2 per data set: ", paste0(T2_DEPTHS$label, " ", T2_DEPTHS$t2, collapse = ", "), "), not t = 9.")
T2_NOTE2 <- paste0("It carries no statistical-tie marks (its mark, like the others, is the highest printed value of the row), and it is lower than ", LEAD_NAME, " on several non-pooled metrics (TF-edge AUPR, per-TF AUROC, per-regulon AP, RPE1, BEELINE medians).")
AUDIT <- list(); AUDIT_DIR <- Sys.getenv("SB_AUDIT", "results/tcaclr2_scoreboard/figures"); dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)   # audit csv files (marks with and without tcaCLR2)

txt <- function(label, x, y, size = 10, col = ink[["primary"]], face = "plain", just = "centre", lineheight = 1.0)
  grid.text(label, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face, lineheight = lineheight))
seg <- function(x0, y0, x1, y1, col = ink[["muted"]], lwd = 1, lty = 1)
  grid.segments(x0, y0, x1, y1, default.units = "native", gp = gpar(col = col, lwd = lwd, lty = lty, lineend = "butt"))

# ---- columns (PORTRAIT: canvas W_CAN x H_CAN inches, printed at 6.9 in wide, scale 6.9 / W_CAN) ----------------------------------------
W_CAN <- 8.3                    # canvas width in inches; printed width 6.9 in, so every font size here prints at 0.83x (>= 8.6 pt -> >= 7.1 pt)
meths <- c("CLR", "Attention", LEAD_KEY, T2_KEY, "GENIE3", "PIDC", "GRNBoost2", "abs_r")
head1 <- c("CLR", ACLR_NAME, LEAD_NAME, T2_NAME, "GENIE3", "PIDC", "GRNBoost2", "|r|")
head2 <- c("historical", "adopted", "lead\ncandidate", "depth set\nby data", "", "", "", "")
EXPL <- "#bfe6d8"
fillc <- c(CLR = pal[["orange"]], Attention = pal[["aqua"]], EXPL, EXPL, GENIE3 = ink[["muted"]], PIDC = ink[["muted"]],
           GRNBoost2 = ink[["muted"]], abs_r = ink[["muted"]])
names(fillc)[3] <- LEAD_KEY; names(fillc)[4] <- T2_KEY
EXPL_KEYS <- c(LEAD_KEY, T2_KEY)
X_LAB <- 0.12                   # row labels start
X_R <- W_CAN - 0.12             # right edge of the table
X_C0 <- 2.45                    # first method column starts
CW <- (X_R - X_C0) / 8          # column width
NUM_X <- 0.09                   # number, left edge inside a cell
ORG_X <- 0.44                   # bar origin (the 'random' tick) inside a cell
BAR_W <- CW - ORG_X - 0.03      # longest bar
cx0 <- function(j) X_C0 + (j - 1) * CW
FS_ROW <- 8.8; FS_NUM <- 9; FS_GRP <- 9.3; FS_HEAD <- 9.2; FS_FOOT <- 8.6; LAB_W <- X_C0 - X_LAB - 0.1 - 0.08   # font sizes (canvas pt) and label width
text_w <- function(s, size, face = "plain") convertWidth(grobWidth(textGrob(s, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face))), "inches", TRUE)
wrap_to <- function(s, size, width, face = "plain") {     # word wrap (never inside a word) to a width in inches
  words <- strsplit(s, " ", fixed = TRUE)[[1]]; lines <- character(); cur <- ""
  for (w in words) { cand <- if (nzchar(cur)) paste(cur, w) else w
    if (!nzchar(cur) || text_w(cand, size, face) <= width) cur <- cand else { lines <- c(lines, cur); cur <- w } }
  c(lines, cur) }

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
  list(group = "Human HepG2 Perturb-seq pseudobulk (later preregistered run)"),
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

# ---- one page (portrait) ------------------------------------------------------------------------------------------------
# Two passes: pass 1 (H = 40) only measures the height used, pass 2 draws on a canvas of exactly that height. y runs upward; T = top of the canvas.
draw_page <- function(rows, sub, footer, H, title = TITLE, tally = FALSE) {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, W_CAN), yscale = c(0, H)))
  T <- H
  txt(title, X_LAB, T - 0.3, 17, face = "bold", just = "left")
  sl <- wrap_to(sub, 9.6, X_R - X_LAB)
  for (k in seq_along(sl)) txt(sl[k], X_LAB, T - 0.62 - (k - 1) * 0.17, 9.6, ink[["secondary"]], just = "left")
  y_hdr <- T - 0.62 - (length(sl) - 1) * 0.17 - 0.52          # column header baseline row (name); sub-labels below it
  nr <- sum(vapply(rows, function(r) is.null(r$group), NA)); ng <- length(rows) - nr
  nrows <- nr
  firsts <- stats::setNames(integer(length(meths)), meths)
  RH <- 0.19; GHH <- 0.185                                       # row height (one label line), group-header height
  LH <- 0.14                                                   # line height of a wrapped row label
  rh <- function(r) { if (!is.null(r$group)) return(GHH); nl <- length(wrap_to(r$label, FS_ROW, LAB_W)); max(RH, 0.06 + LH * nl) }
  y_top <- y_hdr - 0.36
  tot <- sum(vapply(rows, rh, 0)) + (if (tally) max(RH, 0.06 + LH * length(wrap_to(paste0("Rows with the highest value or tied for it (of ", nrows, ")"), FS_ROW, LAB_W, "bold"))) else 0)
  y_bot <- y_top - tot
  # column bands (adopted and exploratory columns) and headers
  grid.rect(cx0(2), y_bot - 0.04, CW, y_hdr + 0.4 - y_bot + 0.04, default.units = "native", just = c("left", "bottom"), gp = gpar(fill = "#eaf7f1", col = NA))
  grid.rect(cx0(3), y_bot - 0.04, 2 * CW, y_hdr + 0.4 - y_bot + 0.04, default.units = "native", just = c("left", "bottom"), gp = gpar(fill = "#f6fcf9", col = pal[["aqua"]], lty = 2, lwd = 1))
  for (j in seq_along(meths)) {
    x <- cx0(j) + CW / 2
    txt(head1[j], x, y_hdr + 0.2, FS_HEAD, face = "bold")
    if (nzchar(head2[j])) txt(head2[j], x, y_hdr - 0.04, 8.6, if (j == 2) "#0d6e4c" else ink[["secondary"]], "italic", lineheight = 0.95)
  }
  seg(X_LAB, y_top + 0.03, X_R, y_top + 0.03, ink[["secondary"]], 1)
  y <- y_top
  for (r in rows) {
    h <- rh(r)
    if (!is.null(r$group)) {
      yc <- y - h * 0.5 - 0.02
      txt(r$group, X_LAB, yc, FS_GRP, ink[["primary"]], "bold", just = "left")
      if (grepl("^Human|^Single cells", r$group)) {   # human data: new in this work
        gw <- text_w(r$group, FS_GRP, "bold")
        grid.roundrect(X_LAB - 0.04, yc, gw + 0.08, 0.17, default.units = "native", just = c("left", "centre"),
                       r = unit(0.03, "inches"), gp = gpar(col = GREEN, lty = "11", lwd = 1.8, fill = NA)) }
      seg(X_LAB, y - h + 0.005, X_R, y - h + 0.005, ink[["grid"]], 0.8)
      y <- y - h; next }
    yc <- y - h / 2
    v <- cell_vals(r)
    # first or shared first: highest value, or tied with it at the printed precision
    best <- if (all(is.na(v))) integer(0) else which(!is.na(v) & vapply(v, function(z) !is.na(z) && fmt(z, r$digits) == fmt(max(v, na.rm = TRUE), r$digits), NA))
    firsts[best] <- firsts[best] + 1L
    if (is.null(r$label)) NULL else if (TRUE) {      # audit: the mark with and without the tcaCLR2 column (the original columns only)
      v0 <- v; v0[T2_KEY] <- NA_real_
      best0 <- if (all(is.na(v0))) integer(0) else which(!is.na(v0) & vapply(v0, function(z) !is.na(z) && fmt(z, r$digits) == fmt(max(v0, na.rm = TRUE), r$digits), NA))
      AUDIT[[length(AUDIT) + 1L]] <<- data.frame(row = r$label, dataset = r$ds, truth = r$truth, metric = r$metric, scope = r$scope, t(round(v, 6)),
        best_without_tcaCLR2 = paste(names(v)[best0], collapse = "+"), best_with_tcaCLR2 = paste(names(v)[best], collapse = "+"),
        tcaCLR2_in_best = T2_KEY %in% names(v)[best], tcaCLR2_takes_mark_from_tcaCLR = (LEAD_KEY %in% names(v)[best0]) && !(LEAD_KEY %in% names(v)[best]),
        tcaCLR2_shares_mark_with_tcaCLR = (LEAD_KEY %in% names(v)[best0]) && (T2_KEY %in% names(v)[best]) && (LEAD_KEY %in% names(v)[best]),
        tcaCLR2_takes_mark_from_other = (T2_KEY %in% names(v)[best]) && !(LEAD_KEY %in% names(v)[best0]) && !(T2_KEY %in% names(v)[best0]) && length(best0) > 0, stringsAsFactors = FALSE) }
    if (nzchar(Sys.getenv("SB_VERBOSE")) && tally) cat(sprintf("ROW %-62s | %s\n", paste(r$ds, r$truth, r$metric), paste(names(v)[best], collapse = ",")))
    span <- max(v - r$base, na.rm = TRUE); if (!is.finite(span) || span <= 0) span <- 1
    ll <- wrap_to(r$label, FS_ROW, LAB_W)
    for (k in seq_along(ll)) txt(ll[k], X_LAB + 0.1, yc + ((length(ll) - 1) / 2 - (k - 1)) * LH, FS_ROW, just = "left")
    for (j in seq_along(meths)) {
      x0 <- cx0(j)
      seg(x0 + ORG_X, yc - RH * 0.42, x0 + ORG_X, yc + RH * 0.42, ink[["secondary"]], 0.8)          # random
      if (is.na(v[j])) { txt("–", x0 + NUM_X + 0.1, yc, 10, ink[["muted"]]); next }
      frac <- (v[j] - r$base) / span; len <- if (frac >= 0) frac * BAR_W else max(frac, -0.04 / BAR_W) * BAR_W
      expl <- meths[j] %in% EXPL_KEYS
      bh <- RH * 0.5
      grid.rect(x0 + ORG_X, yc, abs(len), bh, default.units = "native", just = c(if (len >= 0) "left" else "right", "centre"),
                gp = gpar(fill = fillc[[meths[j]]], col = if (expl) pal[["aqua"]] else NA, lwd = 0.9, lty = if (expl) 2 else 1))
      isb <- j %in% best
      if (isb) grid.circle(x0 + 0.05, yc, 0.025, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
      txt(fmt(v[j], r$digits), x0 + NUM_X, yc, FS_NUM, face = if (isb) "bold" else "plain", just = "left")
    }
    y <- y - h
  }
  if (tally) {
    seg(X_LAB, y, X_R, y, ink[["secondary"]], 1.1)
    tl_ <- wrap_to(paste0("Rows with the highest value or tied for it (of ", nrows, ")"), FS_ROW, LAB_W, "bold")
    h <- max(RH, 0.06 + LH * length(tl_)); yc <- y - h / 2
    for (k in seq_along(tl_)) txt(tl_[k], X_LAB + 0.1, yc + ((length(tl_) - 1) / 2 - (k - 1)) * LH, FS_ROW, face = "bold", just = "left")
    win <- firsts == max(firsts)
    if (nzchar(Sys.getenv("SB_VERBOSE"))) cat("TALLY", paste(names(firsts), firsts, sep = "=", collapse = " "), "of", nrows, "\n")
    for (j in seq_along(meths)) {
      x0 <- cx0(j)
      txt(sprintf("%d/%d", firsts[[j]], nrows), x0 + NUM_X, yc, FS_NUM, face = if (win[j]) "bold" else "plain", just = "left")
    }
    y <- y - h
  }
  # footer (wrapped to the table width; one gap between notes)
  yf <- y_bot - 0.2
  for (k in seq_along(footer)) {
    fl_ <- wrap_to(footer[[k]]$text, FS_FOOT, X_R - X_LAB)
    for (m in seq_along(fl_)) { txt(fl_[m], X_LAB, yf, FS_FOOT, footer[[k]]$col, footer[[k]]$face, just = "left"); yf <- yf - 0.135 }
    yf <- yf - 0.03
  }
  popViewport()
  H - yf + 0.12            # height used (pass 1 returns it; pass 2 should equal H)
}

fl <- function(t) list(text = t, size = 9, col = ink[["secondary"]], face = "plain")
foot_common <- list(
  fl("Bars start at random (tick; 1.00, or 0.50 for the AUROC row) and are scaled within each row; the value, not the bar, carries the comparison across rows."),
  fl("A bold value with a dot is the highest printed value in its row (equal printed values share the mark and the tally), not a statistical winner; no intervals or tests are shown. A grey dash is a cell that was not run. Dotted green outline = human-data groups, new in this work."),
  fl(paste0("Pooled AUPR / base rate (1.00 = random), raw scores unless marked hub-controlled; ", LEAD_NAME, " was picked among about 20 read-outs by its K562 scores (selected, not fitted); K562 and M3D are tuning sets;")),
  fl(paste0(paste0(ACLR_NAME, " and "), LEAD_NAME, " settings were developed on M3D. ", LEAD_NAME, " is the lead candidate: its preregistered E. coli test did not reach the threshold; its preregistered HepG2 test (a separate result on a new truth) did.")))
foot_raw <- c(foot_common, list(fl(T2_NOTE1), fl(T2_NOTE2),
  fl(paste0("DREAM5 yeast row: one exploratory estimation run (plan fixed in advance), not a confirmatory test; competitors from saved score matrices; ", LEAD_NAME, " was not run on the DREAM5 regulon or edge truths (except a descriptive line).")),
  fl(paste0("HepG2 and Kemmeren rows: later preregistered runs, one pass each; ", ACLR_NAME, " and ", LEAD_NAME, " were not tuned on them (Kemmeren is estimate only). HepG2 shares its laboratory and screen family with K562; there ", LEAD_NAME, " is not shown better than CLR, |r|, GENIE3 or PIDC.")),
  fl("RPE1 row: GENIE3, PIDC and GRNBoost2 were run on its 8,745 genes (two duplicated-symbol rows dropped) and scored on the same truth. E. coli complex rows: competitors re-scored from the saved score matrices of the earlier comparison."),
  fl(paste0(paste0("BEELINE: median over the 44 case-truth cells of the primary truths (14 dataset-and-gene-set cases); ", ACLR_NAME, " and "), LEAD_NAME, paste0(" use the CLR-NS operator (", ACLR_NAME, " at its held-out depth, "), LEAD_NAME, " at t = 9, not tuned; ", T2_NAME, " at its own t2); CLR is plain CLR."))))
foot_hub <- c(foot_common, list(fl(T2_NOTE1), fl(T2_NOTE2),
  fl("The yeast complex row of the raw scoreboard is left out here: the study stored its hub-controlled scores only with the MI row-mean stratifier, not the stratifier used on this page, so they are not comparable."),
  fl(paste0("Hub-controlled: each method's scores are converted to percentiles within 10 strata of pi_i + pi_j (pi = stationary distribution of the ", ACLR_NAME, " walk). Not run for direct edges, per-regulon AP and BEELINE.")),
  fl("The competitors were scored on the complex rows (raw scores only) and not hub-controlled, so those cells are dashes here. The HepG2 and Kemmeren rows of the raw scoreboard (later runs) have no hub-controlled version and are not shown.")))

SUB_RAW <- paste0("Raw scores. CLR (historical), ", ACLR_NAME, " (adopted) and the lead method (dashed), then competitors.")
SUB_HUB <- "Appendix. Hub control: each method's scores are compared only among gene pairs of similar hub degree (how strongly connected the genes are). Raw scores are in the main body."
TITLE <- "Every metric, wins and losses"
measure_dev <- png(tempfile(fileext = ".png"), width = 800, height = 800, type = "cairo")   # for text measurements and the height pass
for (nm in c("", "b")) {
  rows <- if (nm == "") raw_rows else hub_rows
  sub <- paste(if (nm == "") SUB_RAW else SUB_HUB, NAMES_NOTE)
  ft <- if (nm == "") foot_raw else foot_hub
  ttl <- if (nm == "") TITLE else "Hub-controlled scoreboard"
  AUDIT <- list()
  H_CAN <- draw_page(rows, sub, ft, 40, title = ttl, tally = nm == "")      # pass 1: height used
  dev.off()
  AUDIT <- list()
  png(file.path(out_dir, paste0("graphical_abstract_short_p3", nm, ".png")), width = round(W_CAN * 300), height = round(H_CAN * 300), res = 300, type = "cairo", bg = "white")
  draw_page(rows, sub, ft, H_CAN, title = ttl, tally = nm == ""); invisible(dev.off())
  cat(sprintf("p3%s: canvas %.2f x %.2f in; printed %.2f x %.2f in (scale %.3f); min font %.1f pt -> %.1f pt printed\n", nm, W_CAN, H_CAN, 6.9, H_CAN * 6.9 / W_CAN, 6.9 / W_CAN, 8.6, 8.6 * 6.9 / W_CAN))
  au <- do.call(rbind, AUDIT); utils::write.csv(au, file.path(AUDIT_DIR, paste0("audit_best_marks_short_p3", nm, ".csv")), row.names = FALSE)
  cnt <- function(col) vapply(meths, function(m) sum(vapply(strsplit(au[[col]], "+", fixed = TRUE), function(z) m %in% z, NA)), 0L)
  tl <- data.frame(method = meths, rows_marked_without_tcaCLR2 = cnt("best_without_tcaCLR2"), rows_marked_with_tcaCLR2 = cnt("best_with_tcaCLR2")); tl$rows_marked_without_tcaCLR2[meths == T2_KEY] <- NA
  utils::write.csv(tl, file.path(AUDIT_DIR, paste0("audit_tally_short_p3", nm, ".csv")), row.names = FALSE)
  cairo_pdf(file.path(out_dir, paste0("graphical_abstract_short_p3", nm, ".pdf")), width = W_CAN, height = H_CAN, bg = "white")
  draw_page(rows, sub, ft, H_CAN, title = ttl, tally = nm == ""); invisible(dev.off())
  png(tempfile(fileext = ".png"), width = 800, height = 800, type = "cairo")
}
invisible(dev.off())
cat("wrote", file.path(out_dir, "graphical_abstract_short_p3{,b}.{png,pdf}"), "and the audit csv files\n")
