#!/usr/bin/env Rscript
# Summary of stage 1 (plan Revision 2): per dataset, method and criterion, the argmax cell of the median over splits (directions averaged within a
# split), the plateau around it, ties, the comparison with the registered cell, and a DESCRIPTIVE agreement report. The earlier "C1, C2, C3 within one
# grid step" declare-or-not rule was withdrawn by Revision 2 and no new rule is applied here (that is a later step): the agreement table says only how far
# apart the criteria's cells are, the surface selection and one-standard-error rule are not implemented.
#   Rscript analysis/tcaclr_depth/stage1_summary.R --files results/tcaclr_depth/stage1_m3d/stage1_M3D.csv,results/tcaclr_depth/stage1_k562/stage1_K562.csv --out results/tcaclr_depth/stage1_summary
# Criteria summarised: C2_adj (null-adjusted reproducibility; C2 raw and C2_null shown), C3 (blocked held-out likelihood, eps 0.01), C3_e10, C3_bic (second reading),
# and C1 (diagnostic only: similarity to CLR). C3 does not depend on the method or on clipO: it is taken once per operator and depth.
suppressPackageStartupMessages(source("analysis/tcaclr_depth/grid_lib.R"))
opt <- gl_options(commandArgs(TRUE), list(files = "", out = "results/tcaclr_depth/stage1_summary"))
stopifnot(nzchar(opt$files)); dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
D <- do.call(rbind, lapply(strsplit(opt$files, ",")[[1]], utils::read.csv, stringsAsFactors = FALSE))
D$clipS <- as.logical(D$clipS); D$clipO <- ifelse(is.na(D$clipO), "NA", as.character(as.logical(D$clipO)))   # aggregate() drops NA groups: code NA as "NA"
crit <- c("C2_adj", "C3", "C3_e10", "C3_bic", "C1", "C2", "C2_null")
selecting <- c("C2_adj", "C3")                       # the two criteria that carry selection; C3_e10 / C3_bic are readings of C3; C1 is diagnostic only
key <- function(d) paste(d$method, d$keff, d$cap, d$clipS, d$clipO, d$t, sep = "|")
# per split: mean of the two directions; then median over splits
A <- stats::aggregate(D[, crit], by = list(dataset = D$dataset, split = D$split, method = D$method, keff = D$keff, cap = D$cap, clipS = D$clipS, clipO = D$clipO, t = D$t), FUN = mean)
M <- stats::aggregate(A[, crit], by = A[, c("dataset", "method", "keff", "cap", "clipS", "clipO", "t")], FUN = stats::median)
M$clipO <- ifelse(M$clipO == "NA", NA, M$clipO == "TRUE")
nsp <- tapply(A$split, A$dataset, function(x) length(unique(x)))
isref <- M$cap == GL$ref$cap & M$keff == GL$ref$keff & M$clipS & M$t == GL$ref$t & !(M$cap == gl_cap(M$keff))   # reference operator (cap 50) vs the grid's keff 10 (cap 40)
G_ <- M[!isref, ]; ref <- M[isref, ]
tg <- GL$depths; kg <- GL$keffs
half_oct <- function(t, k) c(depth = abs(log2(t / 9)) >= 0.5, width = abs(log2(k / 10)) >= 0.5)
arg <- list(); plat <- list(); cmp <- list(); agree <- list()
for (ds in unique(M$dataset)) for (me in c("tcaCLR", "aCLR")) for (cr in crit) {
  x <- G_[G_$dataset == ds & G_$method == me, ]
  if (cr %in% c("C3", "C3_e10", "C3_bic")) { x <- G_[G_$dataset == ds & G_$method == "aCLR", ]; x$clipO <- NA }   # C3 is method- and clipO-free: one value per operator and depth
  v <- x[[cr]]; mx <- max(v, na.rm = TRUE); ties <- which(abs(v - mx) <= 1e-12 * max(1, abs(mx)))
  b <- x[ties[order(x$t[ties], x$keff[ties], !x$clipS[ties])[1L]], ]       # ties: shallowest, then narrowest, then clipped
  ho <- half_oct(b$t, b$keff)
  arg[[length(arg) + 1L]] <- data.frame(dataset = ds, method = me, criterion = cr, role = if (cr %in% selecting) "selecting" else if (cr == "C1") "diagnostic" else "reading",
    keff = b$keff, cap = b$cap, clipS = b$clipS, clipO = b$clipO, t = b$t, median = b[[cr]], n_ties = length(ties), n_splits = nsp[[ds]],
    differs_from_registered_by_half_octave_or_more_in = paste(c("depth", "width")[ho], collapse = "+"), stringsAsFactors = FALSE)
  # plateau: the cell and its nearest neighbours in each axis (medians), and the width along t / keff of the contiguous region within 1% of the maximum
  same <- function(d) x[x$keff == d$keff & x$clipS == d$clipS & (is.na(d$clipO) | x$clipO %in% d$clipO), ]
  ti <- match(b$t, tg); ki <- match(b$keff, kg); pl <- list()
  add <- function(axis, step, k_, t_, cs, co) { r <- x[x$keff == k_ & x$t == t_ & x$clipS == cs & (is.na(co) | x$clipO %in% co), ]
    if (nrow(r)) pl[[length(pl) + 1L]] <<- data.frame(dataset = ds, method = me, criterion = cr, axis = axis, step = step, keff = k_, t = t_, clipS = cs, clipO = co, median = r[[cr]][1L]) }
  for (s in -1:1) { if (ti + s >= 1 && ti + s <= length(tg)) add("t", s, b$keff, tg[ti + s], b$clipS, b$clipO); if (ki + s >= 1 && ki + s <= length(kg)) add("keff", s, kg[ki + s], b$t, b$clipS, b$clipO) }
  add("clipS", 1, b$keff, b$t, !b$clipS, b$clipO); if (!is.na(b$clipO)) add("clipO", 1, b$keff, b$t, b$clipS, !b$clipO)
  plat <- c(plat, pl)
  thr <- mx - 0.01 * abs(mx); run <- function(vals, i) { lo <- i; while (lo > 1 && !is.na(vals[lo - 1]) && vals[lo - 1] >= thr) lo <- lo - 1; hi <- i; while (hi < length(vals) && !is.na(vals[hi + 1]) && vals[hi + 1] >= thr) hi <- hi + 1; hi - lo + 1 }
  vt <- vapply(tg, function(tt) { r <- x[x$keff == b$keff & x$t == tt & x$clipS == b$clipS & (is.na(b$clipO) | x$clipO %in% b$clipO), cr]; if (length(r)) r[1L] else NA_real_ }, 0)
  vk <- vapply(kg, function(kk) { r <- x[x$keff == kk & x$t == b$t & x$clipS == b$clipS & (is.na(b$clipO) | x$clipO %in% b$clipO), cr]; if (length(r)) r[1L] else NA_real_ }, 0)
  arg[[length(arg)]]$plateau_width_t_steps_within_1pct <- run(vt, ti); arg[[length(arg)]]$plateau_width_keff_steps_within_1pct <- run(vk, ki)
  # comparison with the registered cell (cap 50, keff 10, t 9, clipped both) and its nearest grid cell (cap 40)
  rg <- ref[ref$dataset == ds & ref$method == me, ]; if (!is.na(b$clipO) || me == "aCLR") rg <- rg[rg$clipO %in% c(TRUE, NA), ]
  ng <- x[x$keff == 10 & x$t == 9 & x$clipS & (is.na(x$clipO) | x$clipO %in% TRUE), ]
  cmp[[length(cmp) + 1L]] <- data.frame(dataset = ds, method = me, criterion = cr, argmax_value = b[[cr]], registered_cell_value = if (nrow(rg)) rg[[cr]][1L] else NA_real_,
    nearest_grid_cell_value = if (nrow(ng)) ng[[cr]][1L] else NA_real_, argmax_minus_registered = b[[cr]] - if (nrow(rg)) rg[[cr]][1L] else NA_real_)
}
arg <- do.call(rbind, arg); plat <- do.call(rbind, plat); cmp <- do.call(rbind, cmp)
## descriptive agreement of the two selecting criteria (and the C1 diagnostic) per dataset and method
for (ds in unique(arg$dataset)) for (me in c("tcaCLR", "aCLR")) {
  a <- arg[arg$dataset == ds & arg$method == me & arg$criterion %in% c("C2_adj", "C3", "C1"), ]; if (nrow(a) < 2L) next
  st <- function(x, y) c(t_steps = abs(match(x$t, tg) - match(y$t, tg)), keff_steps = abs(match(x$keff, kg) - match(y$keff, kg)), clipS_differs = x$clipS != y$clipS, clipO_differs = isTRUE(x$clipO != y$clipO))
  for (pr in list(c("C2_adj", "C3"), c("C2_adj", "C1"), c("C3", "C1"))) { x <- a[a$criterion == pr[1], ]; y <- a[a$criterion == pr[2], ]; if (!nrow(x) || !nrow(y)) next
    agree[[length(agree) + 1L]] <- data.frame(dataset = ds, method = me, criterion_a = pr[1], criterion_b = pr[2], t(st(x, y)), cell_a = sprintf("t%d keff%g clipS%s clipO%s", x$t, x$keff, x$clipS, x$clipO), cell_b = sprintf("t%d keff%g clipS%s clipO%s", y$t, y$keff, y$clipS, y$clipO)) }
}
agree <- do.call(rbind, agree)
gl_wcsv(arg, file.path(opt$out, "stage1_argmax.csv")); gl_wcsv(plat, file.path(opt$out, "stage1_plateau.csv")); gl_wcsv(cmp, file.path(opt$out, "stage1_registered_comparison.csv"))
gl_wcsv(agree, file.path(opt$out, "stage1_agreement_descriptive.csv")); gl_wcsv(M, file.path(opt$out, "stage1_cell_medians.csv"))
# the clip switches, read separately: median over cells of (unclipped - clipped) on each criterion, per dataset and method
cl <- list()
for (ds in unique(M$dataset)) for (me in c("tcaCLR", "aCLR")) for (cr in c("C2_adj", "C3", "C1")) {
  x <- G_[G_$dataset == ds & G_$method == me, ]; if (cr == "C3") x <- G_[G_$dataset == ds & G_$method == "aCLR", ]
  mk <- function(x) paste(x$keff, x$t, if (me == "tcaCLR" && cr != "C3") x$clipO else NA)
  a <- x[x$clipS, ]; b <- x[!x$clipS, ]; if (me == "tcaCLR" && cr != "C3") { a <- a[a$clipO %in% TRUE, ]; b <- b[b$clipO %in% TRUE, ] }
  m <- merge(a[, c("keff", "t", cr)], b[, c("keff", "t", cr)], by = c("keff", "t"), suffixes = c("_c", "_u"))
  cl[[length(cl) + 1L]] <- data.frame(dataset = ds, method = me, criterion = cr, switch = "clipS", n_cells = nrow(m), median_unclipped_minus_clipped = stats::median(m[[paste0(cr, "_u")]] - m[[paste0(cr, "_c")]]), share_unclipped_higher = mean(m[[paste0(cr, "_u")]] > m[[paste0(cr, "_c")]]))
  if (me == "tcaCLR" && cr != "C3") { a <- x[x$clipO %in% TRUE & x$clipS, ]; b <- x[x$clipO %in% FALSE & x$clipS, ]
    m <- merge(a[, c("keff", "t", cr)], b[, c("keff", "t", cr)], by = c("keff", "t"), suffixes = c("_c", "_u"))
    cl[[length(cl) + 1L]] <- data.frame(dataset = ds, method = me, criterion = cr, switch = "clipO (clipS clipped)", n_cells = nrow(m), median_unclipped_minus_clipped = stats::median(m[[paste0(cr, "_u")]] - m[[paste0(cr, "_c")]]), share_unclipped_higher = mean(m[[paste0(cr, "_u")]] > m[[paste0(cr, "_c")]])) }
}
cl <- do.call(rbind, cl); gl_wcsv(cl, file.path(opt$out, "stage1_clip_switches.csv"))
op <- options(width = 220, digits = 4)
cat("\n== argmax cell of the median over splits\n"); print(arg, row.names = FALSE)
cat("\n== registered cell (keff 10, cap 50, t 9, clipped) and nearest grid cell (cap 40)\n"); print(cmp, row.names = FALSE)
cat("\n== descriptive agreement of the criteria's cells (no rule applied; the Revision 1 declare-or-not rule is withdrawn)\n"); print(agree, row.names = FALSE)
cat("\n== clip switches (unclipped minus clipped, median over cells)\n"); print(cl, row.names = FALSE); options(op)
