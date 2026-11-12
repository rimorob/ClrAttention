#!/usr/bin/env Rscript
# EXPLORATORY. Tables and literal pass/fail for ideas 4 and 7 from the readout caches
# (docs/plans/explore_ideas47_plan.md sections 1-2). Usage: Rscript ideas47_summarize.R <idea 4|7> <out dir>
args <- commandArgs(TRUE); idea <- as.integer(args[1]); out <- args[2]
cache <- file.path(out, "cache"); fs <- list.files(cache, "^readout__.*\\.rds$", full.names = TRUE)
R <- stats::setNames(lapply(fs, readRDS), sub("^readout__(.*)\\.rds$", "\\1", basename(fs)))
meta <- readRDS(file.path(cache, "meta.rds")); br <- meta$base_rate; TR <- meta$truths
if (idea == 7) { mp <- utils::read.csv(file.path(out, "readout_map.csv"), stringsAsFactors = FALSE) } else {
  mp <- data.frame(name = names(R), canonical = names(R), stringsAsFactors = FALSE) }
mp <- mp[mp$canonical %in% names(R), ]
refs <- if (idea == 4) c("K10", "CLR", "a3r") else c("K10", "CLR", "K10_min", "K10_deg")
scopes <- c(raw = "raw", hub_att = "hub_att", hub_mi = "hub_mi")
lev <- do.call(rbind, lapply(mp$name, function(nm) { r <- R[[mp$canonical[mp$name == nm]]]
  do.call(rbind, lapply(names(scopes), function(sc) data.frame(readout = nm, canonical = r$name, scope = sc, truth = TR,
    aupr = vapply(TR, function(e) r[[sc]][[e]]$aupr, 0), base_rate = br[TR], stringsAsFactors = FALSE)))}))
lev$enrichment <- lev$aupr / lev$base_rate
utils::write.csv(lev, file.path(out, "pooled_levels.csv"), row.names = FALSE)
con <- do.call(rbind, lapply(mp$name, function(nm) { r <- R[[mp$canonical[mp$name == nm]]]
  do.call(rbind, lapply(refs, function(rf) { q <- R[[mp$canonical[mp$name == rf]]]
    do.call(rbind, lapply(names(scopes), function(sc) do.call(rbind, lapply(TR, function(e) {
      d <- r[[sc]][[e]]$aupr - q[[sc]][[e]]$aupr; dj <- r[[sc]][[e]]$jack - q[[sc]][[e]]$jack; s <- stats::sd(dj)
      data.frame(readout = nm, versus = rf, scope = sc, truth = e, d_aupr = d, lo_aupr = d - 1.96 * s, hi_aupr = d + 1.96 * s,
                 d_enr = d / br[[e]], lo_enr = (d - 1.96 * s) / br[[e]], hi_enr = (d + 1.96 * s) / br[[e]], stringsAsFactors = FALSE) }))))})) }))
utils::write.csv(con, file.path(out, "contrasts.csv"), row.names = FALSE)
mac <- do.call(rbind, lapply(mp$name, function(nm) { x <- R[[mp$canonical[mp$name == nm]]]$macro; x$readout <- nm
  x$bin <- cut(x$size, c(0, 9, 24, 49, Inf), labels = c("3-9", "10-24", "25-49", ">=50")); x }))
ms <- do.call(rbind, lapply(split(mac, mac$readout), function(d) { rbind(data.frame(readout = d$readout[1], bin = "all", n = nrow(d), mean_ap = mean(d$ap, na.rm = TRUE)),
  do.call(rbind, lapply(split(d, d$bin), function(z) data.frame(readout = z$readout[1], bin = z$bin[1], n = nrow(z), mean_ap = mean(z$ap, na.rm = TRUE))))) }))
ms <- ms[order(match(ms$readout, mp$name), match(ms$bin, c("all", "3-9", "10-24", "25-49", ">=50"))), ]; rownames(ms) <- NULL
utils::write.csv(ms, file.path(out, "macro_summary.csv"), row.names = FALSE)
get <- function(nm, sc, e) lev$aupr[lev$readout == nm & lev$scope == sc & lev$truth == e]
## ---- literal pass/fail --------------------------------------------------------------------------------------------------
pf <- NULL
if (idea == 4) {
  bn <- grep("^B[0-9]$", mp$name, value = TRUE)
  pf <- do.call(rbind, lapply(bn, function(b) {
    cb <- get(b, "raw", "CORUM_all") / br[["CORUM_all"]]; ca <- get("a3r", "raw", "CORUM_all") / br[["CORUM_all"]]
    db <- get(b, "raw", "DoRothEA_AB"); dk <- get("K10", "raw", "DoRothEA_AB"); da <- get("a3r", "raw", "DoRothEA_AB")
    pass <- cb >= ca && db >= dk; fail <- cb <= ca && db <= da
    data.frame(readout = b, corum_all_enr = cb, a3r_corum_all_enr = ca, dorothea_aupr = db, k10_dorothea_aupr = dk, a3r_dorothea_aupr = da,
               corum_ge_a3r = cb >= ca, dorothea_ge_k10 = db >= dk, corum_gt_a3r = cb > ca, dorothea_gt_a3r = db > da,
               verdict = if (pass) "PASS" else if (fail) "FAIL" else "between") })) 
} else {
  inn <- grep("^IN-", mp$name, value = TRUE)
  pf <- do.call(rbind, lapply(inn, function(b) {
    sm <- get(b, "raw", "CORUM_small"); sk <- get("K10", "raw", "CORUM_small"); smin <- get("K10_min", "raw", "CORUM_small"); sdeg <- get("K10_deg", "raw", "CORUM_small")
    dd <- get(b, "raw", "DoRothEA_AB") - get("K10", "raw", "DoRothEA_AB")
    data.frame(readout = b, small_aupr = sm, k10_small = sk, min_small = smin, deg_small = sdeg, small_over_k10 = sm / sk,
               beats_min = sm > smin, beats_deg = sm > sdeg, d_dorothea_vs_k10 = dd, guard_breached = dd < -0.005,
               verdict = if (!(sm > smin)) "DROP (does not beat reciprocity min on small)" else if (sm > sdeg) "not dropped (beats min and degree-scaled)" else "not dropped (beats min only)") }))
}
utils::write.csv(pf, file.path(out, "passfail.csv"), row.names = FALSE)
## ---- markdown tables ------------------------------------------------------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "f"))
f4 <- function(x) ifelse(is.na(x), "", formatC(x, digits = 4, format = "f"))
mdt <- function(d) { h <- paste0("| ", paste(names(d), collapse = " | "), " |"); s <- paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|")
  paste(c(h, s, apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))), collapse = "\n") }
md <- c(sprintf("# Tables, idea %d (generated by analysis/explore_ideas47/ideas47_summarize.R; EXPLORATORY, in-sample)\n", idea))
nm_order <- mp$name
for (sc in names(scopes)) {
  d <- do.call(rbind, lapply(nm_order, function(nm) { x <- lev[lev$readout == nm & lev$scope == sc, ]; v <- stats::setNames(f3(x$enrichment), x$truth)[TR]; c(readout = nm, v) }))
  md <- c(md, sprintf("## Pooled enrichment (AUPR / base rate), scope %s\n", sc), mdt(as.data.frame(d, check.names = FALSE)), "", sprintf("Base rates: %s\n", paste(sprintf("%s %.4f", TR, br[TR]), collapse = ", ")))
}
for (sc in names(scopes)) for (rf in refs[1:min(2, length(refs))]) {
  d <- do.call(rbind, lapply(setdiff(nm_order, rf), function(nm) { x <- con[con$readout == nm & con$versus == rf & con$scope == sc, ]
    v <- stats::setNames(sprintf("%+.4f [%+.4f, %+.4f]", x$d_aupr, x$lo_aupr, x$hi_aupr), x$truth)[TR]; c(readout = nm, v) }))
  md <- c(md, sprintf("## Contrast in AUPR units, readout minus %s, scope %s (paired delete-half jackknife, d +/- 1.96 SD)\n", rf, sc), mdt(as.data.frame(d, check.names = FALSE)), "")
}
for (rf in setdiff(refs, refs[1:2])) {
  d <- do.call(rbind, lapply(setdiff(nm_order, rf), function(nm) { x <- con[con$readout == nm & con$versus == rf & con$scope == "raw", ]
    v <- stats::setNames(sprintf("%+.4f [%+.4f, %+.4f]", x$d_aupr, x$lo_aupr, x$hi_aupr), x$truth)[TR]; c(readout = nm, v) }))
  md <- c(md, sprintf("## Contrast in AUPR units, readout minus %s, scope raw\n", rf), mdt(as.data.frame(d, check.names = FALSE)), "")
}
md <- c(md, "## Per-complex (macro) AP, raw CORUM complexes with >= 10 pairs, mean by size bin (point estimates)\n",
        mdt(data.frame(readout = ms$readout, bin = ms$bin, n = ms$n, mean_AP = f4(ms$mean_ap))), "")
pfm <- pf; pfm[] <- lapply(pfm, function(x) if (is.numeric(x)) formatC(x, digits = 4, format = "f") else as.character(x))
md <- c(md, "## Literal pass/fail (plan sections 1-2; point estimates)\n", mdt(pfm), "")
writeLines(md, file.path(out, "tables.md")); cat(paste(md, collapse = "\n"), "\n")
