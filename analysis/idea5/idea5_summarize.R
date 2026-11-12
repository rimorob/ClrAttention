#!/usr/bin/env Rscript
# EXPLORATORY summary of analysis/idea5/idea5_readout_score.R caches (docs/plans/explore_idea5_plan.md): point estimates,
# jackknife intervals (difference +/- 1.96 SD of paired replicate differences), the literal pass/fail of plan section 7.
# No p-values, no winner selection. Usage: Rscript analysis/idea5/idea5_summarize.R results/idea5/scores results/idea5/models results/idea5/summary
args <- commandArgs(TRUE); sdir <- args[1]; mdir <- args[2]; odir <- args[3]; dir.create(odir, recursive = TRUE, showWarnings = FALSE)
meta <- readRDS(file.path(sdir, "cache", "meta.rds")); mm <- meta$mm; br <- meta$base_rate; TR <- meta$truths
ld <- function(id) { f <- file.path(sdir, "cache", paste0(id, ".rds")); if (file.exists(f)) readRDS(f) else NULL }
SC <- c(pooled = "raw", hub_att = "hub (pi)", hub_mi = "hub (MI row mean)")
get <- function(id, ro, scope, e) { r <- ld(id); if (is.null(r) || is.null(r[[ro]])) NULL else r[[ro]][[scope]][[e]] }
wcsv <- function(x, n) utils::write.csv(x, file.path(odir, n), row.names = FALSE)

## 1. levels --------------------------------------------------------------------------------------------------------
ids <- meta$tasks; lev <- list()
for (id in ids) { r <- ld(id); if (is.null(r)) next
  for (ro in setdiff(names(r), c("t", "isolated_rows"))) for (sc in names(SC)) for (e in TR) { x <- r[[ro]][[sc]][[e]]
    lev[[length(lev) + 1L]] <- data.frame(model = id, readout = ro, depth = r$t[[ro]], scope = SC[[sc]], truth = e, aupr = x$aupr, enrichment = x$aupr / br[[e]], jack_sd = stats::sd(x$jack)) } }
lev <- do.call(rbind, lev); wcsv(lev, "levels_all_models.csv")
iso <- do.call(rbind, lapply(ids, function(id) { r <- ld(id); if (is.null(r)) NULL else data.frame(model = id, isolated_row_fraction = r$isolated_rows) })); wcsv(iso, "isolated_rows.csv")

## 2. contrasts vs comparators ------------------------------------------------------------------------------------------
contr <- function(a, ra, b, rb, scope, e) { pa <- get(a, ra, scope, e); pb <- get(b, rb, scope, e); if (is.null(pa) || is.null(pb)) return(NULL)
  d <- pa$aupr - pb$aupr; dj <- pa$jack - pb$jack; s <- stats::sd(dj)
  data.frame(model = a, readout = ra, comparator = paste(b, rb), scope = SC[[scope]], truth = e, aupr = pa$aupr, comparator_aupr = pb$aupr, d_aupr = d, lo = d - 1.96 * s, hi = d + 1.96 * s,
             d_enrichment = d / br[[e]], jack_sd = s) }
rows <- list()
for (r in seq_len(nrow(mm))) { k <- mm$kind[r]; id <- mm$tag[r]; rd <- if (k == "half") c("att9", "attT") else if (k == "refit") c("pair", "att1", "att9", "attT") else c("pair", "att1", "att9")
  for (ro in rd) for (sc in names(SC)) for (e in TR) {
    rows[[length(rows) + 1L]] <- cbind(contr(id, ro, "CLR_full", "att9", sc, e), kind = k, init = mm$init[r], d = mm$d[r], seed = mm$seed[r], mode = mm$mode[r], vs = "full-data CLR attention t = 9")
    if (k == "half") rows[[length(rows) + 1L]] <- cbind(contr(id, ro, paste0("CLRhalf_", mm$mode[r]), "att9", sc, e), kind = k, init = mm$init[r], d = mm$d[r], seed = mm$seed[r], mode = mm$mode[r], vs = "matched half-data CLR attention t = 9")
    if (k == "refit" && ro %in% c("att9", "attT")) rows[[length(rows) + 1L]] <- cbind(contr(id, ro, "CLR_full", "pair", sc, e), kind = k, init = mm$init[r], d = mm$d[r], seed = mm$seed[r], mode = mm$mode[r], vs = "full-data CLR pairwise score") } }
ct <- do.call(rbind, rows); wcsv(ct, "contrasts_per_model.csv")

## 3. literal pass/fail (plan section 7) ----------------------------------------------------------------------------------
cls <- function(d, lo) ifelse(d >= 0.005 & lo > 0, "win", ifelse(d < -0.005, "loss", ifelse(abs(d) <= 0.005, "tie", "unresolved gain")))
verdict <- function(x) { n <- table(factor(x, levels = c("win", "loss", "tie", "unresolved gain")))
  v <- if (n[["win"]] >= 3) "learned wins" else if (n[["loss"]] >= 3) "learned loses" else if (n[["tie"]] >= 3) "tie" else "mixed"
  c(verdict = v, counts = paste(names(n), n, sep = " ", collapse = "; ")) }
pf <- list()
sel <- ct[ct$kind == "refit" & ct$vs == "full-data CLR attention t = 9" & ct$readout %in% c("attT", "att9"), ]
sel$class <- cls(sel$d_aupr, sel$lo)
for (key in unique(paste(sel$init, sel$d, sel$readout, sel$scope, sel$truth))) { x <- sel[paste(sel$init, sel$d, sel$readout, sel$scope, sel$truth) == key, ]; v <- verdict(x$class)
  pf[[length(pf) + 1L]] <- data.frame(init = x$init[1], d = x$d[1], readout = x$readout[1], scope = x$scope[1], truth = x$truth[1], n_seeds = nrow(x), mean_d_aupr = mean(x$d_aupr), sd_across_seeds = stats::sd(x$d_aupr),
    min_d = min(x$d_aupr), max_d = max(x$d_aupr), mean_jack_sd = mean(x$jack_sd), verdict = v[["verdict"]], counts = v[["counts"]], primary = x$init[1] == "clr" && x$d[1] == 32 && x$readout[1] == "attT" && x$scope[1] == "hub (pi)" && x$truth[1] == "SC") }
pf <- do.call(rbind, pf); wcsv(pf, "passfail_refit.csv")
sel2 <- ct[ct$kind == "half" & ct$vs == "matched half-data CLR attention t = 9" & ct$readout == "attT", ]; sel2$class <- cls(sel2$d_aupr, sel2$lo)
pf2 <- list(); for (key in unique(paste(sel2$init, sel2$d, sel2$scope, sel2$truth))) { x <- sel2[paste(sel2$init, sel2$d, sel2$scope, sel2$truth) == key, ]; v <- verdict(x$class)
  pf2[[length(pf2) + 1L]] <- data.frame(init = x$init[1], d = x$d[1], scope = x$scope[1], truth = x$truth[1], n_models = nrow(x), mean_d_aupr = mean(x$d_aupr), sd_across_models = stats::sd(x$d_aupr), min_d = min(x$d_aupr), max_d = max(x$d_aupr),
    verdict_rule_applied = paste0(v[["verdict"]], " (thresholds as plan section 7 with 3 of 5 seeds replaced by a majority of the ", nrow(x), " models: informal)"), counts = v[["counts"]]) }
wcsv(do.call(rbind, pf2), "passfail_half_matched.csv")
write.csv(sel[sel$init == "clr" & sel$d == 32 & sel$readout == "attT" & sel$scope == "hub (pi)" & sel$truth == "SC", ], file.path(odir, "primary_per_seed.csv"), row.names = FALSE)

## 4. training summaries -------------------------------------------------------------------------------------------------
js <- do.call(rbind, lapply(list.files(mdir, "\\.json$", full.names = TRUE), function(f) { j <- jsonlite::fromJSON(gsub("NaN", "null", paste(readLines(f), collapse = ""))); data.frame(lapply(j[c("tag", "init", "d", "seed", "mode", "lr", "best_epoch", "epochs_run", "steps", "val_mse_step0", "val_r2_step0", "val_mse_best", "val_r2_best", "seconds", "peak_gpu_mb")], function(v) if (is.null(v)) NA else v)) }))
wcsv(js, "training_runs.csv")
cat("summary written to", odir, "\n")
