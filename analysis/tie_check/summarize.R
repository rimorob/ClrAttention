#!/usr/bin/env Rscript
# Merges the parts written by tie_check.R into csv/markdown summaries
# (results/tie_check/). Table cells and bold/dagger marks follow
# analysis/sota_table.R exactly (pooled AUPR over base rate; bold = best in row
# or tied with it by the paired delete-half jackknife difference, Bonferroni
# over the row's comparisons; dagger = tied only after the correction), applied
# once per tie convention. Validates the old-convention numbers against each
# run's metrics.csv and jack_replicates.csv.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "results/tie_check"
runs <- c(M3D = "baselines_m3d", EcoMAC = "baselines_ecomac", K562 = "baselines_k562ess_truthgenes_genie3",
          `K562-8559g` = "baselines_k562ess_genie3")
methods <- c(attention = "Attention", clr = "CLR", genie3 = "GENIE3", pidc = "PIDC", grnboost2 = "GRNBoost2",
             abs_pearson = "|r|", hub = "Hub only")
hub_methods <- c("attention", "clr", "genie3", "pidc", "grnboost2", "abs_pearson")
convs <- c("old", "expected", "optimistic", "negfirst")

P <- list(); ties <- list(); mat <- list(); chk <- list()
for (d in names(runs)) {
  fs <- list.files(file.path(out, "parts", runs[[d]]), full.names = TRUE)
  if (!length(fs)) next
  met <- utils::read.csv(file.path("results", runs[[d]], "metrics.csv"), stringsAsFactors = FALSE)
  jrep <- utils::read.csv(file.path("results", runs[[d]], "jack_replicates.csv"), stringsAsFactors = FALSE)
  for (f in fs) { x <- readRDS(f)
    for (kind in c("plain", "hubstrat")) { r <- x[[kind]]; if (is.null(r)) next
      mname <- if (kind == "plain") x$method else paste0(x$method, "_hubstrat")
      if (kind == "plain") mat[[length(mat) + 1L]] <- data.frame(dataset = d, run = runs[[d]], method = x$method, t(r$.matrix))
      for (tr in setdiff(names(r), ".matrix")) { z <- r[[tr]]
        for (cv in convs) P[[length(P) + 1L]] <- list(dataset = d, run = runs[[d]], method = mname, truth = tr, conv = cv,
                                                       aupr = z$pooled[[cv]], base_rate = z$base_rate, jack = z$jack[, cv])
        ties[[length(ties) + 1L]] <- data.frame(dataset = d, run = runs[[d]], method = mname, truth = tr, t(z$ties), base_rate = z$base_rate)
        ref <- met$aupr[met$method == mname & met$truth == tr]; jr <- jrep$aupr[jrep$method == mname & jrep$truth == tr]
        chk[[length(chk) + 1L]] <- data.frame(run = runs[[d]], method = mname, truth = tr,
          max_abs_diff_pooled = if (length(ref)) abs(z$pooled[["old"]] - ref) else NA,
          max_abs_diff_jack = if (length(jr) == nrow(z$jack)) max(abs(z$jack[, "old"] - jr[order(jrep$rep[jrep$method == mname & jrep$truth == tr])])) else NA) } } }
}
pooled <- do.call(rbind, lapply(P, function(p) data.frame(p[c("dataset", "run", "method", "truth", "conv", "aupr", "base_rate")], ratio = p$aupr / p$base_rate)))
ties <- do.call(rbind, ties); mat <- do.call(rbind, mat); chk <- do.call(rbind, chk)
ties$frac_tied_pairs <- ties$n_tied_pairs / ties$n_pairs; ties$frac_tied_pos <- ties$n_tied_pos / ties$n_pos
dir.create(out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(pooled, file.path(out, "pooled_aupr_by_convention.csv"), row.names = FALSE)
utils::write.csv(ties, file.path(out, "tie_fractions_benchmark_pairs.csv"), row.names = FALSE)
utils::write.csv(mat, file.path(out, "tie_fractions_all_gene_pairs.csv"), row.names = FALSE)
utils::write.csv(chk, file.path(out, "check_vs_existing_metrics.csv"), row.names = FALSE)
jk <- do.call(rbind, lapply(P, function(p) data.frame(dataset = p$dataset, run = p$run, method = p$method, truth = p$truth, conv = p$conv, rep = seq_along(p$jack), aupr = p$jack)))
utils::write.csv(jk, file.path(out, "jack_replicates_by_convention.csv"), row.names = FALSE)
cat("max |old - metrics.csv| pooled:", max(chk$max_abs_diff_pooled, na.rm = TRUE), " jackknife:", max(chk$max_abs_diff_jack, na.rm = TRUE), "\n")

## ---- Table 2 cells with marks, per convention ------------------------------------------------
getp <- function(d, m, tr, cv) { i <- which(pooled$dataset == d & pooled$method == m & pooled$truth == tr & pooled$conv == cv); if (length(i)) pooled$aupr[i] else NA_real_ }
getj <- function(d, m, tr, cv) { p <- P[[which(vapply(P, function(q) q$dataset == d && q$method == m && q$truth == tr && q$conv == cv, NA))[1]]]; p$jack }
cells <- list()
for (hubstrat in c(FALSE, TRUE)) {
  ms <- if (hubstrat) hub_methods else names(methods)
  for (d in unique(pooled$dataset)) for (tr in unique(pooled$truth[pooled$dataset == d])) for (cv in convs) {
    key <- if (hubstrat) paste0(ms, "_hubstrat") else ms
    a <- stats::setNames(vapply(key, function(k) getp(d, k, tr, cv), 0), ms)
    b <- pooled$base_rate[pooled$dataset == d & pooled$truth == tr][1]
    best <- names(which.max(a)); ncomp <- sum(!is.na(a)) - 1L; zc <- stats::qnorm(1 - 0.025 / max(1L, ncomp))
    z <- vapply(ms, function(m) { if (is.na(a[[m]]) || m == best) return(NA_real_)
      (a[[best]] - a[[m]]) / stats::sd(getj(d, key[ms == best], tr, cv) - getj(d, key[ms == m], tr, cv)) }, 0)
    mark <- vapply(ms, function(m) if (is.na(a[[m]])) "-" else if (m == best || (!is.na(z[[m]]) && z[[m]] < 1.96)) "bold"
                   else if (!is.na(z[[m]]) && z[[m]] < zc) "dagger" else "", "")
    cells[[length(cells) + 1L]] <- data.frame(table = if (hubstrat) "hub-controlled" else "raw", dataset = d, truth = tr, conv = cv,
      method = ms, ratio = a / b, mark = mark, rank = rank(-a, na.last = "keep", ties.method = "min"), best = (ms == best), z_vs_best = z)
  }
}
cells <- do.call(rbind, cells)
utils::write.csv(cells, file.path(out, "table2_cells_by_convention.csv"), row.names = FALSE)
fmt <- function(r, mk) { v <- sprintf("%.2f", r); if (is.na(r)) "-" else if (mk == "bold") paste0("**", v, "**") else if (mk == "dagger") paste0("**", v, "**†") else v }

L <- c("# Tie convention check for pooled AUPR (Table 2)", "",
  "Generated by `analysis/tie_check/summarize.R` from `results/tie_check/parts/`. Pooled co-membership AUPR over base rate; bold/dagger marks recomputed by the rule of `analysis/sota_table.R` (paired delete-half jackknife, same 100 splits) separately under each convention. `old` = `ap_fast` (every pair in a tied block gets the precision at the end of the block). `expected` = exact expectation under uniformly random order within each tied block. `optimistic` = positives first inside tied blocks. `negfirst` = negatives first (strict lower bound).", "")
for (tb in c("raw", "hub-controlled")) for (d in unique(cells$dataset)) {
  cc <- cells[cells$table == tb & cells$dataset == d, ]; ms <- unique(cc$method)
  L <- c(L, sprintf("## %s, %s", tb, d), "", paste0("| Truth | Convention | ", paste(methods[ms], collapse = " | "), " |"),
         paste0("|---|---|", paste(rep("---", length(ms)), collapse = "|"), "|"))
  for (tr in unique(cc$truth)) for (cv in convs) { x <- cc[cc$truth == tr & cc$conv == cv, ]
    L <- c(L, sprintf("| %s | %s | %s |", gsub("_", " ", tr), cv, paste(vapply(ms, function(m) { y <- x[x$method == m, ]; fmt(y$ratio, y$mark) }, ""), collapse = " | "))) }
  L <- c(L, "")
}
# changes in rank / marks relative to old
L <- c(L, "## Changes relative to the existing convention", "")
chg <- list()
for (cv in c("expected", "optimistic", "negfirst")) {
  o <- cells[cells$conv == "old", ]; n <- cells[cells$conv == cv, ]
  mm <- merge(o, n, by = c("table", "dataset", "truth", "method"), suffixes = c(".old", ".new"))
  mm <- mm[!is.na(mm$ratio.old), ]
  ch <- mm[mm$rank.old != mm$rank.new | mm$mark.old != mm$mark.new, ]
  if (nrow(ch)) chg[[cv]] <- data.frame(conv = cv, ch[, c("table", "dataset", "truth", "method", "ratio.old", "ratio.new", "rank.old", "rank.new", "mark.old", "mark.new")])
}
chg <- if (length(chg)) do.call(rbind, chg) else NULL
if (!is.null(chg)) { utils::write.csv(chg, file.path(out, "table2_changes_vs_old.csv"), row.names = FALSE)
  L <- c(L, "| Convention | Table | Dataset | Truth | Method | Ratio old | Ratio new | Rank old->new | Mark old->new |", "|---|---|---|---|---|---|---|---|---|")
  for (k in seq_len(nrow(chg))) L <- c(L, sprintf("| %s | %s | %s | %s | %s | %.3f | %.3f | %d->%d | %s->%s |", chg$conv[k], chg$table[k], chg$dataset[k], chg$truth[k], chg$method[k],
    chg$ratio.old[k], chg$ratio.new[k], chg$rank.old[k], chg$rank.new[k], ifelse(chg$mark.old[k] == "", "none", chg$mark.old[k]), ifelse(chg$mark.new[k] == "", "none", chg$mark.new[k]))) }
L <- c(L, "")

## ---- tie fractions ---------------------------------------------------------------------------
L <- c(L, "## Tied fractions (benchmark pairs: pairs between annotated genes of that truth)", "",
       "| Dataset | Truth | Method | Pairs | Positives | Tied pairs (%) | Tied positives | Tied positives (%) | Largest block (%) |", "|---|---|---|---|---|---|---|---|---|")
tp <- ties[!grepl("_hubstrat$", ties$method), ]
for (k in seq_len(nrow(tp))) L <- c(L, sprintf("| %s | %s | %s | %d | %d | %.1f | %d | %.1f | %.1f |", tp$dataset[k], tp$truth[k], tp$method[k], tp$n_pairs[k], tp$n_pos[k],
  100 * tp$frac_tied_pairs[k], tp$n_tied_pos[k], 100 * tp$frac_tied_pos[k], 100 * tp$max_block[k] / tp$n_pairs[k]))
L <- c(L, "", "Tied fraction over all off-diagonal gene pairs (full matrix; pairs sharing a score with at least one other pair):", "",
       "| Dataset | Method | Gene pairs | Tied (%) | Largest block (%) | Largest block value |", "|---|---|---|---|---|---|")
for (k in seq_len(nrow(mat))) L <- c(L, sprintf("| %s | %s | %d | %.1f | %.1f | %g |", mat$dataset[k], mat$method[k], mat$n_pairs[k], 100 * mat$frac_tied[k], 100 * mat$frac_largest_block[k], mat$largest_block_value[k]))
L <- c(L, "")

## ---- attention vs CLR on the E. coli truths -----------------------------------------------------------
L <- c(L, "## Attention minus CLR, pooled AUPR, paired jackknife (estimate +/- 1.96 SD of the 100 paired replicate differences)", "",
       "| Dataset | Truth | Convention | Attention | CLR | Difference | 95% interval | Ratio att/CLR |", "|---|---|---|---|---|---|---|---|")
for (d in unique(pooled$dataset)) for (tr in unique(pooled$truth[pooled$dataset == d])) for (cv in convs) {
  a <- getp(d, "attention", tr, cv); c0 <- getp(d, "clr", tr, cv); dj <- getj(d, "attention", tr, cv) - getj(d, "clr", tr, cv); s <- stats::sd(dj)
  L <- c(L, sprintf("| %s | %s | %s | %.4f | %.4f | %+.4f | [%+.4f, %+.4f] | %.3f |", d, tr, cv, a, c0, a - c0, a - c0 - 1.96 * s, a - c0 + 1.96 * s, a / c0)) }
writeLines(L, file.path(out, "summary.md"))
cat(L, sep = "\n")
