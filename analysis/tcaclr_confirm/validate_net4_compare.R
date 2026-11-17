#!/usr/bin/env Rscript
# Compares the output of score_confirm.R on the consumed DREAM5 net4 study with the stored study (results/complex_confirm/).
# Usage (repo root): Rscript analysis/tcaclr_confirm/validate_net4_compare.R <score_confirm out dir> [<exp source name of the refit run to compare, ...>]
# Writes <out>/validation_compare.csv and prints a summary; exits 1 if any registered comparison exceeds its tolerance.
args <- commandArgs(TRUE); out <- args[1]; refit_srcs <- args[-1]; stored <- "results/complex_confirm"
rd <- function(f) utils::read.csv(f, stringsAsFactors = FALSE)
rows <- list(); add <- function(item, quantity, old, new, tol, gating = TRUE) rows[[length(rows) + 1L]] <<- data.frame(item = item, quantity = quantity, stored = as.numeric(old), new = as.numeric(new),
  abs_diff = abs(as.numeric(old) - as.numeric(new)), tol = tol, gating = gating, stringsAsFactors = FALSE)
nm <- c(P = "P", C1 = "C1", C2 = "C2")
## 1. all contrasts on every truth and scope (stored truth stage)
sc <- rd(file.path(stored, "net4/truth/contrasts.csv")); sc <- sc[sc$role != "complex_level_jackknife_sensitivity" & sc$role != "descriptive_avg_input" & sc$truth != "dream5_ym", ]
mine <- rd(file.path(out, "contrasts.csv")); mine$scope <- as.character(mine$scope)
for (k in seq_len(nrow(sc))) { x <- sc[k, ]; y <- mine[mine$truth == x$truth & mine$scope == x$scope & mine$row == x$row, ]; if (nrow(y) != 1L) { add(paste(x$truth, x$scope, x$row), "MISSING in new output", 0, 1, 0); next }
  it <- paste(x$truth, x$scope, x$row)
  add(it, "d", x$d, y$d, 1e-12); add(it, "se_gene", x$se_gene, y$se_gene, 1e-12); add(it, "gene_lo", x$gene_lo, y$gene_lo, 1e-12); add(it, "gene_hi", x$gene_hi, y$gene_hi, 1e-12)
  add(it, "pct_lo", x$pct_lo, y$gene_pct_lo, 1e-12); add(it, "pct_hi", x$pct_hi, y$gene_pct_hi, 1e-12); add(it, "skew", x$skew, y$gene_skew, 1e-10) }
## 2. complex-level delete-half sensitivity (main)
cl <- rd(file.path(stored, "net4/truth/contrasts.csv")); cl <- cl[cl$role == "complex_level_jackknife_sensitivity", ]; my <- rd(file.path(out, "complex_level_sensitivity.csv"))
pairs_ <- list(P = c("tcaCLR", "aCLR"), C1 = c("tcaCLR", "CLR"), C2 = c("aCLR", "CLR"))
for (r in names(pairs_)) { x <- cl[cl$row == r, ]; y <- my[my$method_a == pairs_[[r]][1] & my$method_b == pairs_[[r]][2], ]; it <- paste("complex-level", r)
  add(it, "d", x$d, y$d, 1e-12); add(it, "SE_C", x$se_gene, y$se_complex, 1e-12); add(it, "lo", x$gene_lo, y$cl_lo, 1e-12); add(it, "hi", x$gene_hi, y$cl_hi, 1e-12); add(it, "pct_lo", x$pct_lo, y$cl_pct_lo, 1e-12); add(it, "pct_hi", x$pct_hi, y$cl_pct_hi, 1e-12); add(it, "skew", x$skew, y$cl_skew, 1e-10) }
## 3. the stored estimate (experiment components, total intervals) from the stored net4 replicates (series and group)
est <- rd(file.path(stored, "estimate.csv")); alt <- rd(file.path(out, "estimate_by_experiment_source.csv"))
for (r in names(pairs_)) { e <- est[est$row == r, ]
  for (src in c("series", "group")) { y <- alt[alt$exp_source == src & alt$row == r, ]; if (!nrow(y)) next; it <- paste("estimate", r, src, "(stored replicates)")
    add(it, "d", e$d, y$d, 1e-12); add(it, "se_gene", e$se_gene, y$se_gene, 1e-12); add(it, "se_exp", e[[paste0("se_", src)]], y$se_exp, 1e-12)
    add(it, "total_lo", e[[paste0("total_", src, "_lo")]], y$total_lo, 1e-12); add(it, "total_hi", e[[paste0("total_", src, "_hi")]], y$total_hi, 1e-12) } }
## 4. experiment SD of every truth (estimate_final.csv, descriptive_sd_series_refits) is not recomputed by the driver; the primary row is above
## 5. per-complex descriptors against per_complex_summary.csv (a3r minus attention, complexes with >= 10 positives)
pcs <- rd(file.path(stored, "net4/truth/per_complex_summary.csv")); br <- rd(file.path(out, "breadth.csv")); x <- pcs[pcs$comparison == "a3r minus attention", ]
add("per-complex AP a3r - attention", "n_complexes_ge10pos", x$n_complexes_ge10pos, br$n_complexes, 0); add("per-complex AP a3r - attention", "n_improved", x$n_improved, br$n_improved, 0)
add("per-complex AP a3r - attention", "n_worse", x$n_worse, br$n_worse, 0); add("per-complex AP a3r - attention", "mean_d_ap", x$mean_d_ap, br$mean_d_ap, 1e-12); add("per-complex AP a3r - attention", "median_d_ap", x$median_d_ap, br$median_d_ap, 1e-12)
## 6. truth structure: md5 of (i, j, label), pair and positive counts of every truth
tcs <- rd(file.path(stored, "net4/truth/truth_counts.csv")); tcm <- rd(file.path(out, "score/truth_counts.csv"))
for (v in setdiff(tcs$truth, "dream5_ym")) { a <- tcs[tcs$truth == v, ]; b <- tcm[tcm$truth == v, ]
  add(paste("truth", v), "structure md5 identical (1 = yes)", 1, as.numeric(identical(a$structure_md5, b$structure_md5)), 0); add(paste("truth", v), "n_pairs", a$n_pairs, b$n_pairs, 0); add(paste("truth", v), "n_pos", a$n_pos, b$n_pos, 0) }
## 7. n_eff (the design states 15.07) and the descriptors of the stored table (printed, no stored reference)
ne <- rd(file.path(out, "score/neff.csv")); cn <- rd("docs/complex_confirm_selections/yeast_complexes_counts.csv"); cn <- cn[cn$variant == "main", ]
add("n_eff main", "overlap-counted (design 15.07; the gate's definition)", cn$n_eff_overlap_counted, ne$n_eff, 1e-10); add("n_eff main", "fractional assignment (design 14.5)", cn$n_eff, ne$n_eff_fractional, 1e-10)
## 8. fixed sentences: with the net4 names the generator reproduces the stored sentences
e <- est[est$row == "P", ]; h <- rd(file.path(out, "estimate.csv"))
for (s in c("S_zero", "S_min", "S_ecoli", "S_fixed")) add(paste("sentence", s), "identical to estimate.csv (1 = yes)", 1, as.numeric(identical(h[[s]], e[[s]])), 0)
## 9. fresh refits (optional): replicate AUPRs and SE_exp against the stored replicates, when the run used the net4 halves
for (src in refit_srcs) { f <- file.path(out, "expjack", paste0("replicates_", src, ".csv")); stored_f <- file.path(stored, "net4/expjack", paste0("replicates_", sub("_refit$", "", src), ".csv"))
  if (!file.exists(f) || !file.exists(stored_f)) next; a <- rd(stored_f); a$method <- c(clr = "CLR", attention = "aCLR", a3r = "tcaCLR")[a$readout]; b <- rd(f)
  m <- merge(a[, c("rep", "truth", "method", "aupr")], b[, c("rep", "truth", "method", "aupr")], by = c("rep", "truth", "method"), suffixes = c("_stored", "_new")); add(paste("fresh refits", src), "replicate AUPRs compared", nrow(m), nrow(m), 0)
  add(paste("fresh refits", src), "max abs difference of replicate AUPR (all truths, CLR, aCLR, tcaCLR)", 0, max(abs(m$aupr_stored - m$aupr_new)), 1e-12)
  yy <- rd(file.path(out, "estimate_by_experiment_source.csv")); yy <- yy[yy$exp_source == src & yy$row == "P", ]; ee <- est[est$row == "P", ]; sname <- sub("_refit$", "", src)
  add(paste("fresh refits", src, "P"), "se_exp", ee[[paste0("se_", sname)]], yy$se_exp, 1e-12); add(paste("fresh refits", src, "P"), "total_lo", ee[[paste0("total_", sname, "_lo")]], yy$total_lo, 1e-12); add(paste("fresh refits", src, "P"), "total_hi", ee[[paste0("total_", sname, "_hi")]], yy$total_hi, 1e-12)
  hf <- file.path(out, "expjack", paste0("refit_hashes_", src, ".csv")); mf <- if (sname == "series") "results/dream5/net4/refit_series/REFIT_SERIES_MANIFEST.tsv" else "results/dream5/net4/refit/REFIT_MANIFEST.tsv"
  if (file.exists(hf) && file.exists(mf)) { hh <- rd(hf); mm <- utils::read.delim(mf, stringsAsFactors = FALSE, colClasses = c(sha256 = "character")); mm$method <- c(clr = "CLR", attention = "aCLR")[mm$method]; mm <- mm[mm$net == 4 & !is.na(mm$method), ]
    j <- merge(hh, mm, by = c("rep", "method"), suffixes = c("_new", "_manifest")); add(paste("fresh refits", src), "packed-score SHA-256 equal to the DREAM5 manifest, CLR and aCLR (count of mismatches)", 0, sum(j$sha256_new != j$sha256_manifest), 0); add(paste("fresh refits", src), "hashes compared", nrow(j), nrow(j), 0) } }
res <- do.call(rbind, rows); res$pass <- res$abs_diff <= res$tol; utils::write.csv(res, file.path(out, "validation_compare.csv"), row.names = FALSE)
cat(sprintf("%d comparisons, %d pass, %d fail; largest abs difference among numeric quantities %.3g\n", nrow(res), sum(res$pass), sum(!res$pass), max(res$abs_diff[res$tol > 0], na.rm = TRUE)))
if (any(!res$pass)) { print(res[!res$pass, ], row.names = FALSE); quit(save = "no", status = 1L) }
