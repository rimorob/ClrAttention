#!/usr/bin/env Rscript
# Validation of the reconstructed pseudobulk (build_pseudobulk_from_sc.R) on K562, where the authors' normalized-bulk file exists:
# (a) agreement of the mean-z profiles; (b) agreement of the energy-test call p < 1e-3 and of the selected set K of the recipe.
# Label-free: only expression profiles and the shipped energy-test p-values; no truth, no score.
# Usage: Rscript validate_pseudobulk_k562.R <rebuilt h5ad> <shipped h5ad> <outdir>
suppressPackageStartupMessages(library(rhdf5))
a <- commandArgs(trailingOnly = TRUE); fn <- a[1]; fs <- a[2]; OUT <- a[3]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat_of <- function(f, grp, field) { x <- h5read(f, paste0(grp, "/", field)); if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(h5read(f, paste0(grp, "/__categories/", field)), error = function(e) NULL); if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
Xn <- h5read(fn, "X"); pn <- cat_of(fn, "/obs", "gene_transcript"); gn <- cat_of(fn, "/var", "gene_name"); vn <- as.numeric(h5read(fn, "/obs/energy_test_p_value"))
Xs <- h5read(fs, "X"); ps <- cat_of(fs, "/obs", "gene_transcript"); gs <- cat_of(fs, "/var", "gene_name"); vs <- as.numeric(h5read(fs, "/obs/energy_test_p_value"))
nc_all <- as.integer(h5read(fs, "/obs/num_cells_filtered"))
stopifnot(identical(gn, gs)); cat("shipped labels:", length(ps), "; rebuilt labels:", length(pn), "; shipped labels absent from the rebuild (no cells in the deposited file):", sum(!ps %in% pn), "(all non-targeting:", all(grepl("non-targeting", ps[!ps %in% pn])), ")\n")
sel <- ps %in% pn; Xs <- Xs[, sel]; vs <- vs[sel]; nc_all <- nc_all[sel]; ps <- ps[sel]; m <- match(ps, pn); cat("order identical:", identical(ps, pn), "\n")
Xn <- Xn[, m]; vn <- vn[m]; pn <- pn[m]
pc <- vapply(seq_along(ps), function(k) stats::cor(Xn[, k], Xs[, k]), 0); sl <- vapply(seq_along(ps), function(k) unname(stats::coef(stats::lm(Xs[, k] ~ Xn[, k]))[2]), 0)
gc_ <- vapply(seq_len(nrow(Xn)), function(g) stats::cor(Xn[g, ], Xs[g, ]), 0)
prof <- data.frame(n_perts = length(ps), per_pert_cor_median = stats::median(pc), per_pert_cor_min = min(pc), per_pert_cor_q01 = stats::quantile(pc, 0.01),
  slope_median = stats::median(sl), slope_q01 = stats::quantile(sl, 0.01), slope_q99 = stats::quantile(sl, 0.99), per_gene_cor_median = stats::median(gc_, na.rm = TRUE),
  max_abs_diff = max(abs(Xn - Xs)), global_cor = stats::cor(as.vector(Xn), as.vector(Xs)))
nt <- grepl("non-targeting", ps); ok <- !nt & is.finite(vs) & is.finite(vn)
a1 <- vs < 1e-3; b1 <- vn < 1e-3
tab <- table(shipped = a1[ok], rebuilt = b1[ok])
Ks <- ps[!nt & a1 & is.finite(vs)]; Kn <- ps[!nt & b1 & is.finite(vn)]
en <- data.frame(tested = sum(ok), shipped_sig = sum(a1[ok]), rebuilt_sig = sum(b1[ok]), both = sum(a1[ok] & b1[ok]), only_shipped = sum(a1[ok] & !b1[ok]), only_rebuilt = sum(!a1[ok] & b1[ok]),
  agreement = mean(a1[ok] == b1[ok]), K_jaccard = length(intersect(Ks, Kn)) / length(union(Ks, Kn)), spearman_logp = stats::cor(log10(vs[ok]), log10(vn[ok]), method = "spearman"))
utils::write.csv(prof, file.path(OUT, "k562_profile_agreement.csv"), row.names = FALSE); utils::write.csv(en, file.path(OUT, "k562_energy_test_agreement.csv"), row.names = FALSE)
nc <- nc_all; bins <- cut(nc[ok], c(0, 20, 40, 80, 160, Inf))
bin <- do.call(rbind, lapply(levels(bins), function(b) { s <- which(bins == b); data.frame(cells = b, n = length(s), shipped_sig = sum(a1[ok][s]), rebuilt_sig = sum(b1[ok][s]), agree = mean(a1[ok][s] == b1[ok][s])) }))
utils::write.csv(bin, file.path(OUT, "k562_energy_test_agreement_by_cells.csv"), row.names = FALSE)
print(t(prof)); print(tab); print(t(en)); print(bin)
