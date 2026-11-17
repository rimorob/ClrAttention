#!/usr/bin/env Rscript
# Builds the input directory of the consumed-data validation (DREAM5 net4 yeast, Complex Portal yeast) in the formats score_confirm.R reads.
# Consumed data only. Usage (repo root): Rscript analysis/tcaclr_confirm/validate_net4_prep.R <input dir>
suppressPackageStartupMessages({ source("analysis/dream5_lib.R"); source("analysis/tcaclr_confirm/tc_lib.R") })
args <- commandArgs(TRUE); dir <- if (length(args)) args[1] else stop("input directory required"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
nd <- d5_read_net(4L); G <- length(nd$gene_names); stopifnot(G == 5950L)
sc <- file.path(dir, "scores"); dir.create(sc, showWarnings = FALSE); writeLines(nd$gene_names, file.path(sc, "genes.txt"))
lf <- "results/complex_confirm/net4/labelfree/cache"; dd <- "results/dream5/net4"
ln <- function(from, to) { to <- file.path(sc, to); if (file.exists(to)) file.remove(to); file.symlink(normalizePath(from), to) }
ln(file.path(lf, "cc_a3r.scores"), "tcaCLR.scores"); ln(file.path(lf, "cc_attention.scores"), "aCLR.scores"); ln(file.path(lf, "cc_clr.scores"), "CLR.scores"); ln(file.path(lf, "cc_mi_rowmean.scores"), "mi_rowmean.scores")
ln(file.path(dd, "abs_pearson__replicate.scores"), "absr.scores"); ln(file.path(dd, "pidc__replicate.scores"), "PIDC.scores")
gm <- c(GENIE3 = "genie3_all", GRNBoost2 = "grnboost2_all")
for (nm in names(gm)) { S <- d5_read_scores(file.path(dd, paste0(gm[[nm]], "__replicate.scores"))); stopifnot(is.matrix(S), all(dim(S) == G))   # directed, symmetrised by max (baselines.R)
  tc_write_bin(tc_pack(pmax(S, t(S))), file.path(sc, paste0(nm, ".scores"))); rm(S); invisible(gc(verbose = FALSE)) }
X <- nd$X; rownames(X) <- nd$gene_names; saveRDS(X, file.path(dir, "expr.rds"))
utils::write.table(data.frame(chip = colnames(X), unit = nd$series), file.path(dir, "units_series.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(chip = colnames(X), unit = nd$group), file.path(dir, "units_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("net4 validation inputs written to", dir, "; genes", G, "; chips", ncol(X), "; series", length(unique(nd$series)), "; groups", nd$n_groups, "\n")
