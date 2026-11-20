#!/usr/bin/env Rscript
# Merges results/tcaclr2/edge_<dataset>.tsv and results/tcaclr2/score_<dataset>/results.csv into results/tcaclr2/summary.csv (one row per dataset and truth scope).
#   Rscript analysis/tcaclr2/collect_tcaclr2.R [--out results/tcaclr2]
a <- commandArgs(TRUE); out <- if (length(a) >= 2L && a[1L] == "--out") a[2L] else "results/tcaclr2"
rs <- Sys.glob(file.path(out, "score_*", "results.csv")); if (!length(rs)) stop("no score_*/results.csv under ", out)
r <- do.call(rbind, lapply(rs, utils::read.csv, stringsAsFactors = FALSE))
e <- do.call(rbind, lapply(Sys.glob(file.path(out, "edge_*.tsv")), utils::read.delim, stringsAsFactors = FALSE))
m <- merge(r, e[, c("dataset", "edge_min", "edge_max", "perms", "eps")], by.x = "truth", by.y = "dataset", all.x = TRUE, sort = FALSE)
f <- file.path(out, "summary.csv"); utils::write.csv(m, f, row.names = FALSE); print(m[, c("truth", "scope", "t_tcaCLR2", "aupr_tcaCLR", "aupr_tcaCLR2", "ratio", "ratio_gene_lo", "ratio_gene_hi")], row.names = FALSE, digits = 4); cat("wrote", f, "\n")
