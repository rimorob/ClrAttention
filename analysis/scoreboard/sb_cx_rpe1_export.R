#!/usr/bin/env Rscript
# Export of the RPE1 expression for the competitor runs (docs/plans/scoreboard_competitors_plan.md, section 8): the matrix load_expr("rpe1")
# of clr_logpt_test.R gives (8,747 genes x 1,023 profiles), with every row of a duplicated gene symbol dropped (HSPA14, 2 rows), as baselines.R
# does for K562. Writes expression.tsv (genes x samples, first column the gene names; the format of baselines.R) and genes_dropped.txt.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
source("analysis/scoreboard/sb_cx_lib.R"); sb_import_truth_defs()
out <- "results/scoreboard/complex_competitors/rpe1_work"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
X <- load_expr("rpe1")$X; dup <- rownames(X) %in% rownames(X)[duplicated(rownames(X))]
sb_say(nrow(X), " genes, dropping ", sum(dup), " rows of duplicated symbols: ", paste(unique(rownames(X)[dup]), collapse = ", "))
writeLines(rownames(X)[dup], file.path(out, "genes_dropped.txt")); X <- X[!dup, ]
utils::write.table(data.frame(gene = rownames(X), X, check.names = FALSE), file.path(out, "expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sb_say("wrote ", nrow(X), " genes x ", ncol(X), " profiles")
