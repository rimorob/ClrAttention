#!/usr/bin/env Rscript
# Human perturbation compendium: Perturb-seq pseudobulk profiles as
# "experiments" (Replogle et al. 2022). Co-expression network inference
# assumes the samples explore many network states; unperturbed cells of one
# condition do not (analysis/perturbseq.R is that contrast). Here each sample
# is one knockdown's mean expression change against control cells (z-scored
# per gene; the authors' normalized pseudobulk), the human analogue of the
# M3D experiment compendium.
#
# Design:
#   - knockdowns with a significant transcriptional phenotype (energy test
#     p < --p_max) are split at random into halves A and B;
#   - networks are built on A only;
#   - truths: DoRothEA A-B and A-C, CollecTRI, CORUM (literature), and
#     knockdown-response modules built from B only (no shared samples);
#   - methods: |r|, MI, CLR on each, CLR on |r| of normal scores (CLR-NS), and
#     soft10 attention on every CLR variant at a depth chosen on random halves
#     of A (criterion: top 1% of the same CLR matrix on the other half);
#   - intervals: paired delete-half gene jackknife (AUPR) for attention minus
#     its CLR and for CLR minus |r|.
#
# Usage (repo root):
#   Rscript analysis/perturbseq_bulk.R --h5 data/perturbseq/K562_essential_normalized_bulk_01.h5ad
#       [--out results/perturbseq_bulk_k562ess] [--p_max 1e-3] [--workers 6] [--tmax 32]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(h5 = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad", out = "results/perturbseq_bulk",
            p_max = 1e-3, workers = 6L, mem_gb = 12, tmax = 32L, seed = 20261015L, alpha = 0.5,
            threads = NULL, jack = 100L)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% c("workers", "tmax", "seed", "threads", "jack")) as.integer(v) else
    if (key %in% c("p_max", "mem_gb", "alpha")) as.numeric(v) else v
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "perturbseq_bulk.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
set.seed(opt$seed)
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

## ---- compendium -----------------------------------------------------------------------
rd <- function(p) rhdf5::h5read(opt$h5, p)
cat_of <- function(grp, field) {
  x <- rd(paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L]
}
Z <- rd("X"); genes <- cat_of("/var", "gene_name"); rownames(Z) <- genes
pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
pval <- as.numeric(rd("/obs/energy_test_p_value"))
ok <- pval < opt$p_max & !grepl("non-targeting", pert)
K <- pert[ok]; K <- K[sample(length(K))]
A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
X <- Z[, A]; rm(Z); invisible(gc())
# Missing values (present in the genome-wide K562 and RPE1 files): drop genes
# missing in more than 1% of the profiles, set the rest to 0 (no change from
# control cells).
nf <- !is.finite(X); bad_gene <- rowMeans(nf) > 0.01
X <- X[!bad_gene, , drop = FALSE]; n_imp <- sum(!is.finite(X)); X[!is.finite(X)] <- 0
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; genes <- rownames(X); G <- length(genes)
say(sprintf("non-finite values: %d genes dropped (> 1%% missing), %d values set to 0", sum(bad_gene), n_imp))
say(sprintf("%d profiles; %d with a significant phenotype (p < %g): %d for networks (A), %d for response modules (B); %d genes",
            length(pert), length(K), opt$p_max, length(A), length(B), G))


## ======================================================================================================================
## Appended for docs/plans/tcaclr_confirm_plan.md (Arm A, HepG2): everything ABOVE this line is analysis/perturbseq_bulk.R lines 1-71,
## byte-identical (the K562 recipe: significant-phenotype filter, seeded half split, missing and constant-gene filters); run with
## --h5 <HepG2 bulk h5ad> --out <dir>. Below: the duplicate-symbol rule of analysis/baselines.R (drop every row of a duplicated symbol, as for RPE1
## HSPA14) and the export of the network half A. Nothing below reads a truth; nothing is scored.
dup <- rownames(X) %in% rownames(X)[duplicated(rownames(X))]
say(sprintf("duplicated gene names: %d rows (%s); constant genes already removed above", sum(dup), paste(unique(rownames(X)[dup]), collapse = ", ")))
writeLines(rownames(X)[dup], file.path(opt$out, "genes_dropped.txt")); X <- X[!dup, , drop = FALSE]
expr_f <- file.path(opt$out, "expression.tsv")
utils::write.table(data.frame(gene = rownames(X), X, check.names = FALSE), expr_f, sep = "\t", quote = FALSE, row.names = FALSE)
sym_of <- function(p) sub("^[0-9]+_(.+)_P[0-9P]+_ENSG[0-9]+$", "\\1", p)
writeLines(rownames(X), file.path(opt$out, "universe_genes.txt"))
writeLines(sort(unique(sym_of(A))), file.path(opt$out, "perturbed_genes_networkhalf.txt"))
writeLines(sort(unique(sym_of(c(A, B)))), file.path(opt$out, "perturbed_genes_significant.txt"))
writeLines(sort(unique(sym_of(pert[!grepl("non-targeting", pert)]))), file.path(opt$out, "perturbed_genes_panel.txt"))
utils::write.csv(data.frame(profiles_in_file = length(pert), nontargeting_labels = sum(grepl("non-targeting", pert)), significant = length(K), half_A_network = length(A),
  half_B = length(B), genes_in_file = length(rd("/var/gene_name")), genes_after_filters = G, duplicated_rows_dropped = sum(dup), universe = nrow(X),
  profiles_in_matrix = ncol(X), expression_sha256 = digest::digest(file = expr_f, algo = "sha256"),
  perturbed_in_A_unique_genes = length(unique(sym_of(A)))), file.path(opt$out, "export_counts.csv"), row.names = FALSE)
say(sprintf("wrote %s: %d genes x %d profiles", expr_f, nrow(X), ncol(X)))
