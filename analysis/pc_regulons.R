#!/usr/bin/env Rscript
# What do the top principal components of the M3D compendium capture?
# For each of PC1-5 (466 replicate-averaged experiments, genes centred), the
# AUROC with which |gene loading| separates a regulon's members from all other
# genes, with a two-sided Mann-Whitney p-value; sigma-factor regulons, the
# large global TF regulons and the ribosomal-protein and flagellar programs.
# Motivation: removing the top 2-3 PCs improves CLR and attention (Supp.
# Note 7); if those PCs are sigma-factor or growth-rate programs, the removal
# deletes real regulatory biology rather than noise.
# Usage (repo root):  Rscript analysis/pc_regulons.R [--out results/pc_regulons]

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 2 && args[1] == "--out") args[2] else "results/pc_regulons"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("analysis/regulondb.R"); source("analysis/regulons.R")
m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; genes <- rownames(X)
sv <- svd(X - rowMeans(X), nu = 5, nv = 5)
ve <- 100 * sv$d^2 / sum(sv$d^2)
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
reg <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)
sel <- c(grep("^sigma:", names(reg$members), value = TRUE),
         "TF:CRP", "TF:FNR", "TF:ArcA", "TF:Fis", "TF:H-NS", "TF:Lrp", "TF:FlhDC", "TF:RpoS")
sel <- sel[sel %in% names(reg$members)]
# Ribosomal-protein genes (rpl*, rpm*, rps*) as a growth-rate program, via gene symbols.
gt <- utils::read.delim("analysis/ecoli_k12_genes.tsv", stringsAsFactors = FALSE)
symcol <- intersect(c("symbol", "Symbol", "gene", "name"), names(gt))[1]
bcol <- intersect(c("bnumber", "locus_tag", "b_number", "locus"), names(gt))[1]
ribo <- if (!is.na(symcol) && !is.na(bcol)) intersect(gt[[bcol]][grepl("^(rpl|rpm|rps)", gt[[symcol]])], genes) else character(0)
sets <- c(lapply(stats::setNames(sel, sel), function(nm) reg$members[[nm]]),
          if (length(ribo)) list(`ribosomal proteins` = match(ribo, genes)))
rows <- list()
for (nm in names(sets)) {
  y <- seq_along(genes) %in% sets[[nm]]
  for (k in 1:5) {
    L <- abs(sv$u[, k])
    w <- stats::wilcox.test(L[y], L[!y])
    rows[[length(rows) + 1L]] <- data.frame(set = nm, n = sum(y), pc = k, var_explained = ve[k],
      auroc = unname(w$statistic) / (sum(y) * sum(!y)), p = w$p.value)
  }
}
res <- do.call(rbind, rows)
res$p_bh <- stats::p.adjust(res$p, "BH")
utils::write.csv(res, file.path(out, "pc_regulon_auroc.csv"), row.names = FALSE)
# experiment loadings of PC1-3: which conditions drive them (growth phase, media)
desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
d1 <- desc$description[match(colnames(X), desc$experiment_name)]
for (k in 1:3) {
  o <- order(sv$v[, k])
  ends <- c(utils::head(o, 4), utils::tail(o, 4))
  cat(sprintf("PC%d (%.1f%%) extreme experiments:\n", k, ve[k]))
  for (j in ends) cat(sprintf("   %+.3f  %s  %s\n", sv$v[j, k], colnames(X)[j], substr(sub(".*(growth_phase:[^,]*).*", "\\1", d1[j]), 1, 60)))
}
cat("\nAUROC of |loading| (BH-adjusted p < 0.05 marked *):\n")
w <- reshape(res[, c("set", "n", "pc", "auroc", "p_bh")], idvar = c("set", "n"), timevar = "pc", direction = "wide")
for (i in seq_len(nrow(w)))
  cat(sprintf("  %-22s n=%4d  %s\n", w$set[i], w$n[i], paste(vapply(1:5, function(k)
    sprintf("%.2f%s", w[[paste0("auroc.", k)]][i], if (w[[paste0("p_bh.", k)]][i] < 0.05) "*" else " "), ""), collapse = "  ")))
cat(sprintf("variance explained PC1-5: %s\n", paste(sprintf("%.1f%%", ve[1:5]), collapse = ", ")))
