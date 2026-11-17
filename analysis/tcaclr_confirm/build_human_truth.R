#!/usr/bin/env Rscript
# Arm A truth build (docs/plans/tcaclr_confirm_plan.md section 4): EBI Complex Portal human (taxon 9606) on a given gene universe.
# Adapted COPY of analysis/complex_confirm/build_yeast_complex_truth.R: no rule change (recursive flattening, paralog molecule sets
# expanded and their pairs excluded from the pair universe, size 3-200 in-universe genes, overlap merge, n_eff, ribosome subunits kept
# separate as in the source with the xsub_excl / ribo_joined variants, ribosome list = merged complexes with >= 30% r-protein genes by the
# label-free regex ^RP[LS][0-9]+[A-Z]?$). Only the identifier map changes (UniProt reviewed human, primary gene symbol; upper-cased) and the
# four ribosomal-subunit complexes are located by Complex Portal id (40S CPX-5223, 60S CPX-5183, 28S CPX-5225, 39S CPX-5226).
# Builds and counts structure only; reads gene names, never an expression value; scores nothing.
# Usage: Rscript build_human_truth.R <universe.txt> <perturbed_genes.txt> <outdir> [prefix]
Sys.setlocale("LC_COLLATE", "C")
source("analysis/tcaclr_confirm/truth_lib.R")
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) >= 3)
u <- sort(unique(toupper(readLines(args[1])))); pert <- sort(unique(toupper(readLines(args[2])))); OUT <- args[3]
PFX <- if (length(args) >= 4) args[4] else "armA_complex_portal"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- "data/complex_portal_human_raw"
up <- utils::read.delim(file.path(RAW, "uniprot_9606_reviewed.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "")
sym <- toupper(up[[2]]); names(sym) <- up[[1]]; sym <- sym[nzchar(sym)]
cp <- utils::read.delim(file.path(RAW, "9606.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
cid <- sub("\\..*$", "", cp[[1]])
parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\([0-9.]*\\)$", "", v)); names(parts) <- cid
sets_acc <- unique(unlist(lapply(parts, function(v) grep("^\\[", v, value = TRUE))))
parts <- lapply(parts, function(v) unlist(lapply(v, function(p) if (grepl("^\\[", p)) strsplit(gsub("[][]", "", p), ",")[[1]] else p)))
flatten <- function(id, seen = character(0)) {
  if (id %in% seen) return(character(0))
  out <- character(0)
  for (p in parts[[id]]) {
    if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
    else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next
    else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p))
  }
  out
}
acc <- lapply(names(parts), flatten); names(acc) <- names(parts)
gen <- lapply(acc, function(a) sort(unique(sym[a[a %in% names(sym)]]))); gen <- gen[!duplicated(names(gen))]
paralog_sets <- lapply(sets_acc, function(p) sort(unique(sym[intersect(strsplit(gsub("[][]", "", p), ",")[[1]], names(sym))])))
# net4 rev-2 rule (one accession mapping to >= 2 ORFs): inapplicable here (one primary symbol per accession); counted for the record
multi_acc <- character(0)
paralog_sets <- unique(paralog_sets[lengths(paralog_sets) >= 2])
log_src <- c(paralog_sets = length(paralog_sets), complexes_in_file = length(parts), proteins_flattened_total = length(unlist(acc)),
             proteins_without_symbol = sum(!unlist(acc) %in% names(sym)), complexes_with_ge3_symbols_any_universe = sum(lengths(gen) >= 3),
             universe_genes = length(u), universe_genes_in_any_complex_ge3 = length(intersect(u, unique(unlist(gen[lengths(gen) >= 3])))))
RPRE <- "^RP[LS][0-9]+[A-Z]?$"; ribo_frac <- function(g) mean(grepl(RPRE, g))
sub_pairs <- list(cyto = c("(^|\\+)CPX-5223(\\+|$)", "(^|\\+)CPX-5183(\\+|$)"), mito = c("(^|\\+)CPX-5225(\\+|$)", "(^|\\+)CPX-5226(\\+|$)"))
H <- make_variants(gen, u, ribo_frac, sub_pairs, pert, "human_complex_portal")
cat("subunit complexes (merged names):\n"); print(H$sub_pairs)
cat("universe case collisions:", length(readLines(args[1])) - length(unique(toupper(readLines(args[1])))), "\n")
ex <- excluded_pairs(paralog_sets, H$sub_pairs, H$main)
xf <- file.path(OUT, paste0(PFX, "_excluded_pairs.tsv")); utils::write.table(ex, xf, sep = "\t", quote = FALSE, row.names = FALSE)
mem <- H$members; mem <- mem[order(mem$variant, mem$complex, mem$gene), ]
mf <- file.path(OUT, paste0(PFX, "_members.tsv")); utils::write.table(mem, mf, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.csv(H$counts, file.path(OUT, paste0(PFX, "_counts.csv")), row.names = FALSE)
utils::write.csv(H$sizes, file.path(OUT, paste0(PFX, "_sizes.csv")), row.names = FALSE)
utils::write.csv(data.frame(dataset = "human_complex_portal", t(log_src)), file.path(OUT, paste0(PFX, "_source_log.csv")), row.names = FALSE)
cat("source log:\n"); print(log_src); print(H$counts, row.names = FALSE)
cat("largest merged complexes:\n"); print(utils::head(H$sizes, 12), row.names = FALSE)
