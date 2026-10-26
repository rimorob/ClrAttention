#!/usr/bin/env Rscript
# Builds the E. coli protein-complex truth exactly as docs/preregistration_clr_logpt.md ("E. coli protein
# complexes"). BUILD ONLY: reads gene universes (expression row names, no values used beyond the
# constant-row filter of baselines.R), the complex source, the operon map. Scores nothing.
#
# Source actually used: EBI Complex Portal taxon 83333 (preregistered availability fallback; EcoCyc
# proteins.dat requires a BioCyc login/licence, see docs/preregistration_clr_logpt.md change log).
#
# Inputs (raw, not committed; hashes recorded in the change log):
#   $RAW/83333.tsv          Complex Portal complextab
#   $RAW/uniprot_83333.tsv  UniProt accession -> b-number table (taxon 83333)
# Outputs (docs/clr_logpt_selections/):
#   ecoli_complexes_members.tsv   long: variant, dataset, complex, bnumber   (SHA-256 gated)
#   ecoli_complexes_counts.csv    counts per variant x dataset x operon handling
#   ecoli_complexes.sha256        SHA-256 of the member file
# Usage: Rscript analysis/build_ecoli_complex_truth.R <RAW dir>
args <- commandArgs(trailingOnly = TRUE); RAW <- args[1]
suppressPackageStartupMessages(library(Matrix))
source("analysis/regulons.R")
OUT <- "docs/clr_logpt_selections"
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
gt <- utils::read.delim("analysis/ecoli_k12_genes.tsv", stringsAsFactors = FALSE, quote = "", comment.char = "")

## ---- gene universes (rows of the expression matrices, as baselines.R) ----
rng0 <- function(X) apply(X, 1, function(x) diff(range(x))) > 0
m3d <- "data/E_coli_v4_Build_6"
raw <- utils::read.delim(file.path(m3d, "avg_E_coli_v4_Build_6_exps466probes4297.tab"), check.names = FALSE)
X <- as.matrix(raw[, -1]); rn <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
universe <- list(m3d = unique(rn[rng0(X)]))
X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
gn <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
universe$ecomac <- unique(gn[rng0(X)]); rm(X, raw)

## ---- UniProt -> b-number ----
up <- utils::read.delim(file.path(RAW, "uniprot_83333.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "")
bn_of <- lapply(strsplit(up$Gene.Names..ordered.locus., "[ ;]+"), function(v) unique(grep("^b[0-9]{4}$", v, value = TRUE)))
names(bn_of) <- up$Entry

## ---- complexes: flatten recursively ----
cp <- utils::read.delim(file.path(RAW, "83333.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
cid <- sub("\\..*$", "", cp[[1]])
parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\(.*$", "", v))
names(parts) <- cid
flatten <- function(id, seen = character(0)) {
  if (id %in% seen) return(character(0))
  out <- character(0)
  for (p in parts[[id]]) {
    if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
    else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next     # non-protein: dropped
    else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p))       # isoform / chain suffix -> parent accession
  }
  out
}
acc <- lapply(names(parts), flatten); names(acc) <- names(parts)
gen <- lapply(acc, function(a) sort(unique(unlist(bn_of[a[a %in% names(bn_of)]]))))
n_unmapped <- sum(!unlist(acc) %in% names(bn_of))
gen <- gen[!duplicated(names(gen))]
log_src <- c(complexes_in_file = length(parts), proteins_flattened_total = length(unlist(acc)),
             proteins_without_bnumber = n_unmapped)

## ---- size filter, overlap merge, variants ----
filt <- function(sets, u, lo = 3L, hi = 200L) {
  s <- lapply(sets, function(g) intersect(g, u)); s[lengths(s) >= lo & lengths(s) <= hi]
}
merge_sets <- function(s) {
  if (length(s) == 0) return(s)
  s <- s[order(names(s))]
  key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]           # exact duplicates
  drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j)
    length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
  s <- s[!drop]                                                                     # proper subsets
  if (length(s) < 2) return(s)
  all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
  inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
  out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
  names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), "")
  out
}
rp <- gt$bnumber[grepl("^rpl[A-Y]$|^rpm[A-J]$|^rps[A-V]$", gt$symbol)]
ribo_frac <- function(g) mean(g %in% rp)
pairs_count <- function(s, genes, with_operon) {
  if (length(s) == 0) return(c(n_pairs = 0, n_pos = 0, n_complexes_pos = 0))
  ix <- lapply(s, match, genes)
  reg <- list(members = ix, class = stats::setNames(rep("complex", length(ix)), names(ix)))
  noop <- if (with_operon) stats::setNames(character(0), character(0)) else operon_of
  p <- comembership_pairs(reg, noop, genes)
  # complexes contributing >= 1 positive (cross-operon) pair
  contrib <- vapply(ix, function(m) { gm <- genes[m]; o <- operon_of[gm]
    if (!with_operon) length(m) >= 2 && (any(is.na(o)) || length(unique(o)) > 1) else length(m) >= 2 }, TRUE)
  c(n_pairs = length(p$label), n_pos = sum(p$label), n_complexes_pos = sum(contrib))
}

rows <- list(); counts <- list(); fallback <- list()
for (ds in names(universe)) {
  u <- universe[[ds]]
  s3 <- filt(gen, u); main <- merge_sets(s3)
  variants <- list(main = main, size2 = merge_sets(filt(gen, u, lo = 2L)), nomerge = s3,
                   noribo = main[!vapply(main, function(g) ribo_frac(g) >= 0.3, TRUE)],
                   rp_module = list(rp_module = intersect(rp, u)))
  for (v in names(variants)) {
    s <- variants[[v]]
    for (k in names(s)) rows[[length(rows) + 1]] <- data.frame(variant = v, dataset = ds, complex = k, bnumber = s[[k]])
    for (wo in c(FALSE, TRUE)) {
      pc <- pairs_count(s, u, wo)
      counts[[length(counts) + 1]] <- data.frame(dataset = ds, variant = v, same_operon_pairs = if (wo) "included" else "excluded",
        n_complexes = length(s), n_genes = length(unique(unlist(s))), t(pc), row.names = NULL)
    }
  }
  pc <- pairs_count(main, u, FALSE)
  fallback[[ds]] <- !(pc["n_pos"] >= 200 && pc["n_complexes_pos"] >= 10)
}
mem <- do.call(rbind, rows)
mem <- mem[order(mem$variant, mem$dataset, mem$complex, mem$bnumber), ]
mf <- file.path(OUT, "ecoli_complexes_members.tsv")
utils::write.table(mem, mf, sep = "\t", quote = FALSE, row.names = FALSE)
cnt <- do.call(rbind, counts)
cnt$operon_fallback_required <- ifelse(cnt$variant == "main" & cnt$same_operon_pairs == "excluded", unlist(fallback)[cnt$dataset], NA)
utils::write.csv(cnt, file.path(OUT, "ecoli_complexes_counts.csv"), row.names = FALSE)
h <- digest::digest(file = mf, algo = "sha256")
writeLines(paste(h, "ecoli_complexes_members.tsv"), file.path(OUT, "ecoli_complexes.sha256"))
cat("source log:\n"); print(log_src); cat("universe sizes:", lengths(universe), "\n")
print(cnt, row.names = FALSE); cat("fallback required:", paste(names(fallback), unlist(fallback)), "\nSHA-256", h, "\n")
