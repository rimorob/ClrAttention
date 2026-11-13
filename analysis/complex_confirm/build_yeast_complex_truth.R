#!/usr/bin/env Rscript
# Builds the S. cerevisiae protein-complex truth for docs/plans/complex_confirm_design.md.
# BUILD AND COUNT ONLY: reads gene universes (gene names of the expression files; no expression value is
# read), the complex source and the UniProt map. Scores nothing; no readout, network or score is loaded.
#
# Rules copied from analysis/build_ecoli_complex_truth.R (the 26 Oct preregistered E. coli truth):
# recursive flattening of sub-complexes, non-protein participants dropped, isoform/chain suffixes to the
# parent accession, size 3-200 in-universe genes, overlap merge (drop exact duplicates and proper subsets,
# then complete linkage on 1 - Jaccard cut at 0.5, union). No operon rule in yeast.
# Ribosome: merged complexes with >= 30% cytoplasmic r-protein ORFs (SGD standard names ^RP[LS][0-9]+[AB]?$,
# the rule of docs/preregistration_reciprocity_yeast.md) form the pre-declared "ribosome-dominated" list.
#
# Inputs (raw, gitignored; SHA-256 in data/complex_confirm_raw/SHA256SUMS):
#   $RAW/cp/559292.tsv                 EBI Complex Portal complextab, S. cerevisiae S288C (curated)
#   $RAW/uniprot/uniprot_559292.tsv    UniProt accession -> ordered locus (ORF)
#   data/yeast/SGD_features.tab        ORF -> standard name (for the r-protein rule)
# Universes: DREAM5 net4 (data/dream5/synapse/net4/net4_gene_ids.tsv, decoys removed; candidate A) and the
# Kemmeren responsive-700 deleteome rows (data/yeast/deleteome_responsive_700.tsv, first column; candidate B).
# Outputs (docs/complex_confirm_selections/):
#   yeast_complexes_members.tsv   long: variant, dataset, complex, orf   (SHA-256 gated at the run stage)
#   yeast_complexes_counts.csv    counts per variant x dataset (genes, pairs, positives, n_eff, shares)
#   yeast_complexes_sizes.csv     per merged complex: genes, positive-pair weight, r-protein fraction
#   yeast_complexes.sha256
# Usage: Rscript analysis/complex_confirm/build_yeast_complex_truth.R data/complex_confirm_raw
args <- commandArgs(trailingOnly = TRUE); RAW <- if (length(args)) args[1] else "data/complex_confirm_raw"
Sys.setlocale("LC_COLLATE", "C")
OUT <- "docs/complex_confirm_selections"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ORF_RE <- "^(Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4})$"

## ---- universes (gene names only) ----
gid <- utils::read.delim("data/dream5/synapse/net4/net4_gene_ids.tsv", stringsAsFactors = FALSE, comment.char = "", quote = "")
names(gid) <- c("id", "name")
universe <- list(dream5_net4 = sort(unique(gid$name[!grepl("^decoy", gid$name)])))
kem <- sub("\t.*$", "", readLines("data/yeast/deleteome_responsive_700.tsv"))[-1]
universe$kemmeren_700 <- sort(unique(toupper(kem)))

## ---- UniProt -> ORF ----
up <- utils::read.delim(file.path(RAW, "uniprot", "uniprot_559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "")
orf_of <- lapply(strsplit(up[[2]], "[ ;/]+"), function(v) unique(grep(ORF_RE, v, value = TRUE)))
names(orf_of) <- up$Entry

## ---- complexes: flatten recursively (verbatim logic of the E. coli builder) ----
cp <- utils::read.delim(file.path(RAW, "cp", "559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
cid <- sub("\\..*$", "", cp[[1]])
# Participants: "ACC(stoich)" or a paralog molecule set "[ACC1,ACC2](stoich)". Only the two cytoplasmic ribosome
# subunits (CPX-1599, CPX-1601) use sets in this release. Rule (pre-declared, before any score): a set is
# expanded to all its members (interchangeable paralogous subunits are all members), and pairs of ORFs from the
# same set are excluded from the pair universe (near-identical paralogs, trivially co-measured on arrays; the
# yeast analogue of the same-operon exclusion of the E. coli truth). The E. coli file 83333.tsv has no sets.
parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\([0-9.]*\\)$", "", v))
names(parts) <- cid
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
gen <- lapply(acc, function(a) sort(unique(unlist(orf_of[a[a %in% names(orf_of)]]))))
gen <- gen[!duplicated(names(gen))]
paralog_sets <- lapply(sets_acc, function(p) sort(unique(unlist(orf_of[intersect(strsplit(gsub("[][]", "", p), ",")[[1]], names(orf_of))]))))
paralog_sets <- paralog_sets[lengths(paralog_sets) >= 2]
log_src <- c(paralog_sets = length(paralog_sets), complexes_in_file = length(parts), proteins_flattened_total = length(unlist(acc)),
             proteins_without_orf = sum(!unlist(acc) %in% names(orf_of)),
             complexes_with_ge3_orfs_any_universe = sum(lengths(gen) >= 3))

## ---- r-protein ORFs ----
sgd <- utils::read.delim("data/yeast/SGD_features.tab", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
rp <- unique(sgd$V4[sgd$V2 == "ORF" & grepl("^RP[LS][0-9]+[AB]?$", sgd$V5)])
ribo_frac <- function(g) mean(g %in% rp)

## ---- size filter, overlap merge (verbatim) ----
filt <- function(sets, u, lo = 3L, hi = 200L) { s <- lapply(sets, function(g) intersect(g, u)); s[lengths(s) >= lo & lengths(s) <= hi] }
merge_sets <- function(s) {
  if (length(s) == 0) return(s)
  s <- s[order(names(s))]
  key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
  drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j)
    length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
  s <- s[!drop]
  if (length(s) < 2) return(s)
  all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
  inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
  out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
  names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), "")
  out
}

## ---- structure counts (no scores) ----
# positive pairs with fractional assignment to the merged complexes that contain them;
# n_eff = 1 / sum(w_k^2) (inverse Herfindahl of positive-pair shares); components of the gene-sharing graph.
structure_of <- function(s) {
  if (!length(s)) return(list(w = numeric(0), n_pos = 0, n_pairs = 0, n_genes = 0, comp = 0))
  g <- sort(unique(unlist(s))); n <- length(g)
  cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }
  excl <- matrix(FALSE, n, n)
  for (ps in paralog_sets) { m <- stats::na.omit(match(ps, g)); if (length(m) >= 2) excl[m, m] <- TRUE }
  diag(excl) <- FALSE; cnt[excl] <- 0L
  ut <- upper.tri(cnt); pos <- cnt > 0 & ut
  w <- vapply(seq_along(s), function(k) { m <- match(s[[k]], g); sub <- cnt[m, m]; x <- sub[upper.tri(sub)]; sum(1 / x[x > 0]) }, 0)
  adj <- crossprod(vapply(s, function(x) g %in% x, logical(n))) > 0
  cl <- seq_along(s); repeat { new <- vapply(seq_along(s), function(k) min(cl[adj[k, ]]), 1L); if (all(new == cl)) break; cl <- new }
  list(w = stats::setNames(w, names(s)), n_pos = sum(pos), n_pairs = n * (n - 1) / 2 - sum(excl & ut), n_genes = n,
       comp = length(unique(cl)), n_excl = sum(excl & ut))
}

rows <- list(); counts <- list(); sizes <- list()
for (ds in names(universe)) {
  u <- universe[[ds]]
  s3 <- filt(gen, u); main <- merge_sets(s3)
  isrp <- vapply(main, function(g) ribo_frac(g) >= 0.3, TRUE)
  variants <- list(main = main, noribo = main[!isrp], ribo_list = main[isrp], nomerge = s3,
                   size2 = merge_sets(filt(gen, u, lo = 2L)),
                   small_3to6 = main[lengths(main) <= 6], large_ge20 = main[lengths(main) >= 20],
                   rp_module = list(rp_module = intersect(rp, u)))
  for (v in names(variants)) {
    s <- variants[[v]]
    for (k in names(s)) rows[[length(rows) + 1]] <- data.frame(variant = v, dataset = ds, complex = k, orf = s[[k]])
    st <- structure_of(s); w <- st$w / max(st$n_pos, 1)
    counts[[length(counts) + 1]] <- data.frame(dataset = ds, variant = v, n_universe = length(u), n_complexes = length(s),
      n_genes = st$n_genes, n_pairs = st$n_pairs, n_paralog_pairs_excluded = st$n_excl, n_pos = st$n_pos, base_rate = st$n_pos / max(st$n_pairs, 1),
      n_complexes_pos = sum(lengths(s) >= 2), n_components = st$comp, n_eff = if (length(w)) 1 / sum(w^2) else 0,
      largest_share = if (length(w)) max(w) else 0,
      ribo_dominated_share = if (v == "main") sum(w[isrp]) else NA_real_)
    if (v == "main") sizes[[ds]] <- data.frame(dataset = ds, complex = names(s), n_genes = lengths(s),
      pos_weight = st$w, share = w, rp_frac = vapply(s, ribo_frac, 0), ribo_dominated = isrp, row.names = NULL)
  }
}
mem <- do.call(rbind, rows); mem <- mem[order(mem$variant, mem$dataset, mem$complex, mem$orf), ]
mf <- file.path(OUT, "yeast_complexes_members.tsv")
utils::write.table(mem, mf, sep = "\t", quote = FALSE, row.names = FALSE)
cnt <- do.call(rbind, counts); utils::write.csv(cnt, file.path(OUT, "yeast_complexes_counts.csv"), row.names = FALSE)
sz <- do.call(rbind, sizes); sz <- sz[order(sz$dataset, -sz$share), ]
utils::write.csv(sz, file.path(OUT, "yeast_complexes_sizes.csv"), row.names = FALSE)
h <- digest::digest(file = mf, algo = "sha256")
writeLines(paste(h, "yeast_complexes_members.tsv"), file.path(OUT, "yeast_complexes.sha256"))
allorf <- unique(unlist(gen[lengths(gen) >= 3]))
cov <- vapply(universe, function(u) mean(allorf %in% u), 0)
utils::write.csv(data.frame(dataset = names(universe), n_universe = lengths(universe), complex_orfs_size_ge3 = length(allorf),
  coverage = cov, t(log_src)[rep(1, length(universe)), , drop = FALSE], row.names = NULL), file.path(OUT, "yeast_complexes_source_log.csv"), row.names = FALSE)
cat("coverage of complex ORFs (complexes with >= 3 ORFs) by universe:", round(cov, 4), "\n")
cat("source log:\n"); print(log_src); cat("universe sizes:", lengths(universe), "\n"); cat("r-protein ORFs:", length(rp), "\n")
print(cnt, row.names = FALSE); cat("SHA-256", h, "\n")
cat("largest merged complexes (net4):\n"); print(utils::head(sz[sz$dataset == "dream5_net4", ], 12), row.names = FALSE)
