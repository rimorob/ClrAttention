#!/usr/bin/env Rscript
# STRUCTURE COUNTS ONLY: EBI Complex Portal human (9606) on a PROXY gene universe (gene names of an existing consumed
# Perturb-seq matrix; no expression value read), plus overlap with the consumed CORUM-merged truth (analysis/human_truths.R::corum_sets).
# Same flatten / size 3-200 / overlap merge / structure_of code as analysis/complex_confirm/build_yeast_complex_truth.R (copy). Scores nothing.
# Usage: Rscript count_human_complex_portal.R <h5ad with /var/gene_name> <outdir>
Sys.setlocale("LC_COLLATE", "C")
args <- commandArgs(trailingOnly = TRUE); H5 <- args[1]; OUT <- args[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source("analysis/human_truths.R", local = TRUE)  # corum_sets (functions only)
vg <- rhdf5::h5read(H5, "/var/gene_name"); cats <- tryCatch(rhdf5::h5read(H5, "/var/__categories/gene_name"), error = function(e) NULL)
u <- if (is.list(vg)) as.character(vg$categories)[vg$codes + 1L] else if (!is.null(cats)) as.character(cats)[vg + 1L] else as.character(vg)
u <- sort(unique(toupper(u)))
up <- utils::read.delim("data/complex_portal_human_raw/uniprot_9606_reviewed.tsv", stringsAsFactors = FALSE, quote = "", comment.char = "")
sym <- toupper(up[[2]]); names(sym) <- up[[1]]
cp <- utils::read.delim("data/complex_portal_human_raw/9606.tsv", stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
cid <- cp[[1]]
parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\([0-9.]*\\)$", "", v)); names(parts) <- cid
parts <- lapply(parts, function(v) unlist(lapply(v, function(p) if (grepl("^\\[", p)) strsplit(gsub("[][]", "", p), ",")[[1]] else p)))
flatten <- function(id, seen = character(0)) { if (id %in% seen) return(character(0)); out <- character(0)
  for (p in parts[[id]]) { p <- sub("^(CPX-[0-9]+).*$", "\\1", p)
    if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
    else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p)) }; out }
acc <- lapply(names(parts), flatten); names(acc) <- names(parts)
gen <- lapply(acc, function(a) sort(unique(na.omit(sym[a[a %in% names(sym)]])))); gen <- gen[!duplicated(names(gen))]
filt <- function(sets, u, lo = 3L, hi = 200L) { s <- lapply(sets, function(g) intersect(g, u)); s[lengths(s) >= lo & lengths(s) <= hi] }
merge_sets <- function(s) {
  if (length(s) == 0) return(s)
  s <- s[order(names(s))]; key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
  drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j) length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
  s <- s[!drop]; if (length(s) < 2) return(s)
  all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
  inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
  out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
  names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), ""); out }
pairs_of <- function(s) { if (!length(s)) return(character(0)); g <- sort(unique(unlist(s))); n <- length(g); cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }; w <- which(cnt > 0 & upper.tri(cnt), arr.ind = TRUE); paste(g[w[, 1]], g[w[, 2]]) }
structure_of <- function(s) { if (!length(s)) return(list(n_pos = 0, n_genes = 0, n_eff = 0, largest = 0))
  g <- sort(unique(unlist(s))); n <- length(g); cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }; pos <- cnt > 0 & upper.tri(cnt)
  w <- vapply(seq_along(s), function(k) { m <- match(s[[k]], g); sub <- cnt[m, m]; x <- sub[upper.tri(sub)]; sum(1 / x[x > 0]) }, 0)
  np <- sum(pos); ws <- w / max(np, 1); list(n_pos = np, n_genes = n, n_eff = if (np) 1 / sum(ws^2) else 0, largest = if (np) max(ws) else 0) }
cls <- function(s, k) switch(k, all = s, `3to6` = s[lengths(s) <= 6], `7to19` = s[lengths(s) >= 7 & lengths(s) <= 19], `20plus` = s[lengths(s) >= 20])
cps <- merge_sets(filt(gen, u)); cor_raw <- corum_sets(u, 3L, 200L)$members
cor_raw <- lapply(cor_raw, function(ix) u[ix]); corm <- merge_sets(cor_raw)
rows <- list()
for (tn in c("ComplexPortal_9606", "CORUM_consumed_ref")) for (k in c("all", "3to6", "7to19", "20plus")) {
  s <- cls(if (tn == "ComplexPortal_9606") cps else corm, k); st <- structure_of(s)
  rows[[length(rows) + 1]] <- data.frame(truth = tn, size_class = k, n_complexes = length(s), n_genes = st$n_genes, n_pos_pairs = st$n_pos, n_eff = round(st$n_eff, 2), largest_share = round(st$largest, 3)) }
cnt <- do.call(rbind, rows); utils::write.csv(cnt, file.path(OUT, "human_cp_counts.csv"), row.names = FALSE)
P <- pairs_of(cps); R <- pairs_of(corm); gP <- unique(unlist(cps)); gR <- unique(unlist(corm))
jm <- vapply(cps, function(a) max(vapply(corm, function(b) length(intersect(a, b)) / length(union(a, b)), 0)), 0)
ov <- data.frame(universe_proxy_genes = length(u), cp_raw_complexes_in_source = length(gen), cp_ge3_in_universe = length(filt(gen, u)), corum_raw_ge3 = length(cor_raw),
  gene_jaccard = round(length(intersect(gP, gR)) / length(union(gP, gR)), 3), n_pos_cp = length(P), n_pos_corum = length(R), shared = length(intersect(P, R)),
  share_cp_pos_in_corum = round(length(intersect(P, R)) / length(P), 3), share_corum_pos_in_cp = round(length(intersect(P, R)) / length(R), 3),
  pair_jaccard = round(length(intersect(P, R)) / length(union(P, R)), 3), cp_complexes_best_jaccard_ge_0.5 = sum(jm >= 0.5), cp_complexes_best_jaccard_lt_0.5 = sum(jm < 0.5))
utils::write.csv(ov, file.path(OUT, "human_cp_overlap_with_corum.csv"), row.names = FALSE)
print(cnt, row.names = FALSE); print(t(ov))
