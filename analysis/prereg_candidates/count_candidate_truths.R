#!/usr/bin/env Rscript
# STRUCTURE COUNTS ONLY for candidate untouched yeast truths on the Kemmeren deleteome gene universe.
# Adapted COPY of analysis/complex_confirm/build_yeast_complex_truth.R (original untouched): same flatten, paralog-set,
# size 3-200, overlap-merge and structure_of() code. No expression value is read (only the ORF row names of the
# deleteome matrix, passed in as a text file); nothing is scored against any truth.
# Usage: Rscript count_candidate_truths.R <kemmeren_orf_list.txt> <outdir>
Sys.setlocale("LC_COLLATE", "C")
args <- commandArgs(trailingOnly = TRUE); UFILE <- args[1]; OUT <- args[2]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- "data/complex_confirm_raw"
ORF_RE <- "^(Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4})$"
u <- sort(unique(readLines(UFILE))); stopifnot(all(grepl(ORF_RE, u)))
## ---- Complex Portal (verbatim logic of the builder) ----
up <- utils::read.delim(file.path(RAW, "uniprot", "uniprot_559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "")
orf_of <- lapply(strsplit(up[[2]], "[ ;/]+"), function(v) unique(grep(ORF_RE, v, value = TRUE))); names(orf_of) <- up$Entry
cp <- utils::read.delim(file.path(RAW, "cp", "559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
cid <- sub("\\..*$", "", cp[[1]])
parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\([0-9.]*\\)$", "", v)); names(parts) <- cid
sets_acc <- unique(unlist(lapply(parts, function(v) grep("^\\[", v, value = TRUE))))
parts <- lapply(parts, function(v) unlist(lapply(v, function(p) if (grepl("^\\[", p)) strsplit(gsub("[][]", "", p), ",")[[1]] else p)))
flatten <- function(id, seen = character(0)) { if (id %in% seen) return(character(0)); out <- character(0)
  for (p in parts[[id]]) { if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
    else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p)) }; out }
acc <- lapply(names(parts), flatten); names(acc) <- names(parts)
gen_cp <- lapply(acc, function(a) sort(unique(unlist(orf_of[a[a %in% names(orf_of)]])))); gen_cp <- gen_cp[!duplicated(names(gen_cp))]
paralog_sets <- lapply(sets_acc, function(p) sort(unique(unlist(orf_of[intersect(strsplit(gsub("[][]", "", p), ",")[[1]], names(orf_of))]))))
multi_acc <- unique(unlist(acc)); multi_acc <- multi_acc[multi_acc %in% names(orf_of) & lengths(orf_of[multi_acc]) >= 2]
paralog_sets <- c(paralog_sets, unname(lapply(multi_acc, function(a) sort(orf_of[[a]]))))
paralog_sets <- unique(paralog_sets[lengths(paralog_sets) >= 2])
sgd <- utils::read.delim("data/yeast/SGD_features.tab", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
rp <- unique(sgd$V4[sgd$V2 == "ORF" & grepl("^RP[LS][0-9]+[AB]?$", sgd$V5)])
ribo_frac <- function(g) mean(g %in% rp)
## ---- CYC2008 (Pu et al. 2009; ORF column used as is) ----
cy <- utils::read.delim("data/cyc2008_raw/CYC2008_complex.tab", stringsAsFactors = FALSE, quote = "", comment.char = "")
gen_cy <- lapply(split(cy$ORF, cy$Complex), function(v) sort(unique(v[grepl(ORF_RE, v)]))); names(gen_cy) <- paste0("CYC:", names(gen_cy))
## ---- filters, merge (verbatim) ----
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
## ---- structure (builder's structure_of) plus pair-set export ----
pairs_of <- function(s) { if (!length(s)) return(character(0))
  g <- sort(unique(unlist(s))); n <- length(g); cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }
  excl <- matrix(FALSE, n, n); for (ps in paralog_sets) { m <- stats::na.omit(match(ps, g)); if (length(m) >= 2) excl[m, m] <- TRUE }
  diag(excl) <- FALSE; cnt[excl] <- 0L; w <- which(cnt > 0 & upper.tri(cnt), arr.ind = TRUE); paste(g[w[, 1]], g[w[, 2]]) }
structure_of <- function(s) {
  if (!length(s)) return(list(w = numeric(0), n_pos = 0, n_genes = 0, n_eff = 0, largest = 0))
  g <- sort(unique(unlist(s))); n <- length(g); cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }
  excl <- matrix(FALSE, n, n); for (ps in paralog_sets) { m <- stats::na.omit(match(ps, g)); if (length(m) >= 2) excl[m, m] <- TRUE }
  diag(excl) <- FALSE; cnt[excl] <- 0L; pos <- cnt > 0 & upper.tri(cnt)
  w <- vapply(seq_along(s), function(k) { m <- match(s[[k]], g); sub <- cnt[m, m]; x <- sub[upper.tri(sub)]; sum(1 / x[x > 0]) }, 0)
  np <- sum(pos); ws <- w / max(np, 1)
  list(w = ws, n_pos = np, n_genes = n, n_eff = if (np) 1 / sum(ws^2) else 0, largest = if (np) max(ws) else 0) }
truths <- list(
  CP_559292_consumed_ref = merge_sets(filt(gen_cp, u)),
  CYC2008 = merge_sets(filt(gen_cy, u)),
  CP_plus_CYC2008_union = merge_sets(filt(c(gen_cp, gen_cy), u)))
cls <- function(s, k) switch(k, all = s, `3to6` = s[lengths(s) <= 6], `7to19` = s[lengths(s) >= 7 & lengths(s) <= 19], `20plus` = s[lengths(s) >= 20])
rows <- list()
for (tn in names(truths)) for (k in c("all", "3to6", "7to19", "20plus")) {
  s <- cls(truths[[tn]], k); st <- structure_of(s)
  isrp <- vapply(s, function(g) ribo_frac(g) >= 0.3, TRUE)
  rows[[length(rows) + 1]] <- data.frame(truth = tn, size_class = k, n_complexes = length(s), n_genes = st$n_genes, n_pos_pairs = st$n_pos,
    n_eff = round(st$n_eff, 2), largest_share = round(st$largest, 3), ribo_share = round(sum(st$w[isrp]), 3), median_size = if (length(s)) median(lengths(s)) else NA)
}
cnt <- do.call(rbind, rows); utils::write.csv(cnt, file.path(OUT, "candidate_truth_counts.csv"), row.names = FALSE)
## ---- raw (unmerged, size 3-200) counts ----
raw <- data.frame(truth = c("CP_559292", "CYC2008"), n_complexes_in_source = c(length(gen_cp), length(gen_cy)),
  n_ge3_in_universe = c(length(filt(gen_cp, u)), length(filt(gen_cy, u))),
  universe_genes_in_any_complex_ge3 = c(length(unique(unlist(filt(gen_cp, u)))), length(unique(unlist(filt(gen_cy, u))))))
utils::write.csv(raw, file.path(OUT, "candidate_truth_raw_counts.csv"), row.names = FALSE)
## ---- overlap with the consumed Complex Portal truth (same universe) ----
ref <- truths$CP_559292_consumed_ref; Pref <- pairs_of(ref); gref <- sort(unique(unlist(ref)))
jmax <- function(a, B) max(vapply(B, function(b) length(intersect(a, b)) / length(union(a, b)), 0))
ov <- list()
for (tn in c("CYC2008", "CP_plus_CYC2008_union")) {
  s <- truths[[tn]]; P <- pairs_of(s); g <- sort(unique(unlist(s))); jm <- vapply(s, jmax, 0, B = ref); nov <- setdiff(P, Pref)
  # n_eff of the pairs NOT in the consumed truth (fractional assignment to complexes holding them)
  sn <- lapply(s, function(x) x)  # complexes; positives restricted below
  wk <- vapply(s, function(x) { pk <- if (length(x) >= 2) { cb <- utils::combn(sort(x), 2); paste(cb[1, ], cb[2, ]) } else character(0); sum(pk %in% nov) }, 0)
  ov[[tn]] <- data.frame(truth = tn, n_complexes = length(s), genes = length(g), genes_in_CP_ref = length(intersect(g, gref)), gene_jaccard = round(length(intersect(g, gref)) / length(union(g, gref)), 3),
    n_pos = length(P), n_pos_shared_with_CP = length(intersect(P, Pref)), share_of_truth_pos_in_CP = round(length(intersect(P, Pref)) / length(P), 3),
    share_of_CP_pos_in_truth = round(length(intersect(P, Pref)) / length(Pref), 3), pair_jaccard = round(length(intersect(P, Pref)) / length(union(P, Pref)), 3),
    n_pos_not_in_CP = length(nov), complexes_best_jaccard_ge_0.5 = sum(jm >= 0.5), complexes_best_jaccard_ge_0.8 = sum(jm >= 0.8), complexes_best_jaccard_lt_0.5 = sum(jm < 0.5),
    median_best_jaccard = round(median(jm), 3), neff_of_novel_pairs_complexwise = round(sum(wk)^2 / sum(wk^2), 2))
}
ov <- do.call(rbind, ov); utils::write.csv(ov, file.path(OUT, "candidate_truth_overlap_with_consumed.csv"), row.names = FALSE)
# CYC2008-only vs CP: reciprocal overlap on the CYC2008 complexes' side
cat("paralog sets:", length(paralog_sets), " universe:", length(u), "\n")
print(raw, row.names = FALSE); print(cnt, row.names = FALSE); print(ov, row.names = FALSE)
## ---- label-free structure: how many deleted genes (mutant column names only; no expression value) are members of each truth's complexes ----
if (length(args) >= 3) {
  hdr <- strsplit(readLines(args[3], n = 1), "\t")[[1]][-1]
  nm <- toupper(sub("[-.][0-9]+$", "", hdr))
  sg <- unique(sgd[sgd$V2 == "ORF" & sgd$V5 != "", c("V4", "V5")]); map <- setNames(sg$V4, toupper(sg$V5))
  del <- unname(map[nm]); cat("mutant columns:", length(hdr), " mapped to ORF via SGD standard name:", sum(!is.na(del)), "\n")
  res <- do.call(rbind, lapply(names(truths), function(tn) { s <- truths[[tn]]; g <- unique(unlist(s)); big <- unique(unlist(s[lengths(s) >= 20]))
    data.frame(truth = tn, mutants_mapped = sum(!is.na(del)), deleted_gene_in_truth_complex = sum(del %in% g, na.rm = TRUE),
               deleted_gene_in_20plus_complex = sum(del %in% big, na.rm = TRUE), truth_genes = length(g), truth_genes_deleted_in_panel = sum(g %in% del)) }))
  print(res, row.names = FALSE); utils::write.csv(res, file.path(OUT, "candidate_truth_deleted_gene_membership.csv"), row.names = FALSE)
}
