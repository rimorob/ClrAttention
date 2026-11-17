#!/usr/bin/env Rscript
# Arm B truth build (docs/plans/tcaclr_confirm_plan.md section 4): CYC2008 on the Kemmeren deleteome universe.
# Adapted COPY of analysis/complex_confirm/build_yeast_complex_truth.R and analysis/prereg_candidates/count_candidate_truths.R:
# no rule change (flatten, paralog sets and paralog-pair exclusion from the Complex Portal yeast release exactly as net4, size 3-200,
# overlap merge, n_eff, subunit variants, ribosome list at >= 30% r-protein ORFs by ^RP[LS][0-9]+[AB]?$). Builds and counts structure only:
# reads the ORF row names and mutant column names of the deleteome matrix, never an expression value; scores nothing.
# Usage: Rscript analysis/tcaclr_confirm/build_cyc2008_truth.R
Sys.setlocale("LC_COLLATE", "C")
source("analysis/tcaclr_confirm/truth_lib.R")
OUT <- "docs/tcaclr_confirm_selections"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RAW <- "data/complex_confirm_raw"; MAT <- "data/yeast/deleteome_all_1484.tsv"
ORF_RE <- "^(Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4})$"
## universe and panel (names only)
u <- sort(unique(utils::read.delim(MAT, header = TRUE, colClasses = c("character", rep("NULL", 1283)), check.names = FALSE)[[1]]))
stopifnot(length(u) == 6169, all(grepl(ORF_RE, u)))
hdr <- strsplit(readLines(MAT, n = 1), "\t")[[1]][-1]; stopifnot(length(hdr) == 1283)
sgd <- utils::read.delim("data/yeast/SGD_features.tab", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
s <- sgd[grepl(ORF_RE, sgd$V4), ]
nm <- stats::setNames(s$V4, toupper(s$V5)); nm <- nm[nzchar(names(nm))]               # standard name, then systematic, then unique alias
al <- do.call(rbind, lapply(seq_len(nrow(s)), function(k) if (nzchar(s$V6[k])) data.frame(a = strsplit(toupper(s$V6[k]), "|", fixed = TRUE)[[1]], orf = s$V4[k]) else NULL))
amb <- unique(al$a[duplicated(al$a)]); al <- al[!al$a %in% amb, ]
x <- toupper(hdr); del <- unname(nm[x]); i <- is.na(del) & x %in% s$V4; del[i] <- x[i]; i <- is.na(del); del[i] <- al$orf[match(x[i], al$a)]
cat("mutant columns:", length(hdr), " mapped to ORF:", sum(!is.na(del)), " unique deleted ORFs:", length(unique(na.omit(del))), "\n")
perturbed <- sort(unique(na.omit(del)))
utils::write.table(data.frame(column = hdr, deleted_orf = del), file.path(OUT, "armB_panel_deleted_orfs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## Complex Portal yeast (consumed reference; also the source of the net4 paralog sets)
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
stopifnot(length(paralog_sets) == 57)
## CYC2008 (Pu et al. 2009; ORF column as is)
cy <- utils::read.delim("data/cyc2008_raw/CYC2008_complex.tab", stringsAsFactors = FALSE, quote = "", comment.char = "")
gen_cy <- lapply(split(cy$ORF, cy$Complex), function(v) sort(unique(v[grepl(ORF_RE, v)]))); names(gen_cy) <- paste0("CYC:", names(gen_cy))
stopifnot(length(gen_cy) == 408)
rp <- unique(sgd$V4[sgd$V2 == "ORF" & grepl("^RP[LS][0-9]+[AB]?$", sgd$V5)])
ribo_frac <- function(g) mean(g %in% rp)
rx <- function(n) paste0("(^|\\+)CYC:", n, "(\\+|$)")
sub_pairs <- list(cyto = c(rx("cytoplasmic ribosomal large subunit"), rx("cytoplasmic ribosomal small subunit")),
                  mito = c(rx("mitochondrial ribosomal large subunit"), rx("mitochondrial ribosomal small subunit")))
B <- make_variants(gen_cy, u, ribo_frac, sub_pairs, perturbed, "kemmeren_cyc2008")
cat("subunit complexes found:\n"); print(B$sub_pairs)
## consumed Complex Portal yeast truth on the same universe (reference for overlap; same code)
Cp <- make_variants(gen_cp, u, ribo_frac, list(cyto = c("^CPX-1599$", "^CPX-1601$"), mito = c("^CPX-1603$", "^CPX-1602$")), perturbed, "kemmeren_cp_consumed")
## overlap with the consumed truth (counts only)
pairs_of <- function(s) { if (!length(s)) return(character(0))
  g <- sort(unique(unlist(s))); n <- length(g); cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }
  excl <- matrix(FALSE, n, n); for (ps in paralog_sets) { m <- stats::na.omit(match(ps, g)); if (length(m) >= 2) excl[m, m] <- TRUE }
  diag(excl) <- FALSE; cnt[excl] <- 0L; w <- which(cnt > 0 & upper.tri(cnt), arr.ind = TRUE); paste(g[w[, 1]], g[w[, 2]]) }
P <- pairs_of(B$main); R <- pairs_of(Cp$main); gB <- sort(unique(unlist(B$main))); gR <- sort(unique(unlist(Cp$main)))
jmax <- function(a, S) max(vapply(S, function(b) length(intersect(a, b)) / length(union(a, b)), 0))
jm <- vapply(B$main, jmax, 0, S = Cp$main); nov <- setdiff(P, R)
wk <- vapply(B$main, function(x) { pk <- if (length(x) >= 2) { cb <- utils::combn(sort(x), 2); paste(cb[1, ], cb[2, ]) } else character(0); sum(pk %in% nov) }, 0)
ov <- data.frame(truth = "CYC2008", n_complexes = length(B$main), genes = length(gB), genes_in_CP_ref = length(intersect(gB, gR)),
  gene_jaccard = round(length(intersect(gB, gR)) / length(union(gB, gR)), 3), n_pos = length(P), n_pos_shared_with_CP = length(intersect(P, R)),
  share_of_truth_pos_in_CP = round(length(intersect(P, R)) / length(P), 3), share_of_CP_pos_in_truth = round(length(intersect(P, R)) / length(R), 3),
  pair_jaccard = round(length(intersect(P, R)) / length(union(P, R)), 3), n_pos_not_in_CP = length(nov),
  complexes_best_jaccard_ge_0.5 = sum(jm >= 0.5), complexes_best_jaccard_lt_0.5 = sum(jm < 0.5), neff_of_novel_pairs_complexwise = round(sum(wk)^2 / sum(wk^2), 2),
  cp_consumed_complexes = length(Cp$main), cp_consumed_pos = length(R), cp_consumed_n_eff = round(Cp$counts$n_eff[Cp$counts$variant == "main"], 2))
## outputs
ex <- excluded_pairs(paralog_sets, B$sub_pairs, B$main)
xf <- file.path(OUT, "armB_cyc2008_excluded_pairs.tsv"); utils::write.table(ex, xf, sep = "\t", quote = FALSE, row.names = FALSE)
mem <- B$members; mem <- mem[order(mem$variant, mem$complex, mem$orf), ]
mf <- file.path(OUT, "armB_cyc2008_members.tsv"); utils::write.table(mem, mf, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.csv(B$counts, file.path(OUT, "armB_cyc2008_counts.csv"), row.names = FALSE)
utils::write.csv(B$sizes, file.path(OUT, "armB_cyc2008_sizes.csv"), row.names = FALSE)
utils::write.csv(ov, file.path(OUT, "armB_overlap_with_consumed_cp_yeast.csv"), row.names = FALSE)
utils::write.csv(Cp$counts[Cp$counts$variant == "main", ], file.path(OUT, "armB_consumed_cp_yeast_counts.csv"), row.names = FALSE)
cat("source: CYC2008 complexes", length(gen_cy), "; >=3 genes in universe", length(B$s3), "; merged", length(B$main), "\n")
print(B$counts, row.names = FALSE); print(t(ov)); print(Cp$counts[Cp$counts$variant == "main", ], row.names = FALSE)
