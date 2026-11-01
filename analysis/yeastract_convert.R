#!/usr/bin/env Rscript
# YEASTRACT export -> normalized association table for analysis/yeast_prep.R regulon_structure
# (docs/preregistration_reciprocity_yeast.md, "Yeast regulons", review2 SF1 iii). DRAFT rev 4.
# Committed BEFORE any YEASTRACT export is performed. If a real export cannot be parsed by this committed
# version, any fix is a logged deviation, made and committed before any association row is used.
#
# Output columns: tf, target, evidence (binding|expression), pmids (PMIDs YEASTRACT attributes to THIS
# evidence row; ";"-separated), pmids_unattributed (PMIDs given only as one combined list per association),
# value (as exported; matrix entries kept, e.g. 1 / -1).
#
# Input kinds (choose by flags; never by inspecting association content):
#   --table <file>           a delimited association table with a header. Columns are detected by fixed regexes:
#                            tf  ~ ^(tf|regulator|transcription.?factor)  target ~ ^(target|gene|orf)
#                            evidence ~ evidence    refs ~ (reference|pubmed|pmid)    value ~ (regulation|type|value|sign)
#                            Evidence labels map by the fixed table EVID below. A row whose label maps to BOTH
#                            classes yields one binding and one expression row; its references are then a combined
#                            list and go to pmids_unattributed (so the computed mode is binding_only for them).
#   --binding_matrix <file> --expression_matrix <file>
#                            YEASTRACT "Generate Regulation Matrix" exports (";"-separated, first column = gene,
#                            other columns = TFs, non-zero = association). No PMIDs exist -> mode binding_only.
#   --out <file>             normalized TSV.
args <- commandArgs(trailingOnly = TRUE); opt <- list()
i <- 1L; while (i < length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
if (is.null(opt$out)) stop("--out required")
EVID <- list(   # fixed label table (case-insensitive regex -> classes), checked in this order
  list(re = "dna.?binding.*(and|plus|&).*expression|expression.*(and|plus|&).*dna.?binding", cls = c("binding", "expression")),
  list(re = "dna.?binding|direct|chip|emsa|footprint", cls = "binding"),
  list(re = "expression|indirect|microarray|rna.?seq", cls = "expression"))
map_ev <- function(lbl) { for (e in EVID) if (grepl(e$re, lbl, ignore.case = TRUE)) return(e$cls); character(0) }
pmids_of <- function(x) { v <- regmatches(x, gregexpr("[0-9]{6,9}", x)); vapply(v, function(z) paste(unique(z), collapse = ";"), "") }
out <- NULL
if (!is.null(opt$table)) {
  first <- readLines(opt$table, 1L); sep <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else ","
  x <- utils::read.delim(opt$table, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", quote = "\"")
  pick <- function(re) { k <- grep(re, names(x), ignore.case = TRUE); if (!length(k)) stop("column not found: ", re); k[1] }
  ctf <- pick("^(tf|regulator|transcription.?factor)"); ctg <- pick("^(target|gene|orf)")
  cev <- pick("evidence"); crf <- pick("(reference|pubmed|pmid)")
  cvl <- grep("(regulation|type|value|sign)", names(x), ignore.case = TRUE); cvl <- if (length(cvl)) cvl[1] else NA
  rows <- lapply(seq_len(nrow(x)), function(k) { cls <- map_ev(x[[cev]][k]); if (!length(cls)) return(NULL)
    p <- pmids_of(x[[crf]][k]); both <- length(cls) > 1
    data.frame(tf = x[[ctf]][k], target = x[[ctg]][k], evidence = cls, pmids = if (both) "" else p,
               pmids_unattributed = if (both) p else "", value = if (is.na(cvl)) "1" else x[[cvl]][k]) })
  out <- do.call(rbind, rows)
  unmapped <- sort(unique(x[[cev]][vapply(x[[cev]], function(l) !length(map_ev(l)), TRUE)]))
  writeLines(unmapped, paste0(opt$out, ".unmapped_evidence_labels.txt"))
}
mat_rows <- function(f, cls) { m <- utils::read.csv(f, sep = ";", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  v <- as.matrix(m[, -1]); w <- which(v != "0" & nzchar(v) & !is.na(v), arr.ind = TRUE)
  data.frame(tf = colnames(v)[w[, 2]], target = m[[1]][w[, 1]], evidence = cls, pmids = "", pmids_unattributed = "", value = v[w]) }
if (!is.null(opt$binding_matrix)) out <- rbind(out, mat_rows(opt$binding_matrix, "binding"))
if (!is.null(opt$expression_matrix)) out <- rbind(out, mat_rows(opt$expression_matrix, "expression"))
if (is.null(out)) stop("no input given")
utils::write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d rows (%d binding, %d expression) -> %s\n", nrow(out), sum(out$evidence == "binding"), sum(out$evidence == "expression"), opt$out))
