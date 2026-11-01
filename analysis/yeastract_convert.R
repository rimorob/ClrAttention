#!/usr/bin/env Rscript
# YEASTRACT export -> normalized association table for analysis/yeast_prep.R regulon_structure
# (docs/preregistration_reciprocity_yeast.md, "Yeast regulons", review2 SF1 iii, review3). DRAFT rev 5.
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
#                            evidence ~ evidence    refs ~ (reference|pubmed|pmid)
#                            value ~ (regulation|sign|value), never a column whose name contains "evidence"
#                            An evidence label is split into tokens at ; , / | and each token must match exactly one
#                            ANCHORED entry of the fixed table EVID (case-insensitive, whole token). The row's classes are
#                            the union of its tokens' classes ("Expression; ChIP-chip" -> binding + expression). A row
#                            with both classes yields one binding and one expression row whose references form a combined
#                            list (pmids_unattributed), so the computed mode is binding_only for them.
#                            ANY unrecognised token STOPS the converter (no silent drops); fixing it is a logged deviation.
#   --binding_matrix <file> --expression_matrix <file>
#                            YEASTRACT "Generate Regulation Matrix" exports (";"-separated, first column = gene,
#                            other columns = TFs, non-zero = association). No PMIDs exist -> mode binding_only.
#   --out <file>             normalized TSV.
args <- commandArgs(trailingOnly = TRUE); opt <- list()
i <- 1L; while (i < length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
if (is.null(opt$out)) stop("--out required")
EVID <- list(   # fixed, anchored label table: whole token (trimmed, lower case) -> classes
  list(re = "^dna[ -]?binding (and|plus|&) expression( evidence)?$", cls = c("binding", "expression")),
  list(re = "^(dna[ -]?binding( evidence)?|direct( evidence)?|chip(-chip|-seq|-exo|-on-chip)?|emsa|dnase i footprinting|footprinting)$", cls = "binding"),
  list(re = "^(expression( evidence)?|indirect( evidence)?|microarray|rna-?seq)$", cls = "expression"))
map_token <- function(tok) { hits <- Filter(function(e) grepl(e$re, tok), EVID)
  if (length(hits) != 1L) return(NULL); hits[[1]]$cls }
map_ev <- function(lbl) {     # NULL if any token is unrecognised or ambiguous
  toks <- trimws(tolower(strsplit(lbl, "[;,/|]")[[1]])); toks <- toks[nzchar(toks)]
  if (!length(toks)) return(NULL); cl <- lapply(toks, map_token)
  if (any(vapply(cl, is.null, TRUE))) return(NULL); sort(unique(unlist(cl))) }
pmids_of <- function(x) { v <- regmatches(x, gregexpr("[0-9]{6,9}", x)); vapply(v, function(z) paste(unique(z), collapse = ";"), "") }
out <- NULL
if (!is.null(opt$table)) {
  first <- readLines(opt$table, 1L); sep <- if (grepl("\t", first)) "\t" else if (grepl(";", first)) ";" else ","
  x <- utils::read.delim(opt$table, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", quote = "\"")
  pick <- function(re) { k <- grep(re, names(x), ignore.case = TRUE); if (!length(k)) stop("column not found: ", re); k[1] }
  ctf <- pick("^(tf|regulator|transcription.?factor)"); ctg <- pick("^(target|gene|orf)")
  cev <- pick("evidence"); crf <- pick("(reference|pubmed|pmid)")
  cvl <- setdiff(grep("(regulation|sign|value)", names(x), ignore.case = TRUE), grep("evidence", names(x), ignore.case = TRUE))
  cvl <- if (length(cvl)) cvl[1] else NA
  bad <- sort(unique(x[[cev]][vapply(x[[cev]], function(l) is.null(map_ev(l)), TRUE)]))
  if (length(bad)) { writeLines(bad, paste0(opt$out, ".unrecognised_evidence_labels.txt"))
    stop("unrecognised evidence label(s): ", paste(utils::head(bad, 10), collapse = " | "), " -> STOP (log a converter deviation)") }
  rows <- lapply(seq_len(nrow(x)), function(k) { cls <- map_ev(x[[cev]][k])
    p <- pmids_of(x[[crf]][k]); both <- length(cls) > 1
    data.frame(tf = x[[ctf]][k], target = x[[ctg]][k], evidence = cls, pmids = if (both) "" else p,
               pmids_unattributed = if (both) p else "", value = if (is.na(cvl)) "1" else x[[cvl]][k]) })
  out <- do.call(rbind, rows)
}
mat_rows <- function(f, cls) { m <- utils::read.csv(f, sep = ";", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  v <- as.matrix(m[, -1]); w <- which(v != "0" & nzchar(v) & !is.na(v), arr.ind = TRUE)
  data.frame(tf = colnames(v)[w[, 2]], target = m[[1]][w[, 1]], evidence = cls, pmids = "", pmids_unattributed = "", value = v[w]) }
if (!is.null(opt$binding_matrix)) out <- rbind(out, mat_rows(opt$binding_matrix, "binding"))
if (!is.null(opt$expression_matrix)) out <- rbind(out, mat_rows(opt$expression_matrix, "expression"))
if (is.null(out)) stop("no input given")
utils::write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d rows (%d binding, %d expression) -> %s\n", nrow(out), sum(out$evidence == "binding"), sum(out$evidence == "expression"), opt$out))
