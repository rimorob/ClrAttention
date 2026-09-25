#!/usr/bin/env Rscript
# Evidence-filtered STRING truths for BEELINE (decided 2026-09-24). BEELINE's
# STRING network mixes every evidence channel, including co-expression (which
# is correlation by another name, so partly circular) and text mining (which
# does not separate direct from indirect relations). Two derived truths, both
# written as TF -> gene edge lists (Gene1 = TF) in BEELINE's format to
# data/beeline/Networks_derived/<species>/:
#
#   STRING-regulatory.csv  STRING v11.0 "actions" with mode = expression
#                          (transcriptional regulation), directional, TF acting
#                          on target, score >= 400. v11.0 is the last release
#                          that published actions. Directed TF -> target.
#   STRING-curated.csv     STRING v12.0 channel scores: edges with
#                          experimental >= 400 or database >= 400 (medium
#                          confidence), i.e. excluding edges supported only by
#                          co-expression, text mining or genomic context; kept
#                          where one end is a TF, oriented TF -> other gene.
#                          Undirected functional association with the
#                          circular channels removed; still not regulation.
#
# Inputs: data/string (tools/get_string.sh), data/beeline (tools/get_beeline.sh).
# Usage (repo root):  Rscript analysis/beeline_truths.R [--min_score 400]
args <- commandArgs(trailingOnly = TRUE)
min_score <- if (length(args) >= 2 && args[1] == "--min_score") as.integer(args[2]) else 400L
sdir <- "data/string"; bdir <- "data/beeline"
tf_file <- function(sp) {
  cand <- c(file.path(bdir, paste0(sp, "-tfs.csv")), file.path(bdir, "Networks", paste0(sp, "-tfs.csv")))
  cand[file.exists(cand)][1]
}
rd <- function(f, ...) utils::read.delim(gzfile(file.path(sdir, f)), stringsAsFactors = FALSE, ...)
for (sp in c(human = "9606", mouse = "10090")) {
  spn <- names(which(c(human = "9606", mouse = "10090") == sp))
  tfs <- toupper(utils::read.csv(tf_file(spn), stringsAsFactors = FALSE)$TF)
  od <- file.path(bdir, "Networks_derived", spn); dir.create(od, recursive = TRUE, showWarnings = FALSE)

  info11 <- rd(sprintf("%s.protein.info.v11.0.txt.gz", sp), quote = "", comment.char = "")
  nm11 <- stats::setNames(toupper(info11[[2]]), info11[[1]])
  act <- rd(sprintf("%s.protein.actions.v11.0.txt.gz", sp), quote = "")
  act <- act[act$mode == "expression" & act$is_directional == "t" & act$a_is_acting == "t" &
               act$score >= min_score, ]
  reg <- unique(data.frame(Gene1 = nm11[act$item_id_a], Gene2 = nm11[act$item_id_b]))
  reg <- reg[!is.na(reg$Gene1) & !is.na(reg$Gene2) & reg$Gene1 != reg$Gene2 & reg$Gene1 %in% tfs, ]
  utils::write.csv(reg, file.path(od, "STRING-regulatory.csv"), row.names = FALSE, quote = FALSE)

  info12 <- rd(sprintf("%s.protein.info.v12.0.txt.gz", sp), quote = "", comment.char = "")
  nm12 <- stats::setNames(toupper(info12[[2]]), info12[[1]])
  L <- utils::read.table(gzfile(file.path(sdir, sprintf("%s.protein.links.detailed.v12.0.txt.gz", sp))),
                         header = TRUE, stringsAsFactors = FALSE)
  L <- L[L$experimental >= min_score | L$database >= min_score, c("protein1", "protein2")]
  a <- nm12[L$protein1]; b <- nm12[L$protein2]
  cur <- unique(rbind(data.frame(Gene1 = a[a %in% tfs], Gene2 = b[a %in% tfs]),
                      data.frame(Gene1 = b[b %in% tfs], Gene2 = a[b %in% tfs])))
  cur <- cur[!is.na(cur$Gene2) & cur$Gene1 != cur$Gene2, ]
  utils::write.csv(cur, file.path(od, "STRING-curated.csv"), row.names = FALSE, quote = FALSE)
  message(sprintf("%s: STRING-regulatory %d edges (%d TFs); STRING-curated %d edges (%d TFs)",
                  spn, nrow(reg), length(unique(reg$Gene1)), nrow(cur), length(unique(cur$Gene1))))
  rm(L, act); invisible(gc())
}
