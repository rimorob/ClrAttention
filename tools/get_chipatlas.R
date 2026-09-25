#!/usr/bin/env Rscript
# Graded ChIP-seq binding for the BEELINE TFs, from ChIP-Atlas "target genes"
# tables (https://chip-atlas.dbcls.jp/data/<genome>/target/<antigen>.<kb>.tsv).
# Each table gives, for every gene whose TSS lies within +/- kb of a peak, the
# MACS2 score of the strongest such peak in each ChIP-seq experiment, and the
# average over all experiments of that antigen. Columns are "SRX...|cell type".
#
# For every TF that BEELINE selects (TFs + 1000 in any dataset of the species)
# and that ChIP-Atlas covers, this writes data/chipatlas/<genome>/<TF>.tsv with
#   gene                 upper-case symbol
#   all                  ChIP-Atlas average over all experiments
#   <dataset>            mean over experiments in cell types matching that
#                        BEELINE dataset (NA column if there are none)
# and data/chipatlas/<genome>/matched_experiments.csv (TF, dataset, number of
# matched experiments, their cell-type labels). Genes absent from a table have
# no peak in the window (score 0).
#
# Cell-type matching (fixed in advance, case-insensitive, on the label):
#   hESC     hESC_H1 / hESC_H9 / hESC_HUESxx / ES_cells / H1 / H9
#   hHep     Hep_G2 / Hepatocytes / Liver
#   mESC     ES_cells / Embryonic_stem_cells
#   mDC      labels containing "dendritic", excluding progenitors
#   mHSC-*   haematopoietic stem / progenitor labels, HPC-7, LSK, multipotent
#            progenitors (one set for mHSC-E, -GM and -L, as in BEELINE)
#
# Usage (repo root):  Rscript tools/get_chipatlas.R [--kb 1] [--workers 8]
args <- commandArgs(trailingOnly = TRUE)
opt <- list(kb = "1", workers = "8", beeline = "data/beeline", out = "data/chipatlas")
i <- 1L
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
root <- opt$beeline
expr_dir <- file.path(root, "BEELINE-data", "inputs", "scRNA-Seq")
tf_file <- function(sp) {
  cand <- c(file.path(root, paste0(sp, "-tfs.csv")), file.path(root, "Networks", paste0(sp, "-tfs.csv")))
  cand[file.exists(cand)][1]
}
species <- list(hg38 = list(sp = "human", ds = c("hESC", "hHep")),
                mm10 = list(sp = "mouse", ds = c("mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L")))
match_ct <- list(
  hESC = function(l) grepl("^(hESC_(H1|H9|HUES[0-9]+)|ES_cells|H1|H9)$", l, ignore.case = TRUE),
  hHep = function(l) grepl("^(Hep_G2|HepG2|Hepatocytes|Liver)$", l, ignore.case = TRUE),
  mESC = function(l) grepl("^(ES_cells|Embryonic_stem_cells|mESCs?)$", l, ignore.case = TRUE),
  mDC = function(l) grepl("dendritic", l, ignore.case = TRUE) & !grepl("progenitor", l, ignore.case = TRUE),
  mHSC = function(l) grepl("haematopoietic|hematopoietic|HPC-7|LSK|multipotent", l, ignore.case = TRUE))
ds_key <- function(d) if (grepl("^mHSC", d)) "mHSC" else d

al <- utils::read.delim("https://chip-atlas.dbcls.jp/data/metadata/analysisList.tab",
                        header = FALSE, stringsAsFactors = FALSE, quote = "")
names(al)[1:4] <- c("antigen", "coloc", "target", "genome")
for (g in names(species)) {
  sp <- species[[g]]
  tfs <- toupper(utils::read.csv(tf_file(sp$sp), stringsAsFactors = FALSE)$TF)
  sel <- unique(unlist(lapply(sp$ds, function(d) {
    go <- utils::read.csv(file.path(expr_dir, d, "GeneOrdering.csv"), stringsAsFactors = FALSE,
                          check.names = FALSE)
    names(go)[1] <- "gene"
    sig <- go[go$VGAMpValue * nrow(go) < 0.01, ]
    toupper(sig$gene[toupper(sig$gene) %in% tfs])
  })))
  av <- al[al$genome == g & al$target == "+", ]
  av <- av[toupper(av$antigen) %in% sel, ]
  keys <- unique(vapply(sp$ds, ds_key, ""))
  od <- file.path(opt$out, g); dir.create(od, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("%s: %d BEELINE-selected TFs, %d with ChIP-Atlas target tables", g,
                  length(sel), nrow(av)))
  one <- function(a) {
    f <- file.path(od, paste0(toupper(a), ".tsv"))
    if (file.exists(f)) return(NULL)
    url <- sprintf("https://chip-atlas.dbcls.jp/data/%s/target/%s.%s.tsv", g, a, opt$kb)
    x <- tryCatch(utils::read.delim(url, check.names = FALSE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(x) || ncol(x) < 3) return(data.frame(tf = toupper(a), dataset = NA, n = 0, labels = "download failed"))
    lab <- sub("^[^|]*\\|", "", names(x)[-(1:2)])
    M <- as.matrix(x[, -(1:2), drop = FALSE])
    out <- data.frame(gene = toupper(x[[1]]), all = x[[2]])
    info <- list()
    for (k in keys) {
      m <- match_ct[[k]](lab)
      out[[k]] <- if (any(m)) rowMeans(M[, m, drop = FALSE]) else NA_real_
      info[[k]] <- data.frame(tf = toupper(a), dataset = k, n = sum(m),
                              labels = paste(unique(lab[m]), collapse = ";"))
    }
    out <- out[!duplicated(out$gene), ]
    utils::write.table(out, paste0(f, ".part"), sep = "\t", row.names = FALSE, quote = FALSE)
    file.rename(paste0(f, ".part"), f)
    do.call(rbind, info)
  }
  res <- parallel::mclapply(av$antigen, one, mc.cores = as.integer(opt$workers),
                            mc.preschedule = FALSE)
  info <- do.call(rbind, res[!vapply(res, is.null, NA)])
  mf <- file.path(od, "matched_experiments.csv")
  if (!is.null(info)) {
    old <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else NULL
    utils::write.csv(rbind(old, info), mf, row.names = FALSE)
  }
  info <- utils::read.csv(mf, stringsAsFactors = FALSE)
  for (k in keys) message(sprintf("   %-5s TFs with >= 1 matched experiment: %d", k,
                                  sum(info$dataset == k & info$n > 0, na.rm = TRUE)))
}
