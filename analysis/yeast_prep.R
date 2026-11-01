#!/usr/bin/env Rscript
# Data acquisition and label-free preparation for docs/preregistration_reciprocity_yeast.md (DRAFT rev 2).
# Written 2026-11-01; parse-checked; NOT yet run on real data (nothing downloaded).
#
# Stages, in protocol order (outputs in data/yeast/ [not committed; hashes are] or
# docs/reciprocity_yeast_selections/ [committed]):
#   download         expression + wild-type controls + annotation only (no truth). SHA-256 -> data/yeast/SHA256SUMS
#   matrix           700 responsive mutants (GSE42527 published M) -> ORF x mutant matrix; deleted-gene resolution
#                    via SGD names/aliases (unresolvable columns dropped); masking of the deleted ORF and ORFs
#                    overlapping / within 500 bp; divergent-pair list (label-free, for a descriptive sensitivity)
#   check            orientation rule + reproduction of the published M for the first 10 mutants
#   rebuild          only if check failed: sample tables of a fixed random 376 mutants (seed 20261105), rebuilt
#   wt               wild-type controls (GSE45115 sample tables) -> ORF x culture matrix
#   splits           10 random halves (seed 20261101) + label-free mutant clusters (complete linkage, 1 - r, h 0.5)
#   feasibility      gate F (sanity check): chance-corrected cross-half reproducibility R' (mutants vs WT) and
#                    overlap share Q of mutant top-1% pairs in the WT top-1%; pass iff median R'_mut > median R'_wt
#                    and Q <= 0.5
#   truth_structure  GATED: preregistration fixed (tracked, == HEAD, no "STATUS: DRAFT"); F passed and committed;
#                    splits/clusters, selections, R check and the test script committed (blob hash recorded);
#                    excluded_pmids.txt and the YEASTRACT export log committed. Builds member lists with the
#                    exclusion (E) and evidence filters, the contamination audit, evaluability counts, n_eff for
#                    P1/S1/P2, the CYC2008 fallback and the power-gate line. Reads membership only; never loads a
#                    network or expression values (only the universe row names).
# No stage scores anything. Stages other than truth_structure never download or open a truth.
# Usage: Rscript analysis/yeast_prep.R <stage> [--threads 4]
#        truth_structure: --yeastract <normalized TSV> --yeastract_mode refs|binding_only [--cyc2008 <file>]
#   normalized YEASTRACT TSV (produced from the logged export; columns): tf, target, evidence
#   (binding|expression), pmids (";"-separated, may be empty in binding_only mode), value (as exported).

args <- commandArgs(trailingOnly = TRUE)
STAGES <- c("download", "matrix", "check", "rebuild", "wt", "splits", "feasibility", "truth_structure")
if (!length(args) || !args[1] %in% STAGES) stop("usage: Rscript analysis/yeast_prep.R <", paste(STAGES, collapse = "|"), "> [options]")
stage <- args[1]
opt <- list(threads = 4L, yeastract = NA_character_, yeastract_mode = NA_character_, cyc2008 = NA_character_)
i <- 2L
while (i < length(args)) { k <- sub("^--", "", args[i]); v <- args[i + 1L]
  opt[[k]] <- if (k == "threads") as.integer(v) else v; i <- i + 2L }

PREREG <- "docs/preregistration_reciprocity_yeast.md"
TEST_SCRIPT <- "analysis/reciprocity_yeast_test.R"
D <- "data/yeast"; SEL <- "docs/reciprocity_yeast_selections"
GEO <- "https://ftp.ncbi.nlm.nih.gov/geo"
ACC <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sha <- function(f) digest::digest(file = f, algo = "sha256")
mkd <- function() { dir.create(D, showWarnings = FALSE, recursive = TRUE); dir.create(SEL, showWarnings = FALSE, recursive = TRUE) }
fetch <- function(url, dest, tries = 4L) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(invisible(dest))
  for (k in seq_len(tries)) {
    ok <- tryCatch({ utils::download.file(url, dest, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (ok && file.exists(dest) && file.info(dest)$size > 0) return(invisible(dest))
    Sys.sleep(2^k)
  }
  stop("download failed: ", url)
}
record_hashes <- function() {
  fs <- list.files(D, recursive = TRUE, full.names = TRUE)
  fs <- fs[basename(fs) != "SHA256SUMS" & !grepl("/truth_raw/", fs)]
  writeLines(sprintf("%s  %s", vapply(fs, sha, ""), sub(paste0("^", D, "/"), "", fs)), file.path(D, "SHA256SUMS"))
  say("hashes: ", length(fs), " files -> ", file.path(D, "SHA256SUMS"))
}
tracked <- function(f) system2("git", c("ls-files", "--error-unmatch", f), stdout = FALSE, stderr = FALSE) == 0L
ORF_RE <- "^(Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4}|R[0-9]{4}[WC]?)$"
# replicate grouping: mutants by deleted gene (the two cultures may be different isolates, e.g.
# rad18-del-4-b with rad18-del-3-a), wild types by dye-swap pair (wt-SV013-a / -b)
group_key <- function(title) ifelse(grepl("-del-", title), sub("-del-.*$", "", title), sub("-[ab]$", "", title))
rd_mat <- function(f) { x <- utils::read.delim(f, check.names = FALSE); m <- as.matrix(x[, -1]); rownames(m) <- x$orf; m }
wr_mat <- function(m, f) utils::write.table(data.frame(orf = rownames(m), m, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = FALSE)

## ---- annotation helpers -------------------------------------------------------------------------
sgd <- function() {      # SGD_features.tab: 4 ORF name, 5 standard name, 6 aliases (|), 9 chromosome, 10/11 coords, 12 strand
  x <- utils::read.delim(file.path(D, "SGD_features.tab"), header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  x <- x[grepl(ORF_RE, x$V4), ]
  data.frame(orf = x$V4, name = toupper(x$V5), alias = toupper(x$V6), chr = x$V9,
             lo = pmin(as.numeric(x$V10), as.numeric(x$V11)), hi = pmax(as.numeric(x$V10), as.numeric(x$V11)),
             strand = x$V12, stringsAsFactors = FALSE)
}
name_to_orf <- function(s) {   # standard name, then systematic, then alias (only if the alias is unique)
  m <- stats::setNames(s$orf, s$name); m <- m[nzchar(names(m))]
  al <- do.call(rbind, lapply(seq_len(nrow(s)), function(k) if (nzchar(s$alias[k]))
    data.frame(a = strsplit(s$alias[k], "|", fixed = TRUE)[[1]], orf = s$orf[k]) else NULL))
  al <- al[!al$a %in% al$a[duplicated(al$a)], ]
  function(x) { x <- toupper(x); out <- unname(m[x]); out[is.na(out) & x %in% s$orf] <- x[is.na(out) & x %in% s$orf]
    miss <- is.na(out); out[miss] <- al$orf[match(x[miss], al$a)]; out }
}
neighbours <- function(s, orf, win = 500) {   # ORFs overlapping or within `win` bp of `orf` (same chromosome)
  k <- match(orf, s$orf); if (is.na(k)) return(character(0))
  same <- s$chr == s$chr[k]; gap <- pmax(0, pmax(s$lo, s$lo[k]) - pmin(s$hi, s$hi[k]))
  s$orf[same & gap <= win]
}
divergent_pairs <- function(s) {   # adjacent genes on opposite strands with 5' ends facing (shared intergenic region)
  out <- list()
  for (ch in unique(s$chr)) { x <- s[s$chr == ch, ]; x <- x[order(x$lo), ]
    if (nrow(x) < 2) next
    for (k in seq_len(nrow(x) - 1L)) if (x$strand[k] == "C" && x$strand[k + 1L] == "W")
      out[[length(out) + 1L]] <- data.frame(a = x$orf[k], b = x$orf[k + 1L]) }
  do.call(rbind, out)
}

## ---- GEO helpers ------------------------------------------------------------------------------
series_meta <- function(gse) {
  f <- file.path(D, paste0(gse, "_series_matrix.txt.gz"))
  L <- readLines(gzfile(f)); L <- L[startsWith(L, "!Sample_")]
  get <- function(key) { x <- L[startsWith(L, paste0(key, "\t"))]; if (!length(x)) return(NULL)
    gsub('^"|"$', "", strsplit(x[1], "\t")[[1]][-1]) }
  data.frame(title = get("!Sample_title"), gsm = get("!Sample_geo_accession"),
             src1 = get("!Sample_source_name_ch1"), src2 = get("!Sample_source_name_ch2"), stringsAsFactors = FALSE)
}
gsm_table <- function(gsm) {
  f <- file.path(D, "gsm", paste0(gsm, ".txt"))
  if (!file.exists(f)) { dir.create(dirname(f), showWarnings = FALSE)
    fetch(sprintf("%s?acc=%s&targ=self&form=text&view=data", ACC, gsm), f); Sys.sleep(0.4) }
  L <- readLines(f); a <- which(L == "!sample_table_begin"); b <- which(L == "!sample_table_end")
  if (!length(a)) stop("no sample table in ", f)
  tb <- utils::read.delim(text = L[(a + 1L):(if (length(b)) b - 1L else length(L))], stringsAsFactors = FALSE)
  stats::setNames(as.numeric(tb$VALUE), as.character(tb$ID_REF))
}
platform <- function() {
  f <- file.path(D, "GPL11232_full.txt"); L <- readLines(f)
  a <- which(L == "!platform_table_begin"); b <- which(L == "!platform_table_end")
  tb <- utils::read.delim(text = L[(a + 1L):(b - 1L)], stringsAsFactors = FALSE, quote = "", comment.char = "")
  frac <- function(re) vapply(tb, function(x) mean(grepl(re, as.character(x))), 0)
  scab <- names(which.max(frac("^SCAB[0-9]{6}$"))); orf <- names(which.max(frac(ORF_RE)))
  if (max(frac("^SCAB[0-9]{6}$")) < 0.5 || max(frac(ORF_RE)) < 0.3)
    stop("GPL11232: could not detect SCAB / ORF columns; inspect ", f, " and set them explicitly (log the choice)")
  data.frame(id = as.character(tb$ID), scab = tb[[scab]], orf = ifelse(grepl(ORF_RE, tb[[orf]]), tb[[orf]], NA),
             stringsAsFactors = FALSE)
}
by_orf <- function(M, orf) {
  keep <- !is.na(orf); M <- M[keep, , drop = FALSE]; orf <- orf[keep]
  S <- rowsum(ifelse(is.na(M), 0, M), orf); N <- rowsum(1 * !is.na(M), orf); S / ifelse(N == 0, NA, N)
}
read_final <- function(g) {
  x <- utils::read.delim(gzfile(file.path(D, paste0(g, "_final_GeneExpressionMatrix.txt.gz"))), check.names = FALSE,
                         stringsAsFactors = FALSE)
  mcol <- grep("^M \\(", names(x)); pcol <- grep("^p-value \\(", names(x))
  M <- as.matrix(x[, mcol]); P <- as.matrix(x[, pcol]); rownames(M) <- rownames(P) <- x$PROBE_ID
  lab <- sub("^M \\((.*?)\\).*$", "\\1", names(x)[mcol])
  colnames(M) <- colnames(P) <- sub("-del-.*$", "", lab)
  list(M = M, P = P, label = names(x)[mcol])
}
rebuild <- function(meta, rule, pf) {   # mean of the oriented hybridizations per group, then per ORF
  key <- group_key(meta$title)
  vals <- lapply(seq_len(nrow(meta)), function(k) { v <- gsm_table(meta$gsm[k])
    if (rule == "flip_ch1" && !grepl("ref", meta$src1[k], ignore.case = TRUE)) v <- -v; v })
  ids <- Reduce(union, lapply(vals, names))
  V <- vapply(vals, function(v) v[ids], numeric(length(ids))); rownames(V) <- ids
  V <- vapply(split(seq_len(ncol(V)), key), function(ix) rowMeans(V[, ix, drop = FALSE], na.rm = TRUE), numeric(nrow(V)))
  rownames(V) <- ids; by_orf(V, pf$orf[match(rownames(V), pf$id)])
}
finish_matrix <- function(X, mask) {   # missing filter (masked entries not counted), median imputation, constant rows
  miss <- rowMeans(is.na(X) & !mask); keep <- miss <= 0.05; X <- X[keep, , drop = FALSE]; mask <- mask[keep, , drop = FALSE]
  X[mask] <- NA; med <- apply(X, 1, stats::median, na.rm = TRUE); ix <- which(is.na(X), arr.ind = TRUE); X[ix] <- med[ix[, 1]]
  X <- X[apply(X, 1, function(v) diff(range(v))) > 0, , drop = FALSE]
  list(X = X, dropped_missing = sum(!keep))
}

## ---- stage: download ----------------------------------------------------------------------------
if (stage == "download") {
  mkd()
  for (g in c("GSE42527", "GSE42526")) {
    fetch(sprintf("%s/series/GSE42nnn/%s/suppl/%s_final_GeneExpressionMatrix.txt.gz", GEO, g, g),
          file.path(D, paste0(g, "_final_GeneExpressionMatrix.txt.gz")))
    fetch(sprintf("%s/series/GSE42nnn/%s/suppl/%s_protocols.txt.gz", GEO, g, g), file.path(D, paste0(g, "_protocols.txt.gz")))
    fetch(sprintf("%s/series/GSE42nnn/%s/matrix/%s_series_matrix.txt.gz", GEO, g, g), file.path(D, paste0(g, "_series_matrix.txt.gz")))
  }
  fetch(sprintf("%s/series/GSE45nnn/GSE45115/matrix/GSE45115_series_matrix.txt.gz", GEO), file.path(D, "GSE45115_series_matrix.txt.gz"))
  fetch(sprintf("%s?acc=GPL11232&targ=self&form=text&view=data", ACC), file.path(D, "GPL11232_full.txt"))
  fetch("https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab", file.path(D, "SGD_features.tab"))
  tryCatch(fetch("http://deleteome.holstegelab.nl/data/downloads/deleteome_all_mutants_svd_transformed.txt",
                 file.path(D, "deleteome_all_mutants_svd_transformed.txt")),   # descriptive sensitivity only; URL to verify
           error = function(e) say("growth-corrected matrix not fetched (descriptive only): ", conditionMessage(e)))
  wt <- series_meta("GSE45115"); say("GSE45115: ", nrow(wt), " wild-type hybridizations")
  for (g in wt$gsm) invisible(gsm_table(g))
  mu <- series_meta("GSE42527"); key <- group_key(mu$title)
  first10 <- sort(unique(key))[1:10]
  for (g in mu$gsm[key %in% first10]) invisible(gsm_table(g))
  record_hashes(); say("download done; commit data/yeast/SHA256SUMS (data files themselves are not committed)")
}

## ---- stage: matrix ------------------------------------------------------------------------------
if (stage == "matrix") {
  mkd(); pf <- platform(); s <- sgd(); to_orf <- name_to_orf(s)
  build <- function(g) {
    r <- read_final(g); X <- by_orf(r$M, pf$orf[match(rownames(r$M), pf$scab)])
    del <- to_orf(colnames(X)); ok <- !is.na(del) & del %in% rownames(X)
    dropped <- data.frame(series = g, column = r$label[!ok], deleted_gene = colnames(X)[!ok])
    X <- X[, ok, drop = FALSE]; del <- del[ok]
    mask <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X)); nmask <- integer(ncol(X))
    for (j in seq_len(ncol(X))) { nb <- intersect(neighbours(s, del[j], 500), rownames(X)); mask[nb, j] <- TRUE; nmask[j] <- length(nb) }
    f <- finish_matrix(X, mask)
    list(X = f$X, dropped = dropped, dropped_missing = f$dropped_missing, nmask = nmask, del = del)
  }
  b <- build("GSE42527")
  wr_mat(b$X, file.path(D, "deleteome_responsive_700.tsv"))
  utils::write.csv(b$dropped, file.path(SEL, "dropped_columns.csv"), row.names = FALSE)
  utils::write.csv(data.frame(mutant = colnames(b$X), deleted_orf = b$del, n_masked = b$nmask), file.path(SEL, "masking_log.csv"), row.names = FALSE)
  say(sprintf("responsive matrix: %d ORFs x %d mutants; %d columns dropped (unresolvable deleted gene); %d ORFs dropped (>5%% missing); masked per column median %d",
              nrow(b$X), ncol(b$X), nrow(b$dropped), b$dropped_missing, stats::median(b$nmask)))
  a <- build("GSE42526"); common <- intersect(rownames(b$X), rownames(a$X))           # descriptive sensitivity set
  wr_mat(cbind(b$X[common, ], a$X[common, ]), file.path(D, "deleteome_all_1484.tsv"))
  utils::write.csv(divergent_pairs(s), file.path(SEL, "divergent_pairs.csv"), row.names = FALSE)   # label-free annotation
  writeLines(sprintf("%s  %s", sha(file.path(D, "deleteome_responsive_700.tsv")), "deleteome_responsive_700.tsv"),
             file.path(SEL, "deleteome_matrix.sha256"))
  say("commit ", SEL, "/{deleteome_matrix.sha256, dropped_columns.csv, masking_log.csv, divergent_pairs.csv}")
}

## ---- stage: check -------------------------------------------------------------------------------
# Candidate rules (fixed in advance): "asis" = VALUE as deposited for both hybridizations;
# "flip_ch1" = negate VALUE of the hybridization whose channel-1 source is the test strain.
if (stage == "check") {
  pf <- platform(); mu <- series_meta("GSE42527"); key <- group_key(mu$title)
  first10 <- sort(unique(key))[1:10]; meta <- mu[key %in% first10, ]
  pub <- read_final("GSE42527"); Mo <- by_orf(pub$M, pf$orf[match(rownames(pub$M), pf$scab)])
  res <- do.call(rbind, lapply(c("asis", "flip_ch1"), function(rule) { R <- rebuild(meta, rule, pf)
    do.call(rbind, lapply(colnames(R), function(s) { j <- which(tolower(colnames(Mo)) == tolower(s))[1]
      g <- intersect(rownames(R), rownames(Mo)); a <- R[g, s]; b <- Mo[g, j]
      data.frame(rule = rule, mutant = s, r = stats::cor(a, b, use = "complete.obs"), med_abs = stats::median(abs(a - b), na.rm = TRUE)) })) }))
  pass <- tapply(res$r >= 0.99 & res$med_abs <= 0.02, res$rule, all)
  rule <- if (any(pass)) names(pass)[pass][1] else names(which.max(tapply(res$r, res$rule, stats::median)))
  res$chosen <- res$rule == rule; res$passed <- any(pass)
  utils::write.csv(res, file.path(SEL, "implementation_check.csv"), row.names = FALSE)
  writeLines(c(rule, if (any(pass)) "published_M_reproduced" else "NOT_reproduced_run_stage_rebuild"), file.path(SEL, "orientation_rule.txt"))
  say("orientation rule: ", rule, "; reproduction ", if (any(pass)) "PASSED" else "FAILED: run stage rebuild (gate F then uses rebuilt mutants)")
}

## ---- stage: rebuild (only after a failed check) -------------------------------------------------
if (stage == "rebuild") {
  rr <- readLines(file.path(SEL, "orientation_rule.txt"))
  if (rr[2] == "published_M_reproduced") stop("check passed; rebuild is not used")
  pf <- platform(); mu <- series_meta("GSE42527"); key <- group_key(mu$title)
  X <- rd_mat(file.path(D, "deleteome_responsive_700.tsv"))
  cand <- intersect(sort(unique(key)), tolower(colnames(X))); Nw <- length(unique(group_key(series_meta("GSE45115")$title)))
  set.seed(20261105); pick <- sort(sample(cand, min(Nw, length(cand))))
  R <- rebuild(mu[key %in% pick, ], rr[1], pf)
  m <- read.csv(file.path(SEL, "masking_log.csv")); s <- sgd()
  R <- R[intersect(rownames(X), rownames(R)), , drop = FALSE]
  mask <- matrix(FALSE, nrow(R), ncol(R), dimnames = dimnames(R))
  for (j in seq_len(ncol(R))) { d <- m$deleted_orf[match(tolower(colnames(R)[j]), tolower(m$mutant))]
    if (!is.na(d)) mask[intersect(neighbours(s, d, 500), rownames(R)), j] <- TRUE }
  f <- finish_matrix(R, mask); wr_mat(f$X, file.path(D, "mutants_rebuilt_376.tsv"))
  writeLines(pick, file.path(SEL, "rebuilt_mutants.txt"))
  say(sprintf("rebuilt %d mutants (rule %s) for gate F", ncol(f$X), rr[1]))
}

## ---- stage: wt ----------------------------------------------------------------------------------
if (stage == "wt") {
  rf <- file.path(SEL, "orientation_rule.txt"); if (!file.exists(rf)) stop("run stage check first")
  rule <- readLines(rf, 1L); pf <- platform(); W <- rebuild(series_meta("GSE45115"), rule, pf)
  X <- rd_mat(file.path(D, "deleteome_responsive_700.tsv"))
  W <- W[intersect(rownames(X), rownames(W)), , drop = FALSE]
  f <- finish_matrix(W, matrix(FALSE, nrow(W), ncol(W)))
  wr_mat(f$X, file.path(D, "wt_controls_GSE45115.tsv"))
  say(sprintf("wild-type matrix: %d ORFs x %d cultures (rule %s)", nrow(f$X), ncol(f$X), rule))
}

## ---- stage: splits (halves + label-free mutant clusters) ---------------------------------------
if (stage == "splits") {
  X <- rd_mat(file.path(D, "deleteome_responsive_700.tsv")); mut <- colnames(X); set.seed(20261101)
  sp <- vapply(1:10, function(s) { h <- integer(length(mut)); h[sample(length(mut), floor(length(mut) / 2))] <- 1L; h + 1L },
               integer(length(mut)))
  colnames(sp) <- sprintf("split_%02d", 1:10)
  utils::write.csv(data.frame(mutant = mut, sp), file.path(SEL, "splits.csv"), row.names = FALSE)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - stats::cor(X)), "complete"), h = 0.5)
  utils::write.csv(data.frame(mutant = mut, cluster = unname(cl)), file.path(SEL, "mutant_clusters.csv"), row.names = FALSE)
  sz <- table(cl); utils::write.csv(data.frame(n_clusters = length(sz), singletons = sum(sz == 1), largest = max(sz),
    t(stats::quantile(as.integer(sz), c(.5, .9, .99)))), file.path(SEL, "mutant_clusters_summary.csv"), row.names = FALSE)
  say("10 random halves of ", length(mut), " mutants; ", length(sz), " mutant clusters (largest ", max(sz), ") -> commit")
}

## ---- stage: feasibility (gate F, sanity check) --------------------------------------------------
if (stage == "feasibility") {
  suppressPackageStartupMessages(library(clr))
  rr <- readLines(file.path(SEL, "orientation_rule.txt"))
  Xm <- if (rr[2] == "published_M_reproduced") rd_mat(file.path(D, "deleteome_responsive_700.tsv")) else {
    f <- file.path(D, "mutants_rebuilt_376.tsv"); if (!file.exists(f)) stop("check failed: run stage rebuild first"); rd_mat(f) }
  W <- rd_mat(file.path(D, "wt_controls_GSE45115.tsv"))
  g <- intersect(rownames(Xm), rownames(W)); Xm <- Xm[g, ]; W <- W[g, ]
  fitC <- function(M) { rng <- apply(M, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) M[rng <= 0, ] <- M[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(M), sd = 1e-9)
    f <- ClrAttention$new(M)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
    f$calibrate(method = "normal", combine = "stouffer"); C <- f$clr_scores; C[lower.tri(C, diag = TRUE)] <- NA; C }
  top <- function(C) { v <- which(!is.na(C)); k <- ceiling(0.01 * length(v)); v[order(C[v], decreasing = TRUE)[seq_len(k)]] }
  repro <- function(M, cols) { h <- sample(cols, floor(length(cols) / 2)); a <- top(fitC(M[, h])); b <- top(fitC(M[, setdiff(cols, h)]))
    (length(intersect(a, b)) / length(a) - 0.01) / 0.99 }
  Nw <- min(ncol(W), ncol(Xm)); set.seed(20261199)
  Rm <- vapply(1:10, function(s) repro(Xm, sample(ncol(Xm), Nw)), 0)
  Rw <- vapply(1:10, function(s) repro(W, seq_len(ncol(W))), 0)
  tm <- top(fitC(Xm)); tw <- top(fitC(W)); Q <- length(intersect(tm, tw)) / length(tm)
  out <- data.frame(mutant_source = if (rr[2] == "published_M_reproduced") "published_M" else "rebuilt_376",
                    N_matched = Nw, genes = length(g), median_Rprime_mut = stats::median(Rm), median_Rprime_wt = stats::median(Rw),
                    iqr_Rprime_mut = stats::IQR(Rm), iqr_Rprime_wt = stats::IQR(Rw), Q_overlap_share = Q)
  out$pass <- out$median_Rprime_mut > out$median_Rprime_wt && Q <= 0.5
  utils::write.csv(out, file.path(SEL, "feasibility_gate.csv"), row.names = FALSE)
  utils::write.csv(data.frame(draw = 1:10, Rprime_mut = Rm, Rprime_wt = Rw), file.path(SEL, "feasibility_draws.csv"), row.names = FALSE)
  say("gate F: ", if (out$pass) "PASS" else "FAIL (verdict 1; yeast not consumed)",
      sprintf(" (median R' mut %.3f vs wt %.3f; Q %.3f)", out$median_Rprime_mut, out$median_Rprime_wt, Q))
}

## ---- stage: truth_structure (GATED) --------------------------------------------------------------
if (stage == "truth_structure") {
  ## gates
  if (!tracked(PREREG)) stop("gate: ", PREREG, " not tracked")
  if (system2("git", c("diff", "--quiet", "HEAD", "--", PREREG)) != 0L) stop("gate: ", PREREG, " differs from HEAD")
  if (any(grepl("STATUS: DRAFT", readLines(PREREG), fixed = TRUE))) stop("gate: ", PREREG, " is still a DRAFT; no truth may be downloaded")
  need <- file.path(SEL, c("excluded_pmids.txt", "feasibility_gate.csv", "splits.csv", "mutant_clusters.csv",
                           "yeast_selections.csv", "R_check.csv", "yeastract_export_log.md"))
  for (f in need) if (!tracked(f)) stop("gate: ", f, " not committed")
  for (f in need) if (system2("git", c("diff", "--quiet", "HEAD", "--", f)) != 0L) stop("gate: ", f, " differs from HEAD")
  if (!isTRUE(as.logical(utils::read.csv(file.path(SEL, "feasibility_gate.csv"))$pass))) stop("gate F failed: verdict 1, stop")
  if (!all(as.logical(utils::read.csv(file.path(SEL, "R_check.csv"))$passed))) stop("synthetic R check not passed: no verdict, fix first")
  if (!tracked(TEST_SCRIPT) || system2("git", c("diff", "--quiet", "HEAD", "--", TEST_SCRIPT)) != 0L)
    stop("gate: the complete test script ", TEST_SCRIPT, " must be committed before any truth is built")
  if (!opt$yeastract_mode %in% c("refs", "binding_only")) stop("--yeastract_mode refs|binding_only (as logged in yeastract_export_log.md)")
  if (is.na(opt$yeastract) || !file.exists(opt$yeastract)) stop("--yeastract <normalized TSV from the logged export>")
  hash <- system2("git", c("log", "-1", "--format=%H", "--", PREREG), stdout = TRUE)
  test_blob <- system2("git", c("rev-parse", paste0("HEAD:", TEST_SCRIPT)), stdout = TRUE)
  E <- readLines(file.path(SEL, "excluded_pmids.txt")); E <- E[grepl("^[0-9]+$", E)]
  universe <- utils::read.delim(file.path(D, "deleteome_responsive_700.tsv"), check.names = FALSE)$orf   # row names only
  s <- sgd(); to_orf <- name_to_orf(s)
  Tdir <- file.path(D, "truth_raw"); dir.create(Tdir, showWarnings = FALSE)
  neff_of <- function(sizes) { p <- choose(sizes, 2); p <- p[p > 0]; if (!length(p)) return(0); sum(p)^2 / sum(p^2) }

  ## regulons (YEASTRACT, normalized export)
  yx <- utils::read.delim(opt$yeastract, stringsAsFactors = FALSE, colClasses = "character")
  stopifnot(all(c("tf", "target", "evidence", "pmids", "value") %in% names(yx)))
  value_table <- as.data.frame(table(yx$value))
  yx <- yx[yx$value != "0" & nzchar(yx$value), ]                          # any non-zero entry (incl. -1) counts
  yx$tf_orf <- to_orf(yx$tf); yx$target_orf <- to_orf(yx$target)
  pm <- strsplit(yx$pmids, "[;, ]+")
  # reanalysis screen (c): PubMed title/abstract of every expression PMID not already in E
  if (opt$yeastract_mode == "refs") {
    cand <- setdiff(unique(unlist(pm[yx$evidence == "expression"])), c(E, "")); hit <- character(0)
    for (ch in split(cand, ceiling(seq_along(cand) / 150))) {
      f <- tempfile(); fetch(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&rettype=abstract&retmode=xml",
                                     paste(ch, collapse = ",")), f)
      x <- paste(readLines(f, warn = FALSE), collapse = " "); arts <- strsplit(x, "<PubmedArticle>", fixed = TRUE)[[1]][-1]
      for (a in arts) { id <- sub(".*?<PMID[^>]*>([0-9]+)</PMID>.*", "\\1", a)
        if (grepl("deleteome|GSE42526|GSE42527|GSE42528|GSE45115", a, ignore.case = TRUE)) hit <- c(hit, id) }
      Sys.sleep(0.4) }
    E2 <- unique(c(E, hit)); writeLines(sort(hit), file.path(SEL, "excluded_pmids_screen_additions.txt"))
  } else E2 <- E
  key <- paste(yx$tf_orf, yx$target_orf)
  bind <- unique(key[yx$evidence == "binding"])
  expr_any <- unique(key[yx$evidence == "expression"])
  expr_clean <- unique(key[yx$evidence == "expression" & vapply(pm, function(v) any(nzchar(v) & !v %in% E2), TRUE)])
  primary <- if (opt$yeastract_mode == "refs") intersect(bind, expr_clean) else bind
  circular <- intersect(bind, expr_any)
  regs_of <- function(keys) { p <- do.call(rbind, strsplit(keys, " ")); p <- p[!is.na(p[, 1]) & !is.na(p[, 2]) & p[, 1] != "NA" & p[, 2] != "NA", , drop = FALSE]
    r <- lapply(split(p[, 2], p[, 1]), function(v) intersect(unique(v), universe)); r[lengths(r) >= 5 & lengths(r) <= 500] }
  comem_pairs <- function(r) { if (!length(r)) return(0); g <- sort(unique(unlist(r)))
    M <- Matrix::sparseMatrix(i = match(unlist(r), g), j = rep(seq_along(r), lengths(r)), x = 1, dims = c(length(g), length(r)))
    C <- Matrix::tcrossprod(M); (Matrix::nnzero(C) - length(g)) / 2 }
  regs <- regs_of(primary); regs_circ <- regs_of(circular)
  regulators <- intersect(names(regs), universe)
  edges <- unique(primary[vapply(strsplit(primary, " "), function(v) v[1] %in% regulators && v[2] %in% universe && v[1] != v[2], TRUE)])
  audit <- data.frame(mode = opt$yeastract_mode, E_size = length(E2), screen_additions = length(E2) - length(E),
    assoc_binding = length(bind), assoc_expression_any = length(expr_any), assoc_binding_and_any_expression = length(circular),
    primary_assoc = length(primary),
    share_circular_removed = if (length(circular)) 1 - length(intersect(circular, primary)) / length(circular) else NA,
    share_expression_only_E = if (opt$yeastract_mode == "refs" && length(expr_any)) 1 - length(expr_clean) / length(expr_any) else NA,
    share_binding_with_expression = if (length(bind)) length(intersect(bind, expr_any)) / length(bind) else NA,
    comem_pos_primary = comem_pairs(regs), comem_pos_circular = comem_pairs(regs_circ),
    share_comem_pos_removed = if (comem_pairs(regs_circ) > 0) 1 - comem_pairs(regs) / comem_pairs(regs_circ) else NA,
    tf_unmapped = sum(is.na(unique(yx$tf_orf))), target_unmapped = sum(is.na(unique(yx$target_orf))))

  ## complexes (Complex Portal 559292 + evidence filter; CYC2008 fallback)
  fetch("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv", file.path(Tdir, "559292.tsv"))
  fetch("https://rest.uniprot.org/uniprotkb/stream?query=organism_id:559292&format=tsv&fields=accession,gene_oln,gene_primary",
        file.path(Tdir, "uniprot_559292.tsv"))
  up <- utils::read.delim(file.path(Tdir, "uniprot_559292.tsv"), stringsAsFactors = FALSE, quote = "")
  orf_of <- lapply(strsplit(up[[2]], "[ ;]+"), function(v) unique(grep(ORF_RE, v, value = TRUE))); names(orf_of) <- up[[1]]
  cp <- utils::read.delim(file.path(Tdir, "559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
  xcol <- grep("cross.?ref", names(cp), ignore.case = TRUE)[1]; ecol <- grep("evidence code", names(cp), ignore.case = TRUE)[1]
  cid <- sub("\\..*$", "", cp[[1]]); parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\(.*$", "", v)); names(parts) <- cid
  cpm <- lapply(regmatches(cp[[xcol]], gregexpr("pubmed:[0-9]+", cp[[xcol]])), function(v) sub("pubmed:", "", v)); names(cpm) <- cid
  only_E <- vapply(cpm, function(v) length(v) > 0 && all(v %in% E2), TRUE)
  flatten <- function(id, seen = character(0)) { if (id %in% seen) return(character(0)); out <- character(0)
    for (p in parts[[id]]) { if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
      else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p)) }; out }
  acc_of <- lapply(names(parts), flatten); names(acc_of) <- names(parts)
  sets <- lapply(acc_of, function(a) sort(unique(unlist(orf_of[a[a %in% names(orf_of)]]))))
  keep_ids <- names(sets)[!duplicated(names(sets)) & !only_E[names(sets)]]
  filt <- function(x, lo = 3L, hi = 200L) { x <- lapply(x, intersect, universe); x[lengths(x) >= lo & lengths(x) <= hi] }
  merge_sets <- function(x) {                                   # as analysis/build_ecoli_complex_truth.R
    if (length(x) == 0) return(x); x <- x[order(names(x))]
    k <- vapply(x, paste, "", collapse = " "); x <- x[!duplicated(k)]
    drop <- vapply(seq_along(x), function(i) any(vapply(seq_along(x)[-i], function(j) length(x[[j]]) > length(x[[i]]) && all(x[[i]] %in% x[[j]]), TRUE)), TRUE)
    x <- x[!drop]; if (length(x) < 2) return(x)
    all_g <- sort(unique(unlist(x))); M <- vapply(x, function(g) all_g %in% g, logical(length(all_g)))
    inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
    out <- lapply(split(seq_along(x), cl), function(ix) sort(unique(unlist(x[ix]))))
    names(out) <- vapply(split(names(x), cl), function(n) paste(n, collapse = "+"), ""); out }
  summarize <- function(merged, src) {
    small <- merged[lengths(merged) >= 3 & lengths(merged) <= 6]
    ne <- neff_of(lengths(small)); sej <- 0.036 / sqrt(max(ne, 1e-9))
    data.frame(source = src, merged_complexes = length(merged), genes = length(unique(unlist(merged))),
               P1_complexes_with_pos = length(small), P1_positive_pairs = sum(choose(lengths(small), 2)), P1_neff = ne,
               P1_seJ_proj = sej, P1_se_proj_rho044 = sej * sqrt(1 + 0.44^2), P1_se_proj_rho1 = sej * sqrt(2), P1_se_proj_rho2 = sej * sqrt(5),
               S1_complexes_with_pos = length(merged), S1_positive_pairs = sum(choose(lengths(merged), 2)), S1_neff = neff_of(lengths(merged)))
  }
  merged <- merge_sets(filt(sets[keep_ids])); cnt <- summarize(merged, "ComplexPortal_559292"); fallback <- "not needed"
  if (cnt$P1_neff < 100) {
    if (is.na(opt$cyc2008) || !file.exists(opt$cyc2008)) { fallback <- "CYC2008 required but not available"
    } else { cy <- utils::read.delim(opt$cyc2008, stringsAsFactors = FALSE, quote = "")
      orfc <- names(which.max(vapply(cy, function(x) mean(grepl(ORF_RE, x)), 0))); cxc <- setdiff(names(cy), orfc)[1]
      cys <- lapply(split(cy[[orfc]], cy[[cxc]]), function(v) sort(unique(v))); names(cys) <- paste0("CYC2008:", names(cys))
      merged <- merge_sets(filt(c(sets[keep_ids], cys))); cnt <- rbind(cnt, summarize(merged, "ComplexPortal_559292+CYC2008"))
      fallback <- "union used" } }
  fin <- cnt[nrow(cnt), ]
  out <- cbind(fin, fallback = fallback, complexes_dropped_only_E = sum(only_E),
    P2_regulons = length(regs), P2_genes = length(unique(unlist(regs))), P2_positive_pairs = audit$comem_pos_primary,
    P2_neff = neff_of(lengths(regs)), S6_regulators = length(regulators), S6_positive_edges = length(edges),
    evaluable_P1 = fin$P1_positive_pairs >= 200 && fin$P1_complexes_with_pos >= 30,
    evaluable_S1 = fin$S1_positive_pairs >= 200 && fin$S1_complexes_with_pos >= 30,
    evaluable_P2 = audit$comem_pos_primary >= 200, evaluable_S6 = length(edges) >= 200,
    power_gate_pass = fin$P1_neff >= 100, prereg_commit = hash, test_script_blob = test_blob)
  out$verdict1 <- !(out$power_gate_pass && out$evaluable_P1 && out$evaluable_P2)
  ## write
  mem <- do.call(rbind, lapply(names(merged), function(k) data.frame(complex = k, orf = merged[[k]])))
  utils::write.table(mem, file.path(SEL, "yeast_complexes_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rmem <- do.call(rbind, lapply(names(regs), function(k) data.frame(regulator = k, orf = regs[[k]])))
  utils::write.table(rmem, file.path(SEL, "yeast_regulons_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(do.call(rbind, strsplit(edges, " ")), file.path(SEL, "yeast_tf_edges.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = c("tf", "target"))
  utils::write.csv(out, file.path(SEL, "truth_structure_counts.csv"), row.names = FALSE)
  utils::write.csv(audit, file.path(SEL, "contamination_audit.csv"), row.names = FALSE)
  utils::write.csv(value_table, file.path(SEL, "yeastract_value_table.csv"), row.names = FALSE)
  utils::write.csv(data.frame(complex = cid, eco = if (!is.na(ecol)) cp[[ecol]] else NA, pmids = vapply(cpm, paste, "", collapse = ";"),
                              only_E = only_E), file.path(SEL, "complex_evidence.csv"), row.names = FALSE)
  writeLines(c(sprintf("%s  yeast_complexes_members.tsv", sha(file.path(SEL, "yeast_complexes_members.tsv"))),
               sprintf("%s  yeast_regulons_members.tsv", sha(file.path(SEL, "yeast_regulons_members.tsv"))),
               sprintf("%s  yeast_tf_edges.tsv", sha(file.path(SEL, "yeast_tf_edges.tsv"))),
               sprintf("%s  truth_raw/559292.tsv", sha(file.path(Tdir, "559292.tsv"))),
               sprintf("%s  truth_raw/uniprot_559292.tsv", sha(file.path(Tdir, "uniprot_559292.tsv"))),
               sprintf("%s  %s", sha(opt$yeastract), basename(opt$yeastract))), file.path(SEL, "yeast_truths.sha256"))
  say("power gate: ", if (out$power_gate_pass) "PASS" else "FAIL", "; verdict 1: ", out$verdict1,
      sprintf(" (P1 n_eff %.1f, P2 n_eff %.1f)", out$P1_neff, out$P2_neff), "; commit ", SEL, " before any scoring")
}
