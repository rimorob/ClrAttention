#!/usr/bin/env Rscript
# Data acquisition and label-free preparation for docs/preregistration_reciprocity_yeast.md.
# DRAFT companion of a DRAFT preregistration: written 2026-11-01, NOT yet run end to end.
#
# Stages (run in order; each writes into data/yeast/ or docs/reciprocity_yeast_selections/):
#   download         expression + wild-type controls + annotation only (no truth). SHA-256 -> data/yeast/SHA256SUMS
#   matrix           700 responsive mutants (GSE42527 published M) -> ORF x mutant matrix, deleted-gene masking
#   wt               wild-type controls (GSE45115 sample tables) -> ORF x culture matrix (orientation from `check`)
#   check            implementation check: rebuild first 10 GSE42527 mutants from their sample tables vs published M
#   splits           10 random halves of the mutants (seed 20261101) -> docs/reciprocity_yeast_selections/splits.csv
#   feasibility      gate F: cross-half top-1% CLR reproducibility, mutants (N matched) vs wild types
#   truth_structure  GATED (preregistration fixed, tracked, identical to HEAD, no "STATUS: DRAFT"):
#                    downloads Complex Portal 559292 + UniProt mapping, reads a YEASTRACT export, builds member
#                    lists, counts positives / complexes / n_eff, applies the CYC2008 fallback rule and writes the
#                    power-gate line. Reads membership only; never loads expression values or any network.
#
# No stage scores anything. Stages other than truth_structure never download or open a truth.
# Usage: Rscript analysis/yeast_prep.R <stage> [--threads 4] [--yeastract path] [--cyc2008 path]
# Run every stage through jobs_pc/sched (# mem / # cpus headers); download is network-bound and light.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript analysis/yeast_prep.R <download|matrix|wt|check|splits|feasibility|truth_structure> [options]")
stage <- args[1]
opt <- list(threads = 4L, yeastract = NA_character_, cyc2008 = NA_character_)
i <- 2L
while (i < length(args)) { k <- sub("^--", "", args[i]); v <- args[i + 1L]
  opt[[k]] <- if (k == "threads") as.integer(v) else v; i <- i + 2L }

PREREG <- "docs/preregistration_reciprocity_yeast.md"
D <- "data/yeast"; SEL <- "docs/reciprocity_yeast_selections"
dir.create(D, showWarnings = FALSE, recursive = TRUE); dir.create(SEL, showWarnings = FALSE, recursive = TRUE)
GEO <- "https://ftp.ncbi.nlm.nih.gov/geo"
ACC <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sha <- function(f) digest::digest(file = f, algo = "sha256")
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
  fs <- list.files(D, recursive = TRUE, full.names = TRUE); fs <- fs[basename(fs) != "SHA256SUMS"]
  writeLines(sprintf("%s  %s", vapply(fs, sha, ""), sub(paste0("^", D, "/"), "", fs)), file.path(D, "SHA256SUMS"))
  say("hashes: ", length(fs), " files -> ", file.path(D, "SHA256SUMS"))
}
# replicate grouping: mutants by deleted gene (the two cultures may be different isolates, e.g.
# rad18-del-4-b with rad18-del-3-a), wild types by dye-swap pair (wt-SV013-a / -b)
group_key <- function(title) ifelse(grepl("-del-", title), sub("-del-.*$", "", title), sub("-[ab]$", "", title))
ORF_RE <- "^(Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4}|R[0-9]{4}[WC]?)$"

## ---- GEO helpers ------------------------------------------------------------------------------
series_meta <- function(gse) {                      # Sample_title + Sample_geo_accession (+ channel sources)
  f <- file.path(D, paste0(gse, "_series_matrix.txt.gz"))
  L <- readLines(gzfile(f)); L <- L[startsWith(L, "!Sample_")]
  get <- function(key) { x <- L[startsWith(L, paste0(key, "\t"))]; if (!length(x)) return(NULL)
    gsub('^"|"$', "", strsplit(x[1], "\t")[[1]][-1]) }
  data.frame(title = get("!Sample_title"), gsm = get("!Sample_geo_accession"),
             src1 = get("!Sample_source_name_ch1"), src2 = get("!Sample_source_name_ch2"), stringsAsFactors = FALSE)
}
gsm_table <- function(gsm) {                        # VALUE column of a GEO sample table (cached)
  f <- file.path(D, "gsm", paste0(gsm, ".txt"))
  if (!file.exists(f)) { dir.create(dirname(f), showWarnings = FALSE)
    fetch(sprintf("%s?acc=%s&targ=self&form=text&view=data", ACC, gsm), f); Sys.sleep(0.4) }
  L <- readLines(f); a <- which(L == "!sample_table_begin"); b <- which(L == "!sample_table_end")
  if (!length(a)) stop("no sample table in ", f)
  tb <- utils::read.delim(text = L[(a + 1L):(if (length(b)) b - 1L else length(L))], stringsAsFactors = FALSE)
  stats::setNames(as.numeric(tb$VALUE), as.character(tb$ID_REF))
}
platform <- function() {                            # spot ID -> SCAB probe -> ORF, detected by value pattern
  f <- file.path(D, "GPL11232_full.txt"); L <- readLines(f)
  a <- which(L == "!platform_table_begin"); b <- which(L == "!platform_table_end")
  tb <- utils::read.delim(text = L[(a + 1L):(b - 1L)], stringsAsFactors = FALSE, quote = "", comment.char = "")
  frac <- function(re) vapply(tb, function(x) mean(grepl(re, as.character(x))), 0)
  scab <- names(which.max(frac("^SCAB[0-9]{6}$"))); orf <- names(which.max(frac(ORF_RE)))
  if (max(frac("^SCAB[0-9]{6}$")) < 0.5 || max(frac(ORF_RE)) < 0.3)
    stop("GPL11232: could not detect SCAB / ORF columns; inspect ", f, " and set them explicitly")
  data.frame(id = as.character(tb$ID), scab = tb[[scab]], orf = ifelse(grepl(ORF_RE, tb[[orf]]), tb[[orf]], NA),
             stringsAsFactors = FALSE)
}
by_orf <- function(M, orf) {                        # average rows (spots/probes) sharing an ORF
  keep <- !is.na(orf); M <- M[keep, , drop = FALSE]; orf <- orf[keep]
  S <- rowsum(ifelse(is.na(M), 0, M), orf); N <- rowsum(1 * !is.na(M), orf); S / ifelse(N == 0, NA, N)
}

## ---- stage: download ----------------------------------------------------------------------------
if (stage == "download") {
  for (g in c("GSE42527", "GSE42526")) {
    fetch(sprintf("%s/series/GSE42nnn/%s/suppl/%s_final_GeneExpressionMatrix.txt.gz", GEO, g, g),
          file.path(D, paste0(g, "_final_GeneExpressionMatrix.txt.gz")))
    fetch(sprintf("%s/series/GSE42nnn/%s/suppl/%s_protocols.txt.gz", GEO, g, g), file.path(D, paste0(g, "_protocols.txt.gz")))
    fetch(sprintf("%s/series/GSE42nnn/%s/matrix/%s_series_matrix.txt.gz", GEO, g, g), file.path(D, paste0(g, "_series_matrix.txt.gz")))
  }
  fetch(sprintf("%s/series/GSE45nnn/GSE45115/matrix/GSE45115_series_matrix.txt.gz", GEO), file.path(D, "GSE45115_series_matrix.txt.gz"))
  fetch(sprintf("%s?acc=GPL11232&targ=self&form=text&view=data", ACC), file.path(D, "GPL11232_full.txt"))
  fetch("https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab", file.path(D, "SGD_features.tab"))
  wt <- series_meta("GSE45115"); say("GSE45115: ", nrow(wt), " wild-type hybridizations")
  for (g in wt$gsm) invisible(gsm_table(g))
  mu <- series_meta("GSE42527"); key <- group_key(mu$title)
  first10 <- sort(unique(key))[1:10]
  for (g in mu$gsm[key %in% first10]) invisible(gsm_table(g))
  record_hashes(); say("download done; commit data/yeast/SHA256SUMS (data files themselves are not committed)")
}

## ---- stage: matrix ------------------------------------------------------------------------------
read_final <- function(g) {
  x <- utils::read.delim(gzfile(file.path(D, paste0(g, "_final_GeneExpressionMatrix.txt.gz"))), check.names = FALSE,
                         stringsAsFactors = FALSE)
  mcol <- grep("^M \\(", names(x)); pcol <- grep("^p-value \\(", names(x))
  M <- as.matrix(x[, mcol]); P <- as.matrix(x[, pcol]); rownames(M) <- rownames(P) <- x$PROBE_ID
  lab <- sub("^M \\((.*?)\\).*$", "\\1", names(x)[mcol])            # first replicate label, e.g. spo7-del-2-b
  colnames(M) <- colnames(P) <- sub("-del-.*$", "", lab)             # deleted gene name (lower case)
  list(M = M, P = P)
}
if (stage == "matrix") {
  pf <- platform(); sgd <- utils::read.delim(file.path(D, "SGD_features.tab"), header = FALSE, stringsAsFactors = FALSE, quote = "")
  name2orf <- stats::setNames(sgd$V4, toupper(sgd$V5)); name2orf <- c(name2orf, stats::setNames(sgd$V4, toupper(sgd$V4)))
  build <- function(g) {
    r <- read_final(g); X <- by_orf(r$M, pf$orf[match(rownames(r$M), pf$scab)])
    del <- name2orf[toupper(colnames(X))]
    unresolved <- colnames(X)[is.na(del) | !del %in% rownames(X)]
    for (j in which(!is.na(del) & del %in% rownames(X))) X[del[j], j] <- NA      # mask deleted gene's own entry
    miss <- rowMeans(is.na(X)); X <- X[miss <= 0.05, , drop = FALSE]
    med <- apply(X, 1, stats::median, na.rm = TRUE); ix <- which(is.na(X), arr.ind = TRUE); X[ix] <- med[ix[, 1]]
    X <- X[apply(X, 1, function(v) diff(range(v))) > 0, , drop = FALSE]
    list(X = X, unresolved = unresolved, dropped_missing = sum(miss > 0.05))
  }
  b <- build("GSE42527")
  utils::write.table(data.frame(orf = rownames(b$X), b$X, check.names = FALSE), file.path(D, "deleteome_responsive_700.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(b$unresolved, file.path(SEL, "unresolved_deleted_gene_labels.txt"))
  say(sprintf("responsive matrix: %d ORFs x %d mutants; %d labels unresolved; %d ORFs dropped (>5%% missing)",
              nrow(b$X), ncol(b$X), length(b$unresolved), b$dropped_missing))
  a <- build("GSE42526"); common <- intersect(rownames(b$X), rownames(a$X))           # descriptive sensitivity set
  utils::write.table(data.frame(orf = common, cbind(b$X[common, ], a$X[common, ]), check.names = FALSE),
                     file.path(D, "deleteome_all_1484.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(sprintf("%s  %s", sha(file.path(D, "deleteome_responsive_700.tsv")), "deleteome_responsive_700.tsv"),
             file.path(SEL, "deleteome_matrix.sha256"))
  say("commit ", SEL, "/deleteome_matrix.sha256 and the unresolved-label log")
}

## ---- stage: check (orientation rule + reproduction of published M) ----------------------------
# Candidate rules (fixed in advance): "asis" = VALUE as deposited for both hybridizations;
# "flip_ch1" = negate VALUE of the hybridization whose channel-1 source is the test strain.
rebuild <- function(meta, rule, pf) {
  stem <- group_key(meta$title)
  vals <- lapply(seq_len(nrow(meta)), function(k) { v <- gsm_table(meta$gsm[k])
    if (rule == "flip_ch1" && !grepl("ref", meta$src1[k], ignore.case = TRUE)) v <- -v; v })
  ids <- Reduce(union, lapply(vals, names))
  V <- vapply(vals, function(v) v[ids], numeric(length(ids))); rownames(V) <- ids
  V <- vapply(split(seq_len(ncol(V)), stem), function(ix) rowMeans(V[, ix, drop = FALSE], na.rm = TRUE), numeric(nrow(V)))
  rownames(V) <- ids; by_orf(V, pf$orf[match(rownames(V), pf$id)])
}
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
  writeLines(c(rule, if (any(pass)) "published_M_reproduced" else "NOT_reproduced_gate_uses_rebuilt_mutants"),
             file.path(SEL, "orientation_rule.txt"))
  say("orientation rule: ", rule, "; reproduction ", if (any(pass)) "PASSED" else "FAILED (gate F uses rebuilt mutants)")
}

## ---- stage: wt ----------------------------------------------------------------------------------
if (stage == "wt") {
  rf <- file.path(SEL, "orientation_rule.txt"); if (!file.exists(rf)) stop("run stage check first")
  rule <- readLines(rf, 1L); pf <- platform(); W <- rebuild(series_meta("GSE45115"), rule, pf)
  X <- utils::read.delim(file.path(D, "deleteome_responsive_700.tsv"), check.names = FALSE)
  W <- W[intersect(X$orf, rownames(W)), , drop = FALSE]
  W <- W[rowMeans(is.na(W)) <= 0.05, , drop = FALSE]
  med <- apply(W, 1, stats::median, na.rm = TRUE); ix <- which(is.na(W), arr.ind = TRUE); W[ix] <- med[ix[, 1]]
  utils::write.table(data.frame(orf = rownames(W), W, check.names = FALSE), file.path(D, "wt_controls_GSE45115.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  say(sprintf("wild-type matrix: %d ORFs x %d cultures (rule %s)", nrow(W), ncol(W), rule))
}

## ---- stage: splits ------------------------------------------------------------------------------
if (stage == "splits") {
  X <- utils::read.delim(file.path(D, "deleteome_responsive_700.tsv"), check.names = FALSE, nrows = 1)
  mut <- colnames(X)[-1]; set.seed(20261101)
  sp <- vapply(1:10, function(s) { h <- integer(length(mut)); h[sample(length(mut), floor(length(mut) / 2))] <- 1L; h + 1L },
               integer(length(mut)))
  colnames(sp) <- sprintf("split_%02d", 1:10)
  utils::write.csv(data.frame(mutant = mut, sp), file.path(SEL, "splits.csv"), row.names = FALSE)
  say("10 random halves of ", length(mut), " mutants -> ", SEL, "/splits.csv (commit)")
}

## ---- stage: feasibility (gate F) ----------------------------------------------------------------
if (stage == "feasibility") {
  suppressPackageStartupMessages(library(clr))
  rd <- function(f) { x <- utils::read.delim(f, check.names = FALSE); m <- as.matrix(x[, -1]); rownames(m) <- x$orf; m }
  rule_ok <- readLines(file.path(SEL, "orientation_rule.txt"))[2] == "published_M_reproduced"
  if (!rule_ok) stop("implementation check failed: rebuild all 700 mutants from sample tables first (not yet implemented)")
  Xm <- rd(file.path(D, "deleteome_responsive_700.tsv")); W <- rd(file.path(D, "wt_controls_GSE45115.tsv"))
  g <- intersect(rownames(Xm), rownames(W)); Xm <- Xm[g, ]; W <- W[g, ]
  fitC <- function(M) { rng <- apply(M, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) M[rng <= 0, ] <- M[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(M), sd = 1e-9)
    f <- ClrAttention$new(M)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
    f$calibrate(method = "normal", combine = "stouffer"); C <- f$clr_scores; C[lower.tri(C, diag = TRUE)] <- NA; C }
  top <- function(C) { v <- which(!is.na(C)); k <- ceiling(0.01 * length(v)); v[order(C[v], decreasing = TRUE)[seq_len(k)]] }
  repro <- function(M, cols) { h <- sample(cols, floor(length(cols) / 2)); a <- top(fitC(M[, h])); b <- top(fitC(M[, setdiff(cols, h)]))
    length(intersect(a, b)) / length(a) }
  Nw <- ncol(W); set.seed(20261104)
  Rm <- vapply(1:10, function(s) repro(Xm, sample(ncol(Xm), Nw)), 0)
  Rw <- vapply(1:10, function(s) repro(W, seq_len(Nw)), 0)
  tm <- top(fitC(Xm)); tw <- top(fitC(W))
  out <- data.frame(N_matched = Nw, genes = length(g), median_R_mut = stats::median(Rm), median_R_wt = stats::median(Rw),
                    min_R_mut = min(Rm), max_R_wt = max(Rw), share_mut_top_in_wt_top = length(intersect(tm, tw)) / length(tm))
  out$pass <- out$median_R_mut >= 2 * out$median_R_wt && out$min_R_mut > out$max_R_wt
  utils::write.csv(out, file.path(SEL, "feasibility_gate.csv"), row.names = FALSE)
  utils::write.csv(data.frame(draw = 1:10, R_mut = Rm, R_wt = Rw), file.path(SEL, "feasibility_draws.csv"), row.names = FALSE)
  say("gate F: ", if (out$pass) "PASS" else "FAIL (verdict 1)", sprintf(" (median R mut %.3f vs wt %.3f)", out$median_R_mut, out$median_R_wt))
}

## ---- stage: truth_structure (GATED) --------------------------------------------------------------
if (stage == "truth_structure") {
  if (system2("git", c("ls-files", "--error-unmatch", PREREG), stdout = FALSE, stderr = FALSE) != 0L) stop("gate: ", PREREG, " not tracked")
  if (system2("git", c("diff", "--quiet", "HEAD", "--", PREREG)) != 0L) stop("gate: ", PREREG, " differs from HEAD")
  if (any(grepl("STATUS: DRAFT", readLines(PREREG), fixed = TRUE))) stop("gate: ", PREREG, " is still a DRAFT; no truth may be downloaded")
  for (f in c("feasibility_gate.csv", "splits.csv"))
    if (system2("git", c("ls-files", "--error-unmatch", file.path(SEL, f)), stdout = FALSE, stderr = FALSE) != 0L) stop("gate: ", f, " not committed")
  if (!isTRUE(as.logical(utils::read.csv(file.path(SEL, "feasibility_gate.csv"))$pass))) stop("gate F failed: verdict 1, stop")
  if (is.na(opt$yeastract) || !file.exists(opt$yeastract)) stop("pass --yeastract <export file> (manual web export, parameters logged)")
  hash <- system2("git", c("log", "-1", "--format=%H", "--", PREREG), stdout = TRUE)
  source("analysis/regulons.R")
  T <- file.path(D, "truth_raw"); dir.create(T, showWarnings = FALSE)
  fetch("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv", file.path(T, "559292.tsv"))
  fetch("https://rest.uniprot.org/uniprotkb/stream?query=organism_id:559292&format=tsv&fields=accession,gene_oln,gene_primary",
        file.path(T, "uniprot_559292.tsv"))
  universe <- utils::read.delim(file.path(D, "deleteome_responsive_700.tsv"), check.names = FALSE)$orf   # row names only
  up <- utils::read.delim(file.path(T, "uniprot_559292.tsv"), stringsAsFactors = FALSE, quote = "")
  orf_of <- lapply(strsplit(up[[2]], "[ ;]+"), function(v) unique(grep(ORF_RE, v, value = TRUE))); names(orf_of) <- up[[1]]
  cp <- utils::read.delim(file.path(T, "559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
  cid <- sub("\\..*$", "", cp[[1]]); parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\(.*$", "", v)); names(parts) <- cid
  flatten <- function(id, seen = character(0)) { if (id %in% seen) return(character(0)); out <- character(0)
    for (p in parts[[id]]) { if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
      else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p)) }; out }
  sets <- lapply(names(parts), function(id) { a <- flatten(id); sort(unique(unlist(orf_of[a[a %in% names(orf_of)]]))) })
  names(sets) <- names(parts); sets <- sets[!duplicated(names(sets))]
  filt <- function(s, lo = 3L, hi = 200L) { s <- lapply(s, intersect, universe); s[lengths(s) >= lo & lengths(s) <= hi] }
  merge_sets <- function(s) {                                   # as analysis/build_ecoli_complex_truth.R
    if (length(s) == 0) return(s); s <- s[order(names(s))]
    key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
    drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j) length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
    s <- s[!drop]; if (length(s) < 2) return(s)
    all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
    inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
    cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
    out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
    names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), ""); out }
  neff_of <- function(s) { p <- choose(lengths(s), 2); if (!length(p) || sum(p) == 0) return(0); sum(p)^2 / sum(p^2) }
  summarize <- function(merged, src) {
    small <- merged[lengths(merged) >= 3 & lengths(merged) <= 6]
    ne <- neff_of(small); se <- 0.036 / sqrt(max(ne, 1e-9)) * sqrt(1 + 0.44^2)
    data.frame(source = src, merged_complexes = length(merged), genes = length(unique(unlist(merged))),
               P1_small_complexes = length(small), P1_positive_pairs = sum(choose(lengths(small), 2)), P1_neff = ne,
               P1_se_proj = se, P1_mde80 = 2.8016 * se, S1_positive_pairs = sum(choose(lengths(merged), 2)), S1_neff = neff_of(merged))
  }
  merged <- merge_sets(filt(sets)); cnt <- summarize(merged, "ComplexPortal_559292")
  if (cnt$P1_neff < 100) {                                       # fallback fixed in the preregistration
    if (is.na(opt$cyc2008) || !file.exists(opt$cyc2008)) { cnt$fallback <- "CYC2008 required but not available"
    } else { cy <- utils::read.delim(opt$cyc2008, stringsAsFactors = FALSE, quote = "")
      orfc <- names(which.max(vapply(cy, function(x) mean(grepl(ORF_RE, x)), 0))); cxc <- setdiff(names(cy), orfc)[1]
      cys <- lapply(split(cy[[orfc]], cy[[cxc]]), function(v) sort(unique(v))); names(cys) <- paste0("CYC2008:", names(cys))
      merged <- merge_sets(filt(c(sets, cys))); cnt <- rbind(cnt, summarize(merged, "ComplexPortal_559292+CYC2008"))
      cnt$fallback <- "union used" }
  } else cnt$fallback <- "not needed"
  # Regulons: YEASTRACT export (semicolon CSV binary matrix, genes x TFs, or long TF;gene table). FORMAT TO VERIFY.
  yx <- utils::read.csv(opt$yeastract, sep = ";", check.names = FALSE, stringsAsFactors = FALSE)
  pairs <- if (ncol(yx) == 2) yx else { m <- as.matrix(yx[, -1]); w <- which(m == 1, arr.ind = TRUE)
    data.frame(tf = colnames(m)[w[, 2]], gene = yx[[1]][w[, 1]]) }
  names(pairs) <- c("tf", "gene"); sgd <- utils::read.delim(file.path(D, "SGD_features.tab"), header = FALSE, stringsAsFactors = FALSE, quote = "")
  to_orf <- function(x) { x <- toupper(x); ifelse(grepl(ORF_RE, x), x, stats::setNames(sgd$V4, toupper(sgd$V5))[x]) }
  pairs$tf_orf <- to_orf(pairs$tf); pairs$gene_orf <- to_orf(pairs$gene)
  regs <- lapply(split(pairs$gene_orf, pairs$tf_orf), function(v) intersect(unique(v), universe)); regs <- regs[lengths(regs) >= 5 & lengths(regs) <= 500]
  cnt$P2_regulons <- length(regs); cnt$P2_genes <- length(unique(unlist(regs)))
  cnt$P2_tf_unmapped <- sum(is.na(unique(pairs$tf_orf))); cnt$power_gate_pass <- utils::tail(cnt$P1_neff, 1) >= 100
  cnt$prereg_commit <- hash
  mem <- do.call(rbind, lapply(names(merged), function(k) data.frame(complex = k, orf = merged[[k]])))
  utils::write.table(mem, file.path(SEL, "yeast_complexes_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  rmem <- do.call(rbind, lapply(names(regs), function(k) data.frame(regulator = k, orf = regs[[k]])))
  utils::write.table(rmem, file.path(SEL, "yeast_regulons_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.csv(cnt, file.path(SEL, "truth_structure_counts.csv"), row.names = FALSE)
  writeLines(c(sprintf("%s  yeast_complexes_members.tsv", sha(file.path(SEL, "yeast_complexes_members.tsv"))),
               sprintf("%s  yeast_regulons_members.tsv", sha(file.path(SEL, "yeast_regulons_members.tsv"))),
               sprintf("%s  truth_raw/559292.tsv", sha(file.path(T, "559292.tsv"))),
               sprintf("%s  truth_raw/uniprot_559292.tsv", sha(file.path(T, "uniprot_559292.tsv"))),
               sprintf("%s  %s", sha(opt$yeastract), basename(opt$yeastract))), file.path(SEL, "yeast_truths.sha256"))
  say("power gate: ", if (utils::tail(cnt$power_gate_pass, 1)) "PASS" else "FAIL (verdict 1)",
      sprintf(" (n_eff %.1f, SE_proj %.4f)", utils::tail(cnt$P1_neff, 1), utils::tail(cnt$P1_se_proj, 1)), "; commit ", SEL)
}
