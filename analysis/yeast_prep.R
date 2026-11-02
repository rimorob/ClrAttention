#!/usr/bin/env Rscript
# Data acquisition and label-free preparation for docs/preregistration_reciprocity_yeast.md (DRAFT rev 5).
# Written 2026-11-01; parse-checked; self-tests selftest_gate and selftest_git pass; NOT yet run on real data.
#
# Stages, in protocol order (outputs in data/yeast/ [not committed; hashes are] or
# docs/reciprocity_yeast_selections/ [committed]):
#   selftest_gate      reports "fixed" or "draft" for the preregistration (or --prereg <file>) and checks the marker
#                      logic on derived copies (fixed <-> draft, invalid date 2026-99-99, placeholder date)
#   selftest_git       builds a throw-away git repository and checks the commit-order gates (strict ancestry,
#                      adding commit, same-commit refusal, re-committed log refusal, shallow-clone refusal)
#   download           expression + wild-type controls + annotation only (no truth). SHA-256 -> data/yeast/SHA256SUMS
#   matrix             700 responsive mutants (published M) -> ORF x mutant matrix; deleted gene resolved via SGD
#                      (unique aliases), verified by its own M (<= -1 and in the column's lowest 1%), else dropped;
#                      masking of the deleted ORF and ORFs overlapping / within 500 bp; divergent-pair list
#   orient_split       DEVIATION 2: fixed 50/50 split of balanced mutants into fitting / held-out sets (seed 20261106)
#   check              DEVIATION 2: per-sub-series orientation rule fitted on the fitting set (>= 5 mutants, else
#                      "asis"); floor: held-out median r >= 0.8; structural dye check part 1 stops, part 2 = balanced only
#   rebuild            only after a failed (but floor-passing) check: 376 random mutants (seed 20261105) rebuilt
#   wt                 wild-type controls (GSE45115 sample tables) with the same rule and the structural check
#   splits             10 random halves (seed 20261101); label-free mutant clusters (complete linkage, 1 - r,
#                      h 0.5); the 200 experiment-jackknife replicates (seeds 20262000 + b; cluster-stratified
#                      fallback if K < 20 or a cluster holds > 25% of mutants)
#   feasibility        gate F (sanity check): chance-corrected R' (mutants vs WT) and overlap share Q; pass iff
#                      median R'_mut > median R'_wt and Q <= 0.5; WT-percentile diagnostic
#   complex_structure  8a, GATED: Complex Portal (+ CYC2008 fallback) membership, evidence filter, merge,
#                      P1/S1 evaluability, n_eff(P1), power gate; never builds regulon membership
#   regulon_structure  8c, GATED (8a committed and pushed; the commit ADDING the export log strictly after 8a and
#                      steps 4-7): runs the committed converter on the raw export files whose SHA-256 are listed in
#                      the export log -> computed mode, E + reanalysis screen, regulons, TF edges, audit, n_eff(P2)
# No stage scores anything. Stages before complex_structure never download or open a truth. Truth stages read
# membership and the universe's row names only.
# Usage: Rscript analysis/yeast_prep.R <stage> [--threads 4]
#        complex_structure: [--cyc2008 <file>]
#        regulon_structure: --yeastract_mode refs|binding_only (cross-check only); raw export files are read from the
#                           export log lines "raw_export: <sha256> <kind> <path>" (kind: table|binding_matrix|expression_matrix)
# Before any truth stage: git fetch origin (the gates require every gate commit to be on the upstream branch).

invisible(Sys.setlocale("LC_COLLATE", "C"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
args <- commandArgs(trailingOnly = TRUE)
STAGES <- c("selftest_gate", "selftest_git", "download", "matrix", "orient_split", "check", "rebuild", "wt", "splits", "feasibility",
            "complex_structure", "regulon_structure")
if (!length(args) || !args[1] %in% STAGES) stop("usage: Rscript analysis/yeast_prep.R <", paste(STAGES, collapse = "|"), "> [options]")
stage <- args[1]
opt <- list(threads = 4L, yeastract_mode = NA_character_, cyc2008 = NA_character_, prereg = NA_character_)
i <- 2L
while (i < length(args)) { k <- sub("^--", "", args[i]); v <- args[i + 1L]
  opt[[k]] <- if (k == "threads") as.integer(v) else v; i <- i + 2L }

PREREG <- "docs/preregistration_reciprocity_yeast.md"
PREP_SCRIPT <- "analysis/yeast_prep.R"
TEST_SCRIPT <- "analysis/reciprocity_yeast_test.R"
CONVERTER <- "analysis/yeastract_convert.R"
D <- "data/yeast"; SEL <- "docs/reciprocity_yeast_selections"
GEO <- "https://ftp.ncbi.nlm.nih.gov/geo"
ACC <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
FIXED_RE <- "^> \\*Fixed [0-9]{4}-[0-9]{2}-[0-9]{2}\\*$"
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
sha <- function(f) digest::digest(file = f, algo = "sha256")
mkd <- function() { dir.create(D, showWarnings = FALSE, recursive = TRUE); dir.create(SEL, showWarnings = FALSE, recursive = TRUE) }
env_info <- function() c(R = R.version.string, RNG = paste(RNGkind(), collapse = "/"), collate = Sys.getlocale("LC_COLLATE"))
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
GITC <- character(0)            # set to c("-C", dir) by selftest_git
git1 <- function(...) suppressWarnings(system2("git", c(GITC, ...), stdout = TRUE, stderr = FALSE))
gitq <- function(...) system2("git", c(GITC, ...), stdout = FALSE, stderr = FALSE) == 0L
tracked <- function(f) gitq("ls-files", "--error-unmatch", f)
clean <- function(f) tracked(f) && gitq("diff", "--quiet", "HEAD", "--", f)
last_commit <- function(f) { h <- git1("log", "-1", "--format=%H", "--", f); if (length(h)) h[1] else NA_character_ }
adding_commit <- function(f) { h <- git1("log", "--diff-filter=A", "--format=%H", "--", f); if (length(h)) h[length(h)] else NA_character_ }
is_ancestor <- function(a, b) gitq("merge-base", "--is-ancestor", a, b)
strictly_before <- function(a, b) !is.na(a) && !is.na(b) && a != b && is_ancestor(a, b)
refuse_shallow <- function() if (identical(git1("rev-parse", "--is-shallow-repository"), "true")) stop("gate: shallow clone; commit-order gates cannot be verified")
on_upstream <- function(cmt) { up <- git1("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
  if (!length(up) || !nzchar(up[1])) stop("gate: no upstream branch; push gate commits and git fetch first")
  is_ancestor(cmt, up[1]) }
DRAFT_HDR_RE <- "^> \\*\\*STATUS: DRAFT"
prereg_status <- function(lines) {   # "fixed" iff exactly one valid fixed-date marker with a real date and no DRAFT header line
  m <- grep(FIXED_RE, lines, value = TRUE)
  ok <- length(m) == 1L && !any(grepl(DRAFT_HDR_RE, lines)) && {
    d <- sub("^> \\*Fixed ([0-9-]+)\\*$", "\\1", m); dd <- suppressWarnings(as.Date(d, format = "%Y-%m-%d"))
    !is.na(dd) && format(dd, "%Y-%m-%d") == d }
  if (ok) "fixed" else "draft"
}
prereg_is_fixed <- function(lines) prereg_status(lines) == "fixed"
gate_prereg <- function() {
  refuse_shallow()
  if (!clean(PREREG)) stop("gate: ", PREREG, " not tracked or differs from HEAD")
  if (!prereg_is_fixed(readLines(PREREG))) stop("gate: ", PREREG, " carries no valid fixed-date marker; no truth may be touched")
  h <- last_commit(PREREG)
  ex <- file.path(SEL, "excluded_pmids.txt")             # E committed at or before the fixing commit
  if (!clean(ex) || !is_ancestor(adding_commit(ex), h)) stop("gate: ", ex, " must be committed at or before the commit that fixed ", PREREG)
  h
}
need_clean <- function(fs) for (f in fs) if (!clean(f)) stop("gate: ", f, " not committed or differs from HEAD")
ORF_RE <- "^(Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4}|R[0-9]{4}[WC]?)$"
REF_RE <- "ref"
group_key <- function(title) {   # strips GEO prefixes such as "[hs1991] "; mutants by deleted gene, WT by title stem
  t <- sub("^\\[[^]]*\\] *", "", title)
  ifelse(grepl("-del", t), sub("-del.*$", "", t), sub("-[a-z0-9]+$", "", t)) }
rd_mat <- function(f) { x <- utils::read.delim(f, check.names = FALSE); m <- as.matrix(x[, -1]); rownames(m) <- x$orf; m }
wr_mat <- function(m, f) utils::write.table(data.frame(orf = rownames(m), m, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = FALSE)
universe_rows <- function() { L <- readLines(file.path(D, "deleteome_responsive_700.tsv")); sub("\t.*$", "", L[-1]) }  # row names only
SEL_FILES <- file.path(SEL, c("yeast_selections.csv", "hepg2_selections.csv", "jurkat_selections.csv", "k562_selections.csv"))

## ---- stage: selftest_gate (BL1) -----------------------------------------------------------------
if (stage == "selftest_gate") {
  f <- if (is.na(opt$prereg)) PREREG else opt$prereg
  L <- readLines(f); now <- prereg_status(L)
  hdr <- grep(DRAFT_HDR_RE, L); fx <- grep(FIXED_RE, L)
  k <- if (length(hdr)) hdr[1] else if (length(fx)) fx[1] else stop("selftest: neither a DRAFT header nor a fixed marker found")
  mk <- function(line) { x <- L; x[k] <- line; x[setdiff(c(hdr, fx), k)] <- "" ; x }
  res <- c(fixed_copy = prereg_status(mk("> *Fixed 2026-11-02*")) == "fixed",
           draft_copy = prereg_status(mk("> **STATUS: DRAFT (selftest copy).**")) == "draft",
           invalid_date = prereg_status(mk("> *Fixed 2026-99-99*")) == "draft",
           feb30 = prereg_status(mk("> *Fixed 2026-02-30*")) == "draft",
           placeholder = prereg_status(mk("> *Fixed YYYY-MM-DD*")) == "draft",
           marker_plus_draft = prereg_status(c(mk("> *Fixed 2026-11-02*"), "> **STATUS: DRAFT stray**")) == "draft")
  print(res)
  if (!all(res)) stop("selftest FAILED: marker logic is wrong")
  cat(now, "\n"); say("selftest passed; ", f, " is ", toupper(now))
}

## ---- stage: selftest_git (commit-order gates in a throw-away repository) ------------------------
if (stage == "selftest_git") {
  td <- tempfile("gitgate"); dir.create(td); GITC <<- c("-C", td)
  g <- function(...) { r <- system2("git", c(GITC, ...), stdout = TRUE, stderr = TRUE); r }
  g("init", "-q"); g("config", "user.email", "t@t"); g("config", "user.name", "t")
  w <- function(f, txt) writeLines(txt, file.path(td, f))
  w("a.txt", "a"); g("add", "a.txt"); g("commit", "-q", "-m", "a")
  w("log.md", "log"); g("add", "log.md"); g("commit", "-q", "-m", "log_after_a")
  t_strict <- strictly_before(last_commit("a.txt"), adding_commit("log.md"))
  w("b.txt", "b"); w("log2.md", "log2"); g("add", "b.txt", "log2.md"); g("commit", "-q", "-m", "same_commit")
  t_same <- !strictly_before(last_commit("b.txt"), adding_commit("log2.md"))
  w("log3.md", "early"); g("add", "log3.md"); g("commit", "-q", "-m", "log3_early")
  w("c.txt", "c"); g("add", "c.txt"); g("commit", "-q", "-m", "c")
  w("log3.md", "touched later"); g("add", "log3.md"); g("commit", "-q", "-m", "log3_touched")
  t_adding <- !strictly_before(last_commit("c.txt"), adding_commit("log3.md")) && strictly_before(last_commit("c.txt"), last_commit("log3.md"))
  t_recommit <- adding_commit("log3.md") != last_commit("log3.md")
  sh <- tempfile("shallow"); system2("git", c("clone", "-q", "--depth", "1", paste0("file://", td), sh), stdout = FALSE, stderr = FALSE)
  GITC <<- c("-C", sh); t_shallow <- inherits(tryCatch(refuse_shallow(), error = function(e) e), "error")
  GITC <<- c("-C", td); t_deep <- !inherits(tryCatch(refuse_shallow(), error = function(e) e), "error")
  res <- c(strict_order_passes = t_strict, same_commit_refused = t_same, adding_commit_used = t_adding,
           recommitted_log_detected = t_recommit, shallow_refused = t_shallow, full_clone_accepted = t_deep)
  print(res); unlink(c(td, sh), recursive = TRUE)
  if (!all(res)) stop("selftest_git FAILED")
  say("selftest_git passed")
}

## ---- annotation helpers -------------------------------------------------------------------------
sgd <- function() {      # SGD_features.tab: 3 qualifier, 4 ORF name, 5 standard name, 6 aliases (|), 9 chr, 10/11 coords, 12 strand
  x <- utils::read.delim(file.path(D, "SGD_features.tab"), header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  x <- x[grepl(ORF_RE, x$V4), ]
  data.frame(orf = x$V4, name = toupper(x$V5), alias = toupper(x$V6), chr = x$V9, dubious = grepl("Dubious", x$V3),
             lo = pmin(as.numeric(x$V10), as.numeric(x$V11)), hi = pmax(as.numeric(x$V10), as.numeric(x$V11)),
             strand = x$V12, stringsAsFactors = FALSE)
}
name_to_orf <- function(s, log_file = NULL) {   # standard name, then systematic, then alias (unique aliases only)
  m <- stats::setNames(s$orf, s$name); m <- m[nzchar(names(m))]
  al <- do.call(rbind, lapply(seq_len(nrow(s)), function(k) if (nzchar(s$alias[k]))
    data.frame(a = strsplit(s$alias[k], "|", fixed = TRUE)[[1]], orf = s$orf[k]) else NULL))
  amb <- unique(al$a[duplicated(al$a)])
  if (!is.null(log_file)) writeLines(sort(amb), log_file)
  al <- al[!al$a %in% amb, ]
  function(x) { x <- toupper(x); out <- unname(m[x]); out[is.na(out) & x %in% s$orf] <- x[is.na(out) & x %in% s$orf]
    miss <- is.na(out); out[miss] <- al$orf[match(x[miss], al$a)]; out }
}
neighbours <- function(s, orf, win = 500) {
  k <- match(orf, s$orf); if (is.na(k)) return(character(0))
  same <- s$chr == s$chr[k]; gap <- pmax(0, pmax(s$lo, s$lo[k]) - pmin(s$hi, s$hi[k]))
  s$orf[same & gap <= win]
}
divergent_pairs <- function(s) {   # adjacent non-dubious genes on opposite strands with 5' ends facing
  s <- s[!s$dubious, ]; out <- list()
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
dye_check <- function(meta, tag) {   # DEVIATION 1: part 1 stops the stage; part 2 returns the balanced groups only
  r1 <- grepl(REF_RE, meta$src1, ignore.case = TRUE); r2 <- grepl(REF_RE, meta$src2, ignore.case = TRUE)
  bad_hyb <- meta$gsm[xor(r1, r2) == FALSE]
  if (length(bad_hyb)) { writeLines(bad_hyb, file.path(SEL, paste0("dye_check_part1_failures_", tag, ".txt")))
    stop("structural dye check part 1 failed: ", length(bad_hyb), " hybridizations without exactly one reference channel") }
  key <- group_key(meta$title)
  bal <- tapply(r1, key, function(v) length(v) == 2 && sum(v) == 1)
  utils::write.csv(data.frame(group = names(bal), n_hyb = as.integer(table(key)[names(bal)]), balanced = as.logical(bal)),
                   file.path(SEL, paste0("dye_balance_", tag, ".csv")), row.names = FALSE)
  say(sprintf("dye check %s: part 1 ok (%d hybridizations); balanced groups %d of %d", tag, nrow(meta), sum(bal), length(bal)))
  names(bal)[bal]
}
GSM_CACHE <- new.env()
gsm_table <- function(gsm) {
  if (!is.null(GSM_CACHE[[gsm]])) return(GSM_CACHE[[gsm]])
  f <- file.path(D, "gsm", paste0(gsm, ".txt"))
  if (!file.exists(f)) { dir.create(dirname(f), showWarnings = FALSE)
    fetch(sprintf("%s?acc=%s&targ=self&form=text&view=data", ACC, gsm), f); Sys.sleep(0.4) }
  L <- readLines(f); a <- which(L == "!sample_table_begin"); b <- which(L == "!sample_table_end")
  if (!length(a)) stop("no sample table in ", f)
  tb <- utils::read.delim(text = L[(a + 1L):(if (length(b)) b - 1L else length(L))], stringsAsFactors = FALSE)
  v <- stats::setNames(as.numeric(tb$VALUE), as.character(tb$ID_REF)); GSM_CACHE[[gsm]] <- v; v
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
RULES <- c("asis", "flip_test_ch1", "flip_ref_ch1")
subseries <- function(title) ifelse(grepl("^\\[", title), sub("^(\\[[^]]*\\]).*$", "\\1", title), "none")   # DEVIATION 2
rule_of <- function(meta, rule) {     # `rule`: one rule, or a named vector sub-series -> rule (missing -> "asis")
  if (length(rule) == 1L && is.null(names(rule))) return(rep(rule, nrow(meta)))
  r <- unname(rule[subseries(meta$title)]); r[is.na(r)] <- "asis"; r }
orient <- function(v, rule, ref_ch1) if ((rule == "flip_test_ch1" && !ref_ch1) || (rule == "flip_ref_ch1" && ref_ch1)) -v else v
rebuild <- function(meta, rule, pf) {   # mean of the oriented hybridizations per group, then per ORF
  key <- group_key(meta$title); r1 <- grepl(REF_RE, meta$src1, ignore.case = TRUE); rr <- rule_of(meta, rule)
  vals <- lapply(seq_len(nrow(meta)), function(k) orient(gsm_table(meta$gsm[k]), rr[k], r1[k]))
  ids <- sort(Reduce(union, lapply(vals, names)))
  V <- vapply(vals, function(v) v[ids], numeric(length(ids))); rownames(V) <- ids
  V <- vapply(split(seq_len(ncol(V)), key), function(ix) rowMeans(V[, ix, drop = FALSE], na.rm = TRUE), numeric(nrow(V)))
  rownames(V) <- ids; by_orf(V, pf$orf[match(rownames(V), pf$id)])
}
finish_matrix <- function(X, mask) {
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
  mkd(); pf <- platform(); s <- sgd(); to_orf <- name_to_orf(s, file.path(SEL, "ambiguous_aliases.txt"))
  build <- function(g) {
    r <- read_final(g); X <- by_orf(r$M, pf$orf[match(rownames(r$M), pf$scab)])
    del <- to_orf(colnames(X)); res <- rep("ok", ncol(X)); ownM <- rep(NA_real_, ncol(X)); q01 <- rep(NA_real_, ncol(X))
    res[is.na(del) | !del %in% rownames(X)] <- "unresolved"
    for (j in which(res == "ok")) { v <- X[, j]; ownM[j] <- v[del[j]]; q01[j] <- stats::quantile(v, 0.01, na.rm = TRUE)
      if (!is.finite(ownM[j]) || ownM[j] > -1 || ownM[j] > q01[j]) res[j] <- "own_M_not_verified" }   # SF6
    log <- data.frame(series = g, column = r$label, deleted_gene = colnames(X), deleted_orf = del, own_M = ownM,
                      column_q01 = q01, status = res)
    ok <- res == "ok"; X <- X[, ok, drop = FALSE]; del <- del[ok]
    mask <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X)); nmask <- integer(ncol(X))
    for (j in seq_len(ncol(X))) { nb <- intersect(neighbours(s, del[j], 500), rownames(X)); mask[nb, j] <- TRUE; nmask[j] <- length(nb) }
    f <- finish_matrix(X, mask)
    list(X = f$X, log = log, dropped_missing = f$dropped_missing, nmask = nmask, del = del)
  }
  b <- build("GSE42527")
  wr_mat(b$X, file.path(D, "deleteome_responsive_700.tsv"))
  utils::write.csv(b$log[b$log$status != "ok", ], file.path(SEL, "dropped_columns.csv"), row.names = FALSE)
  utils::write.csv(cbind(b$log[b$log$status == "ok", ], n_masked = b$nmask), file.path(SEL, "masking_log.csv"), row.names = FALSE)
  say(sprintf("responsive matrix: %d ORFs x %d mutants; dropped columns: %d unresolved, %d own-M not verified; %d ORFs dropped (>5%% missing)",
              nrow(b$X), ncol(b$X), sum(b$log$status == "unresolved"), sum(b$log$status == "own_M_not_verified"), b$dropped_missing))
  a <- build("GSE42526"); common <- intersect(rownames(b$X), rownames(a$X))           # descriptive sensitivity set
  wr_mat(cbind(b$X[common, ], a$X[common, ]), file.path(D, "deleteome_all_1484.tsv"))
  utils::write.csv(divergent_pairs(s), file.path(SEL, "divergent_pairs.csv"), row.names = FALSE)
  writeLines(sprintf("%s  %s", sha(file.path(D, "deleteome_responsive_700.tsv")), "deleteome_responsive_700.tsv"),
             file.path(SEL, "deleteome_matrix.sha256"))
  record_hashes()
  say("commit ", SEL, "/{deleteome_matrix.sha256, dropped_columns.csv, masking_log.csv, divergent_pairs.csv, ambiguous_aliases.txt}")
}

## ---- stage: orient_split (DEVIATION 2) -----------------------------------------------------------
if (stage == "orient_split") {
  mu <- series_meta("GSE42527"); bal <- dye_check(mu, "GSE42527")
  pub <- read_final("GSE42527"); pool <- sort(intersect(unique(bal), tolower(colnames(pub$M))))
  set.seed(20261106); fit <- pool[sort(sample.int(length(pool), floor(length(pool) / 2)))]
  out <- data.frame(group = pool, set = ifelse(pool %in% fit, "fitting", "held_out"))
  utils::write.csv(out, file.path(SEL, "orientation_split.csv"), row.names = FALSE)
  say(sprintf("orientation split: %d balanced mutants -> %d fitting, %d held-out; commit and push before stage check",
              length(pool), sum(out$set == "fitting"), sum(out$set == "held_out")))
}

## ---- stage: check (DEVIATION 2: per-sub-series rule, held-out floor) -------------------------
if (stage == "check") {
  sp_f <- file.path(SEL, "orientation_split.csv")
  if (!clean(sp_f) || !on_upstream(last_commit(sp_f))) stop("gate: orientation_split.csv must be committed and pushed before fitting")
  sp <- utils::read.csv(sp_f); pf <- platform(); mu <- series_meta("GSE42527"); bal <- dye_check(mu, "GSE42527")
  key <- group_key(mu$title); mu$key <- key; mu$ss <- subseries(mu$title)
  pub <- read_final("GSE42527"); Mo <- by_orf(pub$M, pf$orf[match(rownames(pub$M), pf$scab)])
  r_of <- function(groups, rule) { meta <- mu[mu$key %in% groups, ]; if (!nrow(meta)) return(NULL); R <- rebuild(meta, rule, pf)
    do.call(rbind, lapply(colnames(R), function(g) { j <- which(tolower(colnames(Mo)) == tolower(g))[1]
      o <- intersect(rownames(R), rownames(Mo)); a <- R[o, g]; b <- Mo[o, j]
      data.frame(group = g, r = stats::cor(a, b, use = "complete.obs"), med_abs = stats::median(abs(a - b), na.rm = TRUE)) })) }
  fit <- sp$group[sp$set == "fitting"]; held <- sp$group[sp$set == "held_out"]
  ss_of_group <- tapply(mu$ss, mu$key, function(v) if (length(unique(v)) == 1L) v[1] else "mixed")
  rules <- list(); fitrows <- list()
  for (ss in sort(unique(mu$ss[mu$key %in% c(fit, held)]))) {
    g <- fit[ss_of_group[fit] == ss]
    if (length(g) >= 5) {
      med <- vapply(RULES, function(rl) stats::median(r_of(g, rl)$r), 0); rules[[ss]] <- names(which.max(med))
    } else { med <- stats::setNames(rep(NA_real_, 3), RULES); rules[[ss]] <- "asis" }
    fitrows[[ss]] <- data.frame(subseries = ss, n_fit = length(g), rule = rules[[ss]], t(med), default_used = length(g) < 5) }
  rule_map <- unlist(rules)
  utils::write.csv(do.call(rbind, fitrows), file.path(SEL, "orientation_rules.csv"), row.names = FALSE)
  ho <- r_of(held, rule_map); ho$subseries <- ss_of_group[ho$group]; ho$set <- "held_out"
  utils::write.csv(ho, file.path(SEL, "implementation_check.csv"), row.names = FALSE)
  medr <- stats::median(ho$r, na.rm = TRUE)
  status <- if (all(ho$r >= 0.99 & ho$med_abs <= 0.02)) "published_M_reproduced" else if (medr >= 0.8) "NOT_reproduced_run_stage_rebuild" else "F_NOT_COMPUTABLE_verdict_1"
  writeLines(c("per_subseries:orientation_rules.csv", status, sprintf("heldout_median_r=%.4f n=%d", medr, nrow(ho)), env_info()),
             file.path(SEL, "orientation_rule.txt")); record_hashes()
  print(do.call(rbind, fitrows)); print(tapply(ho$r, ho$subseries, stats::median))
  say(sprintf("held-out median r %.3f over %d mutants; status %s", medr, nrow(ho), status))
}

## ---- stage: rebuild (only after a failed, floor-passing check) ----------------------------------
if (stage == "rebuild") {
  rr <- readLines(file.path(SEL, "orientation_rule.txt"))
  if (rr[2] != "NOT_reproduced_run_stage_rebuild") stop("rebuild is used only after a failed but floor-passing check (status: ", rr[2], ")")
  rmap <- utils::read.csv(file.path(SEL, "orientation_rules.csv")); rmap <- stats::setNames(rmap$rule, rmap$subseries)
  pf <- platform(); mu <- series_meta("GSE42527"); bal <- dye_check(mu, "GSE42527"); key <- group_key(mu$title)
  X <- rd_mat(file.path(D, "deleteome_responsive_700.tsv"))
  cand <- sort(intersect(unique(key[key %in% bal]), tolower(colnames(X)))); Nw <- length(unique(group_key(series_meta("GSE45115")$title)))
  set.seed(20261105); pick <- sort(cand[sample.int(length(cand), min(Nw, length(cand)))])
  meta <- mu[key %in% pick, ]
  R <- rebuild(meta, rmap, pf)
  m <- utils::read.csv(file.path(SEL, "masking_log.csv")); s <- sgd()
  R <- R[intersect(rownames(X), rownames(R)), , drop = FALSE]
  mask <- matrix(FALSE, nrow(R), ncol(R), dimnames = dimnames(R))
  for (j in seq_len(ncol(R))) { d <- m$deleted_orf[match(tolower(colnames(R)[j]), tolower(m$deleted_gene))]
    if (!is.na(d)) mask[intersect(neighbours(s, d, 500), rownames(R)), j] <- TRUE }
  f <- finish_matrix(R, mask); wr_mat(f$X, file.path(D, "mutants_rebuilt_376.tsv"))
  writeLines(pick, file.path(SEL, "rebuilt_mutants.txt")); record_hashes()
  say(sprintf("rebuilt %d mutants (per-sub-series rules) for gate F", ncol(f$X)))
}

## ---- stage: wt (DEVIATION 2: "asis" + label-free pair check r > 0) ----------------------------
if (stage == "wt") {
  rf <- file.path(SEL, "orientation_rule.txt"); if (!file.exists(rf)) stop("run stage check first")
  rr <- readLines(rf); if (rr[2] == "F_NOT_COMPUTABLE_verdict_1") stop("orientation floor failed: gate F not computable (verdict 1)")
  pf <- platform(); meta <- series_meta("GSE45115"); bal <- dye_check(meta, "GSE45115")
  meta <- meta[group_key(meta$title) %in% bal, ]; key <- group_key(meta$title)
  one <- function(k) by_orf(matrix(gsm_table(meta$gsm[k]), ncol = 1, dimnames = list(names(gsm_table(meta$gsm[k])), "v")),
                            pf$orf[match(names(gsm_table(meta$gsm[k])), pf$id)])
  pc <- do.call(rbind, lapply(split(seq_len(nrow(meta)), key), function(ix) { a <- one(ix[1]); b <- one(ix[2])
    o <- intersect(rownames(a), rownames(b)); data.frame(group = key[ix[1]], r_pair = stats::cor(a[o, 1], b[o, 1], use = "complete.obs")) }))
  pc$pass <- is.finite(pc$r_pair) & pc$r_pair > 0
  utils::write.csv(pc, file.path(SEL, "wt_pair_check.csv"), row.names = FALSE)
  say(sprintf("WT pair check (r > 0): %d of %d balanced groups pass (%.1f%%)", sum(pc$pass), nrow(pc), 100 * mean(pc$pass)))
  meta <- meta[key %in% pc$group[pc$pass], ]
  W <- rebuild(meta, "asis", pf)
  W <- W[intersect(universe_rows(), rownames(W)), , drop = FALSE]
  f <- finish_matrix(W, matrix(FALSE, nrow(W), ncol(W)))
  wr_mat(f$X, file.path(D, "wt_controls_GSE45115.tsv")); record_hashes()
  say(sprintf("wild-type matrix: %d ORFs x %d cultures", nrow(f$X), ncol(f$X)))
}

## ---- stage: splits (halves, clusters, experiment-jackknife replicates; SF5) ---------------------
if (stage == "splits") {
  X <- rd_mat(file.path(D, "deleteome_responsive_700.tsv")); mut <- colnames(X); n <- length(mut)
  set.seed(20261101)
  sp <- vapply(1:10, function(s) { h <- integer(n); h[sample.int(n, floor(n / 2))] <- 1L; h + 1L }, integer(n))
  colnames(sp) <- sprintf("split_%02d", 1:10)
  utils::write.csv(data.frame(mutant = mut, sp), file.path(SEL, "splits.csv"), row.names = FALSE)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - stats::cor(X)), "complete"), h = 0.5)
  utils::write.csv(data.frame(mutant = mut, cluster = unname(cl)), file.path(SEL, "mutant_clusters.csv"), row.names = FALSE)
  sz <- table(cl); K <- length(sz); fallback <- K < 20 || max(sz) / n > 0.25
  ids <- sort(unique(unname(cl)))
  reps <- vapply(1:200, function(b) { set.seed(20262000 + b)
    keep <- if (!fallback) cl %in% ids[sample.int(K, floor(K / 2))] else {
      k <- logical(n); for (c in ids) { m <- which(cl == c); h <- if (stats::runif(1) < 0.5) floor(length(m) / 2) else ceiling(length(m) / 2)
        if (h > 0) k[m[sample.int(length(m), h)]] <- TRUE }; k }
    as.integer(keep) }, integer(n))
  colnames(reps) <- sprintf("rep_%03d", 1:200)
  utils::write.csv(data.frame(mutant = mut, reps), file.path(SEL, "experiment_jackknife_replicates.csv"), row.names = FALSE)
  utils::write.csv(data.frame(n_clusters = K, singletons = sum(sz == 1), largest = max(sz), largest_share = max(sz) / n,
    stratified_fallback = fallback, rep_N_min = min(colSums(reps)), rep_N_max = max(colSums(reps)),
    t(stats::quantile(as.integer(sz), c(.5, .9, .99))), t(env_info())),
    file.path(SEL, "mutant_clusters_summary.csv"), row.names = FALSE)
  say(sprintf("10 halves; %d clusters (largest %d = %.0f%%); %s; replicate N %d-%d -> commit", K, max(sz), 100 * max(sz) / n,
              if (fallback) "STRATIFIED FALLBACK" else "cluster delete-half", min(colSums(reps)), max(colSums(reps))))
}

## ---- stage: feasibility (gate F, sanity check) --------------------------------------------------
if (stage == "feasibility") {
  suppressPackageStartupMessages(library(clr))
  rr <- readLines(file.path(SEL, "orientation_rule.txt"))
  if (rr[2] == "F_NOT_COMPUTABLE_verdict_1") {
    utils::write.csv(data.frame(pass = FALSE, reason = "orientation floor failed: F not computable"), file.path(SEL, "feasibility_gate.csv"), row.names = FALSE)
    stop("gate F not computable (verdict 1)") }
  Xm <- if (rr[2] == "published_M_reproduced") rd_mat(file.path(D, "deleteome_responsive_700.tsv")) else {
    f <- file.path(D, "mutants_rebuilt_376.tsv"); if (!file.exists(f)) stop("check failed: run stage rebuild first"); rd_mat(f) }
  W <- rd_mat(file.path(D, "wt_controls_GSE45115.tsv"))
  g <- intersect(rownames(Xm), rownames(W)); Xm <- Xm[g, ]; W <- W[g, ]
  fitC <- function(M) { rng <- apply(M, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) M[rng <= 0, ] <- M[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(M), sd = 1e-9)
    f <- ClrAttention$new(M)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
    f$calibrate(method = "normal", combine = "stouffer"); C <- f$clr_scores; C[lower.tri(C, diag = TRUE)] <- NA; C }
  top <- function(C) { v <- which(!is.na(C)); k <- ceiling(0.01 * length(v)); v[order(C[v], decreasing = TRUE)[seq_len(k)]] }
  repro <- function(M, cols) { h <- cols[sample.int(length(cols), floor(length(cols) / 2))]
    a <- top(fitC(M[, h])); b <- top(fitC(M[, setdiff(cols, h)])); (length(intersect(a, b)) / length(a) - 0.01) / 0.99 }
  Nw <- min(ncol(W), ncol(Xm)); set.seed(20261199)
  Rm <- vapply(1:10, function(s) repro(Xm, sample.int(ncol(Xm), Nw)), 0)
  Rw <- vapply(1:10, function(s) repro(W, seq_len(ncol(W))), 0)
  tm <- top(fitC(Xm)); tw <- top(fitC(W)); Q <- length(intersect(tm, tw)) / length(tm)
  lo <- apply(W, 1, stats::quantile, 0.005); hi <- apply(W, 1, stats::quantile, 0.995)
  beyond <- colSums(Xm < lo | Xm > hi)                                   # genes per mutant outside the WT range
  out <- data.frame(mutant_source = if (rr[2] == "published_M_reproduced") "published_M" else "rebuilt_376",
                    N_matched = Nw, genes = length(g), median_Rprime_mut = stats::median(Rm), median_Rprime_wt = stats::median(Rw),
                    iqr_Rprime_mut = stats::IQR(Rm), iqr_Rprime_wt = stats::IQR(Rw), Q_overlap_share = Q,
                    genes_beyond_wt_median = stats::median(beyond), genes_beyond_wt_q10 = stats::quantile(beyond, 0.1),
                    genes_beyond_wt_q90 = stats::quantile(beyond, 0.9), t(env_info()))
  out$pass <- out$median_Rprime_mut > out$median_Rprime_wt && Q <= 0.5
  utils::write.csv(out, file.path(SEL, "feasibility_gate.csv"), row.names = FALSE)
  utils::write.csv(data.frame(draw = 1:10, Rprime_mut = Rm, Rprime_wt = Rw), file.path(SEL, "feasibility_draws.csv"), row.names = FALSE)
  say("gate F: ", if (out$pass) "PASS" else "FAIL (verdict 1; yeast not consumed)",
      sprintf(" (median R' mut %.3f vs wt %.3f; Q %.3f)", out$median_Rprime_mut, out$median_Rprime_wt, Q))
}

## ---- shared truth-stage helpers ------------------------------------------------------------------
neff_of <- function(sizes) { p <- choose(sizes, 2); p <- p[p > 0]; if (!length(p)) return(0); sum(p)^2 / sum(p^2) }
gate_pre_truth <- function() {    # steps 1 and 4-7 committed; universe matrix and logs verified
  hash <- gate_prereg()
  need_clean(c(file.path(SEL, c("excluded_pmids.txt", "feasibility_gate.csv", "splits.csv", "mutant_clusters.csv",
                                "experiment_jackknife_replicates.csv", "R_check.csv", "deleteome_matrix.sha256",
                                "masking_log.csv", "dropped_columns.csv")), SEL_FILES, TEST_SCRIPT, CONVERTER, PREP_SCRIPT))
  if (!isTRUE(as.logical(utils::read.csv(file.path(SEL, "feasibility_gate.csv"))$pass))) stop("gate F failed: verdict 1, stop")
  if (!all(as.logical(utils::read.csv(file.path(SEL, "R_check.csv"))$passed))) stop("synthetic R check not passed: no verdict, fix first")
  want <- strsplit(readLines(file.path(SEL, "deleteome_matrix.sha256"), 1L), "[[:space:]]+")[[1]][1]
  if (!identical(want, sha(file.path(D, "deleteome_responsive_700.tsv")))) stop("gate: expression matrix differs from committed SHA-256")
  pre <- c(file.path(SEL, c("feasibility_gate.csv", "splits.csv", "mutant_clusters.csv", "experiment_jackknife_replicates.csv", "R_check.csv")),
           SEL_FILES, TEST_SCRIPT, CONVERTER)
  for (f in pre) if (!on_upstream(last_commit(f))) stop("gate: the commit of ", f, " is not on the upstream branch (push and git fetch first)")
  list(prereg = hash, test_blob = git1("rev-parse", paste0("HEAD:", TEST_SCRIPT)), conv_blob = git1("rev-parse", paste0("HEAD:", CONVERTER)),
       prep_blob = git1("rev-parse", paste0("HEAD:", PREP_SCRIPT)),
       push_record = paste(vapply(pre, function(f) paste0(basename(f), "=", substr(last_commit(f), 1, 12)), ""), collapse = ";"))
}
gate_after <- function(later, earlier, allow_touch = FALSE) {   # SF2 (review3): the commit ADDING `later` is strictly after
  ca <- adding_commit(later); if (is.na(ca)) stop("gate: ", later, " not committed")   # the last commit of every `earlier`
  if (!allow_touch && last_commit(later) != ca) stop("gate: ", later, " was changed after the commit that added it (log a deviation instead)")
  for (f in earlier) if (!strictly_before(last_commit(f), ca)) stop("gate: ", later, " must be ADDED in a commit strictly after ", f)
  if (!on_upstream(ca)) stop("gate: the commit adding ", later, " is not on the upstream branch (push and git fetch first)")
  invisible(ca)
}
PRE_TRUTH_FILES <- c(file.path(SEL, c("feasibility_gate.csv", "R_check.csv")), SEL_FILES, TEST_SCRIPT, PREP_SCRIPT)

## ---- stage: complex_structure (8a; GATED; SF3 provisional: before any regulon membership) -------
if (stage == "complex_structure") {
  g <- gate_pre_truth(); E <- readLines(file.path(SEL, "excluded_pmids.txt")); E <- E[grepl("^[0-9]+$", E)]
  cnt_f <- file.path(SEL, "complex_structure_counts.csv")
  if (!is.na(opt$cyc2008)) {        # re-run of 8a after a committed first run that required the fallback (review3 note)
    lg <- file.path(SEL, "cyc2008_download_log.md"); need_clean(c(lg, cnt_f))
    if (!grepl("CYC2008 required", utils::read.csv(cnt_f)$fallback[1])) stop("gate: CYC2008 used only after a committed first 8a run that required it")
    gate_after(lg, c(PRE_TRUTH_FILES, cnt_f))
    if (!any(grepl(sha(opt$cyc2008), readLines(lg), fixed = TRUE))) stop("gate: CYC2008 file SHA-256 not listed in ", lg) }
  universe <- universe_rows(); Tdir <- file.path(D, "truth_raw"); dir.create(Tdir, showWarnings = FALSE)
  fetch("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/559292.tsv", file.path(Tdir, "559292.tsv"))
  fetch("https://rest.uniprot.org/uniprotkb/stream?query=organism_id:559292&format=tsv&fields=accession,gene_oln,gene_primary",
        file.path(Tdir, "uniprot_559292.tsv"))
  up <- utils::read.delim(file.path(Tdir, "uniprot_559292.tsv"), stringsAsFactors = FALSE, quote = "")
  orf_of <- lapply(strsplit(up[[2]], "[ ;]+"), function(v) unique(grep(ORF_RE, v, value = TRUE))); names(orf_of) <- up[[1]]
  cp <- utils::read.delim(file.path(Tdir, "559292.tsv"), stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
  xcol <- grep("cross.?ref", names(cp), ignore.case = TRUE)[1]; ecol <- grep("evidence code", names(cp), ignore.case = TRUE)[1]
  cid <- sub("\\..*$", "", cp[[1]]); parts <- lapply(strsplit(cp[[5]], "|", fixed = TRUE), function(v) sub("\\(.*$", "", v)); names(parts) <- cid
  cpm <- lapply(regmatches(cp[[xcol]], gregexpr("pubmed:[0-9]+", cp[[xcol]])), function(v) sub("pubmed:", "", v)); names(cpm) <- cid
  only_E <- vapply(cpm, function(v) length(v) > 0 && all(v %in% E), TRUE)
  flatten <- function(id, seen = character(0)) { if (id %in% seen) return(character(0)); out <- character(0)
    for (p in parts[[id]]) { if (grepl("^CPX-", p)) { if (p %in% names(parts)) out <- c(out, flatten(p, c(seen, id))) }
      else if (grepl("^(CHEBI|URS|EBI|SIGNOR|CHEMBL)", p)) next else out <- c(out, sub("-[0-9]+$|-PRO_[0-9]+$", "", p)) }; out }
  acc_of <- lapply(names(parts), flatten); names(acc_of) <- names(parts)
  sets <- lapply(acc_of, function(a) sort(unique(unlist(orf_of[a[a %in% names(orf_of)]]))))
  keep_ids <- sort(names(sets)[!duplicated(names(sets)) & !only_E[names(sets)]])
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
    evaluable_P1 = fin$P1_positive_pairs >= 200 && fin$P1_complexes_with_pos >= 30,
    evaluable_S1 = fin$S1_positive_pairs >= 200 && fin$S1_complexes_with_pos >= 30,
    power_gate_pass = fin$P1_neff >= 100, prereg_commit = g$prereg, test_script_blob = g$test_blob, converter_blob = g$conv_blob,
    prep_blob = g$prep_blob, push_record = g$push_record, checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), t(env_info()))
  out$reciprocity_verdict1 <- !(out$power_gate_pass && out$evaluable_P1)
  mem <- do.call(rbind, lapply(names(merged), function(k) data.frame(complex = k, orf = merged[[k]])))
  utils::write.table(mem, file.path(SEL, "yeast_complexes_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.csv(out, file.path(SEL, "complex_structure_counts.csv"), row.names = FALSE)
  utils::write.csv(data.frame(complex = cid, eco = if (!is.na(ecol)) cp[[ecol]] else NA, pmids = vapply(cpm, paste, "", collapse = ";"),
                              only_E = only_E), file.path(SEL, "complex_evidence.csv"), row.names = FALSE)
  writeLines(c(sprintf("%s  yeast_complexes_members.tsv", sha(file.path(SEL, "yeast_complexes_members.tsv"))),
               sprintf("%s  truth_raw/559292.tsv", sha(file.path(Tdir, "559292.tsv"))),
               sprintf("%s  truth_raw/uniprot_559292.tsv", sha(file.path(Tdir, "uniprot_559292.tsv"))),
               if (!is.na(opt$cyc2008)) sprintf("%s  %s", sha(opt$cyc2008), basename(opt$cyc2008))), file.path(SEL, "yeast_complexes.sha256"))
  say("8a power gate: ", if (out$power_gate_pass) "PASS" else "FAIL", "; reciprocity verdict 1: ", out$reciprocity_verdict1,
      sprintf(" (P1 n_eff %.1f)", out$P1_neff), "; commit ", SEL, " BEFORE the YEASTRACT export (8b)")
}

## ---- stage: regulon_structure (8c; GATED; SF1, SF2) ---------------------------------------------
if (stage == "regulon_structure") {
  g <- gate_pre_truth()
  lg <- file.path(SEL, "yeastract_export_log.md")
  need_clean(file.path(SEL, c("complex_structure_counts.csv", "yeast_complexes.sha256", "yeastract_export_log.md")))
  for (f in file.path(SEL, c("complex_structure_counts.csv", "yeast_complexes.sha256")))
    if (!on_upstream(last_commit(f))) stop("gate: 8a commit not on the upstream branch (push and git fetch first)")
  ca <- gate_after(lg, c(PRE_TRUTH_FILES, file.path(SEL, "complex_structure_counts.csv")))
  # converter: committed strictly before the export log, or a committed deviation route (review3 SF2 e)
  dev <- file.path(SEL, "converter_deviation.md")
  if (!strictly_before(last_commit(CONVERTER), ca)) {
    need_clean(dev); cd <- adding_commit(dev)
    if (!strictly_before(ca, cd) || !is_ancestor(last_commit(CONVERTER), cd))
      stop("gate: converter changed after the export without a converter_deviation.md committed after the export log and with the change")
    say("converter deviation route used: ", dev) }
  # raw export files listed in the log with SHA-256; 8c runs the committed converter itself (review3 SF2 f)
  ll <- grep("^raw_export: [0-9a-f]{64} (table|binding_matrix|expression_matrix) ", readLines(lg), value = TRUE)
  if (!length(ll)) stop("gate: ", lg, " lists no raw_export lines")
  cargs <- character(0)
  for (x in ll) { p <- strsplit(x, " ", fixed = TRUE)[[1]]; fpath <- paste(p[-(1:3)], collapse = " ")
    if (!file.exists(fpath) || sha(fpath) != p[2]) stop("gate: raw export ", fpath, " missing or SHA-256 differs from the log")
    cargs <- c(cargs, paste0("--", p[3]), fpath) }
  norm <- file.path(D, "yeastract_normalized.tsv")
  rc <- system2("Rscript", c(CONVERTER, cargs, "--out", norm)); if (rc != 0L) stop("converter failed (rc ", rc, "): log a converter deviation")
  opt$yeastract <- norm
  E <- readLines(file.path(SEL, "excluded_pmids.txt")); E <- E[grepl("^[0-9]+$", E)]
  universe <- universe_rows(); s <- sgd(); to_orf <- name_to_orf(s, file.path(SEL, "ambiguous_aliases_yeastract.txt"))
  yx <- utils::read.delim(opt$yeastract, stringsAsFactors = FALSE, colClasses = "character")
  stopifnot(all(c("tf", "target", "evidence", "pmids", "pmids_unattributed", "value") %in% names(yx)))
  value_table <- as.data.frame(table(yx$value))
  yx <- yx[yx$value != "0" & nzchar(yx$value), ]                          # any non-zero entry (incl. -1) counts
  yx$tf_orf <- to_orf(yx$tf); yx$target_orf <- to_orf(yx$target); yx$key <- paste(yx$tf_orf, yx$target_orf)
  pm <- strsplit(yx$pmids, "[;, ]+"); has_pm <- vapply(pm, function(v) any(grepl("^[0-9]+$", v)), TRUE)
  ex <- yx$evidence == "expression"
  ex_keys <- unique(yx$key[ex]); ex_keys_attr <- unique(yx$key[ex & has_pm])
  share_attr <- if (length(ex_keys)) length(ex_keys_attr) / length(ex_keys) else 0
  mode <- if (share_attr >= 0.95) "refs" else "binding_only"                # SF1 (i)-(ii): computed, never chosen
  if (!identical(mode, opt$yeastract_mode)) stop("computed mode '", mode, "' (attributed share ", round(share_attr, 3),
                                                  ") disagrees with --yeastract_mode '", opt$yeastract_mode, "': stop and log")
  if (mode == "refs") {     # reanalysis screen (c): anywhere in the PubMed XML record
    cand <- sort(setdiff(unique(unlist(pm[ex])), c(E, ""))); cand <- cand[grepl("^[0-9]+$", cand)]; hit <- character(0)
    for (ch in split(cand, ceiling(seq_along(cand) / 150))) {
      f <- tempfile(); fetch(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=xml",
                                     paste(ch, collapse = ",")), f)
      x <- paste(readLines(f, warn = FALSE), collapse = " "); arts <- strsplit(x, "<PubmedArticle>", fixed = TRUE)[[1]][-1]
      for (a in arts) { id <- sub(".*?<PMID[^>]*>([0-9]+)</PMID>.*", "\\1", a)
        if (grepl("deleteome|GSE42526|GSE42527|GSE42528|GSE45115", a, ignore.case = TRUE)) hit <- c(hit, id) }
      Sys.sleep(0.4) }
    E2 <- sort(unique(c(E, hit))); writeLines(sort(hit), file.path(SEL, "excluded_pmids_screen_additions.txt"))
  } else E2 <- E
  bind <- unique(yx$key[yx$evidence == "binding"])
  expr_any <- unique(yx$key[ex])
  expr_clean <- unique(yx$key[ex & vapply(pm, function(v) any(grepl("^[0-9]+$", v) & !v %in% E2), TRUE)])
  primary <- if (mode == "refs") intersect(bind, expr_clean) else bind
  circular <- intersect(bind, expr_any)
  regs_of <- function(keys) { p <- do.call(rbind, strsplit(keys, " ")); p <- p[p[, 1] != "NA" & p[, 2] != "NA", , drop = FALSE]
    r <- lapply(split(p[, 2], p[, 1]), function(v) intersect(unique(v), universe)); r[lengths(r) >= 5 & lengths(r) <= 500] }
  comem_pairs <- function(r) { if (!length(r)) return(0); gg <- sort(unique(unlist(r)))
    M <- Matrix::sparseMatrix(i = match(unlist(r), gg), j = rep(seq_along(r), lengths(r)), x = 1, dims = c(length(gg), length(r)))
    C <- Matrix::tcrossprod(M); (Matrix::nnzero(C) - length(gg)) / 2 }
  regs <- regs_of(primary); regs_circ <- regs_of(circular)
  regulators <- intersect(names(regs), universe)
  edges <- unique(primary[vapply(strsplit(primary, " "), function(v)
    v[1] %in% regulators && v[2] %in% universe && !v[2] %in% regulators, TRUE)])       # SF8 (d): TF-TF outside the universe
  cp_prim <- comem_pairs(regs); cp_circ <- comem_pairs(regs_circ)
  audit <- data.frame(mode = mode, attributed_share = share_attr, E_size = length(E2), screen_additions = length(E2) - length(E),
    assoc_binding = length(bind), assoc_expression_any = length(expr_any), assoc_binding_and_any_expression = length(circular),
    primary_assoc = length(primary),
    share_circular_removed = if (mode == "refs" && length(circular)) 1 - length(intersect(circular, primary)) / length(circular) else NA,
    share_expression_only_E = if (mode == "refs" && length(expr_any)) 1 - length(expr_clean) / length(expr_any) else NA,
    share_comem_pos_removed = if (mode == "refs" && cp_circ > 0) 1 - cp_prim / cp_circ else NA,
    share_binding_with_expression = if (length(bind)) length(intersect(bind, expr_any)) / length(bind) else NA,
    comem_pos_primary = cp_prim, comem_pos_circular = cp_circ,
    tf_unmapped = sum(is.na(unique(yx$tf_orf))), target_unmapped = sum(is.na(unique(yx$target_orf))))
  out <- data.frame(P2_regulons = length(regs), P2_genes = length(unique(unlist(regs))), P2_positive_pairs = cp_prim,
    P2_neff = neff_of(lengths(regs)), S6_regulators = length(regulators), S6_positive_edges = length(edges),
    evaluable_P2 = cp_prim >= 200, evaluable_S6 = length(edges) >= 200, mode = mode,
    prereg_commit = g$prereg, test_script_blob = g$test_blob, converter_blob = g$conv_blob, prep_blob = g$prep_blob,
    export_log_adding_commit = ca, push_record = g$push_record, checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), t(env_info()))
  rmem <- do.call(rbind, lapply(names(regs), function(k) data.frame(regulator = k, orf = regs[[k]])))
  utils::write.table(rmem, file.path(SEL, "yeast_regulons_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(do.call(rbind, strsplit(edges, " ")), file.path(SEL, "yeast_tf_edges.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = c("tf", "target"))
  utils::write.csv(out, file.path(SEL, "regulon_structure_counts.csv"), row.names = FALSE)
  utils::write.csv(audit, file.path(SEL, "contamination_audit.csv"), row.names = FALSE)
  utils::write.csv(value_table, file.path(SEL, "yeastract_value_table.csv"), row.names = FALSE)
  writeLines(c(sprintf("%s  yeast_regulons_members.tsv", sha(file.path(SEL, "yeast_regulons_members.tsv"))),
               sprintf("%s  yeast_tf_edges.tsv", sha(file.path(SEL, "yeast_tf_edges.tsv"))),
               sprintf("%s  %s", sha(opt$yeastract), basename(opt$yeastract)),
               sprintf("%s  %s", sha(CONVERTER), CONVERTER), sprintf("%s  %s", sha(PREP_SCRIPT), PREP_SCRIPT)),
             file.path(SEL, "yeast_regulons.sha256"))
  say("8c: mode ", mode, sprintf(" (attributed share %.3f); P2 evaluable %s, n_eff %.1f; S6 evaluable %s", share_attr,
      out$evaluable_P2, out$P2_neff, out$evaluable_S6), "; commit ", SEL, " before any scoring")
}
