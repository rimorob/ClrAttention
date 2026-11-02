#!/usr/bin/env Rscript
# K562 replication stratum of the preregistered scFM comparison (docs/preregistration_scfm.md, change log E14 and the K562 freeze entry).
# Descriptive only: no p-value, no NI claim. Five subsamples of 1,024 non-targeting control cells (results/scfm_k562/subsamples.tsv), the
# 1,244-gene set of results/scfm_k562/geneset.tsv.
#
# Scoring (as for BEELINE, E14): directed TF-node AUPRC ratio and tie-aware EPR by score_truth() of analysis/beeline_baselines.R (evaluated
# from the file), on TF -> target edge lists of DoRothEA A-B and CollecTRI, and knockdown gene -> module-gene edges of the knockdown-response
# modules (response_sets() of analysis/human_truths.R, defaults). Self-pairs 0, non-finite values 0. Per truth and subsample
# d = log2(AUPRC ratio of the reference / AUPRC ratio of the candidate) on the candidate's G-intersection; reported: mean over the five
# subsamples; 95% percentile bootstrap over regulators (B = 1,000; one regulator resample per truth, shared by the five subsamples and all methods).
# Comparators |r|, CLR, CLR-NS are recomputed on every G-intersection from log1p(count / cell total * 13,714) of the same cells.
#
# Truth files are opened only in the stage "truths" below; with --truth synthetic (smoke test) they are never read.
# Usage (repo root):
#   Rscript analysis/scfm_score_k562.R --out results/scfm_score_k562 [--truth real|synthetic] [--boot 1000] [--workers 6]
args0 <- commandArgs(trailingOnly = TRUE)
sopt <- list(out = "results/scfm_score_k562", truth = "real", exports = "~/scfm_models/step4_k562", boot = 1000L, workers = 6L, seed = 20261021L)
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(sopt))
  sopt[[k]] <- if (k %in% c("boot", "workers", "seed")) as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
sopt$exports <- path.expand(sopt$exports)

src <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src)[1]
args <- c("--out", sopt$out, "--threads", "6")
for (e in parse(text = src[seq_len(cut_at - 1L)])) {
  if (identical(deparse(e)[1], "args <- commandArgs(trailingOnly = TRUE)")) next
  eval(e, globalenv())
}
bb <- parse("analysis/beeline_baselines.R")
for (e in bb) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] == "score_truth") eval(e, globalenv())
stopifnot(exists("score_truth"))
source("analysis/human_truths.R")                 # response_sets() (knockdown-response modules)
say(sprintf("scfm_score_k562: truth=%s out=%s boot=%d workers=%d", sopt$truth, sopt$out, sopt$boot, sopt$workers))

H5 <- "data/perturbseq/K562_essential_raw_singlecell_01.h5ad"
gs <- utils::read.delim("results/scfm_k562/geneset.tsv", stringsAsFactors = FALSE)
gs$is_tf_part <- gs$is_tf_part %in% c("True", "TRUE", TRUE); gs$in_variable_part <- gs$in_variable_part %in% c("True", "TRUE", TRUE)
sub <- utils::read.delim("results/scfm_k562/subsamples.tsv", stringsAsFactors = FALSE)
metatxt <- readLines("results/scfm_k562/geneset_meta.json")
MED <- as.numeric(sub("[^0-9.]*([0-9.]+).*", "\\1", grep("median_control_umi", metatxt, value = TRUE)))
stopifnot(is.finite(MED), nrow(gs) == 1244L)
cat_of <- function(grp, field) { codes <- rhdf5::h5read(H5, paste0(grp, "/", field)); cats <- rhdf5::h5read(H5, paste0(grp, "/__categories/", field))
  as.character(cats)[codes + 1L] }
vg <- cat_of("var", "gene_name"); stopifnot(!anyDuplicated(vg[match(gs$symbol, vg)]), all(gs$symbol %in% vg))
gi <- match(gs$symbol, vg)
bar <- as.character(rhdf5::h5read(H5, "obs/cell_barcode"))
SUBS <- 1:5
expr <- list()                                      # log-normalised expression of the five subsamples, gene-set genes x 1024 cells
for (k in SUBS) {
  rows <- sub$file_row[sub$subsample == k]
  C <- rhdf5::h5read(H5, "X", index = list(NULL, rows + 1L)); storage.mode(C) <- "double"      # file_row is 0-based
  stopifnot(identical(bar[rows + 1L], sub$cell_barcode[sub$subsample == k]))
  L <- log1p(t(t(C) / colSums(C)) * MED)[gi, , drop = FALSE]; rownames(L) <- gs$symbol
  expr[[k]] <- L; rm(C)
}
say(sprintf("expression: %d subsamples, %d genes x %d cells, median UMI scale %.0f", length(expr), nrow(expr[[1]]), ncol(expr[[1]]), MED))

rd_mat <- function(src, case, readout) {
  f <- file.path(sopt$exports, src, sprintf("%s.%s.f32", case, readout)); if (!file.exists(f)) return(NULL)
  g <- readLines(file.path(sopt$exports, src, sprintf("%s.genes.txt", case))); G <- length(g)
  m <- matrix(readBin(f, "numeric", n = G * G, size = 4L), G, G); dimnames(m) <- list(g, g); stopifnot(all(is.finite(m))); m
}
rd_genes <- function(src, case) readLines(file.path(sopt$exports, src, sprintf("%s.genes.txt", case)))
orient <- function(A, o) switch(o, T = t(A), A = A, max = pmax(A, t(A)), sym = A)
clean <- function(s) { s[!is.finite(s)] <- 0; diag(s) <- 0; s }
cs <- function(k) paste0("hK562_s", k)
gA <- rd_genes("scprint_medium_A", cs(1)); gS <- rd_genes("scgpt_A", cs(1)); gAS <- gA[gA %in% gS]
GS <- list(full = gs$symbol, A = gA, S = gS, AS = gAS)
for (k in SUBS) stopifnot(identical(rd_genes("scprint_medium_A", cs(k)), gA), identical(rd_genes("scgpt_A", cs(k)), gS))
say(sprintf("gene sets: full %d, A (scPRINT) %d, S (scGPT) %d, AS %d", length(GS$full), length(gA), length(gS), length(gAS)))

## ---- the networks (all five subsamples), per gene-set key --------------------------------------------------------------------
nets <- list()      # nets[[key]][[label]][[k]]: cleaned G x G score matrix
add_net <- function(key, label, k, S) {
  if (is.null(nets[[key]])) nets[[key]] <<- list(); if (is.null(nets[[key]][[label]])) nets[[key]][[label]] <<- list()
  nets[[key]][[label]][[k]] <<- clean(S) }
for (k in SUBS) {
  for (key in names(GS)) {
    g <- GS[[key]]; Xg <- expr[[k]][g, , drop = FALSE]
    add_net(key, "cmp|abs_pearson|sym", k, abs(stats::cor(t(Xg))))
    f <- clr_fit(Xg); add_net(key, "cmp|clr|sym", k, f$clr_scores); f$release(FALSE)
    add_net(key, "cmp|clr_ns|sym", k, ns_clr(Xg))
  }
  A <- rd_mat("scprint_medium_A", cs(k), "R0"); add_net("A", "scprint_medium_A|R0|T", k, t(A)); add_net("A", "scprint_medium_A|R0|A", k, A)
  add_net("A", "scprint_medium_A|R0|max", k, pmax(A, t(A)))
  ii <- match(gAS, gA); add_net("AS", "scprint_medium_A|R0|T", k, t(A[ii, ii]))
  for (rr in c("R0p", "R2")) add_net("A", sprintf("scprint_medium_A|%s|T", rr), k, t(rd_mat("scprint_medium_A", cs(k), rr)))
  add_net("A", "r1_A|R1|T", k, t(rd_mat("r1_A", cs(k), "R1")))
  for (sd in 1001:1003) add_net("A", sprintf("ctrl_perm_s%d|R0|T", sd), k, t(rd_mat(paste0("ctrl_perm_s", sd), cs(k), "R0")))
  for (v in c("rand_esm", "rand_emb")) for (sd in 1:3) add_net("A", sprintf("ctrl_%s_s%d|R0|T", v, sd), k, t(rd_mat(sprintf("ctrl_%s_s%d", v, sd), cs(k), "R0")))
  for (s in c("scprint_large_A", "scprint_small_A")) add_net("A", sprintf("%s|R0|T", s), k, t(rd_mat(s, cs(k), "R0")))
  jj <- match(gAS, gS)
  for (rr in c("attn_last_mean", "tok_cos", "hid_cos", "tok_cos_ln")) {
    M <- rd_mat("scgpt_A", cs(k), rr); o <- if (rr == "attn_last_mean") "T" else "sym"
    add_net("S", sprintf("scgpt_A|%s|%s", rr, o), k, orient(M, o)); add_net("AS", sprintf("scgpt_A|%s|%s", rr, o), k, orient(M[jj, jj], o))
  }
  say(sprintf("subsample %d: networks built", k))
}

## ---- truths (opened here) -----------------------------------------------------------------------------------------------------
if (sopt$truth == "real") {
  dor <- utils::read.delim("data/human_truth/dorothea.tsv", stringsAsFactors = FALSE)
  lv <- vapply(strsplit(dor$dorothea_level, ";"), function(x) any(x %in% c("A", "B")), NA); dor <- dor[lv, ]
  col <- utils::read.delim("data/human_truth/collectri.tsv", stringsAsFactors = FALSE)
  rs <- response_sets(gs$symbol)
  kd <- do.call(rbind, lapply(names(rs$members), function(nm) data.frame(Gene1 = strsplit(nm, "_")[[1]][2], Gene2 = gs$symbol[rs$members[[nm]]])))
  TRUTHS <- list(DoRothEA_AB = data.frame(Gene1 = dor$source_genesymbol, Gene2 = dor$target_genesymbol),
                 CollecTRI = data.frame(Gene1 = col$source_genesymbol, Gene2 = col$target_genesymbol), KD_response = kd)
  say(sprintf("truths (real): %s", paste(sprintf("%s %d edges", names(TRUTHS), vapply(TRUTHS, nrow, 1L)), collapse = "; ")))
} else {            # smoke test: planted from the networks of subsample 1 (never from a truth file)
  set.seed(7L); G <- length(gA); tfi <- sort(sample(match(intersect(gA, gs$symbol[gs$is_tf_part]), gA), 60L))
  plant <- function(sig) do.call(rbind, lapply(tfi, function(r) { s <- rank(sig[r, ]) / G + rnorm(G, sd = 0.3); s[r] <- -Inf
    data.frame(Gene1 = gA[r], Gene2 = gA[order(s, decreasing = TRUE)[1:25]]) }))
  sR0 <- nets$A[["scprint_medium_A|R0|T"]][[1]]; sNS <- nets$A[["cmp|clr_ns|sym"]][[1]]
  TRUTHS <- list(DoRothEA_AB = plant(sNS), CollecTRI = plant(sR0), KD_response = plant(sR0 / stats::sd(sR0) + sNS / stats::sd(sNS)))
  say("truths (SYNTHETIC, smoke test)")
}

## ---- score everything ------------------------------------------------------------------------------------------------------------
rows <- list()
for (key in names(nets)) for (lab in names(nets[[key]])) for (k in SUBS) {
  p <- strsplit(lab, "|", fixed = TRUE)[[1]]; S <- nets[[key]][[lab]][[k]]; g <- GS[[key]]
  for (tn in names(TRUTHS)) { r <- score_truth(S, g, NULL, TRUTHS[[tn]]); if (is.null(r)) next
    rows[[length(rows) + 1L]] <- data.frame(subsample = k, truth = tn, gset = key, G = length(g), source = p[1], readout = p[2], orient = p[3], t(r)) }
}
res <- do.call(rbind, rows); utils::write.csv(res, file.path(sopt$out, "metrics_long.csv"), row.names = FALSE)
say(sprintf("scored %d rows", nrow(res)))

## ---- bootstrap over regulators for the preregistered contrasts --------------------------------------------------------------------
# structure per (gene set, truth): regulators and the label matrix on the TF-node candidate set of score_truth (self-pairs removed)
struct <- function(g, gt) {
  a <- match(toupper(gt$Gene1), toupper(g)); b <- match(toupper(gt$Gene2), toupper(g)); ok <- !is.na(a) & !is.na(b) & a != b
  regs <- sort(unique(a[ok])); L <- matrix(0L, length(regs), length(g)); L[cbind(match(a[ok], regs), b[ok])] <- 1L
  K <- matrix(TRUE, length(regs), length(g)); K[cbind(seq_along(regs), regs)] <- FALSE
  list(regs = regs, L = L, K = K)
}
ap_of <- function(S, st, sel) { sc <- S[st$regs, , drop = FALSE][sel, , drop = FALSE]; lb <- st$L[sel, , drop = FALSE]; kp <- st$K[sel, , drop = FALSE]
  ap_fast(sc[kp], lb[kp]) }
PAIRS <- list(list("CLR-NS vs scPRINT R0", "A", "cmp|clr_ns|sym", "scprint_medium_A|R0|T"),
              list("CLR-NS vs scGPT attention", "S", "cmp|clr_ns|sym", "scgpt_A|attn_last_mean|T"),
              list("CLR-NS vs scGPT token cosine", "S", "cmp|clr_ns|sym", "scgpt_A|tok_cos|sym"),
              list("CLR-NS vs scGPT hidden cosine", "S", "cmp|clr_ns|sym", "scgpt_A|hid_cos|sym"),
              list("scPRINT R0 vs scGPT attention (three-way set)", "AS", "scprint_medium_A|R0|T", "scgpt_A|attn_last_mean|T"),
              list("scPRINT R0 vs scGPT token cosine (three-way set)", "AS", "scprint_medium_A|R0|T", "scgpt_A|tok_cos|sym"),
              list("scPRINT R0 vs scGPT hidden cosine (three-way set)", "AS", "scprint_medium_A|R0|T", "scgpt_A|hid_cos|sym"))
set.seed(sopt$seed)
boot_rows <- list(); draws <- list()
for (tn in names(TRUTHS)) {
  stA <- lapply(unique(vapply(PAIRS, `[[`, "", 2)), function(key) struct(GS[[key]], TRUTHS[[tn]])); names(stA) <- unique(vapply(PAIRS, `[[`, "", 2))
  for (pr in PAIRS) {
    key <- pr[[2]]; st <- stA[[key]]; m <- length(st$regs); if (m < 5L) { say("  ", tn, ": too few regulators for ", pr[[1]]); next }
    sels <- c(list(seq_len(m)), lapply(seq_len(sopt$boot), function(b) sample.int(m, m, replace = TRUE)))   # [[1]] = the full sample
    f <- function(sel) mean(vapply(SUBS, function(k) log2(ap_of(nets[[key]][[pr[[3]]]][[k]], st, sel) / ap_of(nets[[key]][[pr[[4]]]][[k]], st, sel)), 0))
    v <- unlist(parallel::mclapply(sels, f, mc.cores = sopt$workers))
    d_score <- NA_real_
    rr <- function(lab) { p <- strsplit(lab, "|", fixed = TRUE)[[1]]
      x <- res[res$truth == tn & res$gset == key & res$source == p[1] & res$readout == p[2] & res$orient == p[3], ]; x[order(x$subsample), "auprc_ratio"] }
    d_score <- mean(log2(rr(pr[[3]]) / rr(pr[[4]])))
    boot_rows[[length(boot_rows) + 1L]] <- data.frame(contrast = pr[[1]], truth = tn, gset = key, n_regulators = m, d = v[1], d_score_truth = d_score,
      lo = unname(stats::quantile(v[-1], 0.025)), hi = unname(stats::quantile(v[-1], 0.975)), B = sopt$boot)
    draws[[paste(pr[[1]], tn)]] <- v[-1]
    say(sprintf("  %-52s %-12s regs %3d  d = %+.3f [%+.3f, %+.3f]  (score_truth check %+.3f)", pr[[1]], tn, m, v[1],
                stats::quantile(v[-1], 0.025), stats::quantile(v[-1], 0.975), d_score))
    stopifnot(abs(v[1] - d_score) < 1e-8)
  }
}
bt <- do.call(rbind, boot_rows)
# mean over the three truths (same resample in every truth is not available across truths: truths have different regulators; the
# truth-level resamples are independent, so the pooled interval is the percentile interval of the mean of the truth-level resample draws)
pool <- do.call(rbind, lapply(unique(bt$contrast), function(cn) {
  nm <- paste(cn, names(TRUTHS)); nm <- nm[nm %in% names(draws)]
  if (length(nm) < length(TRUTHS)) return(NULL)
  mm <- Reduce(`+`, draws[nm]) / length(nm); d <- mean(bt$d[bt$contrast == cn])
  data.frame(contrast = cn, truth = "mean of the three truths", gset = bt$gset[bt$contrast == cn][1], n_regulators = NA, d = d,
             d_score_truth = mean(bt$d_score_truth[bt$contrast == cn]), lo = unname(stats::quantile(mm, 0.025)), hi = unname(stats::quantile(mm, 0.975)), B = sopt$boot) }))
bt <- rbind(bt, pool); utils::write.csv(bt, file.path(sopt$out, "bootstrap_contrasts.csv"), row.names = FALSE)
for (k in seq_len(nrow(pool))) say(sprintf("  POOLED %-52s d = %+.3f [%+.3f, %+.3f]", pool$contrast[k], pool$d[k], pool$lo[k], pool$hi[k]))

## ---- descriptive point estimates for every other readout (mean over subsamples; no interval) ---------------------------------------
cmpA <- res[res$source == "cmp" & res$gset == "A" & res$readout == "clr_ns", c("subsample", "truth", "auprc_ratio", "epr")]
names(cmpA)[3:4] <- c("ap_ns", "epr_ns")
pt <- merge(res[res$gset == "A" & res$source != "cmp", ], cmpA, by = c("subsample", "truth"))
pt$d_ap <- log2(pt$ap_ns / pt$auprc_ratio); pt$d_epr <- log2(pmax(pt$epr_ns, 1e-9) / pmax(pt$epr, 1e-9))
agg <- stats::aggregate(cbind(auprc_ratio, epr, d_ap, d_epr) ~ truth + source + readout + orient, pt, mean)
utils::write.csv(agg, file.path(sopt$out, "readouts_vs_clrns_mean_over_subsamples.csv"), row.names = FALSE)
key <- stats::aggregate(cbind(auprc_ratio, epr) ~ truth + gset + source + readout + orient, res[res$source == "cmp" | (res$source == "scprint_medium_A" & res$readout == "R0" & res$orient == "T"), ], mean)
utils::write.csv(key, file.path(sopt$out, "key_methods_mean_over_subsamples.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), "", paste("seed", sopt$seed)), file.path(sopt$out, "sessionInfo.txt"))
say("done")
