#!/usr/bin/env Rscript
# scGPT extension (docs/plans/scgpt_extension_plan.md with Revision 1): scoring of CLR-NS, tcaCLR and tcaCLR2 (and |r|, CLR) against the scGPT whole-human readouts, per dataset, on the
# scGPT gene set (the G-intersection of the plan = all scGPT-covered genes) and the same cells. One loop over datasets (hESC, hHep, K562, HepG2, RPE1, Jurkat; a dataset is a key of
# EXT_REG in ext_lib.R, so adding a dataset is a registry entry and a --datasets argument).
#
# Truths (Revision 1): three PRIMARY truths, equal weight: DoRothEA A-B and CollecTRI (directed TF-node score_truth of beeline_baselines.R, as the earlier scFM arms) and Complex Portal
# human complexes (pooled co-membership AUPR with ap_fast on the pair universe of the primary variant `main`, undirected, built on the unit's gene vector). Secondary, descriptive:
# knockdown-response modules (bulk file in the registry; none for Jurkat). Orientation: scGPT attention S = t(A) for the directed truths, max(A, t(A)) for the undirected complex
# pairs (the cosines are symmetric); the CLR-family scores are symmetric.
# Methods (rows `method`): abs_pearson, clr, clr_ns, tcaCLR (t = 9), tcaCLR2 (t2 from the permuted-data edge rule recomputed on each unit's matrix), scgpt_attn, scgpt_tok, scgpt_hid, scgpt_tokln.
#
# Gates (stop on failure, before anything is written): per unit, tcaCLR at t = 9 through gl_make_P / sk_a3r equals the independent rz_operator path (sk_layer) to 1e-12, gl_make_P equals
# rz_lazy_P(rz_operator) on the first permutation to 1e-12; --gate_stored (default for BEELINE datasets): tcaCLR and tcaCLR2 on the full BEELINE case matrix reproduce the stored scoreboard
# values and edge files to 1e-9; --methods legacy (K562 check only): |r|, CLR, CLR-NS and the four scGPT readouts reproduce results/scfm_score_k562/metrics_long.csv to 1e-9.
# No stored K562 tcaCLR values on the CLR-NS operator of the 1,244-gene sets exist (the stored K562 tcaCLR is the MI-CLR operator on 5,671 genes), so for the single-cell datasets the t = 9
# gate is the independent rz_operator implementation; the stored-value gate is the BEELINE one.
#
# Truth files are opened only with --truth real. --truth synthetic (smoke test) plants all truths from the networks of unit 1 and reads no truth file (complex structure is not built).
# Usage (repo root), nice 19:
#   Rscript analysis/scgpt_ext/score_ext.R --datasets HepG2 --out results/scgpt_ext [--truth real|synthetic] [--methods all|legacy] [--units 1,2] [--prep_dir DIR] [--exports DIR]
#       [--gate_stored hESC_500,hHep_500|none] [--threads 4]
# --prep_dir / --exports override the registry for ONE dataset (smoke tests with scratch inputs).
args0 <- commandArgs(trailingOnly = TRUE)
sopt <- list(datasets = "hESC,hHep,K562,HepG2", out = "results/scgpt_ext", truth = "real", methods = "all", units = "all", prep_dir = "", exports = "", gate_stored = "default", threads = 4L, seed = 20261021L, expr = "auto")
i <- 1L
while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(sopt))
  sopt[[k]] <- if (k %in% c("threads", "seed")) as.integer(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
DSN <- strsplit(sopt$datasets, ",")[[1]]
suppressPackageStartupMessages({ source("analysis/tcaclr2/edge_depth.R"); source("analysis/tcaclr2/tsb_lib.R") })       # adaptive_lib, grid_lib, tc_lib (TC), tsb_walk
# beeline.R set-up and helpers (everything before its main loop), score_truth and case_X of beeline_baselines.R: as analysis/tcaclr2/score_scoreboard.R --family beeline
args <- c("--out", sopt$out, "--threads", as.character(sopt$threads))
src0 <- readLines("analysis/beeline.R"); cut_at <- grep("^## ---- main loop", src0)[1]
for (e in parse(text = src0[seq_len(cut_at - 1L)])) { d1 <- deparse(e)[1]
  if (identical(d1, "args <- commandArgs(trailingOnly = TRUE)")) next
  if (startsWith(d1, "script_dir <-")) { script_dir <- "analysis"; next }
  eval(e, globalenv()) }
bsrc <- readLines("analysis/beeline_baselines.R"); a <- grep("^score_truth <- function", bsrc)[1]; z <- grep("^say\\(\"beeline_baselines: stage", bsrc)[1] - 1L
eval(parse(text = bsrc[a:z]), globalenv())
for (e in parse("analysis/beeline_baselines.R")) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] == "case_X") eval(e, globalenv())
source("analysis/scoreboard/sb_lib.R"); source("analysis/explore_stack/stack_lib.R"); source("analysis/robust_z/robust_z_lib.R"); source("analysis/human_truths.R")
source("analysis/scgpt_ext/ext_lib.R")
dir.create(sopt$out, recursive = TRUE, showWarnings = FALSE)
say <- tsb_say(file.path(sopt$out, sprintf("score_ext_%s.log", gsub(",", "_", sopt$datasets))))
say("options: ", paste(names(sopt), unlist(sopt), sep = "=", collapse = " "))
stopifnot(all(DSN %in% names(EXT_REG)), sopt$methods %in% c("all", "legacy"), sopt$truth %in% c("real", "synthetic"))
if ((nzchar(sopt$prep_dir) || nzchar(sopt$exports)) && length(DSN) != 1L) stop("--prep_dir / --exports need exactly one dataset")
LEGACY_ONLY <- sopt$methods == "legacy"
T_START <- Sys.time()

## ---- gates that do not need a dataset -----------------------------------------------------------------------------------------------
gates <- list()
add_gate <- function(g) { gates[[length(gates) + 1L]] <<- g }
stored_cases <- if (sopt$gate_stored == "default") unlist(lapply(DSN[vapply(DSN, function(d) EXT_REG[[d]]$type == "beeline", NA)], function(d) paste0(d, "_500"))) else
                if (sopt$gate_stored == "none") character(0) else strsplit(sopt$gate_stored, ",")[[1]]
for (cc in stored_cases) { g <- ext_gate_stored(cc, say); add_gate(g); say(sprintf("stored-value gate %s: %d checks, max %.3g, %s", cc, nrow(g), max(g$value), if (all(g$pass)) "pass" else "FAIL")) }
if (length(gates) && !all(do.call(rbind, gates)$pass)) stop("STORED-VALUE GATE FAILED; nothing scored")

## ---- datasets ----------------------------------------------------------------------------------------------------------------------
cat_of <- function(h5, grp, field) { codes <- rhdf5::h5read(h5, paste0(grp, "/", field)); cats <- rhdf5::h5read(h5, paste0(grp, "/__categories/", field)); as.character(cats)[codes + 1L] }
for (dn in DSN) {
  R <- EXT_REG[[dn]]; od <- file.path(sopt$out, dn); dir.create(od, recursive = TRUE, showWarnings = FALSE); t_ds <- Sys.time()
  exports <- path.expand(if (nzchar(sopt$exports)) sopt$exports else R$exports); prep <- if (nzchar(sopt$prep_dir)) sopt$prep_dir else R$dir
  say(sprintf("== %s (%s): exports %s", dn, R$type, exports))
  if (R$type == "sc") {
    gs <- utils::read.delim(file.path(prep, "geneset.tsv"), stringsAsFactors = FALSE); sub <- utils::read.delim(file.path(prep, "subsamples.tsv"), stringsAsFactors = FALSE)
    metatxt <- readLines(file.path(prep, "geneset_meta.json"), warn = FALSE)
    MED <- as.numeric(sub("[^0-9.]*([0-9.]+).*", "\\1", grep("median_control_umi", metatxt, value = TRUE)))
    NGS <- as.integer(sub("[^0-9]*([0-9]+).*", "\\1", grep("gene_set_size", metatxt, value = TRUE))); stopifnot(is.finite(MED), nrow(gs) == NGS)
    gs$is_tf_part <- gs$is_tf_part %in% c("True", "TRUE", TRUE)
    EXPR <- if (sopt$expr == "auto") (if (file.exists(R$h5)) "h5" else "files") else sopt$expr; stopifnot(EXPR %in% c("h5", "files"))     # files: after the raw h5ad was deleted
    if (EXPR == "h5") { vg <- cat_of(R$h5, "var", "gene_name"); gi <- match(gs$symbol, vg); stopifnot(!anyNA(gi), !anyDuplicated(gi)); bar <- as.character(rhdf5::h5read(R$h5, "obs/cell_barcode")) }
    else tots <- utils::read.delim(file.path(prep, "cell_totals.tsv"), stringsAsFactors = FALSE)
    say("expression source: ", EXPR)
    UNITS <- sort(unique(sub$subsample)); unit_case <- function(k) paste0(R$prefix, k)
    get_X <- function(k, genes) { rows <- sub$file_row[sub$subsample == k]
      if (EXPR == "h5") {
        C <- rhdf5::h5read(R$h5, "X", index = list(NULL, rows + 1L)); storage.mode(C) <- "double"                  # file_row is 0-based
        stopifnot(identical(bar[rows + 1L], sub$cell_barcode[sub$subsample == k]))
        L <- log1p(t(t(C) / colSums(C)) * MED)[gi, , drop = FALSE]; rownames(L) <- gs$symbol
      } else {
        tk <- tots[tots$subsample == k, ]; stopifnot(identical(tk$cell_barcode, sub$cell_barcode[sub$subsample == k]), identical(tk$file_row, rows))
        C <- matrix(readBin(file.path(prep, "expr", sprintf("%s%d.counts.f32", R$prefix, k)), "numeric", n = nrow(gs) * length(rows), size = 4L), nrow(gs), length(rows))
        L <- log1p(t(t(C) / tk$total_umi) * MED); rownames(L) <- gs$symbol }
      L[genes, , drop = FALSE] }
    tf_symbols <- gs$symbol[gs$is_tf_part]
  } else {
    UNITS <- R$cases; unit_case <- function(k) k
    get_X <- function(k, genes) { X <- case_X(R$ds, as.integer(sub("^.*_", "", k))); stopifnot(all(genes %in% rownames(X))); X[genes, , drop = FALSE] }
    tf_symbols <- NULL                                         # set per unit below (BEELINE human TF list)
  }
  if (sopt$units != "all") UNITS <- UNITS[as.character(UNITS) %in% strsplit(sopt$units, ",")[[1]]]
  gS_of <- function(k) ext_rd_genes(exports, "scgpt_A", unit_case(k))

  ## truths (opened only here, and only with --truth real)
  TRUTHS <- list(); CXDIR <- file.path(od, "complex_truth"); cx_cache <- list()
  cx_for <- function(genes, id) { key <- paste(length(genes), digest::digest(genes), sep = "_"); if (is.null(cx_cache[[key]])) cx_cache[[key]] <<- ext_complex_truth(genes, CXDIR, id, say); cx_cache[[key]] }
  if (sopt$truth == "real") {
    dor <- utils::read.delim("data/human_truth/dorothea.tsv", stringsAsFactors = FALSE)
    lv <- vapply(strsplit(dor$dorothea_level, ";"), function(x) any(x %in% c("A", "B")), NA); dor <- dor[lv, ]
    col <- utils::read.delim("data/human_truth/collectri.tsv", stringsAsFactors = FALSE)
    TRUTHS$DoRothEA_AB <- data.frame(Gene1 = dor$source_genesymbol, Gene2 = dor$target_genesymbol)
    TRUTHS$CollecTRI <- data.frame(Gene1 = col$source_genesymbol, Gene2 = col$target_genesymbol)
    if (R$type == "sc" && !is.na(R$bulk) && file.exists(R$bulk)) {
      rs <- response_sets(gs$symbol, path = R$bulk)
      TRUTHS$KD_response <- do.call(rbind, lapply(names(rs$members), function(nm) data.frame(Gene1 = strsplit(nm, "_")[[1]][2], Gene2 = gs$symbol[rs$members[[nm]]]))) }
    say(sprintf("truths (real): %s", paste(sprintf("%s %d edges", names(TRUTHS), vapply(TRUTHS, nrow, 1L)), collapse = "; ")))
  }

  rows <- list(); edges <- list(); counts <- list(); unit_info <- list(); synth <- NULL
  add_rows <- function(meta, res) rows[[length(rows) + 1L]] <<- data.frame(dataset = dn, meta, t(res), row.names = NULL)
  for (k in UNITS) {
    t_u <- Sys.time(); case <- unit_case(k); gS <- gS_of(k); X <- get_X(k, gS); G <- length(gS)
    say(sprintf("%s unit %s: %d genes x %d cells", dn, k, G, ncol(X)))
    M <- list()                                                                          # method -> list(dir = directed-orientation matrix, sym = symmetric matrix)
    sym1 <- function(S) list(dir = ext_clean(S), sym = ext_clean(S))
    M$abs_pearson <- sym1(abs(stats::cor(t(X)))); f <- clr_fit(X); M$clr <- sym1(f$clr_scores); f$release(FALSE)
    Cn <- ns_clr(X); M$clr_ns <- sym1(Cn)
    if (!LEGACY_ONLY) {
      e <- ext_edge(X); r <- ext_tcaclr(Cn, e$t2, gate = TRUE); edges[[length(edges) + 1L]] <- data.frame(dataset = dn, unit = k, e, gate_rz = r$gate_rz)
      stopifnot(r$gate_rz < 1e-12)
      say(sprintf("  edge %.5f (%.5f-%.5f), t2 = %d, tcaCLR t=9 vs rz_operator path max diff %.2g, gl_make_P vs rz_lazy_P %.2g", e$edge, e$edge_min, e$edge_max, e$t2, r$gate_rz, e$op_diff))
      M$tcaCLR <- sym1(r$tcaCLR); M$tcaCLR2 <- sym1(r$tcaCLR2); rm(r)
    }
    for (nm in names(EXT_SCGPT)) { A <- ext_rd_mat(exports, "scgpt_A", case, EXT_SCGPT[[nm]]); stopifnot(identical(rownames(A), gS))
      M[[nm]] <- if (nm == "scgpt_attn") list(dir = ext_clean(t(A)), sym = ext_clean(pmax(A, t(A)))) else sym1(A) }
    rm(Cn)
    ## truths of this unit
    if (sopt$truth == "synthetic") {      # planted from the networks of the first unit; never from a truth file
      if (is.null(synth)) {
        set.seed(7L); tfs <- if (!is.null(tf_symbols)) intersect(gS, tf_symbols) else intersect(gS, toupper(utils::read.csv("data/beeline/human-tfs.csv")$TF))
        tfi <- sort(sample(match(tfs, gS), min(60L, length(tfs)))); plant <- function(sig) do.call(rbind, lapply(tfi, function(r) { s <- rank(sig[r, ]) / G + stats::rnorm(G, sd = 0.3); s[r] <- -Inf
          data.frame(Gene1 = gS[r], Gene2 = gS[order(s, decreasing = TRUE)[1:25]]) }))
        sA <- M$scgpt_attn$dir; sN <- M$clr_ns$dir
        synth <- list(dir = list(DoRothEA_AB = plant(sN), CollecTRI = plant(sA), KD_response = plant(sA / stats::sd(sA) + sN / stats::sd(sN))), gS = gS,
                      cx = ext_plant_complexes(gS, if (!LEGACY_ONLY) M$tcaCLR$sym else sN))
        say("truths (SYNTHETIC, smoke test): DoRothEA_AB <- CLR-NS, CollecTRI <- scGPT attention, KD_response <- mixture, Complexes <- tcaCLR neighbours") }
      stopifnot(identical(synth$gS, gS)); dirT <- synth$dir; cxp <- synth$cx
    } else {
      dirT <- TRUTHS; cxp <- if (LEGACY_ONLY) list(sets = list(), label = integer(0)) else cx_for(gS, if (R$type == "sc") dn else case)   # the legacy check never touches complexes
    }
    for (mn in names(M)) {
      for (tn in names(dirT)) { x <- score_truth(M[[mn]]$dir, gS, NULL, dirT[[tn]]); if (is.null(x)) next
        add_rows(data.frame(unit = k, truth = tn, gset = "S", G = G, method = mn, orient = if (mn == "scgpt_attn") "T" else "sym"), x) }
      if (!LEGACY_ONLY) add_rows(data.frame(unit = k, truth = "Complexes", gset = "S", G = G, method = mn, orient = if (mn == "scgpt_attn") "max" else "sym"), ext_score_complex(M[[mn]]$sym, cxp))
    }
    counts[[length(counts) + 1L]] <- data.frame(dataset = dn, unit = k, G = G, cells = ncol(X), complexes = length(cxp$sets), complex_pairs = length(cxp$label), complex_positive_pairs = sum(cxp$label),
                                                complex_base_rate = if (length(cxp$label)) mean(cxp$label) else NA_real_)
    unit_info[[length(unit_info) + 1L]] <- data.frame(dataset = dn, unit = k, seconds = as.numeric(difftime(Sys.time(), t_u, units = "secs")), peak_rss_main_GB = tsb_hwm())
    say(sprintf("  unit %s done: %.1f s, main-process peak RSS %.2f GB", k, unit_info[[length(unit_info)]]$seconds, tsb_hwm())); rm(M, X); invisible(gc(verbose = FALSE))
  }
  res <- do.call(rbind, rows)
  ## legacy gate (K562): |r|, CLR, CLR-NS and the four scGPT readouts reproduce results/scfm_score_k562/metrics_long.csv on gene set S
  if (LEGACY_ONLY && sopt$truth == "real" && dn == "K562") {
    lg <- utils::read.csv("results/scfm_score_k562/metrics_long.csv", stringsAsFactors = FALSE); lg <- lg[lg$gset == "S", ]
    lg$method <- ifelse(lg$source == "cmp", lg$readout, paste0("scgpt_", c(attn_last_mean = "attn", tok_cos = "tok", hid_cos = "hid", tok_cos_ln = "tokln")[lg$readout]))
    j <- match(paste(res$unit, res$truth, res$method), paste(lg$subsample, lg$truth, lg$method)); ok <- !is.na(j) & res$truth != "Complexes"
    d <- max(abs(res$auprc_ratio[ok] - lg$auprc_ratio[j[ok]]), abs(res$epr[ok] - lg$epr[j[ok]]))
    add_gate(data.frame(gate = sprintf("K562 legacy rows (%d) vs results/scfm_score_k562/metrics_long.csv: AUPR ratio and EPR", sum(ok)), value = d, tol = 1e-9, pass = d <= 1e-9)); say(sprintf("legacy gate: %d rows, max abs diff %.3g", sum(ok), d))
    if (d > 1e-9) stop("LEGACY GATE FAILED") }
  utils::write.csv(res, file.path(od, if (LEGACY_ONLY) "metrics_ext_legacy_check.csv" else "metrics_ext.csv"), row.names = FALSE)
  if (length(edges)) utils::write.csv(do.call(rbind, edges), file.path(od, "tcaclr2_edge_t2.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, counts), file.path(od, "complex_truth_counts.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, unit_info), file.path(od, "unit_timing.csv"), row.names = FALSE)
  if (length(gates)) utils::write.csv(do.call(rbind, gates), file.path(od, "gates.csv"), row.names = FALSE)
  say(sprintf("%s done: %d metric rows, %.1f s", dn, nrow(res), as.numeric(difftime(Sys.time(), t_ds, units = "secs"))))
}
writeLines(c(capture.output(sessionInfo()), "", paste("seed", sopt$seed)), file.path(sopt$out, "sessionInfo.txt"))
say(sprintf("done: %.1f s, main-process peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T_START, units = "secs")), tsb_hwm()))
