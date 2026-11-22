#!/usr/bin/env Rscript
# Spectral-embedding exploratory pass (docs/plans/spectral_embedding_pass_plan.md): scores the 30 variants of sp_lib.R on every unit of the development datasets with the metrics code of
# analysis/scgpt_ext/score_ext.R (score_truth for DoRothEA A-B, CollecTRI and KD response; ext_score_complex / ap_fast for Complex Portal complexes; same units, same scGPT gene set).
# Output (never results/scgpt_ext): <out>/<dataset>/{metrics_spectral.csv, variant_info.csv, eig_info.csv, ring_planes.csv, gates.csv, unit_timing.csv, complex_truth/}.
# Truth files are opened only with --truth real. --truth synthetic (smoke test) plants every truth from the unit's own matrices and reads no truth file.
# Usage (repo root, nice 19):
#   Rscript analysis/spectral_pass/sp_score.R --datasets HepG2 --out results/spectral_pass [--truth real|synthetic] [--units 1,2] [--only cos,mclust,...] [--krules E,5,10,20,50] [--variants all|name1,name2,...] [--expr auto|files|h5] [--threads 4]
# Jurkat is refused unless --heldout YES is passed (the lead does that once, after freezing the variants).
a0 <- commandArgs(trailingOnly = TRUE)
so <- list(datasets = "hESC,hHep,K562,HepG2,RPE1", out = "results/spectral_pass", truth = "real", units = "all", only = "cos,cos_clr2,mclust,hddc,pgmm,ring,rmix",
           krules = "E,5,10,20,50", variants = "all", expr = "auto", threads = 4L, heldout = "no")
i <- 1L; while (i <= length(a0)) { k <- sub("^--", "", a0[i]); stopifnot(k %in% names(so)); so[[k]] <- if (k == "threads") as.integer(a0[i + 1L]) else a0[i + 1L]; i <- i + 2L }
SP_OUT <- so$out; SP_THREADS <- so$threads
HELD <- so$heldout == "YES"
source("analysis/spectral_pass/sp_lib.R")
if (HELD) SP_DEV <- c("Jurkat")                       # only the sealed dataset, only on explicit request
VARS <- if (so$variants == "all") NULL else strsplit(so$variants, ",")[[1]]
if (HELD && is.null(VARS)) stop("--heldout YES requires --variants <frozen list>")
DSN <- strsplit(so$datasets, ",")[[1]]; ONLY <- strsplit(so$only, ",")[[1]]; KR <- strsplit(so$krules, ",")[[1]]
stopifnot(all(DSN %in% names(EXT_REG)), all(DSN %in% SP_DEV), so$truth %in% c("real", "synthetic"), all(ONLY %in% SP_READOUTS), all(KR %in% SP_KRULES))
dir.create(so$out, recursive = TRUE, showWarnings = FALSE)
say <- tsb_say(file.path(so$out, sprintf("sp_score_%s.log", gsub(",", "_", so$datasets))))
say("options: ", paste(names(so), unlist(so), sep = "=", collapse = " ")); T0 <- Sys.time()
for (dn in DSN) {
  D <- sp_dataset(dn, so$expr); od <- file.path(so$out, dn); dir.create(od, recursive = TRUE, showWarnings = FALSE); t_ds <- Sys.time(); R <- D$R
  UNITS <- D$units; if (so$units != "all") UNITS <- UNITS[as.character(UNITS) %in% strsplit(so$units, ",")[[1]]]
  say(sprintf("== %s (%s), expression source %s, units %s", dn, R$type, D$expr, paste(UNITS, collapse = ",")))
  ## truths: opened only here and only with --truth real (copy of the block of score_ext.R)
  TRUTHS <- list(); CXDIR <- file.path(od, "complex_truth"); cx_cache <- list()
  cx_for <- function(genes, id) { key <- paste(length(genes), digest::digest(genes), sep = "_"); if (is.null(cx_cache[[key]])) cx_cache[[key]] <<- ext_complex_truth(genes, CXDIR, id, say); cx_cache[[key]] }
  if (so$truth == "real") {
    dor <- utils::read.delim("data/human_truth/dorothea.tsv", stringsAsFactors = FALSE)
    lv <- vapply(strsplit(dor$dorothea_level, ";"), function(x) any(x %in% c("A", "B")), NA); dor <- dor[lv, ]
    col <- utils::read.delim("data/human_truth/collectri.tsv", stringsAsFactors = FALSE)
    TRUTHS$DoRothEA_AB <- data.frame(Gene1 = dor$source_genesymbol, Gene2 = dor$target_genesymbol)
    TRUTHS$CollecTRI <- data.frame(Gene1 = col$source_genesymbol, Gene2 = col$target_genesymbol)
    if (R$type == "sc" && !is.na(R$bulk) && file.exists(R$bulk)) {
      gs <- utils::read.delim(file.path(R$dir, "geneset.tsv"), stringsAsFactors = FALSE)
      rs <- response_sets(gs$symbol, path = R$bulk)
      TRUTHS$KD_response <- do.call(rbind, lapply(names(rs$members), function(nm) data.frame(Gene1 = strsplit(nm, "_")[[1]][2], Gene2 = gs$symbol[rs$members[[nm]]]))) }
    say(sprintf("truths (real): %s", paste(sprintf("%s %d edges", names(TRUTHS), vapply(TRUTHS, nrow, 1L)), collapse = "; ")))
  }
  rows <- list(); vinfo <- list(); einfo <- list(); rinfo <- list(); gates <- list(); utime <- list(); synth <- NULL
  for (k in UNITS) {
    t_u <- Sys.time(); u <- D$unit_X(k); X <- u$X; gS <- u$genes; G <- length(gS); say(sprintf("%s unit %s: %d genes x %d cells", dn, k, G, ncol(X)))
    e <- ext_edge(X); g <- sp_gate_edge(dn, k, e)
    gates[[length(gates) + 1L]] <- data.frame(dataset = dn, unit = k, gate = "permuted-data edge and t2 equal results/scgpt_ext/<dataset>/tcaclr2_edge_t2.csv", value = g$diff, tol = 1e-9, pass = g$pass)
    say(sprintf("  edge %.6f (stored %.6f), t2 %d (stored %d): %s", e$edge, g$stored_edge, e$t2, g$stored_t2, if (is.na(g$pass)) "gate skipped (no stored edge file)" else if (g$pass) "gate pass" else "GATE FAIL")); if (isFALSE(g$pass)) stop("EDGE GATE FAILED; nothing scored")
    Cn <- ns_clr(X); emb <- sp_embed(Cn, e$edge, e$t2); sN <- Cn
    say(sprintf("  eigenvalues above the edge among the first %d non-trivial: %d (k_edge = %d)", SP$neig, emb$n_above, emb$k_edge))
    einfo[[length(einfo) + 1L]] <- data.frame(dataset = dn, unit = k, G = G, cells = ncol(X), edge = e$edge, t2 = e$t2, n_above = emb$n_above, k_edge = emb$k_edge, t(stats::setNames(emb$lam, paste0("lam", seq_along(emb$lam)))))
    sc <- sp_unit_scores(emb, say, only = ONLY, krules = KR, variants = VARS)
    vinfo[[length(vinfo) + 1L]] <- data.frame(dataset = dn, unit = k, sc$info)
    if (!is.null(sc$ring_planes)) {
      nullh1 <- if (SP$tda) sp_tda_null(X, e$t2, e$edge, length(sc$ring_planes)) else rep(NA_real_, length(sc$ring_planes))
      rinfo[[length(rinfo) + 1L]] <- do.call(rbind, lapply(seq_along(sc$ring_planes), function(j) { p <- sc$ring_planes[[j]]
        data.frame(dataset = dn, unit = k, plane = j, A_accepted = p$ring$accepted, A_H = p$ring$H, A_dbic = p$ring$bic_ring - p$ring$bic_cluster, A_mu = paste(round(p$ring$mu, 5), collapse = ";"), A_s = paste(round(p$ring$s, 5), collapse = ";"),
                   B_accepted = p$rmix$accepted, B_H = p$rmix$H, B_dbic = p$rmix$bic_ring - p$rmix$bic_cluster, B_mu = paste(round(p$rmix$mu, 5), collapse = ";"), B_s = paste(round(p$rmix$s, 5), collapse = ";"),
                   B_centre = paste(round(p$rmix$centre, 5), collapse = ";"), h1_max_persistence = p$h1, h1_null_plane = nullh1[j], h1_null_max_all_planes = max(nullh1), h1_exceeds_null_max = p$h1 > max(nullh1)) })) }
    ## truths of this unit
    if (so$truth == "synthetic") {
      if (is.null(synth)) { set.seed(7L); tfs <- if (!is.null(D$tf_symbols)) intersect(gS, D$tf_symbols) else intersect(gS, toupper(utils::read.csv("data/beeline/human-tfs.csv")$TF))
        tfi <- sort(sample(match(tfs, gS), min(60L, length(tfs)))); plant <- function(sig) do.call(rbind, lapply(tfi, function(r) { s <- rank(sig[r, ]) / G + stats::rnorm(G, sd = 0.3); s[r] <- -Inf; data.frame(Gene1 = gS[r], Gene2 = gS[order(s, decreasing = TRUE)[1:25]]) }))
        sE <- sp_cos(emb$Y(10)); synth <- list(dir = list(DoRothEA_AB = plant(sN), CollecTRI = plant(sE), KD_response = plant(sE / stats::sd(sE) + sN / stats::sd(sN))), gS = gS, cx = ext_plant_complexes(gS, sN))
        say("truths (SYNTHETIC, smoke test): DoRothEA_AB <- CLR-NS, CollecTRI <- spectral cosine k = 10, KD_response <- mixture, Complexes <- CLR-NS neighbours") }
      stopifnot(identical(synth$gS, gS)); dirT <- synth$dir; cxp <- synth$cx
    } else { dirT <- TRUTHS; cxp <- cx_for(gS, if (R$type == "sc") dn else D$unit_case(k)) }
    for (vn in names(sc$S)) { S <- ext_clean(sc$S[[vn]]); meta <- function(tn) data.frame(unit = k, truth = tn, gset = "S", G = G, method = vn, orient = "sym")
      for (tn in names(dirT)) { x <- score_truth(S, gS, NULL, dirT[[tn]]); if (is.null(x)) next; rows[[length(rows) + 1L]] <- data.frame(dataset = dn, meta(tn), t(x), row.names = NULL) }
      rows[[length(rows) + 1L]] <- data.frame(dataset = dn, meta("Complexes"), t(ext_score_complex(S, cxp)), row.names = NULL) }
    utime[[length(utime) + 1L]] <- data.frame(dataset = dn, unit = k, seconds = as.numeric(difftime(Sys.time(), t_u, units = "secs")), peak_rss_main_GB = tsb_hwm())
    say(sprintf("  unit %s done: %.1f s, main-process peak RSS %.2f GB", k, utime[[length(utime)]]$seconds, tsb_hwm())); rm(sc, emb, X, Cn); invisible(gc(verbose = FALSE))
  }
  wr <- function(x, f) if (length(x)) utils::write.csv(dplyr_bind(x), file.path(od, f), row.names = FALSE)
  wr(rows, "metrics_spectral.csv"); wr(vinfo, "variant_info.csv"); wr(einfo, "eig_info.csv"); wr(rinfo, "ring_planes.csv"); wr(gates, "gates.csv"); wr(utime, "unit_timing.csv")
  say(sprintf("%s done: %d metric rows, %.1f s", dn, length(rows), as.numeric(difftime(Sys.time(), t_ds, units = "secs"))))
}
writeLines(c(capture.output(sessionInfo()), "", paste("seed", SP$seed)), file.path(so$out, sprintf("sessionInfo_%s.txt", gsub(",", "_", so$datasets))))
say(sprintf("done: %.1f s, main-process peak RSS %.2f GB", as.numeric(difftime(Sys.time(), T0, units = "secs")), tsb_hwm()))
